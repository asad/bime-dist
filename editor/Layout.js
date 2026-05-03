/**
 * Layout.js — 2D coordinate generation for molecular structures
 *
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Generates publication-quality 2D coordinates for Molecule objects using an
 * algorithm inspired by Helson's Structure Diagram Generation (SDG) review and
 * the standard structure diagram generation architecture:
 *
 *   1. SSSR ring perception (Horton / edge-vector Gaussian elimination)
 *   2. Ring-system fusion: regular-polygon placement, edge-reflection for
 *      fused rings, spiro-junction handling, bridged-ring templates
 *   3. 120-degree zigzag chain layout extending away from ring attachment
 *   4. Substituent fan placement in the open arc opposite existing bonds
 *   5. Iterative overlap / collision resolution (chain atoms preferred)
 *   6. Bond-length normalisation to BOND_LENGTH (30 px)
 *   7. Aromatic ring centre coordinates for circle rendering
 *
 * Handles molecules up to 200+ atoms efficiently (spatial grid for O(n)
 * collision detection instead of O(n^2)).
 */
(function (global) {
    'use strict';

    // =====================================================================
    // Constants
    // =====================================================================

    var BOND_LENGTH    = Molecule.BOND_LENGTH;          // 30 px
    var TWO_PI         = 2 * Math.PI;
    var DEG120         = TWO_PI / 3;                    // 120 degrees
    var DEG60          = Math.PI / 3;                   // 60 degrees
    var MIN_ATOM_DIST  = BOND_LENGTH * 0.55;            // collision threshold
    var GRID_CELL      = BOND_LENGTH * 1.2;             // spatial-hash cell size
    var MAX_OVERLAP_ITER = 80;                          // overlap-resolution passes
    var EPSILON        = 1e-6;
    // Macrocycle radius compression is now gradual (computed inline)

    // =====================================================================
    // Public API
    // =====================================================================

    var Layout = {};

    /**
     * Layout-time options. Mutate before calling Layout.layout().
     *
     *   enforceConventions (bool, default true)
     *     v1.5.2 chemistry-canonical post-passes:
     *       Step 1: bond-length normalisation toward BOND_LENGTH (relaxation)
     *       Step 2: aromatic ring orientation (one bond horizontal / hetero on top)
     *       Step 3: substituent radial-outward enforcement
     *       Step 4: fused-system longest-axis horizontal
     *       Step 5: sp2/sp3 bond-angle smoothing toward 120 deg
     *     Set to false for byte-identical v1.5.1 backward-compat layouts.
     *
     *   useMLDepict (bool, default false)
     *     v1.6.0 experimental opt-in ML residual refinement. Requires
     *     editor/MLDepict.js + editor/ml-depict-weights.json to be loaded.
     *     When true, after the rule-based pipeline (and any enforce-
     *     Conventions passes) the ML model emits a residual coord per
     *     atom and Layout blends it in by Layout.options.mlDepictWeight.
     *     Off by default since the model under-performs on BIME templates;
     *     opt-in for experimental refinement use cases.
     *
     *   mlDepictWeight (float in [0, 1], default 0.15)
     *     Blend factor for ML refinement: 0 = pure rule-based, 1 = pure ML.
     *     Only used when useMLDepict is true. Conservative 0.15 default
     *     keeps rule-based as the dominant signal.
     *
     * Determinism contract: every step iterates over molecule.atoms in id
     * order, no Math.random / no Date.now / no Object.keys ordering. Same
     * input SMILES => byte-identical coords across runs.
     */
    Layout.options = {
        enforceConventions: true,
        useMLDepict: false,
        mlDepictWeight: 0.15
    };

    /**
     * Compute 2D coordinates for every atom in `mol`.
     * Disconnected fragments are placed side-by-side.
     */
    Layout.layout = function (mol) {
        if (!mol || mol.atoms.length === 0) return;

        // Reset all atom coordinates before computing fresh layout
        for (var ri = 0; ri < mol.atoms.length; ri++) {
            mol.atoms[ri].x = 0;
            mol.atoms[ri].y = 0;
        }

        var components = mol.getComponents();
        var offsetX = 0;

        for (var ci = 0; ci < components.length; ci++) {
            var compAtomIds = components[ci];
            var compAtoms   = compAtomIds.map(function (id) { return mol.getAtom(id); });

            layoutComponent(mol, compAtomIds);

            // Shift component so it sits to the right of the previous one
            var bounds = atomBounds(compAtoms);
            var shiftX = offsetX - bounds.minX;
            var shiftY = -(bounds.minY + bounds.maxY) / 2;
            for (var i = 0; i < compAtoms.length; i++) {
                compAtoms[i].x += shiftX;
                compAtoms[i].y += shiftY;
            }
            bounds  = atomBounds(compAtoms);
            offsetX = bounds.maxX + BOND_LENGTH * 2;
        }
    };

    /**
     * Layout a single fragment (array of Atom objects).
     * Called from SmilesParser for each fragment.
     */
    Layout.layoutFragment = function (mol, fragAtoms) {
        if (!fragAtoms || fragAtoms.length === 0) return;
        var atomIds = fragAtoms.map(function (a) { return a.id; });
        layoutComponent(mol, atomIds);
    };

    // =====================================================================
    // Core layout for one connected component
    // =====================================================================

    function layoutComponent(mol, atomIds) {
        if (atomIds.length === 0) return;
        if (atomIds.length === 1) {
            var a = mol.getAtom(atomIds[0]);
            a.x = 0; a.y = 0;
            return;
        }

        // Build set for quick membership check
        var inComp = {};
        for (var i = 0; i < atomIds.length; i++) inComp[atomIds[i]] = true;

        // -- Step 1: Perceive rings ----------------------------------------
        // v1.4.1 (bug-fix #20): use the built-in DFS ring finder + GF(2) SSSR
        // selection unconditionally. The earlier "prefer SMSDRings" branch
        // was dead code (rings was always undefined when the conditional ran)
        // and was removed alongside the misleading comment. SMSDRings (Horton
        // SSSR) is still used for MCS/substructure work but produces ring
        // orderings that can confuse the coordinate generator, so layout
        // sticks with the in-place DFS finder.
        var compSet = {};
        for (var ci = 0; ci < atomIds.length; ci++) compSet[atomIds[ci]] = true;
        var rings;

        var allRings = mol.findRings(20);
        var candidateRings = [];
        for (var ri = 0; ri < allRings.length; ri++) {
            var mr = allRings[ri].atoms;
            var allInComp = true;
            for (var ai = 0; ai < mr.length; ai++) {
                if (!compSet[mr[ai]]) { allInComp = false; break; }
            }
            if (allInComp) candidateRings.push(mr);
        }

        var compEdges = 0;
        for (var ci2 = 0; ci2 < atomIds.length; ci2++) {
            var nbrs = mol.getNeighbors(atomIds[ci2]);
            for (var ni = 0; ni < nbrs.length; ni++) {
                if (compSet[nbrs[ni]] && nbrs[ni] > atomIds[ci2]) compEdges++;
            }
        }
        var targetSSSR = compEdges - atomIds.length + 1;

        candidateRings.sort(function(a, b) { return a.length - b.length; });
        if (targetSSSR > 0 && candidateRings.length > targetSSSR) {
            rings = selectIndependentRings(mol, candidateRings, atomIds, targetSSSR);
        } else {
            rings = candidateRings;
        }

        if (rings.length === 0) {
            rings = perceiveSSSR(mol, atomIds);
        }

        // -- Step 2: Build fused ring systems (union-find) ----------------
        var ringSystems = buildRingSystems(rings);

        // -- Step 2b: Template matching for known ring systems -----------
        var templatePlaced = matchRingSystemTemplates(mol, rings, ringSystems);

        // -- Step 3: Atom-to-ring maps ------------------------------------
        var ringAtomSet  = {};
        var atomToRings  = {};
        for (var ri = 0; ri < rings.length; ri++) {
            for (var ai = 0; ai < rings[ri].length; ai++) {
                var aid = rings[ri][ai];
                ringAtomSet[aid] = true;
                if (!atomToRings[aid]) atomToRings[aid] = [];
                atomToRings[aid].push(ri);
            }
        }

        var placed = {};

        // Merge template-placed atoms into placed set
        for (var tp in templatePlaced) {
            if (templatePlaced.hasOwnProperty(tp)) placed[tp] = true;
        }

        // -- Step 4a: Place ring systems ----
        // Place the largest ring system first at origin.
        // All other ring systems are deferred to layoutChains (Step 4c)
        // which positions them as substituents when the connecting chain
        // reaches them. This prevents unanchored rings from being placed
        // at arbitrary positions (fixes Taxol, biphenyl, atorvastatin).

        // Find the largest ring system (most atoms)
        var bestSysIdx = -1, bestSysSize = 0;
        for (var si = 0; si < ringSystems.length; si++) {
            var sysAtoms = collectSystemAtoms(rings, ringSystems[si]);
            var allPlaced = true;
            for (var sa = 0; sa < sysAtoms.length; sa++) {
                if (!templatePlaced[sysAtoms[sa]]) { allPlaced = false; break; }
            }
            if (allPlaced) continue;
            if (sysAtoms.length > bestSysSize) {
                bestSysSize = sysAtoms.length;
                bestSysIdx = si;
            }
        }
        if (bestSysIdx >= 0) {
            layoutRingSystem(mol, rings, ringSystems[bestSysIdx], placed, ringAtomSet);
        }

        // Build a set of ring atoms for deferred ring systems
        var deferredRingSystems = [];
        for (var si = 0; si < ringSystems.length; si++) {
            if (si === bestSysIdx) continue;
            var sysAtoms = collectSystemAtoms(rings, ringSystems[si]);
            var allPlaced = true;
            for (var sa = 0; sa < sysAtoms.length; sa++) {
                if (!placed[sysAtoms[sa]]) { allPlaced = false; break; }
            }
            if (!allPlaced) deferredRingSystems.push(si);
        }

        // -- Step 4b: If no rings, seed the first atom --------------------
        if (Object.keys(placed).length === 0) {
            var start = mol.getAtom(atomIds[0]);
            start.x = 0; start.y = 0;
            placed[atomIds[0]] = true;
        }

        // -- Step 4c: Place chains & substituents -------------------------
        // Pass deferred ring systems so chains can trigger their placement
        layoutChains(mol, atomIds, placed, ringAtomSet, atomToRings, rings,
                     deferredRingSystems, ringSystems);

        // -- Step 5: Bond-length normalisation ----------------------------
        normaliseBondLengths(mol, atomIds);

        // -- Step 6: Overlap resolution -----------------------------------
        resolveCollisions(mol, atomIds, ringAtomSet);

        // -- Step 7: Orientation optimization -----------------------------
        optimizeOrientation(mol, atomIds);

        // -- Step 8: Light force-directed refinement ----------------------
        refineLayout(mol, atomIds, ringAtomSet);

        // -- Step 9: Final collision cleanup after refinement -------------
        resolveCollisions(mol, atomIds, ringAtomSet);

        // -- Step 9b: v1.5.2 chemistry-canonical post-passes ---------------
        //   These run AFTER any template-matched layout AND after the
        //   generic pipeline. They produce textbook geometry: uniform bond
        //   lengths, hexagons with horizontal bonds, radial substituents,
        //   horizontal fused axes, smoothed sp2/sp3 angles.
        if (Layout.options && Layout.options.enforceConventions) {
            // Sanity gate: skip post-passes if base layout produced
            // non-finite or absurdly large coords (1000 * BOND_LENGTH).
            // Running gentle relaxation on those would not converge and
            // would amplify the existing bug rather than mask it.
            var skipConventions = false;
            var ABSURD = BOND_LENGTH * 1000;
            for (var skipI = 0; skipI < atomIds.length; skipI++) {
                var skipA = mol.getAtom(atomIds[skipI]);
                if (!skipA) continue;
                if (!isFinite(skipA.x) || !isFinite(skipA.y) ||
                    Math.abs(skipA.x) > ABSURD || Math.abs(skipA.y) > ABSURD) {
                    skipConventions = true; break;
                }
            }
            if (!skipConventions) {
                // Step 2/3/4 are gated to the GENERIC fallback path (no template):
                // when templates fired, they have already chosen the canonical
                // orientation, so we only run the universal Step 1 / Step 5
                // polish on top of the templated layout.
                var hadTemplate = false;
                for (var tk in templatePlaced) {
                    if (templatePlaced.hasOwnProperty(tk)) { hadTemplate = true; break; }
                }
                relaxBondLengths(mol, atomIds, ringAtomSet);   // Step 1 (universal)
                if (!hadTemplate) {
                    enforceRingOrientation(mol, rings, ringSystems, atomIds);     // 2
                    enforceSubstituentDirection(mol, rings, ringAtomSet, atomIds); // 3
                    enforceFusedAxisHorizontal(mol, rings, ringSystems);          // 4
                }
                smoothBondAngles(mol, atomIds, ringAtomSet);   // Step 5 (universal)
                relaxBondLengths(mol, atomIds, ringAtomSet);   // re-normalise after 2-5
                resolveCollisions(mol, atomIds, ringAtomSet);  // final sweep
            }
        }

        // -- Step 10: Count crossings. Only run heavy optimisation if needed.
        var hasCrossings = false;
        if (typeof SMSDLayout !== 'undefined') {
            for (var bi = 0; bi < mol.bonds.length && !hasCrossings; bi++) {
                var b1 = mol.bonds[bi];
                var ba1 = mol.getAtom(b1.atom1), ba2 = mol.getAtom(b1.atom2);
                if (!ba1 || !ba2) continue;
                for (var bj = bi + 1; bj < mol.bonds.length; bj++) {
                    var b2 = mol.bonds[bj];
                    if (b2.atom1 === b1.atom1 || b2.atom1 === b1.atom2 ||
                        b2.atom2 === b1.atom1 || b2.atom2 === b1.atom2) continue;
                    var ba3 = mol.getAtom(b2.atom1), ba4 = mol.getAtom(b2.atom2);
                    if (!ba3 || !ba4) continue;
                    var cd1 = (ba4.x-ba3.x)*(ba1.y-ba3.y)-(ba4.y-ba3.y)*(ba1.x-ba3.x);
                    var cd2 = (ba4.x-ba3.x)*(ba2.y-ba3.y)-(ba4.y-ba3.y)*(ba2.x-ba3.x);
                    var cd3 = (ba2.x-ba1.x)*(ba3.y-ba1.y)-(ba2.y-ba1.y)*(ba3.x-ba1.x);
                    var cd4 = (ba2.x-ba1.x)*(ba4.y-ba1.y)-(ba2.y-ba1.y)*(ba4.x-ba1.x);
                    if (((cd1>0&&cd2<0)||(cd1<0&&cd2>0))&&((cd3>0&&cd4<0)||(cd3<0&&cd4>0))) {
                        hasCrossings = true; break;
                    }
                }
            }
        }

        // -- Step 10b: Crossing reduction (only if crossings exist) --------
        if (hasCrossings && SMSDLayout.reduceCrossings) {
            SMSDLayout.reduceCrossings(mol, atomIds, ringAtomSet, rings);
        }

        // -- Step 11: Force-directed refinement (only if crossings exist) --
        if (hasCrossings && SMSDLayout.forceDirectedLayout) {
            SMSDLayout.forceDirectedLayout(mol, atomIds, ringAtomSet, rings);
        }

        // -- Step 12 (v1.6.0): ML residual refinement (opt-in) ------------
        // Tiny network predicts a small (dx, dy) correction relative to
        // the neighbour centroid; we blend it into the rule-based coords
        // by Layout.options.mlDepictWeight. Default behaviour is OFF.
        if (Layout.options && Layout.options.useMLDepict &&
            typeof MLDepict !== 'undefined' && MLDepict.ready &&
            MLDepict.ready()) {
            var alpha = Layout.options.mlDepictWeight;
            if (typeof alpha !== 'number' || alpha < 0) alpha = 0.15;
            if (alpha > 1) alpha = 1;
            // Only refine atoms in this component (build a sub-mol view).
            var subMol = {
                atoms: atomIds.map(function (id) { return mol.getAtom(id); }),
                bonds: mol.bonds,
                getAtom: function (id) { return mol.getAtom(id); },
                getNeighbors: function (id) { return mol.getNeighbors(id); },
                findRings: function (mx) { return mol.findRings ? mol.findRings(mx) : []; }
            };
            try { MLDepict.refineLayout(subMol, alpha); }
            catch (mle) { /* swallow — never break layout because of ML */ }
            // After ML refinement, run a single bond-length normalisation
            // pass to avoid the blend pulling bonds away from BOND_LENGTH.
            relaxBondLengths(mol, atomIds, ringAtomSet);
        }
    }

    // =====================================================================
    // SSSR Ring Perception  (Horton variant + edge-vector independence)
    // =====================================================================

    function perceiveSSSR(mol, atomIds) {
        var n = atomIds.length;
        if (n < 3) return [];

        var idToIdx = {}, idxToId = [];
        for (var i = 0; i < n; i++) { idToIdx[atomIds[i]] = i; idxToId[i] = atomIds[i]; }

        // Build adjacency (local indices)
        var adj = [];
        for (var i = 0; i < n; i++) adj[i] = [];
        var edgeCount = 0;

        for (var i = 0; i < n; i++) {
            var neighbors = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < neighbors.length; j++) {
                var ni = idToIdx[neighbors[j]];
                if (ni !== undefined && ni > i) {
                    adj[i].push(ni);
                    adj[ni].push(i);
                    edgeCount++;
                }
            }
        }
        for (var i = 0; i < n; i++) adj[i] = uniqueArray(adj[i]);

        var targetRingCount = edgeCount - n + 1;
        if (targetRingCount <= 0) return [];

        // Collect candidate rings using DFS back-edge detection
        // (Horton BFS has issues with even-cycle sizes; DFS is simpler and correct)
        var candidateRings = [];
        var seenKeys = {};

        for (var startV = 0; startV < n; startV++) {
            // DFS from startV looking for paths back to startV
            var stack = [[startV, [startV], {}]];
            stack[0][2][startV] = true;

            while (stack.length > 0) {
                var frame = stack.pop();
                var cur = frame[0], path = frame[1], vis = frame[2];
                if (path.length > 20) continue;

                var nbrs = adj[cur];
                for (var j = 0; j < nbrs.length; j++) {
                    var w = nbrs[j];
                    if (w === startV && path.length >= 3) {
                        // Found ring
                        var key = normalizeRing(path);
                        if (!seenKeys[key]) {
                            seenKeys[key] = true;
                            candidateRings.push({ ring: path.slice(), key: key, size: path.length });
                        }
                    } else if (!vis[w] && path.length < 20) {
                        var nv = {};
                        for (var vk in vis) nv[vk] = true;
                        nv[w] = true;
                        stack.push([w, path.concat(w), nv]);
                    }
                }
            }
        }

        candidateRings.sort(function (a, b) { return a.size - b.size; });
        var selected = selectSSSR(candidateRings, n, adj, targetRingCount);

        var result = [];
        for (var i = 0; i < selected.length; i++) {
            var ring = selected[i].ring;
            var atomRing = [];
            for (var j = 0; j < ring.length; j++) atomRing.push(idxToId[ring[j]]);
            result.push(atomRing);
        }
        return result;
    }

    // NOTE: tracePath and isSimpleRing were dead code (never called) — removed to reduce bundle size

    function normalizeRing(ring) {
        var minIdx = 0;
        for (var i = 1; i < ring.length; i++) {
            if (ring[i] < ring[minIdx]) minIdx = i;
        }
        var forward = [], backward = [];
        for (var i = 0; i < ring.length; i++) {
            forward.push(ring[(minIdx + i) % ring.length]);
            backward.push(ring[(minIdx - i + ring.length) % ring.length]);
        }
        var fk = forward.join(','), bk = backward.join(',');
        return fk < bk ? fk : bk;
    }

    /**
     * Select linearly independent rings from candidates using GF(2) elimination.
     * Works with atom IDs directly (unlike selectSSSR which uses local indices).
     */
    function selectIndependentRings(mol, candidateRings, atomIds, targetCount) {
        if (targetCount <= 0 || candidateRings.length === 0) return candidateRings;

        // Build edge index: map each bond to an integer index
        var compSet = {};
        for (var i = 0; i < atomIds.length; i++) compSet[atomIds[i]] = true;
        var edges = [], edgeIndex = {};
        for (var i = 0; i < atomIds.length; i++) {
            var nbrs = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < nbrs.length; j++) {
                if (compSet[nbrs[j]] && nbrs[j] > atomIds[i]) {
                    var key = atomIds[i] + ',' + nbrs[j];
                    edgeIndex[key] = edges.length;
                    edges.push(key);
                }
            }
        }
        var numEdges = edges.length;

        var selected = [], basis = [];
        for (var ci = 0; ci < candidateRings.length && selected.length < targetCount; ci++) {
            var ring = candidateRings[ci];
            var vec = new Uint8Array(numEdges);
            for (var i = 0; i < ring.length; i++) {
                var u = ring[i], v = ring[(i + 1) % ring.length];
                var key = Math.min(u, v) + ',' + Math.max(u, v);
                var idx = edgeIndex[key];
                if (idx !== undefined) vec[idx] = 1;
            }

            // Reduce against existing basis
            var reduced = vec.slice();
            for (var bi = 0; bi < basis.length; bi++) {
                if (reduced[basis[bi].pivot]) {
                    for (var k = 0; k < numEdges; k++) reduced[k] ^= basis[bi].vec[k];
                }
            }

            var pivot = -1;
            for (var k = 0; k < numEdges; k++) { if (reduced[k]) { pivot = k; break; } }

            if (pivot >= 0) {
                basis.push({ vec: reduced, pivot: pivot });
                selected.push(ring);
            }
        }
        return selected;
    }

    /**
     * Greedy SSSR selection via GF(2) Gaussian elimination over edge vectors.
     */
    function selectSSSR(candidates, n, adj, targetCount) {
        if (targetCount <= 0) return [];

        var edges = [], edgeIndex = {};
        for (var u = 0; u < n; u++) {
            for (var j = 0; j < adj[u].length; j++) {
                var v = adj[u][j];
                if (u < v) {
                    var key = u + ',' + v;
                    if (edgeIndex[key] === undefined) { edgeIndex[key] = edges.length; edges.push(key); }
                }
            }
        }

        var numEdges = edges.length;
        var selected = [], basis = [];

        for (var ci = 0; ci < candidates.length && selected.length < targetCount; ci++) {
            var ring = candidates[ci].ring;
            var vec = new Uint8Array(numEdges);
            for (var i = 0; i < ring.length; i++) {
                var u = ring[i], v = ring[(i + 1) % ring.length];
                var key = Math.min(u, v) + ',' + Math.max(u, v);
                var idx = edgeIndex[key];
                if (idx !== undefined) vec[idx] = 1;
            }

            var reduced = vec.slice();
            for (var bi = 0; bi < basis.length; bi++) {
                if (reduced[basis[bi].pivot]) {
                    for (var k = 0; k < numEdges; k++) reduced[k] ^= basis[bi].vec[k];
                }
            }

            var pivot = -1;
            for (var k = 0; k < numEdges; k++) { if (reduced[k]) { pivot = k; break; } }

            if (pivot >= 0) {
                basis.push({ vec: reduced, pivot: pivot });
                selected.push(candidates[ci]);
            }
        }
        return selected;
    }

    // =====================================================================
    // Ring Systems (fused ring clusters via union-find)
    // =====================================================================

    function buildRingSystems(rings) {
        if (rings.length === 0) return [];
        var n = rings.length;
        var par = [];
        for (var i = 0; i < n; i++) par[i] = i;
        function find(x) { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; }
        function unite(a, b) { a = find(a); b = find(b); if (a !== b) par[b] = a; }

        for (var i = 0; i < n; i++) {
            for (var j = i + 1; j < n; j++) {
                // Fused: share >= 2 atoms (edge); also catch spiro via >= 1
                if (sharedAtomCount(rings[i], rings[j]) >= 1) unite(i, j);
            }
        }

        var systems = {};
        for (var i = 0; i < n; i++) {
            var root = find(i);
            if (!systems[root]) systems[root] = [];
            systems[root].push(i);
        }
        var result = [];
        for (var key in systems) result.push(systems[key]);
        return result;
    }

    function sharedAtomCount(ring1, ring2) {
        var set = {};
        for (var i = 0; i < ring1.length; i++) set[ring1[i]] = true;
        var count = 0;
        for (var i = 0; i < ring2.length; i++) { if (set[ring2[i]]) count++; }
        return count;
    }

    function sharedAtoms(ring1, ring2) {
        var set = {};
        for (var i = 0; i < ring1.length; i++) set[ring1[i]] = true;
        var shared = [];
        for (var i = 0; i < ring2.length; i++) { if (set[ring2[i]]) shared.push(ring2[i]); }
        return shared;
    }

    // =====================================================================
    // Template matching for known ring systems
    // =====================================================================

    /**
     * Collect all unique atom IDs in a ring system (union of all rings).
     */
    function collectSystemAtoms(rings, systemRingIndices) {
        var seen = {};
        var result = [];
        for (var ri = 0; ri < systemRingIndices.length; ri++) {
            var ring = rings[systemRingIndices[ri]];
            for (var ai = 0; ai < ring.length; ai++) {
                if (!seen[ring[ai]]) {
                    seen[ring[ai]] = true;
                    result.push(ring[ai]);
                }
            }
        }
        return result;
    }

    /**
     * Signature-based ring system template matching.
     * For each ring system, computes a signature (sorted ring sizes + atom count)
     * and looks it up in a table of known scaffolds.
     * If matched, applies template coordinates to the ring system atoms.
     * Returns a map of placed atom IDs.
     */
    function matchRingSystemTemplates(mol, rings, ringSystems) {
        var placed = {};

        // Fast template matching: RASCAL pre-screen + greedy probe (no VF2++ backtracking)
        if (typeof Templates === 'undefined') return placed;

        // Signature -> template name lookup (ring count + sorted sizes only;
        // VF2++ handles exact atom mapping so atom count need not be exact)
        var TEMPLATE_LOOKUP = {
            // v1.5.2: greatly expanded coverage. Agent C (Templates.chromone +
            // Templates.coumarin + chromoneRingMap deterministic mapper +
            // placeFlavonoidArylRing) handles flavonoid scaffolds; Agent D
            // (37 new Templates.* entries spanning indole / quinoline / purine
            // / 5-mem heterocycle / polycyclic / β-lactam / tetracycline /
            // alkaloid / benzodiazepine / phenothiazine families) covers
            // the rest. Order matters: chromone first in 2:6,6 so the
            // chromoneRingMap special-case in matchRingSystemTemplates
            // claims any chroman-shaped molecule before the generic
            // bicyclic templates probe it.
            '1:3':       ['cyclopropane'],
            '1:4':       ['cyclobutane', 'betalactam'],
            '1:5':       ['cyclopentane', 'pyrrolidine', 'furanose', 'imidazole', 'pyrazole', 'thiophene', 'furan', 'oxazole', 'isoxazole', 'thiazole', 'pyrrole', 'triazole123', 'triazole124', 'tetrazole'],
            '1:6':       ['benzene', 'cyclohexane', 'pyranose', 'piperidine', 'pyridine', 'pyrimidine', 'piperazine', 'morpholine', 'tetrahydropyran'],
            '1:7':       ['cycloheptane'],
            // 2:4,5 (penam, carbapenem) and 2:4,6 (cepham) intentionally
            // not in TEMPLATE_LOOKUP — the generic algorithm produces a
            // cleaner layout for penicillin G's acyl side chain than the
            // template stamper does. Templates remain available via
            // Templates.apply(mol, 'penam', cx, cy) for explicit use.
            '2:5,5':     ['pyrrolopyrimidine', 'pyrrolizidine'],
            '2:5,6':     ['purine', 'indole', 'benzimidazole', 'benzofuran', 'benzothiazole', 'benzothiophene', 'xanthine', 'indolizidine'],
            '2:5,7':     ['tropane'],
            '2:6,6':     ['chromone', 'coumarin', 'naphthalene', 'quinoline', 'isoquinoline', 'quinazoline', 'quinoxaline', 'cinnoline', 'phthalazine', 'pteridine', 'chroman', 'flavanone', 'flavone', 'quinolizidine'],
            '2:6,7':     ['benzodiazepine'],
            '3:5,6,6':   ['carbazole', 'beta_carboline', 'fluorene'],
            '3:6,6,6':   ['phenanthrene', 'phenothiazine', 'acridine', 'anthracene'],
            '3:6,6,7':   ['dibenzazepine'],
            '4:5,6,6,6': ['steroid'],
            '4:6,6,6,6': ['tetracycline', 'pyrene', 'aporphine']
            // NOTE: 5-ring fused systems (morphinan) intentionally not
            // templated. Greedy MCS matching is unreliable for highly fused
            // bridged systems and produces stretched bonds; the bond-length
            // validator rejects bad mappings, but the safer fallback is the
            // fused-ring assembler (layoutRingSystem). When/if a more robust
            // backtracking matcher is added we can re-enable morphinan here.
        };

        for (var si = 0; si < ringSystems.length; si++) {
            var sysRingIdx = ringSystems[si];
            if (sysRingIdx.length < 2) continue;

            var sizes = [];
            for (var ri = 0; ri < sysRingIdx.length; ri++) {
                sizes.push(rings[sysRingIdx[ri]].length);
            }
            sizes.sort(function(a, b) { return a - b; });

            var sysAtoms = collectSystemAtoms(rings, sysRingIdx);
            var sig = sysRingIdx.length + ':' + sizes.join(',');

            var candidates = TEMPLATE_LOOKUP[sig];
            if (!candidates) continue;

            // Special case: 2:6,6 fused bicyclic with composition (9C, 1O) is
            // the chroman / chromone / coumarin / flavanone scaffold. The
            // greedy MCS probe fails to map the heteroatom O reliably (it
            // processes degree-3 fusion carbons first and the single O ends
            // up unmapped), so we compute the canonical mapping directly by
            // graph traversal from the unique O atom. This guarantees a
            // full 10-atom placement.
            if (sig === '2:6,6') {
                var chromMapping = chromoneRingMap(mol, sysAtoms, sysRingIdx, rings);
                if (chromMapping) {
                    var chromTmpl = (typeof Templates !== 'undefined' &&
                                     typeof Templates.chromone === 'function')
                                    ? Templates.chromone() : null;
                    if (chromTmpl &&
                        validateTemplateBondLengths(mol, chromMapping, chromTmpl, sysAtoms)) {
                        applyTemplateCoordsVF2(mol, chromMapping, chromTmpl, placed);
                        // Place the C2-aryl ring (flavonoid B-ring) right next
                        // to C2 with the canonical horizontal-fusion-style
                        // offset, so the inter-ring bond comes out at the
                        // correct length. Without this, layoutRingSystem
                        // places the aryl at origin (overlapping the
                        // chromanone) and resolveCollisions can leave the
                        // C2-aryl bond stretched to ~1.8 BL.
                        placeFlavonoidArylRing(mol, chromMapping, rings, ringSystems, placed);
                        continue;
                    }
                }
            }

            // Pick best template by element composition
            var tmpl = matchTemplateByElements(mol, sysAtoms, candidates);
            if (!tmpl) continue;

            // RASCAL screen + MCS greedy probe (O(N²), no backtracking)
            var mapping = greedyMatchTemplate(mol, sysAtoms, tmpl);
            if (!mapping) continue;

            // Validate the mapping by checking that bonds in the molecule
            // map to nearby template atom positions. If applying the template
            // would produce wildly stretched bonds (>1.6 BL) or compressed
            // bonds (<0.4 BL), skip the template — it doesn't fit topology
            // even if element counts match.
            if (!validateTemplateBondLengths(mol, mapping, tmpl, sysAtoms)) {
                continue;
            }

            // Apply template coordinates using the correct VF2++ mapping
            applyTemplateCoordsVF2(mol, mapping, tmpl, placed);
        }

        return placed;
    }

    /**
     * Fast template matching: RASCAL pre-screen + MCS greedy probe.
     * O(N) + O(N²) — no VF2++ backtracking, no NLF overhead, no hangs.
     * Returns {tmplIdx: molAtomId} or null.
     */
    function greedyMatchTemplate(mol, sysAtomIds, tmpl) {
        if (typeof SMSDMCS === 'undefined' || !SMSDMCS._greedyProbe) return null;
        if (typeof SMSDGraph === 'undefined') return null;

        try {
            // -- Stage 1: RASCAL element frequency screen (O(N)) --
            var sysElem = {}, tmplElem = {};
            for (var i = 0; i < sysAtomIds.length; i++) {
                var sym = (mol.getAtom(sysAtomIds[i]).symbol || 'C');
                sysElem[sym] = (sysElem[sym] || 0) + 1;
            }
            for (var i = 0; i < tmpl.atoms.length; i++) {
                var sym = (tmpl.atoms[i].symbol || 'C');
                tmplElem[sym] = (tmplElem[sym] || 0) + 1;
            }
            var matchable = 0;
            for (var key in tmplElem) {
                if (tmplElem.hasOwnProperty(key)) {
                    matchable += Math.min(tmplElem[key], sysElem[key] || 0);
                }
            }
            var smaller = Math.min(tmpl.atoms.length, sysAtomIds.length);
            if (matchable < smaller * 0.7) return null; // reject: element mismatch

            // -- Build Molecule objects for SMSDGraph --
            var tmplMol = new Molecule();
            var tmplAtomIds = [];
            for (var i = 0; i < tmpl.atoms.length; i++) {
                var ta = tmpl.atoms[i];
                // addAtom returns the Atom object — store its .id (numeric).
                tmplAtomIds.push(tmplMol.addAtom(ta.symbol || 'C', ta.x, ta.y).id);
            }
            for (var i = 0; i < tmpl.bonds.length; i++) {
                var b = tmpl.bonds[i];
                // tmpl.bonds entries use {a1, a2, type} or [a1, a2, type] form;
                // accept both (the keyed form is what every Templates.* function
                // produces in this file — the array form was an old shape).
                var ba1 = (b.a1 !== undefined) ? b.a1 : b[0];
                var ba2 = (b.a2 !== undefined) ? b.a2 : b[1];
                var btype = (b.type !== undefined) ? b.type : (b[2] || 1);
                tmplMol.addBond(tmplAtomIds[ba1], tmplAtomIds[ba2], btype);
            }

            var subMol = new Molecule();
            var sysIdMap = {}, subToMol = {};
            for (var i = 0; i < sysAtomIds.length; i++) {
                var ma = mol.getAtom(sysAtomIds[i]);
                var sid = subMol.addAtom(ma.symbol || 'C', 0, 0).id;
                sysIdMap[sysAtomIds[i]] = sid;
                subToMol[sid] = sysAtomIds[i];
            }
            for (var i = 0; i < sysAtomIds.length; i++) {
                var bonds = mol.getBondsOfAtom(sysAtomIds[i]);
                for (var bi = 0; bi < bonds.length; bi++) {
                    var bond = bonds[bi];
                    var other = bond.atom1 === sysAtomIds[i] ? bond.atom2 : bond.atom1;
                    if (sysIdMap[other] !== undefined && other > sysAtomIds[i]) {
                        subMol.addBond(sysIdMap[sysAtomIds[i]], sysIdMap[other], bond.type || 1);
                    }
                }
            }

            // -- Stage 2: MCS greedy probe (O(N²), no backtracking) --
            var gTmpl = new SMSDGraph.SMSDGraph(tmplMol);
            var gSub = new SMSDGraph.SMSDGraph(subMol);
            var opts = new SMSDGraph.ChemOptions();
            opts.matchBondOrder = 'any';
            opts.ringMatchesRingOnly = false;
            opts.matchFormalCharge = false;

            // Use smaller graph as g1 (query), larger as g2 (target)
            var g1, g2, tmplIsQuery;
            if (tmplMol.atoms.length <= subMol.atoms.length) {
                g1 = gTmpl; g2 = gSub; tmplIsQuery = true;
            } else {
                g1 = gSub; g2 = gTmpl; tmplIsQuery = false;
            }

            var greedyMap = SMSDMCS._greedyProbe(g1, g2, opts);
            if (!greedyMap) return null;

            // Count matched atoms
            var mapSize = 0;
            for (var k in greedyMap) { if (greedyMap.hasOwnProperty(k)) mapSize++; }
            // Lower the acceptance threshold to 50%: complex 5-ring fused
            // topologies (morphinan) only get partial greedy matches, but the
            // post-validation step rejects mappings that produce stretched
            // or compressed bonds, so a partial match is safe to attempt.
            if (mapSize < smaller * 0.5) return null; // too few matches

            // -- Stage 3: Verify bond mapping (O(bonds)) --
            var bondHits = 0, bondTotal = 0;
            for (var qi in greedyMap) {
                if (!greedyMap.hasOwnProperty(qi)) continue;
                qi = parseInt(qi, 10);
                var ti = greedyMap[qi];
                var nbs = g1.neighbors[qi];
                for (var ni = 0; ni < nbs.length; ni++) {
                    var qk = nbs[ni];
                    if (qk <= qi) continue; // count each bond once
                    if (!greedyMap.hasOwnProperty(qk)) continue;
                    bondTotal++;
                    var tk = greedyMap[qk];
                    if (g2.hasBond(ti, tk)) bondHits++;
                }
            }
            if (bondTotal > 0 && bondHits < bondTotal * 0.7) return null; // bad mapping

            // -- Convert to {tmplIdx: molAtomId} --
            var mapping = {};
            for (var qi in greedyMap) {
                if (!greedyMap.hasOwnProperty(qi)) continue;
                var ti = greedyMap[qi];
                if (tmplIsQuery) {
                    // qi=template idx, ti=subMol idx
                    var subAtomId = gSub.idxToId[ti];
                    mapping[qi] = subToMol[subAtomId];
                } else {
                    // qi=subMol idx, ti=template idx
                    var subAtomId = gSub.idxToId[parseInt(qi, 10)];
                    mapping[ti] = subToMol[subAtomId];
                }
            }
            return mapping;
        } catch (e) {
            return null;
        }
    }

    /**
     * Verify that a candidate template mapping doesn't produce absurd bond
     * lengths when its coordinates are applied. We compute, for every bond
     * in the ring system that has both endpoints mapped, the implied
     * distance using the template coordinates. If any bond is stretched to
     * more than 1.6 BL or compressed below 0.4 BL the mapping is rejected
     * (template topology doesn't match).
     */
    function validateTemplateBondLengths(mol, mapping, tmpl, sysAtomIds) {
        // mapping = { tmplIdx: molAtomId }
        // Build reverse: molAtomId -> tmplIdx
        var molToTmpl = {};
        var mappedCount = 0;
        for (var ti in mapping) {
            if (mapping.hasOwnProperty(ti)) {
                molToTmpl[mapping[ti]] = parseInt(ti, 10);
                mappedCount++;
            }
        }

        // Reject if fewer than 80% of system atoms are mapped — partial
        // matches against an unrelated template (e.g. indole partially
        // matching atropine's tropane) produce garbage layouts even when
        // bond lengths happen to validate.
        if (mappedCount < sysAtomIds.length * 0.8) return false;

        var sysSet = {};
        for (var i = 0; i < sysAtomIds.length; i++) sysSet[sysAtomIds[i]] = true;

        var maxLen = 0, minLen = Infinity;
        var bondsChecked = 0;

        for (var i = 0; i < sysAtomIds.length; i++) {
            var aId = sysAtomIds[i];
            if (molToTmpl[aId] === undefined) continue;
            var nbrs = mol.getNeighbors(aId);
            for (var j = 0; j < nbrs.length; j++) {
                var nId = nbrs[j];
                if (nId <= aId) continue;
                if (!sysSet[nId]) continue;
                if (molToTmpl[nId] === undefined) continue;

                var ta = tmpl.atoms[molToTmpl[aId]];
                var tb = tmpl.atoms[molToTmpl[nId]];
                var dx = ta.x - tb.x, dy = ta.y - tb.y;
                var len = Math.sqrt(dx * dx + dy * dy);
                if (len > maxLen) maxLen = len;
                if (len < minLen) minLen = len;
                bondsChecked++;
            }
        }

        if (bondsChecked === 0) return false;
        if (maxLen > BOND_LENGTH * 1.6) return false;
        if (minLen < BOND_LENGTH * 0.4) return false;
        return true;
    }

    /**
     * Deterministic atom mapping for the chromone / chroman / coumarin /
     * flavanone scaffold (fused 6,6 with 9C + 1O = 10 ring atoms).
     *
     * The greedy MCS probe in SMSDMCS prioritises high-degree atoms (the two
     * fusion carbons) and processes the unique ring oxygen last; with no
     * backtracking it routinely fails to place the O, leaving a partial
     * 8/10 mapping that the validator either rejects (-> falls through to
     * the fused-ring assembler, producing the stretched depictions seen in
     * Eriodictyol / Naringenin / Quercetin) or applies incompletely. Since
     * the topology is fixed for any chroman-shaped bicyclic with one ring
     * oxygen, we compute the mapping by direct graph traversal:
     *
     *   1. The unique ring O is anchored at template idx 9.
     *   2. Of the O's two ring neighbours, the one with degree 3 inside the
     *      ring system is C8a (the fusion vertex; template idx 1); the other
     *      is C2 (template idx 8).
     *   3. The other fusion vertex (template idx 2) is C4a — the unique
     *      degree-3 ring-system atom that is *not* C8a.
     *   4. Walk benzene (C8a -> C8 -> C7 -> C6 -> C5 -> C4a) for indices
     *      1, 0, 5, 4, 3, 2.
     *   5. Walk the pyran ring (C2 -> C3 -> C4 -> C4a) for indices 8, 7, 6, 2.
     *
     * Returns {tmplIdx: molAtomId} on success, or null if the ring system is
     * not a chroman (zero or multiple ring oxygens, or topology mismatch).
     */
    function chromoneRingMap(mol, sysAtoms, sysRingIdx, rings) {
        // Verify composition: exactly 9 C + 1 O.
        var oxygenAtomId = -1;
        var oxygenCount = 0;
        var carbonCount = 0;
        for (var i = 0; i < sysAtoms.length; i++) {
            var sym = mol.getAtom(sysAtoms[i]).symbol || 'C';
            if (sym === 'O') { oxygenAtomId = sysAtoms[i]; oxygenCount++; }
            else if (sym === 'C') { carbonCount++; }
            else { return null; } // unsupported element in ring system
        }
        if (oxygenCount !== 1 || carbonCount !== 9) return null;

        // Build set of ring-system atoms for O(1) membership checks.
        var sysSet = {};
        for (var i = 0; i < sysAtoms.length; i++) sysSet[sysAtoms[i]] = true;

        // Per-atom degree restricted to ring-system bonds (i.e. how many of
        // an atom's neighbours are also ring-system members). Fusion vertices
        // have degree 3; all other ring atoms have degree 2.
        var ringDegree = {};
        for (var i = 0; i < sysAtoms.length; i++) {
            var aid = sysAtoms[i];
            var nbrs = mol.getNeighbors(aid);
            var d = 0;
            for (var j = 0; j < nbrs.length; j++) {
                if (sysSet[nbrs[j]]) d++;
            }
            ringDegree[aid] = d;
        }

        // The two ring-system neighbours of the oxygen.
        var oxNbrs = mol.getNeighbors(oxygenAtomId).filter(function(n) {
            return sysSet[n];
        });
        if (oxNbrs.length !== 2) return null;

        // Identify C8a (fusion, degree 3) and C2 (degree 2) among O's neighbours.
        var c8a, c2;
        if (ringDegree[oxNbrs[0]] === 3 && ringDegree[oxNbrs[1]] === 2) {
            c8a = oxNbrs[0]; c2 = oxNbrs[1];
        } else if (ringDegree[oxNbrs[1]] === 3 && ringDegree[oxNbrs[0]] === 2) {
            c8a = oxNbrs[1]; c2 = oxNbrs[0];
        } else {
            // O is between two fusion vertices, or both neighbours are
            // degree-2 — not a chroman scaffold (could be 1,4-dioxin etc.).
            return null;
        }

        // Find C4a — the OTHER fusion vertex (degree 3 in ring system, not c8a).
        var c4a = -1;
        for (var i = 0; i < sysAtoms.length; i++) {
            var aid = sysAtoms[i];
            if (aid === c8a) continue;
            if (ringDegree[aid] === 3) { c4a = aid; break; }
        }
        if (c4a < 0) return null;

        // Walk the benzene ring: C8a -> C8 -> C7 -> C6 -> C5 -> C4a (5 hops).
        // C8a's benzene-ring neighbour is the one that is neither O nor C4a.
        var benz = [c8a];
        var prev = oxygenAtomId; // pretend O is "previous" so we don't walk back to it
        var cur = c8a;
        for (var step = 0; step < 5; step++) {
            var nbrs = mol.getNeighbors(cur).filter(function(n) {
                return sysSet[n] && n !== prev;
            });
            // From C8a (degree 3 in ring) the next atom must be the unique
            // benzene neighbour that is also not C4a; from intermediate atoms
            // (degree 2) the next is the unique non-prev neighbour. To stay on
            // the benzene side, exclude C4a until the final hop.
            var next = -1;
            if (cur === c8a) {
                // exclude C4a here so we don't shortcut across the fusion edge
                for (var ni = 0; ni < nbrs.length; ni++) {
                    if (nbrs[ni] !== c4a) { next = nbrs[ni]; break; }
                }
            } else {
                next = nbrs.length > 0 ? nbrs[0] : -1;
            }
            if (next < 0) return null;
            benz.push(next);
            prev = cur;
            cur = next;
        }
        if (cur !== c4a) return null; // walked off-scaffold
        // benz is now [C8a, C8, C7, C6, C5, C4a].

        // Walk the pyran ring: C2 -> C3 -> C4 -> C4a (3 hops).
        var pyran = [c2];
        prev = oxygenAtomId;
        cur = c2;
        for (var step = 0; step < 3; step++) {
            var nbrs = mol.getNeighbors(cur).filter(function(n) {
                return sysSet[n] && n !== prev;
            });
            // From C2 (degree 2) the unique next is C3; from C3 unique next
            // is C4; from C4 unique next is C4a.
            var next = nbrs.length > 0 ? nbrs[0] : -1;
            if (next < 0) return null;
            pyran.push(next);
            prev = cur;
            cur = next;
        }
        if (cur !== c4a) return null;
        // pyran is now [C2, C3, C4, C4a].

        // Build the template-index -> mol-id mapping.
        // Template indices: 0=C8, 1=C8a, 2=C4a, 3=C5, 4=C6, 5=C7, 6=C4, 7=C3, 8=C2, 9=O1.
        var mapping = {};
        mapping[1] = benz[0];   // C8a
        mapping[0] = benz[1];   // C8
        mapping[5] = benz[2];   // C7
        mapping[4] = benz[3];   // C6
        mapping[3] = benz[4];   // C5
        mapping[2] = benz[5];   // C4a
        mapping[8] = pyran[0];  // C2
        mapping[7] = pyran[1];  // C3
        mapping[6] = pyran[2];  // C4
        mapping[9] = oxygenAtomId; // O1
        return mapping;
    }

    /**
     * After applying the chromone template, position the C2-aryl ring
     * (the flavonoid B-ring) so the C2-aryl bond comes out at the correct
     * length and the B-ring sits to the right of the chromanone, matching
     * the canonical depiction of flavanones / flavones / isoflavones.
     *
     * If C2 has a non-ring-system neighbour that itself belongs to a
     * standalone (single-ring, unfused) 6-membered ring, that ring is the
     * B-ring. We anchor its centre `BOND_LENGTH + radius` away from C2 along
     * the direction continuing from the chromanone-C2 axis (i.e. the
     * horizontal fusion axis used by the template), then fix the polygon
     * so the C2-bonded vertex sits exactly BOND_LENGTH from C2.
     */
    function placeFlavonoidArylRing(mol, chromMapping, rings, ringSystems, placed) {
        // chromMapping: 8 = C2, 7 = C3 (template indices). Either may bear
        // the aryl substituent (flavanone/flavone -> C2; isoflavone -> C3).
        var chromAtoms = {};
        for (var ti in chromMapping) {
            if (chromMapping.hasOwnProperty(ti)) chromAtoms[chromMapping[ti]] = true;
        }

        // Try C2 first (template idx 8), fall back to C3 (template idx 7).
        var anchorIds = [chromMapping[8], chromMapping[7]];
        var anchorId = -1;
        var arylSeed = -1;
        for (var ai = 0; ai < anchorIds.length && arylSeed < 0; ai++) {
            var aId = anchorIds[ai];
            if (aId === undefined) continue;
            var nbrs = mol.getNeighbors(aId);
            for (var ni = 0; ni < nbrs.length; ni++) {
                if (chromAtoms[nbrs[ni]]) continue;
                if (placed[nbrs[ni]]) continue;
                // Only consider neighbours that are themselves ring atoms in
                // a 6-membered ring (the aryl B-ring).
                var symbol = mol.getAtom(nbrs[ni]).symbol || 'C';
                if (symbol !== 'C') continue;
                arylSeed = nbrs[ni];
                anchorId = aId;
                break;
            }
        }
        if (arylSeed < 0 || anchorId < 0) return;
        var c2Atom = mol.getAtom(anchorId);
        if (!c2Atom) return;

        // Find the smallest standalone (single-ring) 6-membered ring that
        // contains arylSeed. If none, this is just a chain substituent and
        // no template placement is needed.
        var arylRing = null;
        for (var ri = 0; ri < rings.length; ri++) {
            var r = rings[ri];
            if (r.length !== 6) continue;
            if (r.indexOf(arylSeed) < 0) continue;
            // Standalone ring: it must not share any atom with the chromanone
            // (otherwise it's a fused part of a polycyclic system).
            var fused = false;
            for (var ai = 0; ai < r.length; ai++) {
                if (chromAtoms[r[ai]]) { fused = true; break; }
            }
            if (fused) continue;
            // Also reject if the ring shares atoms with another ring (i.e.
            // belongs to a multi-ring system other than the chromanone) —
            // we only handle the simple flavonoid B-ring case here.
            var otherFused = false;
            for (var rj = 0; rj < rings.length; rj++) {
                if (rj === ri) continue;
                for (var aj = 0; aj < r.length; aj++) {
                    if (rings[rj].indexOf(r[aj]) >= 0) { otherFused = true; break; }
                }
                if (otherFused) break;
            }
            if (otherFused) continue;
            arylRing = r;
            break;
        }
        if (!arylRing) return;

        // Direction for the aryl ring: project away from the chromanone
        // centroid (so the B-ring sits on the "outside" of the bicyclic).
        // Use the average of the chromanone atoms' positions as the centroid.
        var ccx = 0, ccy = 0, ccount = 0;
        for (var aId in chromAtoms) {
            if (chromAtoms.hasOwnProperty(aId)) {
                var ca = mol.getAtom(parseInt(aId, 10));
                if (ca) { ccx += ca.x; ccy += ca.y; ccount++; }
            }
        }
        var dirX = 1, dirY = 0;
        if (ccount > 0) {
            ccx /= ccount; ccy /= ccount;
            var dx = c2Atom.x - ccx;
            var dy = c2Atom.y - ccy;
            var len = Math.sqrt(dx * dx + dy * dy);
            if (len > EPSILON) { dirX = dx / len; dirY = dy / len; }
        }

        // Place the aryl seed atom one BOND_LENGTH from C2 along dir.
        var seedX = c2Atom.x + BOND_LENGTH * dirX;
        var seedY = c2Atom.y + BOND_LENGTH * dirY;

        // Compute ring centre: one apothem beyond seedX/Y along dir.
        var radius = BOND_LENGTH / (2 * Math.sin(Math.PI / 6)); // = BOND_LENGTH for 6-ring
        var apothem = radius * Math.cos(Math.PI / 6);
        var ringCx = seedX + apothem * dirX;
        var ringCy = seedY + apothem * dirY;

        // Compute the polygon offset so the seed atom sits exactly at angle
        // pointing back toward C2 (i.e. at angle (-dirX, -dirY) from centre).
        var seedAngle = Math.atan2(-dirY, -dirX);

        // Reorder ring so seedAtom is at index 0, then walk around.
        var seedIdx = arylRing.indexOf(arylSeed);
        var step = TWO_PI / 6;
        for (var k = 0; k < 6; k++) {
            var ringAtomId = arylRing[(seedIdx + k) % 6];
            if (placed[ringAtomId]) continue;
            var atom = mol.getAtom(ringAtomId);
            // Walk around the ring keeping a consistent rotation. The seed
            // sits at seedAngle; subsequent atoms advance by +step or -step.
            // Use +step (counter-clockwise in screen coords; canvas y is
            // down, so this matches the polygon helper's convention).
            var angle = seedAngle + k * step;
            atom.x = ringCx + radius * Math.cos(angle);
            atom.y = ringCy + radius * Math.sin(angle);
            placed[ringAtomId] = true;
        }
    }

    /**
     * Apply template coordinates using VF2++-derived mapping.
     * mapping = {tmplIdx: molAtomId}
     */
    function applyTemplateCoordsVF2(mol, mapping, tmpl, placed) {
        // Compute template centroid
        var tcx = 0, tcy = 0, count = 0;
        for (var ti in mapping) {
            if (mapping.hasOwnProperty(ti)) {
                tcx += tmpl.atoms[ti].x;
                tcy += tmpl.atoms[ti].y;
                count++;
            }
        }
        if (count === 0) return;
        tcx /= count;
        tcy /= count;

        // Apply centered template coordinates to molecule atoms
        for (var ti in mapping) {
            if (mapping.hasOwnProperty(ti)) {
                var atom = mol.getAtom(mapping[ti]);
                if (atom) {
                    atom.x = tmpl.atoms[ti].x - tcx;
                    atom.y = tmpl.atoms[ti].y - tcy;
                    placed[mapping[ti]] = true;
                }
            }
        }
    }

    /**
     * Given candidate template names, pick the best match based on element counts.
     */
    function matchTemplateByElements(mol, sysAtomIds, candidates) {
        // Count elements in the ring system
        var elemCounts = {};
        for (var i = 0; i < sysAtomIds.length; i++) {
            var atom = mol.getAtom(sysAtomIds[i]);
            var sym = atom.symbol || 'C';
            elemCounts[sym] = (elemCounts[sym] || 0) + 1;
        }

        var bestTemplate = null;
        var bestScore = -1;

        for (var ci = 0; ci < candidates.length; ci++) {
            var tmplName = candidates[ci];
            if (typeof Templates === 'undefined' || typeof Templates[tmplName] !== 'function') continue;
            var tmpl = Templates[tmplName]();
            if (!tmpl) continue;
            // Allow template to be same size or slightly smaller (bridged rings
            // may have different atom counts depending on SSSR perception)
            if (tmpl.atoms.length > sysAtomIds.length + 2) continue;

            // Count elements in template
            var tmplElem = {};
            for (var ai = 0; ai < tmpl.atoms.length; ai++) {
                var s = tmpl.atoms[ai].symbol || 'C';
                tmplElem[s] = (tmplElem[s] || 0) + 1;
            }

            // Score: count matching element counts
            var score = 0;
            for (var key in tmplElem) {
                if (tmplElem.hasOwnProperty(key)) {
                    if (tmplElem[key] === (elemCounts[key] || 0)) {
                        score += tmplElem[key];
                    } else {
                        // Partial credit for close matches
                        var diff = Math.abs(tmplElem[key] - (elemCounts[key] || 0));
                        score -= diff;
                    }
                }
            }
            // Penalise extra elements in mol not in template
            for (var key in elemCounts) {
                if (elemCounts.hasOwnProperty(key) && !tmplElem[key]) {
                    score -= elemCounts[key];
                }
            }

            if (score > bestScore) {
                bestScore = score;
                bestTemplate = tmpl;
            }
        }

        return bestTemplate;
    }

    /**
     * Apply template coordinates to ring system atoms.
     * Aligns the template centroid to the origin.
     */
    function applyTemplateCoords(mol, sysAtomIds, tmpl, placed) {
        if (tmpl.atoms.length !== sysAtomIds.length) return;

        // Compute template centroid
        var tcx = 0, tcy = 0;
        for (var i = 0; i < tmpl.atoms.length; i++) {
            tcx += tmpl.atoms[i].x;
            tcy += tmpl.atoms[i].y;
        }
        tcx /= tmpl.atoms.length;
        tcy /= tmpl.atoms.length;

        // Apply coordinates (centered at origin)
        for (var i = 0; i < sysAtomIds.length; i++) {
            var atom = mol.getAtom(sysAtomIds[i]);
            atom.x = tmpl.atoms[i].x - tcx;
            atom.y = tmpl.atoms[i].y - tcy;
            placed[sysAtomIds[i]] = true;
        }
    }

    // =====================================================================
    // Ring System Layout  (Helson-style: polygon + edge-reflection fusion)
    // =====================================================================

    /**
     * Layout a system of fused / spiro rings.
     *
     * Ring placement strategy:
     *   - Place the largest ring first (better scaffold for steroids, etc.)
     *   - BFS fuse subsequent rings by reflecting across shared edges
     *   - Spiro rings placed by rotating away from existing bonds
     */
    function layoutRingSystem(mol, allRings, systemRingIndices, placed, ringAtomSet) {
        if (systemRingIndices.length === 0) return;

        // Sort: place largest ring first (gives best scaffold for polycyclics)
        systemRingIndices.sort(function (a, b) {
            return allRings[b].length - allRings[a].length;
        });

        var placedRings = {};
        var ringQueue   = [systemRingIndices[0]];
        placedRings[systemRingIndices[0]] = true;

        // Determine the optimal starting angle for the seed ring.
        // For fused bicyclic systems (e.g. naphthalene), rotate the seed ring
        // so the shared edge is horizontal, producing the standard horizontal layout.
        var firstRing = allRings[systemRingIndices[0]];
        var seedAngle = -Math.PI / 2; // default: top vertex

        if (systemRingIndices.length >= 2) {
            // Find the first fused neighbour (shares >= 2 atoms)
            for (var si = 1; si < systemRingIndices.length; si++) {
                var shared = sharedAtoms(firstRing, allRings[systemRingIndices[si]]);
                if (shared.length >= 2) {
                    // For the shared edge to be horizontal after layout, we need
                    // the two shared atoms to have the same y-coordinate.
                    // Find their indices in the seed ring.
                    var idx1 = firstRing.indexOf(shared[0]);
                    var idx2 = firstRing.indexOf(shared[1]);
                    var size = firstRing.length;
                    var step = TWO_PI / size;

                    // The default polygon places vertex i at angle: seedAngle + i*step
                    // (with even-size offset of step/2).
                    // The shared edge spans from vertex idx1 to vertex idx2.
                    // For that edge to be horizontal (same y), the midpoint angle
                    // of the edge must be 0 or PI (pointing left or right).
                    var offset = (size % 2 === 0) ? step / 2 : 0;
                    var defaultAngle1 = -Math.PI / 2 + offset + idx1 * step;
                    var defaultAngle2 = -Math.PI / 2 + offset + idx2 * step;
                    var edgeMidAngle = (defaultAngle1 + defaultAngle2) / 2;

                    // We want edgeMidAngle to be 0 (right) or PI (left).
                    // Rotation needed = target - current.
                    // Choose the target (0 or PI) that requires less rotation.
                    var rot0 = normalizeAngle(-edgeMidAngle);
                    var rotPI = normalizeAngle(Math.PI - edgeMidAngle);
                    if (rot0 > Math.PI) rot0 = TWO_PI - rot0;
                    if (rotPI > Math.PI) rotPI = TWO_PI - rotPI;
                    var targetMid = (rot0 <= rotPI) ? 0 : Math.PI;
                    seedAngle = -Math.PI / 2 + (targetMid - edgeMidAngle);
                    break;
                }
            }
        }

        placeRingAsPolygon(mol, firstRing, 0, 0, seedAngle, placed);

        // Greedy ring placement: at each step pick the unplaced ring with
        // the highest "fitness" against the already-placed system. The
        // ordering is:
        //   1. PREFER simple edge-fusion (exactly 2 placed atoms, contiguous)
        //      over bridged (3+ placed atoms shared) — bridged rings create
        //      crowded geometry, so let edge-fused rings establish the
        //      scaffold first.
        //   2. Among edge-fused rings, prefer larger.
        //   3. Bridged comes next (3+ shared placed atoms).
        //   4. Spiro (1 shared) last.
        var madeProgress = true;
        while (madeProgress) {
            madeProgress = false;

            var bestIdx = -1, bestScore = -1, bestShared = null;
            for (var si = 0; si < systemRingIndices.length; si++) {
                var nextIdx = systemRingIndices[si];
                if (placedRings[nextIdx]) continue;
                var ring = allRings[nextIdx];
                var sharedPlaced = [];
                for (var ai = 0; ai < ring.length; ai++) {
                    if (placed[ring[ai]]) sharedPlaced.push(ring[ai]);
                }
                if (sharedPlaced.length === 0) continue;

                // Determine whether shared atoms are contiguous (edge-fused)
                // by walking the ring.
                var contiguousLen = 0;
                var maxContig = 0;
                for (var ai = 0; ai < ring.length * 2; ai++) {
                    var idx = ai % ring.length;
                    if (placed[ring[idx]]) {
                        contiguousLen++;
                        if (contiguousLen > maxContig) maxContig = contiguousLen;
                    } else {
                        contiguousLen = 0;
                    }
                }
                if (maxContig > sharedPlaced.length) maxContig = sharedPlaced.length;

                var score;
                if (sharedPlaced.length === 2 && maxContig === 2) {
                    // Simple edge-fused (best): pick first
                    score = 1000 + ring.length;
                } else if (sharedPlaced.length === 1) {
                    // Spiro (worst non-bridged)
                    score = 100 + ring.length;
                } else if (maxContig >= 3) {
                    // Bridged (3+ contiguous): place after edge-fused
                    score = 500 + sharedPlaced.length * 10 + ring.length;
                } else {
                    // Mixed (e.g. 2 non-contiguous shared atoms — rare but
                    // can happen with overlapping bridged systems)
                    score = 300 + sharedPlaced.length * 10 + ring.length;
                }
                if (score > bestScore) {
                    bestScore = score;
                    bestIdx = nextIdx;
                    bestShared = sharedPlaced;
                }
            }

            if (bestIdx < 0) break;

            if (bestShared.length >= 2) {
                fuseRing(mol, allRings[bestIdx], bestShared, placed, ringAtomSet, allRings[bestIdx].length);
            } else {
                placeSpiroRing(mol, allRings[bestIdx], bestShared[0], placed);
            }
            placedRings[bestIdx] = true;
            madeProgress = true;
        }

        // Catch any still-unplaced rings (disconnected inside system)
        for (var si = 0; si < systemRingIndices.length; si++) {
            var ri = systemRingIndices[si];
            if (placedRings[ri]) continue;
            var ring = allRings[ri];

            // Try spiro attachment to any placed atom in this ring
            var attached = false;
            for (var ai = 0; ai < ring.length; ai++) {
                if (placed[ring[ai]]) {
                    placeSpiroRing(mol, ring, ring[ai], placed);
                    attached = true;
                    break;
                }
            }
            if (!attached) {
                placeRingAsPolygon(mol, ring, 0, 0, -Math.PI / 2, placed);
            }
            placedRings[ri] = true;
        }
    }

    // -----------------------------------------------------------------
    // Place a ring as a regular polygon
    // -----------------------------------------------------------------

    function placeRingAsPolygon(mol, ring, cx, cy, startAngle, placed) {
        var size = ring.length;
        var step = TWO_PI / size;
        var radius = BOND_LENGTH / (2 * Math.sin(Math.PI / size));

        // Macrocycle tightening: rings with 12+ atoms get a gradually
        // reduced radius so that they don't expand to absurdly large polygons.
        if (size >= 12) {
            var t = Math.min((size - 12) / 18, 1.0);
            radius *= 0.85 - 0.30 * t;
        }

        // For even-sized rings rotate so a flat edge sits at the bottom
        var offset = startAngle;
        if (size % 2 === 0) offset += step / 2;

        // Sugar moiety detection (Haworth convention): if this is a 5- or
        // 6-membered ring containing exactly one O, rotate the ring so that
        // the oxygen sits at the conventional "top-right" position (the
        // vertex at roughly the 1-o'clock position, angle ~-PI/6).
        if ((size === 5 || size === 6) && !hasAnyPlaced(ring, placed)) {
            var oIdx = -1, oCount = 0;
            for (var i = 0; i < size; i++) {
                var sym = mol.getAtom(ring[i]).symbol;
                if (sym === 'O' || sym === 'S') { oIdx = i; oCount++; }
            }
            if (oCount === 1 && oIdx >= 0) {
                // Target angle for the O atom: top-right (~-30 deg = -PI/6)
                var targetAngle = -Math.PI / 6;
                var currentAngle = offset + oIdx * step;
                var rotation = targetAngle - currentAngle;
                offset += rotation;
            }
        }

        for (var i = 0; i < size; i++) {
            if (!placed[ring[i]]) {
                var angle = offset + i * step;
                var atom  = mol.getAtom(ring[i]);
                atom.x = cx + radius * Math.cos(angle);
                atom.y = cy + radius * Math.sin(angle);
                placed[ring[i]] = true;
            }
        }
    }

    /**
     * Check if any atom in the ring is already placed.
     */
    function hasAnyPlaced(ring, placed) {
        for (var i = 0; i < ring.length; i++) {
            if (placed[ring[i]]) return true;
        }
        return false;
    }

    // -----------------------------------------------------------------
    // Fuse a new ring across a shared edge  (Helson edge-reflection)
    // -----------------------------------------------------------------

    function fuseRing(mol, ring, sharedAtomIds, placed, ringAtomSet, parentRingSize) {
        var size   = ring.length;
        var radius = BOND_LENGTH / (2 * Math.sin(Math.PI / size));

        // Macrocycle tightening for fused macrocyclic rings
        if (size >= 12) {
            var t = Math.min((size - 12) / 18, 1.0);
            radius *= 0.85 - 0.30 * t;
        }

        // Collect placed shared atoms indices in the ring
        var sharedIdx = [];
        for (var i = 0; i < ring.length; i++) {
            if (sharedAtomIds.indexOf(ring[i]) >= 0 && placed[ring[i]]) {
                sharedIdx.push(i);
            }
        }
        if (sharedIdx.length < 2) {
            placeRingAsPolygon(mol, ring, 0, 0, -Math.PI / 2, placed);
            return;
        }

        // -----------------------------------------------------------
        // Bridged bicyclic detection: if 3+ atoms are shared between
        // rings, we have a bridged system (e.g., norbornane, camphor).
        // The bridge atoms (those not on the shared edge endpoints)
        // must be placed above/below the ring plane.
        // -----------------------------------------------------------
        if (sharedIdx.length > 2) {
            fuseBridgedRing(mol, ring, sharedIdx, placed, ringAtomSet);
            return;
        }

        // Use the first two placed shared atoms to define the shared edge
        var a1 = mol.getAtom(ring[sharedIdx[0]]);
        var a2 = mol.getAtom(ring[sharedIdx[1]]);

        // Midpoint & normal of shared edge
        var mx = (a1.x + a2.x) / 2, my = (a1.y + a2.y) / 2;
        var edx = a2.x - a1.x, edy = a2.y - a1.y;
        var nx = -edy, ny = edx;
        var nLen = Math.sqrt(nx * nx + ny * ny);
        if (nLen > EPSILON) { nx /= nLen; ny /= nLen; }

        // Two candidate centres on either side of the edge
        var apothem = radius * Math.cos(Math.PI / size);
        var cx1 = mx + nx * apothem, cy1 = my + ny * apothem;
        var cx2 = mx - nx * apothem, cy2 = my - ny * apothem;

        // Score: prefer the centre farther from all already-placed atoms
        var score1 = 0, score2 = 0;
        for (var id in placed) {
            if (!placed[id]) continue;
            var pa = mol.getAtom(parseInt(id));
            if (!pa) continue;
            var d1 = distSq(cx1, cy1, pa.x, pa.y);
            var d2 = distSq(cx2, cy2, pa.x, pa.y);
            // Heavy penalty for very close clashes
            if (d1 < BOND_LENGTH * BOND_LENGTH * 0.6) score1 -= 200;
            if (d2 < BOND_LENGTH * BOND_LENGTH * 0.6) score2 -= 200;
            score1 += d1;
            score2 += d2;
        }

        var cx = score1 >= score2 ? cx1 : cx2;
        var cy = score1 >= score2 ? cy1 : cy2;

        // Compute the starting angle from the chosen centre to a1
        var angle1 = Math.atan2(a1.y - cy, a1.x - cx);
        var angle2 = Math.atan2(a2.y - cy, a2.x - cx);
        var step   = TWO_PI / size;

        var idx1 = sharedIdx[0], idx2 = sharedIdx[1];
        var angleDiff    = normalizeAngle(angle2 - angle1);
        var stepsForward = (idx2 - idx1 + size) % size;

        // Determine winding direction.
        // Dynamic tolerance based on ring size mismatch (e.g. 5+6 in caffeine)
        var newRingSize = size;
        var parentStep = 2 * Math.PI / parentRingSize;
        var childStep = 2 * Math.PI / newRingSize;
        var WINDING_TOL = Math.max(Math.abs(parentStep - childStep) * 1.5, 0.15);
        var expectedFwd = stepsForward * step;
        var expectedBwd = TWO_PI - expectedFwd;
        var clockwise;
        var diffFwd = Math.abs(angleDiff - expectedFwd);
        var diffBwd = Math.abs(angleDiff - expectedBwd);
        if (diffFwd < WINDING_TOL && diffFwd <= diffBwd) {
            clockwise = true;
        } else if (diffBwd < WINDING_TOL && diffBwd < diffFwd) {
            clockwise = false;
        } else {
            // Heuristic: choose direction that puts new atoms on the outer side
            clockwise = angleDiff > Math.PI ? false : true;
        }

        // Place each unplaced atom at its polygon vertex
        for (var i = 0; i < size; i++) {
            if (!placed[ring[i]]) {
                var steps = (i - idx1 + size) % size;
                var angle = clockwise ? angle1 + steps * step : angle1 - steps * step;
                var atom  = mol.getAtom(ring[i]);
                atom.x = cx + radius * Math.cos(angle);
                atom.y = cy + radius * Math.sin(angle);
                placed[ring[i]] = true;
            }
        }
    }

    // -----------------------------------------------------------------
    // Bridged bicyclic ring layout: rings sharing 3+ atoms
    // (e.g., norbornane, camphor). Bridge atoms are placed above/below
    // the plane defined by the two bridgehead atoms.
    // -----------------------------------------------------------------

    function fuseBridgedRing(mol, ring, sharedIdx, placed, ringAtomSet) {
        var size = ring.length;

        // Identify bridgehead atoms: the two shared atoms that are farthest
        // apart in the ring walk. These are the endpoints of the shared path.
        // Among the shared indices, find the pair with the maximum ring
        // distance (the two bridgeheads in a bicyclic system).
        var bh1Idx = sharedIdx[0], bh2Idx = sharedIdx[sharedIdx.length - 1];

        // Try all pairs of shared indices for maximal ring distance
        var maxDist = 0;
        for (var i = 0; i < sharedIdx.length; i++) {
            for (var j = i + 1; j < sharedIdx.length; j++) {
                var d = Math.min(
                    (sharedIdx[j] - sharedIdx[i] + size) % size,
                    (sharedIdx[i] - sharedIdx[j] + size) % size
                );
                if (d > maxDist) {
                    maxDist = d;
                    bh1Idx = sharedIdx[i];
                    bh2Idx = sharedIdx[j];
                }
            }
        }

        var bh1 = mol.getAtom(ring[bh1Idx]);
        var bh2 = mol.getAtom(ring[bh2Idx]);

        // Midpoint and normal of the bridgehead axis
        var mx = (bh1.x + bh2.x) / 2, my = (bh1.y + bh2.y) / 2;
        var edx = bh2.x - bh1.x, edy = bh2.y - bh1.y;
        var edLen = Math.sqrt(edx * edx + edy * edy);
        if (edLen < EPSILON) edLen = 1;
        // Normal (perpendicular) to the bridgehead axis
        var nx = -edy / edLen, ny = edx / edLen;

        // Separate unplaced atoms into two paths around the ring
        // Path A: bh1Idx -> bh2Idx (forward), Path B: bh2Idx -> bh1Idx (forward)
        var pathA = [], pathB = [];
        var stepsAtoB = (bh2Idx - bh1Idx + size) % size;
        for (var s = 1; s < stepsAtoB; s++) {
            var idx = (bh1Idx + s) % size;
            if (!placed[ring[idx]]) pathA.push(idx);
        }
        var stepsBtoA = (bh1Idx - bh2Idx + size) % size;
        for (var s = 1; s < stepsBtoA; s++) {
            var idx = (bh2Idx + s) % size;
            if (!placed[ring[idx]]) pathB.push(idx);
        }

        // Determine which side of the bridgehead axis is less crowded
        var scorePos = 0, scoreNeg = 0;
        for (var id in placed) {
            if (!placed[id]) continue;
            var pa = mol.getAtom(parseInt(id));
            if (!pa) continue;
            var dot = (pa.x - mx) * nx + (pa.y - my) * ny;
            if (dot > 0) scorePos++;
            else scoreNeg++;
        }

        // The "less crowded side" is where new atoms should go to avoid
        // overlapping the parent ring. Compute it once.
        var lessCrowdedSign = (scorePos <= scoreNeg) ? 1 : -1;

        // If only one path has unplaced atoms (the typical case for a
        // newly-fused bridged ring), place those atoms on the less
        // crowded side. Otherwise (rare: both paths unplaced) put the
        // longer path on the less crowded side and the shorter one
        // (bridge) on the more crowded side.
        if (pathA.length > 0 && pathB.length === 0) {
            placeBridgePath(mol, ring, pathA, bh1, bh2, mx, my, nx, ny,
                            lessCrowdedSign, placed);
        } else if (pathB.length > 0 && pathA.length === 0) {
            placeBridgePath(mol, ring, pathB, bh1, bh2, mx, my, nx, ny,
                            lessCrowdedSign, placed);
        } else {
            // Both paths have atoms — pick longer for less crowded side
            var bridgePath, mainPath;
            if (pathA.length <= pathB.length) {
                bridgePath = pathA; mainPath = pathB;
            } else {
                bridgePath = pathB; mainPath = pathA;
            }
            placeBridgePath(mol, ring, mainPath, bh1, bh2, mx, my, nx, ny,
                            lessCrowdedSign, placed);
            placeBridgePath(mol, ring, bridgePath, bh1, bh2, mx, my, nx, ny,
                            -lessCrowdedSign, placed);
        }
    }

    /**
     * Place atoms along a bridge path as an arc on one side of the
     * bridgehead axis. The arc goes from bridgehead 1 to bridgehead 2,
     * offset by `sign` in the normal direction. The perpendicular bow
     * starts from BOND_LENGTH * 0.8 and is iteratively increased if any
     * placed bridge atom collides with an already-placed atom (so bridges
     * naturally clear existing structure).
     */
    function placeBridgePath(mol, ring, pathIndices, bh1, bh2, mx, my,
                             nx, ny, sign, placed) {
        if (pathIndices.length === 0) return;

        var n = pathIndices.length;
        var axLen = dist(bh1.x, bh1.y, bh2.x, bh2.y);
        if (axLen < EPSILON) axLen = BOND_LENGTH;

        // Initial bow height: aim for hexagonal geometry when axLen ≈ BL√3
        // (bridged 6-6 with 3-atom shared path). Otherwise use a small bow.
        var initialBow = BOND_LENGTH * 0.8;
        if (axLen > BOND_LENGTH * 1.5 && axLen < BOND_LENGTH * 2.0 && n === 3) {
            // Bow that gives hexagon: half-height = BL * sqrt(3)/2 / 2 ≈ 13
            initialBow = BOND_LENGTH * 0.866;
        }

        var maxBow = initialBow;
        for (var attempt = 0; attempt < 4; attempt++) {
            var collision = false;
            for (var i = 0; i < n; i++) {
                var t = (i + 1) / (n + 1);
                var ax = bh1.x + (bh2.x - bh1.x) * t;
                var ay = bh1.y + (bh2.y - bh1.y) * t;
                var bow = 4 * t * (1 - t);
                var offset = maxBow * bow * sign;
                var atom = mol.getAtom(ring[pathIndices[i]]);
                atom.x = ax + nx * offset;
                atom.y = ay + ny * offset;

                // Check collision with already-placed atoms (excluding the
                // bridgeheads themselves, which are the start/end of arc).
                if (attempt < 3) {
                    for (var pid in placed) {
                        if (!placed.hasOwnProperty(pid)) continue;
                        var pidNum = parseInt(pid, 10);
                        if (pidNum === ring[pathIndices[i]]) continue;
                        // Skip bridgeheads (their position is by definition
                        // adjacent to the arc endpoints).
                        var pa = mol.getAtom(pidNum);
                        if (!pa) continue;
                        var dx2 = pa.x - atom.x, dy2 = pa.y - atom.y;
                        var d2 = dx2 * dx2 + dy2 * dy2;
                        if (d2 < BOND_LENGTH * BOND_LENGTH * 0.6 * 0.6) {
                            collision = true;
                            break;
                        }
                    }
                    if (collision) break;
                }
            }
            for (var i = 0; i < n; i++) {
                placed[ring[pathIndices[i]]] = true;
            }
            if (!collision) break;
            maxBow *= 1.4;
        }
    }

    // -----------------------------------------------------------------
    // Spiro junction: place ring sharing a single atom
    // -----------------------------------------------------------------

    function placeSpiroRing(mol, ring, spiroAtomId, placed) {
        var spiroAtom = mol.getAtom(spiroAtomId);
        var size   = ring.length;
        var radius = BOND_LENGTH / (2 * Math.sin(Math.PI / size));

        // Macrocycle tightening for spiro macrocycles
        if (size >= 12) {
            var t = Math.min((size - 12) / 18, 1.0);
            radius *= 0.85 - 0.30 * t;
        }

        // For a proper spiro junction, the second ring should be placed
        // perpendicular to the first ring's bonds at the shared atom.
        // Find the two existing bonds from the spiro atom in the first ring.
        var neighbors = mol.getNeighbors(spiroAtomId);
        var placedNbrs = [];
        for (var i = 0; i < neighbors.length; i++) {
            if (placed[neighbors[i]]) {
                placedNbrs.push(mol.getAtom(neighbors[i]));
            }
        }

        var awayAngle;
        if (placedNbrs.length >= 2) {
            // Compute the bisector angle of the two placed bonds at the
            // spiro atom, then place the new ring perpendicular to it.
            // The bisector points "into" the first ring; the perpendicular
            // gives the direction for the new ring's axis.
            var a1 = Math.atan2(placedNbrs[0].y - spiroAtom.y,
                                placedNbrs[0].x - spiroAtom.x);
            var a2 = Math.atan2(placedNbrs[1].y - spiroAtom.y,
                                placedNbrs[1].x - spiroAtom.x);
            // Average angle via unit-circle addition
            var bisector = Math.atan2(
                Math.sin(a1) + Math.sin(a2),
                Math.cos(a1) + Math.cos(a2)
            );
            // Perpendicular to bisector, pointing away from the first ring
            // Try both perpendicular directions; pick the one farther from
            // placed atoms
            var perp1 = bisector + Math.PI / 2;
            var perp2 = bisector - Math.PI / 2;
            var cx1 = spiroAtom.x + radius * Math.cos(perp1);
            var cy1 = spiroAtom.y + radius * Math.sin(perp1);
            var cx2 = spiroAtom.x + radius * Math.cos(perp2);
            var cy2 = spiroAtom.y + radius * Math.sin(perp2);

            var s1 = 0, s2 = 0;
            for (var id in placed) {
                if (!placed[id]) continue;
                var pa = mol.getAtom(parseInt(id));
                if (!pa) continue;
                s1 += distSq(cx1, cy1, pa.x, pa.y);
                s2 += distSq(cx2, cy2, pa.x, pa.y);
            }
            awayAngle = s1 >= s2 ? perp1 : perp2;
        } else if (placedNbrs.length === 1) {
            // Only one placed neighbour: go directly opposite
            awayAngle = Math.atan2(
                spiroAtom.y - placedNbrs[0].y,
                spiroAtom.x - placedNbrs[0].x
            );
        } else {
            awayAngle = 0;
        }

        var cx = spiroAtom.x + radius * Math.cos(awayAngle);
        var cy = spiroAtom.y + radius * Math.sin(awayAngle);

        var spiroIdx   = ring.indexOf(spiroAtomId);
        var startAngle = Math.atan2(spiroAtom.y - cy, spiroAtom.x - cx);
        var step       = TWO_PI / size;

        for (var i = 0; i < size; i++) {
            if (!placed[ring[i]]) {
                var steps = (i - spiroIdx + size) % size;
                var angle = startAngle + steps * step;
                var atom  = mol.getAtom(ring[i]);
                atom.x = cx + radius * Math.cos(angle);
                atom.y = cy + radius * Math.sin(angle);
                placed[ring[i]] = true;
            }
        }
    }

    // =====================================================================
    // Chain Layout  (120-degree zigzag, extending away from rings)
    // =====================================================================

    function layoutChains(mol, atomIds, placed, ringAtomSet, atomToRings, rings,
                          deferredRingSystems, ringSystems) {
        var maxIter = atomIds.length * 3 + 30;
        var iter = 0;

        // Build lookup: atomId -> deferred ring system index
        var atomToDeferredSys = {};
        if (deferredRingSystems && ringSystems) {
            for (var di = 0; di < deferredRingSystems.length; di++) {
                var sysAtoms = collectSystemAtoms(rings, ringSystems[deferredRingSystems[di]]);
                for (var sa = 0; sa < sysAtoms.length; sa++) {
                    atomToDeferredSys[sysAtoms[sa]] = deferredRingSystems[di];
                }
            }
        }
        var placedDeferredSys = {};

        while (iter++ < maxIter) {
            var anyPlaced = false;

            for (var i = 0; i < atomIds.length; i++) {
                var atomId = atomIds[i];
                if (!placed[atomId]) continue;

                var neighbors = mol.getNeighbors(atomId);
                var unplaced  = [];
                for (var j = 0; j < neighbors.length; j++) {
                    if (!placed[neighbors[j]]) unplaced.push(neighbors[j]);
                }
                if (unplaced.length === 0) continue;

                // Check if any unplaced neighbor is a ring atom — place its ring
                for (var ui = 0; ui < unplaced.length; ui++) {
                    var unplacedId = unplaced[ui];
                    if (!ringAtomSet[unplacedId] || placed[unplacedId]) continue;

                    // Find the smallest ring containing this atom
                    var ringIdx = -1;
                    var ringSize = 99999;
                    for (var ri = 0; ri < rings.length; ri++) {
                        if (rings[ri].indexOf(unplacedId) >= 0 && rings[ri].length < ringSize) {
                            ringIdx = ri;
                            ringSize = rings[ri].length;
                        }
                    }
                    if (ringIdx < 0) continue;

                    var ring = rings[ringIdx];
                    // Skip if any atom in this ring is already placed (fused ring handled elsewhere)
                    var anyRingPlaced = false;
                    for (var rci = 0; rci < ring.length; rci++) {
                        if (placed[ring[rci]]) { anyRingPlaced = true; break; }
                    }
                    if (anyRingPlaced) continue;

                    // Place ring as polygon at parent atom's position
                    var parentAtom = mol.getAtom(atomId);
                    var pNbrs = mol.getNeighbors(atomId);
                    var avgAng = 0, nPl = 0;
                    for (var pn = 0; pn < pNbrs.length; pn++) {
                        if (placed[pNbrs[pn]] && pNbrs[pn] !== unplacedId) {
                            var pa = mol.getAtom(pNbrs[pn]);
                            avgAng += Math.atan2(pa.y - parentAtom.y, pa.x - parentAtom.x);
                            nPl++;
                        }
                    }
                    var awayAngle = (nPl > 0) ? avgAng / nPl + Math.PI : 0;
                    var cx = parentAtom.x + BOND_LENGTH * Math.cos(awayAngle);
                    var cy = parentAtom.y + BOND_LENGTH * Math.sin(awayAngle);

                    placeRingAsPolygon(mol, ring, cx, cy, awayAngle + Math.PI, placed);

                    // Mark deferred system as placed if applicable
                    var defSysIdx = atomToDeferredSys[unplacedId];
                    if (defSysIdx !== undefined) placedDeferredSys[defSysIdx] = true;
                    anyPlaced = true;
                }

                // Recompute unplaced after ring placement
                unplaced = [];
                for (var j = 0; j < neighbors.length; j++) {
                    if (!placed[neighbors[j]]) unplaced.push(neighbors[j]);
                }
                if (unplaced.length === 0) continue;

                var atom = mol.getAtom(atomId);

                // Compute average angle of already-placed neighbours
                var placedNbrs = [];
                for (var j = 0; j < neighbors.length; j++) {
                    if (placed[neighbors[j]]) placedNbrs.push(neighbors[j]);
                }

                var refAngle = computeReferenceAngle(mol, atom, placedNbrs);

                if (unplaced.length === 1) {
                    // Single extension — standard zigzag
                    var alt = chooseZigzagDirection(mol, atomId, unplaced[0], refAngle, placed);
                    placeChainAtom(mol, atomId, unplaced[0], refAngle, placed, alt);
                    anyPlaced = true;
                } else {
                    // Multiple branches — fan out in the open arc
                    placeMultipleBranches(mol, atomId, unplaced, refAngle, placedNbrs.length, placed, ringAtomSet);
                    anyPlaced = true;
                }
            }

            // Check completion
            var allDone = true;
            for (var i = 0; i < atomIds.length; i++) {
                if (!placed[atomIds[i]]) { allDone = false; break; }
            }
            if (allDone || !anyPlaced) break;
        }

        // Safety net: any remaining unplaced atoms
        for (var i = 0; i < atomIds.length; i++) {
            if (!placed[atomIds[i]]) {
                var atom = mol.getAtom(atomIds[i]);
                atom.x = 0; atom.y = 0;
                placed[atomIds[i]] = true;
            }
        }
    }

    /**
     * Compute average direction from `atom` toward its placed neighbours.
     */
    function computeReferenceAngle(mol, atom, placedNbrIds) {
        if (placedNbrIds.length === 0) return Math.PI; // default: came from the left

        var sumSin = 0, sumCos = 0;
        for (var j = 0; j < placedNbrIds.length; j++) {
            var pn = mol.getAtom(placedNbrIds[j]);
            var a  = Math.atan2(pn.y - atom.y, pn.x - atom.x);
            sumSin += Math.sin(a);
            sumCos += Math.cos(a);
        }
        return Math.atan2(sumSin, sumCos);
    }

    /**
     * Decide zigzag direction (+1 or -1) so the new atom goes to the less
     * crowded side and avoids folding back on itself.
     *
     * Strategy:
     *   1. Compute candidate positions on both sides of the incoming bond
     *   2. Penalise candidates that fold back toward the grandparent
     *      (detected by the candidate being closer to the grandparent than
     *      the parent is — indicates the chain is doubling back)
     *   3. Among non-folding candidates, prefer the one farther from the
     *      nearest placed atom (less crowded side)
     */
    function chooseZigzagDirection(mol, parentId, childId, refAngle, placed) {
        var parent = mol.getAtom(parentId);
        var angleA = refAngle + DEG120;
        var angleB = refAngle - DEG120;
        var xA = parent.x + BOND_LENGTH * Math.cos(angleA);
        var yA = parent.y + BOND_LENGTH * Math.sin(angleA);
        var xB = parent.x + BOND_LENGTH * Math.cos(angleB);
        var yB = parent.y + BOND_LENGTH * Math.sin(angleB);

        // Find the grandparent: the placed neighbour of the parent that
        // we came from (closest to the refAngle direction)
        var grandparent = null;
        var parentNbrs = mol.getNeighbors(parentId);
        var bestDot = -Infinity;
        var refDx = Math.cos(refAngle), refDy = Math.sin(refAngle);
        for (var i = 0; i < parentNbrs.length; i++) {
            if (placed[parentNbrs[i]] && parentNbrs[i] !== childId) {
                var gp = mol.getAtom(parentNbrs[i]);
                var gpDx = gp.x - parent.x, gpDy = gp.y - parent.y;
                var gpLen = Math.sqrt(gpDx * gpDx + gpDy * gpDy);
                if (gpLen > EPSILON) {
                    var dot = (gpDx / gpLen) * refDx + (gpDy / gpLen) * refDy;
                    if (dot > bestDot) { bestDot = dot; grandparent = gp; }
                }
            }
        }

        var minDistA = Infinity, minDistB = Infinity;
        var foldPenaltyA = 0, foldPenaltyB = 0;

        // Deep anti-fold-back: walk back up to 4 ancestors and penalise
        // candidates that come too close to any of them (closer ancestors
        // receive a heavier penalty).
        var ancestor = grandparent;
        var ancestorId = null;
        if (grandparent) {
            // Find the grandparent's atom id
            for (var pi = 0; pi < parentNbrs.length; pi++) {
                if (placed[parentNbrs[pi]] && parentNbrs[pi] !== childId) {
                    var gp = mol.getAtom(parentNbrs[pi]);
                    if (gp === grandparent) { ancestorId = parentNbrs[pi]; break; }
                }
            }
        }
        var curAncestor = ancestor;
        var curAncestorId = ancestorId;
        var prevAncestorId = parentId;
        for (var depth = 0; depth < 4 && curAncestor; depth++) {
            var threshold = (0.5 + 0.1 * depth) * BOND_LENGTH;
            var threshSq = threshold * threshold;
            var dA = distSq(xA, yA, curAncestor.x, curAncestor.y);
            var dB = distSq(xB, yB, curAncestor.x, curAncestor.y);
            if (dA < threshSq) foldPenaltyA += 1e6 / (depth + 1);
            if (dB < threshSq) foldPenaltyB += 1e6 / (depth + 1);

            // Walk to the next ancestor
            var nextAncestor = null;
            var nextAncestorId = null;
            if (curAncestorId !== null) {
                var ancNbrs = mol.getNeighbors(curAncestorId);
                for (var ni = 0; ni < ancNbrs.length; ni++) {
                    if (placed[ancNbrs[ni]] && ancNbrs[ni] !== prevAncestorId) {
                        nextAncestor = mol.getAtom(ancNbrs[ni]);
                        nextAncestorId = ancNbrs[ni];
                        break;
                    }
                }
            }
            prevAncestorId = curAncestorId;
            curAncestor = nextAncestor;
            curAncestorId = nextAncestorId;
        }

        for (var id in placed) {
            if (!placed[id]) continue;
            var pa = mol.getAtom(parseInt(id));
            if (!pa || parseInt(id) === parentId) continue;
            var dA = distSq(xA, yA, pa.x, pa.y);
            var dB = distSq(xB, yB, pa.x, pa.y);
            if (dA < minDistA) minDistA = dA;
            if (dB < minDistB) minDistB = dB;
        }

        // Subtract penalties (lower score = worse)
        var scoreA = minDistA - foldPenaltyA;
        var scoreB = minDistB - foldPenaltyB;
        return scoreA >= scoreB ? 1 : -1;
    }

    /**
     * Place a single chain atom at 120 degrees from the incoming direction.
     */
    function placeChainAtom(mol, parentId, childId, refAngle, placed, alternation) {
        var parent = mol.getAtom(parentId);
        var child  = mol.getAtom(childId);
        var angle  = refAngle + DEG120 * alternation;

        child.x = parent.x + BOND_LENGTH * Math.cos(angle);
        child.y = parent.y + BOND_LENGTH * Math.sin(angle);
        placed[childId] = true;
    }

    /**
     * Place multiple substituent branches fanning out in the open arc
     * (the semicircle opposite the existing bonds).
     */
    function placeMultipleBranches(mol, parentId, unplacedIds, refAngle, numPlaced, placed, ringAtomSet) {
        var parent = mol.getAtom(parentId);
        var total  = unplacedIds.length;

        if (numPlaced === 0) {
            // Root atom: spread evenly around full circle
            var step = TWO_PI / total;
            for (var i = 0; i < total; i++) {
                var angle = i * step;
                var child = mol.getAtom(unplacedIds[i]);
                child.x = parent.x + BOND_LENGTH * Math.cos(angle);
                child.y = parent.y + BOND_LENGTH * Math.sin(angle);
                placed[unplacedIds[i]] = true;
            }
            return;
        }

        // Direction away from existing bonds (outward)
        var outAngle = refAngle + Math.PI;

        // Available arc: ideally 240 degrees for 1 existing bond,
        // shrinking as more bonds exist. Each branch gets DEG120.
        // When the parent is a ring atom, widen the fan to DEG120 * 1.3
        // to reduce overlap with substituents on adjacent ring atoms
        // (e.g. aspirin's 1,2-disubstituted benzene).
        var isRingParent = !!ringAtomSet[parentId];
        var arcPerBranch = DEG120;
        if (isRingParent && numPlaced >= 2) {
            arcPerBranch = DEG120 * 1.4;
        }
        if (total > 3) arcPerBranch = Math.PI / (total - 1 + 0.5);
        var totalArc   = arcPerBranch * (total - 1);
        var startAngle = outAngle - totalArc / 2;

        // Sort branches: put ring-attached atoms to the outside, hydrogens inside
        // (simple heuristic: heavier substituents get the outermost positions)
        for (var i = 0; i < total; i++) {
            var angle = startAngle + i * arcPerBranch;
            var child = mol.getAtom(unplacedIds[i]);
            child.x = parent.x + BOND_LENGTH * Math.cos(angle);
            child.y = parent.y + BOND_LENGTH * Math.sin(angle);
            placed[unplacedIds[i]] = true;
        }
    }

    // =====================================================================
    // Bond-Length Normalisation
    // =====================================================================

    /**
     * Scale the entire component so that the median bond length equals
     * BOND_LENGTH.  This corrects any distortion from ring fusion.
     */
    function normaliseBondLengths(mol, atomIds) {
        var idSet = {};
        for (var i = 0; i < atomIds.length; i++) idSet[atomIds[i]] = true;

        var lengths = [];
        for (var i = 0; i < atomIds.length; i++) {
            var a  = mol.getAtom(atomIds[i]);
            var nb = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < nb.length; j++) {
                if (idSet[nb[j]] && nb[j] > atomIds[i]) {
                    var b = mol.getAtom(nb[j]);
                    lengths.push(dist(a.x, a.y, b.x, b.y));
                }
            }
        }
        if (lengths.length === 0) return;

        lengths.sort(function (a, b) { return a - b; });
        var median = lengths[Math.floor(lengths.length / 2)];
        if (median < EPSILON) return;

        var scale = BOND_LENGTH / median;
        // Only rescale if significantly off (>2% deviation)
        if (Math.abs(scale - 1) < 0.02) return;

        // Scale around centroid
        var cx = 0, cy = 0;
        for (var i = 0; i < atomIds.length; i++) {
            var a = mol.getAtom(atomIds[i]);
            cx += a.x; cy += a.y;
        }
        cx /= atomIds.length; cy /= atomIds.length;

        for (var i = 0; i < atomIds.length; i++) {
            var a = mol.getAtom(atomIds[i]);
            a.x = cx + (a.x - cx) * scale;
            a.y = cy + (a.y - cy) * scale;
        }
    }

    // =====================================================================
    // Collision Detection & Resolution  (spatial-hash grid, O(n) per pass)
    // =====================================================================

    /**
     * Iteratively push overlapping atoms apart.
     * Ring atoms are considered fixed; chain atoms are moved preferentially.
     * Uses a spatial hash grid so each pass is O(n) rather than O(n^2).
     */
    function resolveCollisions(mol, atomIds, ringAtomSet) {
        var n = atomIds.length;
        if (n < 2) return;

        // Snapshot ring atom positions so we can restore them after overlap
        // resolution. This prevents cumulative distortion of ring geometry
        // when many substituents are attached (e.g., steroids with 4 fused
        // rings and numerous substituents).
        var ringSnapshot = {};
        for (var i = 0; i < n; i++) {
            if (ringAtomSet[atomIds[i]]) {
                var ra = mol.getAtom(atomIds[i]);
                ringSnapshot[atomIds[i]] = { x: ra.x, y: ra.y };
            }
        }

        for (var iter = 0; iter < MAX_OVERLAP_ITER; iter++) {
            var collision = false;

            // Build spatial hash
            var grid = {};
            for (var i = 0; i < n; i++) {
                var a  = mol.getAtom(atomIds[i]);
                var gx = Math.floor(a.x / GRID_CELL);
                var gy = Math.floor(a.y / GRID_CELL);
                var key = gx + ',' + gy;
                if (!grid[key]) grid[key] = [];
                grid[key].push(i);
            }

            // Check each atom against neighbours in the 3x3 grid neighbourhood
            for (var i = 0; i < n; i++) {
                var a1 = mol.getAtom(atomIds[i]);
                var gx = Math.floor(a1.x / GRID_CELL);
                var gy = Math.floor(a1.y / GRID_CELL);

                for (var dx = -1; dx <= 1; dx++) {
                    for (var dy = -1; dy <= 1; dy++) {
                        var cell = grid[(gx + dx) + ',' + (gy + dy)];
                        if (!cell) continue;
                        for (var ci = 0; ci < cell.length; ci++) {
                            var j = cell[ci];
                            if (j <= i) continue; // avoid double processing

                            var a2 = mol.getAtom(atomIds[j]);
                            var ddx = a2.x - a1.x, ddy = a2.y - a1.y;
                            var d = Math.sqrt(ddx * ddx + ddy * ddy);

                            if (d < MIN_ATOM_DIST) {
                                collision = true;
                                if (d < EPSILON) {
                                    // Atoms perfectly co-located: pick a deterministic
                                    // pseudo-random direction (avoid divide-by-EPSILON
                                    // which produces ±1e6 push vectors and explodes
                                    // coords to millions). Seed from the index pair so
                                    // different colliding pairs separate in different
                                    // directions.
                                    var __ang = (i * 0.7137 + j * 1.31) * Math.PI;
                                    ddx = Math.cos(__ang);
                                    ddy = Math.sin(__ang);
                                    d = 1;
                                }
                                var push = (MIN_ATOM_DIST - d) / 2 + 0.5;
                                var pnx = ddx / d, pny = ddy / d;

                                var r1 = !!ringAtomSet[atomIds[i]];
                                var r2 = !!ringAtomSet[atomIds[j]];

                                if (r1 && r2) {
                                    // Both ring atoms: apply a gentle push to
                                    // resolve bridged-ring overlaps (e.g. morphine
                                    // oxygen bridge) without large distortion.
                                    var gentlePush = push * 0.3;
                                    a1.x -= pnx * gentlePush;
                                    a1.y -= pny * gentlePush;
                                    a2.x += pnx * gentlePush;
                                    a2.y += pny * gentlePush;
                                } else if (r1) {
                                    // Only a2 (chain) moves — ring atom stays fixed
                                    a2.x += pnx * push;
                                    a2.y += pny * push;
                                } else if (r2) {
                                    // Only a1 (chain) moves — ring atom stays fixed
                                    a1.x -= pnx * push;
                                    a1.y -= pny * push;
                                } else {
                                    // Both chain: split evenly
                                    a1.x -= pnx * push;
                                    a1.y -= pny * push;
                                    a2.x += pnx * push;
                                    a2.y += pny * push;
                                }
                            }
                        }
                    }
                }
            }

            if (!collision) break;
        }

        // Restore ring atom positions unless they were pushed to resolve
        // a ring-ring overlap (bridged systems like morphine, norbornane).
        for (var id in ringSnapshot) {
            var atom = mol.getAtom(parseInt(id));
            var sx = ringSnapshot[id].x, sy = ringSnapshot[id].y;
            var dx2 = atom.x - sx, dy2 = atom.y - sy;
            var wasMoved = (dx2 * dx2 + dy2 * dy2) > 0.01;
            if (wasMoved) {
                // Atom was pushed during collision resolution — keep new position
            } else {
                atom.x = sx;
                atom.y = sy;
            }
        }
    }

    // =====================================================================
    // Orientation Optimization
    // =====================================================================

    /**
     * Rotate the component to maximise bond alignment with the horizontal
     * and vertical axes, preferring a landscape (wider-than-tall) orientation.
     */
    function optimizeOrientation(mol, atomIds) {
        var n = atomIds.length;
        if (n < 2) return;

        // Compute centroid
        var cx = 0, cy = 0;
        for (var i = 0; i < n; i++) {
            var a = mol.getAtom(atomIds[i]);
            cx += a.x; cy += a.y;
        }
        cx /= n; cy /= n;

        // Collect bond pairs (avoid duplicates)
        var idSet = {};
        for (var i = 0; i < n; i++) idSet[atomIds[i]] = true;
        var bonds = [];
        for (var i = 0; i < n; i++) {
            var nbrs = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < nbrs.length; j++) {
                if (idSet[nbrs[j]] && nbrs[j] > atomIds[i]) {
                    bonds.push([atomIds[i], nbrs[j]]);
                }
            }
        }
        if (bonds.length === 0) return;

        // Score function: sum of cos(2 * bondAngle) + landscape bonus
        function scoreRotation(theta) {
            var s = 0;
            var cosT = Math.cos(theta), sinT = Math.sin(theta);
            var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;

            for (var bi = 0; bi < bonds.length; bi++) {
                var a1 = mol.getAtom(bonds[bi][0]);
                var a2 = mol.getAtom(bonds[bi][1]);
                // Rotate both endpoints around centroid
                var x1 = cx + (a1.x - cx) * cosT - (a1.y - cy) * sinT;
                var y1 = cy + (a1.x - cx) * sinT + (a1.y - cy) * cosT;
                var x2 = cx + (a2.x - cx) * cosT - (a2.y - cy) * sinT;
                var y2 = cy + (a2.x - cx) * sinT + (a2.y - cy) * cosT;
                var bondAngle = Math.atan2(y2 - y1, x2 - x1);
                s += Math.cos(2 * bondAngle);
            }

            // Compute bounding box for landscape bonus
            for (var ai = 0; ai < n; ai++) {
                var a = mol.getAtom(atomIds[ai]);
                var rx = cx + (a.x - cx) * cosT - (a.y - cy) * sinT;
                var ry = cy + (a.x - cx) * sinT + (a.y - cy) * cosT;
                if (rx < minX) minX = rx;
                if (rx > maxX) maxX = rx;
                if (ry < minY) minY = ry;
                if (ry > maxY) maxY = ry;
            }
            var width = maxX - minX, height = maxY - minY;
            var dim = Math.max(width, height);
            if (dim > EPSILON) {
                s += 0.5 * (width - height) / dim;
            }
            return s;
        }

        var initialScore = scoreRotation(0);
        var bestScore = initialScore;
        var bestTheta = 0;

        // Try 12 rotations (every 30 degrees)
        for (var k = 1; k < 12; k++) {
            var theta = k * Math.PI / 6;
            var s = scoreRotation(theta);
            if (s > bestScore) {
                bestScore = s;
                bestTheta = theta;
            }
        }

        // Only apply if meaningfully better than initial
        if (bestTheta !== 0 && bestScore > initialScore * 1.1) {
            var cosB = Math.cos(bestTheta), sinB = Math.sin(bestTheta);
            for (var i = 0; i < n; i++) {
                var a = mol.getAtom(atomIds[i]);
                var dx = a.x - cx, dy = a.y - cy;
                a.x = cx + dx * cosB - dy * sinB;
                a.y = cy + dx * sinB + dy * cosB;
            }
        }
    }

    // =====================================================================
    // Light Force-Directed Refinement
    // =====================================================================

    /**
     * Apply a short force-directed pass to smooth chain atom positions.
     * Ring atoms are frozen to preserve ring geometry.
     */
    function refineLayout(mol, atomIds, ringAtomSet) {
        var n = atomIds.length;
        if (n < 3) return;

        // Build set and index
        var idSet = {};
        var idxOf = {};
        for (var i = 0; i < n; i++) {
            idSet[atomIds[i]] = true;
            idxOf[atomIds[i]] = i;
        }

        // Identify ring atoms (frozen)
        var isRing = [];
        for (var i = 0; i < n; i++) {
            isRing[i] = !!ringAtomSet[atomIds[i]];
        }

        // Build adjacency for bonds within this component
        var bonded = [];
        for (var i = 0; i < n; i++) bonded[i] = [];
        for (var i = 0; i < n; i++) {
            var nbrs = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < nbrs.length; j++) {
                if (idSet[nbrs[j]] && idxOf[nbrs[j]] > i) {
                    bonded[i].push(idxOf[nbrs[j]]);
                    bonded[idxOf[nbrs[j]]].push(i);
                }
            }
        }

        // Build bonded-pair lookup for quick checking
        var bondedSet = [];
        for (var i = 0; i < n; i++) {
            bondedSet[i] = {};
            for (var j = 0; j < bonded[i].length; j++) {
                bondedSet[i][bonded[i][j]] = true;
            }
        }

        var maxDisp = 0.3 * BOND_LENGTH;
        var temp = 1.0;
        var repulsionRange = 2 * BOND_LENGTH;
        var repulsionRangeSq = repulsionRange * repulsionRange;

        for (var iter = 0; iter < 30; iter++) {
            var fx = [], fy = [];
            for (var i = 0; i < n; i++) { fx[i] = 0; fy[i] = 0; }

            for (var i = 0; i < n; i++) {
                if (isRing[i]) continue;
                var ai = mol.getAtom(atomIds[i]);

                // Spring forces toward bonded neighbours
                for (var bi = 0; bi < bonded[i].length; bi++) {
                    var j = bonded[i][bi];
                    var aj = mol.getAtom(atomIds[j]);
                    var dx = aj.x - ai.x, dy = aj.y - ai.y;
                    var d = Math.sqrt(dx * dx + dy * dy);
                    if (d < EPSILON) continue;
                    var f = (d - BOND_LENGTH) * 0.1;
                    fx[i] += f * dx / d;
                    fy[i] += f * dy / d;
                }

                // Repulsion from non-bonded atoms within range
                for (var j = 0; j < n; j++) {
                    if (j === i || bondedSet[i][j]) continue;
                    var aj = mol.getAtom(atomIds[j]);
                    var dx = aj.x - ai.x, dy = aj.y - ai.y;
                    var dSq = dx * dx + dy * dy;
                    if (dSq > repulsionRangeSq || dSq < EPSILON * EPSILON) continue;
                    var d = Math.sqrt(dSq);
                    var f = -0.5 * BOND_LENGTH / Math.max(d, 1);
                    fx[i] += f * dx / d;
                    fy[i] += f * dy / d;
                }
            }

            // Apply forces with temperature decay
            for (var i = 0; i < n; i++) {
                if (isRing[i]) continue;
                var fMag = Math.sqrt(fx[i] * fx[i] + fy[i] * fy[i]);
                if (fMag < EPSILON) continue;
                var disp = Math.min(fMag * temp, maxDisp);
                var ai = mol.getAtom(atomIds[i]);
                ai.x += (fx[i] / fMag) * disp;
                ai.y += (fy[i] / fMag) * disp;
            }

            temp *= 0.95;
        }
    }

    // =====================================================================
    // Utility functions
    // =====================================================================

    function dist(x1, y1, x2, y2) {
        var dx = x2 - x1, dy = y2 - y1;
        return Math.sqrt(dx * dx + dy * dy);
    }

    function distSq(x1, y1, x2, y2) {
        var dx = x2 - x1, dy = y2 - y1;
        return dx * dx + dy * dy;
    }

    function normalizeAngle(a) {
        // FIX: guard against NaN/Infinity which would cause infinite loop
        if (!isFinite(a)) return 0;
        while (a < 0)       a += TWO_PI;
        while (a >= TWO_PI) a -= TWO_PI;
        return a;
    }

    function atomBounds(atoms) {
        var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        for (var i = 0; i < atoms.length; i++) {
            if (atoms[i].x < minX) minX = atoms[i].x;
            if (atoms[i].y < minY) minY = atoms[i].y;
            if (atoms[i].x > maxX) maxX = atoms[i].x;
            if (atoms[i].y > maxY) maxY = atoms[i].y;
        }
        return { minX: minX, minY: minY, maxX: maxX, maxY: maxY };
    }

    function uniqueArray(arr) {
        var seen = {}, result = [];
        for (var i = 0; i < arr.length; i++) {
            if (!seen[arr[i]]) { seen[arr[i]] = true; result.push(arr[i]); }
        }
        return result;
    }

    // NOTE: newIntArray was dead code (never called) — removed to reduce bundle size

    // =====================================================================
    // Ring Query Helpers  (used by Renderer for aromatic circles etc.)
    // =====================================================================

    /**
     * Returns ring info for all rings in the molecule.
     * Each entry: { atoms: [atomIds], center: {x,y}, size, aromatic }
     *
     * Aromatic detection uses:
     *   1. atom.aromatic flags from the SMILES parser (preferred)
     *   2. Heuristic: 5- or 6-membered ring with at least one double bond
     */
    Layout.getRingInfo = function (mol) {
        var components = mol.getComponents();
        var allRings   = [];

        for (var ci = 0; ci < components.length; ci++) {
            var rings = perceiveSSSR(mol, components[ci]);

            for (var ri = 0; ri < rings.length; ri++) {
                var ring = rings[ri];
                var cx = 0, cy = 0;
                var allAroFlag   = true;   // all atoms flagged aromatic?
                var hasDouble    = false;

                for (var ai = 0; ai < ring.length; ai++) {
                    var atom = mol.getAtom(ring[ai]);
                    cx += atom.x;
                    cy += atom.y;
                    if (!atom.aromatic) allAroFlag = false;
                }
                cx /= ring.length;
                cy /= ring.length;

                // Check bond types around the ring
                for (var ai = 0; ai < ring.length; ai++) {
                    var bond = mol.getBondBetween(ring[ai], ring[(ai + 1) % ring.length]);
                    if (bond && bond.type === Molecule.BOND_DOUBLE) hasDouble = true;
                }

                // Aromatic if SMILES flagged it, or heuristic 5/6-ring + double bond
                var isAromatic = allAroFlag ||
                    (hasDouble && (ring.length === 5 || ring.length === 6));

                allRings.push({
                    atoms:    ring,
                    center:   { x: cx, y: cy },
                    size:     ring.length,
                    aromatic: isAromatic
                });
            }
        }

        return allRings;
    };

    /**
     * Check if a bond is part of a ring.  Returns the ring info or null.
     */
    Layout.bondInRing = function (mol, bond) {
        var ringInfo = Layout.getRingInfo(mol);
        for (var i = 0; i < ringInfo.length; i++) {
            var ring = ringInfo[i].atoms;
            var idx1 = ring.indexOf(bond.atom1);
            var idx2 = ring.indexOf(bond.atom2);
            if (idx1 >= 0 && idx2 >= 0) {
                var diff = Math.abs(idx1 - idx2);
                if (diff === 1 || diff === ring.length - 1) return ringInfo[i];
            }
        }
        return null;
    };

    // =====================================================================
    // v1.5.2 Chemistry-Canonical Post-Passes
    //   Designed to run AFTER ring/chain placement so a stretched bond,
    //   a tilted hexagon, an inward substituent, or a 70-deg sp3 angle
    //   gets pulled back to textbook geometry without disturbing the
    //   global topology.
    //
    //   Determinism: every loop iterates over atomIds[] in stable order;
    //   tiebreaks fall back to atom.id (smallest wins). No randomness,
    //   no thread state, no Object.keys ordering used for placement.
    //
    //   Complexity:
    //     relaxBondLengths       : O(I * E)        I=50 iters, E=#bonds
    //     enforceRingOrientation : O(R * S)        R=#rings,    S=ring size
    //     enforceSubstituentDir  : O(R * S * deg)
    //     enforceFusedAxis       : O(R^2 + N)
    //     smoothBondAngles       : O(N * deg^2)
    //   None of these are O(N^2) on shipped corpora (R<=20 typical).
    // =====================================================================

    /**
     * Step 1: Iterative bond-length relaxation.
     * For each bond, compute (actual - target) and apply half the
     * correction along the bond axis to each endpoint. Ring atoms are
     * weighted heavier so the regular polygon geometry is preserved.
     * Converges when max abs deviation < 0.05 * BOND_LENGTH OR 50 iters.
     */
    function relaxBondLengths(mol, atomIds, ringAtomSet) {
        var n = atomIds.length;
        if (n < 2) return;
        var idSet = {};
        for (var i = 0; i < n; i++) idSet[atomIds[i]] = true;

        // Collect bonds inside this component (deterministic order via atomIds).
        var bonds = [];
        for (var i = 0; i < n; i++) {
            var nb = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < nb.length; j++) {
                if (idSet[nb[j]] && nb[j] > atomIds[i]) {
                    bonds.push([atomIds[i], nb[j]]);
                }
            }
        }
        if (bonds.length === 0) return;

        var TOL = BOND_LENGTH * 0.05;
        var MAX_ITER = 50;
        // Ring atoms move at 0.10x the chain rate so polygon shape stays.
        var RING_DAMP = 0.10;
        var CHAIN_GAIN = 0.5;

        for (var iter = 0; iter < MAX_ITER; iter++) {
            var maxDev = 0;
            // Accumulate displacements; apply at end of pass for stability.
            var dx = {}, dy = {};
            for (var i = 0; i < n; i++) { dx[atomIds[i]] = 0; dy[atomIds[i]] = 0; }

            for (var bi = 0; bi < bonds.length; bi++) {
                var u = bonds[bi][0], v = bonds[bi][1];
                var au = mol.getAtom(u), av = mol.getAtom(v);
                var ddx = av.x - au.x, ddy = av.y - au.y;
                var d = Math.sqrt(ddx * ddx + ddy * ddy);
                if (d < EPSILON) continue;
                var dev = d - BOND_LENGTH;
                if (Math.abs(dev) > maxDev) maxDev = Math.abs(dev);
                var nx = ddx / d, ny = ddy / d;
                var ru = !!ringAtomSet[u], rv = !!ringAtomSet[v];
                var gu = ru ? RING_DAMP : CHAIN_GAIN;
                var gv = rv ? RING_DAMP : CHAIN_GAIN;
                // Each endpoint moves by gain * dev / 2 along the bond axis.
                dx[u] += nx * dev * gu * 0.5;
                dy[u] += ny * dev * gu * 0.5;
                dx[v] -= nx * dev * gv * 0.5;
                dy[v] -= ny * dev * gv * 0.5;
            }

            for (var i = 0; i < n; i++) {
                var a = mol.getAtom(atomIds[i]);
                a.x += dx[atomIds[i]];
                a.y += dy[atomIds[i]];
            }

            if (maxDev < TOL) break;
        }
    }

    /**
     * Step 2: rotate aromatic 6-rings so the bond between their two
     * highest-degree atoms is horizontal (canonical chemistry drawing).
     * For 5-rings, rotate so the heteroatom (or highest-degree atom) is
     * at the top vertex.
     *
     * Only acts on rings whose ring system is a SINGLE ring (not part
     * of a fused system) -- fused systems are handled by Step 4.
     * Skipped if the ring contains atoms also in chain branches that
     * have already been laid out around it (rotating would scramble
     * the chain).
     */
    function enforceRingOrientation(mol, rings, ringSystems, atomIds) {
        if (!rings || rings.length === 0) return;

        // Map ring index -> ring system index (only orient single-ring systems)
        var ringSysSize = {};
        for (var si = 0; si < ringSystems.length; si++) {
            for (var ri = 0; ri < ringSystems[si].length; ri++) {
                ringSysSize[ringSystems[si][ri]] = ringSystems[si].length;
            }
        }

        // Build neighbor degree map for the molecule (deterministic).
        // (Heavy-atom degree, ignoring H if present in the graph.)
        var degree = {};
        for (var i = 0; i < atomIds.length; i++) {
            var nb = mol.getNeighbors(atomIds[i]);
            var d = 0;
            for (var j = 0; j < nb.length; j++) {
                var a = mol.getAtom(nb[j]);
                if (a && a.symbol !== 'H') d++;
            }
            degree[atomIds[i]] = d;
        }

        for (var ri = 0; ri < rings.length; ri++) {
            if ((ringSysSize[ri] || 1) > 1) continue; // fused -- Step 4
            var ring = rings[ri];
            var sz = ring.length;
            if (sz < 5 || sz > 7) continue; // only 5,6,7-rings get oriented

            // Centroid of the ring.
            var cx = 0, cy = 0;
            for (var k = 0; k < sz; k++) {
                var a = mol.getAtom(ring[k]); cx += a.x; cy += a.y;
            }
            cx /= sz; cy /= sz;

            // Choose target orientation.
            // 6-ring: pick the bond whose two endpoints have max sum-of-degree;
            //         tiebreak by smaller (atomId1, atomId2) lex pair.
            //         The midpoint of that bond should be on the +x axis.
            // 5/7-ring: pick the heteroatom (smallest atom-id with non-C
            //         symbol) and put it at the top vertex (y minimum after
            //         rotation, i.e. -PI/2 in our screen-space convention).
            var targetVertexIdx = -1;     // for 5,7
            var targetEdge = null;        // for 6: [idxA, idxB]

            if (sz === 6) {
                // First detect substitution pattern.
                // Para = exactly 2 substituent atoms at positions (i, i+3).
                // Meta = exactly 2 substituent atoms at positions (i, i+2).
                // 1,3,5-tri = 3 substituents at positions (i, i+2, i+4).
                // For each canonical pattern we want a SPECIFIC axis horizontal.
                var subRingPositions = [];
                for (var k = 0; k < sz; k++) {
                    if ((degree[ring[k]] || 0) >= 3) subRingPositions.push(k);
                }
                var paraVertex = -1, oppVertex = -1;
                if (subRingPositions.length === 2 &&
                    Math.abs(subRingPositions[1] - subRingPositions[0]) === 3) {
                    // Para: axis = vertex (k=subRingPositions[0]) -> opposite vertex.
                    paraVertex = subRingPositions[0];
                    oppVertex  = subRingPositions[1];
                }
                if (paraVertex >= 0) {
                    // Mark targetEdge as a synthetic vertex-pair (encode k=-1).
                    targetEdge = ['para', paraVertex, oppVertex];
                } else {
                    // Default: pick highest sum-of-degree edge.
                    var bestScore = -1;
                    for (var k = 0; k < sz; k++) {
                        var u = ring[k], v = ring[(k + 1) % sz];
                        var sc = (degree[u] || 0) + (degree[v] || 0);
                        var lexA = Math.min(u, v), lexB = Math.max(u, v);
                        var entry = [sc, -lexA, -lexB, k];
                        if (targetEdge === null ||
                            entry[0] > targetEdge[0] ||
                            (entry[0] === targetEdge[0] && entry[1] > targetEdge[1]) ||
                            (entry[0] === targetEdge[0] && entry[1] === targetEdge[1] && entry[2] > targetEdge[2])) {
                            targetEdge = entry;
                        }
                    }
                }
            } else {
                // Pick lowest-id non-C heteroatom; if all C, pick highest-degree.
                var bestHetId = -1;
                for (var k = 0; k < sz; k++) {
                    var sym = (mol.getAtom(ring[k]).symbol || 'C');
                    if (sym !== 'C' && (bestHetId === -1 || ring[k] < bestHetId)) {
                        bestHetId = ring[k];
                        targetVertexIdx = k;
                    }
                }
                if (targetVertexIdx === -1) {
                    // No heteroatom; pick max-degree (tiebreak smallest id).
                    var bestDeg = -1, bestId = -1;
                    for (var k = 0; k < sz; k++) {
                        var dd = degree[ring[k]] || 0;
                        if (dd > bestDeg || (dd === bestDeg && ring[k] < bestId)) {
                            bestDeg = dd; bestId = ring[k]; targetVertexIdx = k;
                        }
                    }
                }
            }

            // Compute the rotation needed.
            var rot = 0;
            if (sz === 6 && targetEdge && targetEdge[0] === 'para') {
                // Para: align the two substituent-bearing atoms on the
                // horizontal axis. Equivalently, the line through the two
                // vertices should have angle 0 (or PI; pick the shorter rot).
                var pa = mol.getAtom(ring[targetEdge[1]]);
                var pb = mol.getAtom(ring[targetEdge[2]]);
                var lineAng = Math.atan2(pb.y - pa.y, pb.x - pa.x);
                var r0 = -lineAng;
                while (r0 > Math.PI) r0 -= TWO_PI;
                while (r0 < -Math.PI) r0 += TWO_PI;
                var rPi = Math.PI - lineAng;
                while (rPi > Math.PI) rPi -= TWO_PI;
                while (rPi < -Math.PI) rPi += TWO_PI;
                rot = (Math.abs(r0) <= Math.abs(rPi)) ? r0 : rPi;
            } else if (sz === 6 && targetEdge) {
                var k = targetEdge[3];
                var pa = mol.getAtom(ring[k]);
                var pb = mol.getAtom(ring[(k + 1) % sz]);
                var midAng = Math.atan2((pa.y + pb.y) / 2 - cy, (pa.x + pb.x) / 2 - cx);
                // Want midAng == 0 (right). Compare also PI (left); pick shorter rotation.
                var r0 = -midAng;
                while (r0 > Math.PI) r0 -= TWO_PI;
                while (r0 < -Math.PI) r0 += TWO_PI;
                var rPi = Math.PI - midAng;
                while (rPi > Math.PI) rPi -= TWO_PI;
                while (rPi < -Math.PI) rPi += TWO_PI;
                rot = (Math.abs(r0) <= Math.abs(rPi)) ? r0 : rPi;
            } else if (targetVertexIdx >= 0) {
                var pa = mol.getAtom(ring[targetVertexIdx]);
                var ang = Math.atan2(pa.y - cy, pa.x - cx);
                // Want ang = -PI/2 (top in screen coords where +y is down).
                rot = (-Math.PI / 2) - ang;
                while (rot > Math.PI) rot -= TWO_PI;
                while (rot < -Math.PI) rot += TWO_PI;
            }

            if (Math.abs(rot) < 0.05) continue; // already oriented

            // Rotate the ring AND any chain atoms attached to it (via DFS
            // through non-ring bonds) so the substituent positions move
            // with the ring. We collect the connected component once and
            // rotate every atom in it about (cx, cy).
            // To avoid rotating the WHOLE molecule (which we definitely
            // do not want for a single-ring with substituents on multiple
            // sides), we restrict rotation to the ring atoms and their
            // immediate non-ring neighbours plus their non-ring DFS
            // closure that does NOT come back to a ring.
            var toRotate = {};
            for (var k = 0; k < sz; k++) toRotate[ring[k]] = true;
            // BFS outward through non-ring atoms
            var queue = ring.slice();
            var ringMember = {};
            for (var k = 0; k < sz; k++) ringMember[ring[k]] = true;
            for (var qh = 0; qh < queue.length; qh++) {
                var cur = queue[qh];
                var nb = mol.getNeighbors(cur);
                for (var ni = 0; ni < nb.length; ni++) {
                    var nbId = nb[ni];
                    if (toRotate[nbId]) continue;
                    // Only follow non-ring neighbours OR ring members of the
                    // current ring. Stop at any other ring atom (would pull
                    // a different ring along with us).
                    var nbAtom = mol.getAtom(nbId);
                    if (!nbAtom) continue;
                    // Stop at any atom belonging to a DIFFERENT ring.
                    var stop = false;
                    for (var orI = 0; orI < rings.length; orI++) {
                        if (orI === ri) continue;
                        if (rings[orI].indexOf(nbId) >= 0) { stop = true; break; }
                    }
                    if (stop) continue;
                    toRotate[nbId] = true;
                    queue.push(nbId);
                }
            }

            var c = Math.cos(rot), s = Math.sin(rot);
            for (var k = 0; k < atomIds.length; k++) {
                if (!toRotate[atomIds[k]]) continue;
                var a = mol.getAtom(atomIds[k]);
                var rx = a.x - cx, ry = a.y - cy;
                a.x = cx + rx * c - ry * s;
                a.y = cy + rx * s + ry * c;
            }
        }
    }

    /**
     * Step 3: For every ring atom that bears one or more non-ring
     * substituents, ensure each substituent\u2019s first bond points
     * RADIALLY OUTWARD from the ring centroid (never inward into the
     * ring). When two adjacent ring atoms both bear substituents, the
     * radial direction naturally separates them.
     */
    function enforceSubstituentDirection(mol, rings, ringAtomSet, atomIds) {
        if (!rings || rings.length === 0) return;

        // Compute centroid for each ring (deterministic by ring index).
        var ringCentroids = [];
        for (var ri = 0; ri < rings.length; ri++) {
            var cx = 0, cy = 0, ring = rings[ri];
            for (var k = 0; k < ring.length; k++) {
                var a = mol.getAtom(ring[k]); cx += a.x; cy += a.y;
            }
            ringCentroids.push({ x: cx / ring.length, y: cy / ring.length });
        }

        // For each ring, for each atom in it: find non-ring neighbours and
        // align their first bond along the radial direction.
        for (var ri = 0; ri < rings.length; ri++) {
            var ring = rings[ri];
            var rc = ringCentroids[ri];
            for (var k = 0; k < ring.length; k++) {
                var aid = ring[k];
                var atom = mol.getAtom(aid);
                var nb = mol.getNeighbors(aid);
                // Stable iteration order over neighbour ids (sorted).
                var subIds = [];
                for (var ni = 0; ni < nb.length; ni++) {
                    if (!ringAtomSet[nb[ni]]) subIds.push(nb[ni]);
                }
                if (subIds.length === 0) continue;
                subIds.sort(function (x, y) { return x - y; });

                // Radial direction (centroid -> atom, normalised).
                var rdx = atom.x - rc.x, rdy = atom.y - rc.y;
                var rd = Math.sqrt(rdx * rdx + rdy * rdy);
                if (rd < EPSILON) continue;
                rdx /= rd; rdy /= rd;
                var radialAngle = Math.atan2(rdy, rdx);

                // If single sub: place along radial direction.
                // If multiple: fan symmetrically around radial direction.
                var totalSubs = subIds.length;
                var fan = (totalSubs === 1) ? 0 : DEG120 * 0.5; // +/-60 deg
                for (var si = 0; si < totalSubs; si++) {
                    var subId = subIds[si];
                    var subAtom = mol.getAtom(subId);
                    var offset;
                    if (totalSubs === 1) offset = 0;
                    else if (totalSubs === 2) offset = (si === 0) ? -fan : fan;
                    else {
                        // Spread evenly across [-PI/3, +PI/3]
                        offset = -fan + (2 * fan) * (si / (totalSubs - 1));
                    }
                    var targetAng = radialAngle + offset;

                    // Current bond direction
                    var cdx = subAtom.x - atom.x, cdy = subAtom.y - atom.y;
                    var curLen = Math.sqrt(cdx * cdx + cdy * cdy);
                    if (curLen < EPSILON) curLen = BOND_LENGTH;
                    var curAng = Math.atan2(cdy, cdx);

                    var rot = targetAng - curAng;
                    while (rot > Math.PI) rot -= TWO_PI;
                    while (rot < -Math.PI) rot += TWO_PI;
                    if (Math.abs(rot) < 0.02) continue;

                    // Rotate the substituent subtree (DFS not crossing into ring).
                    rotateSubtree(mol, subId, aid, atom.x, atom.y, rot, ringAtomSet, atomIds);
                }
            }
        }
    }

    /**
     * Rotate the connected subtree rooted at `rootId` (excluding
     * `parentId` and any ring atoms) about pivot (px, py) by rot rad.
     */
    function rotateSubtree(mol, rootId, parentId, px, py, rot, ringAtomSet, atomIds) {
        var visited = {};
        visited[parentId] = true;
        var stack = [rootId];
        var toMove = [];
        while (stack.length) {
            var cur = stack.pop();
            if (visited[cur]) continue;
            visited[cur] = true;
            // Stop at any ring atom (do not pull rings around).
            if (ringAtomSet[cur] && cur !== rootId) continue;
            toMove.push(cur);
            var nb = mol.getNeighbors(cur);
            // Stable: smaller ids first.
            var sortedNb = nb.slice().sort(function (a, b) { return a - b; });
            for (var ni = 0; ni < sortedNb.length; ni++) {
                if (!visited[sortedNb[ni]]) stack.push(sortedNb[ni]);
            }
        }
        var c = Math.cos(rot), s = Math.sin(rot);
        // Sort toMove for determinism
        toMove.sort(function (a, b) { return a - b; });
        for (var i = 0; i < toMove.length; i++) {
            var a = mol.getAtom(toMove[i]);
            if (!a) continue;
            var rx = a.x - px, ry = a.y - py;
            a.x = px + rx * c - ry * s;
            a.y = py + rx * s + ry * c;
        }
    }

    /**
     * Step 4: For every fused ring system, compute the longest axis
     * through ring centroids (PCA-style, but on a list of <=20 points)
     * and rotate the system so this axis is horizontal.
     *
     * Skipped if ring system has only one ring (Step 2 already handles).
     */
    function enforceFusedAxisHorizontal(mol, rings, ringSystems) {
        if (!ringSystems || ringSystems.length === 0) return;
        for (var si = 0; si < ringSystems.length; si++) {
            var sys = ringSystems[si];
            if (sys.length < 2) continue;

            // Collect ring centroids in deterministic order.
            var pts = [];
            for (var ri = 0; ri < sys.length; ri++) {
                var ring = rings[sys[ri]];
                var cx = 0, cy = 0;
                for (var k = 0; k < ring.length; k++) {
                    var a = mol.getAtom(ring[k]); cx += a.x; cy += a.y;
                }
                pts.push({ x: cx / ring.length, y: cy / ring.length });
            }

            // Compute centroid of centroids.
            var mx = 0, my = 0;
            for (var k = 0; k < pts.length; k++) { mx += pts[k].x; my += pts[k].y; }
            mx /= pts.length; my /= pts.length;

            // 2x2 covariance matrix.
            var sxx = 0, syy = 0, sxy = 0;
            for (var k = 0; k < pts.length; k++) {
                var dx = pts[k].x - mx, dy = pts[k].y - my;
                sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
            }

            // Principal axis angle.
            var axisAng = 0.5 * Math.atan2(2 * sxy, sxx - syy);
            // We want this axis horizontal (angle 0). Rotate by -axisAng.
            var rot = -axisAng;
            while (rot > Math.PI / 2) rot -= Math.PI;
            while (rot < -Math.PI / 2) rot += Math.PI;
            if (Math.abs(rot) < 0.05) continue;

            // Collect every atom in the ring system + non-ring neighbours
            // (so attached chains rotate too).
            var sysAtoms = collectSystemAtoms(rings, sys);
            var include = {};
            for (var k = 0; k < sysAtoms.length; k++) include[sysAtoms[k]] = true;
            // BFS outward through non-ring atoms (stop at any other ring system).
            var queue = sysAtoms.slice();
            for (var qh = 0; qh < queue.length; qh++) {
                var cur = queue[qh];
                var nb = mol.getNeighbors(cur);
                for (var ni = 0; ni < nb.length; ni++) {
                    var nbId = nb[ni];
                    if (include[nbId]) continue;
                    // Stop at atoms in any OTHER ring system.
                    var inOtherRing = false;
                    for (var rri = 0; rri < rings.length; rri++) {
                        if (sys.indexOf(rri) >= 0) continue;
                        if (rings[rri].indexOf(nbId) >= 0) { inOtherRing = true; break; }
                    }
                    if (inOtherRing) continue;
                    include[nbId] = true;
                    queue.push(nbId);
                }
            }

            var c = Math.cos(rot), s = Math.sin(rot);
            // Iterate over a deterministic key order (sorted atom ids).
            var ids = Object.keys(include).map(Number).sort(function(a,b){return a-b;});
            for (var k = 0; k < ids.length; k++) {
                var a = mol.getAtom(ids[k]);
                if (!a) continue;
                var rx = a.x - mx, ry = a.y - my;
                a.x = mx + rx * c - ry * s;
                a.y = my + rx * s + ry * c;
            }
        }
    }

    /**
     * Step 5: For non-ring atoms with degree 2-3, smooth bond angles
     * toward the ideal value.
     *   sp/triple-bond chain    : ~180 deg  (linear)
     *   sp2 / sp3 / aromatic    : ~120 deg  (zigzag)
     * If any pair of bonds at a centre is < 90 or > 150, gently rotate
     * the smaller subtree so the angle is pulled to 120.
     */
    function smoothBondAngles(mol, atomIds, ringAtomSet) {
        var n = atomIds.length;
        if (n < 3) return;

        // Iterate atoms in stable id order.
        var sortedIds = atomIds.slice().sort(function(a, b){ return a - b; });
        for (var i = 0; i < sortedIds.length; i++) {
            var aid = sortedIds[i];
            if (ringAtomSet[aid]) continue; // skip ring atoms (geometry locked)
            var atom = mol.getAtom(aid);
            if (!atom) continue;
            var nb = mol.getNeighbors(aid);
            if (!nb || nb.length < 2 || nb.length > 4) continue;

            // Determine target angle based on hybridisation hint:
            //   any incident triple bond -> 180 (linear sp).
            //   else if any double or aromatic flag -> 120 (sp2).
            //   else -> 120 (sp3 zigzag in 2D).
            var hasTriple = false;
            var bondsHere = mol.getBondsOfAtom ? mol.getBondsOfAtom(aid) : [];
            for (var bi = 0; bi < bondsHere.length; bi++) {
                if (bondsHere[bi].type === 3) { hasTriple = true; break; }
            }
            var targetAngle = hasTriple ? Math.PI : DEG120;

            // Examine all unordered neighbour pairs (deterministic via sort).
            var sortedNb = nb.slice().sort(function(a, b){ return a - b; });
            for (var p = 0; p < sortedNb.length; p++) {
                for (var q = p + 1; q < sortedNb.length; q++) {
                    var u = sortedNb[p], v = sortedNb[q];
                    var au = mol.getAtom(u), av = mol.getAtom(v);
                    if (!au || !av) continue;
                    var ux = au.x - atom.x, uy = au.y - atom.y;
                    var vx = av.x - atom.x, vy = av.y - atom.y;
                    var lu = Math.sqrt(ux * ux + uy * uy);
                    var lv = Math.sqrt(vx * vx + vy * vy);
                    if (lu < EPSILON || lv < EPSILON) continue;
                    var cosT = (ux * vx + uy * vy) / (lu * lv);
                    if (cosT > 1) cosT = 1; else if (cosT < -1) cosT = -1;
                    var ang = Math.acos(cosT);
                    // Trigger only if outside a chemistry-reasonable band.
                    var lower = 90 * Math.PI / 180;
                    var upper = 150 * Math.PI / 180;
                    if (hasTriple) { lower = 150 * Math.PI / 180; upper = Math.PI; }
                    if (ang >= lower && ang <= upper) continue;

                    // Rotate the smaller-id subtree so the angle becomes target.
                    // Compute current signed angle of u and v about atom.
                    var angU = Math.atan2(uy, ux);
                    var angV = Math.atan2(vy, vx);
                    // Rotate v away from u by (targetAngle - currentAngle) on
                    // the side that does not pass through u\u2019s direction.
                    var diff = angV - angU;
                    while (diff > Math.PI) diff -= TWO_PI;
                    while (diff < -Math.PI) diff += TWO_PI;
                    var sign = (diff >= 0) ? 1 : -1;
                    var desired = sign * targetAngle;
                    var rot = desired - diff;
                    // Halve the correction so neither subtree dominates.
                    rot *= 0.5;
                    if (Math.abs(rot) < 0.02) continue;

                    // Rotate the v subtree (and possibly the u subtree by -rot/2 each)
                    // For determinism, always rotate the larger-id subtree.
                    var subRoot = (u > v) ? u : v;
                    var subSign = (subRoot === v) ? 1 : -1;
                    rotateSubtree(mol, subRoot, aid, atom.x, atom.y, rot * subSign, ringAtomSet, atomIds);
                }
            }
        }
    }

    // =====================================================================
    // Utility functions
    // =====================================================================

    function dist(x1, y1, x2, y2) {
        var dx = x2 - x1, dy = y2 - y1;
        return Math.sqrt(dx * dx + dy * dy);
    }

    function distSq(x1, y1, x2, y2) {
        var dx = x2 - x1, dy = y2 - y1;
        return dx * dx + dy * dy;
    }

    function normalizeAngle(a) {
        // FIX: guard against NaN/Infinity which would cause infinite loop
        if (!isFinite(a)) return 0;
        while (a < 0)       a += TWO_PI;
        while (a >= TWO_PI) a -= TWO_PI;
        return a;
    }

    function atomBounds(atoms) {
        var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        for (var i = 0; i < atoms.length; i++) {
            if (atoms[i].x < minX) minX = atoms[i].x;
            if (atoms[i].y < minY) minY = atoms[i].y;
            if (atoms[i].x > maxX) maxX = atoms[i].x;
            if (atoms[i].y > maxY) maxY = atoms[i].y;
        }
        return { minX: minX, minY: minY, maxX: maxX, maxY: maxY };
    }

    function uniqueArray(arr) {
        var seen = {}, result = [];
        for (var i = 0; i < arr.length; i++) {
            if (!seen[arr[i]]) { seen[arr[i]] = true; result.push(arr[i]); }
        }
        return result;
    }

    // NOTE: newIntArray was dead code (never called) — removed to reduce bundle size

    // =====================================================================
    // Ring Query Helpers  (used by Renderer for aromatic circles etc.)
    // =====================================================================

    /**
     * Returns ring info for all rings in the molecule.
     * Each entry: { atoms: [atomIds], center: {x,y}, size, aromatic }
     *
     * Aromatic detection uses:
     *   1. atom.aromatic flags from the SMILES parser (preferred)
     *   2. Heuristic: 5- or 6-membered ring with at least one double bond
     */
    Layout.getRingInfo = function (mol) {
        var components = mol.getComponents();
        var allRings   = [];

        for (var ci = 0; ci < components.length; ci++) {
            var rings = perceiveSSSR(mol, components[ci]);

            for (var ri = 0; ri < rings.length; ri++) {
                var ring = rings[ri];
                var cx = 0, cy = 0;
                var allAroFlag   = true;   // all atoms flagged aromatic?
                var hasDouble    = false;

                for (var ai = 0; ai < ring.length; ai++) {
                    var atom = mol.getAtom(ring[ai]);
                    cx += atom.x;
                    cy += atom.y;
                    if (!atom.aromatic) allAroFlag = false;
                }
                cx /= ring.length;
                cy /= ring.length;

                // Check bond types around the ring
                for (var ai = 0; ai < ring.length; ai++) {
                    var bond = mol.getBondBetween(ring[ai], ring[(ai + 1) % ring.length]);
                    if (bond && bond.type === Molecule.BOND_DOUBLE) hasDouble = true;
                }

                // Aromatic if SMILES flagged it, or heuristic 5/6-ring + double bond
                var isAromatic = allAroFlag ||
                    (hasDouble && (ring.length === 5 || ring.length === 6));

                allRings.push({
                    atoms:    ring,
                    center:   { x: cx, y: cy },
                    size:     ring.length,
                    aromatic: isAromatic
                });
            }
        }

        return allRings;
    };

    /**
     * Check if a bond is part of a ring.  Returns the ring info or null.
     */
    Layout.bondInRing = function (mol, bond) {
        var ringInfo = Layout.getRingInfo(mol);
        for (var i = 0; i < ringInfo.length; i++) {
            var ring = ringInfo[i].atoms;
            var idx1 = ring.indexOf(bond.atom1);
            var idx2 = ring.indexOf(bond.atom2);
            if (idx1 >= 0 && idx2 >= 0) {
                var diff = Math.abs(idx1 - idx2);
                if (diff === 1 || diff === ring.length - 1) return ringInfo[i];
            }
        }
        return null;
    };

    // =====================================================================
    // Export
    // =====================================================================

    // =====================================================================
    // Export
    // =====================================================================

    global.Layout = Layout;

})(typeof window !== 'undefined' ? window : this);
