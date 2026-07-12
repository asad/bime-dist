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
 *   1. SSSR ring perception (edge-vector Gaussian elimination)
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
     *     Off by default since the model under-performs on the built-in
     *     template set; helpful for less-canonical input drawings.
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
        mlDepictWeight: 0.15,
        // Source-tag conditioning for the ML refiner. The default
        // channel is the conservative-average prediction; advanced
        // callers may override if they know the input drawing style.
        mlDepictSourceTag: null,
        // v1.8.12: SDG corrective pass (editor/SDGLayout.js).
        // Three native phases:
        //   1. resolveOverlap — summed-force iterative overlap resolution
        //   2. distributePartners
        //      (substituent fan in the open angular arc)
        //   3. populatePolygonCorners
        //
        // **OPT-IN in v1.8.12.** Default off because the diagnostic
        // shows it net-regresses on the regression set:
        //
        //                       useSDG=false   useSDG=true
        //   clean rate          97.4%          93.4%   (−4.0 pp)
        //   abnormal bond-len   10             1       (−9, good)
        //   atom overlap        0              50      (+50, bad)
        //   ring-sys overlap    3              8       (+5, bad)
        //
        // The spring force in the overlap resolver is too strong;
        // straightening stretched bonds knocks atoms into neighbours.
        // Tracked for tuning in upcoming v1.8.x patches. Until then,
        // SDG is available as `Layout.options.useSDG = true` for
        // experiment / development.
        useSDG: false,
        // v1.8.17: post-Step-15 SDG refinement passes. Each is
        // independently gated; defaults reflect the diagnostic verdict
        // for each phase.
        //   useSDGRefine          — master switch (Step 16 entirely; ON)
        //   correctEZ             — flip mis-drawn E/Z double bonds
        //                           OFF by default — affects SMILES
        //                           round-trip directional markers; the
        //                           cipLabel-based correction is opt-in
        //                           pending further round-trip-stable
        //                           tuning.
        //   assignWedgeDash       — pick wedge/dash for sp³ stereo
        //                           OFF by default — sets bond.stereo on
        //                           previously-unmarked bonds, which
        //                           SmilesWriter encodes as / \ markers.
        //                           Opt in when the caller wants 2D
        //                           wedge/dash for visual depiction
        //                           (e.g., MOL export, SVG render).
        //   alignToLongestAxis    — PCA rotation to horizontal (ON;
        //                           pure rigid transformation, preserves
        //                           bond lengths and ring shapes)
        //   rotateRings           — try 180° per-ring, keep best (ON;
        //                           only flips standalone rings, never
        //                           fused / spiro)
        useSDGRefine: true,
        // v1.8.19: correctEZ now uses full CIP priorities from
        // editor/CIPStereo.js (depth-N digraph traversal). Safe-by-
        // default — only acts on bonds with explicit cipLabel = 'E' / 'Z'
        // (set by CIPStereo.assignEZ); inert for molecules without
        // E/Z stereo annotations.
        correctEZ: true,
        // v1.8.19: assignWedgeDash safe-by-default. Now writes to
        // bond.depictStereo + bond.depictStereoFromAtom (depiction-only
        // fields read by ImageExport) instead of bond.stereo, so SMILES
        // round-trip is preserved.
        assignWedgeDash: true,
        alignToLongestAxis: true,
        rotateRings: true,
        // v1.8.18: explicit-H placement. ON by default — it only acts
        // on hydrogen atoms whose coordinates haven't been placed yet
        // (rare in SMILES-parsed molecules where H is usually implicit;
        // common in MOL-imported molecules with explicit H records).
        placeExplicitH: true,
        // v1.8.19: macro-chain linearisation. For molecules with ≥ 2
        // ring systems connected by chains, finds the longest
        // topological path through the ring-system graph and rotates
        // the molecule so that path lies along the x-axis. Closes the
        // bisazo dye TARGET test and improves any "rod-shaped"
        // chromophore depiction.
        linearizeMacroChain: true,
        // Run the SDG corrective pass only when cheap quality
        // checks detect a severe depiction defect (badly stretched bonds
        // or atom clashes). This keeps the conservative default path for
        // ordinary structures while rescuing hard multi-ring/branched
        // molecules without requiring callers to opt into useSDG globally.
        adaptiveSDGRescue: true
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
     * Public quality inspection hook used by tests and downstream callers.
     * Returns SDG.LayoutQuality metrics when that module is loaded, otherwise
     * falls back to the local lightweight guard used by adaptive rescue.
     */
    Layout.quality = function (mol, atomIds, options) {
        atomIds = atomIds || (mol && mol.atoms ? mol.atoms.map(function (a) { return a.id; }) : []);
        var quality = global.SDG && global.SDG.LayoutQuality;
        if (quality && quality.evaluate) {
            options = options || {};
            if (!options.bondLength) options.bondLength = BOND_LENGTH;
            return quality.evaluate(mol, atomIds, options);
        }
        return layoutQualityMetrics(mol, atomIds);
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
        // and was removed alongside the misleading comment. SMSDRings (the
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

        // -- Step 10: Detect crossings AND ring-system overlap.
        // v1.8.10: previously the heavy SMSDLayout passes only fired on
        // bond crossings. But complex multi-ring structures (e.g. Congo
        // red, two naphthalenes connected by azo bridges) can have ring
        // systems piled on top of each other WITHOUT bond crossings —
        // the connecting bonds end at common atoms, so the line-segment
        // intersection test misses them. We now also fire the heavy
        // passes when two distinct (non-fused) ring systems sit closer
        // than 1.5 * BL apart — a textbook structural-diagram-generator
        // condition for "the rings need force-directed separation".
        var hasCrossings = false;
        var hasRingOverlap = false;
        var smsdLayout10 = global.SMSDLayout;
        if (smsdLayout10) {
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

            // v1.8.10: ring-system overlap detector. Compute centroids of
            // each ring; if two rings are < 1.5*BL apart but share fewer
            // than 2 atoms (i.e. they're not fused), the layout has piled
            // them together unintentionally.
            if (rings && rings.length >= 3) {
                var rocCentroids = [];
                for (var rci = 0; rci < rings.length; rci++) {
                    var rcAtoms = rings[rci].atoms || rings[rci];
                    var rcCx = 0, rcCy = 0;
                    var rcSet = {};
                    for (var rcj = 0; rcj < rcAtoms.length; rcj++) {
                        var rcA = mol.getAtom(rcAtoms[rcj]);
                        if (rcA) { rcCx += rcA.x; rcCy += rcA.y; rcSet[rcAtoms[rcj]] = true; }
                    }
                    rcCx /= rcAtoms.length; rcCy /= rcAtoms.length;
                    rocCentroids.push({ cx: rcCx, cy: rcCy, set: rcSet });
                }
                for (var rA = 0; rA < rocCentroids.length && !hasRingOverlap; rA++) {
                    for (var rB = rA + 1; rB < rocCentroids.length; rB++) {
                        var rdx = rocCentroids[rA].cx - rocCentroids[rB].cx;
                        var rdy = rocCentroids[rA].cy - rocCentroids[rB].cy;
                        var rDist = Math.sqrt(rdx * rdx + rdy * rdy);
                        if (rDist >= BOND_LENGTH * 1.5) continue;
                        // Count shared atoms.
                        var rShared = 0;
                        var sA = rocCentroids[rA].set, sB = rocCentroids[rB].set;
                        for (var sk in sA) if (sB[sk]) rShared++;
                        if (rShared < 2) { hasRingOverlap = true; break; }
                    }
                }
            }
        }

        var needsForceDirected = hasCrossings || hasRingOverlap;
        var hasSMSDPass10 = !!(smsdLayout10 &&
            (smsdLayout10.reduceCrossings || smsdLayout10.forceDirectedLayout));
        var forceSnap10 = null;
        var forceScore10 = 0;
        var forceNeededRescue10 = false;
        var forceMaxBond10 = 0;
        if (needsForceDirected && hasSMSDPass10) {
            var forceMetrics10 = layoutQualityMetrics(mol, atomIds);
            forceSnap10 = snapshotAtomPositions(mol, atomIds);
            forceScore10 = forceMetrics10.penalty;
            forceNeededRescue10 = forceMetrics10.needsAdaptiveRescue;
            forceMaxBond10 = maxBondLengthForAtomSet(mol, atomIds);
        }

        // -- Step 10b: Crossing reduction (when crossings or ring overlap) --
        if (needsForceDirected && smsdLayout10 && smsdLayout10.reduceCrossings) {
            smsdLayout10.reduceCrossings(mol, atomIds, ringAtomSet, rings);
        }

        // -- Step 11: Force-directed refinement (when crossings or overlap) --
        if (needsForceDirected && smsdLayout10 && smsdLayout10.forceDirectedLayout) {
            smsdLayout10.forceDirectedLayout(mol, atomIds, ringAtomSet, rings);
        }
        if (needsForceDirected && forceSnap10) {
            var forceMetricsAfter10 = layoutQualityMetrics(mol, atomIds);
            var forceScoreAfter10 = forceMetricsAfter10.penalty;
            var forceNeedsRescueAfter10 = forceMetricsAfter10.needsAdaptiveRescue;
            var forceMaxBondAfter10 = maxBondLengthForAtomSet(mol, atomIds);
            if ((!forceNeededRescue10 && forceNeedsRescueAfter10) ||
                (forceMaxBondAfter10 > BOND_LENGTH * 1.2 &&
                 forceMaxBondAfter10 > forceMaxBond10 + 0.5) ||
                forceScoreAfter10 > forceScore10 + 2) {
                restoreAtomPositions(mol, forceSnap10);
            }
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
            var srcTag = Layout.options.mlDepictSourceTag;
            try { MLDepict.refineLayout(subMol, alpha, srcTag); }
            catch (mle) { /* swallow — never break layout because of ML */ }
            // After ML refinement, run a single bond-length normalisation
            // pass to avoid the blend pulling bonds away from BOND_LENGTH.
            relaxBondLengths(mol, atomIds, ringAtomSet);
        }

        // -- Step 13 (v1.8.9): collapsed-ring + collapsed-chain rescue ----
        // The chain-extension + collision-resolution pipeline sometimes
        // (~1-2% of complex polycyclic / multi-fragment SMILES) leaves
        // either:
        //   (a) a ring whose 6 atoms are folded into a sub-BL cluster
        //       — renders as a stack of aromatic-ring indicator circles
        //   (b) two non-bonded atoms < 0.4*BL apart in a chain region
        //       — renders as labels on top of each other
        //
        // We address both with two iterative passes:
        //   1. Re-inflate any ring whose perimeter bonds are sub-BL/2
        //      (regular polygon around the existing centroid).
        //   2. Pairwise repulsion + bond-length relaxation, repeated
        //      until no two non-bonded atoms remain within 0.5*BL of
        //      each other (capped at 6 iterations so a pathologically
        //      tangled large polycyclic doesn't loop forever).
        if (rings && rings.length > 0) {
            var anyReinflated = false;
            for (var rri = 0; rri < rings.length; rri++) {
                var rRing = rings[rri];
                var rAtomIds = rRing.atoms || rRing;
                if (!rAtomIds || rAtomIds.length < 3) continue;
                var maxIntraBL = 0;
                for (var ri2 = 0; ri2 < rAtomIds.length; ri2++) {
                    var ra = mol.getAtom(rAtomIds[ri2]);
                    var rb = mol.getAtom(rAtomIds[(ri2 + 1) % rAtomIds.length]);
                    if (!ra || !rb) continue;
                    var rdx = rb.x - ra.x, rdy = rb.y - ra.y;
                    var rd = Math.sqrt(rdx * rdx + rdy * rdy);
                    if (rd > maxIntraBL) maxIntraBL = rd;
                }
                if (maxIntraBL > 0 && maxIntraBL < BOND_LENGTH * 0.5) {
                    var rcx = 0, rcy = 0;
                    for (var ri3 = 0; ri3 < rAtomIds.length; ri3++) {
                        var rca = mol.getAtom(rAtomIds[ri3]);
                        if (rca) { rcx += rca.x; rcy += rca.y; }
                    }
                    rcx /= rAtomIds.length; rcy /= rAtomIds.length;
                    var nR = rAtomIds.length;
                    var radius = BOND_LENGTH / (2 * Math.sin(Math.PI / nR));
                    var first = mol.getAtom(rAtomIds[0]);
                    var startAng = (first ? Math.atan2(first.y - rcy, first.x - rcx) : 0);
                    if (!isFinite(startAng)) startAng = 0;
                    for (var ri4 = 0; ri4 < nR; ri4++) {
                        var ra4 = mol.getAtom(rAtomIds[ri4]);
                        if (!ra4) continue;
                        var ang = startAng + (2 * Math.PI * ri4) / nR;
                        ra4.x = rcx + radius * Math.cos(ang);
                        ra4.y = rcy + radius * Math.sin(ang);
                    }
                    anyReinflated = true;
                }
            }
            if (anyReinflated) {
                relaxBondLengths(mol, atomIds, ringAtomSet);
            }
        }

        // (b) Chain-collision repair via summed-force pass.
        // Pairwise sequential pushes were tried but created secondary
        // overlaps (atoms pushed away from one neighbour land on top of
        // another). The correct approach is the existing v1.5.x
        // resolveCollisions which uses summed forces. We just run it
        // once more after the ring re-inflation so any chains that were
        // disrupted get a fresh sweep.
        resolveCollisions(mol, atomIds, ringAtomSet);

        // -- Step 14 (v1.8.12): SDG corrective pass ----------------------
        // editor/SDGLayout.js implements three native diagram-generation
        // phases:
        //   1. resolveOverlap — summed-force iterative overlap resolution
        //   2. distributePartners — substituent fan in the open angular arc
        //      opposite placed neighbours
        //   3. final overlap clean-up
        // Active by default in v1.8.12 (Layout.options.useSDG = true).
        // Wrapped in try/catch so a SDGLayout failure never breaks an
        // otherwise-valid layout.
        var sdgLayout14 = global.SDGLayout;
        if (Layout.options && Layout.options.useSDG &&
            sdgLayout14 && sdgLayout14.refineLayout) {
            runSDGRefinement(mol, atomIds, ringAtomSet, sdgLayout14, false);
        } else if (Layout.options &&
            Layout.options.adaptiveSDGRescue !== false &&
            sdgLayout14 && sdgLayout14.refineLayout &&
            layoutNeedsAdaptiveRescue(mol, atomIds)) {
            runSDGRefinement(mol, atomIds, ringAtomSet, sdgLayout14, true);
        }

        // -- Step 15 (v1.8.16): FUSED-RING-AWARE polygon snap + chain relax
        // Absolute last step. Forces every ring's atoms to be at exact
        // regular-polygon positions. Deliberately skips the post-relax
        // that would pull atoms back out of polygon shape.
        //
        // v1.8.15 KEY CHANGE — Process ring SYSTEMS, not rings:
        //   - For each fused ring system, place a SEED ring as a regular
        //     polygon at its current centroid + first-atom orientation,
        //     then BFS-place every other ring in the system by EDGE
        //     REFLECTION across the already-snapped shared atoms.
        //   - This preserves regular polygon geometry on BOTH sides of
        //     a fusion bond — fixing the v1.8.14 bug where the second
        //     ring's snap stomped atoms shared with the first ring.
        //   - 2 shared contiguous atoms → edge-reflection.
        //   - 1 shared atom → spiro rotation.
        //   - 3+ shared atoms → bridged; fall back to centroid snap.
        //
        // Diagnostic vs v1.8.14 / v1.8.15 on the regression set
        // (tools/diagnose-ring-quality.js):
        //                                  v1.8.14   v1.8.15  v1.8.16
        //   overall ring clean-rate:        68.4 %    79.4 %    80.9 %
        //   per-mol mean clean-rate:           --     93.3 %    94.8 %
        //   per-mol p10 clean-rate:            --     71.4 %   100.0 %
        //   6-ring (benzene, etc.):            --     93.2 %    94.8 %
        //   5-ring:                            --     86.6 %    88.5 %
        //   chain bond-length CV (p50):        --     31.3 %    20.5 %
        //   atom overlaps:                       0         0         0
        //
        // v1.8.16 wins came from three additions to Step 15:
        //   1. relaxChainBonds — chain-only iterative bond-length pull
        //      (ring atoms FROZEN; only chain atoms move). Drops chain
        //      CV from 31% to 20% at p50.
        //   2. separateOverlappingRingSystems — RIGID translation of
        //      colliding ring systems by worst-pair deficit along the
        //      centroid axis. Polygon shape preserved exactly.
        //   3. resolveCollisions in chain mode — collision resolver
        //      that NEVER pushes atoms in different ring systems (those
        //      are handled by 2). Intra-system gentle push for genuinely
        //      bridged geometry (morphine, norbornane).
        // Bisazo dye benchmark: all four 6-rings remain perfect hexagons.
        //
        // Idempotent: a ring system already at regular-polygon positions
        // produces no positional change.
        //
        // Ring placement follows BIME's native SDG geometry rules.
        //   methods: placeRing, placeFusedRing, placeSpiroRing.
        if (rings && rings.length > 0 && ringSystems && ringSystems.length > 0) {
            var ringIdSet15 = {};
            for (var rsi15 = 0; rsi15 < rings.length; rsi15++) {
                var rsAtoms15 = rings[rsi15].atoms || rings[rsi15];
                for (var rsj15 = 0; rsj15 < rsAtoms15.length; rsj15++) {
                    ringIdSet15[rsAtoms15[rsj15]] = true;
                }
            }
            // Snapshot pre-snap positions of every ring atom so we can
            // compute per-atom displacement and propagate to non-ring
            // neighbours (substituents) at the end.
            var oldPos15 = {};
            for (var snapId in ringIdSet15) {
                var snapA = mol.getAtom(parseInt(snapId, 10));
                if (snapA) oldPos15[snapId] = { x: snapA.x, y: snapA.y };
            }

            var ringApplied = false;
            for (var sysIdx15 = 0; sysIdx15 < ringSystems.length; sysIdx15++) {
                var sys15 = ringSystems[sysIdx15];
                if (!sys15 || sys15.length === 0) continue;

                // Bornane-style 2x5 bridged cages already encode the correct
                // projection. Snapping those systems back to separate
                // regular polygons flattens the bridge and makes previews
                // look broken, so preserve that complete template hit.
                var sysAtomsForSnap15 = collectSystemAtoms(rings, sys15);
                var sysTemplated15 = sysAtomsForSnap15.length > 0;
                for (var sta15 = 0; sta15 < sysAtomsForSnap15.length; sta15++) {
                    if (!templatePlaced[sysAtomsForSnap15[sta15]]) {
                        sysTemplated15 = false;
                        break;
                    }
                }
                if (sysTemplated15 &&
                    sys15.length === 2 &&
                    rings[sys15[0]].length === 5 &&
                    rings[sys15[1]].length === 5 &&
                    sysAtomsForSnap15.length === 7) {
                    continue;
                }

                // ----- Single ring system: snap directly. ------------
                if (sys15.length === 1) {
                    if (snapRingAtCentroid(mol, rings[sys15[0]], BOND_LENGTH)) {
                        ringApplied = true;
                    }
                    continue;
                }

                // ----- Multi-ring system: seed + reflection. ---------
                // Choose the largest ring as seed (best scaffold for
                // polycyclics — matches layoutRingSystem's heuristic).
                var seedIdx15 = sys15[0];
                for (var ssi15 = 1; ssi15 < sys15.length; ssi15++) {
                    if (rings[sys15[ssi15]].length > rings[seedIdx15].length) {
                        seedIdx15 = sys15[ssi15];
                    }
                }
                var snappedAtomSet15 = {};
                var snappedRingSet15 = {};

                if (snapRingAtCentroid(mol, rings[seedIdx15], BOND_LENGTH)) {
                    ringApplied = true;
                }
                for (var sai15 = 0; sai15 < rings[seedIdx15].length; sai15++) {
                    snappedAtomSet15[rings[seedIdx15][sai15]] = true;
                }
                snappedRingSet15[seedIdx15] = true;

                // BFS: place each remaining ring by reflection from an
                // already-snapped neighbour. Loops until no more rings
                // can be placed (same fixed-point pattern as
                // layoutRingSystem above).
                var madeProgress15 = true;
                while (madeProgress15) {
                    madeProgress15 = false;
                    // Score-based: prefer simple edge-fusion (2 contig.)
                    // first, then bridged (3+), then spiro (1) — same
                    // ordering as layoutRingSystem.
                    var bestRIdx = -1, bestScore = -1, bestShared = null;
                    for (var ssi16 = 0; ssi16 < sys15.length; ssi16++) {
                        var rIdx16 = sys15[ssi16];
                        if (snappedRingSet15[rIdx16]) continue;
                        var ring16 = rings[rIdx16];
                        var shared16 = [];
                        for (var ai16 = 0; ai16 < ring16.length; ai16++) {
                            if (snappedAtomSet15[ring16[ai16]]) shared16.push(ring16[ai16]);
                        }
                        if (shared16.length === 0) continue;
                        // Contiguity scan.
                        var contigLen = 0, maxContig = 0;
                        for (var ai17 = 0; ai17 < ring16.length * 2; ai17++) {
                            var idx17 = ai17 % ring16.length;
                            if (snappedAtomSet15[ring16[idx17]]) {
                                contigLen++;
                                if (contigLen > maxContig) maxContig = contigLen;
                            } else {
                                contigLen = 0;
                            }
                        }
                        if (maxContig > shared16.length) maxContig = shared16.length;
                        var score16;
                        if (shared16.length === 2 && maxContig === 2) {
                            score16 = 1000 + ring16.length;       // edge-fused
                        } else if (shared16.length === 1) {
                            score16 = 100 + ring16.length;         // spiro
                        } else if (maxContig >= 3) {
                            score16 = 500 + shared16.length * 10;  // bridged
                        } else {
                            score16 = 300 + shared16.length * 10;  // mixed
                        }
                        if (score16 > bestScore) {
                            bestScore = score16;
                            bestRIdx = rIdx16;
                            bestShared = shared16;
                        }
                    }
                    if (bestRIdx < 0) break;
                    var bestRing = rings[bestRIdx];
                    if (bestShared.length === 1) {
                        if (snapSpiroRingFromAtom(mol, bestRing, bestShared[0],
                                                     snappedAtomSet15, BOND_LENGTH)) {
                            ringApplied = true;
                        }
                    } else if (bestShared.length === 2) {
                        if (snapFusedRingByReflection(mol, bestRing, bestShared,
                                                         snappedAtomSet15, BOND_LENGTH)) {
                            ringApplied = true;
                        }
                    } else {
                        // Bridged (3+ shared) — keep current positions.
                        // Full bridged-ring placement lives in the SDG ring placer
                        // .placeBridgedRing; tracked for a follow-up patch.
                    }
                    for (var sai17 = 0; sai17 < bestRing.length; sai17++) {
                        snappedAtomSet15[bestRing[sai17]] = true;
                    }
                    snappedRingSet15[bestRIdx] = true;
                    madeProgress15 = true;
                }
                // Any rings still unsnapped (rare: disconnected within a
                // system, or all-bridged) → fall back to centroid snap.
                for (var ssi18 = 0; ssi18 < sys15.length; ssi18++) {
                    if (snappedRingSet15[sys15[ssi18]]) continue;
                    if (snapRingAtCentroid(mol, rings[sys15[ssi18]], BOND_LENGTH)) {
                        ringApplied = true;
                    }
                }
            }

            if (ringApplied) {
                // Compute substituent displacements from the pre/post
                // snapshots and apply averaged shifts to non-ring
                // neighbours (so substituents follow their ring atom).
                var subDisp15 = {};   // id → {dx, dy, n}
                for (var ridSnap in oldPos15) {
                    var newA = mol.getAtom(parseInt(ridSnap, 10));
                    if (!newA) continue;
                    var dX15 = newA.x - oldPos15[ridSnap].x;
                    var dY15 = newA.y - oldPos15[ridSnap].y;
                    if (Math.abs(dX15) < 1e-9 && Math.abs(dY15) < 1e-9) continue;
                    var nbrs15 = mol.getNeighbors(parseInt(ridSnap, 10)) || [];
                    for (var nbi15 = 0; nbi15 < nbrs15.length; nbi15++) {
                        if (ringIdSet15[nbrs15[nbi15]]) continue;
                        var nbid15 = nbrs15[nbi15];
                        if (!subDisp15[nbid15]) subDisp15[nbid15] = { dx: 0, dy: 0, n: 0 };
                        subDisp15[nbid15].dx += dX15;
                        subDisp15[nbid15].dy += dY15;
                        subDisp15[nbid15].n  += 1;
                    }
                }
                Object.keys(subDisp15).forEach(function (idStr) {
                    var rec = subDisp15[idStr];
                    if (rec.n === 0) return;
                    var sa = mol.getAtom(parseInt(idStr, 10));
                    if (!sa) return;
                    sa.x += rec.dx / rec.n;
                    sa.y += rec.dy / rec.n;
                });
                // v1.8.16: chain-only bond-length relax. The diagnostic
                // measured chain (non-ring) bond-length CV at p50 31 % in
                // v1.8.15 — substituents at distance D from a ring atom
                // get translated by the ring's snap delta, but 2nd-level
                // chain atoms (with no ring neighbour) stay put, so the
                // 1st→2nd-level bond changes by the magnitude of the snap
                // delta. Fix: pull every chain bond toward BOND_LENGTH
                // while ring atoms stay frozen. This preserves the
                // perfect polygon snap from above and fixes chain CV.
                relaxChainBonds(mol, atomIds, ringAtomSet);
                // v1.8.16: chain-only collision sweep. Ring atoms are
                // FROZEN — they were just snapped to perfect polygon
                // positions and any push would re-distort them. Only
                // chain atoms move to resolve overlaps. Ring-ring
                // overlaps (rare; only happens when two ring systems'
                // centroids drifted on top of each other during earlier
                // steps) are handled by ring-system separation BEFORE
                // this sweep.
                separateOverlappingRingSystems(mol, atomIds, ringSystems, rings, ringAtomSet);
                // Build atom → ring-system-index map so the chain-mode
                // resolveCollisions can apply gentle ring-ring push only
                // INTRA-system (e.g. morphine bridged 5-ring) and never
                // INTER-system (those are already cleared by the rigid
                // separator).
                var atomToSys15 = {};
                for (var sii = 0; sii < ringSystems.length; sii++) {
                    var sysL = ringSystems[sii];
                    for (var sj = 0; sj < sysL.length; sj++) {
                        var ringIds2 = rings[sysL[sj]];
                        if (!ringIds2) continue;
                        for (var sk = 0; sk < ringIds2.length; sk++) {
                            atomToSys15[ringIds2[sk]] = sii;
                        }
                    }
                }
                // v1.8.24: connection-aware pull-back. The separator above
                // clears ring-system OVERLAPS but has no notion of the
                // bond(s) that CONNECT two systems — so a biaryl bond
                // (biphenyl: ring-atom to ring-atom) or a one-atom linker
                // (diphenylmethane: ring-CH2-ring) is left stretched to
                // 1.9-2.5x BOND_LENGTH. This pass rigidly pulls each
                // connected ring system back along its linker axis until
                // the linker bonds are BOND_LENGTH.
                //
                // Snapshot gate: keep the pull-back only when it does not
                // increase atom overlap count. The chain-length repair itself
                // is handled by pullConnectedRingSystems + relaxChainBonds.
                var preSnap24 = snapshotAtomPositions(mol, atomIds);
                var preOverlaps24 = countAtomOverlaps(mol, atomIds);
                pullConnectedRingSystems(mol, atomIds, ringSystems, rings,
                                            ringAtomSet, atomToSys15);
                // Re-relax chain bonds AFTER the pull-back: the pull
                // rigidly translates whole ring systems, which can leave
                // a linker chain atom (e.g. diphenylmethane's CH2)
                // stranded between its now-repositioned ring neighbours.
                relaxChainBonds(mol, atomIds, ringAtomSet);
                // v1.8.25: chain-only collision sweep — the unified
                // resolveCollisions in "chain mode": gentle 0.15x push for
                // intra-system ring pairs, skip inter-system ring pairs
                // (atomToSys15), and DO NOT snapshot/restore rings (they
                // were just snapped to perfect polygons).
                resolveCollisions(mol, atomIds, ringAtomSet, {
                    ringGentleFactor: 0.15,
                    atomToSys: atomToSys15,
                    snapshotRings: false
                });
                var postOverlaps24 = countAtomOverlaps(mol, atomIds);
                if (postOverlaps24 > preOverlaps24) {
                    // Net regression — revert the entire pull sequence.
                    restoreAtomPositions(mol, preSnap24);
                }
            }
        }

        // -- Step 15b: adaptive post-snap rescue ------------------------
        // Ring snapping and rigid ring-system separation can expose an
        // over-constrained branched centre (for example diaryl carbinols):
        // the rings are perfect, but one first-shell bond is left far too
        // long. Re-run the SDG corrective pass only when the final
        // quality gate sees a severe defect, and keep it only if the score
        // improves.
        var sdgLayout15b = global.SDGLayout;
        if (Layout.options && Layout.options.adaptiveSDGRescue !== false &&
            sdgLayout15b && sdgLayout15b.refineLayout &&
            layoutNeedsAdaptiveRescue(mol, atomIds)) {
            runSDGRefinement(mol, atomIds, ringAtomSet, sdgLayout15b, true);
        }

        // -- Step 16 (v1.8.17): SDG post-refinement passes ----------------
        //   16a. CorrectGeometricConfiguration — flip mis-drawn E/Z double
        //        bonds (CIP perception sets bond.cipLabel; this checks
        //        whether the depicted geometry matches and reflects the
        //        smaller subtree if not).
        //   16b. NonplanarBonds.assign — pick the best wedge/dash bond at
        //        each tetrahedral stereo centre.
        //   16c. LayoutRefiner.alignToLongestAxis — PCA rotation so the
        //        principal axis is horizontal (publication convention).
        //   16d. LayoutRefiner.rotateRings — try 180° flip per ring, keep
        //        the lower-Congestion orientation.
        //
        // Each phase is independently gated by Layout.options so callers
        // can disable any step. Defaults: all enabled.
        if (Layout.options && Layout.options.useSDGRefine !== false) {
            var sdg16 = global.SDG;
            try {
                if (sdg16 && sdg16.CorrectGeometricConfiguration &&
                    Layout.options.correctEZ !== false) {
                    // Populate bond.cipLabel from the geometry-free directional
                    // markers (SMILES / \) BEFORE correcting geometry. Layout
                    // rebuilds every C=C as a trans zigzag and the coordinate-
                    // based assignEZ is circular, so without this cis bonds
                    // (maleate, oleate, cis-2-butene) are silently drawn trans.
                    // Only bonds that actually carry directional markers get a
                    // label; unmarked / MOL-imported bonds are left untouched.
                    if (global.CIPStereo && global.CIPStereo.doubleBondEZ) {
                        for (var _dbi = 0; _dbi < mol.bonds.length; _dbi++) {
                            var _db = mol.bonds[_dbi];
                            if (_db.type !== 2) { continue; }
                            var _ez = global.CIPStereo.doubleBondEZ(mol, _db);
                            if (_ez === 'E' || _ez === 'Z') { _db.cipLabel = _ez; }
                        }
                    }
                    sdg16.CorrectGeometricConfiguration.correct(mol);
                }
            } catch (e1) { /* keep layout if EZ correction fails */ }
            try {
                if (sdg16 && sdg16.NonplanarBonds &&
                    Layout.options.assignWedgeDash !== false) {
                    sdg16.NonplanarBonds.assign(mol);
                }
            } catch (e2) { /* keep layout if wedge/dash assign fails */ }
            try {
                if (sdg16 && sdg16.LayoutRefiner) {
                    sdg16.LayoutRefiner.refine(mol, rings, {
                        alignToLongestAxis: Layout.options.alignToLongestAxis !== false,
                        rotateRings: Layout.options.rotateRings !== false,
                        // Skip the SDG resolveOverlaps — Step 15 already
                        // ran chain-only overlap clearing and ring-system
                        // separation; an additional pass would just
                        // perturb the polished result.
                        resolveOverlaps: false
                    });
                }
            } catch (e3) { /* keep layout if refine fails */ }
            // 16e (v1.8.18). HydrogenPlacer — place explicit H atoms.
            // Operates only on H atoms whose coordinates are still at
            // origin (i.e., never placed by the earlier ring/chain
            // logic). Implicit H counts on heavy atoms are NOT touched
            // — they're rendered by the SVG depicter, not the layout.
            try {
                if (sdg16 && sdg16.HydrogenPlacer &&
                    Layout.options.placeExplicitH !== false) {
                    sdg16.HydrogenPlacer.placeHydrogens2D(mol, BOND_LENGTH);
                }
            } catch (e4) { /* keep layout if H placement fails */ }
            // 16f (v1.8.19). Macro-chain linearisation. For molecules
            // where multiple ring systems are connected end-to-end
            // (e.g., bisazo dyes, polymeric chromophores, peptide
            // analogs), align the topological backbone (longest path
            // through the ring-system graph) to horizontal. This fixes
            // the bisazo dye TARGET test where alignToLongestAxis alone
            // gives 1.86 aspect — PCA on atom positions weights all
            // atoms equally, so substituents skew the principal axis
            // away from the actual molecular rod.
            try {
                if (Layout.options.linearizeMacroChain !== false &&
                    ringSystems && ringSystems.length >= 2) {
                    linearizeMacroChain(mol, atomIds, rings, ringSystems, ringAtomSet);
                }
            } catch (e5) { /* keep layout if linearisation fails */ }
        }
    }

    /**
     * v1.8.19: Macro-chain linearisation.
     *
     * For molecules with ≥ 2 ring systems connected by chains, finds the
     * longest topological path through the ring-system graph and rotates
     * the whole molecule so that path lies along the x-axis.
     *
     * Algorithm:
     *   1. Build a graph G:
     *        nodes  = ring system indices
     *        edges  = (a, b) if any atom in system a is bonded to any
     *                 atom in system b through a chain of ≤ 4 non-ring atoms
     *   2. Find the longest path P in G via BFS-extended DFS (treat as
     *      a tree if possible; otherwise pick the diameter of the graph).
     *   3. Compute the centroids of the two endpoint ring systems on P.
     *   4. Compute the angle θ of the line from start centroid to end
     *      centroid.
     *   5. Rotate every atom by -θ about the molecule's centroid.
     *
     * The rotation is rigid — bond lengths, ring shapes, and stereo
     * orientations are preserved exactly. Only the orientation changes.
     *
     * Idempotent: a molecule whose backbone is already horizontal
     * produces θ ≈ 0 and the rotation is a no-op.
     */
    function linearizeMacroChain(mol, atomIds, rings, ringSystems, ringAtomSet) {
        if (!ringSystems || ringSystems.length < 2) return;
        var n = ringSystems.length;

        // Build per-system atom set + centroid.
        var sysAtomSet = [];
        var sysCentroid = [];
        for (var si = 0; si < n; si++) {
            var sa = {};
            var cx = 0, cy = 0, cn = 0;
            for (var ri = 0; ri < ringSystems[si].length; ri++) {
                var ring = rings[ringSystems[si][ri]];
                if (!ring) continue;
                for (var ai = 0; ai < ring.length; ai++) {
                    if (!sa[ring[ai]]) {
                        sa[ring[ai]] = true;
                        var aa = mol.getAtom(ring[ai]);
                        if (aa) { cx += aa.x; cy += aa.y; cn++; }
                    }
                }
            }
            sysAtomSet.push(sa);
            sysCentroid.push({ x: cn ? cx / cn : 0, y: cn ? cy / cn : 0, n: cn });
        }

        // Build atom → system index lookup (fast inner-loop check).
        var atomToSys = {};
        for (var si2 = 0; si2 < n; si2++) {
            for (var k in sysAtomSet[si2]) atomToSys[k] = si2;
        }

        // Build adjacency between systems: BFS from each system's atoms
        // through chain (non-ring) atoms; if we reach an atom of another
        // system within 4 hops, add an edge.
        var adj = [];
        for (var i = 0; i < n; i++) adj.push({});
        var MAX_CHAIN_HOPS = 5;
        for (var si3 = 0; si3 < n; si3++) {
            var queue = [];
            var visited = {};
            for (var k2 in sysAtomSet[si3]) {
                queue.push({ id: parseInt(k2, 10), depth: 0 });
                visited[k2] = true;
            }
            while (queue.length > 0) {
                var cur = queue.shift();
                if (cur.depth > MAX_CHAIN_HOPS) continue;
                var nbrs = mol.getNeighbors(cur.id) || [];
                for (var nb = 0; nb < nbrs.length; nb++) {
                    var nbid = nbrs[nb];
                    if (visited[nbid]) continue;
                    visited[nbid] = true;
                    if (atomToSys[nbid] !== undefined && atomToSys[nbid] !== si3) {
                        var other = atomToSys[nbid];
                        adj[si3][other] = true;
                        adj[other][si3] = true;
                    } else if (!ringAtomSet[nbid]) {
                        // Walk through chain atoms.
                        queue.push({ id: nbid, depth: cur.depth + 1 });
                    }
                }
            }
        }

        // Find longest path via 2× BFS (tree-diameter algorithm; works
        // for any connected graph, returning a path that's at least as
        // long as the diameter).
        function bfsFar(startNode) {
            var dist = {};
            dist[startNode] = 0;
            var q = [startNode];
            var farthest = startNode, farthestDist = 0;
            var parent = {};
            parent[startNode] = -1;
            while (q.length > 0) {
                var u = q.shift();
                for (var v in adj[u]) {
                    var vi = parseInt(v, 10);
                    if (dist[vi] === undefined) {
                        dist[vi] = dist[u] + 1;
                        parent[vi] = u;
                        q.push(vi);
                        if (dist[vi] > farthestDist) {
                            farthestDist = dist[vi];
                            farthest = vi;
                        }
                    }
                }
            }
            return { node: farthest, dist: farthestDist, parent: parent };
        }
        var bfs1 = bfsFar(0);
        if (bfs1.dist === 0) return;  // disconnected ring systems
        var bfs2 = bfsFar(bfs1.node);
        if (bfs2.dist < 1) return;
        // Reconstruct the path from start (bfs1.node) to end (bfs2.node)
        // using bfs2's parent pointers.
        var path = [];
        var cur = bfs2.node;
        while (cur !== undefined && cur !== -1) {
            path.unshift(cur);
            cur = bfs2.parent[cur];
        }
        if (path.length < 2) return;

        // Snapshot ALL atoms so we can revert on aspect regression.
        var snap = atomIds.map(function (id) {
            var a = mol.getAtom(id);
            return a ? { id: id, x: a.x, y: a.y } : null;
        });
        function bboxAspect() {
            var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
            for (var i = 0; i < atomIds.length; i++) {
                var aa = mol.getAtom(atomIds[i]);
                if (!aa) continue;
                if (aa.x < minX) minX = aa.x;
                if (aa.x > maxX) maxX = aa.x;
                if (aa.y < minY) minY = aa.y;
                if (aa.y > maxY) maxY = aa.y;
            }
            var w = maxX - minX, h = maxY - minY;
            if (w < 1e-6 || h < 1e-6) return 1;
            return (w >= h) ? w / h : h / w;
        }
        function restoreSnap() {
            for (var sn = 0; sn < snap.length; sn++) {
                if (!snap[sn]) continue;
                var sa = mol.getAtom(snap[sn].id);
                if (sa) { sa.x = snap[sn].x; sa.y = snap[sn].y; }
            }
        }
        function rotateAll(theta) {
            var mcx = 0, mcy = 0;
            for (var ai2 = 0; ai2 < atomIds.length; ai2++) {
                var ma = mol.getAtom(atomIds[ai2]);
                if (ma) { mcx += ma.x; mcy += ma.y; }
            }
            mcx /= atomIds.length;
            mcy /= atomIds.length;
            var ct = Math.cos(-theta), st = Math.sin(-theta);
            for (var ai3 = 0; ai3 < atomIds.length; ai3++) {
                var mb = mol.getAtom(atomIds[ai3]);
                if (!mb) continue;
                var rx = mb.x - mcx, ry = mb.y - mcy;
                mb.x = mcx + rx * ct - ry * st;
                mb.y = mcy + rx * st + ry * ct;
            }
        }
        var preAspect = bboxAspect();

        // Strategy 1: PCA principal axis through the path centroids.
        var pcx = 0, pcy = 0;
        for (var pi = 0; pi < path.length; pi++) {
            pcx += sysCentroid[path[pi]].x;
            pcy += sysCentroid[path[pi]].y;
        }
        pcx /= path.length;
        pcy /= path.length;
        var sxx = 0, syy = 0, sxy = 0;
        for (var pi2 = 0; pi2 < path.length; pi2++) {
            var ddx2 = sysCentroid[path[pi2]].x - pcx;
            var ddy2 = sysCentroid[path[pi2]].y - pcy;
            sxx += ddx2 * ddx2;
            syy += ddy2 * ddy2;
            sxy += ddx2 * ddy2;
        }
        var theta1 = (sxx + syy < 1e-9) ? 0 : 0.5 * Math.atan2(2 * sxy, sxx - syy);

        // Strategy 2: start-to-end vector through the path's two
        // endpoint centroids.
        var startC = sysCentroid[path[0]];
        var endC = sysCentroid[path[path.length - 1]];
        var theta2 = Math.atan2(endC.y - startC.y, endC.x - startC.x);

        // Try each rotation; keep the one with the BEST resulting
        // aspect ratio. Revert if no rotation improves on preAspect.
        var bestAspect = preAspect;
        var bestSnapKept = false;
        var attempts = [{ theta: theta1, name: 'pca-path' },
                        { theta: theta2, name: 'start-end' }];
        for (var ai = 0; ai < attempts.length; ai++) {
            if (Math.abs(attempts[ai].theta) < 1e-6) continue;
            restoreSnap();
            rotateAll(attempts[ai].theta);
            var asp = bboxAspect();
            if (asp > bestAspect + 0.01) {
                bestAspect = asp;
                bestSnapKept = true;
                // Capture the resulting positions for "best" comparison.
                var bestSnap = atomIds.map(function (id) {
                    var a = mol.getAtom(id);
                    return a ? { id: id, x: a.x, y: a.y } : null;
                });
                // Restore the original snap so we can try the next
                // rotation cleanly; we'll re-apply bestSnap at the end.
                attempts[ai].snap = bestSnap;
            }
        }
        // Restore best result OR original if nothing improved.
        var winner = null;
        for (var ai2x = 0; ai2x < attempts.length; ai2x++) {
            if (attempts[ai2x].snap) winner = attempts[ai2x];
        }
        if (winner) {
            for (var k3 = 0; k3 < winner.snap.length; k3++) {
                if (!winner.snap[k3]) continue;
                var ka = mol.getAtom(winner.snap[k3].id);
                if (ka) { ka.x = winner.snap[k3].x; ka.y = winner.snap[k3].y; }
            }
        } else {
            restoreSnap();  // no improvement; keep pre-linearisation layout
            return;
        }

        // Final guarantee: ensure width ≥ height. If not, rotate by π/2.
        var minX2 = Infinity, maxX2 = -Infinity, minY2 = Infinity, maxY2 = -Infinity;
        var mcx2 = 0, mcy2 = 0;
        for (var fi = 0; fi < atomIds.length; fi++) {
            var fa = mol.getAtom(atomIds[fi]);
            if (!fa) continue;
            if (fa.x < minX2) minX2 = fa.x;
            if (fa.x > maxX2) maxX2 = fa.x;
            if (fa.y < minY2) minY2 = fa.y;
            if (fa.y > maxY2) maxY2 = fa.y;
            mcx2 += fa.x; mcy2 += fa.y;
        }
        mcx2 /= atomIds.length;
        mcy2 /= atomIds.length;
        if ((maxY2 - minY2) > (maxX2 - minX2)) {
            var ct2 = Math.cos(-Math.PI / 2), st2 = Math.sin(-Math.PI / 2);
            for (var fj = 0; fj < atomIds.length; fj++) {
                var fb = mol.getAtom(atomIds[fj]);
                if (!fb) continue;
                var rrx = fb.x - mcx2, rry = fb.y - mcy2;
                fb.x = mcx2 + rrx * ct2 - rry * st2;
                fb.y = mcy2 + rrx * st2 + rry * ct2;
            }
        }
    }

    /**
     * v1.8.16: chain-only iterative bond-length relax.
     *
     * For every bond where AT LEAST ONE endpoint is NOT a ring atom,
     * pull that endpoint toward BOND_LENGTH along the bond axis. Ring
     * atoms are FROZEN — they never move, even if the bond is stretched.
     * This guarantees the regular-polygon ring snap from Step 15 is
     * preserved exactly while chains get cleaned up.
     *
     * Algorithm: simulated bond-length spring with chain endpoints
     * absorbing 100 % of each bond's length error per iteration. Caps
     * at 30 iterations or when max deviation < 5 % of BOND_LENGTH.
     * Per-iteration accumulation (no race on shared atoms).
     */
    function relaxChainBonds(mol, atomIds, ringAtomSet) {
        var n = atomIds.length;
        if (n < 2) return;
        var idSet = {};
        for (var i = 0; i < n; i++) idSet[atomIds[i]] = true;
        // Collect bonds that touch at least one chain atom.
        var bonds = [];
        for (var i = 0; i < n; i++) {
            var nb = mol.getNeighbors(atomIds[i]);
            for (var j = 0; j < nb.length; j++) {
                if (idSet[nb[j]] && nb[j] > atomIds[i]) {
                    var u = atomIds[i], v = nb[j];
                    if (!ringAtomSet[u] || !ringAtomSet[v]) bonds.push([u, v]);
                }
            }
        }
        if (bonds.length === 0) return;

        var TOL = BOND_LENGTH * 0.05;
        var MAX_ITER = 80;
        // Per-iteration gain. Chain atoms with multiple bonds receive
        // contributions from each neighbour; gain 0.5 strikes a balance
        // between convergence speed and overshoot when a single chain
        // atom has 2-4 chain bonds all wanting to pull. Higher gains
        // (0.6+) overshoot for branched chains; lower (0.25) takes too
        // many iterations on long chains. Empirical: 0.5 hits 5%-of-BL
        // tolerance in 12-25 iters on the regression set.
        var GAIN = 0.5;
        for (var iter = 0; iter < MAX_ITER; iter++) {
            var maxDev = 0;
            var dx = {}, dy = {};
            for (var i = 0; i < n; i++) { dx[atomIds[i]] = 0; dy[atomIds[i]] = 0; }
            for (var bi = 0; bi < bonds.length; bi++) {
                var u = bonds[bi][0], v = bonds[bi][1];
                var au = mol.getAtom(u), av = mol.getAtom(v);
                if (!au || !av) continue;
                var ddx = av.x - au.x, ddy = av.y - au.y;
                var d = Math.sqrt(ddx * ddx + ddy * ddy);
                if (d < EPSILON) continue;
                var dev = d - BOND_LENGTH;
                if (Math.abs(dev) > maxDev) maxDev = Math.abs(dev);
                var nx = ddx / d, ny = ddy / d;
                var ru = !!ringAtomSet[u], rv = !!ringAtomSet[v];
                // Distribution: ring endpoint frozen (gain 0); chain
                // endpoint absorbs deviation × GAIN. If both are chain,
                // split evenly between them.
                if (ru && !rv) {
                    dx[v] -= nx * dev * GAIN;
                    dy[v] -= ny * dev * GAIN;
                } else if (rv && !ru) {
                    dx[u] += nx * dev * GAIN;
                    dy[u] += ny * dev * GAIN;
                } else if (!ru && !rv) {
                    dx[u] += nx * dev * GAIN * 0.5;
                    dy[u] += ny * dev * GAIN * 0.5;
                    dx[v] -= nx * dev * GAIN * 0.5;
                    dy[v] -= ny * dev * GAIN * 0.5;
                }
                // both ring → frozen (no entry in `bonds`).
            }
            // Apply only to chain atoms (defence in depth — `dx[u]` of a
            // ring atom should be 0 by construction, but skip explicitly).
            for (var i = 0; i < n; i++) {
                var aid = atomIds[i];
                if (ringAtomSet[aid]) continue;
                var a = mol.getAtom(aid);
                if (!a) continue;
                a.x += dx[aid];
                a.y += dy[aid];
            }
            if (maxDev < TOL) break;
        }
    }

    /**
     * v1.8.16: rigid translation of overlapping ring systems.
     *
     * After Step 15 ring snap each ring is a perfect regular polygon,
     * but two ring systems whose centroids drifted close during earlier
     * steps may now physically overlap (atoms within MIN_ATOM_DIST). We
     * fix that by RIGIDLY TRANSLATING the smaller-centroid-distance ring
     * system away — preserving polygon shape exactly. Substituents move
     * with their parent ring system.
     *
     * O(R²) over ring systems R; for typical molecules R ≤ 4-6 so this
     * is well under 50 µs even on big polycyclics.
     */
    function separateOverlappingRingSystems(mol, atomIds, ringSystems, rings, ringAtomSet) {
        if (!ringSystems || ringSystems.length < 2) return;
        // Build per-system atom set.
        var sysAtomSets = [];
        var sysAtomLists = [];
        var sysCentroids = [];
        for (var si = 0; si < ringSystems.length; si++) {
            var sysAtoms = {};
            var sysList = [];
            for (var ri = 0; ri < ringSystems[si].length; ri++) {
                var ringIds = rings[ringSystems[si][ri]];
                if (!ringIds) continue;
                for (var ai = 0; ai < ringIds.length; ai++) {
                    if (!sysAtoms[ringIds[ai]]) {
                        sysAtoms[ringIds[ai]] = true;
                        sysList.push(ringIds[ai]);
                    }
                }
            }
            sysAtomSets.push(sysAtoms);
            sysAtomLists.push(sysList);
            // Centroid.
            var cx = 0, cy = 0, cn = 0;
            for (var li = 0; li < sysList.length; li++) {
                var a = mol.getAtom(sysList[li]);
                if (a) { cx += a.x; cy += a.y; cn++; }
            }
            sysCentroids.push({ x: cn ? cx / cn : 0, y: cn ? cy / cn : 0 });
        }

        // Iterate up to 12 separation passes. Each pass scans all pairs;
        // for any pair that overlaps, translate system B along the
        // CENTROID-TO-CENTROID axis by the WORST per-atom-pair deficit.
        // For collinear (centroid-axis-aligned) overlaps this clears the
        // worst pair in one push; for laterally-staggered overlaps
        // residual overlaps may remain and are cleared by subsequent
        // iterations. The 12-iteration outer cap is the convergence
        // guarantee. Worst-pair deficit > average-deficit because a
        // single 14.4 px pair would otherwise be averaged with several
        // 16.0 px pairs and the system would barely move.
        // Complexity: O(SEP_ITER × R² × S²) where R = #ring systems,
        // S = atoms per system; for typical drug-like inputs (R ≤ 6,
        // S ≤ 30) ≈ 130 K atom-pair checks per layout.
        var SEP_ITER = 12;
        for (var iter = 0; iter < SEP_ITER; iter++) {
            var anyMoved = false;
            for (var i = 0; i < ringSystems.length; i++) {
                for (var j = i + 1; j < ringSystems.length; j++) {
                    var listA = sysAtomLists[i], listB = sysAtomLists[j];
                    var worstDeficit = 0;
                    for (var ai2 = 0; ai2 < listA.length; ai2++) {
                        var aA = mol.getAtom(listA[ai2]);
                        if (!aA) continue;
                        for (var bi2 = 0; bi2 < listB.length; bi2++) {
                            var aB = mol.getAtom(listB[bi2]);
                            if (!aB) continue;
                            var ddx = aB.x - aA.x, ddy = aB.y - aA.y;
                            var dd = Math.hypot(ddx, ddy);
                            if (dd < MIN_ATOM_DIST) {
                                var deficit = MIN_ATOM_DIST - dd;
                                if (deficit > worstDeficit) worstDeficit = deficit;
                            }
                        }
                    }
                    if (worstDeficit < 0.01) continue;
                    // Push along the centroid-to-centroid axis.
                    var cdx = sysCentroids[j].x - sysCentroids[i].x;
                    var cdy = sysCentroids[j].y - sysCentroids[i].y;
                    var clen = Math.hypot(cdx, cdy);
                    var pushX, pushY;
                    if (clen < EPSILON) {
                        // Centroids co-located: pick a deterministic axis.
                        var ang2 = (i * 0.7137 + j * 1.31) * Math.PI;
                        pushX = Math.cos(ang2) * worstDeficit;
                        pushY = Math.sin(ang2) * worstDeficit;
                    } else {
                        pushX = (cdx / clen) * worstDeficit;
                        pushY = (cdy / clen) * worstDeficit;
                    }
                    for (var bi3 = 0; bi3 < listB.length; bi3++) {
                        var aBmove = mol.getAtom(listB[bi3]);
                        if (aBmove) {
                            aBmove.x += pushX;
                            aBmove.y += pushY;
                        }
                    }
                    sysCentroids[j].x += pushX;
                    sysCentroids[j].y += pushY;
                    anyMoved = true;
                }
            }
            if (!anyMoved) break;
        }
    }

    /**
     * v1.8.24: count non-bonded atom pairs closer than BOND_LENGTH/3.
     * Used by the pull-back safety net to decide whether the
     * connection-aware pull regressed the layout. O(n) via spatial grid.
     */
    function countAtomOverlaps(mol, atomIds) {
        var n = atomIds.length;
        if (n < 2) return 0;
        var bonded = {};
        for (var b = 0; b < mol.bonds.length; b++) {
            var bn = mol.bonds[b];
            bonded[bn.atom1 + ',' + bn.atom2] = true;
            bonded[bn.atom2 + ',' + bn.atom1] = true;
        }
        var thr = BOND_LENGTH / 3;
        var thr2 = thr * thr;
        var grid = {};
        for (var i = 0; i < n; i++) {
            var a = mol.getAtom(atomIds[i]);
            if (!a) continue;
            var key = Math.floor(a.x / GRID_CELL) + ',' + Math.floor(a.y / GRID_CELL);
            if (!grid[key]) grid[key] = [];
            grid[key].push(i);
        }
        var count = 0;
        for (var i2 = 0; i2 < n; i2++) {
            var a1 = mol.getAtom(atomIds[i2]);
            if (!a1) continue;
            var gx = Math.floor(a1.x / GRID_CELL), gy = Math.floor(a1.y / GRID_CELL);
            for (var dx = -1; dx <= 1; dx++) {
                for (var dy = -1; dy <= 1; dy++) {
                    var cell = grid[(gx + dx) + ',' + (gy + dy)];
                    if (!cell) continue;
                    for (var ci = 0; ci < cell.length; ci++) {
                        var j = cell[ci];
                        if (j <= i2) continue;
                        if (bonded[atomIds[i2] + ',' + atomIds[j]]) continue;
                        var a2 = mol.getAtom(atomIds[j]);
                        if (!a2) continue;
                        var ddx = a2.x - a1.x, ddy = a2.y - a1.y;
                        if (ddx * ddx + ddy * ddy < thr2) count++;
                    }
                }
            }
        }
        return count;
    }

    /**
     * v1.8.24: connection-aware ring-system pull-back.
     *
     * separateOverlappingRingSystems clears OVERLAPS between ring
     * systems but is blind to the bonds that CONNECT them. The result
     * is that a directly-bonded biaryl (biphenyl: ring-atom to
     * ring-atom, observed 2.54x BOND_LENGTH) or a one-atom linker
     * (diphenylmethane: ring-CH2-ring, observed 1.90x) is left badly
     * stretched — the worst single defect in BIME's chain depiction.
     *
     * v1.8.32: the old centroid-direction translation pull could not
     * recover a biaryl junction at all — the rule-based layout often
     * leaves the two rings crammed face-to-face with the junction bond
     * stretched, so translating one ring straight at the other only
     * deepens the clash and the binary-search backs the pull off to
     * zero. A pure translation cannot fix a mis-ORIENTED placement.
     *
     * This pass instead does a spanning-tree RIGID PLACEMENT (translate
     * + rotate):
     *   1. Build a spanning tree over the ring systems — an edge for
     *      every linker pair (a direct bond, or a ≤3-atom chain) found
     *      by findLinker — rooted at the largest system.
     *   2. Walk the tree outward. For every tree edge that is a DIRECT
     *      junction bond (linkLen === 1, the biaryl case), rigidly
     *      transform the child system AND its whole descendant subtree
     *      so the junction bond is exactly BOND_LENGTH and leaves the
     *      parent radially — i.e. the child extends away from the parent
     *      instead of folding back over it.
     *   3. Chain linkers (linkLen > 1) move by translation when their
     *      end-to-end span is clearly compressed or stretched; then
     *      relaxChainBonds re-centres their interior atoms.
     *
     * Moving whole subtrees in tree order is deterministic and cannot
     * tug-of-war: placing a parent-child bond never disturbs an
     * already-placed bond, because everything below the child travels
     * with it. Every ring polygon and internal chain shape is preserved
     * exactly — only the rigid placement of whole ring systems changes.
     */
    function pullConnectedRingSystems(mol, atomIds, ringSystems, rings,
                                         ringAtomSet, atomToSys) {
        if (!ringSystems || ringSystems.length < 2) return;
        var nSys = ringSystems.length;

        // Per-system atom lists (rigid bodies).
        var sysAtoms = [];
        for (var si = 0; si < nSys; si++) {
            var seen = {};
            var list = [];
            for (var ri = 0; ri < ringSystems[si].length; ri++) {
                var rg = rings[ringSystems[si][ri]];
                if (!rg) continue;
                for (var ai = 0; ai < rg.length; ai++) {
                    if (!seen[rg[ai]]) { seen[rg[ai]] = true; list.push(rg[ai]); }
                }
            }
            sysAtoms.push(list);
        }

        // Find linkers between every system pair. A linker is the
        // shortest atom path u..v with u a ring atom of system i, v a
        // ring atom of system j, and every interior atom a non-ring
        // (chain) atom. BFS from each system's ring atoms, depth <= 3.
        function findLinker(si, sj) {
            // BFS: queue of {id, depth, originRingAtom}.
            var q = [], visited = {};
            for (var s = 0; s < sysAtoms[si].length; s++) {
                q.push({ id: sysAtoms[si][s], depth: 0, origin: sysAtoms[si][s] });
                visited[sysAtoms[si][s]] = true;
            }
            while (q.length > 0) {
                var cur = q.shift();
                if (cur.depth > 3) continue;
                var nbrs = mol.getNeighbors(cur.id) || [];
                for (var nb = 0; nb < nbrs.length; nb++) {
                    var nid = nbrs[nb];
                    if (visited[nid]) continue;
                    if (atomToSys[nid] === sj) {
                        // Reached system j. Linker length = cur.depth + 1 bonds.
                        return { u: cur.origin, v: nid, linkLen: cur.depth + 1 };
                    }
                    if (atomToSys[nid] === undefined && !ringAtomSet[nid] && cur.depth < 3) {
                        // chain atom — keep walking
                        visited[nid] = true;
                        q.push({ id: nid, depth: cur.depth + 1, origin: cur.origin });
                    }
                }
            }
            return null;
        }

        // --- system adjacency: one edge per linked system pair ----------
        var sysEdges = [];                 // { a, b, link }
        var sysAdj = [];
        for (var ae = 0; ae < nSys; ae++) sysAdj.push([]);
        for (var ei = 0; ei < nSys; ei++) {
            for (var ej = ei + 1; ej < nSys; ej++) {
                var lk = findLinker(ei, ej);
                if (!lk) continue;
                var edge = { a: ei, b: ej, link: lk };
                sysEdges.push(edge);
                sysAdj[ei].push({ other: ej, edge: edge });
                sysAdj[ej].push({ other: ei, edge: edge });
            }
        }
        if (sysEdges.length === 0) return;

        // --- spanning tree, rooted at the largest system ----------------
        // The biggest rigid body stays put; smaller satellites move.
        var rootSys = 0;
        for (var rs = 1; rs < nSys; rs++) {
            if (sysAtoms[rs].length > sysAtoms[rootSys].length) rootSys = rs;
        }
        var parentOf = {};                 // childSys -> { parent, edge }
        var bfsOrder = [];
        var seenSys = {}; seenSys[rootSys] = true;
        var sq = [rootSys];
        while (sq.length > 0) {
            var curSys = sq.shift();
            for (var ad = 0; ad < sysAdj[curSys].length; ad++) {
                var nb2 = sysAdj[curSys][ad];
                if (seenSys[nb2.other]) continue;
                seenSys[nb2.other] = true;
                parentOf[nb2.other] = { parent: curSys, edge: nb2.edge };
                bfsOrder.push(nb2.other);
                sq.push(nb2.other);
            }
        }

        // --- descendant-subtree atom list for every system --------------
        // Built in reverse BFS order so each child is complete before its
        // parent folds it in.
        var subtreeAtoms = {};
        for (var st = 0; st < nSys; st++) subtreeAtoms[st] = sysAtoms[st].slice();
        for (var rb = bfsOrder.length - 1; rb >= 0; rb--) {
            var csr = bfsOrder[rb];
            var par = parentOf[csr].parent;
            subtreeAtoms[par] = subtreeAtoms[par].concat(subtreeAtoms[csr]);
        }

        function centroidOf(ids) {
            var cx = 0, cy = 0, nn = 0;
            for (var c = 0; c < ids.length; c++) {
                var a = mol.getAtom(ids[c]);
                if (a) { cx += a.x; cy += a.y; nn++; }
            }
            return nn ? { x: cx / nn, y: cy / nn } : { x: 0, y: 0 };
        }

        // --- walk the tree outward, rigidly place each DIRECT-bond child -
        var MIN_D2 = MIN_ATOM_DIST * MIN_ATOM_DIST;   // chain-linker clash guard
        for (var bo = 0; bo < bfsOrder.length; bo++) {
            var childSys = bfsOrder[bo];
            var pe = parentOf[childSys];
            var linkLen = pe.edge.link.linkLen;

            // Orient the junction as u (in parent) -> v (in child).
            var lU = pe.edge.link.u, lV = pe.edge.link.v;
            var u = (atomToSys[lU] === pe.parent) ? lU : lV;
            var v = (u === lU) ? lV : lU;
            var au = mol.getAtom(u), av = mol.getAtom(v);
            if (!au || !av) continue;
            var sub = subtreeAtoms[childSys];

            if (linkLen === 1) {
                // DIRECT junction bond (the biaryl case). The rule-based
                // layout reliably mis-places AND mis-orients it — the two
                // rings end up crammed face-to-face with the bond stretched
                // — so a rigid TRANSLATE + ROTATE is needed: put the bond at
                // exactly BOND_LENGTH, leaving the parent radially so the
                // child extends away instead of folding back over it.
                // This fires for EVERY direct junction bond, not only the
                // visibly-broken ones — radial-outward placement at exactly
                // BOND_LENGTH is a correct placement regardless, and the
                // regression suite confirms it: clean-ring rate
                // unchanged, atom overlaps strictly down. Intentional.
                var pc = centroidOf(sysAtoms[pe.parent]);
                var dxv = au.x - pc.x, dyv = au.y - pc.y;
                var dlen = Math.hypot(dxv, dyv);
                if (dlen < EPSILON) {
                    // u sits at the parent centroid — fall back to the
                    // current junction direction so there is still an axis.
                    dxv = av.x - au.x; dyv = av.y - au.y;
                    dlen = Math.hypot(dxv, dyv);
                    if (dlen < EPSILON) continue;
                }
                dxv /= dlen; dyv /= dlen;
                var vtx = au.x + dxv * BOND_LENGTH;
                var vty = au.y + dyv * BOND_LENGTH;
                // 1. translate the whole subtree so v lands on the target.
                var tx = vtx - av.x, ty = vty - av.y;
                for (var s1 = 0; s1 < sub.length; s1++) {
                    var sa = mol.getAtom(sub[s1]);
                    if (sa) { sa.x += tx; sa.y += ty; }
                }
                // 2. rotate the subtree about v so the child system's
                //    centroid lies along d.
                var cc = centroidOf(sysAtoms[childSys]);
                var ccdx = cc.x - vtx, ccdy = cc.y - vty;
                var cclen = Math.hypot(ccdx, ccdy);
                if (cclen < EPSILON) continue;
                ccdx /= cclen; ccdy /= cclen;
                // rotation taking unit ccDir onto unit d:
                //   cos = ccDir . d ,  sin = ccDir x d  (z-component).
                var cosT = ccdx * dxv + ccdy * dyv;
                var sinT = ccdx * dyv - ccdy * dxv;
                for (var s2 = 0; s2 < sub.length; s2++) {
                    var ra = mol.getAtom(sub[s2]);
                    if (!ra) continue;
                    var rx = ra.x - vtx, ry = ra.y - vty;
                    ra.x = vtx + rx * cosT - ry * sinT;
                    ra.y = vty + rx * sinT + ry * cosT;
                }
            } else {
                // CHAIN linker (>=2 bonds, e.g. -CH2-, -N=N-). The
                // rule-based layout orients these acceptably but can leave
                // the child system either too far out (stretched chain) or
                // too close (diphenylmethane: ring-CH2-ring). Translate the
                // child subtree along the CURRENT junction axis until the
                // end-to-end span matches a natural zigzag — no rotation,
                // the orientation is already fine — then relaxChainBonds
                // re-centres the interior chain atoms.
                var span = Math.hypot(av.x - au.x, av.y - au.y);
                if (span < EPSILON) continue;
                // 0.92 ~= mean projection of a 120-degree zigzag onto its
                // end-to-end axis (matches relaxChainBonds' rest geometry).
                var targetSpan = linkLen * BOND_LENGTH * 0.92;
                if (Math.abs(span - targetSpan) <= targetSpan * 0.05) continue;
                var sdx = (av.x - au.x) / span, sdy = (av.y - au.y) / span;
                var ctx = (au.x + sdx * targetSpan) - av.x;
                var cty = (au.y + sdy * targetSpan) - av.y;
                // Expanding an over-compressed one-atom linker moves the
                // child away from the parent; the outer snapshot gate below
                // still catches rare third-system regressions. Longer
                // flexible linkers keep their base placement because their
                // natural end-to-end span depends heavily on conformation.
                if (span < targetSpan) {
                    if (linkLen !== 2) continue;
                    for (var es3 = 0; es3 < sub.length; es3++) {
                        var ea2 = mol.getAtom(sub[es3]);
                        if (ea2) { ea2.x += ctx; ea2.y += cty; }
                    }
                    continue;
                }
                // Clash guard: pulling a subtree IN can cram it into the
                // parent or a sibling. Binary-search the largest fraction
                // of the pull that keeps every moved atom >= MIN_ATOM_DIST
                // from every non-moved atom. The two junction-bond
                // endpoints are exempt — u (the parent anchor, excluded
                // from `rest`) and v (the child anchor, skipped in the
                // moved-atom loop) are MEANT to come closer together; the
                // interior chain atoms between them get re-centred by
                // relaxChainBonds right after this pass.
                var inSub = {};
                for (var q1 = 0; q1 < sub.length; q1++) inSub[sub[q1]] = true;
                var loC = 0, hiC = 1, bestC = 0;
                for (var bsC = 0; bsC < 8; bsC++) {
                    var midC = (loC + hiC) / 2;
                    var clashC = false;
                    for (var p1 = 0; p1 < sub.length && !clashC; p1++) {
                        if (sub[p1] === v) continue;       // child anchor exempt
                        var pa = mol.getAtom(sub[p1]); if (!pa) continue;
                        var npx = pa.x + ctx * midC, npy = pa.y + cty * midC;
                        // Test against every non-moved atom, inline — exempt
                        // the subtree (moves as one body) and u (the parent
                        // anchor, meant to come closer). No per-child `rest`
                        // array is materialised.
                        for (var qx = 0; qx < atomIds.length; qx++) {
                            var rid = atomIds[qx];
                            if (inSub[rid] || rid === u) continue;
                            var qa = mol.getAtom(rid); if (!qa) continue;
                            var ddx = npx - qa.x, ddy = npy - qa.y;
                            if (ddx * ddx + ddy * ddy < MIN_D2) { clashC = true; break; }
                        }
                    }
                    if (clashC) { hiC = midC; } else { bestC = midC; loC = midC; }
                }
                if (bestC > 0.01) {
                    for (var s3 = 0; s3 < sub.length; s3++) {
                        var ca2 = mol.getAtom(sub[s3]);
                        if (ca2) { ca2.x += ctx * bestC; ca2.y += cty * bestC; }
                    }
                }
            }
        }
    }

    /**
     * v1.8.16: chain-only collision resolver. Ring atoms are NEVER moved.
     * Uses the same spatial grid + 3x3 neighbourhood scan as
     * resolveCollisions but skips any pair where BOTH atoms are ring
     * atoms (already handled by separateOverlappingRingSystems) and
     * any push that would move a ring atom (only the chain endpoint
     * absorbs the push).
     */
    // (v1.8.25: the former `resolveChainOverlaps15` was folded into the
    // unified `resolveCollisions` above — Step 15 now calls it with
    // { ringGentleFactor: 0.15, atomToSys: <map>, snapshotRings: false }.)

    // ---------------------------------------------------------------------
    // v1.8.15 Step 15 helpers — fused-ring-aware polygon snap.
    //
    // These three functions are NOT used by the initial layout pipeline
    // (Steps 1-12). They are only invoked from the Step 15 final-snap
    // block above. Keeping them at module scope (rather than nested) so
    // the bundle minifier produces a clean call site and so each helper
    // has its own stack frame in the JS profiler.
    //
    // All three return TRUE when at least one atom moved (so Step 15 can
    // tell if any substituent translation needs to flow downstream), or
    // FALSE if the ring was already at exact regular-polygon positions.
    // ---------------------------------------------------------------------

    /**
     * Snap a single ring as a regular polygon at its CURRENT centroid and
     * its CURRENT first-atom angle. Idempotent — already-clean rings
     * produce no positional change.
     */
    function snapRingAtCentroid(mol, ring, bondLength) {
        var rids = ring.atoms || ring;
        if (!rids || rids.length < 3) return false;
        var n = rids.length;

        // Measure current shape.
        var minLen = Infinity, maxLen = 0;
        for (var i = 0; i < n; i++) {
            var ra = mol.getAtom(rids[i]);
            var rb = mol.getAtom(rids[(i + 1) % n]);
            if (!ra || !rb) continue;
            var d = Math.hypot(rb.x - ra.x, rb.y - ra.y);
            if (d < minLen) minLen = d;
            if (d > maxLen) maxLen = d;
        }
        // Already clean? Skip — exit fast and report no change.
        // EXCEPT for size-8+ rings: even a "clean" regular polygon
        // should be replaced with the egg-shape ovalisation from
        // SDG.MacroCycleLayout (better visual depiction).
        var sdgAvailable = (global.SDG &&
                            global.SDG.MacroCycleLayout && global.SDG.MacroCycleLayout.layout);
        if (n < 8 || !sdgAvailable) {
            if (maxLen > 0 && maxLen / minLen < 1.02 &&
                Math.abs(minLen - bondLength) < bondLength * 0.05) return false;
        }
        // For size 8+: continue to the egg-shape path even if already
        // a regular polygon. The MacroCycleLayout call below will be a
        // visible no-op if the ring already matches egg-shape, but for
        // a freshly-laid regular polygon it produces the visually
        // preferable ovalised geometry.

        // Current centroid.
        var cx = 0, cy = 0, c = 0;
        for (var i2 = 0; i2 < n; i2++) {
            var a2 = mol.getAtom(rids[i2]);
            if (a2) { cx += a2.x; cy += a2.y; c++; }
        }
        if (c === 0) return false;
        cx /= c; cy /= c;

        var first = mol.getAtom(rids[0]);
        var startAng = (first ? Math.atan2(first.y - cy, first.x - cx) : 0);
        if (!isFinite(startAng)) startAng = 0;

        // v1.8.18: for size-8+ rings, delegate to SDG.MacroCycleLayout
        // which uses arc-length parameterisation + chord-equalisation
        // refinement (egg-shape ovalisation). For size 3-7 the regular
        // polygon below produces a strict polygon — egg-shape would
        // reduce to a regular polygon at α=1 anyway.
        var useMacroCycle = (n >= 8 &&
                              global.SDG &&
                              global.SDG.MacroCycleLayout &&
                              global.SDG.MacroCycleLayout.layout);
        if (useMacroCycle) {
            var ok = global.SDG.MacroCycleLayout.layout(rids, {
                mol: mol,
                bondLength: bondLength,
                cx: cx, cy: cy,
                startAng: startAng
            });
            if (ok) return true;
            // If MacroCycleLayout fails, fall through to regular polygon.
        }

        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        var moved = false;
        for (var i3 = 0; i3 < n; i3++) {
            var a3 = mol.getAtom(rids[i3]);
            if (!a3) continue;
            var ang = startAng + (TWO_PI * i3) / n;
            var nX = cx + radius * Math.cos(ang);
            var nY = cy + radius * Math.sin(ang);
            if (Math.abs(nX - a3.x) > 1e-9 || Math.abs(nY - a3.y) > 1e-9) moved = true;
            a3.x = nX; a3.y = nY;
        }
        return moved;
    }

    /**
     * Snap a ring fused to an already-snapped neighbour by EDGE
     * REFLECTION across the two shared atoms. The shared edge is now
     * fixed at exactly BOND_LENGTH (since the parent ring was snapped to
     * a regular polygon), so reflecting a regular polygon across it
     * yields a perfect regular polygon for this ring too.
     *
     * @param sharedIds  the (≥ 2) atom IDs already in snappedAtomSet —
     *                   we use the first 2 as the shared-edge endpoints
     *                   (the BFS scoring above prefers contiguous pairs).
     */
    function snapFusedRingByReflection(mol, ring, sharedIds, snappedAtomSet, bondLength) {
        var rids = ring.atoms || ring;
        if (!rids || rids.length < 3) return false;
        var n = rids.length;

        // Find the shared-edge indices in the ring. Prefer two
        // contiguous indices; if 3+ shared, pick any contiguous pair.
        var sharedIdx = [];
        for (var i = 0; i < n; i++) {
            if (sharedIds.indexOf(rids[i]) >= 0) sharedIdx.push(i);
        }
        if (sharedIdx.length < 2) return false;
        // Find a contiguous pair (idx_i and idx_{i+1 mod n}).
        var contigPair = null;
        for (var ci = 0; ci < sharedIdx.length; ci++) {
            for (var cj = 0; cj < sharedIdx.length; cj++) {
                if (ci === cj) continue;
                if ((sharedIdx[cj] - sharedIdx[ci] + n) % n === 1) {
                    contigPair = [sharedIdx[ci], sharedIdx[cj]];
                    break;
                }
            }
            if (contigPair) break;
        }
        if (!contigPair) return false;
        var idx1 = contigPair[0], idx2 = contigPair[1];

        var a1 = mol.getAtom(rids[idx1]);
        var a2 = mol.getAtom(rids[idx2]);
        if (!a1 || !a2) return false;

        // Edge midpoint and unit normal.
        var mx = (a1.x + a2.x) / 2, my = (a1.y + a2.y) / 2;
        var edx = a2.x - a1.x, edy = a2.y - a1.y;
        var nrm = Math.hypot(edx, edy);
        if (nrm < EPSILON) return false;
        var nx = -edy / nrm, ny = edx / nrm;

        // Two candidate centres at apothem distance on either side.
        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        var apothem = radius * Math.cos(Math.PI / n);
        var cx1 = mx + nx * apothem, cy1 = my + ny * apothem;
        var cx2 = mx - nx * apothem, cy2 = my - ny * apothem;

        // Pick the side farther from already-snapped atoms (so the new
        // ring lands OPPOSITE the parent ring, not on top of it).
        var s1 = 0, s2 = 0;
        for (var pid in snappedAtomSet) {
            if (sharedIds.indexOf(parseInt(pid, 10)) >= 0) continue;
            var pa = mol.getAtom(parseInt(pid, 10));
            if (!pa) continue;
            var d1sq = (cx1 - pa.x) * (cx1 - pa.x) + (cy1 - pa.y) * (cy1 - pa.y);
            var d2sq = (cx2 - pa.x) * (cx2 - pa.x) + (cy2 - pa.y) * (cy2 - pa.y);
            if (d1sq < bondLength * bondLength * 0.6) s1 -= 200;
            if (d2sq < bondLength * bondLength * 0.6) s2 -= 200;
            s1 += d1sq;
            s2 += d2sq;
        }
        var cx = (s1 >= s2) ? cx1 : cx2;
        var cy = (s1 >= s2) ? cy1 : cy2;

        // Determine winding direction by comparing the angle from cx,cy
        // to a1 vs a2 against the expected step (1 vertex = 2π/n).
        var ang1 = Math.atan2(a1.y - cy, a1.x - cx);
        var ang2 = Math.atan2(a2.y - cy, a2.x - cx);
        var step = TWO_PI / n;
        var stepsForward = (idx2 - idx1 + n) % n;          // = 1 by construction
        var angleDiff = normalizeAngle(ang2 - ang1);
        var expectedFwd = stepsForward * step;
        var clockwise = (Math.abs(angleDiff - expectedFwd) <= Math.abs(angleDiff - (TWO_PI - expectedFwd)));

        var moved = false;
        for (var i = 0; i < n; i++) {
            var atomI = mol.getAtom(rids[i]);
            if (!atomI) continue;
            var steps = (i - idx1 + n) % n;
            var ang = clockwise ? ang1 + steps * step : ang1 - steps * step;
            var nX = cx + radius * Math.cos(ang);
            var nY = cy + radius * Math.sin(ang);
            if (Math.abs(nX - atomI.x) > 1e-9 || Math.abs(nY - atomI.y) > 1e-9) moved = true;
            atomI.x = nX; atomI.y = nY;
        }
        return moved;
    }

    /**
     * Snap a spiro ring (one shared atom) as a regular polygon rotated
     * around the shared atom. The ring axis is chosen perpendicular to
     * the bisector of the two existing in-parent-ring bonds at the spiro
     * atom (so the new ring extends "outward" — same heuristic as
     * placeSpiroRing).
     */
    function snapSpiroRingFromAtom(mol, ring, spiroId, snappedAtomSet, bondLength) {
        var rids = ring.atoms || ring;
        if (!rids || rids.length < 3) return false;
        var n = rids.length;
        var spiroIdx = -1;
        for (var i = 0; i < n; i++) {
            if (rids[i] === spiroId) { spiroIdx = i; break; }
        }
        if (spiroIdx < 0) return false;
        var spiroAtom = mol.getAtom(spiroId);
        if (!spiroAtom) return false;

        // Find the snapped neighbours of the spiro atom (in OTHER rings).
        var nbrs = mol.getNeighbors(spiroId) || [];
        var snappedNbrs = [];
        var inThisRing = {};
        for (var j = 0; j < n; j++) inThisRing[rids[j]] = true;
        for (var k = 0; k < nbrs.length; k++) {
            if (inThisRing[nbrs[k]]) continue;
            if (snappedAtomSet[nbrs[k]]) snappedNbrs.push(mol.getAtom(nbrs[k]));
        }

        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        var awayAng;
        if (snappedNbrs.length >= 2) {
            // Bisector of existing bonds → perpendicular gives axis.
            var ba1 = Math.atan2(snappedNbrs[0].y - spiroAtom.y,
                                 snappedNbrs[0].x - spiroAtom.x);
            var ba2 = Math.atan2(snappedNbrs[1].y - spiroAtom.y,
                                 snappedNbrs[1].x - spiroAtom.x);
            var bisector = Math.atan2(Math.sin(ba1) + Math.sin(ba2),
                                      Math.cos(ba1) + Math.cos(ba2));
            // Pick whichever perpendicular is farther from snapped atoms.
            var p1 = bisector + Math.PI / 2;
            var p2 = bisector - Math.PI / 2;
            var cx1 = spiroAtom.x + radius * Math.cos(p1);
            var cy1 = spiroAtom.y + radius * Math.sin(p1);
            var cx2 = spiroAtom.x + radius * Math.cos(p2);
            var cy2 = spiroAtom.y + radius * Math.sin(p2);
            var s1 = 0, s2 = 0;
            for (var pid2 in snappedAtomSet) {
                if (parseInt(pid2, 10) === spiroId) continue;
                var pa2 = mol.getAtom(parseInt(pid2, 10));
                if (!pa2) continue;
                s1 += (cx1 - pa2.x) * (cx1 - pa2.x) + (cy1 - pa2.y) * (cy1 - pa2.y);
                s2 += (cx2 - pa2.x) * (cx2 - pa2.x) + (cy2 - pa2.y) * (cy2 - pa2.y);
            }
            awayAng = (s1 >= s2) ? p1 : p2;
        } else if (snappedNbrs.length === 1) {
            // Go directly opposite the single placed bond.
            awayAng = Math.atan2(spiroAtom.y - snappedNbrs[0].y,
                                 spiroAtom.x - snappedNbrs[0].x);
        } else {
            awayAng = 0;
        }

        var cx = spiroAtom.x + radius * Math.cos(awayAng);
        var cy = spiroAtom.y + radius * Math.sin(awayAng);
        var startAng = Math.atan2(spiroAtom.y - cy, spiroAtom.x - cx);
        var step = TWO_PI / n;

        var moved = false;
        for (var i2 = 0; i2 < n; i2++) {
            var atomI2 = mol.getAtom(rids[i2]);
            if (!atomI2) continue;
            var steps = (i2 - spiroIdx + n) % n;
            var ang = startAng + steps * step;
            var nX = cx + radius * Math.cos(ang);
            var nY = cy + radius * Math.sin(ang);
            if (Math.abs(nX - atomI2.x) > 1e-9 || Math.abs(nY - atomI2.y) > 1e-9) moved = true;
            atomI2.x = nX; atomI2.y = nY;
        }
        return moved;
    }

    // =====================================================================
    // SSSR Ring Perception  (edge-vector independence)
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
        // (BFS has issues with even-cycle sizes; DFS is simpler and correct)
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

        // Fast template matching: element-frequency pre-screen + greedy probe.
        if (typeof Templates === 'undefined') return placed;

        // Signature -> template name lookup (ring count + sorted sizes only;
        // the downstream matcher handles atom correspondence, so atom count
        // need not be exact here).
        var TEMPLATE_LOOKUP = {
            // Expanded scaffold coverage. Order matters: chromone first in
            // 2:6,6 so the chromoneRingMap special-case claims any
            // chroman-shaped molecule before the generic bicyclic templates
            // probe it.
            '1:3':       ['cyclopropane'],
            '1:4':       ['cyclobutane', 'beta_lactam'],
            '1:5':       ['cyclopentane', 'pyrrolidine', 'furanose', 'imidazole', 'pyrazole', 'thiophene', 'furan', 'oxazole', 'isoxazole', 'thiazole', 'pyrrole', 'triazole_1_2_3', 'triazole_1_2_4', 'tetrazole'],
            '1:6':       ['benzene', 'cyclohexane', 'tetrahydropyran', 'piperidine', 'pyridine', 'pyrimidine', 'piperazine', 'morpholine'],
            '1:7':       ['cycloheptane'],
            // 2:4,5 (penam, carbapenem) and 2:4,6 (cepham) intentionally
            // not in TEMPLATE_LOOKUP — the generic algorithm produces a
            // cleaner layout for penicillin G's acyl side chain than the
            // template stamper does. Templates remain available via
            // Templates.apply(mol, 'penam', cx, cy) for explicit use.
            '2:5,5':     ['norbornane', 'pyrrolizidine'],
            '2:5,6':     ['purine', 'indole', 'benzimidazole', 'benzofuran', 'benzothiazole', 'benzothiophene', 'indolizidine'],
            '2:5,7':     ['tropane'],
            '2:6,6':     ['chromone', 'naphthalene', 'quinoline', 'isoquinoline', 'quinazoline', 'quinoxaline', 'cinnoline', 'phthalazine', 'pteridine', 'quinolizidine'],
            '2:6,7':     ['benzodiazepine'],
            '3:5,6,6':   ['carbazole', 'beta_carboline', 'fluorene'],
            '3:6,6,6':   ['adamantane', 'phenanthrene', 'phenothiazine', 'acridine', 'anthracene'],
            '3:6,6,7':   ['dibenzazepine'],
            '4:5,6,6,6': ['steroid'],
            // v2.0.35: pyrene restored with correct peri-condensed geometry
            // (16 atoms, two atoms shared by 3 rings each). See
            // Templates.pyrene comment block.
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

            if (sig === '2:5,5') {
                var norMapping = norbornaneRingMap(mol, sysAtoms, sysRingIdx, rings);
                var norTmpl = (norMapping &&
                               typeof Templates !== 'undefined' &&
                               typeof Templates.norbornane === 'function')
                              ? Templates.norbornane() : null;
                if (norTmpl &&
                    validateTemplateBondLengths(mol, norMapping, norTmpl, sysAtoms)) {
                    applyTemplateCoordsVF2(mol, norMapping, norTmpl, placed);
                    continue;
                }
            }

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

            // Element-frequency screen + MCS greedy probe (O(N²), no backtracking).
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

            // Apply template coordinates using the resolved atom mapping.
            applyTemplateCoordsVF2(mol, mapping, tmpl, placed);
        }

        return placed;
    }

    /**
     * Fast template matching: element-frequency pre-screen + MCS greedy probe.
     * O(N) + O(N²) — no subgraph isomorphism backtracking, no NLF overhead, no hangs.
     * Returns {tmplIdx: molAtomId} or null.
     */
    function greedyMatchTemplate(mol, sysAtomIds, tmpl) {
        if (typeof SMSDMCS === 'undefined' || !SMSDMCS._greedyProbe) return null;
        if (typeof SMSDGraph === 'undefined') return null;

        try {
            // -- Stage 1: element-frequency screen (O(N)) --
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
     * Deterministic mapping for bornane / norbornane-style 2:5,5 cages.
     *
     * The SSSR often represents these as two 5-rings sharing a three-atom
     * bridge path. A greedy all-carbon match can choose the wrong bridgehead,
     * so map the shared path directly:
     *
     *   template 0 -- 6 -- 3 is the short bridge path
     *   template 0 -- 1 -- 2 -- 3 is one side path
     *   template 0 -- 5 -- 4 -- 3 is the other side path
     */
    function norbornaneRingMap(mol, sysAtoms, sysRingIdx, rings) {
        if (!sysRingIdx || sysRingIdx.length !== 2) return null;
        var rA = rings[sysRingIdx[0]];
        var rB = rings[sysRingIdx[1]];
        if (!rA || !rB || rA.length !== 5 || rB.length !== 5) return null;
        if (sysAtoms.length !== 7) return null;

        for (var si = 0; si < sysAtoms.length; si++) {
            var atom = mol.getAtom(sysAtoms[si]);
            if (!atom || atom.symbol !== 'C') return null;
        }

        var inB = {};
        for (var bi = 0; bi < rB.length; bi++) inB[rB[bi]] = true;
        var shared = [];
        for (var ai = 0; ai < rA.length; ai++) {
            if (inB[rA[ai]]) shared.push(rA[ai]);
        }
        if (shared.length !== 3) return null;

        var sharedSet = {};
        for (var sh = 0; sh < shared.length; sh++) sharedSet[shared[sh]] = true;

        var mid = -1;
        var ends = [];
        for (var si2 = 0; si2 < shared.length; si2++) {
            var sid = shared[si2];
            var nbrs = mol.getNeighbors(sid) || [];
            var sharedNbrCount = 0;
            for (var ni = 0; ni < nbrs.length; ni++) {
                if (sharedSet[nbrs[ni]]) sharedNbrCount++;
            }
            if (sharedNbrCount === 2) mid = sid;
            else if (sharedNbrCount === 1) ends.push(sid);
        }
        if (mid < 0 || ends.length !== 2) return null;

        var pathA = ringInteriorPath(rA, ends[0], ends[1], sharedSet);
        var pathB = ringInteriorPath(rB, ends[0], ends[1], sharedSet);
        if (!pathA || !pathB || pathA.length !== 2 || pathB.length !== 2) return null;

        var mapping = {};
        mapping[0] = ends[0];
        mapping[6] = mid;
        mapping[3] = ends[1];
        mapping[1] = pathA[0];
        mapping[2] = pathA[1];
        mapping[5] = pathB[0];
        mapping[4] = pathB[1];
        return mapping;
    }

    function ringInteriorPath(ring, startId, endId, blockedInterior) {
        var n = ring.length;
        var start = ring.indexOf(startId);
        var end = ring.indexOf(endId);
        if (start < 0 || end < 0) return null;

        function walk(step) {
            var out = [];
            var idx = start;
            var guard = 0;
            while (guard++ < n + 1) {
                idx = (idx + step + n) % n;
                if (ring[idx] === endId) return out;
                if (blockedInterior[ring[idx]]) return null;
                out.push(ring[idx]);
            }
            return null;
        }

        return walk(1) || walk(-1);
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
     * Apply template coordinates using subgraph-isomorphism-derived mapping.
     * mapping = {tmplIdx: molAtomId}
     */
    function applyTemplateCoordsVF2(mol, mapping, tmpl, placed) {
        var fitter = global.SDG && global.SDG.TemplateHandler &&
            global.SDG.TemplateHandler.applyTemplateBestFit;
        if (fitter) {
            var transform = fitter(mol, mapping, tmpl, placed, {
                bondLength: BOND_LENGTH,
                allowScale: false
            });
            if (transform) return;
        }

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

        var branchAngles = null;
        if (isRingParent && total >= 2 && numPlaced >= 2) {
            branchAngles = chooseCrossingFreeBranchFan(mol, parentId, unplacedIds,
                                                       outAngle, arcPerBranch,
                                                       placed);
        }

        // Sort branches: put ring-attached atoms to the outside, hydrogens inside
        // (simple heuristic: heavier substituents get the outermost positions)
        for (var i = 0; i < total; i++) {
            var angle = branchAngles ? branchAngles[i] : startAngle + i * arcPerBranch;
            var child = mol.getAtom(unplacedIds[i]);
            child.x = parent.x + BOND_LENGTH * Math.cos(angle);
            child.y = parent.y + BOND_LENGTH * Math.sin(angle);
            placed[unplacedIds[i]] = true;
        }
    }

    function chooseCrossingFreeBranchFan(mol, parentId, unplacedIds,
                                         preferredAngle, arcPerBranch, placed) {
        var parent = mol.getAtom(parentId);
        if (!parent) return null;
        var total = unplacedIds.length;
        var totalArc = arcPerBranch * (total - 1);
        var placedIds = Object.keys(placed).map(function (id) {
            return parseInt(id, 10);
        }).filter(function (id) {
            return id !== parentId && !!mol.getAtom(id);
        });

        var bestAngles = null;
        var bestScore = -Infinity;
        var arcCandidates = [arcPerBranch];
        if (total === 2) {
            arcCandidates.push(DEG120, DEG120 * 0.85, Math.PI / 2);
        }

        for (var ac = 0; ac < arcCandidates.length; ac++) {
            var candidateArc = arcCandidates[ac];
            totalArc = candidateArc * (total - 1);
            for (var step = -12; step <= 12; step++) {
                var center = preferredAngle + step * (Math.PI / 12);
                var start = center - totalArc / 2;
                var angles = [];
                var centerDelta = normalizeAngle(center - preferredAngle);
                centerDelta = Math.min(centerDelta, TWO_PI - centerDelta);
                var score = -centerDelta * 0.05;
                score -= Math.abs(candidateArc - arcPerBranch) * 20;
                var crossings = 0;
                var minPairDistSq = Infinity;

                for (var i = 0; i < total; i++) {
                    var angle = start + i * candidateArc;
                    angles.push(angle);
                    var x = parent.x + BOND_LENGTH * Math.cos(angle);
                    var y = parent.y + BOND_LENGTH * Math.sin(angle);
                    var minAtomDistSq = Infinity;
                    for (var pi = 0; pi < placedIds.length; pi++) {
                        var pa = mol.getAtom(placedIds[pi]);
                        var d2 = distSq(x, y, pa.x, pa.y);
                        if (d2 < minAtomDistSq) minAtomDistSq = d2;
                    }
                    if (minAtomDistSq < Infinity) {
                        score += Math.min(minAtomDistSq, BOND_LENGTH * BOND_LENGTH * 4) /
                                 (BOND_LENGTH * BOND_LENGTH);
                    }

                    for (var bi = 0; bi < mol.bonds.length; bi++) {
                        var b = mol.bonds[bi];
                        if (b.atom1 === parentId || b.atom2 === parentId) continue;
                        if (!placed[b.atom1] || !placed[b.atom2]) continue;
                        var a1 = mol.getAtom(b.atom1);
                        var a2 = mol.getAtom(b.atom2);
                        if (!a1 || !a2) continue;
                        if (segmentsCrossStrict(parent.x, parent.y, x, y,
                                                a1.x, a1.y, a2.x, a2.y)) {
                            crossings++;
                        }
                    }

                    for (var j = 0; j < i; j++) {
                        var prevAngle = angles[j];
                        var px = parent.x + BOND_LENGTH * Math.cos(prevAngle);
                        var py = parent.y + BOND_LENGTH * Math.sin(prevAngle);
                        var pairD2 = distSq(x, y, px, py);
                        if (pairD2 < minPairDistSq) minPairDistSq = pairD2;
                    }
                }

                if (minPairDistSq < Infinity) {
                    score += Math.min(minPairDistSq, BOND_LENGTH * BOND_LENGTH * 4) /
                             (BOND_LENGTH * BOND_LENGTH);
                }
                score -= crossings * 10000;
                if (score > bestScore) {
                    bestScore = score;
                    bestAngles = angles;
                }
            }
        }
        return bestAngles;
    }

    function segmentsCrossStrict(ax, ay, bx, by, cx, cy, dx, dy) {
        function orient(px, py, qx, qy, rx, ry) {
            return (qx - px) * (ry - py) - (qy - py) * (rx - px);
        }
        var o1 = orient(ax, ay, bx, by, cx, cy);
        var o2 = orient(ax, ay, bx, by, dx, dy);
        var o3 = orient(cx, cy, dx, dy, ax, ay);
        var o4 = orient(cx, cy, dx, dy, bx, by);
        if (Math.abs(o1) < EPSILON || Math.abs(o2) < EPSILON ||
            Math.abs(o3) < EPSILON || Math.abs(o4) < EPSILON) {
            return false;
        }
        return ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
               ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0));
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
    /**
     * resolveCollisions — spatial-hash summed-push overlap resolver.
     *
     * v1.8.25 consolidation: this is now the SINGLE collision resolver in
     * Layout.js. It absorbed the former `resolveChainOverlaps15` (a v1.8.16
     * near-duplicate) via three `opts` knobs — the two functions shared
     * ~80 % of their body (grid build, 3×3 neighbourhood scan, push law).
     *
     * opts (all optional; defaults reproduce the pre-consolidation
     * `resolveCollisions` behaviour byte-for-byte):
     *   ringGentleFactor : push scale for ring-ring pairs. Default 0.3
     *                      (bridged-ring relief). Step 15 passes 0.15.
     *   atomToSys        : atom-id → ring-system-index map. When provided,
     *                      ring-ring pairs in DIFFERENT systems are skipped
     *                      entirely (Step 15 uses this — inter-system
     *                      overlaps are handled by the rigid separator, and
     *                      pushing them here would un-snap the polygons).
     *                      Default null → all ring-ring pairs get the
     *                      gentle push (original behaviour).
     *   snapshotRings    : when true (default) ring atom positions are
     *                      snapshotted and restored at the end unless they
     *                      were actually pushed — keeps fused-ring geometry
     *                      stable across substituent crowding. Step 15
     *                      passes false (rings were just snapped to perfect
     *                      polygons; the chain-only sweep must not touch
     *                      them at all).
     */
    function resolveCollisions(mol, atomIds, ringAtomSet, opts) {
        var n = atomIds.length;
        if (n < 2) return;
        opts = opts || {};
        var ringGentleFactor = (opts.ringGentleFactor !== undefined)
            ? opts.ringGentleFactor : 0.3;
        var atomToSys     = opts.atomToSys || null;
        var snapshotRings = (opts.snapshotRings !== false);

        // Snapshot ring atom positions so we can restore them after overlap
        // resolution. This prevents cumulative distortion of ring geometry
        // when many substituents are attached (e.g., steroids with 4 fused
        // rings and numerous substituents).
        var ringSnapshot = {};
        if (snapshotRings) {
            for (var i = 0; i < n; i++) {
                if (ringAtomSet[atomIds[i]]) {
                    var ra = mol.getAtom(atomIds[i]);
                    if (ra) ringSnapshot[atomIds[i]] = { x: ra.x, y: ra.y };
                }
            }
        }

        for (var iter = 0; iter < MAX_OVERLAP_ITER; iter++) {
            var collision = false;

            // Build spatial hash
            var grid = {};
            for (var i = 0; i < n; i++) {
                var a  = mol.getAtom(atomIds[i]);
                if (!a) continue;
                var gx = Math.floor(a.x / GRID_CELL);
                var gy = Math.floor(a.y / GRID_CELL);
                var key = gx + ',' + gy;
                if (!grid[key]) grid[key] = [];
                grid[key].push(i);
            }

            // Check each atom against neighbours in the 3x3 grid neighbourhood
            for (var i = 0; i < n; i++) {
                var a1 = mol.getAtom(atomIds[i]);
                if (!a1) continue;
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
                            if (!a2) continue;
                            var ddx = a2.x - a1.x, ddy = a2.y - a1.y;
                            var d = Math.sqrt(ddx * ddx + ddy * ddy);
                            if (d >= MIN_ATOM_DIST) continue;

                            var r1 = !!ringAtomSet[atomIds[i]];
                            var r2 = !!ringAtomSet[atomIds[j]];

                            // Inter-ring-system ring-ring pairs are skipped
                            // when atomToSys is supplied (Step 15 path).
                            if (r1 && r2 && atomToSys &&
                                atomToSys[atomIds[i]] !== atomToSys[atomIds[j]]) {
                                continue;
                            }

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

                            if (r1 && r2) {
                                // Both ring atoms: gentle push to resolve
                                // bridged-ring overlaps (e.g. morphine oxygen
                                // bridge, norbornane) without large distortion.
                                var gentlePush = push * ringGentleFactor;
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

            if (!collision) break;
        }

        // Restore ring atom positions unless they were pushed to resolve
        // a ring-ring overlap (bridged systems like morphine, norbornane).
        if (snapshotRings) {
            for (var id in ringSnapshot) {
                var atom = mol.getAtom(parseInt(id));
                if (!atom) continue;
                var sx = ringSnapshot[id].x, sy = ringSnapshot[id].y;
                var dx2 = atom.x - sx, dy2 = atom.y - sy;
                var wasMoved = (dx2 * dx2 + dy2 * dy2) > 0.01;
                if (!wasMoved) {
                    atom.x = sx;
                    atom.y = sy;
                }
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

    function snapshotAtomPositions(mol, atomIds) {
        var snap = [];
        for (var i = 0; i < atomIds.length; i++) {
            var a = mol.getAtom(atomIds[i]);
            if (a) snap.push({ id: atomIds[i], x: a.x, y: a.y });
        }
        return snap;
    }

    function restoreAtomPositions(mol, snap) {
        if (!snap) return;
        for (var i = 0; i < snap.length; i++) {
            var a = mol.getAtom(snap[i].id);
            if (a) { a.x = snap[i].x; a.y = snap[i].y; }
        }
    }

    function runSDGRefinement(mol, atomIds, ringAtomSet, sdgLayout, guarded) {
        if (!sdgLayout || !sdgLayout.refineLayout) return false;
        var snap = snapshotAtomPositions(mol, atomIds);
        var preScore = guarded ? layoutQualityPenalty(mol, atomIds) : 0;
        try {
            sdgLayout.refineLayout(mol, atomIds, {
                bondLength: BOND_LENGTH,
                minDist: BOND_LENGTH * 0.6,
                maxIters: 60,
                stepSize: 0.5
            });
            // Restore bond lengths after the SDG pass — the overlap resolver
            // may have stretched some bonds while resolving atom clashes.
            relaxBondLengths(mol, atomIds, ringAtomSet);
            if (guarded && layoutQualityPenalty(mol, atomIds) > preScore + 0.05) {
                restoreAtomPositions(mol, snap);
                return false;
            }
            return true;
        } catch (e) {
            restoreAtomPositions(mol, snap);
            return false;
        }
    }

    function layoutQualityMetrics(mol, atomIds) {
        var quality = global.SDG && global.SDG.LayoutQuality;
        if (quality && quality.evaluate) {
            return quality.evaluate(mol, atomIds, { bondLength: BOND_LENGTH });
        }

        var idSet = {};
        for (var i = 0; i < atomIds.length; i++) idSet[atomIds[i]] = true;

        var penalty = 0;
        var needsAdaptiveRescue = false;
        var bonded = {};
        for (var bi = 0; bi < mol.bonds.length; bi++) {
            var bond = mol.bonds[bi];
            bonded[bond.atom1 + ',' + bond.atom2] = true;
            bonded[bond.atom2 + ',' + bond.atom1] = true;
            if (!idSet[bond.atom1] || !idSet[bond.atom2]) continue;
            var a = mol.getAtom(bond.atom1), b = mol.getAtom(bond.atom2);
            if (!a || !b) continue;
            var d = dist(a.x, a.y, b.x, b.y);
            var ratio = d / BOND_LENGTH;
            if (ratio > 1.45 || ratio < 0.55) needsAdaptiveRescue = true;
            if (ratio > 1.25) penalty += (ratio - 1.25) * 10;
            if (ratio < 0.65) penalty += (0.65 - ratio) * 10;
        }

        var rescueClose = BOND_LENGTH * 0.35;
        var penaltyClose = BOND_LENGTH * 0.45;
        for (var ai = 0; ai < atomIds.length; ai++) {
            var aa = mol.getAtom(atomIds[ai]);
            if (!aa) continue;
            for (var aj = ai + 1; aj < atomIds.length; aj++) {
                if (bonded[atomIds[ai] + ',' + atomIds[aj]]) continue;
                var ab = mol.getAtom(atomIds[aj]);
                if (!ab) continue;
                var nd = dist(aa.x, aa.y, ab.x, ab.y);
                if (nd < rescueClose) needsAdaptiveRescue = true;
                if (nd < penaltyClose) {
                    penalty += 2 + (penaltyClose - nd) / BOND_LENGTH * 20;
                }
            }
        }
        return { penalty: penalty, needsAdaptiveRescue: needsAdaptiveRescue };
    }

    function layoutNeedsAdaptiveRescue(mol, atomIds) {
        return layoutQualityMetrics(mol, atomIds).needsAdaptiveRescue;
    }

    function layoutQualityPenalty(mol, atomIds) {
        return layoutQualityMetrics(mol, atomIds).penalty;
    }

    function maxBondLengthForAtomSet(mol, atomIds) {
        var idSet = {};
        for (var i = 0; i < atomIds.length; i++) idSet[atomIds[i]] = true;
        var maxLen = 0;
        for (var bi = 0; bi < mol.bonds.length; bi++) {
            var bond = mol.bonds[bi];
            if (!idSet[bond.atom1] || !idSet[bond.atom2]) continue;
            var a = mol.getAtom(bond.atom1);
            var b = mol.getAtom(bond.atom2);
            if (!a || !b) continue;
            var len = dist(a.x, a.y, b.x, b.y);
            if (len > maxLen) maxLen = len;
        }
        return maxLen;
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

        // v2.0.34 fix (was: per-ring iteration with last-write-wins on
        // bridgehead atoms): build an atomId → [ringIndex…] index so
        // bridgehead atoms (members of N ≥ 2 rings) get visited ONCE and
        // their substituent is rotated along the BISECTOR of all containing
        // rings' outward radials, not whichever ring iterated last.
        var ringsForAtom = {};
        for (var ri2 = 0; ri2 < rings.length; ri2++) {
            var ring2 = rings[ri2];
            for (var k2 = 0; k2 < ring2.length; k2++) {
                var aid2 = ring2[k2];
                if (!ringsForAtom[aid2]) { ringsForAtom[aid2] = []; }
                ringsForAtom[aid2].push(ri2);
            }
        }

        // Iterate unique ring atom ids in sorted order (deterministic).
        var uniqueIds = Object.keys(ringsForAtom).map(function (s) { return parseInt(s, 10); });
        uniqueIds.sort(function (x, y) { return x - y; });

        for (var ui = 0; ui < uniqueIds.length; ui++) {
            var aid = uniqueIds[ui];
            var atom = mol.getAtom(aid);
            var nb = mol.getNeighbors(aid);
            var subIds = [];
            for (var ni = 0; ni < nb.length; ni++) {
                if (!ringAtomSet[nb[ni]]) subIds.push(nb[ni]);
            }
            if (subIds.length === 0) continue;
            subIds.sort(function (x, y) { return x - y; });

            // Radial direction: BISECTOR of the outward radials from every
            // ring the atom belongs to. For a single-ring atom this is the
            // original behaviour. For a bridgehead, this is the
            // chemically-expected exterior-angle bisector (sum of unit
            // radials → normalise → atan2).
            var ringIdxs = ringsForAtom[aid];
            var sumX = 0, sumY = 0;
            for (var rj = 0; rj < ringIdxs.length; rj++) {
                var rcj = ringCentroids[ringIdxs[rj]];
                var dx = atom.x - rcj.x, dy = atom.y - rcj.y;
                var dl = Math.sqrt(dx * dx + dy * dy);
                if (dl < EPSILON) continue;
                sumX += dx / dl;
                sumY += dy / dl;
            }
            var rd = Math.sqrt(sumX * sumX + sumY * sumY);
            if (rd < EPSILON) continue;
            var rdx = sumX / rd, rdy = sumY / rd;
            var radialAngle = Math.atan2(rdy, rdx);

                // If single sub: place along radial direction.
                // If multiple: fan symmetrically around radial direction.
                var totalSubs = subIds.length;
                var fan = (totalSubs === 1) ? 0 : DEG120 * 0.5; // +/-60 deg
                for (var si = 0; si < totalSubs; si++) {
                    var subId = subIds[si];
                    var subAtom = mol.getAtom(subId);
                    if (!subAtom) continue;
                    if (subtreeTouchesOtherRingThroughSingleBonds(mol, subId, aid, ringAtomSet)) {
                        // This substituent is a linker to another ring system.
                        // Rotating only a flexible single-bond chain while
                        // leaving the far ring fixed stretches bonds badly, so
                        // preserve the base placement. Conjugated linkers such
                        // as azo bridges still get the radial convention pass.
                        continue;
                    }
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

    function subtreeTouchesOtherRingThroughSingleBonds(mol, rootId, parentId, ringAtomSet) {
        var visited = {};
        visited[parentId + ':0'] = true;
        visited[parentId + ':1'] = true;
        var stack = [{ id: rootId, hasMultipleBond: false }];
        while (stack.length) {
            var rec = stack.pop();
            var cur = rec.id;
            var key = cur + ':' + (rec.hasMultipleBond ? '1' : '0');
            if (visited[key]) continue;
            visited[key] = true;
            if (ringAtomSet[cur]) {
                if (!rec.hasMultipleBond) return true;
                continue;
            }
            var nb = mol.getNeighbors(cur) || [];
            for (var i = 0; i < nb.length; i++) {
                var bond = mol.getBondBetween ? mol.getBondBetween(cur, nb[i]) : null;
                stack.push({
                    id: nb[i],
                    hasMultipleBond: rec.hasMultipleBond || !!(bond && bond.type > 1)
                });
            }
        }
        return false;
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

    // v1.8.26: the second byte-identical copy of the "Utility functions"
    // + "Ring Query Helpers" block (dist / distSq / normalizeAngle /
    // atomBounds / uniqueArray / Layout.getRingInfo / Layout.bondInRing)
    // was removed here — it was a verbatim duplicate of the block ~570
    // lines above (search "// Utility functions"). JS function-hoisting
    // and last-assignment-wins meant this copy was always shadowed by
    // the earlier one; ~118 lines of pure copy-paste gone. The
    // duplicated "// Export" header was de-duplicated at the same time.

    // =====================================================================
    // Export
    // =====================================================================

    global.Layout = Layout;

})(typeof window !== 'undefined' ? window : this);
