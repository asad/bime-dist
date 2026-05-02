/**
 * SMSDLayout.js — Layout support functions for 2D coordinate generation
 *
 * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK.
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * Ported from SMSD upstream, layout.hpp / ring_finder.hpp
 *
 * Provides:
 *   - computeSSSR(graph)     — minimum cycle basis sorted ascending by size
 *   - layoutSSSR(graph)      — rings ordered for 2D layout (largest system first,
 *                               fused rings by shared-edge adjacency)
 *   - reduceCrossings(mol, atomIds, ringAtomSet, rings)
 *                             — simulated annealing on ring system orientations
 *                               to minimise bond crossings
 *   - forceDirectedLayout(mol, atomIds, ringAtomSet, rings, maxIter)
 *                             — PrEd-inspired force-directed refinement that
 *                               preserves edge crossing zones (Bertault 2000)
 *   - smacofLayout(mol, atomIds, maxIter, nInit)
 *                             — SMACOF stress majorisation (de Leeuw) with
 *                               multi-init random starts
 *   - SCAFFOLD_TEMPLATES      — 45 canonical scaffold coordinate arrays
 *
 * Dependencies: window.Molecule (for BOND_LENGTH constant)
 */
(function(global) {
    'use strict';

    // =====================================================================
    // Point2D helper
    // =====================================================================

    function Point2D(x, y) {
        this.x = x || 0;
        this.y = y || 0;
    }

    Point2D.prototype.distTo = function(other) {
        var dx = this.x - other.x;
        var dy = this.y - other.y;
        return Math.sqrt(dx * dx + dy * dy);
    };

    Point2D.prototype.add = function(other) {
        return new Point2D(this.x + other.x, this.y + other.y);
    };

    Point2D.prototype.sub = function(other) {
        return new Point2D(this.x - other.x, this.y - other.y);
    };

    Point2D.prototype.scale = function(s) {
        return new Point2D(this.x * s, this.y * s);
    };

    // =====================================================================
    // Segment intersection test
    // =====================================================================

    /**
     * Test whether segment (p1,p2) intersects segment (p3,p4).
     * Returns true if the segments properly intersect (not counting shared
     * endpoints or collinear overlaps).
     */
    function segmentsIntersect(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y) {
        var d1x = p2x - p1x, d1y = p2y - p1y;
        var d2x = p4x - p3x, d2y = p4y - p3y;
        var denom = d1x * d2y - d1y * d2x;
        if (Math.abs(denom) < 1e-12) { return false; } // parallel
        var dx = p3x - p1x, dy = p3y - p1y;
        var t = (dx * d2y - dy * d2x) / denom;
        var u = (dx * d1y - dy * d1x) / denom;
        // Proper intersection: exclude endpoints to avoid counting shared atoms
        var eps = 0.001;
        return t > eps && t < (1 - eps) && u > eps && u < (1 - eps);
    }

    // =====================================================================
    // computeSSSR — Minimum Cycle Basis via Horton + GF(2) elimination
    // =====================================================================

    /**
     * Compute the Smallest Set of Smallest Rings (SSSR / Minimum Cycle Basis)
     * for a graph described by neighbor lists.
     *
     * @param {Object} graph — { n: number, neighbors: number[][] }
     *        graph.n            = number of atoms
     *        graph.neighbors[i] = array of neighbor indices for atom i
     * @returns {Array<Array<number>>} — rings sorted ascending by size,
     *          each ring is an ordered array of atom indices
     */
    function computeSSSR(graph) {
        var n = graph.n;
        if (n < 3) { return []; }

        var neighbors = graph.neighbors;

        // Count edges
        var edgeCount = 0;
        var edgeIndex = {}; // 'u,v' -> index  (u < v)
        var edgeList = [];  // [[u,v], ...]
        var i, j, u, v, nbs;
        for (i = 0; i < n; i++) {
            nbs = neighbors[i];
            for (j = 0; j < nbs.length; j++) {
                v = nbs[j];
                if (v > i) {
                    edgeIndex[i + ',' + v] = edgeCount;
                    edgeList.push([i, v]);
                    edgeCount++;
                }
            }
        }

        var cycleRank = edgeCount - n;
        // Count connected components to adjust cycle rank
        var visited = new Array(n);
        for (i = 0; i < n; i++) { visited[i] = false; }
        var numComponents = 0;
        for (i = 0; i < n; i++) {
            if (visited[i]) { continue; }
            numComponents++;
            var stack = [i];
            visited[i] = true;
            while (stack.length > 0) {
                u = stack.pop();
                nbs = neighbors[u];
                for (j = 0; j < nbs.length; j++) {
                    if (!visited[nbs[j]]) {
                        visited[nbs[j]] = true;
                        stack.push(nbs[j]);
                    }
                }
            }
        }
        cycleRank += numComponents;
        if (cycleRank <= 0) { return []; }

        // BFS shortest-path trees from every atom
        var dist = new Array(n);
        var parent = new Array(n);
        for (i = 0; i < n; i++) {
            dist[i] = new Array(n);
            parent[i] = new Array(n);
            for (j = 0; j < n; j++) {
                dist[i][j] = -1;
                parent[i][j] = -1;
            }
        }

        for (var src = 0; src < n; src++) {
            dist[src][src] = 0;
            var queue = [src];
            var head = 0;
            while (head < queue.length) {
                u = queue[head++];
                nbs = neighbors[u];
                for (j = 0; j < nbs.length; j++) {
                    v = nbs[j];
                    if (dist[src][v] < 0) {
                        dist[src][v] = dist[src][u] + 1;
                        parent[src][v] = u;
                        queue.push(v);
                    }
                }
            }
        }

        // Horton candidate generation: for each edge (u,v) and each vertex w,
        // if dist[w][u] + 1 + dist[w][v] == ring size, we have a candidate
        var candidates = [];
        var words = Math.ceil(edgeCount / 32);

        for (var ei = 0; ei < edgeCount; ei++) {
            var eu = edgeList[ei][0], ev = edgeList[ei][1];
            for (var w = 0; w < n; w++) {
                if (dist[w][eu] < 0 || dist[w][ev] < 0) { continue; }
                var ringSize = dist[w][eu] + dist[w][ev] + 1;
                if (ringSize < 3) { continue; }

                // Trace paths from w to eu and from w to ev
                var pathU = tracePath(parent, w, eu);
                var pathV = tracePath(parent, w, ev);

                // Check that paths don't share internal vertices
                if (pathU === null || pathV === null) { continue; }
                var overlap = false;
                var inPathU = {};
                for (j = 1; j < pathU.length; j++) { inPathU[pathU[j]] = true; }
                for (j = 1; j < pathV.length; j++) {
                    if (inPathU[pathV[j]]) { overlap = true; break; }
                }
                if (overlap) { continue; }

                // Build ring: pathU + reverse(pathV[1..])
                var ring = pathU.slice();
                for (j = pathV.length - 1; j >= 1; j--) {
                    ring.push(pathV[j]);
                }

                // Build edge vector (GF(2) bit vector over edges)
                var vec = new Array(words);
                for (j = 0; j < words; j++) { vec[j] = 0; }
                for (j = 0; j < ring.length; j++) {
                    var a = ring[j], b = ring[(j + 1) % ring.length];
                    var lo = a < b ? a : b, hi = a < b ? b : a;
                    var eidx = edgeIndex[lo + ',' + hi];
                    if (eidx !== undefined) {
                        vec[eidx >>> 5] |= (1 << (eidx & 31));
                    }
                }

                candidates.push({
                    ring: ring,
                    size: ring.length,
                    vec: vec
                });
            }
        }

        // Sort candidates by size (ascending)
        candidates.sort(function(a, b) { return a.size - b.size; });

        // Deduplicate by edge vector
        var uniqueCands = [];
        var seenVecs = {};
        for (i = 0; i < candidates.length; i++) {
            var key = candidates[i].vec.join(',');
            if (!seenVecs[key]) {
                seenVecs[key] = true;
                uniqueCands.push(candidates[i]);
            }
        }

        // GF(2) Gaussian elimination to find linearly independent set
        var basis = [];
        var basisVecs = [];

        for (i = 0; i < uniqueCands.length && basis.length < cycleRank; i++) {
            var cv = uniqueCands[i].vec.slice();

            // Reduce against existing basis
            for (j = 0; j < basisVecs.length; j++) {
                var pivot = findPivot(basisVecs[j], words);
                if (pivot >= 0 && (cv[pivot >>> 5] & (1 << (pivot & 31)))) {
                    for (var wi = 0; wi < words; wi++) {
                        cv[wi] ^= basisVecs[j][wi];
                    }
                }
            }

            // Check if still non-zero (linearly independent)
            var nonzero = false;
            for (j = 0; j < words; j++) {
                if (cv[j] !== 0) { nonzero = true; break; }
            }
            if (nonzero) {
                basisVecs.push(cv);
                basis.push(uniqueCands[i].ring);
            }
        }

        // Sort result ascending by ring size
        basis.sort(function(a, b) { return a.length - b.length; });
        return basis;
    }

    /** Trace BFS parent pointers from src to dst. Returns path [src,...,dst] or null. */
    function tracePath(parentMatrix, src, dst) {
        if (src === dst) { return [src]; }
        var path = [];
        var cur = dst;
        var maxSteps = parentMatrix.length + 1;
        while (cur !== src && maxSteps-- > 0) {
            path.push(cur);
            cur = parentMatrix[src][cur];
            if (cur < 0) { return null; }
        }
        if (cur !== src) { return null; }
        path.push(src);
        path.reverse();
        return path;
    }

    /** Find lowest set bit position in a GF(2) vector. */
    function findPivot(vec, words) {
        for (var w = 0; w < words; w++) {
            if (vec[w] !== 0) {
                // Count trailing zeros
                var x = vec[w];
                var bit = 0;
                while ((x & 1) === 0) { x >>>= 1; bit++; }
                return (w << 5) | bit;
            }
        }
        return -1;
    }

    // =====================================================================
    // layoutSSSR — rings ordered for 2D layout
    // =====================================================================

    /**
     * Order SSSR rings for 2D coordinate generation:
     *   1. Group into fused ring systems (two rings sharing >= 2 atoms)
     *   2. Sort systems: largest system first
     *   3. Within each system: BFS from the largest ring outward via
     *      shared-edge adjacency
     *
     * @param {Object} graph — same format as computeSSSR input
     * @returns {Array<Array<number>>} — rings ordered for layout
     */
    function layoutSSSR(graph) {
        var rings = computeSSSR(graph);
        if (rings.length <= 1) { return rings; }

        var nRings = rings.length;

        // Build ring adjacency (shared >= 2 atoms = shared edge)
        var ringAdj = new Array(nRings);
        var i, j, k;
        for (i = 0; i < nRings; i++) { ringAdj[i] = []; }

        for (i = 0; i < nRings; i++) {
            var setI = {};
            for (k = 0; k < rings[i].length; k++) { setI[rings[i][k]] = true; }
            for (j = i + 1; j < nRings; j++) {
                var shared = 0;
                for (k = 0; k < rings[j].length; k++) {
                    if (setI[rings[j][k]]) { shared++; }
                }
                if (shared >= 2) {
                    ringAdj[i].push(j);
                    ringAdj[j].push(i);
                }
            }
        }

        // Union-find to group into systems
        var parent = new Array(nRings);
        for (i = 0; i < nRings; i++) { parent[i] = i; }
        function find(x) {
            while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        }
        function unite(a, b) {
            a = find(a); b = find(b);
            if (a !== b) { parent[a] = b; }
        }

        for (i = 0; i < nRings; i++) {
            for (j = 0; j < ringAdj[i].length; j++) {
                unite(i, ringAdj[i][j]);
            }
        }

        // Group rings by system
        var systemMap = {};
        for (i = 0; i < nRings; i++) {
            var root = find(i);
            if (!systemMap[root]) { systemMap[root] = []; }
            systemMap[root].push(i);
        }

        // Collect systems, sort by total atom count descending
        var systems = [];
        for (k in systemMap) {
            if (systemMap.hasOwnProperty(k)) {
                systems.push(systemMap[k]);
            }
        }
        systems.sort(function(a, b) {
            var totalA = 0, totalB = 0;
            for (var si = 0; si < a.length; si++) { totalA += rings[a[si]].length; }
            for (var si2 = 0; si2 < b.length; si2++) { totalB += rings[b[si2]].length; }
            return totalB - totalA;
        });

        // BFS within each system from the largest ring
        var ordered = [];
        for (var si = 0; si < systems.length; si++) {
            var sys = systems[si];
            if (sys.length === 1) {
                ordered.push(rings[sys[0]]);
                continue;
            }

            // Find largest ring in system
            var bestIdx = sys[0], bestLen = rings[sys[0]].length;
            for (j = 1; j < sys.length; j++) {
                if (rings[sys[j]].length > bestLen) {
                    bestLen = rings[sys[j]].length;
                    bestIdx = sys[j];
                }
            }

            // BFS from largest ring
            var visited = {};
            var queue = [bestIdx];
            visited[bestIdx] = true;
            var head = 0;
            while (head < queue.length) {
                var cur = queue[head++];
                ordered.push(rings[cur]);
                var adj = ringAdj[cur];
                for (j = 0; j < adj.length; j++) {
                    if (!visited[adj[j]]) {
                        visited[adj[j]] = true;
                        queue.push(adj[j]);
                    }
                }
            }
        }

        return ordered;
    }

    // =====================================================================
    // reduceCrossings — Simulated annealing on ring system orientations
    // =====================================================================

    /**
     * Post-processing step: reduce bond crossings using simulated annealing
     * on ring system orientations. Each ring system can be reflected (flipped)
     * about its principal axis.
     *
     * @param {Molecule} mol       — molecule with existing 2D coordinates
     * @param {Array<number>} atomIds — atom IDs in this component
     * @param {Object} ringAtomSet — { atomId: true } for ring atoms
     * @param {Array<Array<number>>} rings — SSSR rings as atom ID arrays
     * @param {number} [maxIter]   — SA iterations (default 200)
     */
    function reduceCrossings(mol, atomIds, ringAtomSet, rings, maxIter) {
        if (!rings || rings.length < 2) { return; }
        if (!atomIds || atomIds.length < 4) { return; }
        if (maxIter === undefined) { maxIter = 200; }

        // Build ring systems (union-find on shared edges)
        var nRings = rings.length;
        var rParent = new Array(nRings);
        var i, j, k;
        for (i = 0; i < nRings; i++) { rParent[i] = i; }
        function rFind(x) {
            while (rParent[x] !== x) { rParent[x] = rParent[rParent[x]]; x = rParent[x]; }
            return x;
        }
        function rUnite(a, b) {
            a = rFind(a); b = rFind(b);
            if (a !== b) { rParent[a] = b; }
        }

        for (i = 0; i < nRings; i++) {
            var setI = {};
            for (k = 0; k < rings[i].length; k++) { setI[rings[i][k]] = true; }
            for (j = i + 1; j < nRings; j++) {
                var shared = 0;
                for (k = 0; k < rings[j].length; k++) {
                    if (setI[rings[j][k]]) { shared++; }
                }
                if (shared >= 2) { rUnite(i, j); }
            }
        }

        // Collect ring systems with their atom sets
        var sysMap = {};
        for (i = 0; i < nRings; i++) {
            var root = rFind(i);
            if (!sysMap[root]) { sysMap[root] = {}; }
            for (k = 0; k < rings[i].length; k++) {
                sysMap[root][rings[i][k]] = true;
            }
        }

        var systemList = [];
        for (k in sysMap) {
            if (sysMap.hasOwnProperty(k)) {
                var atoms = [];
                for (var aid in sysMap[k]) {
                    if (sysMap[k].hasOwnProperty(aid)) { atoms.push(+aid); }
                }
                if (atoms.length >= 3) { systemList.push(atoms); }
            }
        }

        if (systemList.length < 1) { return; }

        // Collect all bonds in this component for crossing checks
        var bonds = [];
        var inComp = {};
        for (i = 0; i < atomIds.length; i++) { inComp[atomIds[i]] = true; }
        for (i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            if (inComp[b.atom1] && inComp[b.atom2]) {
                bonds.push(b);
            }
        }

        if (bonds.length < 2) { return; }

        /** Count total bond crossings */
        function countCrossings() {
            var count = 0;
            for (var bi = 0; bi < bonds.length; bi++) {
                var b1 = bonds[bi];
                var a1 = mol.getAtom(b1.atom1), a2 = mol.getAtom(b1.atom2);
                for (var bj = bi + 1; bj < bonds.length; bj++) {
                    var b2 = bonds[bj];
                    // Skip if bonds share an atom
                    if (b1.atom1 === b2.atom1 || b1.atom1 === b2.atom2 ||
                        b1.atom2 === b2.atom1 || b1.atom2 === b2.atom2) {
                        continue;
                    }
                    var a3 = mol.getAtom(b2.atom1), a4 = mol.getAtom(b2.atom2);
                    if (segmentsIntersect(a1.x, a1.y, a2.x, a2.y,
                                          a3.x, a3.y, a4.x, a4.y)) {
                        count++;
                    }
                }
            }
            return count;
        }

        /** Flip a ring system: reflect all atoms about the system centroid's axis */
        function flipSystem(sysAtoms) {
            // Compute centroid
            var cx = 0, cy = 0;
            for (var si = 0; si < sysAtoms.length; si++) {
                var atom = mol.getAtom(sysAtoms[si]);
                cx += atom.x;
                cy += atom.y;
            }
            cx /= sysAtoms.length;
            cy /= sysAtoms.length;

            // Find principal axis via covariance
            var cxx = 0, cxy = 0, cyy = 0;
            for (var si2 = 0; si2 < sysAtoms.length; si2++) {
                var a = mol.getAtom(sysAtoms[si2]);
                var dx = a.x - cx, dy = a.y - cy;
                cxx += dx * dx;
                cxy += dx * dy;
                cyy += dy * dy;
            }

            // Principal axis angle
            var theta = 0.5 * Math.atan2(2 * cxy, cxx - cyy);

            // Reflect about this axis (through centroid)
            var cosT = Math.cos(2 * theta);
            var sinT = Math.sin(2 * theta);
            for (var si3 = 0; si3 < sysAtoms.length; si3++) {
                var atom2 = mol.getAtom(sysAtoms[si3]);
                var dx2 = atom2.x - cx, dy2 = atom2.y - cy;
                // Reflection matrix about axis at angle theta:
                // [cos2t  sin2t] [dx]
                // [sin2t -cos2t] [dy]
                atom2.x = cx + cosT * dx2 + sinT * dy2;
                atom2.y = cy + sinT * dx2 - cosT * dy2;
            }
        }

        /** Reverse a flip (same operation applied again = identity) */
        function unflipSystem(sysAtoms) {
            flipSystem(sysAtoms); // reflection is self-inverse
        }

        // Phase 1 candidates: whole ring systems
        var phase1Candidates = systemList.slice();

        // Phase 2 candidates: individual rings with exclusive/shared split
        // Exclusive atoms get flipped, shared atoms stay as pivots
        var phase2Candidates = [];
        var allRingAtomSets = [];
        for (i = 0; i < rings.length; i++) {
            var rset = {};
            for (k = 0; k < rings[i].length; k++) rset[rings[i][k]] = true;
            allRingAtomSets.push(rset);
        }
        for (i = 0; i < rings.length; i++) {
            var exclusive = [], shared = [];
            for (k = 0; k < rings[i].length; k++) {
                var aid = rings[i][k];
                var inOther = false;
                for (j = 0; j < rings.length; j++) {
                    if (j !== i && allRingAtomSets[j][aid]) { inOther = true; break; }
                }
                if (inOther) shared.push(aid);
                else exclusive.push(aid);
            }
            if (exclusive.length >= 1 && shared.length >= 2) {
                phase2Candidates.push({ exclusive: exclusive, shared: shared });
            }
        }

        /** Flip only exclusive atoms of a ring, keeping shared as pivots */
        function flipRingExclusive(cand) {
            if (cand.shared.length < 2) return;
            var p1 = mol.getAtom(cand.shared[0]), p2 = mol.getAtom(cand.shared[1]);
            var dx = p2.x - p1.x, dy = p2.y - p1.y;
            var len2 = dx * dx + dy * dy;
            if (len2 < 0.001) return;
            for (var ei = 0; ei < cand.exclusive.length; ei++) {
                var atom = mol.getAtom(cand.exclusive[ei]);
                var px = atom.x - p1.x, py = atom.y - p1.y;
                var dot = (px * dx + py * dy) / len2;
                atom.x = p1.x + 2 * dot * dx - px;
                atom.y = p1.y + 2 * dot * dy - py;
            }
        }

        // Two-phase simulated annealing
        var currentCrossings = countCrossings();
        if (currentCrossings === 0) { return; }

        var rng = 12345;
        function nextRandom() {
            rng = (rng * 1103515245 + 12345) & 0x7fffffff;
            return rng / 0x7fffffff;
        }

        var phase1Iters = Math.floor(maxIter * 0.4);
        var phase2Iters = maxIter - phase1Iters;

        // Phase 1: flip whole ring systems
        for (var iter = 0; iter < phase1Iters && currentCrossings > 0; iter++) {
            if (phase1Candidates.length === 0) break;
            var idx = Math.floor(nextRandom() * phase1Candidates.length) % phase1Candidates.length;
            flipSystem(phase1Candidates[idx]);
            var nc = countCrossings();
            if (nc < currentCrossings) { currentCrossings = nc; }
            else { unflipSystem(phase1Candidates[idx]); }
        }

        // Phase 2: flip individual rings (exclusive atoms only)
        for (var iter2 = 0; iter2 < phase2Iters && currentCrossings > 0; iter2++) {
            if (phase2Candidates.length === 0) break;
            var idx2 = Math.floor(nextRandom() * phase2Candidates.length) % phase2Candidates.length;
            var cand = phase2Candidates[idx2];
            flipRingExclusive(cand);
            var nc2 = countCrossings();
            if (nc2 < currentCrossings) { currentCrossings = nc2; }
            else { flipRingExclusive(cand); } // self-inverse
        }
    }

    // =====================================================================
    // Scaffold Templates — 10 canonical ring coordinate arrays
    // Each template is { name, coords: [[x,y], ...] } with BOND_LENGTH spacing
    // =====================================================================

    var BL = (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) ? Molecule.BOND_LENGTH : 30;

    var SCAFFOLD_TEMPLATES = [
        // 1. Benzene (6-ring, regular hexagon)
        { name: 'benzene', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            return c;
        })() },
        // 2. Naphthalene (10 atoms, 2 fused 6-rings)
        { name: 'naphthalene', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 4; i++) {
                var a = -Math.PI / 6 + i * Math.PI / 3;
                c.push([dx + r * Math.cos(a), r * Math.sin(a)]);
            }
            return c;
        })() },
        // 3. Indole (9 atoms, 5+6 fused)
        { name: 'indole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            var ax = c[0][0]; var ay = c[0][1];
            var bx = c[1][0]; var by = c[1][1];
            var mx = (ax + bx) / 2; var my = (ay + by) / 2;
            var nx = -(by - ay); var ny = bx - ax;
            var nl = Math.sqrt(nx * nx + ny * ny);
            nx /= nl; ny /= nl;
            var d5 = BL * 0.8507;
            c.push([mx + nx * d5, my + ny * d5]);
            var ang1 = Math.atan2(c[6][1] - ay, c[6][0] - ax);
            var ang2 = Math.atan2(c[6][1] - by, c[6][0] - bx);
            c.push([c[6][0] + BL * 0.5 * Math.cos(ang1 + 0.6), c[6][1] + BL * 0.5 * Math.sin(ang1 + 0.6)]);
            c.push([c[6][0] + BL * 0.5 * Math.cos(ang2 - 0.6), c[6][1] + BL * 0.5 * Math.sin(ang2 - 0.6)]);
            return c;
        })() },
        // 4. Purine (9 atoms, 5+6 fused nitrogen heterocycle)
        { name: 'purine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            var mx = (c[0][0] + c[1][0]) / 2;
            var my = (c[0][1] + c[1][1]) / 2;
            var nx = -(c[1][1] - c[0][1]);
            var ny = c[1][0] - c[0][0];
            var nl = Math.sqrt(nx * nx + ny * ny);
            c.push([mx + nx / nl * BL * 0.85, my + ny / nl * BL * 0.85]);
            c.push([(c[0][0] + c[6][0]) / 2 + nx / nl * BL * 0.3,
                    (c[0][1] + c[6][1]) / 2 + ny / nl * BL * 0.3]);
            c.push([(c[1][0] + c[6][0]) / 2 + nx / nl * BL * 0.3,
                    (c[1][1] + c[6][1]) / 2 + ny / nl * BL * 0.3]);
            return c;
        })() },
        // 5. Steroid (17 atoms, 4-ring ABCD skeleton)
        { name: 'steroid', coords: (function() {
            var c = []; var r = BL; var s3 = Math.sqrt(3);
            // Ring A (6)
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            // Ring B (6) fused to A at atoms 4,5
            var ox = BL * s3; var oy = 0;
            for (var i = 0; i < 4; i++) {
                var a = -Math.PI / 6 + i * Math.PI / 3;
                c.push([ox + r * Math.cos(a), oy + r * Math.sin(a)]);
            }
            // Ring C (6) fused to B
            var ox2 = 2 * BL * s3;
            for (var i = 0; i < 4; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([ox2 + r * Math.cos(a), r * Math.sin(a)]);
            }
            // Ring D (5) fused to C
            var ox3 = 3 * BL * s3 * 0.95;
            for (var i = 0; i < 3; i++) {
                var a = -Math.PI / 5 + i * 2 * Math.PI / 5;
                c.push([ox3 + BL * 0.85 * Math.cos(a), BL * 0.85 * Math.sin(a)]);
            }
            return c;
        })() },
        // 6. Morphinan (16 atoms, pentacyclic)
        { name: 'morphinan', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            var s3 = Math.sqrt(3);
            var ox = BL * s3;
            for (var i = 0; i < 4; i++) {
                var a = -Math.PI / 6 + i * Math.PI / 3;
                c.push([ox + r * Math.cos(a), r * Math.sin(a)]);
            }
            var ox2 = 2 * BL * s3;
            for (var i = 0; i < 4; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([ox2 + r * Math.cos(a), r * Math.sin(a)]);
            }
            c.push([BL * s3, -BL * 1.2]);
            c.push([BL * s3 * 0.5, -BL * 0.9]);
            return c;
        })() },
        // 7. Quinoline (9 atoms, 6+6 fused with N)
        { name: 'quinoline', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            var s3 = Math.sqrt(3);
            for (var i = 0; i < 4; i++) {
                var a = -Math.PI / 6 + i * Math.PI / 3;
                c.push([BL * s3 + r * Math.cos(a), r * Math.sin(a)]);
            }
            // Collapse shared atoms (naphthalene-like but 9 unique atoms)
            return c.slice(0, 9);
        })() },
        // 8. Biphenyl (12 atoms, 2 separate 6-rings)
        { name: 'biphenyl', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            var gap = BL * 2.5;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([gap + r * Math.cos(a), r * Math.sin(a)]);
            }
            return c;
        })() },
        // 9. Cyclohexane (6 atoms, chair-like flat)
        { name: 'cyclohexane', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            return c;
        })() },
        // 10. Piperidine (6 atoms, N-containing 6-ring)
        { name: 'piperidine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) {
                var a = Math.PI / 2 + i * Math.PI / 3;
                c.push([r * Math.cos(a), r * Math.sin(a)]);
            }
            return c;
        })() },
        // 11. Pyridine (6-ring with N)
        { name: 'pyridine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 12. Pyrimidine (6-ring with 2N)
        { name: 'pyrimidine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 13. Imidazole (5-ring with 2N)
        { name: 'imidazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 14. Pyrazole (5-ring with 2 adjacent N)
        { name: 'pyrazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 15. Thiophene (5-ring with S)
        { name: 'thiophene', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 16. Furan (5-ring with O)
        { name: 'furan', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 17. Pyrrolidine (5-ring saturated N)
        { name: 'pyrrolidine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 18. Piperazine (6-ring with 2N)
        { name: 'piperazine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 19. Morpholine (6-ring with N and O)
        { name: 'morpholine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 20. Tetrahydropyran (6-ring with O)
        { name: 'tetrahydropyran', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 21. Isoquinoline (6+6 fused, N in second ring)
        { name: 'isoquinoline', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            return c;
        })() },
        // 22. Quinazoline (6+6 fused, 2N in second ring)
        { name: 'quinazoline', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            return c;
        })() },
        // 23. Benzimidazole (6+5 fused)
        { name: 'benzimidazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            var cx = (c[0][0] + c[1][0]) / 2, cy = (c[0][1] + c[1][1]) / 2;
            var nx = -(c[1][1] - c[0][1]), ny = c[1][0] - c[0][0];
            var nl = Math.sqrt(nx*nx+ny*ny); nx /= nl; ny /= nl;
            var pr = BL * 0.85;
            c.push([cx + nx*pr*0.5, cy + ny*pr*0.5]);
            c.push([cx + nx*pr, cy + ny*pr]);
            c.push([cx + nx*pr*0.5, cy + ny*pr*1.5]);
            return c;
        })() },
        // 24. Benzofuran (6+5 fused with O)
        { name: 'benzofuran', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            var cx = (c[0][0] + c[1][0]) / 2, cy = (c[0][1] + c[1][1]) / 2;
            var nx = -(c[1][1] - c[0][1]), ny = c[1][0] - c[0][0];
            var nl = Math.sqrt(nx*nx+ny*ny); nx /= nl; ny /= nl;
            var pr = BL * 0.85;
            c.push([cx + nx*pr*0.5, cy + ny*pr*0.5]);
            c.push([cx + nx*pr, cy + ny*pr]);
            c.push([cx + nx*pr*0.5, cy + ny*pr*1.5]);
            return c;
        })() },
        // 25. Benzothiazole (6+5 fused with N,S)
        { name: 'benzothiazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            var cx = (c[0][0] + c[1][0]) / 2, cy = (c[0][1] + c[1][1]) / 2;
            var nx = -(c[1][1] - c[0][1]), ny = c[1][0] - c[0][0];
            var nl = Math.sqrt(nx*nx+ny*ny); nx /= nl; ny /= nl;
            var pr = BL * 0.85;
            c.push([cx + nx*pr*0.5, cy + ny*pr*0.5]);
            c.push([cx + nx*pr, cy + ny*pr]);
            c.push([cx + nx*pr*0.5, cy + ny*pr*1.5]);
            return c;
        })() },
        // 26. Phenanthrene (3 fused 6-rings, angular)
        { name: 'phenanthrene', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([2*dx + r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 27. Acridine (3 fused 6-rings with central N)
        { name: 'acridine', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([2*dx + r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 28. Triazole 1,2,3 (5-ring with 3N)
        { name: 'triazole123', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 29. Triazole 1,2,4 (5-ring with 3N)
        { name: 'triazole124', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 30. Tetrazole (5-ring with 4N)
        { name: 'tetrazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 31. Oxazole (5-ring with N,O)
        { name: 'oxazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 32. Isoxazole (5-ring with N,O adjacent)
        { name: 'isoxazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 33. Thiazole (5-ring with N,S)
        { name: 'thiazole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 34. Pyrrole (5-ring with NH)
        { name: 'pyrrole', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 35. Carbazole (indole + extra 6-ring, 3 fused)
        { name: 'carbazole', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 3; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([2*dx + r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() },
        // 36. Pyrrolopyrimidine (5+6 fused, kinase scaffold)
        { name: 'pyrrolopyrimidine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            var cx = (c[0][0] + c[4][0]) / 2, cy = (c[0][1] + c[4][1]) / 2;
            var nx = (c[4][1] - c[0][1]), ny = -(c[4][0] - c[0][0]);
            var nl = Math.sqrt(nx*nx+ny*ny); nx /= nl; ny /= nl;
            for (var i = 0; i < 4; i++) { c.push([cx + nx*BL*(0.5+i*0.4), cy + ny*BL*(0.5+i*0.4)]); }
            return c;
        })() },
        // 37. Benzodiazepine (6+7 fused, anxiolytics)
        { name: 'benzodiazepine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            var cx = (c[0][0] + c[1][0]) / 2, cy = (c[0][1] + c[1][1]) / 2;
            var nx = -(c[1][1] - c[0][1]), ny = c[1][0] - c[0][0];
            var nl = Math.sqrt(nx*nx+ny*ny); nx /= nl; ny /= nl;
            for (var i = 0; i < 5; i++) {
                var ang = -Math.PI/3 + i * Math.PI / 4;
                c.push([cx + nx*BL*1.2*Math.cos(ang), cy + ny*BL*1.2*Math.sin(ang)]);
            }
            return c;
        })() },
        // 38. Beta-lactam (4-ring, penicillin/cephalosporin core)
        { name: 'betalactam', coords: (function() {
            var c = []; var r = BL * 0.7;
            for (var i = 0; i < 4; i++) { c.push([r * Math.cos(Math.PI/4 + i*Math.PI/2), r * Math.sin(Math.PI/4 + i*Math.PI/2)]); }
            return c;
        })() },
        // 39. Cyclopropane (3-ring)
        { name: 'cyclopropane', coords: (function() {
            var c = []; var r = BL * 0.6;
            for (var i = 0; i < 3; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/3), r * Math.sin(Math.PI/2 + i*2*Math.PI/3)]); }
            return c;
        })() },
        // 40. Cyclopentane (5-ring saturated)
        { name: 'cyclopentane', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            return c;
        })() },
        // 41. Cyclobutane (4-ring saturated)
        { name: 'cyclobutane', coords: (function() {
            var c = []; var r = BL * 0.7;
            for (var i = 0; i < 4; i++) { c.push([r * Math.cos(Math.PI/4 + i*Math.PI/2), r * Math.sin(Math.PI/4 + i*Math.PI/2)]); }
            return c;
        })() },
        // 42. Cycloheptane (7-ring)
        { name: 'cycloheptane', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 7; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/7), r * Math.sin(Math.PI/2 + i*2*Math.PI/7)]); }
            return c;
        })() },
        // 43. Pteridine (purine variant, 2 fused 6-rings with 4N)
        { name: 'pteridine', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            return c;
        })() },
        // 44. Xanthine (purine + carbonyl, caffeine/theophylline core)
        { name: 'xanthine', coords: (function() {
            var c = []; var r = BL;
            for (var i = 0; i < 5; i++) { c.push([r * Math.cos(Math.PI/2 + i*2*Math.PI/5), r * Math.sin(Math.PI/2 + i*2*Math.PI/5)]); }
            var cx = (c[3][0] + c[4][0]) / 2, cy = (c[3][1] + c[4][1]) / 2;
            var nx = (c[4][1] - c[3][1]), ny = -(c[4][0] - c[3][0]);
            var nl = Math.sqrt(nx*nx+ny*ny); nx /= nl; ny /= nl;
            for (var i = 0; i < 4; i++) { var a = i * Math.PI/3; c.push([cx + nx*BL*Math.cos(a), cy + ny*BL*Math.sin(a)]); }
            return c;
        })() },
        // 45. Phenothiazine (3 fused 6-rings with N,S, antipsychotics)
        { name: 'phenothiazine', coords: (function() {
            var c = []; var r = BL; var dx = BL * Math.sqrt(3);
            for (var i = 0; i < 6; i++) { c.push([r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([dx + r * Math.cos(-Math.PI/6 + i*Math.PI/3), r * Math.sin(-Math.PI/6 + i*Math.PI/3)]); }
            for (var i = 0; i < 4; i++) { c.push([2*dx + r * Math.cos(Math.PI/2 + i*Math.PI/3), r * Math.sin(Math.PI/2 + i*Math.PI/3)]); }
            return c;
        })() }
    ];

    // =====================================================================
    // BFS shortest-path distance matrix (for SMACOF target distances)
    // =====================================================================

    function bfsDistMatrix(mol, atomIds) {
        var n = atomIds.length;
        var idToIdx = {};
        var i, j;
        for (i = 0; i < n; i++) { idToIdx[atomIds[i]] = i; }
        var dist = new Array(n);
        for (i = 0; i < n; i++) {
            dist[i] = new Array(n);
            for (j = 0; j < n; j++) { dist[i][j] = (i === j) ? 0 : -1; }
        }
        for (var src = 0; src < n; src++) {
            var queue = [src];
            var head = 0;
            while (head < queue.length) {
                var u = queue[head++];
                var uid = atomIds[u];
                var nbs = mol.getNeighbors(uid);
                for (j = 0; j < nbs.length; j++) {
                    var vi = idToIdx[nbs[j]];
                    if (vi !== undefined && dist[src][vi] < 0) {
                        dist[src][vi] = dist[src][u] + 1;
                        queue.push(vi);
                    }
                }
            }
            // Fill unreachable with large value
            for (j = 0; j < n; j++) {
                if (dist[src][j] < 0) { dist[src][j] = n; }
            }
        }
        return dist;
    }

    // =====================================================================
    // smacofLayout — SMACOF stress majorisation (de Leeuw)
    // =====================================================================

    /**
     * SMACOF stress majorisation with multi-init random starts.
     * Computes target distance matrix from graph shortest paths * BOND_LENGTH,
     * then iteratively minimises stress. Keeps lowest-stress result.
     *
     * @param {Molecule} mol     — molecule with atom coordinates
     * @param {Array<number>} atomIds — atom indices to layout
     * @param {number} [maxIter] — iterations per init (default 100)
     * @param {number} [nInit]   — number of random starts (default 3)
     * @returns {number} — final stress value
     */
    function smacofLayout(mol, atomIds, maxIter, nInit) {
        if (!atomIds || atomIds.length < 3) { return 0; }
        maxIter = maxIter || 100;
        nInit = nInit || 3;

        var n = atomIds.length;
        var bondLen = (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) ? Molecule.BOND_LENGTH : 30;

        // Build target distance matrix from graph shortest paths
        var graphDist = bfsDistMatrix(mol, atomIds);
        var target = new Array(n);
        var weights = new Array(n);
        var i, j;
        for (i = 0; i < n; i++) {
            target[i] = new Array(n);
            weights[i] = new Array(n);
            for (j = 0; j < n; j++) {
                target[i][j] = graphDist[i][j] * bondLen;
                // Inverse-square distance weights (de Leeuw standard)
                weights[i][j] = (i !== j && target[i][j] > 0)
                    ? 1.0 / (target[i][j] * target[i][j]) : 0;
            }
        }

        // Simple LCG for reproducible random inits
        var rngState = 42;
        function nextRng() {
            rngState = (rngState * 1103515245 + 12345) & 0x7fffffff;
            return rngState / 0x7fffffff;
        }

        var bestStress = Infinity;
        var bestX = new Array(n);
        var bestY = new Array(n);

        for (var init = 0; init < nInit; init++) {
            // Initialise: first init uses current coords, rest use random
            var x = new Array(n);
            var y = new Array(n);
            if (init === 0) {
                for (i = 0; i < n; i++) {
                    var atom = mol.getAtom(atomIds[i]);
                    x[i] = atom.x;
                    y[i] = atom.y;
                }
            } else {
                for (i = 0; i < n; i++) {
                    x[i] = (nextRng() - 0.5) * bondLen * n * 0.5;
                    y[i] = (nextRng() - 0.5) * bondLen * n * 0.5;
                }
            }

            // SMACOF iterations
            for (var iter = 0; iter < maxIter; iter++) {
                // Compute current distances
                var curDist = new Array(n);
                for (i = 0; i < n; i++) {
                    curDist[i] = new Array(n);
                    for (j = 0; j < n; j++) {
                        if (i === j) { curDist[i][j] = 0; continue; }
                        var dx = x[i] - x[j], dy = y[i] - y[j];
                        curDist[i][j] = Math.sqrt(dx * dx + dy * dy);
                    }
                }

                // Guttman transform: X_new = V^+ B(X) X
                var newX = new Array(n);
                var newY = new Array(n);
                for (i = 0; i < n; i++) {
                    var sumWi = 0;
                    var bx = 0, by = 0;
                    for (j = 0; j < n; j++) {
                        if (i === j) { continue; }
                        sumWi += weights[i][j];
                        var dij = curDist[i][j];
                        if (dij > 1e-12) {
                            var ratio = weights[i][j] * target[i][j] / dij;
                            bx += ratio * (x[i] - x[j]);
                            by += ratio * (y[i] - y[j]);
                        }
                    }
                    if (sumWi > 1e-12) {
                        newX[i] = (1.0 / sumWi) * bx;
                        newY[i] = (1.0 / sumWi) * by;
                    } else {
                        newX[i] = x[i];
                        newY[i] = y[i];
                    }
                }

                // Centre to origin
                var cx = 0, cy = 0;
                for (i = 0; i < n; i++) { cx += newX[i]; cy += newY[i]; }
                cx /= n; cy /= n;
                for (i = 0; i < n; i++) { newX[i] -= cx; newY[i] -= cy; }

                x = newX;
                y = newY;
            }

            // Compute stress
            var stress = 0;
            for (i = 0; i < n; i++) {
                for (j = i + 1; j < n; j++) {
                    var dx = x[i] - x[j], dy = y[i] - y[j];
                    var d = Math.sqrt(dx * dx + dy * dy);
                    var diff = d - target[i][j];
                    stress += weights[i][j] * diff * diff;
                }
            }

            if (stress < bestStress) {
                bestStress = stress;
                for (i = 0; i < n; i++) { bestX[i] = x[i]; bestY[i] = y[i]; }
            }
        }

        // Apply best coordinates
        for (i = 0; i < n; i++) {
            var atom = mol.getAtom(atomIds[i]);
            atom.x = bestX[i];
            atom.y = bestY[i];
        }
        return bestStress;
    }

    // =====================================================================
    // forceDirectedLayout — PrEd-inspired (Bertault 2000) crossing-aware
    // =====================================================================

    /**
     * Force-directed refinement that preserves edge crossing zones.
     * Attractive springs between bonded atoms, repulsive forces between
     * non-bonded, and crossing penalty forces perpendicular to crossing
     * segments (Bertault 2000 PrEd approach).
     *
     * @param {Molecule} mol          — molecule with 2D coordinates
     * @param {Array<number>} atomIds — atoms in this component
     * @param {Object} ringAtomSet    — { atomId: true } for ring atoms
     * @param {Array<Array<number>>} rings — SSSR rings
     * @param {number} [maxIter]      — iterations (default 80)
     */
    function forceDirectedLayout(mol, atomIds, ringAtomSet, rings, maxIter) {
        if (!atomIds || atomIds.length < 3) { return; }
        maxIter = maxIter || 80;

        var n = atomIds.length;
        var bondLen = (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) ? Molecule.BOND_LENGTH : 30;
        var i, j, k;

        // Build lookup
        var idToIdx = {};
        for (i = 0; i < n; i++) { idToIdx[atomIds[i]] = i; }

        // Collect bonds within this component
        var bonds = [];
        var bondAdj = new Array(n);
        for (i = 0; i < n; i++) { bondAdj[i] = []; }

        for (i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            var ui = idToIdx[b.atom1];
            var vi = idToIdx[b.atom2];
            if (ui !== undefined && vi !== undefined) {
                bonds.push([ui, vi]);
                bondAdj[ui].push(vi);
                bondAdj[vi].push(ui);
            }
        }

        // Ring atom flags (local indices)
        var isRing = new Array(n);
        for (i = 0; i < n; i++) {
            isRing[i] = ringAtomSet && ringAtomSet[atomIds[i]] ? true : false;
        }

        // Force parameters
        var kSpring = 0.3;          // attractive spring constant
        var kRepulse = bondLen * bondLen * 2.0; // repulsive constant
        var kCrossing = 0.5;        // crossing penalty strength
        var damping = 0.85;         // velocity damping per iteration
        var maxDisp = bondLen * 0.4; // max displacement per step
        var ringDamping = 0.3;      // ring atoms move less

        // Position and velocity arrays
        var px = new Array(n), py = new Array(n);
        var vx = new Array(n), vy = new Array(n);
        for (i = 0; i < n; i++) {
            var atom = mol.getAtom(atomIds[i]);
            px[i] = atom.x;
            py[i] = atom.y;
            vx[i] = 0;
            vy[i] = 0;
        }

        // Pre-compute non-bonded pairs mask
        var isBonded = new Array(n);
        for (i = 0; i < n; i++) {
            isBonded[i] = {};
            for (j = 0; j < bondAdj[i].length; j++) {
                isBonded[i][bondAdj[i][j]] = true;
            }
        }

        for (var iter = 0; iter < maxIter; iter++) {
            var fx = new Array(n), fy = new Array(n);
            for (i = 0; i < n; i++) { fx[i] = 0; fy[i] = 0; }

            // 1. Attractive forces (springs along bonds)
            for (k = 0; k < bonds.length; k++) {
                var u = bonds[k][0], v = bonds[k][1];
                var dx = px[v] - px[u], dy = py[v] - py[u];
                var d = Math.sqrt(dx * dx + dy * dy);
                if (d < 1e-6) { continue; }
                var stretch = (d - bondLen) / d;
                var fAtt = kSpring * stretch;
                fx[u] += fAtt * dx;
                fy[u] += fAtt * dy;
                fx[v] -= fAtt * dx;
                fy[v] -= fAtt * dy;
            }

            // 2. Repulsive forces (non-bonded pairs, Coulomb-like)
            for (i = 0; i < n; i++) {
                for (j = i + 1; j < n; j++) {
                    if (isBonded[i][j]) { continue; }
                    var dx = px[j] - px[i], dy = py[j] - py[i];
                    var d2 = dx * dx + dy * dy;
                    if (d2 < 1e-6) { d2 = 1e-6; }
                    var fRep = -kRepulse / d2;
                    var d = Math.sqrt(d2);
                    fx[i] += fRep * dx / d;
                    fy[i] += fRep * dy / d;
                    fx[j] -= fRep * dx / d;
                    fy[j] -= fRep * dy / d;
                }
            }

            // 3. Crossing penalty forces (Bertault PrEd)
            // For each pair of crossing edges, push atoms perpendicular
            // to the crossing segments to reduce overlap
            for (var bi = 0; bi < bonds.length; bi++) {
                var a1 = bonds[bi][0], a2 = bonds[bi][1];
                for (var bj = bi + 1; bj < bonds.length; bj++) {
                    var a3 = bonds[bj][0], a4 = bonds[bj][1];
                    // Skip shared-endpoint bonds
                    if (a1 === a3 || a1 === a4 || a2 === a3 || a2 === a4) { continue; }
                    if (!segmentsIntersect(px[a1], py[a1], px[a2], py[a2],
                                           px[a3], py[a3], px[a4], py[a4])) {
                        continue;
                    }
                    // Compute perpendicular push for crossing resolution
                    // Push endpoints of each edge away from the other edge
                    var mx1 = (px[a1] + px[a2]) / 2, my1 = (py[a1] + py[a2]) / 2;
                    var mx2 = (px[a3] + px[a4]) / 2, my2 = (py[a3] + py[a4]) / 2;

                    // Direction perpendicular to first edge
                    var e1x = px[a2] - px[a1], e1y = py[a2] - py[a1];
                    var e1len = Math.sqrt(e1x * e1x + e1y * e1y);
                    if (e1len < 1e-6) { continue; }
                    var perp1x = -e1y / e1len, perp1y = e1x / e1len;

                    // Direction perpendicular to second edge
                    var e2x = px[a4] - px[a3], e2y = py[a4] - py[a3];
                    var e2len = Math.sqrt(e2x * e2x + e2y * e2y);
                    if (e2len < 1e-6) { continue; }
                    var perp2x = -e2y / e2len, perp2y = e2x / e2len;

                    // Determine which side the midpoint of the other edge is on
                    var dot1 = (mx2 - mx1) * perp1x + (my2 - my1) * perp1y;
                    var sign1 = dot1 >= 0 ? 1 : -1;
                    var dot2 = (mx1 - mx2) * perp2x + (my1 - my2) * perp2y;
                    var sign2 = dot2 >= 0 ? 1 : -1;

                    // Push first edge endpoints away from crossing
                    var push = kCrossing * bondLen * 0.1;
                    fx[a1] -= sign1 * push * perp1x;
                    fy[a1] -= sign1 * push * perp1y;
                    fx[a2] -= sign1 * push * perp1x;
                    fy[a2] -= sign1 * push * perp1y;

                    // Push second edge endpoints away
                    fx[a3] -= sign2 * push * perp2x;
                    fy[a3] -= sign2 * push * perp2y;
                    fx[a4] -= sign2 * push * perp2x;
                    fy[a4] -= sign2 * push * perp2y;
                }
            }

            // 4. Update velocities and positions with damping
            var temperature = 1.0 - iter / maxIter; // cooling
            var maxF = maxDisp * temperature;
            for (i = 0; i < n; i++) {
                var scale = isRing[i] ? ringDamping : 1.0;
                vx[i] = (vx[i] + fx[i] * scale) * damping;
                vy[i] = (vy[i] + fy[i] * scale) * damping;

                // Clamp displacement
                var vLen = Math.sqrt(vx[i] * vx[i] + vy[i] * vy[i]);
                if (vLen > maxF) {
                    vx[i] = vx[i] / vLen * maxF;
                    vy[i] = vy[i] / vLen * maxF;
                }
                px[i] += vx[i];
                py[i] += vy[i];
            }
        }

        // Write back positions
        for (i = 0; i < n; i++) {
            var atom = mol.getAtom(atomIds[i]);
            atom.x = px[i];
            atom.y = py[i];
        }
    }

    // =====================================================================
    // Public API
    // =====================================================================

    var SMSDLayout = {
        computeSSSR: computeSSSR,
        layoutSSSR: layoutSSSR,
        reduceCrossings: reduceCrossings,
        forceDirectedLayout: forceDirectedLayout,
        smacofLayout: smacofLayout,
        SCAFFOLD_TEMPLATES: SCAFFOLD_TEMPLATES,
        Point2D: Point2D,
        _segmentsIntersect: segmentsIntersect,
        _bfsDistMatrix: bfsDistMatrix
    };

    global.SMSDLayout = SMSDLayout;

})(typeof window !== 'undefined' ? window : this);
