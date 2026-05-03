/**
 * SMSDRings.js — SSSR, Relevant Cycles and Unique Ring Families
 * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * Licensed under the Apache License, Version 2.0
 *
 * Ported from SMSD C++ ring_finder.hpp (header-only).
 *
 * Implements:
 *   - Horton candidate generation
 *   - 2-phase GF(2) elimination (Vismara's rule)
 *   - SSSR (Minimum Cycle Basis)  — computeRings()
 *   - Relevant Cycle Basis        — computeRelevantCycles()
 *   - Unique Ring Families         — computeURFs()
 *
 * Depends on: SMSDGraph (window.SMSDGraph)
 */
(function() {
    'use strict';

    // ========================================================================
    // Portable bit helpers (operating on arrays of 32-bit words for GF(2))
    // ========================================================================

    /**
     * Count trailing zeros in a 32-bit integer.
     * @param {number} x
     * @return {number}
     */
    function ctz32(x) {
        if (x === 0) return 32;
        var n = 0;
        if ((x & 0x0000FFFF) === 0) { n += 16; x >>>= 16; }
        if ((x & 0x000000FF) === 0) { n += 8;  x >>>= 8;  }
        if ((x & 0x0000000F) === 0) { n += 4;  x >>>= 4;  }
        if ((x & 0x00000003) === 0) { n += 2;  x >>>= 2;  }
        if ((x & 0x00000001) === 0) { n += 1; }
        return n;
    }

    // ========================================================================
    // Edge key helpers
    // ========================================================================

    /**
     * Canonical edge key string "lo:hi".
     * @param {number} i
     * @param {number} j
     * @return {string}
     */
    function edgeKey(i, j) {
        var lo = i < j ? i : j;
        var hi = i < j ? j : i;
        return lo + ':' + hi;
    }

    // ========================================================================
    // BFS shortest path avoiding a specific edge
    // ========================================================================

    /**
     * BFS shortest path from start to target, excluding the undirected
     * edge (avoidU, avoidV).
     *
     * @param {Object} g          Graph with .n and .neighbors[]
     * @param {number} start
     * @param {number} target
     * @param {number} avoidU
     * @param {number} avoidV
     * @return {Array<number>}     Atom path [start, ..., target] or []
     */
    function bfsPath(g, start, target, avoidU, avoidV) {
        if (start === target) return [start];
        var n = g.n;
        var parent = new Int32Array(n);
        var vis = new Uint8Array(n);
        var i;
        for (i = 0; i < n; i++) parent[i] = -1;
        vis[start] = 1;
        var queue = [start];
        var head = 0;
        while (head < queue.length) {
            var u = queue[head++];
            var nbs = g.neighbors[u];
            for (var k = 0; k < nbs.length; k++) {
                var v = nbs[k];
                if ((u === avoidU && v === avoidV) || (u === avoidV && v === avoidU))
                    continue;
                if (vis[v]) continue;
                vis[v] = 1;
                parent[v] = u;
                if (v === target) {
                    var path = [];
                    for (var c = target; c !== -1; c = parent[c]) path.push(c);
                    if (path[path.length - 1] !== start) path.push(start);
                    path.reverse();
                    return path;
                }
                queue.push(v);
            }
        }
        return [];
    }

    // ========================================================================
    // Graph topology helpers
    // ========================================================================

    /**
     * Build edge index map: edgeKey -> sequential integer index.
     * @param {Object} g
     * @return {Object}  { map: {string: number}, list: Array<[number,number]> }
     */
    function buildEdgeIndex(g) {
        var map = {};
        var list = [];
        var idx = 0;
        for (var i = 0; i < g.n; i++) {
            var nbs = g.neighbors[i];
            for (var k = 0; k < nbs.length; k++) {
                var j = nbs[k];
                if (j > i) {
                    var key = edgeKey(i, j);
                    if (!(key in map)) {
                        map[key] = idx;
                        list.push([i, j]);
                        idx++;
                    }
                }
            }
        }
        return { map: map, list: list, count: idx };
    }

    /**
     * Count connected components via BFS.
     * @param {Object} g
     * @return {number}
     */
    function countComponents(g) {
        var n = g.n;
        var vis = new Uint8Array(n);
        var comp = 0;
        for (var i = 0; i < n; i++) {
            if (vis[i]) continue;
            comp++;
            var queue = [i];
            vis[i] = 1;
            var head = 0;
            while (head < queue.length) {
                var u = queue[head++];
                var nbs = g.neighbors[u];
                for (var k = 0; k < nbs.length; k++) {
                    var v = nbs[k];
                    if (!vis[v]) { vis[v] = 1; queue.push(v); }
                }
            }
        }
        return comp;
    }

    /**
     * Count edges in the graph.
     * @param {Object} g
     * @return {number}
     */
    function countEdges(g) {
        var total = 0;
        for (var i = 0; i < g.n; i++) total += g.neighbors[i].length;
        return total / 2;
    }

    /**
     * Cycle rank = |E| - |V| + components.
     * @param {Object} g
     * @return {number}
     */
    function cycleRank(g) {
        return countEdges(g) - g.n + countComponents(g);
    }

    // ========================================================================
    // GF(2) vector helpers (using Int32Array words, 32 bits per word)
    // ========================================================================

    /**
     * Create a zero GF(2) vector with the given number of 32-bit words.
     * @param {number} numWords
     * @return {Int32Array}
     */
    function gf2zero(numWords) {
        return new Int32Array(numWords);
    }

    /**
     * XOR vector b into vector a (in place).
     * @param {Int32Array} a
     * @param {Int32Array} b
     */
    function gf2xor(a, b) {
        for (var i = 0; i < a.length; i++) a[i] ^= b[i];
    }

    /**
     * Clone a GF(2) vector.
     * @param {Int32Array} v
     * @return {Int32Array}
     */
    function gf2clone(v) {
        var c = new Int32Array(v.length);
        for (var i = 0; i < v.length; i++) c[i] = v[i];
        return c;
    }

    /**
     * Check if a GF(2) vector is all zeros.
     * @param {Int32Array} v
     * @return {boolean}
     */
    function gf2isZero(v) {
        for (var i = 0; i < v.length; i++) if (v[i] !== 0) return false;
        return true;
    }

    /**
     * Find the lowest set bit position in a GF(2) vector, or -1 if zero.
     * @param {Int32Array} v
     * @return {number}
     */
    function gf2findPivot(v) {
        for (var w = 0; w < v.length; w++) {
            if (v[w] !== 0) return (w << 5) + ctz32(v[w] >>> 0);
        }
        return -1;
    }

    /**
     * Build edge-incidence GF(2) vector for a cycle.
     * @param {Array<number>} cycle    Atom indices forming the cycle
     * @param {Object}        edgeMap  edgeKey -> index
     * @param {number}        numWords
     * @return {Int32Array}
     */
    function cycleToEdgeVec(cycle, edgeMap, numWords) {
        var vec = gf2zero(numWords);
        for (var i = 0; i < cycle.length; i++) {
            var a = cycle[i];
            var b = cycle[(i + 1) % cycle.length];
            var key = edgeKey(a, b);
            if (key in edgeMap) {
                var idx = edgeMap[key];
                vec[idx >> 5] ^= (1 << (idx & 31));
            }
        }
        return vec;
    }

    /**
     * Canonical edge-set key for cycle deduplication.
     * @param {Array<number>} cycle
     * @return {string}
     */
    function canonicalCycleKey(cycle) {
        var keys = [];
        for (var i = 0; i < cycle.length; i++) {
            var a = cycle[i];
            var b = cycle[(i + 1) % cycle.length];
            var lo = a < b ? a : b;
            var hi = a < b ? b : a;
            keys.push(lo * 100000 + hi);
        }
        keys.sort(function(x, y) { return x - y; });
        return keys.join(',');
    }

    // ========================================================================
    // Union-Find
    // ========================================================================

    function ufFind(parent, i) {
        while (parent[i] !== i) {
            parent[i] = parent[parent[i]]; // path halving
            i = parent[i];
        }
        return i;
    }

    function ufUnion(parent, a, b) {
        a = ufFind(parent, a);
        b = ufFind(parent, b);
        if (a !== b) parent[a] = b;
    }

    // ========================================================================
    // Horton + GF(2) Pipeline (shared between SSSR and RCB)
    // ========================================================================

    /**
     * Run the full Horton candidate generation and 2-phase GF(2) pipeline.
     *
     * @param {Object} g  Graph with .n and .neighbors[]
     * @return {{ sssr: Array<Array<number>>, rcb: Array<Array<number>> }}
     */
    function runHortonPipeline(g) {
        var result = { sssr: [], rcb: [] };
        if (!g || g.n === 0) return result;

        var cRank = cycleRank(g);
        if (cRank <= 0) return result;

        var edgeInfo = buildEdgeIndex(g);
        var edgeMap = edgeInfo.map;
        var edgeList = edgeInfo.list;
        var totalEdges = edgeInfo.count;
        var numWords = (totalEdges + 31) >> 5;

        // --- Horton candidate generation ---
        var seenMasks = {};
        var candidateVecs = [];
        var candidateAtoms = [];
        var candidateSizes = [];

        var v, ei, eU, eV, path1, path2;
        var p1Interior, disjoint, ringSize, atoms, vec;
        var key, maskKey;

        for (v = 0; v < g.n; v++) {
            for (ei = 0; ei < edgeList.length; ei++) {
                eU = edgeList[ei][0];
                eV = edgeList[ei][1];

                path1 = bfsPath(g, v, eU, eU, eV);
                if (path1.length === 0) continue;
                path2 = bfsPath(g, v, eV, eU, eV);
                if (path2.length === 0) continue;

                // Check that interior vertices of path1 and path2 are disjoint
                p1Interior = {};
                for (var pi = 1; pi < path1.length; pi++) p1Interior[path1[pi]] = 1;
                disjoint = true;
                for (var pj = 1; pj < path2.length - 1; pj++) {
                    if (path2[pj] in p1Interior) { disjoint = false; break; }
                }
                if (!disjoint) continue;

                ringSize = path1.length + path2.length - 1;
                if (ringSize < 3 || ringSize > 20) continue;

                // Build atom list for this candidate ring
                atoms = [];
                var ai;
                for (ai = 0; ai < path1.length; ai++) atoms.push(path1[ai]);
                for (ai = path2.length - 1; ai >= 1; ai--) atoms.push(path2[ai]);

                // Build GF(2) edge-incidence vector
                vec = gf2zero(numWords);
                // Add the avoided edge
                key = edgeKey(eU, eV);
                if (key in edgeMap) {
                    var eidx = edgeMap[key];
                    vec[eidx >> 5] ^= (1 << (eidx & 31));
                }
                var valid = true;
                for (ai = 0; ai + 1 < path1.length && valid; ai++) {
                    key = edgeKey(path1[ai], path1[ai + 1]);
                    if (!(key in edgeMap)) { valid = false; }
                    else {
                        eidx = edgeMap[key];
                        vec[eidx >> 5] ^= (1 << (eidx & 31));
                    }
                }
                for (ai = 0; ai + 1 < path2.length && valid; ai++) {
                    key = edgeKey(path2[ai], path2[ai + 1]);
                    if (!(key in edgeMap)) { valid = false; }
                    else {
                        eidx = edgeMap[key];
                        vec[eidx >> 5] ^= (1 << (eidx & 31));
                    }
                }
                if (!valid) continue;

                // Deduplicate by canonical edge set
                maskKey = canonicalCycleKey(atoms);
                if (maskKey in seenMasks) continue;
                seenMasks[maskKey] = 1;

                candidateVecs.push(vec);
                candidateAtoms.push(atoms);
                candidateSizes.push(ringSize);
            }
        }

        // Sort candidates by size (smallest first)
        var order = [];
        for (var oi = 0; oi < candidateVecs.length; oi++) order.push(oi);
        order.sort(function(a, b) { return candidateSizes[a] - candidateSizes[b]; });

        // --- 2-phase GF(2) interleaved pipeline ---
        var rcbSeen = {};
        var basisMatrix = [];
        var basisSizes = [];
        var basisByLeadingBit = new Int32Array(totalEdges);
        for (var bi = 0; bi < totalEdges; bi++) basisByLeadingBit[bi] = -1;

        for (var oo = 0; oo < order.length; oo++) {
            var idx = order[oo];
            var mask = candidateVecs[idx];
            var cSize = candidateSizes[idx];

            // Phase 1: reduce against STRICTLY SHORTER basis vectors only
            var phase1Mask = gf2clone(mask);
            for (var b = 0; b < totalEdges; b++) {
                if ((phase1Mask[b >> 5] >>> (b & 31)) & 1) {
                    var basisIdx = basisByLeadingBit[b];
                    if (basisIdx !== -1 && basisSizes[basisIdx] < cSize) {
                        gf2xor(phase1Mask, basisMatrix[basisIdx]);
                    }
                }
            }

            var relevant = !gf2isZero(phase1Mask);
            if (!relevant) continue;

            // Add to RCB (deduplicated)
            var rcbKey = canonicalCycleKey(candidateAtoms[idx]);
            if (!(rcbKey in rcbSeen)) {
                rcbSeen[rcbKey] = 1;
                result.rcb.push(candidateAtoms[idx]);
            }

            // Phase 2: reduce against ALL basis vectors for SSSR
            var phase2Mask = gf2clone(phase1Mask);
            for (b = 0; b < totalEdges; b++) {
                if ((phase2Mask[b >> 5] >>> (b & 31)) & 1) {
                    basisIdx = basisByLeadingBit[b];
                    if (basisIdx !== -1) {
                        gf2xor(phase2Mask, basisMatrix[basisIdx]);
                    }
                }
            }

            var independent = false;
            var newLeadingBit = -1;
            for (b = 0; b < totalEdges; b++) {
                if ((phase2Mask[b >> 5] >>> (b & 31)) & 1) {
                    independent = true;
                    newLeadingBit = b;
                    break;
                }
            }

            if (independent && result.sssr.length < cRank) {
                result.sssr.push(candidateAtoms[idx]);
                basisMatrix.push(phase2Mask);
                basisSizes.push(cSize);
                basisByLeadingBit[newLeadingBit] = basisMatrix.length - 1;
            }
        }

        return result;
    }

    // ========================================================================
    // Public API: computeRings (SSSR)
    // ========================================================================

    /**
     * Compute the Smallest Set of Smallest Rings (SSSR / Minimum Cycle Basis).
     *
     * Algorithm:
     *   1. Horton candidate generation: for every vertex v and every edge (eU,eV),
     *      BFS from v to eU and v to eV avoiding edge (eU,eV), form cycle from
     *      the two vertex-disjoint paths + the edge.
     *   2. Sort candidates by size.
     *   3. GF(2) Gaussian elimination to select a linearly independent set.
     *
     * @param {Object} g  Graph with .n and .neighbors[]
     * @return {Array<Array<number>>}  SSSR rings (each = ordered atom indices)
     */
    function computeRings(g) {
        return runHortonPipeline(g).sssr;
    }

    // ========================================================================
    // Public API: computeRelevantCycles
    // ========================================================================

    /**
     * Compute the Relevant Cycle Basis (union of all minimum cycle bases).
     * A cycle is relevant if it cannot be expressed as the XOR of strictly
     * shorter cycles.
     *
     * @param {Object} g  Graph with .n and .neighbors[]
     * @return {Array<Array<number>>}  Relevant cycles
     */
    function computeRelevantCycles(g) {
        return runHortonPipeline(g).rcb;
    }

    // ========================================================================
    // Public API: computeURFs (Unique Ring Families)
    // ========================================================================

    /**
     * Compute Unique Ring Families (Kolodzik et al. 2012).
     *
     * Two relevant cycles belong to the same URF if they have the same
     * orbit signature (lexmin rotation of orbit labels in both directions)
     * and all atoms are vertex-transitive (same orbit).
     *
     * @param {Object} g  Graph with .n, .neighbors[], .orbit[]
     * @return {Array<Array<Array<number>>>}  URFs, each = array of ring arrays
     */
    function computeURFs(g) {
        // Ensure canonical labelling / orbit info is available
        if (g.ensureCanonical) g.ensureCanonical();

        var rc = computeRelevantCycles(g);
        if (rc.length === 0) return [];

        var nrc = rc.length;
        var orbit = g.orbit || [];

        // Generate orbit-based signature for a cycle:
        // lexmin sequence of orbit[atom] values, considering both directions
        function getCycleOrbitSig(c) {
            var len = c.length;
            var best = null;
            for (var dir = 0; dir <= 1; dir++) {
                for (var start = 0; start < len; start++) {
                    var seq = [];
                    for (var i = 0; i < len; i++) {
                        var ci = dir === 0
                            ? (start + i) % len
                            : (start - i + len) % len;
                        seq.push(orbit[c[ci]] || 0);
                    }
                    if (best === null) {
                        best = seq;
                    } else {
                        // Lexicographic compare
                        for (var x = 0; x < len; x++) {
                            if (seq[x] < best[x]) { best = seq; break; }
                            if (seq[x] > best[x]) break;
                        }
                    }
                }
            }
            return best;
        }

        // Compute signatures
        var sigs = [];
        for (var si = 0; si < nrc; si++) sigs.push(getCycleOrbitSig(rc[si]));

        // Union-Find family array
        var family = [];
        for (var fi = 0; fi < nrc; fi++) family.push(fi);

        // Compare all pairs
        for (var i = 0; i < nrc; i++) {
            for (var j = i + 1; j < nrc; j++) {
                if (rc[i].length !== rc[j].length) continue;

                // Orbit-signature match, guarded: apply only when all atoms
                // in both cycles belong to the same orbit (vertex-transitive)
                var sigMatch = true;
                var len = sigs[i].length;
                for (var x = 0; x < len; x++) {
                    if (sigs[i][x] !== sigs[j][x]) { sigMatch = false; break; }
                }

                if (sigMatch) {
                    var orb0 = orbit[rc[i][0]] || 0;
                    var allSame = true;
                    var a;
                    for (a = 0; a < rc[i].length; a++) {
                        if ((orbit[rc[i][a]] || 0) !== orb0) { allSame = false; break; }
                    }
                    if (allSame) {
                        for (a = 0; a < rc[j].length; a++) {
                            if ((orbit[rc[j][a]] || 0) !== orb0) { allSame = false; break; }
                        }
                    }
                    if (allSame) {
                        ufUnion(family, i, j);
                    }
                }
            }
        }

        // Group by family
        var groups = {};
        for (var gi = 0; gi < nrc; gi++) {
            var rep = ufFind(family, gi);
            if (!(rep in groups)) groups[rep] = [];
            groups[rep].push(gi);
        }

        var result = [];
        for (var key in groups) {
            if (!groups.hasOwnProperty(key)) continue;
            var members = groups[key];
            var fam = [];
            for (var m = 0; m < members.length; m++) fam.push(rc[members[m]]);
            result.push(fam);
        }

        return result;
    }

    // ========================================================================
    // Export
    // ========================================================================

    window.SMSDRings = {
        computeRings: computeRings,
        computeRelevantCycles: computeRelevantCycles,
        computeURFs: computeURFs,
        // Expose internals for testing
        _bfsPath: bfsPath,
        _cycleRank: cycleRank
    };

})();
