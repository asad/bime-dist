/**
 * SMSDBatch.js — Batch substructure/MCS search, fingerprints and screening
 * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * Licensed under the Apache License, Version 2.0
 *
 * Ported from SMSD C++ batch.hpp (header-only, OpenMP sections skipped).
 * Single-threaded ES5 JavaScript port.
 *
 * Implements:
 *   - batchSubstructure:  1 query vs N targets (VF2PP)
 *   - batchMCS:           1 query vs N targets (MCS funnel)
 *   - screenAndMatch:     RASCAL pre-screen + exact MCS on hits
 *   - generateFingerprint: path-based molecular fingerprint
 *   - generateMcsFingerprint: MCS-aware path fingerprint
 *   - ecfpCounts:         count-based ECFP circular fingerprint
 *   - fcfpCounts:         count-based FCFP pharmacophore fingerprint
 *   - topologicalTorsion: 4-atom path torsion fingerprint (count-based)
 *   - tanimotoFingerprint: Tanimoto coefficient between fingerprints
 *   - diceSimilarity:     Dice similarity for binary fingerprints
 *   - cosineSimilarity:   Cosine similarity for binary fingerprints
 *   - soergelDistance:     Soergel distance (1 - Tanimoto)
 *   - countTanimoto/countDice/countCosine: count-vector similarity
 *   - fingerprintSubset:  subset check for substructure screening
 *   - fingerprintFromSmiles: one-call SMILES-to-fingerprint pipeline
 *   - countsToArray:      sparse count map to dense Float64Array
 *   - batchSimilarity:    query vs targets with sorted results
 *
 * Depends on: SMSDVF2 (window.SMSDVF2), SMSDMCS (window.SMSDMCS)
 */
(function() {
    'use strict';

    // ========================================================================
    // FNV-1a constants (32-bit variant for ES5 safe integer arithmetic)
    // ========================================================================

    var FNV1A_SEED  = 0x811C9DC5;
    var FNV1A_PRIME = 0x01000193;

    // ========================================================================
    // Bit-counting helper (popcount for 32-bit integer)
    // ========================================================================

    /**
     * Population count (Hamming weight) for a 32-bit integer.
     * @param {number} x
     * @return {number}
     */
    function popcount32(x) {
        x = x - ((x >>> 1) & 0x55555555);
        x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
        return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24;
    }

    // ========================================================================
    // FNV-1a hash (32-bit)
    // ========================================================================

    /**
     * FNV-1a 32-bit hash over an array of integers.
     * @param {Array<number>} data
     * @param {number}        len
     * @return {number}
     */
    function fnv1a(data, len) {
        var h = FNV1A_SEED;
        for (var i = 0; i < len; i++) {
            h ^= (data[i] & 0xFFFFFFFF);
            h = Math.imul ? Math.imul(h, FNV1A_PRIME) : imul32(h, FNV1A_PRIME);
            h = h >>> 0; // keep unsigned
        }
        return h >>> 0;
    }

    /**
     * 32-bit integer multiply fallback for environments without Math.imul.
     * @param {number} a
     * @param {number} b
     * @return {number}
     */
    function imul32(a, b) {
        var ah = (a >>> 16) & 0xFFFF;
        var al = a & 0xFFFF;
        var bh = (b >>> 16) & 0xFFFF;
        var bl = b & 0xFFFF;
        return ((al * bl) + (((ah * bl + al * bh) << 16) >>> 0)) | 0;
    }

    // ========================================================================
    // Charge-aware implicit hydrogen count with aromatic pi contribution.
    // Multi-valence elements (S, P, N) expand their valence shell based on
    // charge and bonding context. Aromatic atoms contribute one electron to
    // the pi system, reducing the effective valence available for H.
    // ========================================================================

    /**
     * Standard valence table: lists allowed valences per atomic number.
     * Multi-valence elements have multiple entries (e.g. S: 2, 4, 6).
     * @type {Object<number, Array<number>>}
     */
    var STD_VALENCE = {
        1: [1], 5: [3], 6: [4], 7: [3, 5], 8: [2], 9: [1],
        14: [4], 15: [3, 5], 16: [2, 4, 6], 17: [1], 33: [3, 5],
        34: [2, 4, 6], 35: [1], 53: [1, 3, 5, 7]
    };

    /**
     * Compute implicit hydrogen count for a single atom.
     *
     * @param {Object} graph  MolGraph
     * @param {number} idx    Atom index
     * @return {number}       Implicit H count (>= 0)
     */
    function computeImplicitH(graph, idx) {
        var z = graph.atomicNum ? graph.atomicNum[idx] : (graph.label ? graph.label[idx] : 6);
        var deg = graph.degree ? graph.degree[idx] : (graph.neighbors ? graph.neighbors[idx].length : 0);
        var chg = graph.formalCharge ? graph.formalCharge[idx] : 0;
        var arom = graph.aromatic ? graph.aromatic[idx] : false;

        var vals = STD_VALENCE[z];
        if (!vals) { return 0; }

        // Sum of bond orders (explicit valence)
        var bondSum = 0;
        if (graph.neighbors && graph.bondOrder) {
            var nbs = graph.neighbors[idx];
            for (var k = 0; k < nbs.length; k++) {
                var bo = graph.bondOrder(idx, nbs[k]);
                bondSum += (bo > 0) ? bo : 1;
            }
        } else {
            bondSum = deg;
        }

        // Aromatic atoms contribute 1 electron to pi system
        if (arom) { bondSum += 1; }

        // Charge shifts the target valence: positive charge reduces electrons
        // available, negative charge adds electrons
        var target = bondSum - chg;

        // Find the smallest standard valence >= target
        var bestVal = vals[vals.length - 1];
        for (var vi = 0; vi < vals.length; vi++) {
            if (vals[vi] >= target) { bestVal = vals[vi]; break; }
        }

        var implH = bestVal - bondSum + chg;
        return Math.max(0, implH);
    }

    // ========================================================================
    // Path-based fingerprint generation
    // ========================================================================

    /**
     * Generate a path-based molecular fingerprint.
     *
     * Enumerates all paths up to maxLength bonds via iterative DFS.
     * Each path is hashed (atom labels + bond orders interleaved) using FNV-1a.
     * Both forward and reverse paths are hashed for canonical invariance.
     *
     * @param {Object}  graph      MolGraph with .n, .neighbors[], .label[], .bondOrder()
     * @param {number}  [maxLength=7]  Maximum path length in bonds
     * @param {number}  [fpSize=1024]  Fingerprint size in bits (rounded to multiple of 32)
     * @return {Uint32Array}        Fingerprint as array of 32-bit words
     */
    function generateFingerprint(graph, maxLength, fpSize) {
        if (!graph || typeof graph.n !== 'number') throw new Error('graph is required');
        if (maxLength === undefined || maxLength === null) maxLength = 7;
        if (fpSize === undefined || fpSize === null) fpSize = 1024;
        if (typeof fpSize !== 'number' || fpSize <= 0) throw new Error('fpSize must be a positive number');
        if (typeof maxLength !== 'number' || maxLength < 1) throw new Error('maxLength must be >= 1');

        // Clamp path length to 7 (matching C++ buffer constraints)
        if (maxLength > 7) maxLength = 7;

        // Round fpSize to multiple of 32
        fpSize = ((fpSize + 31) >> 5) << 5;
        var numWords = fpSize >> 5;
        var fp = new Uint32Array(numWords);

        var n = graph.n;
        if (n === 0) return fp;

        var label = graph.label;
        var neighbors = graph.neighbors;

        // Set a bit in the fingerprint
        function setBit(hash) {
            var bit = hash % fpSize;
            fp[bit >> 5] |= (1 << (bit & 31));
        }

        // Iterative DFS stack frames
        // Each frame: { atom, depth, nbIdx, path[], bonds[] }
        var visited = new Uint8Array(n);
        var stackAtom = new Int32Array(16);
        var stackDepth = new Int32Array(16);
        var stackNbIdx = new Int32Array(16);
        // path and bonds stored per-depth level
        var pathLabels = new Int32Array(8);  // atom labels along path
        var pathBonds = new Int32Array(7);   // bond orders between consecutive atoms

        var pathData = new Int32Array(16);  // interleaved: [label0, bond01, label1, ...]
        var revData = new Int32Array(16);

        for (var start = 0; start < n; start++) {
            // Single atom path
            pathData[0] = label[start];
            setBit(fnv1a(pathData, 1));

            // Multi-atom paths via iterative DFS
            var top = 0;
            stackAtom[0] = start;
            stackDepth[0] = 0;
            stackNbIdx[0] = 0;
            pathLabels[0] = label[start];
            visited[start] = 1;

            while (top >= 0) {
                var fAtom = stackAtom[top];
                var fDepth = stackDepth[top];
                var expanded = false;

                if (fDepth < maxLength) {
                    var nbs = neighbors[fAtom];
                    var nbCount = nbs.length;
                    for (var ni = stackNbIdx[top]; ni < nbCount; ni++) {
                        var nb = nbs[ni];
                        if (visited[nb]) continue;

                        var newDepth = fDepth + 1;
                        var bondOrd = graph.bondOrder ? graph.bondOrder(fAtom, nb) : 1;

                        // Build interleaved path: [label0, bond01, label1, bond12, label2, ...]
                        var pathLen = 0;
                        for (var d = 0; d <= fDepth; d++) {
                            pathData[pathLen++] = pathLabels[d];
                            if (d < fDepth) {
                                pathData[pathLen++] = pathBonds[d] + 1000;
                            }
                        }
                        pathData[pathLen++] = bondOrd + 1000;
                        pathData[pathLen++] = label[nb];

                        // Hash forward
                        setBit(fnv1a(pathData, pathLen));

                        // Hash reversed for canonical invariance
                        for (var r = 0; r < pathLen; r++) {
                            revData[r] = pathData[pathLen - 1 - r];
                        }
                        setBit(fnv1a(revData, pathLen));

                        // Push deeper if possible
                        if (newDepth < maxLength && top + 1 < 15) {
                            visited[nb] = 1;
                            stackNbIdx[top] = ni + 1; // resume at next neighbor on backtrack

                            top++;
                            stackAtom[top] = nb;
                            stackDepth[top] = newDepth;
                            stackNbIdx[top] = 0;

                            // Copy path state to new depth level
                            // pathLabels[0..fDepth] already in place from parent
                            pathLabels[newDepth] = label[nb];
                            pathBonds[fDepth] = bondOrd;

                            expanded = true;
                            break;
                        }
                    }
                }

                if (!expanded) {
                    visited[stackAtom[top]] = 0;
                    top--;
                }
            }

            visited[start] = 0;
        }

        return fp;
    }

    // ========================================================================
    // MCS-aware path fingerprint
    // ========================================================================

    /**
     * Hash a single atom for MCS fingerprint.
     * @param {Object}  g
     * @param {boolean} tautAware
     * @param {number}  atom
     * @return {number}
     */
    function mcsAtomHash(g, tautAware, atom) {
        var h = 17;
        var atomicLabel = g.atomicNum ? g.atomicNum[atom] : (g.label ? g.label[atom] : 6);
        var hasTautClass = g.tautomerClass && g.tautomerClass.length > 0 && g.tautomerClass[atom] >= 0;
        if (tautAware && hasTautClass) atomicLabel = 999;
        h = Math.imul ? Math.imul(h, 37) + atomicLabel : imul32(h, 37) + atomicLabel;
        h = h | 0;
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + (g.ring && g.ring[atom] ? 1 : 0);
        h = h | 0;
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + (g.aromatic && g.aromatic[atom] ? 1 : 0);
        h = h | 0;
        var tc = (g.tautomerClass && g.tautomerClass.length > 0) ? g.tautomerClass[atom] : -1;
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + (tautAware ? 0 : tc);
        h = h | 0;
        var deg = g.degree ? g.degree[atom] : (g.neighbors ? g.neighbors[atom].length : 0);
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + deg;
        return h | 0;
    }

    /**
     * Hash a bond for MCS fingerprint.
     * @param {Object} g
     * @param {number} from
     * @param {number} to
     * @return {number}
     */
    function mcsBondHash(g, from, to) {
        var h = 17;
        var bo = g.bondOrder ? g.bondOrder(from, to) : 1;
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + bo;
        h = h | 0;
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + (g.bondInRing && g.bondInRing(from, to) ? 1 : 0);
        h = h | 0;
        h = (Math.imul ? Math.imul(h, 37) : imul32(h, 37)) + (g.bondAromatic && g.bondAromatic(from, to) ? 1 : 0);
        return h | 0;
    }

    /**
     * Hash a path canonically (direction-independent).
     * @param {Object}  g
     * @param {boolean} tautAware
     * @param {Array<number>} path  Atom indices
     * @param {number}  len
     * @return {number}
     */
    function mcsHashPath(g, tautAware, path, len) {
        var fwd = 17, rev = 17;
        for (var i = 0; i < len; i++) {
            fwd = (Math.imul ? Math.imul(fwd, 31) : imul32(fwd, 31)) + mcsAtomHash(g, tautAware, path[i]);
            fwd = fwd | 0;
            rev = (Math.imul ? Math.imul(rev, 31) : imul32(rev, 31)) + mcsAtomHash(g, tautAware, path[len - 1 - i]);
            rev = rev | 0;
            if (i < len - 1) {
                fwd = (Math.imul ? Math.imul(fwd, 31) : imul32(fwd, 31)) + mcsBondHash(g, path[i], path[i + 1]);
                fwd = fwd | 0;
                rev = (Math.imul ? Math.imul(rev, 31) : imul32(rev, 31)) + mcsBondHash(g, path[len - 1 - i], path[len - 2 - i]);
                rev = rev | 0;
            }
        }
        fwd = fwd & 0x7FFFFFFF;
        rev = rev & 0x7FFFFFFF;
        return fwd < rev ? fwd : rev;
    }

    /**
     * Generate an MCS-aware path fingerprint.
     *
     * Encodes element, ring, aromatic, tautomer class, degree per atom;
     * bond order, ring, aromatic per bond. Paths are hashed canonically.
     *
     * @param {Object}  graph
     * @param {number}  [pathLength=7]
     * @param {number}  [fpSize=2048]
     * @return {Uint32Array}
     */
    function generateMcsFingerprint(graph, pathLength, fpSize) {
        if (!graph || typeof graph.n !== 'number') throw new Error('graph is required');
        if (pathLength === undefined || pathLength === null) pathLength = 7;
        if (fpSize === undefined || fpSize === null) fpSize = 2048;
        if (typeof fpSize !== 'number' || fpSize <= 0) throw new Error('fpSize must be a positive number');
        if (typeof pathLength !== 'number' || pathLength < 1) throw new Error('pathLength must be >= 1');

        fpSize = ((fpSize + 31) >> 5) << 5;
        var numWords = fpSize >> 5;
        var fp = new Uint32Array(numWords);
        var n = graph.n;
        if (n === 0) return fp;

        // Identify heavy atoms
        var isHeavy = new Uint8Array(n);
        var heavyCount = 0;
        for (var i = 0; i < n; i++) {
            var z = graph.atomicNum ? graph.atomicNum[i] : (graph.label ? graph.label[i] : 6);
            isHeavy[i] = (z !== 1) ? 1 : 0;
            if (isHeavy[i]) heavyCount++;
        }
        if (heavyCount === 0) return fp;

        var tautAware = graph.tautomerClass && graph.tautomerClass.length > 0;

        function mcSetBit(hash) {
            var bit = ((hash % fpSize) + fpSize) % fpSize; // handle negatives
            fp[bit >> 5] |= (1 << (bit & 31));
        }

        var pathBuf = new Int32Array(pathLength + 1);
        var visited = new Uint8Array(n);

        // Recursive path enumeration
        function enumerate(depth, maxDepth) {
            var cur = pathBuf[depth - 1];
            var nbs = graph.neighbors[cur];
            for (var k = 0; k < nbs.length; k++) {
                var nb = nbs[k];
                if (visited[nb] || !isHeavy[nb]) continue;
                pathBuf[depth] = nb;
                visited[nb] = 1;

                // Convert typed array slice to regular array for hashing
                var pathArr = [];
                for (var pi = 0; pi <= depth; pi++) pathArr.push(pathBuf[pi]);
                mcSetBit(mcsHashPath(graph, tautAware, pathArr, depth + 1));

                if (depth < maxDepth) {
                    enumerate(depth + 1, maxDepth);
                }
                visited[nb] = 0;
            }
        }

        for (var start = 0; start < n; start++) {
            if (!isHeavy[start]) continue;
            pathBuf[0] = start;
            visited[start] = 1;
            mcSetBit(mcsHashPath(graph, tautAware, [start], 1));
            enumerate(1, pathLength);
            visited[start] = 0;
        }

        return fp;
    }

    // ========================================================================
    // Tanimoto coefficient
    // ========================================================================

    /**
     * Compute Tanimoto coefficient between two fingerprints.
     *
     * @param {Uint32Array} fp1
     * @param {Uint32Array} fp2
     * @return {number}  Float in [0, 1]
     */
    function tanimotoFingerprint(fp1, fp2) {
        if (!fp1 || !fp2 || fp1.length === 0 || fp2.length === 0) return 0.0;

        var minWords = fp1.length < fp2.length ? fp1.length : fp2.length;
        var maxWords = fp1.length > fp2.length ? fp1.length : fp2.length;
        var longer = fp1.length > fp2.length ? fp1 : fp2;

        var andBits = 0, orBits = 0;
        for (var i = 0; i < minWords; i++) {
            andBits += popcount32((fp1[i] & fp2[i]) >>> 0);
            orBits  += popcount32((fp1[i] | fp2[i]) >>> 0);
        }
        // Extra words in the longer fingerprint contribute only to OR
        for (i = minWords; i < maxWords; i++) {
            orBits += popcount32(longer[i] >>> 0);
        }

        return orBits === 0 ? 1.0 : andBits / orBits;
    }

    /**
     * Check if fingerprint a is a subset of fingerprint b.
     * Used for substructure pre-screening: if query is substructure of target,
     * then query FP must be a subset of target FP.
     *
     * @param {Uint32Array} query
     * @param {Uint32Array} target
     * @return {boolean}
     */
    function fingerprintSubset(query, target) {
        if (!query || query.length === 0) return true;
        var qWords = query.length;
        var tWords = target ? target.length : 0;
        for (var i = 0; i < qWords; i++) {
            var qw = query[i];
            var tw = (i < tWords) ? target[i] : 0;
            if ((qw & tw) !== qw) return false;
        }
        return true;
    }

    // ========================================================================
    // Batch Substructure: 1 query vs N targets
    // ========================================================================

    /**
     * Check if query is a substructure of each target molecule.
     *
     * @param {Object}        query    MolGraph
     * @param {Array<Object>} targets  Array of MolGraphs
     * @param {Object}        [opts]   ChemOptions for SMSDVF2.isSubstructure
     * @return {Array<number>}         Indices of targets where query is a substructure
     */
    function batchSubstructure(query, targets, opts) {
        if (!query || typeof query.n !== 'number') throw new Error('query graph is required');
        if (!targets || !targets.length) return [];
        if (query.n === 0) return [];
        if (!opts) opts = {};

        var VF2 = window.SMSDVF2;
        if (!VF2) throw new Error('SMSDVF2 not loaded');

        var hits = [];
        for (var i = 0; i < targets.length; i++) {
            if (targets[i].n >= query.n) {
                if (VF2.isSubstructure(query, targets[i], opts)) {
                    hits.push(i);
                }
            }
        }
        return hits;
    }

    // ========================================================================
    // Batch MCS: 1 query vs N targets
    // ========================================================================

    /**
     * Find MCS between query and each target. Returns results sorted by
     * Tanimoto score descending.
     *
     * @param {Object}        query      MolGraph
     * @param {Array<Object>} targets    Array of MolGraphs
     * @param {Object}        [chemOpts] ChemOptions
     * @param {Object}        [mcsOpts]  McsOptions
     * @return {Array<Object>} Array of { targetIdx, mapping, size, tanimoto }
     */
    function batchMCS(query, targets, chemOpts, mcsOpts) {
        if (!query || typeof query.n !== 'number') throw new Error('query graph is required');
        if (!targets || !targets.length) return [];
        if (query.n === 0) return [];
        if (!chemOpts) chemOpts = {};
        if (!mcsOpts) mcsOpts = {};
        if (mcsOpts.timeoutMs !== undefined && (typeof mcsOpts.timeoutMs !== 'number' || mcsOpts.timeoutMs <= 0)) {
            throw new Error('timeoutMs must be a positive number');
        }
        if (mcsOpts.maxResults !== undefined && (typeof mcsOpts.maxResults !== 'number' || mcsOpts.maxResults < 1)) {
            throw new Error('maxResults must be >= 1');
        }

        var MCS = window.SMSDMCS;
        if (!MCS) throw new Error('SMSDMCS not loaded');

        var results = [];
        for (var i = 0; i < targets.length; i++) {
            if (targets[i].n === 0) continue;
            var res = MCS.findMCS(query, targets[i], chemOpts, mcsOpts);
            // findMCS returns { mapping, size, tanimoto, overlap }
            if (res && res.mapping && res.size > 0) {
                var denom = (query.n + targets[i].n - res.size);
                var tanimoto = denom > 0 ? (res.size / denom) : 0.0;
                results.push({
                    targetIdx: i,
                    mapping: res.mapping,
                    size: res.size,
                    tanimoto: tanimoto
                });
            }
        }

        // Sort by tanimoto descending
        results.sort(function(a, b) { return b.tanimoto - a.tanimoto; });
        return results;
    }

    // ========================================================================
    // Screen and Match: RASCAL pre-screen + exact MCS on hits
    // ========================================================================

    /**
     * Two-phase pipeline:
     *   Phase 1: RASCAL similarity upper-bound screening
     *   Phase 2: Exact MCS only on targets passing threshold
     *
     * @param {Object}        query        MolGraph
     * @param {Array<Object>} targets      Array of MolGraphs
     * @param {number}        [threshold=0.3] RASCAL threshold
     * @param {Object}        [chemOpts]   ChemOptions
     * @param {Object}        [mcsOpts]    McsOptions
     * @return {Array<Object>} Array of { targetIdx, mapping, size, tanimoto }
     */
    function screenAndMatch(query, targets, threshold, chemOpts, mcsOpts) {
        if (!query || typeof query.n !== 'number') throw new Error('query graph is required');
        if (!targets || !targets.length) return [];
        if (query.n === 0) return [];
        if (threshold === undefined || threshold === null) threshold = 0.3;
        if (typeof threshold !== 'number' || threshold < 0 || threshold > 1) {
            throw new Error('threshold must be a number between 0 and 1');
        }
        if (!chemOpts) chemOpts = {};
        if (!mcsOpts) mcsOpts = {};

        var MCS = window.SMSDMCS;
        if (!MCS) throw new Error('SMSDMCS not loaded');

        var N = targets.length;

        // Phase 1: RASCAL screening
        var hits = [];
        for (var i = 0; i < N; i++) {
            if (targets[i].n === 0) continue;
            var score = 0.0;
            if (MCS.similarityUpperBound) {
                score = MCS.similarityUpperBound(query, targets[i]);
            } else {
                // Fallback: no screening, pass everything
                score = 1.0;
            }
            if (score >= threshold) {
                hits.push(i);
            }
        }

        if (hits.length === 0) return [];

        // Phase 2: Exact MCS on hits only
        var results = [];
        for (var h = 0; h < hits.length; h++) {
            var idx = hits[h];
            var res = MCS.findMCS(query, targets[idx], chemOpts, mcsOpts);
            // findMCS returns { mapping, size, tanimoto, overlap }
            if (res && res.mapping && res.size > 0) {
                var denom = (query.n + targets[idx].n - res.size);
                var tanimoto = denom > 0 ? (res.size / denom) : 0.0;
                results.push({
                    targetIdx: idx,
                    mapping: res.mapping,
                    size: res.size,
                    tanimoto: tanimoto
                });
            }
        }

        // Sort by tanimoto descending
        results.sort(function(a, b) { return b.tanimoto - a.tanimoto; });
        return results;
    }

    // ========================================================================
    // Fingerprint format conversions
    // ========================================================================

    /**
     * Convert a fingerprint to a hexadecimal string (lowercase, zero-padded).
     * @param {Uint32Array} fp
     * @return {string}
     */
    function toHex(fp) {
        var out = '';
        for (var i = 0; i < fp.length; i++) {
            var word = fp[i] >>> 0;
            var hex = word.toString(16);
            while (hex.length < 8) hex = '0' + hex;
            out += hex;
        }
        return out;
    }

    /**
     * Parse a hexadecimal string back into a fingerprint.
     * @param {string} hex
     * @return {Uint32Array}
     */
    function fromHex(hex) {
        var words = Math.ceil(hex.length / 8);
        var fp = new Uint32Array(words);
        for (var i = 0; i < words; i++) {
            var start = i * 8;
            var chunk = hex.substring(start, start + 8);
            fp[i] = parseInt(chunk, 16) >>> 0;
        }
        return fp;
    }

    /**
     * Convert a fingerprint to a binary string of '0' and '1' characters.
     * @param {Uint32Array} fp
     * @param {number}      fpSize  Fingerprint size in bits
     * @return {string}
     */
    function toBinaryString(fp, fpSize) {
        var out = '';
        for (var i = 0; i < fpSize; i++) {
            out += (fp[i >> 5] >>> (i & 31)) & 1 ? '1' : '0';
        }
        return out;
    }

    // ========================================================================
    // ECFP (Extended Connectivity Fingerprint) — count-based
    //
    // Morgan-style iterative radius expansion.  Each atom starts with an
    // initial invariant (atomic number, degree, charge, aromatic, ring, H
    // count).  At each iteration the neighbourhood hashes are folded in.
    // The result is an Object mapping identifier -> count.
    // ========================================================================

    /**
     * Generate a count-based ECFP fingerprint.
     *
     * @param {Object}  graph      MolGraph with .n, .neighbors[], .atomicNum[], etc.
     * @param {number}  [radius=2] Morgan radius (2 = ECFP4)
     * @return {Object}            Map of {identifier: count}
     */
    function ecfpCounts(graph, radius) {
        if (!graph || typeof graph.n !== 'number') throw new Error('graph is required');
        if (radius === undefined || radius === null) radius = 2;
        if (typeof radius !== 'number' || radius < 0) throw new Error('radius must be >= 0');
        var n = graph.n;
        if (n === 0) return {};

        // Initial atom invariants
        var current = new Int32Array(n);
        var data = new Int32Array(8);
        for (var i = 0; i < n; i++) {
            data[0] = graph.atomicNum ? graph.atomicNum[i] : (graph.label ? graph.label[i] : 6);
            data[1] = graph.degree ? graph.degree[i] : (graph.neighbors ? graph.neighbors[i].length : 0);
            data[2] = graph.formalCharge ? graph.formalCharge[i] : 0;
            data[3] = graph.aromatic ? (graph.aromatic[i] ? 1 : 0) : 0;
            data[4] = graph.ring ? (graph.ring[i] ? 1 : 0) : 0;
            // Charge-aware implicit H with aromatic pi contribution
            data[5] = computeImplicitH(graph, i);
            current[i] = fnv1a(data, 6);
        }

        var counts = {};
        // Collect radius-0 identifiers
        for (i = 0; i < n; i++) {
            var id0 = current[i] >>> 0;
            counts[id0] = (counts[id0] || 0) + 1;
        }

        // Iterative expansion
        var next = new Int32Array(n);
        var nbHash = new Int32Array(128);
        for (var r = 1; r <= radius; r++) {
            for (i = 0; i < n; i++) {
                var nbs = graph.neighbors[i];
                var deg2 = nbs.length;
                // Sort neighbour hashes for canonical invariance
                for (var k = 0; k < deg2; k++) {
                    nbHash[k] = current[nbs[k]];
                }
                // Insertion sort
                for (var a = 1; a < deg2; a++) {
                    var tmp = nbHash[a], b = a - 1;
                    while (b >= 0 && nbHash[b] > tmp) { nbHash[b + 1] = nbHash[b]; b--; }
                    nbHash[b + 1] = tmp;
                }
                // Hash: current atom + sorted neighbours + radius
                data[0] = current[i];
                data[1] = r;
                for (k = 0; k < deg2; k++) data[k + 2] = nbHash[k];
                next[i] = fnv1a(data, deg2 + 2);
            }
            // Swap current <-> next (temp var swap, avoids element-by-element copy)
            var swapTmp = current;
            current = next;
            next = swapTmp;
            // Collect counts from new current
            for (i = 0; i < n; i++) {
                var id = current[i] >>> 0;
                counts[id] = (counts[id] || 0) + 1;
            }
        }

        return counts;
    }

    // ========================================================================
    // FCFP (Functional Class Fingerprint) — count-based
    //
    // Like ECFP but initial invariants are pharmacophore feature classes:
    // donor, acceptor, positive ionisable, negative ionisable, aromatic,
    // hydrophobic — instead of atomic properties.
    // ========================================================================

    /**
     * Pharmacophore feature classification for FCFP.
     * Returns a 6-bit feature mask per atom.
     * @param {Object} graph
     * @return {Uint8Array}
     */
    function pharmacophoreFeatures(graph) {
        var n = graph.n;
        var feat = new Uint8Array(n);
        var SG = window.SMSDGraph;
        for (var i = 0; i < n; i++) {
            var z = graph.atomicNum ? graph.atomicNum[i] : 6;
            var deg = graph.degree ? graph.degree[i] : (graph.neighbors ? graph.neighbors[i].length : 0);
            var arom = graph.aromatic ? graph.aromatic[i] : false;
            var chg = graph.formalCharge ? graph.formalCharge[i] : 0;
            var inRing = graph.ring ? graph.ring[i] : false;
            var f = 0;
            // Donor: N-H, O-H, S-H (N or O or S with implicit H)
            // Charge-aware implicit H with aromatic pi contribution
            var implicitH = computeImplicitH(graph, i);
            if ((z === 7 || z === 8 || z === 16) && implicitH > 0) f |= 1;
            // Acceptor: N or O with lone pair (pyridine N is acceptor)
            if (z === 7 || z === 8) {
                // Pyridine N: aromatic ring N with degree 2 — H-bond acceptor
                if (SG && SG.isPyridineN && SG.isPyridineN(graph, i)) {
                    f |= 2;
                } else if (z === 8 || (z === 7 && deg < 4)) {
                    f |= 2;
                }
            }
            // Positive ionisable: N with positive charge or basic N (exclude aniline)
            if (z === 7 && (chg > 0 || (!arom && implicitH > 0))) {
                var isAniline = SG && SG.isAnilineN && SG.isAnilineN(graph, i);
                if (!isAniline) f |= 4;
            }
            // Negative ionisable: O with negative charge or acidic O-H (exclude ester O)
            if ((z === 8 && chg < 0) || (z === 8 && implicitH > 0)) {
                var isEster = SG && SG.isEsterO && SG.isEsterO(graph, i);
                if (!isEster) f |= 8;
            }
            // Aromatic
            if (arom) f |= 16;
            // Hydrophobic: C or S in non-polar context
            if ((z === 6 && !arom && deg <= 2) || (z === 16 && deg <= 2 && !chg)) f |= 32;
            feat[i] = f;
        }
        return feat;
    }

    /**
     * Generate a count-based FCFP fingerprint.
     *
     * @param {Object}  graph      MolGraph
     * @param {number}  [radius=2] Morgan radius (2 = FCFP4)
     * @return {Object}            Map of {identifier: count}
     */
    function fcfpCounts(graph, radius) {
        if (!graph || typeof graph.n !== 'number') throw new Error('graph is required');
        if (radius === undefined || radius === null) radius = 2;
        if (typeof radius !== 'number' || radius < 0) throw new Error('radius must be >= 0');
        var n = graph.n;
        if (n === 0) return {};

        var feat = pharmacophoreFeatures(graph);

        // Initial invariant from pharmacophore class
        var current = new Int32Array(n);
        var data = new Int32Array(8);
        for (var i = 0; i < n; i++) {
            data[0] = feat[i];
            current[i] = fnv1a(data, 1);
        }

        var counts = {};
        for (i = 0; i < n; i++) {
            var id0 = current[i] >>> 0;
            counts[id0] = (counts[id0] || 0) + 1;
        }

        // Iterative expansion (same as ECFP)
        var next = new Int32Array(n);
        var nbHash = new Int32Array(128);
        for (var r = 1; r <= radius; r++) {
            for (i = 0; i < n; i++) {
                var nbs = graph.neighbors[i];
                var deg = nbs.length;
                for (var k = 0; k < deg; k++) nbHash[k] = current[nbs[k]];
                for (var a = 1; a < deg; a++) {
                    var tmp = nbHash[a], b2 = a - 1;
                    while (b2 >= 0 && nbHash[b2] > tmp) { nbHash[b2 + 1] = nbHash[b2]; b2--; }
                    nbHash[b2 + 1] = tmp;
                }
                data[0] = current[i];
                data[1] = r;
                for (k = 0; k < deg; k++) data[k + 2] = nbHash[k];
                next[i] = fnv1a(data, deg + 2);
            }
            // Swap current <-> next (temp var swap, avoids element-by-element copy)
            var swapTmp2 = current;
            current = next;
            next = swapTmp2;
            for (i = 0; i < n; i++) {
                var id = current[i] >>> 0;
                counts[id] = (counts[id] || 0) + 1;
            }
        }

        return counts;
    }

    // ========================================================================
    // Topological Torsion fingerprint
    //
    // 4-atom linear path fingerprint.  Each atom is encoded as
    // (atomicNum, degree, piElectrons, hCount).  All 4-atom paths are
    // enumerated and hashed.  Result is count-based.
    // ========================================================================

    /**
     * Generate a topological torsion fingerprint (count-based).
     *
     * @param {Object}  graph  MolGraph
     * @return {Object}        Map of {identifier: count}
     */
    function topologicalTorsion(graph) {
        if (!graph || typeof graph.n !== 'number') throw new Error('graph is required');
        var n = graph.n;
        if (n === 0) return {};

        // Precompute atom typing: (atomicNum, degree, piElectrons, hCount)
        var atomType = new Int32Array(n * 4);
        for (var i = 0; i < n; i++) {
            var z = graph.atomicNum ? graph.atomicNum[i] : 6;
            var deg = graph.degree ? graph.degree[i] : graph.neighbors[i].length;
            var hCount = computeImplicitH(graph, i);
            // Pi electrons: sum of (bondOrder - 1) for all bonds
            var piE = 0;
            var nbs = graph.neighbors[i];
            for (var k = 0; k < nbs.length; k++) {
                var bo = graph.bondOrder ? graph.bondOrder(i, nbs[k]) : 1;
                if (bo > 1) piE += (bo - 1);
            }
            atomType[i * 4] = z;
            atomType[i * 4 + 1] = deg;
            atomType[i * 4 + 2] = piE;
            atomType[i * 4 + 3] = hCount;
        }

        var counts = {};
        var pathData = new Int32Array(16); // 4 atoms x 4 features

        // Enumerate all 4-atom paths
        for (var a0 = 0; a0 < n; a0++) {
            var nbs0 = graph.neighbors[a0];
            for (var j1 = 0; j1 < nbs0.length; j1++) {
                var a1 = nbs0[j1];
                if (a1 === a0) continue;
                var nbs1 = graph.neighbors[a1];
                for (var j2 = 0; j2 < nbs1.length; j2++) {
                    var a2 = nbs1[j2];
                    if (a2 === a0 || a2 === a1) continue;
                    var nbs2 = graph.neighbors[a2];
                    for (var j3 = 0; j3 < nbs2.length; j3++) {
                        var a3 = nbs2[j3];
                        if (a3 === a0 || a3 === a1 || a3 === a2) continue;

                        // Canonical direction: use smaller endpoint first
                        var forward = true;
                        if (atomType[a0 * 4] > atomType[a3 * 4]) {
                            forward = false;
                        } else if (atomType[a0 * 4] === atomType[a3 * 4]) {
                            if (atomType[a0 * 4 + 1] > atomType[a3 * 4 + 1]) {
                                forward = false;
                            } else if (atomType[a0 * 4 + 1] === atomType[a3 * 4 + 1]) {
                                if (a0 > a3) forward = false;
                            }
                        }
                        if (!forward) continue; // skip reverse to avoid double-counting

                        // Pack atom types
                        var idx = 0;
                        var atoms4 = [a0, a1, a2, a3];
                        for (var q = 0; q < 4; q++) {
                            var ai = atoms4[q];
                            pathData[idx++] = atomType[ai * 4];
                            pathData[idx++] = atomType[ai * 4 + 1];
                            pathData[idx++] = atomType[ai * 4 + 2];
                            pathData[idx++] = atomType[ai * 4 + 3];
                        }

                        var h = fnv1a(pathData, 16) >>> 0;
                        counts[h] = (counts[h] || 0) + 1;
                    }
                }
            }
        }

        return counts;
    }

    // ========================================================================
    // Dice similarity — binary fingerprints
    // ========================================================================

    /**
     * Compute Dice similarity between two binary fingerprints.
     * Dice = 2 * |A AND B| / (|A| + |B|)
     *
     * @param {Uint32Array} fp1
     * @param {Uint32Array} fp2
     * @return {number}
     */
    function diceSimilarity(fp1, fp2) {
        if (!fp1 || !fp2 || fp1.length === 0 || fp2.length === 0) return 0.0;

        var minWords = fp1.length < fp2.length ? fp1.length : fp2.length;
        var maxWords = fp1.length > fp2.length ? fp1.length : fp2.length;
        var longer = fp1.length > fp2.length ? fp1 : fp2;

        var andBits = 0, aBits = 0, bBits = 0;
        for (var i = 0; i < minWords; i++) {
            andBits += popcount32((fp1[i] & fp2[i]) >>> 0);
            aBits   += popcount32(fp1[i] >>> 0);
            bBits   += popcount32(fp2[i] >>> 0);
        }
        for (i = minWords; i < maxWords; i++) {
            if (fp1.length > fp2.length) {
                aBits += popcount32(longer[i] >>> 0);
            } else {
                bBits += popcount32(longer[i] >>> 0);
            }
        }

        var denom = aBits + bBits;
        return denom === 0 ? 1.0 : (2 * andBits) / denom;
    }

    // ========================================================================
    // Cosine similarity — binary fingerprints
    // ========================================================================

    /**
     * Compute Cosine similarity between two binary fingerprints.
     * Cosine = |A AND B| / sqrt(|A| * |B|)
     *
     * @param {Uint32Array} fp1
     * @param {Uint32Array} fp2
     * @return {number}
     */
    function cosineSimilarity(fp1, fp2) {
        if (!fp1 || !fp2 || fp1.length === 0 || fp2.length === 0) return 0.0;

        var minWords = fp1.length < fp2.length ? fp1.length : fp2.length;
        var maxWords = fp1.length > fp2.length ? fp1.length : fp2.length;
        var longer = fp1.length > fp2.length ? fp1 : fp2;

        var andBits = 0, aBits = 0, bBits = 0;
        for (var i = 0; i < minWords; i++) {
            andBits += popcount32((fp1[i] & fp2[i]) >>> 0);
            aBits   += popcount32(fp1[i] >>> 0);
            bBits   += popcount32(fp2[i] >>> 0);
        }
        for (i = minWords; i < maxWords; i++) {
            if (fp1.length > fp2.length) {
                aBits += popcount32(longer[i] >>> 0);
            } else {
                bBits += popcount32(longer[i] >>> 0);
            }
        }

        var denom = Math.sqrt(aBits * bBits);
        return denom === 0 ? 1.0 : andBits / denom;
    }

    // ========================================================================
    // Soergel distance — 1 - Tanimoto
    // ========================================================================

    /**
     * Compute Soergel distance between two binary fingerprints.
     * Soergel = 1 - Tanimoto
     *
     * @param {Uint32Array} fp1
     * @param {Uint32Array} fp2
     * @return {number}
     */
    function soergelDistance(fp1, fp2) {
        return 1.0 - tanimotoFingerprint(fp1, fp2);
    }

    // ========================================================================
    // Count-vector similarity functions
    //
    // For count-based fingerprints (Objects mapping id -> count).
    // ========================================================================

    /**
     * Tanimoto similarity for count vectors.
     * T = sum(min(a,b)) / (sum(a) + sum(b) - sum(min(a,b)))
     *
     * @param {Object} cv1  Map {id: count}
     * @param {Object} cv2  Map {id: count}
     * @return {number}
     */
    function countTanimoto(cv1, cv2) {
        if (!cv1 || !cv2) return 0.0;
        var sumMin = 0, sumA = 0, sumB = 0;
        var key;
        for (key in cv1) {
            if (cv1.hasOwnProperty(key)) {
                var a = cv1[key];
                sumA += a;
                if (cv2.hasOwnProperty(key)) {
                    sumMin += (a < cv2[key] ? a : cv2[key]);
                }
            }
        }
        for (key in cv2) {
            if (cv2.hasOwnProperty(key)) {
                sumB += cv2[key];
            }
        }
        var denom = sumA + sumB - sumMin;
        return denom === 0 ? 1.0 : sumMin / denom;
    }

    /**
     * Dice similarity for count vectors.
     * D = 2 * sum(min(a,b)) / (sum(a) + sum(b))
     *
     * @param {Object} cv1
     * @param {Object} cv2
     * @return {number}
     */
    function countDice(cv1, cv2) {
        if (!cv1 || !cv2) return 0.0;
        var sumMin = 0, sumA = 0, sumB = 0;
        var key;
        for (key in cv1) {
            if (cv1.hasOwnProperty(key)) {
                var a = cv1[key];
                sumA += a;
                if (cv2.hasOwnProperty(key)) {
                    sumMin += (a < cv2[key] ? a : cv2[key]);
                }
            }
        }
        for (key in cv2) {
            if (cv2.hasOwnProperty(key)) {
                sumB += cv2[key];
            }
        }
        var denom = sumA + sumB;
        return denom === 0 ? 1.0 : (2 * sumMin) / denom;
    }

    /**
     * Cosine similarity for count vectors.
     * C = sum(a[i]*b[i]) / (sqrt(sum(a[i]^2)) * sqrt(sum(b[i]^2)))
     *
     * @param {Object} cv1
     * @param {Object} cv2
     * @return {number}
     */
    function countCosine(cv1, cv2) {
        if (!cv1 || !cv2) return 0.0;
        var dot = 0, normA = 0, normB = 0;
        var key;
        for (key in cv1) {
            if (cv1.hasOwnProperty(key)) {
                var a = cv1[key];
                normA += a * a;
                if (cv2.hasOwnProperty(key)) {
                    dot += a * cv2[key];
                }
            }
        }
        for (key in cv2) {
            if (cv2.hasOwnProperty(key)) {
                var bv = cv2[key];
                normB += bv * bv;
            }
        }
        var denom = Math.sqrt(normA) * Math.sqrt(normB);
        return denom === 0 ? 1.0 : dot / denom;
    }

    // ========================================================================
    // fingerprintFromSmiles — parse SMILES, build graph, generate fingerprint
    // ========================================================================

    /**
     * One-call fingerprint generation from a SMILES string.
     * @param {string} smiles  The SMILES string
     * @param {string} [type='path']  Fingerprint type: 'path'|'ecfp'|'fcfp'|'torsion'|'mcs'
     * @param {Object} [opts]  Options: { radius, maxLength, fpSize }
     * @returns {Uint32Array|Object}  Binary FP (path/mcs) or count map (ecfp/fcfp/torsion)
     */
    function fingerprintFromSmiles(smiles, type, opts) {
        if (!window.SmilesParser || !window.SmilesParser.parse) {
            throw new Error('SmilesParser not available');
        }
        if (!window.SMSDGraph || !window.SMSDGraph.SMSDGraph) {
            throw new Error('SMSDGraph not available');
        }

        var mol = window.SmilesParser.parse(smiles);
        if (!mol) { throw new Error('Failed to parse SMILES'); }
        if (mol.parseErrors && mol.parseErrors.length > 0) {
            throw new Error('SMILES parse error: ' + mol.parseErrors[0]);
        }

        var SG = window.SMSDGraph.SMSDGraph;
        var graph = new SG(mol);

        if (!type) { type = 'path'; }
        if (!opts) { opts = {}; }

        switch (type) {
            case 'path':
                return generateFingerprint(graph, opts.maxLength, opts.fpSize);
            case 'ecfp':
                return ecfpCounts(graph, opts.radius);
            case 'fcfp':
                return fcfpCounts(graph, opts.radius);
            case 'torsion':
                return topologicalTorsion(graph);
            case 'mcs':
                return generateMcsFingerprint(graph, opts.maxLength, opts.fpSize);
            default:
                throw new Error('Unknown fingerprint type: ' + type);
        }
    }

    // ========================================================================
    // countsToArray — convert sparse count map to dense Float64Array
    // ========================================================================

    /**
     * Convert a sparse count map {hash: count} to a dense Float64Array.
     * Hashes are folded into the target size via modulo.
     * @param {Object} counts   Sparse count map from ecfpCounts/fcfpCounts/topologicalTorsion
     * @param {number} fpSize   Target array size (default 2048)
     * @returns {Float64Array}  Dense count vector
     */
    function countsToArray(counts, fpSize) {
        if (!fpSize || fpSize <= 0) { fpSize = 2048; }
        var arr = new Float64Array(fpSize);
        var key;
        for (key in counts) {
            if (counts.hasOwnProperty(key)) {
                // Use absolute value of hash for folding
                var hash = +key;
                var idx = ((hash % fpSize) + fpSize) % fpSize;
                arr[idx] += counts[key];
            }
        }
        return arr;
    }

    // ========================================================================
    // batchSimilarity — compare query fingerprint against array of targets
    // ========================================================================

    /**
     * Compare a query fingerprint against an array of target fingerprints.
     * Returns sorted (descending) array of { index, tanimoto } pairs.
     * @param {Uint32Array} queryFP     Query binary fingerprint
     * @param {Uint32Array[]} targetFPs Array of target binary fingerprints
     * @returns {Array<{index: number, tanimoto: number}>}  Sorted results
     */
    function batchSimilarity(queryFP, targetFPs) {
        if (!queryFP || !targetFPs) { return []; }
        var results = [];
        for (var i = 0; i < targetFPs.length; i++) {
            var sim = tanimotoFingerprint(queryFP, targetFPs[i]);
            results.push({ index: i, tanimoto: sim });
        }
        results.sort(function(a, b) { return b.tanimoto - a.tanimoto; });
        return results;
    }

    // ========================================================================
    // Export
    // ========================================================================

    window.SMSDBatch = {
        // Batch operations
        batchSubstructure: batchSubstructure,
        batchMCS: batchMCS,
        screenAndMatch: screenAndMatch,

        // Fingerprint generation
        generateFingerprint: generateFingerprint,
        generateMcsFingerprint: generateMcsFingerprint,
        fingerprintFromSmiles: fingerprintFromSmiles,

        // Circular fingerprints (count-based)
        ecfpCounts: ecfpCounts,
        fcfpCounts: fcfpCounts,

        // Topological Torsion fingerprint (count-based)
        topologicalTorsion: topologicalTorsion,

        // Binary fingerprint comparison
        tanimotoFingerprint: tanimotoFingerprint,
        diceSimilarity: diceSimilarity,
        cosineSimilarity: cosineSimilarity,
        soergelDistance: soergelDistance,
        fingerprintSubset: fingerprintSubset,
        batchSimilarity: batchSimilarity,

        // Count-vector similarity
        countTanimoto: countTanimoto,
        countDice: countDice,
        countCosine: countCosine,

        // Format conversions
        toHex: toHex,
        fromHex: fromHex,
        toBinaryString: toBinaryString,
        countsToArray: countsToArray,

        // Expose for testing
        _popcount32: popcount32,
        _fnv1a: fnv1a,
        _pharmacophoreFeatures: pharmacophoreFeatures
    };

})();
