/**
 * SMSDMCS.js — Maximum Common Substructure (MCS) engine
 *
 * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK.
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * Ported from SMSD upstream, mcs.hpp
 *
 * Coverage-driven funnel with 7+ levels:
 *   L0    Label-frequency + degree-sequence upper bound
 *   L0.25 Linear chain fast-path (PEG, polyethylene)
 *   L0.5  Tree fast-path (branched polymers, dendrimers, glycogen)
 *   L0.75 Greedy probe
 *   L1    Substructure containment check (via SMSDVF2)
 *   L1.25 Augmenting path refinement
 *   L1.5  Seed-and-extend (bond growth)
 *   L2    McSplit partition refinement (bidirectional, RRSplit maximality)
 *   L3    Bron-Kerbosch + Tomita pivoting + k-core reduction
 *   L4    McGregor extension (bond-grow + atom-frontier DFS + forced assignment)
 *   L5    Extra seeds (only when McSplit and BK disagree)
 *
 * Dependencies: window.SMSDGraph, window.SMSDVF2
 */
(function(global) {
    'use strict';

    // -----------------------------------------------------------------------
    // McsOptions defaults
    // -----------------------------------------------------------------------
    var McsOptions = {
        induced: false,
        connectedOnly: true,
        disconnectedMCS: false,
        maximizeBonds: false,
        minFragmentSize: 1,
        maxFragments: 999999,
        timeoutMs: 10000,  // default; overridden by adaptive timeout in findMCS
        extraSeeds: true,
        seedNeighborhoodRadius: 2,
        seedMaxAnchors: 12,
        useTwoHopNLF: true,
        useThreeHopNLF: false,
        templateFuzzyAtoms: 2,
        reactionAware: false  // v6.4: reaction-aware MCS with charge relaxation
    };

    var MAX_NODE_LIMIT = 500000;
    var SEED_EXTEND_NODE_LIMIT = 100000;
    var GREEDY_PROBE_MAX_SIZE = 40;
    var SEED_EXTEND_MAX_ATOMS = 50;
    var BK_SKIP_RATIO = 0.8;

    // -----------------------------------------------------------------------
    // Utility helpers
    // -----------------------------------------------------------------------
    function mergeOpts(defaults, user) {
        var result = {};
        var k;
        for (k in defaults) {
            if (defaults.hasOwnProperty(k)) {
                result[k] = defaults[k];
            }
        }
        if (user) {
            for (k in user) {
                if (user.hasOwnProperty(k)) {
                    result[k] = user[k];
                }
            }
        }
        return result;
    }

    function mapSize(m) {
        var n = 0, k;
        for (k in m) {
            if (m.hasOwnProperty(k)) { n++; }
        }
        return n;
    }

    function mapCopy(m) {
        var r = {}, k;
        for (k in m) {
            if (m.hasOwnProperty(k)) { r[k] = m[k]; }
        }
        return r;
    }

    function mapKeys(m) {
        var arr = [], k;
        for (k in m) {
            if (m.hasOwnProperty(k)) { arr.push(+k); }
        }
        return arr;
    }

    function arrFill(n, v) {
        var a = new Array(n);
        for (var i = 0; i < n; i++) { a[i] = v; }
        return a;
    }

    function arrCopy(src) {
        var a = new Array(src.length);
        for (var i = 0; i < src.length; i++) { a[i] = src[i]; }
        return a;
    }

    function boolArr(n) { return arrFill(n, false); }

    /**
     * Count heteroatoms (non-carbon, non-hydrogen) in a mapping.
     * Used as tiebreaker: prefer MCS with more heteroatoms when sizes are equal.
     */
    function countHeteroatomsInMapping(g1, mapping) {
        var count = 0;
        var k;
        for (k in mapping) {
            if (mapping.hasOwnProperty(k)) {
                var anum = g1.atomicNum[+k];
                if (anum !== 6 && anum !== 1) { count++; }
            }
        }
        return count;
    }

    /**
     * Compare two mappings: returns true if candidate is better than current best.
     * When sizes are equal, prefer the mapping with more heteroatoms.
     */
    function isBetterMapping(g1, candidate, candidateSize, bestMapping, bestSize, weightMode, candidateScore, bestScore) {
        if (weightMode) { return candidateScore > bestScore; }
        if (candidateSize > bestSize) { return true; }
        if (candidateSize === bestSize && candidateSize > 0) {
            return countHeteroatomsInMapping(g1, candidate) >
                   countHeteroatomsInMapping(g1, bestMapping);
        }
        return false;
    }

    // -----------------------------------------------------------------------
    // TimeBudget
    // -----------------------------------------------------------------------
    function TimeBudget(ms) {
        this.deadline = Date.now() + ms;
    }
    TimeBudget.prototype.expired = function() {
        return Date.now() >= this.deadline;
    };
    TimeBudget.prototype.remainingMs = function() {
        return Math.max(0, this.deadline - Date.now());
    };

    // -----------------------------------------------------------------------
    // Graph accessors (delegate to SMSDGraph)
    // We expect g to have: n, atomicNum[], aromatic[], ring[], degree[],
    //   neighbors[][], label[], orbit[], morganRank[], formalCharge[],
    //   canonicalHash, canonicalLabel[],
    //   bondOrder(i,j), hasBond(i,j)
    // And SMSDGraph to provide: atomsCompatFast(g1,i,g2,j,C),
    //   bondsCompatible(g1,i,j,g2,k,l,C)
    // -----------------------------------------------------------------------
    function atomsCompat(g1, qi, g2, tj, C) {
        return global.SMSDGraph.atomsCompat(g1, qi, g2, tj, C);
    }

    function bondsCompat(g1, qi, qk, g2, tj, tk, C) {
        return global.SMSDGraph.bondsCompat(g1, qi, qk, g2, tj, tk, C);
    }

    // -----------------------------------------------------------------------
    // NLF (Neighbour Label Frequency) helpers
    // -----------------------------------------------------------------------
    function mcsNlfLabel(g, idx) {
        return (g.atomicNum[idx] << 1) | (g.aromatic[idx] ? 1 : 0);
    }

    function buildNLF(g, collector) {
        var nlf = new Array(g.n);
        var i, freq, lab, arr;
        for (i = 0; i < g.n; i++) {
            freq = {};
            collector(g, i, freq);
            arr = [];
            var sorted = [];
            for (lab in freq) {
                if (freq.hasOwnProperty(lab)) { sorted.push(+lab); }
            }
            sorted.sort(function(a, b) { return a - b; });
            for (var s = 0; s < sorted.length; s++) {
                arr.push(sorted[s]);
                arr.push(freq[sorted[s]]);
            }
            nlf[i] = arr;
        }
        return nlf;
    }

    function buildNLF1(g) {
        return buildNLF(g, function(g, idx, freq) {
            var nbs = g.neighbors[idx];
            for (var k = 0; k < nbs.length; k++) {
                var lab = mcsNlfLabel(g, nbs[k]);
                freq[lab] = (freq[lab] || 0) + 1;
            }
        });
    }

    function buildNLF2(g) {
        return buildNLF(g, function(g, idx, freq) {
            var direct = boolArr(g.n);
            direct[idx] = true;
            var nbs = g.neighbors[idx], i, j, nb;
            for (i = 0; i < nbs.length; i++) { direct[nbs[i]] = true; }
            var seen = boolArr(g.n);
            for (i = 0; i < nbs.length; i++) {
                nb = nbs[i];
                var nbs2 = g.neighbors[nb];
                for (j = 0; j < nbs2.length; j++) {
                    var v = nbs2[j];
                    if (direct[v] || seen[v]) { continue; }
                    seen[v] = true;
                    var lab = mcsNlfLabel(g, v);
                    freq[lab] = (freq[lab] || 0) + 1;
                }
            }
        });
    }

    function buildNLF3(g) {
        return buildNLF(g, function(g, idx, freq) {
            var level1 = boolArr(g.n);
            level1[idx] = true;
            var nbs = g.neighbors[idx], i, j, nb;
            for (i = 0; i < nbs.length; i++) { level1[nbs[i]] = true; }
            var level2 = boolArr(g.n);
            for (i = 0; i < nbs.length; i++) {
                nb = nbs[i];
                var nbs2 = g.neighbors[nb];
                for (j = 0; j < nbs2.length; j++) {
                    if (!level1[nbs2[j]]) { level2[nbs2[j]] = true; }
                }
            }
            for (var v = 0; v < g.n; v++) {
                if (!level2[v]) { continue; }
                var nbs3 = g.neighbors[v];
                for (j = 0; j < nbs3.length; j++) {
                    if (!level1[nbs3[j]] && !level2[nbs3[j]]) {
                        var lab = mcsNlfLabel(g, nbs3[j]);
                        freq[lab] = (freq[lab] || 0) + 1;
                    }
                }
            }
        });
    }

    function nlfOk(fq, ft) {
        var fi = 0, ti = 0;
        var fqsz = fq.length, ftsz = ft.length;
        while (fi < fqsz) {
            var fLabel = fq[fi], fFreq = fq[fi + 1];
            while (ti < ftsz && ft[ti] < fLabel) { ti += 2; }
            if (ti >= ftsz || ft[ti] !== fLabel || ft[ti + 1] < fFreq) { return false; }
            fi += 2;
        }
        return true;
    }

    function nlfCheckOk(qi, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                         useTwoHop, useThreeHop) {
        if (!nlfOk(qNLF1[qi], tNLF1[tj])) { return false; }
        if (useTwoHop && qNLF2 && qNLF2.length > 0 && !nlfOk(qNLF2[qi], tNLF2[tj])) { return false; }
        if (useThreeHop && qNLF3 && qNLF3.length > 0 && !nlfOk(qNLF3[qi], tNLF3[tj])) { return false; }
        return true;
    }

    // -----------------------------------------------------------------------
    // L0: Label-frequency upper bound
    // -----------------------------------------------------------------------
    function ubLabel(g, i, C) {
        if (C && !C.matchAtomType) { return 0; }
        var key = g.atomicNum[i] << 2;
        if (C && C.aromaticityMode === 'strict') {
            key |= g.aromatic[i] ? 2 : 0;
        }
        if (C && C.ringMatchesRingOnly) {
            key |= g.ring[i] ? 1 : 0;
        }
        return key;
    }

    function labelFrequencyUpperBound(g1, g2, C) {
        if (g1.n === 0 || g2.n === 0) { return 0; }
        var freq1 = {}, freq2 = {}, i, lab, ub = 0;
        for (i = 0; i < g1.n; i++) {
            lab = ubLabel(g1, i, C);
            freq1[lab] = (freq1[lab] || 0) + 1;
        }
        for (i = 0; i < g2.n; i++) {
            lab = ubLabel(g2, i, C);
            freq2[lab] = (freq2[lab] || 0) + 1;
        }
        for (lab in freq1) {
            if (freq1.hasOwnProperty(lab) && freq2.hasOwnProperty(lab)) {
                ub += Math.min(freq1[lab], freq2[lab]);
            }
        }
        return Math.min(ub, Math.min(g1.n, g2.n));
    }

    // -----------------------------------------------------------------------
    // L0 (alternate): Degree-sequence upper bound
    // Sorts degree sequences of both graphs in descending order, then
    // counts how many positions have min(degQ[i], degT[i]) >= 1.
    // Tighter than label-frequency alone for degree-heterogeneous graphs.
    // -----------------------------------------------------------------------
    function degreeSequenceUpperBound(g1, g2, C) {
        if (g1.n === 0 || g2.n === 0) { return 0; }

        // Build per-label buckets of sorted degrees
        var buckets1 = {}, buckets2 = {}, i, lab;
        for (i = 0; i < g1.n; i++) {
            lab = ubLabel(g1, i, C);
            if (!buckets1[lab]) { buckets1[lab] = []; }
            buckets1[lab].push(g1.degree[i]);
        }
        for (i = 0; i < g2.n; i++) {
            lab = ubLabel(g2, i, C);
            if (!buckets2[lab]) { buckets2[lab] = []; }
            buckets2[lab].push(g2.degree[i]);
        }

        // Sort each bucket descending
        var descSort = function(a, b) { return b - a; };
        for (lab in buckets1) {
            if (buckets1.hasOwnProperty(lab)) { buckets1[lab].sort(descSort); }
        }
        for (lab in buckets2) {
            if (buckets2.hasOwnProperty(lab)) { buckets2[lab].sort(descSort); }
        }

        // For each label, pair up atoms greedily by degree
        var ub = 0;
        for (lab in buckets1) {
            if (!buckets1.hasOwnProperty(lab) || !buckets2.hasOwnProperty(lab)) { continue; }
            var dq = buckets1[lab], dt = buckets2[lab];
            var pairs = Math.min(dq.length, dt.length);
            for (i = 0; i < pairs; i++) {
                // Both must have degree >= 1 (or we allow degree-0 isolated atoms)
                if (dq[i] >= 0 && dt[i] >= 0) { ub++; }
            }
        }
        return Math.min(ub, Math.min(g1.n, g2.n));
    }

    // -----------------------------------------------------------------------
    // Mapped bond counter
    // -----------------------------------------------------------------------
    function countMappedBonds(g1, m) {
        var count = 0, qi, qk, nbs;
        for (qi in m) {
            if (!m.hasOwnProperty(qi)) { continue; }
            qi = +qi;
            nbs = g1.neighbors[qi];
            for (var k = 0; k < nbs.length; k++) {
                qk = nbs[k];
                if (qk > qi && m.hasOwnProperty(qk)) { count++; }
            }
        }
        return count;
    }

    // -----------------------------------------------------------------------
    // MCS scoring
    // -----------------------------------------------------------------------
    function mcsScore(g1, m, M) {
        return M.maximizeBonds ? countMappedBonds(g1, m) : mapSize(m);
    }

    // -----------------------------------------------------------------------
    // Largest connected component
    // -----------------------------------------------------------------------
    function largestConnected(g1, m) {
        if (mapSize(m) === 0) { return m; }
        var mapped = {}, qi;
        for (qi in m) {
            if (m.hasOwnProperty(qi)) { mapped[+qi] = true; }
        }
        var seen = {};
        var comps = [];
        for (qi in m) {
            if (!m.hasOwnProperty(qi)) { continue; }
            qi = +qi;
            if (seen[qi]) { continue; }
            var comp = [];
            var dq = [qi];
            seen[qi] = true;
            while (dq.length > 0) {
                var u = dq.pop();
                comp.push(u);
                var nbs = g1.neighbors[u];
                for (var k = 0; k < nbs.length; k++) {
                    var v = nbs[k];
                    if (mapped[v] && !seen[v] && g1.hasBond(u, v)) {
                        seen[v] = true;
                        dq.push(v);
                    }
                }
            }
            comps.push(comp);
        }
        var bestIdx = 0;
        for (var i = 1; i < comps.length; i++) {
            if (comps[i].length > comps[bestIdx].length) { bestIdx = i; }
        }
        if (comps[bestIdx].length === mapSize(m)) { return m; }
        var result = {};
        for (var c = 0; c < comps[bestIdx].length; c++) {
            var a = comps[bestIdx][c];
            result[a] = m[a];
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Fragment constraints (dMCS)
    // -----------------------------------------------------------------------
    function applyFragmentConstraints(g1, m, minFragSize, maxFrags) {
        var sz = mapSize(m);
        if (sz === 0 || (minFragSize <= 1 && maxFrags >= sz)) { return m; }
        var mapped = {}, qi;
        for (qi in m) {
            if (m.hasOwnProperty(qi)) { mapped[+qi] = true; }
        }
        var seen = {};
        var frags = [];
        for (qi in m) {
            if (!m.hasOwnProperty(qi)) { continue; }
            qi = +qi;
            if (seen[qi]) { continue; }
            var frag = [];
            var dq = [qi];
            seen[qi] = true;
            while (dq.length > 0) {
                var u = dq.pop();
                frag.push(u);
                var nbs = g1.neighbors[u];
                for (var k = 0; k < nbs.length; k++) {
                    var v = nbs[k];
                    if (mapped[v] && !seen[v] && g1.hasBond(u, v)) {
                        seen[v] = true;
                        dq.push(v);
                    }
                }
            }
            frags.push(frag);
        }
        // Remove small fragments
        var filtered = [];
        for (var i = 0; i < frags.length; i++) {
            if (frags[i].length >= minFragSize) { filtered.push(frags[i]); }
        }
        // Sort descending by size
        filtered.sort(function(a, b) { return b.length - a.length; });
        if (filtered.length > maxFrags) { filtered.length = maxFrags; }
        var result = {};
        for (var f = 0; f < filtered.length; f++) {
            for (var j = 0; j < filtered[f].length; j++) {
                var atom = filtered[f][j];
                result[atom] = m[atom];
            }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Ring anchor guard
    // -----------------------------------------------------------------------
    function applyRingAnchorGuard(g1, g2, m, C) {
        if (!(C && C.ringMatchesRingOnly) || mapSize(m) === 0) { return m; }
        var qHasRing = false, tHasRing = false, i;
        for (i = 0; i < g1.n; i++) { if (g1.ring[i]) { qHasRing = true; break; } }
        for (i = 0; i < g2.n; i++) { if (g2.ring[i]) { tHasRing = true; break; } }
        if (qHasRing && !tHasRing) {
            var mappedRing = 0, qi;
            for (qi in m) {
                if (m.hasOwnProperty(qi) && g1.ring[+qi]) { mappedRing++; }
            }
            if (mappedRing === 0) { return {}; }
        }
        return m;
    }

    // -----------------------------------------------------------------------
    // Prune to induced subgraph
    // -----------------------------------------------------------------------
    function pruneToInduced(g1, g2, m, C) {
        if (mapSize(m) === 0) { return m; }
        var M = mapCopy(m);
        var changed;
        do {
            changed = false;
            var keys = mapKeys(M);
            for (var i = 0; i < keys.length && !changed; i++) {
                for (var j = i + 1; j < keys.length && !changed; j++) {
                    var qi = keys[i], qj = keys[j];
                    var ti = M[qi], tj = M[qj];
                    var qOrd = g1.bondOrder(qi, qj), tOrd = g2.bondOrder(ti, tj);
                    var qHas = qOrd !== 0, tHas = tOrd !== 0;
                    var ok = (qHas === tHas);
                    if (ok && qHas) { ok = bondsCompat(g1, qi, qj, g2, ti, tj, C); }
                    if (!ok) {
                        if (g1.degree[qi] >= g1.degree[qj]) { delete M[qi]; }
                        else { delete M[qj]; }
                        changed = true;
                    }
                }
            }
        } while (changed);
        return M;
    }

    // -----------------------------------------------------------------------
    // Post-process MCS (ppx)
    // -----------------------------------------------------------------------
    function ppx(g1, g2, ext, C, M) {
        var changed = true;
        while (changed) {
            var startSize = mapSize(ext);
            if (M.induced) { ext = pruneToInduced(g1, g2, ext, C); }
            if (!M.disconnectedMCS && M.connectedOnly) { ext = largestConnected(g1, ext); }
            changed = mapSize(ext) < startSize;
        }
        ext = applyRingAnchorGuard(g1, g2, ext, C);
        if (M.disconnectedMCS && (M.minFragmentSize > 1 || M.maxFragments < 999999)) {
            ext = applyFragmentConstraints(g1, ext, M.minFragmentSize, M.maxFragments);
        }
        return ext;
    }

    // -----------------------------------------------------------------------
    // Label-frequency pruning check
    // -----------------------------------------------------------------------
    function isPruned(curSize, bestSize, qLF, tLF, freqSize) {
        var potential = curSize;
        for (var lbl = 0; lbl < freqSize; lbl++) {
            if (qLF[lbl] > 0) { potential += Math.min(qLF[lbl], tLF[lbl]); }
        }
        return potential <= bestSize;
    }

    // -----------------------------------------------------------------------
    // L0.75: Greedy probe
    // -----------------------------------------------------------------------
    function greedyProbe(g1, g2, C, fuzzyTol) {
        var n1 = g1.n, n2 = g2.n;
        if (n1 === 0 || n2 === 0) { return {}; }
        if (fuzzyTol === undefined) { fuzzyTol = 0; }

        // Sort query atoms: ring first, then descending degree
        var order = [];
        for (var idx = 0; idx < n1; idx++) { order.push(idx); }
        order.sort(function(a, b) {
            var sa = (g1.ring[a] ? 1000 : 0) + g1.degree[a];
            var sb = (g1.ring[b] ? 1000 : 0) + g1.degree[b];
            return sb - sa;
        });

        var usedT = boolArr(n2);
        var mapping = {};
        var fuzzyUsed = 0;

        for (var oi = 0; oi < n1; oi++) {
            var qi = order[oi];
            var bestTj = -1, bestScore = -1, bestIsFuzzy = false;
            for (var tj = 0; tj < n2; tj++) {
                if (usedT[tj]) { continue; }
                var exactMatch = atomsCompat(g1, qi, g2, tj, C);
                var isFuzzy = false;
                if (!exactMatch) {
                    // Fuzzy atom tolerance: allow mismatch if within budget
                    if (fuzzyTol > 0 && fuzzyUsed < fuzzyTol &&
                        g1.degree[qi] > 0 && g2.degree[tj] > 0) {
                        isFuzzy = true;
                    } else {
                        continue;
                    }
                }
                var ok = true;
                var nbs = g1.neighbors[qi];
                for (var nb = 0; nb < nbs.length; nb++) {
                    var qk = nbs[nb];
                    if (!mapping.hasOwnProperty(qk)) { continue; }
                    var tk = mapping[qk];
                    var qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                    if ((qOrd === 0) !== (tOrd === 0)) { ok = false; break; }
                    if (qOrd !== 0 && !bondsCompat(g1, qi, qk, g2, tj, tk, C)) { ok = false; break; }
                }
                if (!ok) { continue; }
                var score = (isFuzzy ? -500 : 0)
                          + (g2.ring[tj] === g1.ring[qi] ? 100 : 0)
                          + (g2.degree[tj] === g1.degree[qi] ? 50 : 0)
                          + Math.min(g1.degree[qi], g2.degree[tj]);
                if (score > bestScore) { bestScore = score; bestTj = tj; bestIsFuzzy = isFuzzy; }
            }
            if (bestTj >= 0) {
                mapping[qi] = bestTj; usedT[bestTj] = true;
                if (bestIsFuzzy) { fuzzyUsed++; }
            }
        }
        return mapping;
    }

    // -----------------------------------------------------------------------
    // Greedy atom extend
    // -----------------------------------------------------------------------
    function greedyAtomExtend(g1, g2, seed, C, M) {
        var n1 = g1.n, n2 = g2.n;
        var q2t = arrFill(n1, -1), t2q = arrFill(n2, -1);
        var k;
        for (k in seed) {
            if (seed.hasOwnProperty(k)) { q2t[+k] = seed[k]; t2q[seed[k]] = +k; }
        }
        var progress = true;
        while (progress) {
            progress = false;
            for (var qi = 0; qi < n1; qi++) {
                if (q2t[qi] >= 0) { continue; }
                var onFrontier = false;
                var nbs = g1.neighbors[qi];
                for (var nb = 0; nb < nbs.length; nb++) {
                    if (q2t[nbs[nb]] >= 0) { onFrontier = true; break; }
                }
                if (!onFrontier) { continue; }
                var bestTj = -1, bestScore = -1;
                for (var tj = 0; tj < n2; tj++) {
                    if (t2q[tj] >= 0) { continue; }
                    if (!atomsCompat(g1, qi, g2, tj, C)) { continue; }
                    var consistent = true;
                    for (var nb2 = 0; nb2 < nbs.length; nb2++) {
                        var qk = nbs[nb2];
                        if (q2t[qk] < 0) { continue; }
                        var tk = q2t[qk];
                        var qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                        if (qOrd !== 0 && tOrd !== 0) {
                            if (!bondsCompat(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
                        } else if (M.induced && ((qOrd !== 0) !== (tOrd !== 0))) {
                            consistent = false; break;
                        }
                    }
                    if (!consistent) { continue; }
                    var score = (g2.ring[tj] && g1.ring[qi] ? 50 : 0)
                              + Math.min(g1.degree[qi], g2.degree[tj]);
                    if (score > bestScore) { bestScore = score; bestTj = tj; }
                }
                if (bestTj >= 0) {
                    q2t[qi] = bestTj; t2q[bestTj] = qi;
                    progress = true;
                }
            }
        }
        var result = {};
        for (var i = 0; i < n1; i++) {
            if (q2t[i] >= 0) { result[i] = q2t[i]; }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Build frontier for McGregor DFS
    // -----------------------------------------------------------------------
    function buildFrontier(g1, cur, usedQ, inFrontier, frontierBuf, connectedOnly) {
        var count = 0, qi, nbs, k, qn;
        if (connectedOnly === undefined) { connectedOnly = true; }
        for (qi in cur) {
            if (!cur.hasOwnProperty(qi)) { continue; }
            qi = +qi;
            nbs = g1.neighbors[qi];
            for (k = 0; k < nbs.length; k++) {
                qn = nbs[k];
                if (!usedQ[qn] && !inFrontier[qn] && !cur.hasOwnProperty(qn)) {
                    inFrontier[qn] = true;
                    frontierBuf[count++] = qn;
                }
            }
        }
        if (count === 0 && !connectedOnly) {
            for (var i = 0; i < g1.n; i++) {
                if (!usedQ[i] && !cur.hasOwnProperty(i)) {
                    frontierBuf[count++] = i;
                }
            }
        }
        for (var f = 0; f < count; f++) { inFrontier[frontierBuf[f]] = false; }
        return count;
    }

    // -----------------------------------------------------------------------
    // Find best candidate for McGregor
    // -----------------------------------------------------------------------
    function findBestCandidate(g1, g2, C, cur, usedT,
                               qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                               useTwoHopNLF, useThreeHopNLF,
                               frontierCount, frontierBuf,
                               candBuf, bestCandBuf) {
        var bestQi = -1, bestCandSize = 999999999, bestCandCount = 0;
        for (var fi = 0; fi < frontierCount; fi++) {
            var qi = frontierBuf[fi];
            var candCount = 0;
            for (var tj = 0; tj < g2.n; tj++) {
                if (usedT[tj]) { continue; }
                if (!atomsCompat(g1, qi, g2, tj, C)) { continue; }
                if (!nlfCheckOk(qi, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                useTwoHopNLF, useThreeHopNLF)) { continue; }
                var ok = true;
                var nbs = g1.neighbors[qi];
                for (var nb = 0; nb < nbs.length; nb++) {
                    var qk = nbs[nb];
                    if (!cur.hasOwnProperty(qk)) { continue; }
                    var tl = cur[qk];
                    var qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tl);
                    if ((qOrd === 0) !== (tOrd === 0)) { ok = false; break; }
                    if (qOrd !== 0 && !bondsCompat(g1, qi, qk, g2, tj, tl, C)) { ok = false; break; }
                }
                if (ok) { candBuf[candCount++] = tj; }
            }
            if (candCount === 0) { continue; }
            if (candCount < bestCandSize) {
                bestCandSize = candCount;
                bestQi = qi;
                bestCandCount = candCount;
                for (var cc = 0; cc < candCount; cc++) { bestCandBuf[cc] = candBuf[cc]; }
            }
        }
        return { qi: bestQi, count: bestCandCount };
    }

    // -----------------------------------------------------------------------
    // McGregor DFS (Level 4 — atom-frontier)
    // -----------------------------------------------------------------------
    function mcGregorDFS(g1, g2, C, cur, bestRef,
                         qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                         useTwoHopNLF, useThreeHopNLF,
                         tb, localDeadline, depth,
                         usedQ, usedT,
                         candBuf, bestCandBuf,
                         qLabelFreq, tLabelFreq, freqSize,
                         inFrontier, frontierBuf, jointQ, jointT) {

        if (Date.now() >= localDeadline || tb.expired()) { return; }
        var curSz = mapSize(cur);
        var bestSz = mapSize(bestRef.best);
        if (curSz > bestSz) { bestRef.best = mapCopy(cur); }
        if (isPruned(curSz, mapSize(bestRef.best), qLabelFreq, tLabelFreq, freqSize)) { return; }

        var frontierCount = buildFrontier(g1, cur, usedQ, inFrontier, frontierBuf);
        if (frontierCount === 0) { return; }

        var res = findBestCandidate(g1, g2, C, cur, usedT,
                                    qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                    useTwoHopNLF, useThreeHopNLF,
                                    frontierCount, frontierBuf,
                                    candBuf, bestCandBuf);
        var bestQi = res.qi, bestCandCount = res.count;
        if (bestQi === -1) { return; }

        // Unit propagation (forced assignment)
        var forcedQ = [], forcedT = [];
        while (bestCandCount === 1 && !(Date.now() >= localDeadline || tb.expired())) {
            var fq = bestQi, ft = bestCandBuf[0];
            cur[fq] = ft; usedQ[fq] = true; usedT[ft] = true;
            qLabelFreq[jointQ[fq]]--; tLabelFreq[jointT[ft]]--;
            forcedQ.push(fq); forcedT.push(ft);
            depth++;
            curSz = mapSize(cur);
            if (curSz > mapSize(bestRef.best)) { bestRef.best = mapCopy(cur); }
            if (isPruned(curSz, mapSize(bestRef.best), qLabelFreq, tLabelFreq, freqSize)) { break; }

            frontierCount = buildFrontier(g1, cur, usedQ, inFrontier, frontierBuf);
            if (frontierCount === 0) { break; }
            res = findBestCandidate(g1, g2, C, cur, usedT,
                                    qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                    useTwoHopNLF, useThreeHopNLF,
                                    frontierCount, frontierBuf,
                                    candBuf, bestCandBuf);
            bestQi = res.qi; bestCandCount = res.count;
            if (bestQi === -1) { break; }
        }

        if (bestQi !== -1 && bestCandCount > 1) {
            var branchLimit = depth < 5 ? bestCandCount : Math.min(bestCandCount, 16);
            var savedCands = arrCopy(bestCandBuf);
            for (var i = 0; i < branchLimit; i++) {
                if (Date.now() >= localDeadline || tb.expired()) { break; }
                var bestTj = savedCands[i];
                cur[bestQi] = bestTj; usedQ[bestQi] = true; usedT[bestTj] = true;
                qLabelFreq[jointQ[bestQi]]--; tLabelFreq[jointT[bestTj]]--;
                mcGregorDFS(g1, g2, C, cur, bestRef,
                            qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                            useTwoHopNLF, useThreeHopNLF,
                            tb, localDeadline, depth + 1,
                            usedQ, usedT, candBuf, bestCandBuf,
                            qLabelFreq, tLabelFreq, freqSize,
                            inFrontier, frontierBuf, jointQ, jointT);
                qLabelFreq[jointQ[bestQi]]++; tLabelFreq[jointT[bestTj]]++;
                delete cur[bestQi]; usedQ[bestQi] = false; usedT[bestTj] = false;
            }
        }

        // Undo forced assignments
        for (var f = forcedQ.length - 1; f >= 0; f--) {
            var ufq = forcedQ[f], uft = forcedT[f];
            qLabelFreq[jointQ[ufq]]++; tLabelFreq[jointT[uft]]++;
            delete cur[ufq]; usedQ[ufq] = false; usedT[uft] = false;
        }
    }

    // -----------------------------------------------------------------------
    // McGregor bond-grow (Level 4 — bond-oriented)
    // -----------------------------------------------------------------------
    function mcGregorBondGrow(g1, g2, C, cur, bestRef,
                              qNLF1, tNLF1, useTwoHopNLF, useThreeHopNLF,
                              qNLF2, tNLF2, qNLF3, tNLF3,
                              tb, localDeadline, depth,
                              usedQ, usedT,
                              qLabelFreq, tLabelFreq, freqSize,
                              q2tMap, inFrontier, frontierBuf,
                              candBuf, bestCandBuf, jointQ, jointT) {

        if (Date.now() >= localDeadline || tb.expired()) { return; }
        var curSz = mapSize(cur);
        if (curSz > mapSize(bestRef.best)) { bestRef.best = mapCopy(cur); }
        if (isPruned(curSz, mapSize(bestRef.best), qLabelFreq, tLabelFreq, freqSize)) { return; }

        // Find best frontier bond
        var bestQk = -1, bestQi_local = -1, bestCandSize = 999999999, bestCandCount = 0;
        var qi, nbs, k, qk, candCount, tk, tj, tOrd, qOrd;

        for (qi in cur) {
            if (!cur.hasOwnProperty(qi)) { continue; }
            qi = +qi;
            var mappedTi = q2tMap[qi];
            nbs = g1.neighbors[qi];
            for (k = 0; k < nbs.length; k++) {
                qk = nbs[k];
                if (usedQ[qk] || inFrontier[qk]) { continue; }
                candCount = 0;
                var tNbs = g2.neighbors[mappedTi];
                for (var t = 0; t < tNbs.length; t++) {
                    tk = tNbs[t];
                    if (usedT[tk]) { continue; }
                    qOrd = g1.bondOrder(qi, qk); tOrd = g2.bondOrder(mappedTi, tk);
                    if ((qOrd === 0) !== (tOrd === 0)) { continue; }
                    if (qOrd !== 0 && !bondsCompat(g1, qi, qk, g2, mappedTi, tk, C)) { continue; }
                    if (!atomsCompat(g1, qk, g2, tk, C)) { continue; }
                    if (!nlfCheckOk(qk, tk, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                    useTwoHopNLF, useThreeHopNLF)) { continue; }
                    var ok = true;
                    var qkNbs = g1.neighbors[qk];
                    for (var nn = 0; nn < qkNbs.length; nn++) {
                        var qn = qkNbs[nn];
                        if (qn === qi || !usedQ[qn]) { continue; }
                        var tn = q2tMap[qn];
                        var qOrd2 = g1.bondOrder(qk, qn), tOrd2 = g2.bondOrder(tk, tn);
                        if ((qOrd2 === 0) !== (tOrd2 === 0)) { ok = false; break; }
                        if (qOrd2 !== 0 && !bondsCompat(g1, qk, qn, g2, tk, tn, C)) { ok = false; break; }
                    }
                    if (ok) { candBuf[candCount++] = tk; }
                }
                if (candCount > 0 && candCount < bestCandSize) {
                    bestCandSize = candCount; bestQk = qk; bestQi_local = qi;
                    bestCandCount = candCount;
                    for (var cc = 0; cc < candCount; cc++) { bestCandBuf[cc] = candBuf[cc]; }
                }
            }
        }

        // Fallback: disconnected extension
        if (bestQk === -1) {
            for (var i = 0; i < g1.n; i++) {
                if (usedQ[i]) { continue; }
                candCount = 0;
                for (tj = 0; tj < g2.n; tj++) {
                    if (usedT[tj]) { continue; }
                    if (!atomsCompat(g1, i, g2, tj, C)) { continue; }
                    if (!nlfCheckOk(i, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                    useTwoHopNLF, useThreeHopNLF)) { continue; }
                    candBuf[candCount++] = tj;
                }
                if (candCount > 0 && candCount < bestCandSize) {
                    bestCandSize = candCount; bestQk = i; bestQi_local = -1;
                    bestCandCount = candCount;
                    for (var cc2 = 0; cc2 < candCount; cc2++) { bestCandBuf[cc2] = candBuf[cc2]; }
                }
            }
        }
        if (bestQk === -1) { return; }

        // Unit propagation
        var forcedQ = [], forcedT = [];
        while (bestCandCount === 1 && !(Date.now() >= localDeadline || tb.expired())) {
            var fq = bestQk, ft = bestCandBuf[0];
            cur[fq] = ft; usedQ[fq] = true; usedT[ft] = true; q2tMap[fq] = ft;
            qLabelFreq[jointQ[fq]]--; tLabelFreq[jointT[ft]]--;
            forcedQ.push(fq); forcedT.push(ft);
            depth++;
            if (mapSize(cur) > mapSize(bestRef.best)) { bestRef.best = mapCopy(cur); }
            if (isPruned(mapSize(cur), mapSize(bestRef.best), qLabelFreq, tLabelFreq, freqSize)) { break; }

            bestQk = -1; bestCandSize = 999999999; bestCandCount = 0;
            for (qi in cur) {
                if (!cur.hasOwnProperty(qi)) { continue; }
                qi = +qi;
                var mTi = q2tMap[qi];
                nbs = g1.neighbors[qi];
                for (k = 0; k < nbs.length; k++) {
                    qk = nbs[k];
                    if (usedQ[qk]) { continue; }
                    candCount = 0;
                    var tNbs2 = g2.neighbors[mTi];
                    for (var t2 = 0; t2 < tNbs2.length; t2++) {
                        tk = tNbs2[t2];
                        if (usedT[tk]) { continue; }
                        qOrd = g1.bondOrder(qi, qk); tOrd = g2.bondOrder(mTi, tk);
                        if ((qOrd === 0) !== (tOrd === 0)) { continue; }
                        if (qOrd !== 0 && !bondsCompat(g1, qi, qk, g2, mTi, tk, C)) { continue; }
                        if (!atomsCompat(g1, qk, g2, tk, C)) { continue; }
                        if (!nlfCheckOk(qk, tk, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                        useTwoHopNLF, useThreeHopNLF)) { continue; }
                        var ok2 = true;
                        var qkNbs2 = g1.neighbors[qk];
                        for (var nn2 = 0; nn2 < qkNbs2.length; nn2++) {
                            var qn2 = qkNbs2[nn2];
                            if (qn2 === qi || !usedQ[qn2]) { continue; }
                            var tn2 = q2tMap[qn2];
                            var qOrd3 = g1.bondOrder(qk, qn2), tOrd3 = g2.bondOrder(tk, tn2);
                            if ((qOrd3 === 0) !== (tOrd3 === 0)) { ok2 = false; break; }
                            if (qOrd3 !== 0 && !bondsCompat(g1, qk, qn2, g2, tk, tn2, C)) { ok2 = false; break; }
                        }
                        if (ok2) { candBuf[candCount++] = tk; }
                    }
                    if (candCount > 0 && candCount < bestCandSize) {
                        bestCandSize = candCount; bestQk = qk; bestQi_local = qi;
                        bestCandCount = candCount;
                        for (var cc3 = 0; cc3 < candCount; cc3++) { bestCandBuf[cc3] = candBuf[cc3]; }
                    }
                }
            }
            if (bestQk === -1) { break; }
        }

        if (bestQk !== -1 && bestCandCount > 1) {
            var branchLimit = depth < 5 ? bestCandCount : Math.min(bestCandCount, 16);
            var savedCands = [];
            for (var sc = 0; sc < bestCandCount; sc++) { savedCands.push(bestCandBuf[sc]); }
            for (var bi = 0; bi < branchLimit; bi++) {
                if (Date.now() >= localDeadline || tb.expired()) { break; }
                var btj = savedCands[bi];
                cur[bestQk] = btj; usedQ[bestQk] = true; usedT[btj] = true; q2tMap[bestQk] = btj;
                qLabelFreq[jointQ[bestQk]]--; tLabelFreq[jointT[btj]]--;
                mcGregorBondGrow(g1, g2, C, cur, bestRef,
                                 qNLF1, tNLF1, useTwoHopNLF, useThreeHopNLF,
                                 qNLF2, tNLF2, qNLF3, tNLF3,
                                 tb, localDeadline, depth + 1,
                                 usedQ, usedT, qLabelFreq, tLabelFreq, freqSize,
                                 q2tMap, inFrontier, frontierBuf,
                                 candBuf, bestCandBuf, jointQ, jointT);
                qLabelFreq[jointQ[bestQk]]++; tLabelFreq[jointT[btj]]++;
                delete cur[bestQk]; usedQ[bestQk] = false; usedT[btj] = false; q2tMap[bestQk] = -1;
            }
        }

        // Undo forced
        for (var uf = forcedQ.length - 1; uf >= 0; uf--) {
            var ufq = forcedQ[uf], uft = forcedT[uf];
            qLabelFreq[jointQ[ufq]]++; tLabelFreq[jointT[uft]]++;
            delete cur[ufq]; usedQ[ufq] = false; usedT[uft] = false; q2tMap[ufq] = -1;
        }
    }

    // -----------------------------------------------------------------------
    // McGregor extend entry point (Level 4)
    // -----------------------------------------------------------------------
    function mcGregorExtend(g1, g2, seed, C, tb, localMillis,
                            useTwoHopNLF, useThreeHopNLF, connectedOnly) {
        var localDeadline = Date.now() + localMillis;
        var bestRef = { best: mapCopy(seed) };

        var qNLF1 = buildNLF1(g1), tNLF1 = buildNLF1(g2);
        var emptyNLF = [];
        var qNLF2 = useTwoHopNLF ? buildNLF2(g1) : emptyNLF;
        var tNLF2 = useTwoHopNLF ? buildNLF2(g2) : emptyNLF;
        var qNLF3 = useThreeHopNLF ? buildNLF3(g1) : emptyNLF;
        var tNLF3 = useThreeHopNLF ? buildNLF3(g2) : emptyNLF;

        var n1 = g1.n, n2 = g2.n;
        var usedQ = boolArr(n1), usedT = boolArr(n2);
        var candBuf = new Array(n2), bestCandBuf = new Array(n2);
        var k;
        for (k in seed) {
            if (seed.hasOwnProperty(k)) { usedQ[+k] = true; usedT[seed[k]] = true; }
        }

        var maxLabel = 0, i;
        for (i = 0; i < n1; i++) { if (g1.label[i] > maxLabel) { maxLabel = g1.label[i]; } }
        for (i = 0; i < n2; i++) { if (g2.label[i] > maxLabel) { maxLabel = g2.label[i]; } }
        var freqSize = maxLabel + 1;
        var qLabelFreq = arrFill(freqSize, 0), tLabelFreq = arrFill(freqSize, 0);
        var jointQ = new Array(n1), jointT = new Array(n2);
        for (i = 0; i < n1; i++) { jointQ[i] = g1.label[i]; if (!usedQ[i]) { qLabelFreq[jointQ[i]]++; } }
        for (i = 0; i < n2; i++) { jointT[i] = g2.label[i]; if (!usedT[i]) { tLabelFreq[jointT[i]]++; } }

        var inFrontier = boolArr(n1);
        var frontierBuf = new Array(n1);

        // Fast bail-out
        if (mapSize(seed) > 0) {
            var hasExtensible = false;
            for (k in seed) {
                if (!seed.hasOwnProperty(k)) { continue; }
                var nbs = g1.neighbors[+k];
                for (var nb = 0; nb < nbs.length; nb++) {
                    if (usedQ[nbs[nb]]) { continue; }
                    for (var tj = 0; tj < n2; tj++) {
                        if (usedT[tj]) { continue; }
                        if (atomsCompat(g1, nbs[nb], g2, tj, C)) {
                            hasExtensible = true;
                            break;
                        }
                    }
                    if (hasExtensible) { break; }
                }
                if (hasExtensible) { break; }
            }
            if (!hasExtensible && connectedOnly) { return bestRef.best; }
        }

        // Bond-growth first pass for larger molecules
        if (mapSize(seed) > 0 && n2 >= 20) {
            var q2tMap = arrFill(n1, -1);
            for (k in seed) {
                if (seed.hasOwnProperty(k)) { q2tMap[+k] = seed[k]; }
            }
            var bondBestRef = { best: mapCopy(seed) };
            var bondDeadline = Date.now() + Math.max(1, Math.min(localMillis / 10, 20));

            var usedQCopy = arrCopy(usedQ);
            var usedTCopy = arrCopy(usedT);
            var qLFCopy = arrCopy(qLabelFreq);
            var tLFCopy = arrCopy(tLabelFreq);
            var inFCopy = arrCopy(inFrontier);
            var curCopy = mapCopy(seed);

            mcGregorBondGrow(g1, g2, C, curCopy, bondBestRef,
                             qNLF1, tNLF1, useTwoHopNLF, useThreeHopNLF,
                             qNLF2, tNLF2, qNLF3, tNLF3,
                             tb, bondDeadline, 0,
                             usedQCopy, usedTCopy, qLFCopy, tLFCopy, freqSize,
                             q2tMap, inFCopy, frontierBuf,
                             candBuf, bestCandBuf, jointQ, jointT);
            for (i = 0; i < n1; i++) { inFrontier[i] = false; }
            if (mapSize(bondBestRef.best) > mapSize(bestRef.best)) {
                bestRef.best = bondBestRef.best;
            }
        }

        // Atom-frontier DFS
        var curDFS = mapCopy(seed);
        mcGregorDFS(g1, g2, C, curDFS, bestRef,
                    qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                    useTwoHopNLF, useThreeHopNLF,
                    tb, localDeadline, 0,
                    usedQ, usedT, candBuf, bestCandBuf,
                    qLabelFreq, tLabelFreq, freqSize,
                    inFrontier, frontierBuf, jointQ, jointT);

        return bestRef.best;
    }

    // =======================================================================
    // L1.5: Seed-and-extend
    // =======================================================================
    function greedyBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced) {
        var mapSz = 0, i;
        for (i = 0; i < n1; i++) { if (q2t[i] >= 0) { mapSz++; } }
        var progress = true;
        while (progress) {
            progress = false;
            for (var qi = 0; qi < n1; qi++) {
                if (q2t[qi] >= 0) { continue; }
                var bestTj = -1, bestScore = -1;
                var nbs = g1.neighbors[qi];
                for (var nb = 0; nb < nbs.length; nb++) {
                    if (q2t[nbs[nb]] < 0) { continue; }
                    var tNb = q2t[nbs[nb]];
                    var tNbs = g2.neighbors[tNb];
                    for (var t = 0; t < tNbs.length; t++) {
                        var tj = tNbs[t];
                        if (t2q[tj] >= 0) { continue; }
                        if (!atomsCompat(g1, qi, g2, tj, C)) { continue; }
                        if (!bondsCompat(g1, qi, nbs[nb], g2, tj, tNb, C)) { continue; }
                        var consistent = true;
                        for (var nk = 0; nk < nbs.length; nk++) {
                            var qk = nbs[nk];
                            if (qk === nbs[nb] || q2t[qk] < 0) { continue; }
                            var tk = q2t[qk];
                            var qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                            if (qOrd !== 0 && tOrd !== 0) {
                                if (!bondsCompat(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
                            } else if (induced && ((qOrd !== 0) !== (tOrd !== 0))) {
                                consistent = false; break;
                            }
                        }
                        if (!consistent) { continue; }
                        var score = (g2.ring[tj] && g1.ring[qi] ? 50 : 0)
                                  + Math.min(g1.degree[qi], g2.degree[tj]);
                        if (score > bestScore) { bestScore = score; bestTj = tj; }
                    }
                }
                if (bestTj >= 0) {
                    q2t[qi] = bestTj; t2q[bestTj] = qi; mapSz++; progress = true;
                }
            }
        }
        return mapSz;
    }

    function backtrackBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced,
                                  bestQ2T, bestSizeRef, tb, nodeCountRef, nodeLimit) {
        var curSize = 0, i;
        for (i = 0; i < n1; i++) { if (q2t[i] >= 0) { curSize++; } }
        if (curSize > bestSizeRef.val) {
            bestSizeRef.val = curSize;
            for (i = 0; i < n1; i++) { bestQ2T[i] = q2t[i]; }
        }
        if (nodeCountRef.val > nodeLimit || tb.expired()) { return; }

        var frontier = 0;
        for (var qi = 0; qi < n1; qi++) {
            if (q2t[qi] >= 0) { continue; }
            var nbs = g1.neighbors[qi];
            for (var nb = 0; nb < nbs.length; nb++) {
                if (q2t[nbs[nb]] >= 0) { frontier++; break; }
            }
        }
        if (curSize + frontier <= bestSizeRef.val) { return; }

        var bestQi = -1, bestConstraint = -1;
        for (qi = 0; qi < n1; qi++) {
            if (q2t[qi] >= 0) { continue; }
            var mappedNb = 0;
            var qnbs = g1.neighbors[qi];
            for (var mnb = 0; mnb < qnbs.length; mnb++) {
                if (q2t[qnbs[mnb]] >= 0) { mappedNb++; }
            }
            if (mappedNb === 0) { continue; }
            var constraint = mappedNb * 100 + (g1.ring[qi] ? 50 : 0) + g1.degree[qi];
            if (constraint > bestConstraint) { bestConstraint = constraint; bestQi = qi; }
        }
        if (bestQi < 0) { return; }

        qi = bestQi;
        var qiNbs = g1.neighbors[qi];
        for (var ni = 0; ni < qiNbs.length; ni++) {
            if (q2t[qiNbs[ni]] < 0) { continue; }
            var tNbs = g2.neighbors[q2t[qiNbs[ni]]];
            for (var ti = 0; ti < tNbs.length; ti++) {
                var tj = tNbs[ti];
                if (t2q[tj] >= 0) { continue; }
                nodeCountRef.val++;
                if (nodeCountRef.val > nodeLimit || tb.expired()) { return; }
                if (!atomsCompat(g1, qi, g2, tj, C)) { continue; }
                var consistent = true;
                for (var ck = 0; ck < qiNbs.length; ck++) {
                    var qk = qiNbs[ck];
                    if (q2t[qk] < 0) { continue; }
                    var tk = q2t[qk];
                    var qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                    if (qOrd !== 0 && tOrd !== 0) {
                        if (!bondsCompat(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
                    } else if (induced && ((qOrd !== 0) !== (tOrd !== 0))) {
                        consistent = false; break;
                    }
                }
                if (!consistent) { continue; }
                q2t[qi] = tj; t2q[tj] = qi;
                backtrackBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced,
                                    bestQ2T, bestSizeRef, tb, nodeCountRef, nodeLimit);
                q2t[qi] = -1; t2q[tj] = -1;
                if (nodeCountRef.val > nodeLimit || tb.expired()) { return; }
            }
        }
    }

    function seedExtendMCS(g1, g2, C, tb, upperBound, induced) {
        var n1 = g1.n, n2 = g2.n;
        if (n1 === 0 || n2 === 0) { return {}; }

        var maxLabel = 0, i, j;
        for (i = 0; i < n1; i++) { if (g1.label[i] > maxLabel) { maxLabel = g1.label[i]; } }
        for (j = 0; j < n2; j++) { if (g2.label[j] > maxLabel) { maxLabel = g2.label[j]; } }
        var labelFreqG1 = arrFill(maxLabel + 1, 0), labelFreqG2 = arrFill(maxLabel + 1, 0);
        for (i = 0; i < n1; i++) { labelFreqG1[g1.label[i]]++; }
        for (j = 0; j < n2; j++) { labelFreqG2[g2.label[j]]++; }

        // Compatible targets per query atom
        var compatT = new Array(n1);
        for (i = 0; i < n1; i++) {
            compatT[i] = [];
            for (j = 0; j < n2; j++) {
                if (atomsCompat(g1, i, g2, j, C)) { compatT[i].push(j); }
            }
        }

        // Collect bonds
        var g1Bonds = [];
        for (i = 0; i < n1; i++) {
            var nbs = g1.neighbors[i];
            for (var nb = 0; nb < nbs.length; nb++) {
                if (nbs[nb] > i) { g1Bonds.push({ u: i, v: nbs[nb] }); }
            }
        }
        var g1BondCount = g1Bonds.length;

        // Score and sort seeds
        var seedOrder = new Array(g1BondCount);
        for (var b = 0; b < g1BondCount; b++) {
            var u = g1Bonds[b].u, v = g1Bonds[b].v;
            var totalAtoms = n1 + n2;
            var uFreq = labelFreqG1[g1.label[u]]
                + (g1.label[u] <= maxLabel ? (labelFreqG2[g1.label[u]] || 0) : 0);
            var vFreq = labelFreqG1[g1.label[v]]
                + (g1.label[v] <= maxLabel ? (labelFreqG2[g1.label[v]] || 0) : 0);
            var score = Math.floor((totalAtoms * 2) / Math.max(1, uFreq))
                      + Math.floor((totalAtoms * 2) / Math.max(1, vFreq))
                      + Math.floor((n2 * 2) / Math.max(1, compatT[u].length))
                      + Math.floor((n2 * 2) / Math.max(1, compatT[v].length))
                      + g1.degree[u] + g1.degree[v]
                      + (g1.ring[u] ? 10 : 0) + (g1.ring[v] ? 10 : 0);
            seedOrder[b] = { idx: b, score: score };
        }
        seedOrder.sort(function(a, b2) { return b2.score - a.score; });

        var bestQ2T = arrFill(n1, -1);
        var bestSize = 0;
        var nodeCount = { val: 0 };
        var minN = Math.min(n1, n2);
        var maxSeeds = Math.min(g1BondCount, 10);
        var atomUsedAsSeed = boolArr(n1);
        var triedOrbitPairs = {};
        var seedsTried = 0;

        for (var si = 0; si < g1BondCount && seedsTried < maxSeeds; si++) {
            if (tb.expired() || nodeCount.val > MAX_NODE_LIMIT) { break; }
            var bondIdx = seedOrder[si].idx;
            var qu = g1Bonds[bondIdx].u, qv = g1Bonds[bondIdx].v;
            var oA = Math.min(g1.orbit[qu], g1.orbit[qv]);
            var oB = Math.max(g1.orbit[qu], g1.orbit[qv]);
            var pairKey = oA + '_' + oB;
            if (triedOrbitPairs[pairKey]) { continue; }
            triedOrbitPairs[pairKey] = true;
            if (atomUsedAsSeed[qu] && atomUsedAsSeed[qv]) { continue; }
            atomUsedAsSeed[qu] = true; atomUsedAsSeed[qv] = true; seedsTried++;

            // Forward: qu->ta, qv->tb2
            for (var ca = 0; ca < compatT[qu].length; ca++) {
                var ta = compatT[qu][ca];
                if (tb.expired() || nodeCount.val > MAX_NODE_LIMIT) { break; }
                var taNbs = g2.neighbors[ta];
                for (var tb2i = 0; tb2i < taNbs.length; tb2i++) {
                    var tb2 = taNbs[tb2i];
                    nodeCount.val++;
                    if (nodeCount.val > MAX_NODE_LIMIT || tb.expired()) { break; }
                    if (!atomsCompat(g1, qv, g2, tb2, C)) { continue; }
                    if (!bondsCompat(g1, qu, qv, g2, ta, tb2, C)) { continue; }
                    var q2t = arrFill(n1, -1), t2q = arrFill(n2, -1);
                    q2t[qu] = ta; t2q[ta] = qu; q2t[qv] = tb2; t2q[tb2] = qv;
                    var ms = greedyBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced);
                    nodeCount.val += ms;
                    if (ms > bestSize) {
                        bestSize = ms;
                        for (var ci = 0; ci < n1; ci++) { bestQ2T[ci] = q2t[ci]; }
                    }
                }
            }
            // Reverse: qv->ta, qu->tb2
            for (var ca2 = 0; ca2 < compatT[qv].length; ca2++) {
                var ta2 = compatT[qv][ca2];
                if (tb.expired() || nodeCount.val > MAX_NODE_LIMIT) { break; }
                var ta2Nbs = g2.neighbors[ta2];
                for (var tb3i = 0; tb3i < ta2Nbs.length; tb3i++) {
                    var tb3 = ta2Nbs[tb3i];
                    nodeCount.val++;
                    if (nodeCount.val > MAX_NODE_LIMIT || tb.expired()) { break; }
                    if (!atomsCompat(g1, qu, g2, tb3, C)) { continue; }
                    if (!bondsCompat(g1, qu, qv, g2, tb3, ta2, C)) { continue; }
                    var q2t2 = arrFill(n1, -1), t2q2 = arrFill(n2, -1);
                    q2t2[qu] = tb3; t2q2[tb3] = qu; q2t2[qv] = ta2; t2q2[ta2] = qv;
                    var ms2 = greedyBondExtend(g1, g2, C, q2t2, t2q2, n1, n2, induced);
                    nodeCount.val += ms2;
                    if (ms2 > bestSize) {
                        bestSize = ms2;
                        for (var ci2 = 0; ci2 < n1; ci2++) { bestQ2T[ci2] = q2t2[ci2]; }
                    }
                }
            }
            if (bestSize >= upperBound || bestSize >= minN) { break; }
            if (bestSize >= (upperBound * 19 / 20) && nodeCount.val > SEED_EXTEND_NODE_LIMIT / 2) { break; }
        }

        // Backtrack refinement
        if (bestSize > 0 && bestSize < minN && !tb.expired() && nodeCount.val < MAX_NODE_LIMIT) {
            var q2tBT = arrFill(n1, -1), t2qBT = arrFill(n2, -1);
            for (i = 0; i < n1; i++) {
                if (bestQ2T[i] >= 0) { q2tBT[i] = bestQ2T[i]; t2qBT[bestQ2T[i]] = i; }
            }
            var bestSizeRef = { val: bestSize };
            backtrackBondExtend(g1, g2, C, q2tBT, t2qBT, n1, n2, induced,
                                bestQ2T, bestSizeRef, tb, nodeCount, MAX_NODE_LIMIT);
            bestSize = 0;
            for (i = 0; i < n1; i++) { if (bestQ2T[i] >= 0) { bestSize++; } }
        }

        var result = {};
        for (i = 0; i < n1; i++) { if (bestQ2T[i] >= 0) { result[i] = bestQ2T[i]; } }
        return result;
    }

    // =======================================================================
    // L2: McSplit partition refinement
    // =======================================================================

    // BitClass — bitset for McSplit partition classes
    function BitClass(n) {
        this.n_bits = n;
        this.wordCount = ((n + 31) >>> 5); // 32-bit words for JS (no 64-bit int)
        this.words = arrFill(this.wordCount, 0);
        this.card = 0;
    }
    BitClass.prototype.set = function(i) {
        var w = i >>> 5, bit = 1 << (i & 31);
        if (!(this.words[w] & bit)) { this.words[w] |= bit; this.card++; }
    };
    BitClass.prototype.clear = function(i) {
        var w = i >>> 5, bit = 1 << (i & 31);
        if (this.words[w] & bit) { this.words[w] &= ~bit; this.card--; }
    };
    BitClass.prototype.get = function(i) {
        return (this.words[i >>> 5] & (1 << (i & 31))) !== 0;
    };
    BitClass.prototype.empty = function() { return this.card === 0; };
    BitClass.prototype.cardinality = function() { return this.card; };
    BitClass.prototype.clone = function() {
        var c = new BitClass(this.n_bits);
        c.card = this.card;
        for (var w = 0; w < this.wordCount; w++) { c.words[w] = this.words[w]; }
        return c;
    };
    BitClass.prototype.nextSetBit = function(from) {
        if (from >= this.n_bits) { return -1; }
        var w = from >>> 5;
        var mask = ~0 << (from & 31);
        var bits = this.words[w] & mask;
        if (bits) { return (w << 5) | ctz32(bits); }
        for (w++; w < this.wordCount; w++) {
            if (this.words[w]) { return (w << 5) | ctz32(this.words[w]); }
        }
        return -1;
    };
    BitClass.prototype.andWith = function(adj, n) {
        var r = new BitClass(n);
        for (var i = this.nextSetBit(0); i >= 0; i = this.nextSetBit(i + 1)) {
            if (adj[i]) { r.set(i); }
        }
        return r;
    };
    BitClass.prototype.andNotWith = function(adj, n) {
        var r = new BitClass(n);
        for (var i = this.nextSetBit(0); i >= 0; i = this.nextSetBit(i + 1)) {
            if (!adj[i]) { r.set(i); }
        }
        return r;
    };

    // Count trailing zeros for 32-bit int
    function ctz32(x) {
        if (x === 0) { return 32; }
        var n = 0;
        if ((x & 0xFFFF) === 0) { n += 16; x >>>= 16; }
        if ((x & 0xFF) === 0) { n += 8; x >>>= 8; }
        if ((x & 0xF) === 0) { n += 4; x >>>= 4; }
        if ((x & 0x3) === 0) { n += 2; x >>>= 2; }
        if ((x & 0x1) === 0) { n += 1; }
        return n;
    }

    function popcount32(x) {
        x = x - ((x >>> 1) & 0x55555555);
        x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
        return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24;
    }

    function checkBondCompat(g1, g2, C, induced, qi, tj, q2t) {
        var nbs = g1.neighbors[qi];
        for (var k = 0; k < nbs.length; k++) {
            var qk = nbs[k];
            if (q2t[qk] < 0) { continue; }
            var tk = q2t[qk];
            var qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
            if (qOrd !== 0 && tOrd !== 0) {
                if (!bondsCompat(g1, qi, qk, g2, tj, tk, C)) { return false; }
            } else if (induced && ((qOrd !== 0) !== (tOrd !== 0))) { return false; }
        }
        return true;
    }

    function mcSplitRecurse(g1, g2, C, induced, tb,
                            qSets, tSets, numClasses,
                            q2t, t2q, curSize, upperBound,
                            bestQ2T, bestSizeRef, nodeCountRef,
                            n1, n2, depth, maxDepth) {
        nodeCountRef.val++;
        if (nodeCountRef.val > MAX_NODE_LIMIT) { return; }
        if ((nodeCountRef.val & 1023) === 0 && tb.expired()) { return; }
        if (curSize > bestSizeRef.val) {
            bestSizeRef.val = curSize;
            for (var cp = 0; cp < n1; cp++) { bestQ2T[cp] = q2t[cp]; }
        }
        if (curSize + upperBound <= bestSizeRef.val || depth >= maxDepth) { return; }

        if (curSize > 0) {
            var canExtend = false;
            for (var ch = 0; ch < numClasses && !canExtend; ch++) {
                if (!qSets[ch].empty() && !tSets[ch].empty()) { canExtend = true; }
            }
            if (!canExtend) { return; }
        }

        // Select most constrained class
        var bestClass = -1, bestMin = 999999999, bestClassQC = 0, bestClassTC = 0;
        for (var c = 0; c < numClasses; c++) {
            var qc = qSets[c].cardinality(), tc = tSets[c].cardinality();
            if (qc === 0 || tc === 0) { continue; }
            var m = Math.min(qc, tc);
            if (m < bestMin) { bestMin = m; bestClass = c; bestClassQC = qc; bestClassTC = tc; }
        }
        if (bestClass === -1) { return; }

        var branchFromTarget = bestClassTC < bestClassQC;

        if (!branchFromTarget) {
            // Select qi with connectivity-aware ordering
            var qi = -1, qiBestScore = -1;
            for (var v = qSets[bestClass].nextSetBit(0); v >= 0; v = qSets[bestClass].nextSetBit(v + 1)) {
                var conn = false;
                var vNbs = g1.neighbors[v];
                for (var vn = 0; vn < vNbs.length; vn++) {
                    if (q2t[vNbs[vn]] >= 0) { conn = true; break; }
                }
                var sc = (conn ? 1000 : 0) + (g1.ring[v] ? 100 : 0) + g1.degree[v] * 10;
                if (sc > qiBestScore) { qiBestScore = sc; qi = v; }
            }

            // Check if qi has any mapped neighbor — if so, orbit symmetry
            // is broken and we must not skip orbit-equivalent targets
            var qiHasMappedNb = false;
            var qiNbsOrbit = g1.neighbors[qi];
            for (var on = 0; on < qiNbsOrbit.length; on++) {
                if (q2t[qiNbsOrbit[on]] >= 0) { qiHasMappedNb = true; break; }
            }

            // RRSplit: try each orbit of target
            var triedOrbits = {};
            for (var tj = tSets[bestClass].nextSetBit(0); tj >= 0; tj = tSets[bestClass].nextSetBit(tj + 1)) {
                nodeCountRef.val++;
                if (nodeCountRef.val > MAX_NODE_LIMIT || ((nodeCountRef.val & 1023) === 0 && tb.expired())) { return; }
                // Only apply orbit pruning when qi has no mapped neighbor
                if (!qiHasMappedNb && triedOrbits[g2.orbit[tj]]) { continue; }
                triedOrbits[g2.orbit[tj]] = true;
                // Fast atom compatibility pre-check before full bond compat
                if (!atomsCompat(g1, qi, g2, tj, C)) { continue; }
                if (!checkBondCompat(g1, g2, C, induced, qi, tj, q2t)) { continue; }

                q2t[qi] = tj; t2q[tj] = qi;
                refineAndRecurse(g1, g2, C, induced, tb, qSets, tSets, numClasses,
                                 q2t, t2q, curSize, bestQ2T, bestSizeRef, nodeCountRef,
                                 n1, n2, depth, maxDepth, qi, tj);
                q2t[qi] = -1; t2q[tj] = -1;
            }

            // Skip qi branch
            mcSplitSkipVertex(g1, g2, C, induced, tb, qSets, tSets, numClasses,
                              q2t, t2q, curSize, bestQ2T, bestSizeRef, nodeCountRef,
                              n1, n2, depth, maxDepth, bestClass, qi, true);
        } else {
            var tj2 = -1, tjBestScore = -1;
            for (var v2 = tSets[bestClass].nextSetBit(0); v2 >= 0; v2 = tSets[bestClass].nextSetBit(v2 + 1)) {
                var conn2 = false;
                var v2Nbs = g2.neighbors[v2];
                for (var vn2 = 0; vn2 < v2Nbs.length; vn2++) {
                    if (t2q[v2Nbs[vn2]] >= 0) { conn2 = true; break; }
                }
                var sc2 = (conn2 ? 1000 : 0) + (g2.ring[v2] ? 100 : 0) + g2.degree[v2] * 10;
                if (sc2 > tjBestScore) { tjBestScore = sc2; tj2 = v2; }
            }

            // Check if tj2 has any mapped neighbor — if so, orbit symmetry
            // on the query side is broken and we must not skip orbit-equivalent qi
            var tj2HasMappedNb = false;
            var tj2NbsOrbit = g2.neighbors[tj2];
            for (var on2 = 0; on2 < tj2NbsOrbit.length; on2++) {
                if (t2q[tj2NbsOrbit[on2]] >= 0) { tj2HasMappedNb = true; break; }
            }

            var triedOrbits2 = {};
            for (var qi2 = qSets[bestClass].nextSetBit(0); qi2 >= 0; qi2 = qSets[bestClass].nextSetBit(qi2 + 1)) {
                nodeCountRef.val++;
                if (nodeCountRef.val > MAX_NODE_LIMIT || ((nodeCountRef.val & 1023) === 0 && tb.expired())) { return; }
                // Only apply orbit pruning when tj2 has no mapped neighbor
                if (!tj2HasMappedNb && triedOrbits2[g1.orbit[qi2]]) { continue; }
                triedOrbits2[g1.orbit[qi2]] = true;
                // Fast atom compatibility pre-check
                if (!atomsCompat(g1, qi2, g2, tj2, C)) { continue; }

                var bondOk = true;
                var tjNbs = g2.neighbors[tj2];
                for (var tn = 0; tn < tjNbs.length; tn++) {
                    var tkn = tjNbs[tn];
                    if (t2q[tkn] < 0) { continue; }
                    var qkn = t2q[tkn];
                    var qOrdN = g1.bondOrder(qi2, qkn), tOrdN = g2.bondOrder(tj2, tkn);
                    if (qOrdN !== 0 && tOrdN !== 0) {
                        if (!bondsCompat(g1, qi2, qkn, g2, tj2, tkn, C)) { bondOk = false; break; }
                    } else if (induced && ((qOrdN !== 0) !== (tOrdN !== 0))) { bondOk = false; break; }
                }
                if (!bondOk) { continue; }

                q2t[qi2] = tj2; t2q[tj2] = qi2;
                refineAndRecurse(g1, g2, C, induced, tb, qSets, tSets, numClasses,
                                 q2t, t2q, curSize, bestQ2T, bestSizeRef, nodeCountRef,
                                 n1, n2, depth, maxDepth, qi2, tj2);
                q2t[qi2] = -1; t2q[tj2] = -1;
            }

            mcSplitSkipVertex(g1, g2, C, induced, tb, qSets, tSets, numClasses,
                              q2t, t2q, curSize, bestQ2T, bestSizeRef, nodeCountRef,
                              n1, n2, depth, maxDepth, bestClass, tj2, false);
        }
    }

    function refineAndRecurse(g1, g2, C, induced, tb,
                              qSets, tSets, numClasses,
                              q2t, t2q, curSize,
                              bestQ2T, bestSizeRef, nodeCountRef,
                              n1, n2, depth, maxDepth, qi, tj) {
        var qiAdj = boolArr(n1), tjAdj = boolArr(n2);
        var qNbs = g1.neighbors[qi], tNbs = g2.neighbors[tj];
        var k;
        for (k = 0; k < qNbs.length; k++) { qiAdj[qNbs[k]] = true; }
        for (k = 0; k < tNbs.length; k++) { tjAdj[tNbs[k]] = true; }

        var newQ = [], newT = [], newUB = 0;
        for (var c = 0; c < numClasses; c++) {
            if (qSets[c].empty() || tSets[c].empty()) { continue; }
            var qAdj = qSets[c].andWith(qiAdj, n1);
            qAdj.clear(qi);
            var qNon = qSets[c].andNotWith(qiAdj, n1);
            qNon.clear(qi);
            var tAdj = tSets[c].andWith(tjAdj, n2);
            tAdj.clear(tj);
            var tNon = tSets[c].andNotWith(tjAdj, n2);
            tNon.clear(tj);

            if (qAdj.cardinality() > 0 && tAdj.cardinality() > 0) {
                newUB += Math.min(qAdj.cardinality(), tAdj.cardinality());
                newQ.push(qAdj); newT.push(tAdj);
            }
            if (qNon.cardinality() > 0 && tNon.cardinality() > 0) {
                newUB += Math.min(qNon.cardinality(), tNon.cardinality());
                newQ.push(qNon); newT.push(tNon);
            }
        }

        if (curSize + 1 + newUB > bestSizeRef.val) {
            mcSplitRecurse(g1, g2, C, induced, tb, newQ, newT, newQ.length,
                           q2t, t2q, curSize + 1, newUB,
                           bestQ2T, bestSizeRef, nodeCountRef,
                           n1, n2, depth + 1, maxDepth);
        }
    }

    function mcSplitSkipVertex(g1, g2, C, induced, tb,
                               qSets, tSets, numClasses,
                               q2t, t2q, curSize,
                               bestQ2T, bestSizeRef, nodeCountRef,
                               n1, n2, depth, maxDepth,
                               bestClass, vertex, isQuery) {
        var skipQ = [], skipT = [], skipUB = 0;
        for (var c = 0; c < numClasses; c++) {
            if (qSets[c].empty() || tSets[c].empty()) { continue; }
            if (c === bestClass) {
                var reduced, reducedT;
                if (isQuery) {
                    reduced = qSets[c].clone(); reduced.clear(vertex);
                    reducedT = tSets[c];
                } else {
                    reduced = qSets[c];
                    reducedT = tSets[c].clone(); reducedT.clear(vertex);
                }
                var checkSet = isQuery ? reduced : reducedT;
                if (!checkSet.empty()) {
                    var qCard = isQuery ? reduced.cardinality() : qSets[c].cardinality();
                    var tCard = isQuery ? tSets[c].cardinality() : reducedT.cardinality();
                    skipUB += Math.min(qCard, tCard);
                    skipQ.push(isQuery ? reduced : qSets[c]);
                    skipT.push(isQuery ? tSets[c] : reducedT);
                }
            } else {
                skipUB += Math.min(qSets[c].cardinality(), tSets[c].cardinality());
                skipQ.push(qSets[c]);
                skipT.push(tSets[c]);
            }
        }
        if (curSize + skipUB > bestSizeRef.val) {
            mcSplitRecurse(g1, g2, C, induced, tb, skipQ, skipT, skipQ.length,
                           q2t, t2q, curSize, skipUB,
                           bestQ2T, bestSizeRef, nodeCountRef,
                           n1, n2, depth + 1, maxDepth);
        }
    }

    function mcSplitSeed(g1, g2, C, induced, tb) {
        var n1 = g1.n, n2 = g2.n;
        if (n1 === 0 || n2 === 0) { return { seed: {}, nodeCount: 0 }; }

        // Build initial partition classes by label
        var qGroups = {}, tGroups = {}, i, lab;
        for (i = 0; i < n1; i++) {
            lab = g1.label[i];
            if (!qGroups[lab]) { qGroups[lab] = []; }
            qGroups[lab].push(i);
        }
        for (i = 0; i < n2; i++) {
            lab = g2.label[i];
            if (!tGroups[lab]) { tGroups[lab] = []; }
            tGroups[lab].push(i);
        }

        // 2-hop orbit refinement: split label groups by 2-hop NLF signature
        // to distinguish atoms that Morgan cannot tell apart
        function twoHopOrbitKey(g, idx) {
            var nbs = g.neighbors[idx];
            var hop2Sig = [];
            for (var h = 0; h < nbs.length; h++) {
                var nb2 = g.neighbors[nbs[h]];
                var cnt = 0;
                for (var h2 = 0; h2 < nb2.length; h2++) {
                    if (nb2[h2] !== idx) { cnt++; }
                }
                hop2Sig.push(mcsNlfLabel(g, nbs[h]) * 100 + cnt);
            }
            hop2Sig.sort(function(a, b) { return a - b; });
            return hop2Sig.join('_');
        }

        function refineGroupByTwoHop(g, members) {
            var buckets = {}, mi;
            for (mi = 0; mi < members.length; mi++) {
                var key = twoHopOrbitKey(g, members[mi]);
                if (!buckets[key]) { buckets[key] = []; }
                buckets[key].push(members[mi]);
            }
            var result = [];
            for (var bk in buckets) {
                if (buckets.hasOwnProperty(bk)) { result.push(buckets[bk]); }
            }
            return result;
        }

        var initQSets = [], initTSets = [];
        var qlabel, tlabel;
        for (qlabel in qGroups) {
            if (!qGroups.hasOwnProperty(qlabel)) { continue; }
            for (tlabel in tGroups) {
                if (!tGroups.hasOwnProperty(tlabel)) { continue; }
                var qList = qGroups[qlabel], tList = tGroups[tlabel];
                if (qList.length === 0 || tList.length === 0) { continue; }
                if (!atomsCompat(g1, qList[0], g2, tList[0], C)) { continue; }

                // Refine by 2-hop orbit signature for tighter initial partitions
                var qSubGroups = refineGroupByTwoHop(g1, qList);
                var tSubGroups = refineGroupByTwoHop(g2, tList);

                for (var qsi = 0; qsi < qSubGroups.length; qsi++) {
                    for (var tsi = 0; tsi < tSubGroups.length; tsi++) {
                        var qSub = qSubGroups[qsi], tSub = tSubGroups[tsi];
                        if (qSub.length === 0 || tSub.length === 0) { continue; }
                        if (!atomsCompat(g1, qSub[0], g2, tSub[0], C)) { continue; }
                        var qbs = new BitClass(n1), tbs = new BitClass(n2);
                        for (var qi = 0; qi < qSub.length; qi++) { qbs.set(qSub[qi]); }
                        for (var ti = 0; ti < tSub.length; ti++) { tbs.set(tSub[ti]); }
                        initQSets.push(qbs); initTSets.push(tbs);
                    }
                }
            }
        }
        if (initQSets.length === 0) { return { seed: {}, nodeCount: 0 }; }

        var numClasses = initQSets.length;
        var initUB = 0;
        for (var c = 0; c < numClasses; c++) {
            initUB += Math.min(initQSets[c].cardinality(), initTSets[c].cardinality());
        }

        var q2t = arrFill(n1, -1), t2q = arrFill(n2, -1);
        var bestQ2T = arrFill(n1, -1);
        var bestSizeRef = { val: 0 };
        var nodeCountRef = { val: 0 };

        mcSplitRecurse(g1, g2, C, induced, tb,
                       initQSets, initTSets, numClasses,
                       q2t, t2q, 0, initUB,
                       bestQ2T, bestSizeRef, nodeCountRef,
                       n1, n2, 0, Math.min(n1, n2) + 1);

        var seed = {};
        for (i = 0; i < n1; i++) {
            if (bestQ2T[i] >= 0) { seed[i] = bestQ2T[i]; }
        }
        return { seed: seed, nodeCount: nodeCountRef.val };
    }

    // =======================================================================
    // L3: Bron-Kerbosch + Tomita pivoting + k-core reduction
    // =======================================================================

    // 32-bit adjacency bitset operations for the product graph
    function pgPopcount(words) {
        var total = 0;
        for (var w = 0; w < words.length; w++) { total += popcount32(words[w]); }
        return total;
    }

    function pgIterateBits(words, N) {
        var result = [];
        for (var w = 0; w < words.length; w++) {
            var bits = words[w];
            while (bits) {
                var bit = ctz32(bits);
                var v = (w << 5) | bit;
                if (v < N) { result.push(v); }
                bits &= (bits - 1);
            }
        }
        return result;
    }

    function greedyClique(adj, N, words) {
        var bestStart = 0, bestDeg = 0, i, w, deg;
        for (i = 0; i < N; i++) {
            deg = pgPopcount(adj[i]);
            if (deg > bestDeg) { bestDeg = deg; bestStart = i; }
        }
        var clique = [bestStart];
        var cand = arrCopy(adj[bestStart]);
        while (true) {
            var bestV = -1, bestConn = -1;
            for (w = 0; w < words; w++) {
                var bits = cand[w];
                while (bits) {
                    var bit = ctz32(bits);
                    var v = (w << 5) | bit;
                    var conn = 0;
                    for (var cw = 0; cw < words; cw++) { conn += popcount32(adj[v][cw] & cand[cw]); }
                    if (conn > bestConn) { bestConn = conn; bestV = v; }
                    bits &= (bits - 1);
                }
            }
            if (bestV === -1) { break; }
            clique.push(bestV);
            for (w = 0; w < words; w++) { cand[w] &= adj[bestV][w]; }
        }
        return clique;
    }

    function colorBound(P, adj, words, N) {
        var color = arrFill(N, -1);
        var maxColor = 0;
        for (var w = 0; w < words; w++) {
            var bits = P[w];
            while (bits) {
                var bit = ctz32(bits);
                var v = (w << 5) | bit;
                var usedColors = 0;
                for (var nw = 0; nw < words; nw++) {
                    var nbits = adj[v][nw] & P[nw];
                    while (nbits) {
                        var nbit = ctz32(nbits);
                        var u = (nw << 5) | nbit;
                        if (color[u] >= 0 && color[u] < 32) { usedColors |= (1 << color[u]); }
                        nbits &= (nbits - 1);
                    }
                }
                color[v] = ctz32(~usedColors);
                if (color[v] > maxColor) { maxColor = color[v]; }
                bits &= (bits - 1);
            }
        }
        return maxColor + 1;
    }

    function bronKerboschPivot(R, P, X, adj, currentBestRef, N, words,
                               equivClass, tb, nodes, g1, g2, depth) {
        if (tb.expired()) { return; }
        var rSize = 0, pSize = 0, xSize = 0, w;
        for (w = 0; w < words; w++) {
            rSize += popcount32(R[w]);
            pSize += popcount32(P[w]);
            xSize += popcount32(X[w]);
        }
        if (pSize === 0 && xSize === 0) {
            if (rSize > currentBestRef.best.length) {
                currentBestRef.best = pgIterateBits(R, N);
            }
            return;
        }
        if (rSize + pSize <= currentBestRef.best.length) { return; }

        // Partition bound
        if (pSize > 2) {
            var qiSeen = boolArr(g1.n), tjSeen = boolArr(g2.n);
            for (w = 0; w < words; w++) {
                var pbits = P[w];
                while (pbits) {
                    var pbit = ctz32(pbits);
                    var pv = (w << 5) | pbit;
                    if (pv < N) { qiSeen[nodes[pv].qi] = true; tjSeen[nodes[pv].tj] = true; }
                    pbits &= (pbits - 1);
                }
            }
            var qiCard = 0, tjCard = 0;
            for (var gi = 0; gi < g1.n; gi++) { if (qiSeen[gi]) { qiCard++; } }
            for (var gj = 0; gj < g2.n; gj++) { if (tjSeen[gj]) { tjCard++; } }
            if (rSize + Math.min(qiCard, tjCard) <= currentBestRef.best.length) { return; }
        }

        if (pSize > 0 && rSize + colorBound(P, adj, words, N) <= currentBestRef.best.length) { return; }

        // Choose pivot
        var pivot = -1, pivotConn = -1;
        for (w = 0; w < words; w++) {
            var pbits2 = P[w] | X[w];
            while (pbits2) {
                var pbit2 = ctz32(pbits2);
                var pu = (w << 5) | pbit2;
                var pconn = 0;
                for (var cw = 0; cw < words; cw++) { pconn += popcount32(P[cw] & adj[pu][cw]); }
                if (pconn > pivotConn) { pivotConn = pconn; pivot = pu; }
                pbits2 &= (pbits2 - 1);
            }
        }

        var candidates = arrFill(words, 0);
        if (pivot >= 0) {
            for (w = 0; w < words; w++) { candidates[w] = P[w] & ~adj[pivot][w]; }
        } else {
            for (w = 0; w < words; w++) { candidates[w] = P[w]; }
        }

        var triedClasses = {};
        for (w = 0; w < words; w++) {
            var cbits = candidates[w];
            while (cbits) {
                if (tb.expired()) { return; }
                var cbit = ctz32(cbits);
                var cv = (w << 5) | cbit;
                cbits &= (cbits - 1);
                if (triedClasses[equivClass[cv]]) { continue; }
                triedClasses[equivClass[cv]] = true;

                var newR = arrFill(words, 0), newP = arrFill(words, 0), newX = arrFill(words, 0);
                for (var cw2 = 0; cw2 < words; cw2++) {
                    newR[cw2] = R[cw2];
                    newP[cw2] = P[cw2] & adj[cv][cw2];
                    newX[cw2] = X[cw2] & adj[cv][cw2];
                }
                newR[cv >>> 5] |= (1 << (cv & 31));
                bronKerboschPivot(newR, newP, newX, adj, currentBestRef, N, words,
                                  equivClass, tb, nodes, g1, g2, depth + 1);
                P[cv >>> 5] &= ~(1 << (cv & 31));
                X[cv >>> 5] |= (1 << (cv & 31));
            }
        }
    }

    function maximumCliqueSeed(g1, g2, C, induced, tb) {
        var n1 = g1.n, n2 = g2.n;
        var nodes = [];
        for (var i = 0; i < n1; i++) {
            for (var j = 0; j < n2; j++) {
                if (!atomsCompat(g1, i, g2, j, C)) { continue; }
                if (!induced && g1.degree[i] > g2.degree[j]) { continue; }
                nodes.push({ qi: i, tj: j });
            }
        }
        var N = nodes.length;
        if (N === 0) { return {}; }

        var words = ((N + 31) >>> 5);
        var adj = new Array(N);
        for (var u = 0; u < N; u++) { adj[u] = arrFill(words, 0); }

        for (u = 0; u < N; u++) {
            var nu = nodes[u];
            for (var v = u + 1; v < N; v++) {
                if (tb.expired()) { break; }
                var nv = nodes[v];
                if (nu.qi === nv.qi || nu.tj === nv.tj) { continue; }
                var qOrd = g1.bondOrder(nu.qi, nv.qi), tOrd = g2.bondOrder(nu.tj, nv.tj);
                var ok;
                if (qOrd !== 0 && tOrd !== 0) {
                    ok = bondsCompat(g1, nu.qi, nv.qi, g2, nu.tj, nv.tj, C);
                } else if (induced) {
                    ok = (qOrd === 0 && tOrd === 0);
                } else {
                    continue;
                }
                if (ok) {
                    adj[u][v >>> 5] |= (1 << (v & 31));
                    adj[v][u >>> 5] |= (1 << (u & 31));
                }
            }
        }

        // Product graph size cap: skip BK entirely for very large graphs
        if (N > 1500) { return {}; }

        // Greedy clique
        var bestClique = greedyClique(adj, N, words);

        // K-core reduction: skip for large product graphs (pgSize > 2500)
        var pgSize = g1.n * g2.n;
        if (bestClique.length > 1 && !tb.expired() && pgSize <= 2500) {
            var minDeg = Math.max(2, bestClique.length);
            var changed = true;
            while (changed && !tb.expired()) {
                changed = false;
                for (u = 0; u < N; u++) {
                    var deg = pgPopcount(adj[u]);
                    if (deg === 0 || deg >= minDeg) { continue; }
                    for (var w = 0; w < words; w++) {
                        var bits = adj[u][w];
                        while (bits) {
                            var bit = ctz32(bits);
                            var vv = (w << 5) | bit;
                            adj[vv][u >>> 5] &= ~(1 << (u & 31));
                            bits &= (bits - 1);
                        }
                        adj[u][w] = 0;
                    }
                    changed = true;
                }
            }
        }

        // Equivalence classes
        var equivClass = new Array(N);
        var sig2class = {};
        var nextClass = 0;
        for (i = 0; i < N; i++) {
            var nd = nodes[i];
            var ndDeg = pgPopcount(adj[i]);
            var sig = g1.orbit[nd.qi] + '_' + g2.orbit[nd.tj] + '_'
                    + (g1.morganRank ? (g1.morganRank[nd.qi] & 0xFF) : 0) + '_'
                    + g1.degree[nd.qi] + '_' + ndDeg;
            if (sig2class[sig] !== undefined) {
                equivClass[i] = sig2class[sig];
            } else {
                sig2class[sig] = nextClass;
                equivClass[i] = nextClass++;
            }
        }

        // BK search
        var P = arrFill(words, 0), X = arrFill(words, 0), R = arrFill(words, 0);
        for (i = 0; i < N; i++) {
            var ideg = pgPopcount(adj[i]);
            if (ideg > 0) { P[i >>> 5] |= (1 << (i & 31)); }
        }

        var currentBestRef = { best: bestClique };
        bronKerboschPivot(R, P, X, adj, currentBestRef, N, words,
                          equivClass, tb, nodes, g1, g2, 0);

        var seed = {};
        var usedQSet = {}, usedTSet = {};
        for (var ci = 0; ci < currentBestRef.best.length; ci++) {
            var ndv = nodes[currentBestRef.best[ci]];
            if (!usedQSet[ndv.qi] && !usedTSet[ndv.tj]) {
                seed[ndv.qi] = ndv.tj;
                usedQSet[ndv.qi] = true;
                usedTSet[ndv.tj] = true;
            }
        }
        return seed;
    }

    // =======================================================================
    // L5: Extra seed generators
    // =======================================================================

    function ringAnchorSeed(g1, g2, C, tb) {
        var seed = {};
        var R1 = [], R2 = [], i;
        for (i = 0; i < g1.n; i++) { if (g1.ring[i]) { R1.push(i); } }
        for (i = 0; i < g2.n; i++) { if (g2.ring[i]) { R2.push(i); } }
        var used = boolArr(g2.n);
        for (var ri = 0; ri < R1.length; ri++) {
            var qi = R1[ri];
            var bestJ = -1, bestScore = -999999;
            for (var rj = 0; rj < R2.length; rj++) {
                var j = R2[rj];
                if (used[j]) { continue; }
                if (!atomsCompat(g1, qi, g2, j, C)) { continue; }
                var score = (g1.aromatic[qi] && g2.aromatic[j] ? 10 : 0)
                          + Math.min(g1.degree[qi], g2.degree[j]);
                if (score > bestScore) { bestScore = score; bestJ = j; }
            }
            if (bestJ >= 0) { seed[qi] = bestJ; used[bestJ] = true; }
            if (tb.expired()) { break; }
        }
        return seed;
    }

    function labelDegreeAnchorSeed(g1, g2, C, tb) {
        var seed = {};
        var freq = {};
        var i;
        for (i = 0; i < g1.n; i++) {
            freq[g1.label[i]] = (freq[g1.label[i]] || 0) + 1;
        }
        var Q = [];
        for (i = 0; i < g1.n; i++) { Q.push(i); }
        Q.sort(function(a, b) {
            var ra = freq[g1.label[a]] || 0, rb = freq[g1.label[b]] || 0;
            return ra !== rb ? ra - rb : g1.degree[b] - g1.degree[a];
        });
        var used = boolArr(g2.n);
        for (var qi = 0; qi < Q.length; qi++) {
            var q = Q[qi];
            var bestJ = -1, bestScore = -999999;
            for (var tj = 0; tj < g2.n; tj++) {
                if (used[tj]) { continue; }
                if (!atomsCompat(g1, q, g2, tj, C)) { continue; }
                var score = 10 - Math.abs(g1.degree[q] - g2.degree[tj]);
                if (score > bestScore) { bestScore = score; bestJ = tj; }
            }
            if (bestJ >= 0) { seed[q] = bestJ; used[bestJ] = true; }
            if (tb.expired()) { break; }
        }
        return seed;
    }

    // =======================================================================
    // L0.25: Linear chain fast-path
    // =======================================================================
    function linearChainFastPath(g1, g2, C) {
        var g1Chain = true, g2Chain = true, i;
        for (i = 0; i < g1.n && g1Chain; i++) { if (g1.degree[i] > 2) { g1Chain = false; } }
        for (i = 0; i < g2.n && g2Chain; i++) { if (g2.degree[i] > 2) { g2Chain = false; } }
        if (!g1Chain || !g2Chain || g1.n < 2 || g2.n < 2) { return null; }

        var walkChain = function(g) {
            var start = -1;
            for (var s = 0; s < g.n; s++) {
                if (g.degree[s] <= 1) { start = s; break; }
            }
            if (start < 0) { start = 0; }
            var seq = [];
            var visited = boolArr(g.n);
            var cur = start;
            while (cur >= 0) {
                visited[cur] = true;
                seq.push(cur);
                var next = -1;
                var nbs = g.neighbors[cur];
                for (var nb = 0; nb < nbs.length; nb++) {
                    if (!visited[nbs[nb]]) { next = nbs[nb]; break; }
                }
                cur = next;
            }
            return seq;
        };

        var seq1 = walkChain(g1);
        var seq2 = walkChain(g2);
        var n1 = seq1.length, n2 = seq2.length;

        // DP: longest common subpath
        var dp = new Array(n1 + 1);
        for (i = 0; i <= n1; i++) { dp[i] = arrFill(n2 + 1, 0); }
        var bestLen = 0, bestI = 0, bestJ = 0;
        for (i = 1; i <= n1; i++) {
            var a1 = seq1[i - 1];
            for (var j = 1; j <= n2; j++) {
                var a2 = seq2[j - 1];
                if (!atomsCompat(g1, a1, g2, a2, C)) { continue; }
                if (i > 1 && j > 1) {
                    var b1prev = seq1[i - 2], b2prev = seq2[j - 2];
                    if (dp[i - 1][j - 1] > 0 &&
                        bondsCompat(g1, b1prev, a1, g2, b2prev, a2, C)) {
                        dp[i][j] = dp[i - 1][j - 1] + 1;
                    } else {
                        dp[i][j] = 1;
                    }
                } else {
                    dp[i][j] = 1;
                }
                if (dp[i][j] > bestLen) { bestLen = dp[i][j]; bestI = i; bestJ = j; }
            }
        }
        if (bestLen <= 0) { return null; }
        var chainMCS = {};
        for (var k = 0; k < bestLen; k++) {
            chainMCS[seq1[bestI - bestLen + k]] = seq2[bestJ - bestLen + k];
        }
        return chainMCS;
    }

    // =======================================================================
    // L0.5: Tree fast-path
    // =======================================================================
    function treeFastPath(g1, g2, C) {
        if (g1.n < 10 || g2.n < 10) { return null; }

        var isTree = function(g) {
            var edgeCount = 0;
            for (var i = 0; i < g.n; i++) { edgeCount += g.neighbors[i].length; }
            edgeCount = edgeCount / 2;
            if (edgeCount !== g.n - 1) { return false; }
            var vis = boolArr(g.n);
            var bfs = [0]; vis[0] = true;
            var cnt = 1, head = 0;
            while (head < bfs.length) {
                var u = bfs[head++];
                var nbs = g.neighbors[u];
                for (var k = 0; k < nbs.length; k++) {
                    if (!vis[nbs[k]]) { vis[nbs[k]] = true; cnt++; bfs.push(nbs[k]); }
                }
            }
            return cnt === g.n;
        };

        if (!isTree(g1) || !isTree(g2)) { return null; }

        var findCentroid = function(g) {
            if (g.n <= 2) { return 0; }
            var deg = new Array(g.n);
            var leaves = [];
            for (var i = 0; i < g.n; i++) {
                deg[i] = g.neighbors[i].length;
                if (deg[i] <= 1) { leaves.push(i); }
            }
            var remaining = g.n;
            while (remaining > 2) {
                var sz = leaves.length;
                remaining -= sz;
                var newLeaves = [];
                for (var li = 0; li < sz; li++) {
                    var u = leaves[li];
                    var nbs = g.neighbors[u];
                    for (var k = 0; k < nbs.length; k++) {
                        deg[nbs[k]]--;
                        if (deg[nbs[k]] === 1) { newLeaves.push(nbs[k]); }
                    }
                }
                leaves = newLeaves;
            }
            return leaves[0];
        };

        var rootTree = function(g, root) {
            var parent = arrFill(g.n, -1);
            var children = new Array(g.n);
            for (var ci = 0; ci < g.n; ci++) { children[ci] = []; }
            var vis = boolArr(g.n);
            var bfs = [root]; vis[root] = true;
            var bfsOrder = [];
            var head = 0;
            while (head < bfs.length) {
                var u = bfs[head++];
                bfsOrder.push(u);
                var nbs = g.neighbors[u];
                for (var k = 0; k < nbs.length; k++) {
                    var v = nbs[k];
                    if (!vis[v]) {
                        vis[v] = true;
                        parent[v] = u;
                        children[u].push(v);
                        bfs.push(v);
                    }
                }
            }
            var postorder = [];
            for (var pi = bfsOrder.length - 1; pi >= 0; pi--) { postorder.push(bfsOrder[pi]); }
            return { parent: parent, children: children, postorder: postorder, root: root };
        };

        var root1 = findCentroid(g1);
        var root2 = findCentroid(g2);
        var rt1 = rootTree(g1, root1);
        var rt2 = rootTree(g2, root2);

        var n1t = g1.n, n2t = g2.n;
        // dp[u][v] = max common subtree size
        var dpTree = new Array(n1t);
        var btTree = new Array(n1t);
        for (var di = 0; di < n1t; di++) {
            dpTree[di] = arrFill(n2t, 0);
            btTree[di] = new Array(n2t);
            for (var dj = 0; dj < n2t; dj++) { btTree[di][dj] = []; }
        }

        for (var po1 = 0; po1 < rt1.postorder.length; po1++) {
            var u = rt1.postorder[po1];
            for (var po2 = 0; po2 < rt2.postorder.length; po2++) {
                var v = rt2.postorder[po2];
                if (!atomsCompat(g1, u, g2, v, C)) {
                    dpTree[u][v] = 0;
                    continue;
                }
                var ch1 = rt1.children[u];
                var ch2 = rt2.children[v];
                if (ch1.length === 0 || ch2.length === 0) {
                    dpTree[u][v] = 1;
                    continue;
                }
                // Greedy bipartite matching by descending dp value
                var d1 = ch1.length, d2 = ch2.length;
                var edges = [];
                for (var ci1 = 0; ci1 < d1; ci1++) {
                    for (var ci2 = 0; ci2 < d2; ci2++) {
                        var w = dpTree[ch1[ci1]][ch2[ci2]];
                        if (w > 0 && bondsCompat(g1, u, ch1[ci1], g2, v, ch2[ci2], C)) {
                            edges.push({ weight: w, i: ci1, j: ci2 });
                        }
                    }
                }
                edges.sort(function(a, b) { return b.weight - a.weight; });
                var usedI = boolArr(d1), usedJ = boolArr(d2);
                var childSum = 0;
                var matching = [];
                for (var ei = 0; ei < edges.length; ei++) {
                    var e = edges[ei];
                    if (usedI[e.i] || usedJ[e.j]) { continue; }
                    usedI[e.i] = true; usedJ[e.j] = true;
                    childSum += e.weight;
                    matching.push([e.i, e.j]);
                }
                dpTree[u][v] = 1 + childSum;
                btTree[u][v] = matching;
            }
        }

        // Find best root pair
        var treeBest = 0, bestU = -1, bestV = -1;
        for (u = 0; u < n1t; u++) {
            for (v = 0; v < n2t; v++) {
                if (dpTree[u][v] > treeBest) { treeBest = dpTree[u][v]; bestU = u; bestV = v; }
            }
        }
        if (treeBest <= 0) { return null; }

        var treeMap = {};
        var reconstruct = function(u, v) {
            treeMap[u] = v;
            var bt = btTree[u][v];
            for (var bi = 0; bi < bt.length; bi++) {
                reconstruct(rt1.children[u][bt[bi][0]], rt2.children[v][bt[bi][1]]);
            }
        };
        reconstruct(bestU, bestV);
        return treeMap;
    }

    // =======================================================================
    // RASCAL screening
    // =======================================================================
    function similarityUpperBound(g1, g2) {
        if (g1.n === 0 || g2.n === 0) { return 0.0; }
        var count1 = {}, count2 = {}, i;
        for (i = 0; i < g1.n; i++) { count1[g1.label[i]] = (count1[g1.label[i]] || 0) + 1; }
        for (i = 0; i < g2.n; i++) { count2[g2.label[i]] = (count2[g2.label[i]] || 0) + 1; }
        var atomOverlap = 0, lab;
        for (lab in count1) {
            if (count1.hasOwnProperty(lab) && count2.hasOwnProperty(lab)) {
                atomOverlap += Math.min(count1[lab], count2[lab]);
            }
        }
        if (atomOverlap === 0) { return 0.0; }
        var total = g1.n + g2.n;
        return atomOverlap / (total - atomOverlap);
    }

    // =======================================================================
    // Main: findMCS — 7-level coverage-driven funnel
    // =======================================================================
    function findMCS(g1, g2, chemOpts, mcsOpts) {
        var C = chemOpts || { matchAtomType: true, matchBondType: true };
        var M = mergeOpts(McsOptions, mcsOpts);

        // Progressive callback support
        var onProgress = (mcsOpts && typeof mcsOpts.onProgress === 'function')
            ? mcsOpts.onProgress : null;
        var progressPhases = [0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00];
        var phaseIdx = 0;

        function reportProgress(best, bestSize, n1, n2, phaseName) {
            if (!onProgress) { return; }
            if (phaseIdx >= progressPhases.length) { return; }
            var pct = progressPhases[phaseIdx];
            phaseIdx++;
            onProgress({
                phase: phaseName,
                budget: pct,
                bestSoFar: {
                    mapping: best,
                    size: bestSize,
                    tanimoto: tanimoto(bestSize, n1, n2),
                    overlap: simpson(bestSize, n1, n2)
                }
            });
        }

        // Adaptive timeout: scale with molecular complexity, capped at 30s
        var adaptiveMs = Math.min(30000, 500 + g1.n * g2.n * 2);
        var timeoutMs = (mcsOpts && mcsOpts.timeoutMs !== undefined)
            ? mcsOpts.timeoutMs
            : adaptiveMs;
        M.timeoutMs = timeoutMs;
        var tb = new TimeBudget(timeoutMs);

        // Ensure orbit arrays exist (identity orbits by default)
        if (!g1.orbit) {
            g1.orbit = new Array(g1.n);
            for (var i = 0; i < g1.n; i++) g1.orbit[i] = i;
        }
        if (!g2.orbit) {
            g2.orbit = new Array(g2.n);
            for (var i = 0; i < g2.n; i++) g2.orbit[i] = i;
        }

        // Edge case: empty molecules
        if (g1.n === 0 || g2.n === 0) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
        }

        // Edge case: single-atom molecules
        if (g1.n === 1 && g2.n === 1) {
            if (atomsCompat(g1, 0, g2, 0, C)) {
                return { mapping: { 0: 0 }, size: 1, tanimoto: 1.0, overlap: 1.0 };
            }
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
        }

        // Ensure canonical data if available
        if (g1.ensureCanonical) { g1.ensureCanonical(); }
        if (g2.ensureCanonical) { g2.ensureCanonical(); }

        // DSB (degree-sequence bound) is only safe for induced MCS;
        // non-induced uses the looser label-frequency upper bound (LFUB)
        var lfub = labelFrequencyUpperBound(g1, g2, C);
        var upperBound = M.induced
            ? Math.min(lfub, degreeSequenceUpperBound(g1, g2, C))
            : lfub;
        var best = {};
        var bestSize = 0, bestScore = 0;
        var weightMode = M.maximizeBonds;

        // Identity check
        if (g1 === g2 || (g1.n === g2.n && g1.n > 0 &&
            g1.canonicalHash !== undefined && g1.canonicalHash === g2.canonicalHash &&
            g1.canonicalLabel && g2.canonicalLabel)) {
            var id = {};
            var idValid = true;
            if (g1 === g2) {
                for (var ii = 0; ii < g1.n; ii++) { id[ii] = ii; }
            } else {
                // Map via canonical labels — guard against hash collision
                var head = arrFill(g2.n, -1), nxt = arrFill(g2.n, -1);
                for (var j = g2.n - 1; j >= 0; j--) {
                    nxt[j] = head[g2.canonicalLabel[j]];
                    head[g2.canonicalLabel[j]] = j;
                }
                for (var ii2 = 0; ii2 < g1.n; ii2++) {
                    var cl = g1.canonicalLabel[ii2];
                    var tgt = head[cl];
                    if (tgt < 0) {
                        // Canonical label mismatch — hash collision; abort
                        idValid = false;
                        break;
                    }
                    // Atom-compatibility cross-check (catches collisions where
                    // canonicalLabel coincides but atoms differ, e.g. charge)
                    if (!atomsCompat(g1, ii2, g2, tgt, C)) {
                        idValid = false;
                        break;
                    }
                    id[ii2] = tgt;
                    head[cl] = nxt[head[cl]];
                }
            }
            if (idValid) {
                if (M.disconnectedMCS && (M.minFragmentSize > 1 || M.maxFragments < 999999)) {
                    id = applyFragmentConstraints(g1, id, M.minFragmentSize, M.maxFragments);
                }
                var idSize = mapSize(id);
                if (idSize === g1.n || idSize > 0) {
                    return { mapping: id, size: idSize, tanimoto: tanimoto(idSize, g1.n, g2.n), overlap: simpson(idSize, g1.n, g2.n) };
                }
            }
            // fall through to normal MCS pipeline on collision
        }

        var minN = Math.min(g1.n, g2.n);

        // --- L0.25: Linear chain fast-path ---
        if (!weightMode) {
            var chainResult = linearChainFastPath(g1, g2, C);
            if (chainResult) {
                var chainSize = mapSize(chainResult);
                var chainScore = mcsScore(g1, chainResult, M);
                if (!weightMode ? chainSize >= bestSize : chainScore > bestScore) {
                    best = chainResult; bestSize = chainSize; bestScore = chainScore;
                }
                if (bestSize >= upperBound) {
                    return formatResult(best, g1.n, g2.n);
                }
            }
        }

        // --- L0.5: Tree fast-path ---
        if (!weightMode) {
            var treeResult = treeFastPath(g1, g2, C);
            if (treeResult) {
                var treeSize = mapSize(treeResult);
                var treeScore = mcsScore(g1, treeResult, M);
                if (treeSize > bestSize) {
                    best = treeResult; bestSize = treeSize; bestScore = treeScore;
                }
                if (bestSize >= upperBound) {
                    return formatResult(best, g1.n, g2.n);
                }
            }
        }

        // --- L0.75: Greedy probe ---
        if (!weightMode && minN > 2 && minN < GREEDY_PROBE_MAX_SIZE && upperBound >= minN) {
            var greedy = greedyProbe(g1, g2, C, M.templateFuzzyAtoms);
            var greedySize = mapSize(greedy);
            if (greedySize >= upperBound) {
                var greedyC = applyRingAnchorGuard(g1, g2,
                    M.connectedOnly ? largestConnected(g1, greedy) : greedy, C);
                if (mapSize(greedyC) >= upperBound) {
                    return formatResult(greedyC, g1.n, g2.n);
                }
            }
            if (greedySize > upperBound * 0.6 && greedySize > bestSize) {
                var greedyC2 = applyRingAnchorGuard(g1, g2,
                    M.connectedOnly ? largestConnected(g1, greedy) : greedy, C);
                if (mapSize(greedyC2) > bestSize) {
                    best = greedyC2; bestSize = mapSize(greedyC2);
                    bestScore = mcsScore(g1, best, M);
                }
            }
        }

        // Report after greedy phase
        reportProgress(best, bestSize, g1.n, g2.n, 'greedy');

        // --- L1.25: Augmenting path refinement ---
        if (bestSize > 0 && bestSize < upperBound && !tb.expired()) {
            var augmented = mapCopy(best);
            var grew = true;
            while (grew && !tb.expired()) {
                grew = false;
                var augUsedQ = boolArr(g1.n), augUsedT = boolArr(g2.n);
                var aq;
                for (aq in augmented) {
                    if (augmented.hasOwnProperty(aq)) {
                        augUsedQ[+aq] = true; augUsedT[augmented[aq]] = true;
                    }
                }
                for (aq in augmented) {
                    if (!augmented.hasOwnProperty(aq)) { continue; }
                    var aqi = +aq, ati = augmented[aq];
                    var aqNbs = g1.neighbors[aqi];
                    for (var an = 0; an < aqNbs.length; an++) {
                        var aqk = aqNbs[an];
                        if (augUsedQ[aqk]) { continue; }
                        var candidate = -1, count = 0;
                        var atNbs = g2.neighbors[ati];
                        for (var at = 0; at < atNbs.length; at++) {
                            var atk = atNbs[at];
                            if (augUsedT[atk]) { continue; }
                            if (atomsCompat(g1, aqk, g2, atk, C) &&
                                bondsCompat(g1, aqi, aqk, g2, ati, atk, C)) {
                                candidate = atk; count++;
                            }
                        }
                        if (count === 1) {
                            var consistent = true;
                            var aqkNbs = g1.neighbors[aqk];
                            for (var acn = 0; acn < aqkNbs.length; acn++) {
                                var aqm = aqkNbs[acn];
                                if (!augmented.hasOwnProperty(aqm)) { continue; }
                                var atm = augmented[aqm];
                                var aqOrd = g1.bondOrder(aqk, aqm), atOrd = g2.bondOrder(candidate, atm);
                                if ((aqOrd !== 0) !== (atOrd !== 0)) { consistent = false; break; }
                                if (aqOrd !== 0 && !bondsCompat(g1, aqk, aqm, g2, candidate, atm, C)) {
                                    consistent = false; break;
                                }
                            }
                            if (consistent) {
                                augmented[aqk] = candidate;
                                augUsedQ[aqk] = true; augUsedT[candidate] = true;
                                grew = true;
                            }
                        }
                    }
                }
            }
            var augSize = mapSize(augmented);
            if (augSize > bestSize) {
                best = augmented; bestSize = augSize;
                bestScore = mcsScore(g1, best, M);
            }
            if (!weightMode && bestSize >= upperBound) {
                return formatResult(best, g1.n, g2.n);
            }
        }

        // --- L1: Substructure containment ---
        if (!weightMode && g1.n > 0 && g2.n > 0 &&
            Math.min(g1.n, g2.n) <= Math.max(g1.n, g2.n) * 3 / 4 &&
            global.SMSDVF2 && !tb.expired()) {
            var sml = g1.n <= g2.n ? g1 : g2;
            var lrg = g1.n <= g2.n ? g2 : g1;
            var swapped = g1.n > g2.n;
            var subMap = global.SMSDVF2.findSubstructure(sml, lrg, C, tb.remainingMs());
            if (subMap && mapSize(subMap) > 0) {
                var result;
                if (!swapped) {
                    result = subMap;
                } else {
                    result = {};
                    for (var sk in subMap) {
                        if (subMap.hasOwnProperty(sk)) { result[subMap[sk]] = +sk; }
                    }
                }
                // Run through ppx so connectedOnly / induced / fragment
                // constraints are honoured (fixes disconnected-query bug)
                result = ppx(g1, g2, result, C, M);
                if (mapSize(result) > 0) {
                    return formatResult(result, g1.n, g2.n);
                }
            }
        }

        // Report after seed-extend phase
        reportProgress(best, bestSize, g1.n, g2.n, 'seed-extend');

        // --- L1.5: Seed-and-extend ---
        if (minN >= 4 && Math.max(g1.n, g2.n) <= SEED_EXTEND_MAX_ATOMS && !tb.expired()) {
            var seSeed = seedExtendMCS(g1, g2, C, tb, upperBound, M.induced);
            var seScore = mcsScore(g1, seSeed, M);
            var seSize = mapSize(seSeed);
            if (isBetterMapping(g1, seSeed, seSize, best, bestSize, weightMode, seScore, bestScore)) {
                best = seSeed; bestSize = seSize; bestScore = seScore;
            }
            if (!weightMode && bestSize >= upperBound) {
                return formatResult(ppx(g1, g2, best, C, M), g1.n, g2.n);
            }
        }

        // Report after seed-extend actual
        reportProgress(best, bestSize, g1.n, g2.n, 'seed-extend-done');

        // --- L2: McSplit ---
        var mcSplitExhaustive = false;
        var mcSplitSize = 0;
        if (!tb.expired()) {
            var mcResult = mcSplitSeed(g1, g2, C, M.induced, tb);
            var mcSeed = mcResult.seed;
            mcSplitSize = mapSize(mcSeed);
            var mcScore = mcsScore(g1, mcSeed, M);
            if (isBetterMapping(g1, mcSeed, mcSplitSize, best, bestSize, weightMode, mcScore, bestScore)) {
                best = mcSeed; bestSize = mcSplitSize; bestScore = mcScore;
            }
            mcSplitExhaustive = mcResult.nodeCount < 200000;
        }
        if (!weightMode && mcSplitExhaustive && bestSize >= upperBound) {
            return formatResult(ppx(g1, g2, best, C, M), g1.n, g2.n);
        }

        // Near-optimal fast path: greedy atom extension
        if (bestSize >= upperBound - 2 && bestSize > 0 && bestSize < upperBound && !tb.expired()) {
            var ext = ppx(g1, g2, greedyAtomExtend(g1, g2, ppx(g1, g2, best, C, M), C, M), C, M);
            var extScore = mcsScore(g1, ext, M);
            var extSize = mapSize(ext);
            if (isBetterMapping(g1, ext, extSize, best, bestSize, weightMode, extScore, bestScore)) {
                best = ext; bestSize = extSize; bestScore = extScore;
            }
            if (!weightMode && bestSize >= upperBound) {
                return formatResult(best, g1.n, g2.n);
            }
        }

        // Report after McSplit phase
        reportProgress(best, bestSize, g1.n, g2.n, 'mcsplit');

        // --- L3: Bron-Kerbosch ---
        var bkSize = 0;
        if (bestSize < Math.floor(upperBound * BK_SKIP_RATIO) && !tb.expired()) {
            var cliqueSeed = maximumCliqueSeed(g1, g2, C, M.induced, tb);
            bkSize = mapSize(cliqueSeed);
            var cScore = mcsScore(g1, cliqueSeed, M);
            if (isBetterMapping(g1, cliqueSeed, bkSize, best, bestSize, weightMode, cScore, bestScore)) {
                best = cliqueSeed; bestSize = bkSize; bestScore = cScore;
            }
        }
        if (!weightMode && bestSize >= upperBound) {
            return formatResult(ppx(g1, g2, best, C, M), g1.n, g2.n);
        }

        // Report after BK phase
        reportProgress(best, bestSize, g1.n, g2.n, 'bk');

        // --- L4: McGregor extension with seeds ---
        var seeds = [];
        if (bestSize > 0) { seeds.push(mapCopy(best)); }

        // --- L5: Extra seeds (only when McSplit and BK disagree) ---
        var skipExtraSeeds = mcSplitSize > 0 && bkSize > 0 && mcSplitSize === bkSize;
        if (M.extraSeeds && !skipExtraSeeds && bestSize < upperBound && !tb.expired()) {
            var s = ringAnchorSeed(g1, g2, C, tb);
            if (mapSize(s) > 0) { seeds.push(s); }
            if (!tb.expired()) {
                s = labelDegreeAnchorSeed(g1, g2, C, tb);
                if (mapSize(s) > 0) { seeds.push(s); }
            }
        }

        var perSeedMs = Math.max(1, Math.floor(tb.remainingMs() / Math.max(1, seeds.length)));
        for (var si = 0; si < seeds.length; si++) {
            if (tb.expired()) { break; }
            var extSeed = ppx(g1, g2,
                mcGregorExtend(g1, g2, seeds[si], C, tb, perSeedMs,
                               M.useTwoHopNLF, M.useThreeHopNLF, M.connectedOnly),
                C, M);
            var seedScore = mcsScore(g1, extSeed, M);
            var seedSize = mapSize(extSeed);
            if (isBetterMapping(g1, extSeed, seedSize, best, bestSize, weightMode, seedScore, bestScore)) {
                best = extSeed; bestSize = seedSize; bestScore = seedScore;
            }
            if (!weightMode && bestSize >= upperBound) {
                return formatResult(best, g1.n, g2.n);
            }
        }

        // Report after McGregor phase
        reportProgress(best, bestSize, g1.n, g2.n, 'mcgregor');

        // Last resort: start from empty seed
        if (bestScore <= 0 && !tb.expired()) {
            best = ppx(g1, g2,
                mcGregorExtend(g1, g2, {}, C, tb, tb.remainingMs(),
                               M.useTwoHopNLF, M.useThreeHopNLF, M.connectedOnly),
                C, M);
        }

        return formatResult(best, g1.n, g2.n);
    }

    // -----------------------------------------------------------------------
    // Tanimoto and formatting
    // -----------------------------------------------------------------------
    function tanimoto(mcsSize, n1, n2) {
        if (n1 === 0 && n2 === 0) { return 1.0; }
        if (mcsSize === 0) { return 0.0; }
        var denom = n1 + n2 - mcsSize;
        if (denom <= 0) { return 0.0; }
        return mcsSize / denom;
    }

    function simpson(mcsSize, n1, n2) {
        // Both empty: vacuously fully overlapping (consistent with tanimoto).
        if (n1 === 0 && n2 === 0) { return 1.0; }
        var minN = Math.min(n1, n2);
        if (minN === 0) { return 0.0; }
        return mcsSize / minN;
    }

    function formatResult(mapping, n1, n2) {
        var sz = mapSize(mapping);
        return {
            mapping: mapping,
            size: sz,
            tanimoto: tanimoto(sz, n1, n2),
            overlap: simpson(sz, n1, n2)
        };
    }

    // =======================================================================
    // extractSubgraph — extract a subgraph from a SMSDGraph-like object
    // given an array of atom indices. Re-perceives ring membership and
    // aromaticity on the extracted subgraph to prevent aromatic flag leakage.
    // =======================================================================
    function extractSubgraph(g, atomIndices) {
        var subN = atomIndices.length;
        if (subN === 0) { return null; }

        // Build old-to-new index mapping
        var oldToNew = {};
        var k;
        for (k = 0; k < subN; k++) {
            oldToNew[atomIndices[k]] = k;
        }

        // Copy atom properties
        var sub = {
            n: subN,
            atomicNum: new Array(subN),
            formalCharge: new Array(subN),
            aromatic: new Array(subN),
            ring: new Array(subN),
            degree: new Array(subN),
            label: new Array(subN),
            neighbors: new Array(subN),
            bondOrd: {},
            bondRing: {},
            bondArom: {},
            orbit: new Array(subN),
            morganRank: new Array(subN)
        };

        for (k = 0; k < subN; k++) {
            var oldIdx = atomIndices[k];
            sub.atomicNum[k] = g.atomicNum[oldIdx];
            sub.formalCharge[k] = g.formalCharge ? g.formalCharge[oldIdx] : 0;
            sub.aromatic[k] = false; // will re-perceive below
            sub.ring[k] = false;     // will re-perceive below
            sub.label[k] = g.label ? g.label[oldIdx] : sub.atomicNum[k];
            sub.neighbors[k] = [];
            sub.orbit[k] = k;
            sub.morganRank[k] = g.morganRank ? g.morganRank[oldIdx] : 0;
        }

        // Copy bonds where both endpoints are in the subgraph
        for (k = 0; k < subN; k++) {
            var oldI = atomIndices[k];
            var nbs = g.neighbors[oldI];
            for (var ni = 0; ni < nbs.length; ni++) {
                var oldJ = nbs[ni];
                if (oldToNew.hasOwnProperty(oldJ)) {
                    var newJ = oldToNew[oldJ];
                    sub.neighbors[k].push(newJ);
                    var bkey = k + ',' + newJ;
                    sub.bondOrd[bkey] = g.bondOrder(oldI, oldJ);
                    if (g.bondInRing) { sub.bondRing[bkey] = g.bondInRing(oldI, oldJ); }
                    if (g.bondAromatic) { sub.bondArom[bkey] = g.bondAromatic(oldI, oldJ); }
                }
            }
            sub.degree[k] = sub.neighbors[k].length;
        }

        // Attach bondOrder/hasBond/bondInRing/bondAromatic methods
        sub.bondOrder = function(i, j) { return this.bondOrd[i + ',' + j] || 0; };
        sub.hasBond = function(i, j) { return this.bondOrder(i, j) !== 0; };
        sub.bondInRing = function(i, j) { return !!this.bondRing[i + ',' + j]; };
        sub.bondAromatic = function(i, j) { return !!this.bondArom[i + ',' + j]; };

        // Re-perceive ring membership via simple cycle detection (DFS back-edges)
        var inRing = boolArr(subN);
        var ringBonds = {};
        var visited = arrFill(subN, -1);
        var parent = arrFill(subN, -1);
        var stack = [];

        for (var root = 0; root < subN; root++) {
            if (visited[root] >= 0) { continue; }
            stack.push({ node: root, parentNode: -1, nbIdx: 0 });
            visited[root] = 0;
            while (stack.length > 0) {
                var frame = stack[stack.length - 1];
                var u = frame.node;
                var uNbs = sub.neighbors[u];
                if (frame.nbIdx < uNbs.length) {
                    var w = uNbs[frame.nbIdx++];
                    if (visited[w] < 0) {
                        visited[w] = visited[u] + 1;
                        parent[w] = u;
                        stack.push({ node: w, parentNode: u, nbIdx: 0 });
                    } else if (w !== frame.parentNode && visited[w] < visited[u]) {
                        // Back-edge found: mark all atoms on the cycle
                        var cur = u;
                        while (cur !== w) {
                            inRing[cur] = true;
                            var p = parent[cur];
                            var rbk = Math.min(cur, p) + ',' + Math.max(cur, p);
                            ringBonds[rbk] = true;
                            cur = p;
                        }
                        inRing[w] = true;
                        var rbkClose = Math.min(u, w) + ',' + Math.max(u, w);
                        ringBonds[rbkClose] = true;
                    }
                } else {
                    stack.pop();
                }
            }
        }

        for (k = 0; k < subN; k++) { sub.ring[k] = inRing[k]; }
        // Update bond ring flags based on re-perception
        for (var rk in ringBonds) {
            if (ringBonds.hasOwnProperty(rk)) {
                sub.bondRing[rk] = true;
                // Reverse key
                var parts = rk.split(',');
                sub.bondRing[parts[1] + ',' + parts[0]] = true;
            }
        }

        // Re-perceive aromaticity: a ring is aromatic if all atoms in it
        // were aromatic in the original graph. Only set aromatic if the atom
        // remains in a ring in the subgraph.
        for (k = 0; k < subN; k++) {
            if (sub.ring[k] && g.aromatic[atomIndices[k]]) {
                sub.aromatic[k] = true;
            }
        }

        // Kekulisation fallback: atoms that were aromatic in the original
        // graph but are no longer in a complete ring in the subgraph need
        // their aromatic flag cleared and explicit single/double bonds assigned.
        // This prevents invalid aromatic SMILES from partial ring extraction.
        for (k = 0; k < subN; k++) {
            if (g.aromatic[atomIndices[k]] && !sub.ring[k]) {
                sub.aromatic[k] = false;
                // Fix bonds: convert aromatic (order 4/1.5) to explicit kekulised
                var kNbs = sub.neighbors[k];
                for (var ki = 0; ki < kNbs.length; ki++) {
                    var kj = kNbs[ki];
                    var bk1 = k + ',' + kj;
                    var bk2 = kj + ',' + k;
                    // Clear aromatic bond flags
                    sub.bondArom[bk1] = false;
                    sub.bondArom[bk2] = false;
                    // If bond order was inherited as aromatic (0 or unset),
                    // assign based on original graph bond order
                    var origOrd = g.bondOrder(atomIndices[k], atomIndices[kj]);
                    if (!sub.bondOrd[bk1] || sub.bondOrd[bk1] === 0) {
                        var assignOrd = (origOrd === 2) ? 2 : 1;
                        sub.bondOrd[bk1] = assignOrd;
                        sub.bondOrd[bk2] = assignOrd;
                    }
                }
            }
        }

        return sub;
    }

    // =======================================================================
    // mcsToSmiles — extract the MCS as a canonical SMILES string.
    // Requires global.SmilesWriter to be available.
    // =======================================================================
    function mcsToSmiles(g1, mapping) {
        if (!mapping || mapSize(mapping) === 0) { return ''; }
        if (!global.SmilesWriter) { return ''; }

        var atomIndices = mapKeys(mapping);
        atomIndices.sort(function(a, b) { return a - b; });
        var sub = extractSubgraph(g1, atomIndices);
        if (!sub || sub.n === 0) { return ''; }

        // Build a minimal molecule-like object for SmilesWriter
        var atoms = [];
        var bonds = [];
        var SYMBOLS = {};
        // Reverse lookup: atomic number -> symbol
        if (global.SMSDGraph && global.SMSDGraph.ATOMIC_NUMBERS) {
            var AN = global.SMSDGraph.ATOMIC_NUMBERS;
            for (var sym in AN) {
                if (AN.hasOwnProperty(sym)) { SYMBOLS[AN[sym]] = sym; }
            }
        }
        // Common fallbacks
        if (!SYMBOLS[6]) { SYMBOLS[6] = 'C'; }
        if (!SYMBOLS[7]) { SYMBOLS[7] = 'N'; }
        if (!SYMBOLS[8]) { SYMBOLS[8] = 'O'; }
        if (!SYMBOLS[1]) { SYMBOLS[1] = 'H'; }
        if (!SYMBOLS[16]) { SYMBOLS[16] = 'S'; }
        if (!SYMBOLS[15]) { SYMBOLS[15] = 'P'; }
        if (!SYMBOLS[9]) { SYMBOLS[9] = 'F'; }
        if (!SYMBOLS[17]) { SYMBOLS[17] = 'Cl'; }
        if (!SYMBOLS[35]) { SYMBOLS[35] = 'Br'; }
        if (!SYMBOLS[53]) { SYMBOLS[53] = 'I'; }

        for (var i = 0; i < sub.n; i++) {
            atoms.push({
                id: i,
                symbol: SYMBOLS[sub.atomicNum[i]] || ('*'),
                charge: sub.formalCharge[i] || 0,
                aromatic: sub.aromatic[i]
            });
        }

        // Build bonds from neighbors (only add each bond once)
        var addedBonds = {};
        for (var ai = 0; ai < sub.n; ai++) {
            var nbsSub = sub.neighbors[ai];
            for (var ni2 = 0; ni2 < nbsSub.length; ni2++) {
                var aj = nbsSub[ni2];
                if (aj <= ai) { continue; }
                var bondKey = ai + ',' + aj;
                if (addedBonds[bondKey]) { continue; }
                addedBonds[bondKey] = true;
                bonds.push({
                    atom1: ai,
                    atom2: aj,
                    type: sub.bondOrder(ai, aj) || 1
                });
            }
        }

        // Use a real Molecule object so SmilesWriter has all required methods
        if (typeof global.Molecule !== 'function') { return ''; }
        var mol = new global.Molecule();
        var idMap = {};
        for (var i = 0; i < atoms.length; i++) {
            var a = atoms[i];
            var newAtom = mol.addAtom(a.symbol, 0, 0);
            idMap[i] = newAtom.id;
            newAtom.charge = a.charge;
            newAtom.aromatic = a.aromatic;
        }
        for (var bi = 0; bi < bonds.length; bi++) {
            var b = bonds[bi];
            mol.addBond(idMap[b.atom1], idMap[b.atom2], b.type);
        }
        try {
            return global.SmilesWriter.write(mol);
        } catch (e) {
            return '';
        }
    }

    // =======================================================================
    // findAllMCS — return top-N MCS mappings using multi-seed perturbation.
    // Each seed starts from a different initial atom pair or shuffled order,
    // then extends via the full findMCS pipeline. Duplicate mappings (by
    // sorted query-atom set) are filtered out.
    // =======================================================================
    function findAllMCS(g1, g2, chemOpts, mcsOpts, topN) {
        topN = topN || 5;
        var C = chemOpts || { matchAtomType: true, matchBondType: true };
        var M = mergeOpts(McsOptions, mcsOpts);
        var perSeedTimeout = Math.max(100, Math.floor(M.timeoutMs / (topN + 2)));

        // Canonical SMILES dedup key: use mcsToSmiles for dedup when available,
        // fall back to sorted query-atom set
        function dedupKey(mapping) {
            var smi = mcsToSmiles(g1, mapping);
            if (smi && smi.length > 0) { return 'smi:' + smi; }
            return 'idx:' + mapKeys(mapping).sort(function(a, b) { return a - b; }).join(',');
        }

        // First, get the best MCS via normal findMCS
        var results = [];
        var seenKeys = {};
        var best = findMCS(g1, g2, C, M);
        if (best && best.size > 0) {
            var key0 = dedupKey(best.mapping);
            results.push(best);
            seenKeys[key0] = true;
        }

        if (results.length >= topN || g1.n === 0 || g2.n === 0) { return results; }

        // Generate seed perturbations by choosing different anchor atoms
        var n1 = g1.n, n2 = g2.n;
        var minN = Math.min(n1, n2);
        var seedPairs = [];
        var si, sj;

        // Collect compatible pairs sorted by degree product (diverse seeds)
        for (si = 0; si < n1 && seedPairs.length < topN * 4; si++) {
            for (sj = 0; sj < n2 && seedPairs.length < topN * 4; sj++) {
                if (atomsCompat(g1, si, g2, sj, C)) {
                    seedPairs.push({
                        qi: si, tj: sj,
                        score: g1.degree[si] + g2.degree[sj]
                            + (g1.ring[si] && g2.ring[sj] ? 20 : 0)
                    });
                }
            }
        }
        // Sort by score descending, then take diverse subset
        seedPairs.sort(function(a, b) { return b.score - a.score; });

        // Try each seed pair with a shorter timeout
        for (var sp = 0; sp < seedPairs.length && results.length < topN; sp++) {
            var pair = seedPairs[sp];
            var seedMap = {};
            seedMap[pair.qi] = pair.tj;

            // Extend the seed via greedy + McGregor
            var extMopts = mergeOpts(M, { timeoutMs: perSeedTimeout });
            var tb2 = new TimeBudget(perSeedTimeout);
            var extended = greedyAtomExtend(g1, g2, seedMap, C, extMopts);
            extended = ppx(g1, g2, extended, C, extMopts);

            // Try McGregor if the seed produced something
            if (mapSize(extended) > 0 && !tb2.expired()) {
                var mcgResult = mcGregorExtend(g1, g2, extended, C, tb2,
                    Math.max(1, perSeedTimeout / 2),
                    extMopts.useTwoHopNLF, extMopts.useThreeHopNLF, extMopts.connectedOnly);
                mcgResult = ppx(g1, g2, mcgResult, C, extMopts);
                if (mapSize(mcgResult) > mapSize(extended)) { extended = mcgResult; }
            }

            var extSize = mapSize(extended);
            if (extSize === 0) { continue; }
            var extKey = dedupKey(extended);
            if (seenKeys[extKey]) { continue; }
            seenKeys[extKey] = true;
            results.push(formatResult(extended, n1, n2));
        }

        // Sort results by size descending
        results.sort(function(a, b) { return b.size - a.size; });
        if (results.length > topN) { results.length = topN; }
        return results;
    }

    // =======================================================================
    // mcs_size — lightweight function returning just the MCS size without
    // computing the full mapping. Uses the label-frequency and degree-sequence
    // upper bounds, then a fast greedy probe as a lower bound. If the bounds
    // match, returns immediately. Otherwise falls back to the full algorithm
    // but with a reduced timeout.
    // =======================================================================
    function mcs_size(g1, g2, chemOpts, mcsOpts) {
        if (g1.n === 0 || g2.n === 0) { return 0; }
        var C = chemOpts || { matchAtomType: true, matchBondType: true };
        var M = mergeOpts(McsOptions, mcsOpts);

        // Quick upper bound
        var ub = Math.min(
            labelFrequencyUpperBound(g1, g2, C),
            degreeSequenceUpperBound(g1, g2, C)
        );
        if (ub === 0) { return 0; }

        // Quick lower bound via greedy probe
        var greedy = greedyProbe(g1, g2, C, M.templateFuzzyAtoms);
        var lb = mapSize(greedy);

        // If bounds match, we know the exact answer
        if (lb >= ub) { return ub; }

        // Single-atom trivial check
        if (g1.n === 1 || g2.n === 1) {
            if (lb > 0) { return lb; }
            for (var i = 0; i < g1.n; i++) {
                for (var j = 0; j < g2.n; j++) {
                    if (atomsCompat(g1, i, g2, j, C)) { return 1; }
                }
            }
            return 0;
        }

        // Fall back to full MCS with reduced timeout for screening
        var screenTimeout = Math.min(M.timeoutMs, 2000);
        var result = findMCS(g1, g2, C, mergeOpts(M, { timeoutMs: screenTimeout }));
        return result ? result.size : lb;
    }

    // =======================================================================
    // findMcsSmarts — find the largest MCS substructure that matches a
    // given SMARTS pattern. Runs findAllMCS and filters results through
    // a SMARTS substructure check, returning the largest match.
    // Requires global.SmartsParser and global.SMSDVF2 to be available.
    // =======================================================================
    function findMcsSmarts(g1, g2, smartsPattern, chemOpts, mcsOpts, topN) {
        topN = topN || 10;
        if (!smartsPattern || !global.SmartsParser || !global.SMSDVF2) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
        }

        var allMCS = findAllMCS(g1, g2, chemOpts, mcsOpts, topN);
        if (!allMCS || allMCS.length === 0) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
        }

        var queryGraph;
        try {
            queryGraph = global.SmartsParser.parse(smartsPattern);
        } catch (e) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
        }
        if (!queryGraph || queryGraph.n === 0) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
        }

        // Check each MCS result: extract subgraph, test SMARTS containment
        for (var i = 0; i < allMCS.length; i++) {
            var mcs = allMCS[i];
            if (!mcs || mcs.size === 0) { continue; }
            var atomIndices = mapKeys(mcs.mapping);
            atomIndices.sort(function(a, b) { return a - b; });
            var sub = extractSubgraph(g1, atomIndices);
            if (!sub || sub.n === 0) { continue; }

            // Check if the SMARTS pattern is a substructure of the MCS subgraph
            var subMatch = global.SMSDVF2.findSubstructure(
                queryGraph, sub, chemOpts, 2000
            );
            if (subMatch && mapSize(subMatch) > 0) {
                return mcs;
            }
        }

        return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0 };
    }

    // -----------------------------------------------------------------------
    // translateToAtomIds — convert SMSDGraph indices to Molecule atom IDs
    // -----------------------------------------------------------------------
    function translateToAtomIds(mapping, g1, g2) {
        var result = {};
        var k;
        for (k in mapping) {
            if (mapping.hasOwnProperty(k)) {
                var idx1 = +k;
                var idx2 = mapping[k];
                // Bounds checking
                if (idx1 < 0 || idx1 >= g1.n) { continue; }
                if (idx2 < 0 || idx2 >= g2.n) { continue; }
                var id1 = g1.idxToId ? g1.idxToId[idx1] : idx1;
                var id2 = g2.idxToId ? g2.idxToId[idx2] : idx2;
                result[id1] = id2;
            }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // findMCSFromSmiles — convenience: parse SMILES, build graphs, run MCS
    // -----------------------------------------------------------------------
    function findMCSFromSmiles(smi1, smi2, chemOpts, mcsOpts) {
        if (!global.SmilesParser || !global.SmilesParser.parse) {
            throw new Error('SmilesParser not available');
        }
        if (!global.SMSDGraph || !global.SMSDGraph.SMSDGraph) {
            throw new Error('SMSDGraph not available');
        }

        var mol1 = global.SmilesParser.parse(smi1);
        if (!mol1 || (mol1.parseErrors && mol1.parseErrors.length > 0)) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0, mcsSmiles: '', error: 'Failed to parse first SMILES' };
        }
        var mol2 = global.SmilesParser.parse(smi2);
        if (!mol2 || (mol2.parseErrors && mol2.parseErrors.length > 0)) {
            return { mapping: {}, size: 0, tanimoto: 0.0, overlap: 0.0, mcsSmiles: '', error: 'Failed to parse second SMILES' };
        }

        var SG = global.SMSDGraph.SMSDGraph;
        var g1 = new SG(mol1);
        var g2 = new SG(mol2);

        var result = findMCS(g1, g2, chemOpts, mcsOpts);
        var smilesOut = mcsToSmiles(g1, result.mapping);

        return {
            mapping: result.mapping,
            size: result.size,
            tanimoto: result.tanimoto,
            overlap: result.overlap,
            mcsSmiles: smilesOut
        };
    }

    // =======================================================================
    // ChemOptions.copyOf — safe deep copy of chemistry matching options
    // =======================================================================

    /**
     * Create an independent deep copy of a ChemOptions-like object.
     * Prevents mutation of shared option objects across concurrent searches.
     *
     * @param {Object} opts — chemistry options (matchAtomType, matchBondType, etc.)
     * @returns {Object} — deep copy with all nested values duplicated
     */
    function chemOptsCopyOf(opts) {
        if (!opts) { return { matchAtomType: true, matchBondType: true }; }
        var copy = {};
        var k;
        for (k in opts) {
            if (opts.hasOwnProperty(k)) {
                var v = opts[k];
                if (v && typeof v === 'object' && !Array.isArray(v)) {
                    // Shallow clone nested objects (one level deep is sufficient)
                    var inner = {};
                    for (var ik in v) {
                        if (v.hasOwnProperty(ik)) { inner[ik] = v[ik]; }
                    }
                    copy[k] = inner;
                } else if (Array.isArray(v)) {
                    copy[k] = v.slice();
                } else {
                    copy[k] = v;
                }
            }
        }
        return copy;
    }

    // =======================================================================
    // Charge relaxation for reaction centres (v6.4)
    // =======================================================================

    /**
     * Relax formal charge matching at reaction centre atoms.
     * In reaction-aware MCS, formal charge differences between reactant and
     * product sides at mapped centres should not prevent matching, since
     * charge transfer is a normal outcome of bond-breaking/formation.
     *
     * Returns a modified ChemOptions copy with matchFormalCharge = false
     * when reactionAware is enabled.
     *
     * @param {Object} chemOpts — original chemistry options
     * @param {Object} mcsOpts  — MCS options (checked for reactionAware flag)
     * @returns {Object} — potentially relaxed chemistry options (safe copy)
     */
    function relaxChargeForReaction(chemOpts, mcsOpts) {
        var C = chemOptsCopyOf(chemOpts);
        if (mcsOpts && mcsOpts.reactionAware) {
            C.matchFormalCharge = false;
        }
        return C;
    }

    // =======================================================================
    // Composite MCS scoring (v6.4)
    // Weights: size(0.25) + hetero(0.40) + rarity(0.25) + conn(0.10)
    // =======================================================================

    /**
     * Element rarity score: how uncommon are the matched atoms?
     * Common elements (C, N, O, S) get low rarity; less common get higher.
     */
    var RARITY_TABLE = { 6: 0.1, 7: 0.3, 8: 0.25, 16: 0.35, 15: 0.6,
        9: 0.5, 17: 0.55, 35: 0.65, 53: 0.7, 5: 0.7, 14: 0.6, 34: 0.7 };

    function atomRarity(atomicNum) {
        return RARITY_TABLE[atomicNum] || 0.8;
    }

    /**
     * Connectivity score: ratio of bonds preserved in the MCS mapping
     * relative to total bonds in the smaller graph.
     */
    function connectivityScore(g1, g2, mapping) {
        var bondCount = 0;
        var totalBonds = 0;
        var k, k2;
        for (k in mapping) {
            if (!mapping.hasOwnProperty(k)) { continue; }
            var qi = +k;
            var tj = mapping[k];
            var qNbs = g1.neighbors[qi];
            for (var ni = 0; ni < qNbs.length; ni++) {
                var qk = qNbs[ni];
                if (qk > qi && mapping.hasOwnProperty(qk)) {
                    var tk = mapping[qk];
                    if (g2.hasBond(tj, tk)) { bondCount++; }
                }
                if (qk > qi) { totalBonds++; }
            }
        }
        return totalBonds > 0 ? bondCount / totalBonds : 0;
    }

    /**
     * Composite score for ranking MCS candidates.
     * Weights: size(0.25) + hetero(0.40) + rarity(0.25) + conn(0.10)
     *
     * @param {Object} g1      — query graph
     * @param {Object} g2      — target graph
     * @param {Object} mapping — MCS mapping { queryIdx: targetIdx }
     * @returns {number} — composite score in [0, 1]
     */
    function compositeScore(g1, g2, mapping) {
        var sz = mapSize(mapping);
        if (sz === 0) { return 0; }
        var maxPossible = Math.min(g1.n, g2.n);
        var sizeNorm = sz / maxPossible;

        // Heteroatom fraction
        var heteroCount = countHeteroatomsInMapping(g1, mapping);
        var heteroNorm = sz > 0 ? heteroCount / sz : 0;

        // Rarity: average rarity of matched atoms
        var raritySum = 0;
        var k;
        for (k in mapping) {
            if (mapping.hasOwnProperty(k)) {
                raritySum += atomRarity(g1.atomicNum[+k]);
            }
        }
        var rarityNorm = sz > 0 ? raritySum / sz : 0;

        // Connectivity
        var connNorm = connectivityScore(g1, g2, mapping);

        return 0.25 * sizeNorm + 0.40 * heteroNorm + 0.25 * rarityNorm + 0.10 * connNorm;
    }

    // =======================================================================
    // findNearMCS — K-1 and K-2 near-optimal candidates (v6.4)
    // =======================================================================

    /**
     * Generate near-MCS candidates by removing 1 or 2 atoms (K-1, K-2)
     * from the best MCS mapping. Returns an array of candidate mappings
     * ranked by composite score.
     *
     * Useful for reaction-aware analysis where the optimal MCS may miss
     * reaction centre atoms due to charge/bond-order changes.
     *
     * @param {Object} g1       — query graph
     * @param {Object} g2       — target graph
     * @param {Object} chemOpts — chemistry matching options
     * @param {Object} mcsOpts  — MCS options (reactionAware checked)
     * @param {number} [topN]   — max candidates to return (default 10)
     * @returns {Array<Object>} — ranked candidates with mapping, size, score
     */
    function findNearMCS(g1, g2, chemOpts, mcsOpts, topN) {
        topN = topN || 10;
        var M = mergeOpts(McsOptions, mcsOpts);
        var C = M.reactionAware ? relaxChargeForReaction(chemOpts, M) : (chemOpts || { matchAtomType: true, matchBondType: true });

        // Get all MCS candidates as starting points
        var allMCS = findAllMCS(g1, g2, C, M, Math.min(topN, 5));
        if (!allMCS || allMCS.length === 0) { return []; }

        var candidates = [];
        var seenKeys = {};
        var ci, k, k2;

        // Helper: canonical key for dedup
        function candKey(m) {
            var keys = mapKeys(m).sort(function(a, b) { return a - b; });
            return keys.join(',');
        }

        // Add base MCS results
        for (ci = 0; ci < allMCS.length; ci++) {
            var base = allMCS[ci].mapping;
            var bKey = candKey(base);
            if (!seenKeys[bKey]) {
                seenKeys[bKey] = true;
                candidates.push({
                    mapping: base,
                    size: mapSize(base),
                    score: compositeScore(g1, g2, base),
                    level: 'K-0'
                });
            }
        }

        // Generate K-1 candidates (remove one atom)
        for (ci = 0; ci < allMCS.length; ci++) {
            var baseMap = allMCS[ci].mapping;
            var baseKeys = mapKeys(baseMap);
            for (var ri = 0; ri < baseKeys.length; ri++) {
                var reduced = mapCopy(baseMap);
                delete reduced[baseKeys[ri]];
                var rKey = candKey(reduced);
                if (seenKeys[rKey]) { continue; }
                seenKeys[rKey] = true;

                // Verify remaining mapping is still valid
                var valid = true;
                for (k in reduced) {
                    if (!reduced.hasOwnProperty(k)) { continue; }
                    if (!atomsCompat(g1, +k, g2, reduced[k], C)) { valid = false; break; }
                }
                if (!valid) { continue; }

                candidates.push({
                    mapping: reduced,
                    size: mapSize(reduced),
                    score: compositeScore(g1, g2, reduced),
                    level: 'K-1'
                });
            }
        }

        // Generate K-2 candidates (remove two atoms) — only from best MCS
        if (allMCS.length > 0) {
            var bestBase = allMCS[0].mapping;
            var bestKeys = mapKeys(bestBase);
            if (bestKeys.length > 4) { // Only worth it for larger MCS
                for (var ri1 = 0; ri1 < bestKeys.length && candidates.length < topN * 3; ri1++) {
                    for (var ri2 = ri1 + 1; ri2 < bestKeys.length && candidates.length < topN * 3; ri2++) {
                        var reduced2 = mapCopy(bestBase);
                        delete reduced2[bestKeys[ri1]];
                        delete reduced2[bestKeys[ri2]];
                        var rKey2 = candKey(reduced2);
                        if (seenKeys[rKey2]) { continue; }
                        seenKeys[rKey2] = true;

                        var valid2 = true;
                        for (k2 in reduced2) {
                            if (!reduced2.hasOwnProperty(k2)) { continue; }
                            if (!atomsCompat(g1, +k2, g2, reduced2[k2], C)) { valid2 = false; break; }
                        }
                        if (!valid2) { continue; }

                        candidates.push({
                            mapping: reduced2,
                            size: mapSize(reduced2),
                            score: compositeScore(g1, g2, reduced2),
                            level: 'K-2'
                        });
                    }
                }
            }
        }

        // Sort by composite score descending
        candidates.sort(function(a, b) { return b.score - a.score; });
        if (candidates.length > topN) { candidates.length = topN; }
        return candidates;
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------
    var SMSDMCS = {};

    SMSDMCS.McsOptions = McsOptions;

    SMSDMCS.findMCS = findMCS;
    SMSDMCS.findMCSFromSmiles = findMCSFromSmiles;
    SMSDMCS.findAllMCS = findAllMCS;
    SMSDMCS.findNearMCS = findNearMCS;
    SMSDMCS.findMcsSmarts = findMcsSmarts;
    SMSDMCS.mcs_size = mcs_size;
    SMSDMCS.extractSubgraph = extractSubgraph;
    SMSDMCS.mcsToSmiles = mcsToSmiles;
    SMSDMCS.translateToAtomIds = translateToAtomIds;

    SMSDMCS.compositeScore = compositeScore;
    SMSDMCS.chemOptsCopyOf = chemOptsCopyOf;
    SMSDMCS.similarityUpperBound = similarityUpperBound;

    // Expose internals for testing
    SMSDMCS._labelFrequencyUpperBound = labelFrequencyUpperBound;
    SMSDMCS._degreeSequenceUpperBound = degreeSequenceUpperBound;
    SMSDMCS._greedyProbe = greedyProbe;
    SMSDMCS._largestConnected = largestConnected;
    SMSDMCS._linearChainFastPath = linearChainFastPath;
    SMSDMCS._treeFastPath = treeFastPath;
    SMSDMCS._relaxChargeForReaction = relaxChargeForReaction;
    SMSDMCS._connectivityScore = connectivityScore;
    SMSDMCS._atomRarity = atomRarity;

    global.SMSDMCS = SMSDMCS;

})(typeof window !== 'undefined' ? window : this);
