/**
 * SMSDVF2.js — VF2++ Substructure Search Engine
 *
 * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * Ported from SMSD upstream, vf2pp.hpp
 *
 * Algorithms:
 *   - Greedy probe (O(N) fast-path)
 *   - FASTiso ordering (forward-checking, smallest domain first, < 30 atoms)
 *   - VF3-Light ordering (label rarity sort, >= 30 atoms)
 *   - VF2++ backtracking with NLF-1/2 pruning
 *   - Full feasibility check: degree, NLF, mapped-neighbor adjacency, bond compat
 */
(function() {
    'use strict';

    var SG = window.SMSDGraph;
    var atomsCompat = SG.atomsCompat;
    var bondsCompat = SG.bondsCompat;
    var nlfOk       = SG.nlfOk;

    // ========================================================================
    // TimeBudget — check Date.now() every 1024 iterations
    // ========================================================================

    function TimeBudget(ms) {
        this.deadline = Date.now() + Math.max(1, ms);
        this.counter = 0;
    }

    TimeBudget.prototype.expired = function() {
        this.counter++;
        if ((this.counter & 1023) !== 0) return false;
        return Date.now() >= this.deadline;
    };

    // ========================================================================
    // popcount — count set bits in a 32-bit integer
    // ========================================================================

    function popcount32(x) {
        x = x - ((x >>> 1) & 0x55555555);
        x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
        return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24;
    }

    // ES5-safe count trailing zeros using de Bruijn sequence
    var DE_BRUIJN = [
        0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
        31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    ];
    function ctz32(v) {
        // v must be non-zero and a power of 2 (already isolated lsb)
        return DE_BRUIJN[((v * 0x077CB531) >>> 27) & 31];
    }

    // ========================================================================
    // Bit-set helpers for candidate domains (array of 32-bit integers)
    // We use 32-bit words since JS bitwise ops are 32-bit
    // ========================================================================

    function wordCount(n) {
        return ((n + 31) >>> 5);
    }

    function setBit(words, idx) {
        words[idx >>> 5] |= (1 << (idx & 31));
    }

    function clearBit(words, idx) {
        words[idx >>> 5] &= ~(1 << (idx & 31));
    }

    function testBit(words, idx) {
        return (words[idx >>> 5] & (1 << (idx & 31))) !== 0;
    }

    function iterateBits(words, nWords, out) {
        var count = 0;
        for (var w = 0; w < nWords; w++) {
            var bits = words[w];
            while (bits !== 0) {
                var lsb = bits & (-bits);
                var bitIdx = ctz32(lsb);
                out[count++] = (w << 5) | bitIdx;
                bits ^= lsb;
            }
        }
        return count;
    }

    function popcountWords(words, nWords) {
        var c = 0;
        for (var w = 0; w < nWords; w++) c += popcount32(words[w] | 0);
        return c;
    }

    // ========================================================================
    // VF2PPMatcher — main backtracking engine
    // ========================================================================

    function VF2PPMatcher(gq, gt, opts, timeoutMs) {
        this.gq = gq;
        this.gt = gt;
        this.opts = opts;
        this.tb = new TimeBudget(timeoutMs || 10000);

        this.Nq = gq.n;
        this.Nt = gt.n;
        this.tWords = wordCount(gt.n);

        // Mapping arrays
        this.q2t = new Array(gq.n);
        this.t2q = new Array(gt.n);
        var i, j;
        for (i = 0; i < gq.n; i++) this.q2t[i] = -1;
        for (i = 0; i < gt.n; i++) this.t2q[i] = -1;

        // Used mask (bit-set of matched target atoms)
        this.usedMask = new Array(this.tWords);
        for (i = 0; i < this.tWords; i++) this.usedMask[i] = 0;

        // NLF caches
        this.qNLF1 = gq.getNLF1();
        this.tNLF1 = gt.getNLF1();
        this.useTwoHop = opts.useTwoHopNLF !== false;
        // Adaptive: disable 2-hop for small molecules
        if (this.Nq <= 12 || this.Nt <= 12) this.useTwoHop = false;
        this.qNLF2 = null;
        this.tNLF2 = null;
        if (this.useTwoHop) {
            this.qNLF2 = gq.getNLF2();
            this.tNLF2 = gt.getNLF2();
        }

        // Ensure tautomer classes if needed
        if (opts.tautomerAware) {
            gq.ensureTautomers(opts.pH);
            gt.ensureTautomers(opts.pH);
        }

        // Build query neighbor lists sorted by descending degree
        this.qNeighborsByDegDesc = new Array(gq.n);
        for (i = 0; i < gq.n; i++) {
            var nbCopy = gq.neighbors[i].slice();
            // Insertion sort descending by degree
            for (var a = 1; a < nbCopy.length; a++) {
                var key = nbCopy[a], keyDeg = gq.degree[key];
                var b = a - 1;
                while (b >= 0 && gq.degree[nbCopy[b]] < keyDeg) {
                    nbCopy[b + 1] = nbCopy[b];
                    b--;
                }
                nbCopy[b + 1] = key;
            }
            this.qNeighborsByDegDesc[i] = nbCopy;
        }

        // Target adjacency sets (for fast adjacency check)
        // adjSet[i] is an object where adjSet[i][j] = true if i-j are bonded
        this.adjSet = new Array(gt.n);
        for (i = 0; i < gt.n; i++) {
            var s = {};
            var gtNb = gt.neighbors[i];
            for (j = 0; j < gtNb.length; j++) s[gtNb[j]] = true;
            this.adjSet[i] = s;
        }

        // Build candidate domains
        this.domain = new Array(gq.n);
        for (i = 0; i < gq.n; i++) {
            var dom = new Array(this.tWords);
            for (j = 0; j < this.tWords; j++) dom[j] = 0;
            for (j = 0; j < gt.n; j++) {
                if (atomsCompat(gq, i, gt, j, opts)) {
                    setBit(dom, j);
                }
            }
            this.domain[i] = dom;
        }

        // Frontier mask for VF2++ candidate selection
        this.frontierMask = new Array(this.tWords);
        for (i = 0; i < this.tWords; i++) this.frontierMask[i] = 0;
        this.frontierStack = new Array(gq.n);
        for (i = 0; i < gq.n; i++) {
            this.frontierStack[i] = new Array(this.tWords);
        }

        // Candidate buffers (per-depth)
        this.candBuf = new Array(Math.max(1, gq.n));
        for (i = 0; i < Math.max(1, gq.n); i++) {
            this.candBuf[i] = new Array(gt.n + 64);
        }
        this.probeCandBuf = new Array(gt.n + 64);

        // Stats
        this.found = false;
        this.timedOut = false;
        this.nodesVisited = 0;
    }

    // ========================================================================
    // Domain candidates: available & compatible target atoms for query atom qi
    // ========================================================================

    VF2PPMatcher.prototype.domainCandidates = function(qi, out) {
        var count = 0;
        var tw = this.tWords;
        var dom = this.domain[qi];
        var used = this.usedMask;
        for (var w = 0; w < tw; w++) {
            var bits = dom[w] & ~used[w];
            while (bits !== 0) {
                var lsb = bits & (-bits);
                var bitIdx = ctz32(lsb);
                out[count++] = (w << 5) | bitIdx;
                bits ^= lsb;
            }
        }
        return count;
    };

    // ========================================================================
    // Feasibility check
    // ========================================================================

    VF2PPMatcher.prototype.feasible = function(qi, tj) {
        // Degree check
        if (this.gt.degree[tj] < this.gq.degree[qi]) return false;

        // Ring-only check
        if (this.opts.ringMatchesRingOnly && this.gq.ring[qi] && !this.gt.ring[tj]) return false;

        // NLF-1 check
        if (!nlfOk(this.qNLF1[qi], this.tNLF1[tj])) return false;

        // NLF-2 check
        if (this.useTwoHop && !nlfOk(this.qNLF2[qi], this.tNLF2[tj])) return false;

        // Mapped-neighbor adjacency and bond compatibility
        var nbSorted = this.qNeighborsByDegDesc[qi];
        var adjSetTj = this.adjSet[tj];
        for (var k = 0; k < nbSorted.length; k++) {
            var qk = nbSorted[k];
            var tk = this.q2t[qk];
            if (tk !== -1) {
                // Must be adjacent in target
                if (!adjSetTj[tk]) return false;
                // Bond compatibility
                if (!bondsCompat(this.gq, qi, qk, this.gt, tj, tk, this.opts)) return false;
            }
        }

        return true;
    };

    // ========================================================================
    // Candidate selection (VF2++ with frontier mask)
    // ========================================================================

    VF2PPMatcher.prototype.selectCandidates = function(qi, out) {
        var tw = this.tWords;
        var dom = this.domain[qi];
        var used = this.usedMask;
        var fm = this.frontierMask;
        var w, bits, count, lsb, bitIdx;

        // Check if frontier is non-empty AND qi is connected to mapped subgraph
        var hasFrontier = false;
        for (w = 0; w < tw; w++) {
            if (fm[w] !== 0) { hasFrontier = true; break; }
        }
        if (hasFrontier) {
            var qiConnected = false;
            var qNb = this.gq.neighbors[qi];
            for (var k = 0; k < qNb.length; k++) {
                if (this.q2t[qNb[k]] !== -1) { qiConnected = true; break; }
            }
            if (!qiConnected) hasFrontier = false;
        }

        if (hasFrontier) {
            var hasAny = false;
            for (w = 0; w < tw; w++) {
                if ((dom[w] & ~used[w] & fm[w]) !== 0) { hasAny = true; break; }
            }
            if (hasAny) {
                count = 0;
                for (w = 0; w < tw; w++) {
                    bits = dom[w] & ~used[w] & fm[w];
                    while (bits !== 0) {
                        lsb = bits & (-bits);
                        bitIdx = ctz32(lsb);
                        out[count++] = (w << 5) | bitIdx;
                        bits ^= lsb;
                    }
                }
                this.sortByDegreeProximity(out, count, qi);
                return count;
            }
        }

        // Fallback: all domain candidates
        count = 0;
        for (w = 0; w < tw; w++) {
            bits = dom[w] & ~used[w];
            while (bits !== 0) {
                lsb = bits & (-bits);
                bitIdx = ctz32(lsb);
                out[count++] = (w << 5) | bitIdx;
                bits ^= lsb;
            }
        }
        this.sortByDegreeProximity(out, count, qi);
        return count;
    };

    // ========================================================================
    // Sort candidates by degree proximity to query atom
    // ========================================================================

    VF2PPMatcher.prototype.sortByDegreeProximity = function(cands, n, qi) {
        var qd = this.gq.degree[qi];
        var tdeg = this.gt.degree;
        for (var i = 1; i < n; i++) {
            var key = cands[i];
            var keyDist = tdeg[key] - qd;
            if (keyDist < 0) keyDist = -keyDist;
            var j = i - 1;
            while (j >= 0) {
                var d = tdeg[cands[j]] - qd;
                if (d < 0) d = -d;
                if (d <= keyDist) break;
                cands[j + 1] = cands[j];
                j--;
            }
            cands[j + 1] = key;
        }
    };

    // ========================================================================
    // Frontier management
    // ========================================================================

    VF2PPMatcher.prototype.onMatch = function(pos, tj) {
        var tw = this.tWords;
        var fm = this.frontierMask;
        var stack = this.frontierStack[pos];
        var w;
        // Save frontier state
        for (w = 0; w < tw; w++) stack[w] = fm[w];
        // Add unmatched neighbors of tj to frontier
        var gtNb = this.gt.neighbors[tj];
        for (var k = 0; k < gtNb.length; k++) {
            var u = gtNb[k];
            if (this.t2q[u] === -1) setBit(fm, u);
        }
        // Remove tj itself from frontier
        clearBit(fm, tj);
    };

    VF2PPMatcher.prototype.onUnmatch = function(pos) {
        var tw = this.tWords;
        var fm = this.frontierMask;
        var stack = this.frontierStack[pos];
        for (var w = 0; w < tw; w++) fm[w] = stack[w];
    };

    // ========================================================================
    // Greedy probe — O(N) fast-path
    // ========================================================================

    VF2PPMatcher.prototype.greedyProbe = function(order) {
        var matched = 0;
        var i, qi, nCands, bestTj, bestDist, c, tj, dist;
        for (i = 0; i < this.Nq; i++) {
            // Timeout check — prevents greedy probe from hanging on
            // complex ring systems (e.g., morphinan 5-ring bridged)
            if (this.tb && this.tb.expired()) return false;
            qi = order[i];
            nCands = this.domainCandidates(qi, this.probeCandBuf);
            bestTj = -1;
            bestDist = 0x7FFFFFFF;
            for (c = 0; c < nCands; c++) {
                tj = this.probeCandBuf[c];
                dist = this.gt.degree[tj] - this.gq.degree[qi];
                if (dist < 0) dist = -dist;
                if (dist < bestDist && this.feasible(qi, tj)) {
                    bestTj = tj;
                    bestDist = dist;
                    if (dist === 0) break;
                }
            }
            if (bestTj === -1) {
                // Undo all matched so far
                for (var k = matched - 1; k >= 0; k--) {
                    var qk = order[k], tk = this.q2t[qk];
                    this.q2t[qk] = -1;
                    this.t2q[tk] = -1;
                    clearBit(this.usedMask, tk);
                }
                return false;
            }
            this.q2t[qi] = bestTj;
            this.t2q[bestTj] = qi;
            setBit(this.usedMask, bestTj);
            matched++;
        }
        return true;
    };

    // ========================================================================
    // FASTiso ordering — forward-checking, smallest domain first
    // ========================================================================

    VF2PPMatcher.prototype.fastisoOrder = function(order) {
        var Nq = this.Nq;
        var tw = this.tWords;
        var picked = new Array(Nq);
        var i, pos, bestQ, bestSize, bestDeg, size, w;
        for (i = 0; i < Nq; i++) picked[i] = false;

        // Connected components via BFS
        var compId = new Array(Nq);
        for (i = 0; i < Nq; i++) compId[i] = -1;
        var nComp = 0;
        for (var s = 0; s < Nq; s++) {
            if (compId[s] >= 0) continue;
            var cid = nComp++;
            compId[s] = cid;
            var bfs = [s], head = 0;
            while (head < bfs.length) {
                var u = bfs[head++];
                var nb = this.gq.neighbors[u];
                for (var k = 0; k < nb.length; k++) {
                    if (compId[nb[k]] < 0) {
                        compId[nb[k]] = cid;
                        bfs.push(nb[k]);
                    }
                }
            }
        }

        // Working copy of domains for forward checking
        var workDomain = new Array(Nq);
        for (i = 0; i < Nq; i++) {
            workDomain[i] = this.domain[i].slice();
        }

        var reachable = new Array(tw);
        var activeComp = -1;

        for (pos = 0; pos < Nq; pos++) {
            bestQ = -1; bestSize = 0x7FFFFFFF; bestDeg = -1;

            // Prefer atoms from active component
            if (activeComp >= 0) {
                for (i = 0; i < Nq; i++) {
                    if (picked[i] || compId[i] !== activeComp) continue;
                    size = popcountWords(workDomain[i], tw);
                    if (size < bestSize || (size === bestSize && this.gq.degree[i] > bestDeg)) {
                        bestQ = i; bestSize = size; bestDeg = this.gq.degree[i];
                    }
                }
            }

            // Global fallback
            if (bestQ < 0) {
                for (i = 0; i < Nq; i++) {
                    if (picked[i]) continue;
                    size = popcountWords(workDomain[i], tw);
                    if (size < bestSize || (size === bestSize && this.gq.degree[i] > bestDeg)) {
                        bestQ = i; bestSize = size; bestDeg = this.gq.degree[i];
                    }
                }
                activeComp = compId[bestQ];
            }

            order[pos] = bestQ;
            picked[bestQ] = true;

            // Forward-check: restrict neighbor domains to reachable from bestQ's domain
            var bestNb = this.gq.neighbors[bestQ];
            for (var nk = 0; nk < bestNb.length; nk++) {
                var qk2 = bestNb[nk];
                if (picked[qk2]) continue;
                // Check if qk2's domain is already empty
                var hasAny = false;
                for (w = 0; w < tw && !hasAny; w++) {
                    if (workDomain[qk2][w]) hasAny = true;
                }
                if (!hasAny) continue;
                // Build reachable mask: union of neighbors of all target atoms in bestQ's domain
                for (w = 0; w < tw; w++) reachable[w] = 0;
                for (w = 0; w < tw; w++) {
                    var bits = workDomain[bestQ][w];
                    while (bits !== 0) {
                        var lsb = bits & (-bits);
                        var bitIdx = ctz32(lsb);
                        var tj2 = (w << 5) | bitIdx;
                        bits ^= lsb;
                        if (tj2 < this.Nt) {
                            var tNb = this.gt.neighbors[tj2];
                            for (var rr = 0; rr < tNb.length; rr++) {
                                setBit(reachable, tNb[rr]);
                            }
                        }
                    }
                }
                for (w = 0; w < tw; w++) {
                    workDomain[qk2][w] &= reachable[w];
                }
            }
        }
    };

    // ========================================================================
    // VF3-Light ordering — sort by label rarity, then degree, then domain size
    // ========================================================================

    VF2PPMatcher.prototype.vf3LightOrder = function(order) {
        var Nq = this.Nq;
        var tw = this.tWords;
        var gq = this.gq;
        var domain = this.domain;
        var i;

        // Count label frequency in query
        var labelFreq = {};
        for (i = 0; i < Nq; i++) {
            var lbl = gq.label[i];
            labelFreq[lbl] = (labelFreq[lbl] || 0) + 1;
        }

        // Initialize identity
        for (i = 0; i < Nq; i++) order[i] = i;

        // Sort using three-key comparator
        // We need a stable-ish sort; Array.sort is fine in practice
        var domSize = new Array(Nq);
        for (i = 0; i < Nq; i++) {
            domSize[i] = popcountWords(domain[i], tw);
        }

        order.sort(function(a, b) {
            var fA = labelFreq[gq.label[a]] || 0;
            var fB = labelFreq[gq.label[b]] || 0;
            if (fA !== fB) return fA - fB;
            if (gq.degree[a] !== gq.degree[b]) return gq.degree[b] - gq.degree[a];
            return domSize[a] - domSize[b];
        });
    };

    // ========================================================================
    // Check if any query atom has empty domain
    // ========================================================================

    VF2PPMatcher.prototype.anyEmptyDomain = function() {
        var tw = this.tWords;
        for (var i = 0; i < this.Nq; i++) {
            var has = false;
            for (var w = 0; w < tw; w++) {
                if (this.domain[i][w] !== 0) { has = true; break; }
            }
            if (!has) return true;
        }
        return false;
    };

    // ========================================================================
    // Backtracking (exists — find first match)
    // ========================================================================

    VF2PPMatcher.prototype.backtrack = function(order, pos) {
        if (this.tb.expired()) { this.timedOut = true; return; }
        this.nodesVisited++;
        if (this.found) return;
        if (pos === this.Nq) { this.found = true; return; }

        var qi = order[pos];
        var buf = this.candBuf[pos];
        var nCands = this.selectCandidates(qi, buf);

        for (var c = 0; c < nCands; c++) {
            var tj = buf[c];
            if (!this.feasible(qi, tj)) continue;

            this.q2t[qi] = tj;
            this.t2q[tj] = qi;
            setBit(this.usedMask, tj);
            this.onMatch(pos, tj);

            this.backtrack(order, pos + 1);

            if (this.found || this.timedOut) return;

            this.q2t[qi] = -1;
            this.t2q[tj] = -1;
            clearBit(this.usedMask, tj);
            this.onUnmatch(pos);
        }
    };

    // ========================================================================
    // Backtracking (enumerate all matches)
    // ========================================================================

    VF2PPMatcher.prototype.enumerateRec = function(order, pos, out, maxSolutions) {
        if (this.tb.expired()) { this.timedOut = true; return; }
        this.nodesVisited++;
        if (out.length >= maxSolutions) return;

        if (pos === this.Nq) {
            // Collect mapping
            var mapping = {};
            for (var k = 0; k < this.Nq; k++) {
                mapping[order[k]] = this.q2t[order[k]];
            }
            out.push(mapping);
            return;
        }

        var qi = order[pos];
        var buf = this.candBuf[pos];
        var nCands = this.selectCandidates(qi, buf);

        for (var c = 0; c < nCands; c++) {
            if (this.tb.expired()) { this.timedOut = true; return; }
            var tj = buf[c];
            if (!this.feasible(qi, tj)) continue;

            this.q2t[qi] = tj;
            this.t2q[tj] = qi;
            setBit(this.usedMask, tj);
            this.onMatch(pos, tj);

            this.enumerateRec(order, pos + 1, out, maxSolutions);

            this.q2t[qi] = -1;
            this.t2q[tj] = -1;
            clearBit(this.usedMask, tj);
            this.onUnmatch(pos);

            if (out.length >= maxSolutions || this.timedOut) return;
        }
    };

    // ========================================================================
    // exists — does a substructure match exist?
    // ========================================================================

    VF2PPMatcher.prototype.exists = function() {
        if (this.Nq === 0) return true;
        if (this.anyEmptyDomain()) return false;

        var order = new Array(this.Nq);
        if (this.Nq > 30) {
            this.vf3LightOrder(order);
        } else {
            this.fastisoOrder(order);
        }

        // Greedy probe for small molecules
        if (this.Nq <= 50 && this.greedyProbe(order)) {
            return true;
        }

        // Reset frontier
        var w;
        for (w = 0; w < this.tWords; w++) this.frontierMask[w] = 0;

        this.backtrack(order, 0);
        return this.found;
    };

    // ========================================================================
    // enumerate — find up to maxSolutions mappings
    // ========================================================================

    VF2PPMatcher.prototype.enumerate = function(maxSolutions) {
        var out = [];
        if (this.Nq === 0) { out.push({}); return out; }

        var order = new Array(this.Nq);
        if (this.Nq > 30) {
            this.vf3LightOrder(order);
        } else {
            this.fastisoOrder(order);
        }

        // Reset frontier
        var w;
        for (w = 0; w < this.tWords; w++) this.frontierMask[w] = 0;

        this.enumerateRec(order, 0, out, maxSolutions);
        return out;
    };

    // ========================================================================
    // Public API
    // ========================================================================

    /**
     * isSubstructure(query, target, opts) -> boolean
     *
     * @param {SMSDGraph} query  — query molecular graph
     * @param {SMSDGraph} target — target molecular graph
     * @param {ChemOptions} opts — matching options (optional)
     * @returns {boolean}
     */
    function isSubstructure(query, target, opts, timeoutMs) {
        opts = opts || new SG.ChemOptions();
        if (query.n === 0) return true;
        if (query.n > target.n) return false;
        var tmo = (typeof timeoutMs === 'number' && timeoutMs > 0) ? timeoutMs
                  : (opts.timeoutMs || 10000);
        var matcher = new VF2PPMatcher(query, target, opts, tmo);
        return matcher.exists();
    }

    /**
     * findSubstructure(query, target, opts, timeoutMs) -> mapping or null
     *
     * Returns a single mapping {queryIdx: targetIdx} or null if no match.
     * Optional timeoutMs overrides opts.timeoutMs (used by MCS L1 to honour
     * the parent search budget).
     */
    function findSubstructure(query, target, opts, timeoutMs) {
        opts = opts || new SG.ChemOptions();
        if (query.n === 0) return {};
        if (query.n > target.n) return null;
        var tmo = (typeof timeoutMs === 'number' && timeoutMs > 0) ? timeoutMs
                  : (opts.timeoutMs || 10000);
        var matcher = new VF2PPMatcher(query, target, opts, tmo);
        var all = matcher.enumerate(1);
        return all.length > 0 ? all[0] : null;
    }

    /**
     * findAllSubstructures(query, target, opts, maxMappings, timeoutMs) -> array of mappings
     */
    function findAllSubstructures(query, target, opts, maxMappings, timeoutMs) {
        opts = opts || new SG.ChemOptions();
        maxMappings = maxMappings || 10000;
        if (query.n === 0) return [{}];
        if (query.n > target.n) return [];
        var tmo = (typeof timeoutMs === 'number' && timeoutMs > 0) ? timeoutMs
                  : (opts.timeoutMs || 10000);
        var matcher = new VF2PPMatcher(query, target, opts, tmo);
        return matcher.enumerate(maxMappings);
    }

    // ========================================================================
    // Export
    // ========================================================================

    window.SMSDVF2 = {
        isSubstructure: isSubstructure,
        findSubstructure: findSubstructure,
        findAllSubstructures: findAllSubstructures,
        VF2PPMatcher: VF2PPMatcher
    };

})();
