/**
 * RDT.js - Reaction Decoder Tool: pure-JavaScript port for BIME v1.2.0.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Atom-atom mapping (AAM) and bond-change annotation for reactions.
 * Faithful port of the published Reaction Decoder Tool algorithm:
 *   Rahman SA, Torrance G, Baldacci L, Cuesta SM, Fenninger F, Gopal N,
 *   Choudhary S, May JW, Holliday GL, Steinbeck C, Thornton JM.
 *   Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions.
 *   Bioinformatics 2016, 32(13):2065-66.
 *
 * Methodology reference (EC-BLAST, similarity-driven AAM):
 *   Rahman SA, Cuesta SM, Furnham N, Holliday GL, Thornton JM.
 *   EC-BLAST: A Tool to Automatically Search and Compare Enzyme Reactions.
 *   Nature Methods 2014, 11(2):171-74.
 *
 * Algorithm summary
 * -----------------
 *   1. Standardise: perceive aromaticity (via SMSDGraph), strip incoming maps.
 *   2. Build similarity matrix: MCS for every (reactant_i, product_j) pair.
 *   3. Run four pairing strategies: MIN, MAX, MIXTURE, RING.
 *   4. (Optional, v1.2.0) Refine component pairing with max-weight bipartite
 *      matching (Munkres-Kuhn) over the strategy's MCS-size matrix.
 *   5. Annotate bond-change events (formed / cleaved / orderChange /
 *      stereoChange / hydrogenChange).
 *   6. Score each candidate via the fitness function (lower = better).
 *   7. Select winner by fitness, with deterministic tie-break.
 *
 * v1.2.0 additions (all opt-in or default-true and additive)
 * ---------------------------------------------------------
 *   - options.includeHydrogens (default true): emits a 'hydrogenChange' event
 *     for every mapped atom whose implicit-H count differs across the
 *     reaction (e.g. ethanol oxidation: -1 H on the alpha-carbon, -1 H on O).
 *   - options.includeStereo (default true): runs CIP perception on each side
 *     before annotating bond changes so 'stereoChange' events fire when a
 *     mapped atom flips R<->S during the reaction (e.g. SN2 inversion).
 *   - options.useBipartitePostPass (default true since v1.3.0): after the
 *     four strategies run, re-pair reactant/product components using
 *     max-weight bipartite matching when this strictly improves total MCS
 *     coverage versus the greedy strategy pick. See
 *     `_bipartiteComponentPairing` below. To restore v1.2.x default-off
 *     behaviour, pass `{ useBipartitePostPass: false }`.
 *
 * Public API
 * ----------
 *   RDT.mapReaction(reaction, options) -> { mapping, bondChanges, score, strategy, ... }
 *   RDT.runMinPairwise(reaction, options)
 *   RDT.runMaxPairwise(reaction, options)
 *   RDT.runMixturePairwise(reaction, options)
 *   RDT.runRingPairwise(reaction, options)
 *   RDT.annotateBondChanges(reactantSide, productSide, atomMap, opts) -> bondChanges[]
 *   RDT.fitness(reaction, atomMap, bondChanges) -> number
 *   RDT._bipartiteComponentPairing(rN, pN, weights) -> [{reactantIdx,productIdx},...]
 */
(function(global) {
    'use strict';

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    var STRATEGIES = ['MIN', 'MAX', 'MIXTURE', 'RING'];

    // Tie-break order when two strategies post equal fitness:
    //   Lex order chosen to be deterministic (alphabetical).
    var TIEBREAK_ORDER = { MAX: 0, MIN: 1, MIXTURE: 2, RING: 3 };

    var DEFAULT_TIMEOUT_MS = 10000;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    function nowMs() {
        return (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    }

    function shallowCopyOpts(opts) {
        var out = {};
        if (!opts) { return out; }
        for (var k in opts) {
            if (opts.hasOwnProperty(k)) { out[k] = opts[k]; }
        }
        return out;
    }

    // Split a Molecule with a reaction arrow into reactant and product side
    // arrays. Each side is a list of independent Molecule objects, one per
    // connected component on that side of the arrow midpoint.
    function splitReactionSides(reaction) {
        if (!reaction) {
            return { reactants: [], products: [], reactantAtomIds: [], productAtomIds: [] };
        }
        if (!reaction.reactionArrow) {
            // No arrow: treat whole molecule as a "reactant" side, no products.
            return {
                reactants: extractComponents(reaction, function() { return true; }),
                products: [],
                reactantAtomIds: reaction.atoms.map(function(a) { return a.id; }),
                productAtomIds: []
            };
        }
        var arrow = reaction.reactionArrow;
        var midX = (arrow.x1 + arrow.x2) / 2;

        var reactantAtomIds = [];
        var productAtomIds = [];
        var reactantSet = {};
        var productSet = {};

        // Determine each *connected component*'s side once (centroid) so that
        // all atoms of a component end up on the same side, even if individual
        // atom x-coordinates straddle the midpoint slightly.
        var components = reaction.getComponents();
        for (var ci = 0; ci < components.length; ci++) {
            var comp = components[ci];
            var sumX = 0;
            for (var ai = 0; ai < comp.length; ai++) {
                var atom = reaction.getAtom(comp[ai]);
                sumX += atom ? atom.x : 0;
            }
            var cx = sumX / comp.length;
            var isReactant = cx < midX;
            for (var bi = 0; bi < comp.length; bi++) {
                if (isReactant) {
                    reactantAtomIds.push(comp[bi]);
                    reactantSet[comp[bi]] = true;
                } else {
                    productAtomIds.push(comp[bi]);
                    productSet[comp[bi]] = true;
                }
            }
        }

        var reactants = extractComponents(reaction, function(atomId) { return !!reactantSet[atomId]; });
        var products = extractComponents(reaction, function(atomId) { return !!productSet[atomId]; });

        return {
            reactants: reactants,
            products: products,
            reactantAtomIds: reactantAtomIds,
            productAtomIds: productAtomIds
        };
    }

    // Build new Molecule objects per connected component for atoms that pass
    // the predicate. The new molecule re-uses the ORIGINAL atom and bond IDs
    // by replaying them via direct field assignment (so that downstream code
    // can reference back to the source reaction by atom/bond id).
    function extractComponents(source, atomPredicate) {
        var Mol = global.Molecule;
        var visited = {};
        var componentsOut = [];
        var i, j;

        // Group atoms into connected components within the predicate-filtered subgraph.
        function neighborsInSubgraph(atomId) {
            var bonds = source.getBondsOfAtom(atomId);
            var out = [];
            for (var bi = 0; bi < bonds.length; bi++) {
                var other = bonds[bi].otherAtom(atomId);
                if (atomPredicate(other)) { out.push(other); }
            }
            return out;
        }

        var allAtomIds = [];
        for (i = 0; i < source.atoms.length; i++) {
            if (atomPredicate(source.atoms[i].id)) { allAtomIds.push(source.atoms[i].id); }
        }
        // Sort for deterministic component order
        allAtomIds.sort(function(a, b) { return a - b; });

        for (i = 0; i < allAtomIds.length; i++) {
            var startId = allAtomIds[i];
            if (visited[startId]) { continue; }
            var stack = [startId];
            var compIds = [];
            visited[startId] = true;
            while (stack.length > 0) {
                var cur = stack.pop();
                compIds.push(cur);
                var nbs = neighborsInSubgraph(cur);
                for (var ni = 0; ni < nbs.length; ni++) {
                    if (!visited[nbs[ni]]) {
                        visited[nbs[ni]] = true;
                        stack.push(nbs[ni]);
                    }
                }
            }
            compIds.sort(function(a, b) { return a - b; });

            // Build a Molecule for this component preserving original ids
            var sub = new Mol();
            var idSet = {};
            for (j = 0; j < compIds.length; j++) {
                var src = source.getAtom(compIds[j]);
                if (!src) { continue; }
                // Manually push to keep original id (addAtom would mint a new one)
                var atom = src.clone();
                atom.id = src.id;
                sub.atoms.push(atom);
                sub._atomMap[atom.id] = atom;
                sub._adjacency[atom.id] = [];
                idSet[atom.id] = true;
            }
            // Add bonds whose endpoints are both inside the component
            var seenBond = {};
            for (j = 0; j < compIds.length; j++) {
                var bonds = source.getBondsOfAtom(compIds[j]);
                for (var bi2 = 0; bi2 < bonds.length; bi2++) {
                    var bond = bonds[bi2];
                    if (seenBond[bond.id]) { continue; }
                    if (!idSet[bond.atom1] || !idSet[bond.atom2]) { continue; }
                    seenBond[bond.id] = true;
                    var nb = bond.clone();
                    nb.id = bond.id;
                    nb.atom1 = bond.atom1;
                    nb.atom2 = bond.atom2;
                    sub.bonds.push(nb);
                    sub._bondMap[nb.id] = nb;
                    sub._adjacency[nb.atom1].push(nb.id);
                    sub._adjacency[nb.atom2].push(nb.id);
                }
            }
            componentsOut.push(sub);
        }
        return componentsOut;
    }

    // Tally element symbol counts in a list of Molecule components (heavy atoms only).
    function elementCounts(components) {
        var counts = {};
        for (var i = 0; i < components.length; i++) {
            var atoms = components[i].atoms;
            for (var j = 0; j < atoms.length; j++) {
                var sym = atoms[j].symbol;
                counts[sym] = (counts[sym] || 0) + 1;
            }
        }
        return counts;
    }

    function elementCountsEqual(a, b) {
        var k;
        for (k in a) { if (a.hasOwnProperty(k) && (a[k] || 0) !== (b[k] || 0)) { return false; } }
        for (k in b) { if (b.hasOwnProperty(k) && (a[k] || 0) !== (b[k] || 0)) { return false; } }
        return true;
    }

    // Extract atoms with mapNumber > 0 from a side; returns { mapNum -> [atomId, ...] }
    function collectPreMaps(side) {
        var out = {};
        for (var i = 0; i < side.length; i++) {
            var atoms = side[i].atoms;
            for (var j = 0; j < atoms.length; j++) {
                var n = atoms[j].mapNumber || 0;
                if (n > 0) {
                    if (!out[n]) { out[n] = []; }
                    out[n].push(atoms[j].id);
                }
            }
        }
        return out;
    }

    // Strip incoming map numbers on the working reaction (we recompute them).
    function stripMapNumbers(reaction) {
        for (var i = 0; i < reaction.atoms.length; i++) { reaction.atoms[i].mapNumber = 0; }
    }

    // -----------------------------------------------------------------------
    // Cache of MCS results between (reactant_i, product_j)
    // -----------------------------------------------------------------------

    function buildSimilarityCache(reactants, products, chemOpts, mcsOpts, ringOnly) {
        var cache = [];
        var SG = global.SMSDGraph && global.SMSDGraph.SMSDGraph;
        if (!SG || !global.SMSDMCS || !global.SMSDMCS.findMCS) {
            throw new Error('RDT requires SMSDGraph and SMSDMCS modules');
        }

        var cOpts = shallowCopyOpts(chemOpts);
        cOpts.matchAtomType = (cOpts.matchAtomType !== undefined) ? cOpts.matchAtomType : true;
        cOpts.matchBondType = (cOpts.matchBondType !== undefined) ? cOpts.matchBondType : false;
        // Charge changes are normal across reactions (ionisation, protonation):
        cOpts.matchFormalCharge = false;
        if (ringOnly) { cOpts.ringMatchesRingOnly = true; }

        var mOpts = shallowCopyOpts(mcsOpts);
        if (mOpts.timeoutMs === undefined) { mOpts.timeoutMs = 2000; }

        for (var i = 0; i < reactants.length; i++) {
            cache[i] = [];
            var g1 = new SG(reactants[i]);
            for (var j = 0; j < products.length; j++) {
                var g2 = new SG(products[j]);
                var mcs;
                try {
                    mcs = global.SMSDMCS.findMCS(g1, g2, cOpts, mOpts);
                } catch (e) {
                    mcs = { mapping: {}, size: 0 };
                }
                var atomMapping = global.SMSDMCS.translateToAtomIds(mcs.mapping || {}, g1, g2);
                cache[i][j] = {
                    mcsSize: mcs.size || 0,
                    mapping: atomMapping,
                    rIdx: i,
                    pIdx: j,
                    rSize: g1.n,
                    pSize: g2.n,
                    ringSize: countRingMappedAtoms(reactants[i], atomMapping, true),
                    chainSize: (mcs.size || 0) - countRingMappedAtoms(reactants[i], atomMapping, true)
                };
            }
        }
        return cache;
    }

    // Count atoms in the mapping that lie in a ring on the reactant side.
    function countRingMappedAtoms(reactantMol, atomMapping, _unused) {
        if (!atomMapping) { return 0; }
        var ringSet = {};
        var rings = reactantMol.findRings ? reactantMol.findRings(8) : [];
        for (var i = 0; i < rings.length; i++) {
            for (var j = 0; j < rings[i].atoms.length; j++) {
                ringSet[rings[i].atoms[j]] = true;
            }
        }
        var n = 0;
        for (var aid in atomMapping) {
            if (atomMapping.hasOwnProperty(aid) && ringSet[+aid]) { n++; }
        }
        return n;
    }

    // Whether either side of a (reactant_i, product_j) pair contains a ring.
    function pairHasRing(reactants, products, i, j) {
        var rRings = reactants[i].findRings ? reactants[i].findRings(8) : [];
        if (rRings.length > 0) { return true; }
        var pRings = products[j].findRings ? products[j].findRings(8) : [];
        return pRings.length > 0;
    }

    // -----------------------------------------------------------------------
    // Strategy implementations: each picks (reactant_i, product_j) pairs in
    // a deterministic order driven by the cache, removing chosen indices
    // from the available pool until exhausted.
    // -----------------------------------------------------------------------

    // Generic greedy pair selection driven by a comparator(a, b) function.
    // comparator returns negative if a should be picked before b.
    function selectPairs(cache, reactants, products, comparator, filter) {
        var rUsed = {};
        var pUsed = {};
        var pairs = [];
        var candidates = [];

        for (var i = 0; i < reactants.length; i++) {
            for (var j = 0; j < products.length; j++) {
                if (filter && !filter(i, j, cache[i][j])) { continue; }
                candidates.push(cache[i][j]);
            }
        }
        // Stable, deterministic sort
        candidates.sort(function(a, b) {
            var c = comparator(a, b);
            if (c !== 0) { return c; }
            // Tie-break by indices for full determinism
            if (a.rIdx !== b.rIdx) { return a.rIdx - b.rIdx; }
            return a.pIdx - b.pIdx;
        });

        for (var k = 0; k < candidates.length; k++) {
            var c = candidates[k];
            if (rUsed[c.rIdx] || pUsed[c.pIdx]) { continue; }
            if (c.mcsSize === 0) { continue; } // nothing to map
            rUsed[c.rIdx] = true;
            pUsed[c.pIdx] = true;
            pairs.push(c);
        }
        return pairs;
    }

    // Construct a final atom-id mapping object (reactantAtomId -> productAtomId)
    // by unioning every chosen pair's mapping.
    function pairsToAtomMap(pairs) {
        var atomMap = {};
        for (var k = 0; k < pairs.length; k++) {
            var pm = pairs[k].mapping || {};
            for (var aid in pm) {
                if (pm.hasOwnProperty(aid)) {
                    // First wins: deterministic via candidate sort order
                    if (!atomMap.hasOwnProperty(aid)) { atomMap[aid] = pm[aid]; }
                }
            }
        }
        return atomMap;
    }

    // -----------------------------------------------------------------------
    // Max-weight bipartite component pairing — refinement of the published
    // RDT 2016 strategies' greedy per-component selection. Standard textbook
    // Munkres-Kuhn assignment (Kuhn HW. The Hungarian method for the
    // assignment problem. Naval Research Logistics Quarterly 1955;2:83-97.
    // Munkres J. Algorithms for the assignment and transportation problems.
    // SIAM J 1957;5:32-38). Cost matrix entries are negative MCS sizes from
    // the SMSDMCS pair cache, so a min-cost assignment maximises MCS coverage.
    //
    // ON by default since v1.3.0 (was opt-in in v1.2.0). The dispatcher
    // re-pairs components iff the bipartite assignment strictly exceeds the
    // strategy-greedy total MCS size — so on reactions where greedy is
    // already optimal, this is a no-op. Pass `{ useBipartitePostPass: false }`
    // to restore the v1.2.x default-off behaviour.
    //
    // Implementation: square the cost matrix to size N = max(rN, pN) by
    // padding with zero-weight dummy rows/columns, then run the textbook
    // O(N^3) Hungarian algorithm. For typical reactions N <= 8 so this
    // completes in microseconds.
    // -----------------------------------------------------------------------
    function _bipartiteComponentPairing(rN, pN, weights) {
        if (rN === 0 || pN === 0) { return []; }
        // Build square N x N cost matrix; pad with zeroes for unmatched slots.
        var N = Math.max(rN, pN);
        var INF = 1e9;
        var cost = new Array(N);
        var i, j;
        for (i = 0; i < N; i++) {
            cost[i] = new Array(N);
            for (j = 0; j < N; j++) {
                if (i < rN && j < pN) {
                    var w = (weights[i] && typeof weights[i][j] === 'number') ? weights[i][j] : 0;
                    cost[i][j] = -w;            // maximise w  =>  minimise -w
                } else {
                    cost[i][j] = 0;             // dummy padding (free assignment)
                }
            }
        }

        // Textbook O(N^3) Hungarian (Munkres-Kuhn) — minimise total cost.
        // Uses 1-indexed potentials u/v and column slack vectors as in
        // standard references (e.g. Cormen et al., e-maxx tutorial).
        var u = new Array(N + 1);
        var v = new Array(N + 1);
        var p = new Array(N + 1);  // p[j] = row assigned to column j
        var way = new Array(N + 1);
        for (i = 0; i <= N; i++) { u[i] = 0; v[i] = 0; p[i] = 0; way[i] = 0; }

        for (i = 1; i <= N; i++) {
            p[0] = i;
            var j0 = 0;
            var minv = new Array(N + 1);
            var used = new Array(N + 1);
            for (j = 0; j <= N; j++) { minv[j] = INF; used[j] = false; }
            do {
                used[j0] = true;
                var i0 = p[j0];
                var delta = INF;
                var j1 = 0;
                for (j = 1; j <= N; j++) {
                    if (!used[j]) {
                        var cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                        if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                        if (minv[j] < delta) { delta = minv[j]; j1 = j; }
                    }
                }
                for (j = 0; j <= N; j++) {
                    if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                    else         { minv[j] -= delta; }
                }
                j0 = j1;
            } while (p[j0] !== 0);
            do {
                var j2 = way[j0];
                p[j0] = p[j2];
                j0 = j2;
            } while (j0 !== 0);
        }

        // p[j] = row index assigned to column j (1-based). Flip back to
        // 0-based and drop dummy assignments (where i >= rN or j >= pN).
        var pairs = [];
        for (j = 1; j <= N; j++) {
            var rIdx = p[j] - 1;
            var pIdx = j - 1;
            if (rIdx < 0 || rIdx >= rN || pIdx >= pN) { continue; }
            // Skip pairs with zero MCS — they represent "no real match".
            var ww = (weights[rIdx] && typeof weights[rIdx][pIdx] === 'number') ? weights[rIdx][pIdx] : 0;
            if (ww <= 0) { continue; }
            pairs.push({ reactantIdx: rIdx, productIdx: pIdx, mcsSize: ww });
        }
        // Deterministic ordering by reactant index, then product index.
        pairs.sort(function(a, b) {
            if (a.reactantIdx !== b.reactantIdx) { return a.reactantIdx - b.reactantIdx; }
            return a.productIdx - b.productIdx;
        });
        return pairs;
    }

    // Apply the bipartite post-pass to a strategy result, mutating in place
    // ONLY if total MCS strictly improves. Returns true iff the result was
    // re-paired. Pulls component->component MCS atom maps from result.cache.
    function applyBipartitePostPass(result, options) {
        if (!result || !result.cache || !result.sides) { return false; }
        var cache = result.cache;
        var rN = result.sides.reactants.length;
        var pN = result.sides.products.length;
        if (rN === 0 || pN === 0) { return false; }

        // Greedy total = sum of MCS sizes from the strategy's chosen pairs.
        var greedyTotal = 0;
        var i;
        for (i = 0; i < result.pairs.length; i++) { greedyTotal += (result.pairs[i].mcsSize || 0); }

        // Build weight matrix (MCS sizes) from the cache.
        var weights = [];
        for (i = 0; i < rN; i++) {
            weights[i] = [];
            for (var j = 0; j < pN; j++) {
                weights[i][j] = (cache[i] && cache[i][j]) ? (cache[i][j].mcsSize || 0) : 0;
            }
        }

        var bipPairs = _bipartiteComponentPairing(rN, pN, weights);
        var bipTotal = 0;
        for (i = 0; i < bipPairs.length; i++) { bipTotal += (bipPairs[i].mcsSize || 0); }

        if (bipTotal <= greedyTotal) { return false; }

        // Convert bipartite (rIdx,pIdx) pairs back into the cache cells the
        // rest of the pipeline expects.
        var newPairs = [];
        for (i = 0; i < bipPairs.length; i++) {
            var bp = bipPairs[i];
            var cell = cache[bp.reactantIdx][bp.productIdx];
            if (cell && cell.mcsSize > 0) { newPairs.push(cell); }
        }
        result.pairs = newPairs;
        result.mapping = pairsToAtomMap(newPairs);
        result.mapping = mergePreMaps(result.mapping, options._preMaps);
        var bcOpts = {
            includeHydrogens: (options.includeHydrogens !== false),
            includeStereo: (options.includeStereo !== false)
        };
        result.bondChanges = annotateBondChanges(result.sides.reactants, result.sides.products, result.mapping, bcOpts);
        result.score = fitness({ reactants: result.sides.reactants, products: result.sides.products },
            result.mapping, result.bondChanges);
        result.bipartiteApplied = true;
        result.bipartiteGreedyTotal = greedyTotal;
        result.bipartiteTotal = bipTotal;
        return true;
    }

    // -----------------------------------------------------------------------
    // Bond-change annotation
    // -----------------------------------------------------------------------

    function buildBondLookupBySide(side) {
        // bondMap[atomA + ',' + atomB] = bond (with smaller atom id first)
        var lookup = {};
        for (var i = 0; i < side.length; i++) {
            var bonds = side[i].bonds;
            for (var j = 0; j < bonds.length; j++) {
                var b = bonds[j];
                var a = b.atom1 < b.atom2 ? b.atom1 : b.atom2;
                var c = b.atom1 < b.atom2 ? b.atom2 : b.atom1;
                lookup[a + ',' + c] = b;
            }
        }
        return lookup;
    }

    // Given an atomMap (reactantAtomId -> productAtomId), produce a list of
    // bond-change events covering both sides. Atoms not present in atomMap
    // are unmapped; we record events conservatively.
    //
    // Optional `opts`:
    //   includeHydrogens (default true) — emit 'hydrogenChange' events when a
    //     mapped atom's implicit-H count changes (oxidation, reduction, etc).
    //   includeStereo (default true) — emit 'stereoChange' events when a
    //     mapped atom's CIP label (R/S) flips. CIP labels must already be
    //     assigned on the source molecules; mapReaction() does this for you.
    function annotateBondChanges(reactantSide, productSide, atomMap, opts) {
        opts = opts || {};
        var includeHydrogens = (opts.includeHydrogens !== false);
        var includeStereo = (opts.includeStereo !== false);

        var rBonds = buildBondLookupBySide(reactantSide);
        var pBonds = buildBondLookupBySide(productSide);
        var events = [];
        var seenProductBonds = {};

        // Atom + parent-molecule lookup for stereo / H info.
        var rAtoms = {};
        var rAtomMol = {};
        var i;
        for (i = 0; i < reactantSide.length; i++) {
            for (var ji = 0; ji < reactantSide[i].atoms.length; ji++) {
                rAtoms[reactantSide[i].atoms[ji].id] = reactantSide[i].atoms[ji];
                rAtomMol[reactantSide[i].atoms[ji].id] = reactantSide[i];
            }
        }
        var pAtoms = {};
        var pAtomMol = {};
        for (i = 0; i < productSide.length; i++) {
            for (var ji2 = 0; ji2 < productSide[i].atoms.length; ji2++) {
                pAtoms[productSide[i].atoms[ji2].id] = productSide[i].atoms[ji2];
                pAtomMol[productSide[i].atoms[ji2].id] = productSide[i];
            }
        }

        // Pass 1: walk reactant bonds, classify formed/cleaved/orderChange
        var rBondIds = Object.keys(rBonds);
        rBondIds.sort();
        for (i = 0; i < rBondIds.length; i++) {
            var rb = rBonds[rBondIds[i]];
            var pa1 = atomMap.hasOwnProperty(rb.atom1) ? atomMap[rb.atom1] : null;
            var pa2 = atomMap.hasOwnProperty(rb.atom2) ? atomMap[rb.atom2] : null;
            if (pa1 == null || pa2 == null) {
                // At least one endpoint is unmapped: bond is effectively cleaved.
                events.push({
                    type: 'cleaved',
                    reactantBond: rb.id,
                    productBond: null,
                    atoms: [rb.atom1, null],
                    beforeOrder: rb.type,
                    afterOrder: null
                });
                continue;
            }
            var key = pa1 < pa2 ? (pa1 + ',' + pa2) : (pa2 + ',' + pa1);
            var pb = pBonds[key];
            if (!pb) {
                events.push({
                    type: 'cleaved',
                    reactantBond: rb.id,
                    productBond: null,
                    atoms: [rb.atom1, pa1],
                    beforeOrder: rb.type,
                    afterOrder: null
                });
                continue;
            }
            seenProductBonds[pb.id] = true;
            if (pb.type !== rb.type) {
                events.push({
                    type: 'orderChange',
                    reactantBond: rb.id,
                    productBond: pb.id,
                    atoms: [rb.atom1, pa1],
                    beforeOrder: rb.type,
                    afterOrder: pb.type
                });
            }
        }

        // Pass 2: walk product bonds; any not covered = formed.
        var pBondIds = Object.keys(pBonds);
        pBondIds.sort();
        for (i = 0; i < pBondIds.length; i++) {
            var pb2 = pBonds[pBondIds[i]];
            if (seenProductBonds[pb2.id]) { continue; }
            // Need both endpoints to be in the *image* of atomMap (i.e. mapped from reactant side).
            // If both endpoints are mapped from the reactant side then it's a formed bond.
            // Otherwise, it's a bond between unmapped product atoms (e.g. unbalanced reaction)
            // and we skip it to avoid spurious events.
            var pAtomsMappedFromR = false;
            var p1Mapped = false, p2Mapped = false;
            for (var rk in atomMap) {
                if (atomMap.hasOwnProperty(rk)) {
                    if (atomMap[rk] === pb2.atom1) { p1Mapped = true; }
                    if (atomMap[rk] === pb2.atom2) { p2Mapped = true; }
                }
            }
            if (p1Mapped && p2Mapped) { pAtomsMappedFromR = true; }
            if (!pAtomsMappedFromR) { continue; }

            events.push({
                type: 'formed',
                reactantBond: null,
                productBond: pb2.id,
                atoms: [null, pb2.atom1],
                beforeOrder: null,
                afterOrder: pb2.type
            });
        }

        // Pass 3: stereo changes at mapped atoms (CIP R<->S flip).
        // CIP perception is performed by mapReaction() on each side before
        // this pass runs; if callers invoke annotateBondChanges directly they
        // are responsible for CipStereo.assign() if they want stereo events.
        if (includeStereo) {
            var rIds = Object.keys(atomMap);
            rIds.sort(function(a, b) { return (+a) - (+b); });
            for (i = 0; i < rIds.length; i++) {
                var rid = +rIds[i];
                var pid = atomMap[rid];
                var ra = rAtoms[rid];
                var pa = pAtoms[pid];
                if (!ra || !pa) { continue; }
                var beforeRS = ra.cipLabel || '';
                var afterRS = pa.cipLabel || '';
                if (beforeRS !== afterRS && (beforeRS || afterRS)) {
                    events.push({
                        type: 'stereoChange',
                        reactantBond: null,
                        productBond: null,
                        atoms: [rid, pid],
                        beforeOrder: null,
                        afterOrder: null,
                        beforeStereo: beforeRS,
                        afterStereo: afterRS,
                        before: beforeRS,
                        after: afterRS,
                        atom: rid,
                        productAtom: pid,
                        mapNumber: ra.mapNumber || 0
                    });
                }
            }
        }

        // Pass 4 (v1.2.0): hydrogen-count changes at mapped atoms.
        // For oxidation / reduction / dehydration the heavy-atom skeleton is
        // unchanged but each mapped atom's implicit-H count shifts. We emit a
        // 'hydrogenChange' event whose deltaH = hCountProduct - hCountReactant.
        // Molecule.calcHydrogens(atomId) returns the implicit-H count derived
        // from the standard valence rules already in BIME.
        if (includeHydrogens) {
            var hIds = Object.keys(atomMap);
            hIds.sort(function(a, b) { return (+a) - (+b); });
            for (i = 0; i < hIds.length; i++) {
                var rhid = +hIds[i];
                var phid = atomMap[rhid];
                var rha = rAtoms[rhid];
                var pha = pAtoms[phid];
                if (!rha || !pha) { continue; }
                var rmol = rAtomMol[rhid];
                var pmol = pAtomMol[phid];
                if (!rmol || !pmol) { continue; }
                if (typeof rmol.calcHydrogens !== 'function' ||
                    typeof pmol.calcHydrogens !== 'function') { continue; }
                var hR, hP;
                try {
                    hR = rmol.calcHydrogens(rhid);
                    hP = pmol.calcHydrogens(phid);
                } catch (e) { continue; }
                if (typeof hR !== 'number' || typeof hP !== 'number') { continue; }
                var dH = hP - hR;
                if (dH !== 0) {
                    events.push({
                        type: 'hydrogenChange',
                        reactantBond: null,
                        productBond: null,
                        atoms: [rhid, phid],
                        beforeOrder: null,
                        afterOrder: null,
                        atom: rhid,
                        productAtom: phid,
                        deltaH: dH,
                        beforeH: hR,
                        afterH: hP,
                        mapNumber: rha.mapNumber || 0
                    });
                }
            }
        }

        return events;
    }

    // -----------------------------------------------------------------------
    // Fitness function
    // -----------------------------------------------------------------------

    // Compute the fitness score for a candidate. Lower is better.
    function fitness(reactionOrSides, atomMap, bondChanges) {
        var sides = reactionOrSides;
        if (sides && sides.atoms && sides.bonds) {
            sides = splitReactionSides(sides);
        }
        var rSide = sides.reactants || [];
        var pSide = sides.products || [];

        var changeCount = 0;
        for (var i = 0; i < bondChanges.length; i++) {
            var ev = bondChanges[i];
            if (ev.type === 'formed' || ev.type === 'cleaved') { changeCount += 1; }
            else if (ev.type === 'orderChange') { changeCount += 0.5; }
            else if (ev.type === 'stereoChange') { changeCount += 0.25; }
        }

        // Ring-preservation bonus: +0.5 per ring-atom kept in a ring.
        var rRingSet = {};
        var pRingSet = {};
        var k, ringIdx, atomsArr;
        for (k = 0; k < rSide.length; k++) {
            var rRings = rSide[k].findRings ? rSide[k].findRings(8) : [];
            for (ringIdx = 0; ringIdx < rRings.length; ringIdx++) {
                atomsArr = rRings[ringIdx].atoms;
                for (var ai = 0; ai < atomsArr.length; ai++) { rRingSet[atomsArr[ai]] = true; }
            }
        }
        for (k = 0; k < pSide.length; k++) {
            var pRings = pSide[k].findRings ? pSide[k].findRings(8) : [];
            for (ringIdx = 0; ringIdx < pRings.length; ringIdx++) {
                atomsArr = pRings[ringIdx].atoms;
                for (var aj = 0; aj < atomsArr.length; aj++) { pRingSet[atomsArr[aj]] = true; }
            }
        }
        var ringPreserved = 0;
        for (var rid in atomMap) {
            if (atomMap.hasOwnProperty(rid)) {
                if (rRingSet[+rid] && pRingSet[atomMap[rid]]) { ringPreserved += 1; }
            }
        }

        // Chemical-filter penalty: aromatic mapped to non-aromatic.
        var rAtoms = {};
        var pAtoms = {};
        for (k = 0; k < rSide.length; k++) {
            for (var ki = 0; ki < rSide[k].atoms.length; ki++) {
                rAtoms[rSide[k].atoms[ki].id] = rSide[k].atoms[ki];
            }
        }
        for (k = 0; k < pSide.length; k++) {
            for (var kj = 0; kj < pSide[k].atoms.length; kj++) {
                pAtoms[pSide[k].atoms[kj].id] = pSide[k].atoms[kj];
            }
        }
        var aromaticMismatches = 0;
        for (var rid2 in atomMap) {
            if (atomMap.hasOwnProperty(rid2)) {
                var ra = rAtoms[+rid2];
                var pa = pAtoms[atomMap[rid2]];
                if (ra && pa && ra.aromatic && !pa.aromatic) { aromaticMismatches += 1; }
                if (ra && pa && !ra.aromatic && pa.aromatic) { aromaticMismatches += 1; }
            }
        }

        return changeCount - 0.5 * ringPreserved + 1.0 * aromaticMismatches;
    }

    // -----------------------------------------------------------------------
    // Strategy runners
    // -----------------------------------------------------------------------

    // Build a candidate {atomMap, pairs} from a comparator-based strategy.
    function runStrategy(reaction, options, strategyName) {
        var sides = splitReactionSides(reaction);
        var reactants = sides.reactants;
        var products = sides.products;

        // Empty sides: graceful no-op
        if (reactants.length === 0 || products.length === 0) {
            return {
                strategy: strategyName,
                mapping: {},
                bondChanges: [],
                score: 0,
                pairs: [],
                sides: sides
            };
        }

        var chemOpts = options.chemOpts || {};
        var mcsOpts = options.mcsOpts || {};
        var ringOnly = (strategyName === 'RING');
        var cache;
        try {
            cache = buildSimilarityCache(reactants, products, chemOpts, mcsOpts, ringOnly);
        } catch (e) {
            return {
                strategy: strategyName,
                mapping: {},
                bondChanges: [],
                score: Number.POSITIVE_INFINITY,
                pairs: [],
                sides: sides,
                error: e.message
            };
        }

        var comparator;
        var filter = null;
        switch (strategyName) {
            case 'MIN':
                // Smallest non-zero MCS first
                comparator = function(a, b) { return a.mcsSize - b.mcsSize; };
                filter = function(_i, _j, cell) { return cell.mcsSize > 0; };
                break;
            case 'MAX':
                // Largest MCS first
                comparator = function(a, b) { return b.mcsSize - a.mcsSize; };
                break;
            case 'MIXTURE':
                // Largest MCS-relative-to-min(|R_i|,|P_j|) first
                comparator = function(a, b) {
                    var ar = a.mcsSize / Math.max(1, Math.min(a.rSize, a.pSize));
                    var br = b.mcsSize / Math.max(1, Math.min(b.rSize, b.pSize));
                    if (br !== ar) { return br - ar; }
                    return b.mcsSize - a.mcsSize;
                };
                break;
            case 'RING':
                // Ring-constrained pairs only; largest MCS first.
                filter = function(i, j, _cell) { return pairHasRing(reactants, products, i, j); };
                comparator = function(a, b) {
                    if (b.ringSize !== a.ringSize) { return b.ringSize - a.ringSize; }
                    return b.mcsSize - a.mcsSize;
                };
                break;
            default:
                comparator = function(a, b) { return b.mcsSize - a.mcsSize; };
        }

        var pairs = selectPairs(cache, reactants, products, comparator, filter);

        // After the strategy's primary picks, greedy-fill any remaining (i, j)
        // slots so each component participates in the mapping. The fill pass
        // uses MAX-style ranking on whatever cells are still available.
        // This matches the published RDT spirit: each strategy commits to its
        // primary choice first, then assigns remaining components by overlap.
        if (pairs.length < Math.min(reactants.length, products.length)) {
            var rUsed2 = {}, pUsed2 = {};
            for (var pi = 0; pi < pairs.length; pi++) {
                rUsed2[pairs[pi].rIdx] = true;
                pUsed2[pairs[pi].pIdx] = true;
            }
            var fill = selectPairs(cache, reactants, products,
                function(a, b) { return b.mcsSize - a.mcsSize; },
                function(i, j, cell) { return !rUsed2[i] && !pUsed2[j] && cell.mcsSize > 0; });
            pairs = pairs.concat(fill);
        }

        var atomMap = pairsToAtomMap(pairs);

        // Honor pre-existing atom maps in the input reaction:
        // if user pre-mapped any (rAtomId -> mapNumber) and a product atom carries
        // the same mapNumber, prefer that pairing.
        atomMap = mergePreMaps(atomMap, options._preMaps);

        var bcOpts = {
            includeHydrogens: (options.includeHydrogens !== false),
            includeStereo: (options.includeStereo !== false)
        };
        var bondChanges = annotateBondChanges(sides.reactants, sides.products, atomMap, bcOpts);
        var score = fitness({ reactants: sides.reactants, products: sides.products }, atomMap, bondChanges);

        return {
            strategy: strategyName,
            mapping: atomMap,
            bondChanges: bondChanges,
            score: score,
            pairs: pairs,
            sides: sides,
            cache: cache
        };
    }

    function mergePreMaps(atomMap, preMaps) {
        if (!preMaps) { return atomMap; }
        var merged = {};
        for (var k in atomMap) { if (atomMap.hasOwnProperty(k)) { merged[k] = atomMap[k]; } }
        // preMaps is a list of [rAtomId, pAtomId] pairs derived from user-supplied mapNumbers.
        for (var i = 0; i < preMaps.length; i++) {
            var pair = preMaps[i];
            // Override any conflicting auto-mapping with the user's pre-map.
            // (Overwrite previous reactant->product entry if it conflicts.)
            // Remove any entries pointing to pair[1] from a different reactant.
            for (var rk in merged) {
                if (merged.hasOwnProperty(rk) && merged[rk] === pair[1] && (+rk) !== pair[0]) {
                    delete merged[rk];
                }
            }
            merged[pair[0]] = pair[1];
        }
        return merged;
    }

    // -----------------------------------------------------------------------
    // Public entry points
    // -----------------------------------------------------------------------

    function runMinPairwise(reaction, options) { return runStrategy(reaction, options || {}, 'MIN'); }
    function runMaxPairwise(reaction, options) { return runStrategy(reaction, options || {}, 'MAX'); }
    function runMixturePairwise(reaction, options) { return runStrategy(reaction, options || {}, 'MIXTURE'); }
    function runRingPairwise(reaction, options) { return runStrategy(reaction, options || {}, 'RING'); }

    function pickWinner(results) {
        // Lower fitness wins. Ties broken by:
        //   1. more atoms mapped (completeness)
        //   2. fewer total bond-change events (formed+cleaved+orderChange)
        //   3. more ring-preservation
        //   4. fixed strategy name order (TIEBREAK_ORDER)
        var best = null;
        for (var i = 0; i < results.length; i++) {
            var r = results[i];
            if (!isFinite(r.score)) { continue; }
            if (!best) { best = r; continue; }
            if (r.score < best.score - 1e-9) { best = r; continue; }
            if (Math.abs(r.score - best.score) <= 1e-9) {
                var rMapped = Object.keys(r.mapping).length;
                var bMapped = Object.keys(best.mapping).length;
                if (rMapped > bMapped) { best = r; continue; }
                if (rMapped === bMapped) {
                    var rEvents = countEvents(r.bondChanges);
                    var bEvents = countEvents(best.bondChanges);
                    if (rEvents < bEvents) { best = r; continue; }
                    if (rEvents === bEvents) {
                        var rRing = countRingPreserved(r);
                        var bRing = countRingPreserved(best);
                        if (rRing > bRing) { best = r; continue; }
                        if (rRing === bRing) {
                            if (TIEBREAK_ORDER[r.strategy] < TIEBREAK_ORDER[best.strategy]) { best = r; }
                        }
                    }
                }
            }
        }
        return best;
    }

    function countEvents(events) {
        var n = 0;
        for (var i = 0; i < events.length; i++) {
            var t = events[i].type;
            if (t === 'formed' || t === 'cleaved' || t === 'orderChange') { n++; }
        }
        return n;
    }

    function countRingPreserved(result) {
        var sides = result.sides;
        var rRingSet = {};
        var pRingSet = {};
        var rings, j;
        for (var k = 0; k < sides.reactants.length; k++) {
            rings = sides.reactants[k].findRings ? sides.reactants[k].findRings(8) : [];
            for (j = 0; j < rings.length; j++) {
                for (var a = 0; a < rings[j].atoms.length; a++) { rRingSet[rings[j].atoms[a]] = true; }
            }
        }
        for (var kp = 0; kp < sides.products.length; kp++) {
            rings = sides.products[kp].findRings ? sides.products[kp].findRings(8) : [];
            for (j = 0; j < rings.length; j++) {
                for (var ap = 0; ap < rings[j].atoms.length; ap++) { pRingSet[rings[j].atoms[ap]] = true; }
            }
        }
        var n = 0;
        for (var rid in result.mapping) {
            if (result.mapping.hasOwnProperty(rid) && rRingSet[+rid] && pRingSet[result.mapping[rid]]) { n++; }
        }
        return n;
    }

    function applyMappingToReaction(reaction, atomMap) {
        // Strip and re-apply: assign a new contiguous mapNumber starting at 1,
        // sorted by reactant atom id for determinism.
        for (var i = 0; i < reaction.atoms.length; i++) { reaction.atoms[i].mapNumber = 0; }
        var rIds = Object.keys(atomMap).map(function(k) { return +k; });
        rIds.sort(function(a, b) { return a - b; });
        var nextN = 1;
        for (var k = 0; k < rIds.length; k++) {
            var rid = rIds[k];
            var pid = atomMap[rid];
            var ra = reaction.getAtom(rid);
            var pa = reaction.getAtom(pid);
            if (ra && pa) {
                ra.mapNumber = nextN;
                pa.mapNumber = nextN;
                nextN += 1;
            }
        }
    }

    // Validate atom-balance: returns { balanced: bool, reactantCounts, productCounts }.
    function checkBalance(sides) {
        var rc = elementCounts(sides.reactants);
        var pc = elementCounts(sides.products);
        return { balanced: elementCountsEqual(rc, pc), reactantCounts: rc, productCounts: pc };
    }

    // Run CIP perception on every component of every side. Best-effort: any
    // exception inside CipStereo.assign is swallowed so a CIP failure on one
    // pathological component never breaks the AAM run as a whole.
    function runCipOnSides(sides) {
        var Cip = (typeof global.CipStereo !== 'undefined') ? global.CipStereo : null;
        if (!Cip || typeof Cip.assign !== 'function') { return false; }
        var k;
        for (k = 0; k < sides.reactants.length; k++) {
            try { Cip.assign(sides.reactants[k]); } catch (e) { /* best-effort */ }
        }
        for (k = 0; k < sides.products.length; k++) {
            try { Cip.assign(sides.products[k]); } catch (e) { /* best-effort */ }
        }
        return true;
    }

    function mapReaction(reaction, options) {
        var opts = shallowCopyOpts(options);
        var t0 = nowMs();
        var timeoutMs = opts.timeoutMs || DEFAULT_TIMEOUT_MS;
        var debug = !!opts.debug;
        var enabledStrategies = (opts.strategies && opts.strategies.length > 0) ? opts.strategies : STRATEGIES;
        var includeStereo = (opts.includeStereo !== false);
        // v1.3.0: bipartite post-pass is now ON by default. Strictly improves
        // total MCS coverage (or is a no-op when greedy is already optimal).
        // Set `useBipartitePostPass: false` to restore v1.2.x default-off
        // behaviour.
        var useBipartite = (opts.useBipartitePostPass !== false);

        if (!reaction || !reaction.atoms || reaction.atoms.length === 0) {
            return {
                mapping: {},
                bondChanges: [],
                score: 0,
                strategy: null,
                timedOut: false,
                warnings: ['empty reaction'],
                strategyResults: []
            };
        }

        // Snapshot pre-existing maps as user constraints (rAtomId -> pAtomId pairs).
        var preMapPairs = [];
        var sides0 = splitReactionSides(reaction);
        var rPreMaps = collectPreMaps(sides0.reactants);
        var pPreMaps = collectPreMaps(sides0.products);
        for (var n in rPreMaps) {
            if (rPreMaps.hasOwnProperty(n) && pPreMaps.hasOwnProperty(n)) {
                var rIds = rPreMaps[n], pIds = pPreMaps[n];
                if (rIds.length === 1 && pIds.length === 1) {
                    preMapPairs.push([rIds[0], pIds[0]]);
                }
            }
        }
        opts._preMaps = preMapPairs;

        // Strip incoming maps before running strategies (we recompute them).
        stripMapNumbers(reaction);

        var balance = checkBalance(sides0);
        var warnings = [];
        if (!balance.balanced) {
            warnings.push('atom-balance mismatch (best-effort mapping)');
        }
        if (sides0.reactants.length === 0 || sides0.products.length === 0) {
            warnings.push('one side empty');
            // Restore pre-maps (if any) onto the reaction so user input is preserved.
            for (var pi = 0; pi < preMapPairs.length; pi++) {
                var ra = reaction.getAtom(preMapPairs[pi][0]);
                var pa = reaction.getAtom(preMapPairs[pi][1]);
                if (ra && pa) { ra.mapNumber = pa.mapNumber = (pi + 1); }
            }
            return {
                mapping: {},
                bondChanges: [],
                score: 0,
                strategy: null,
                timedOut: false,
                warnings: warnings,
                strategyResults: []
            };
        }

        var results = [];
        var timedOut = false;
        for (var i = 0; i < enabledStrategies.length; i++) {
            if ((nowMs() - t0) > timeoutMs) { timedOut = true; break; }
            var s = enabledStrategies[i];
            if (STRATEGIES.indexOf(s) < 0) { continue; }
            var r = runStrategy(reaction, opts, s);

            // (v1.2.0) CIP perception is run AFTER strategies extract per-side
            // sub-molecules so each side has its own correctly-perceived R/S
            // labels. Then re-annotate bond changes so stereoChange events
            // reflect the just-computed cipLabel deltas (e.g. SN2 inversion).
            // CIP failures are best-effort and never abort the mapping.
            if (includeStereo && r && r.sides) {
                runCipOnSides(r.sides);
                var bcOpts = {
                    includeHydrogens: (opts.includeHydrogens !== false),
                    includeStereo: true
                };
                r.bondChanges = annotateBondChanges(r.sides.reactants, r.sides.products, r.mapping, bcOpts);
                r.score = fitness({ reactants: r.sides.reactants, products: r.sides.products }, r.mapping, r.bondChanges);
            }

            // (v1.2.0) Optional bipartite post-pass. Only re-pair if it
            // strictly improves total MCS coverage versus the strategy greedy.
            if (useBipartite && r && r.sides) {
                applyBipartitePostPass(r, opts);
            }

            results.push(r);
            if (debug) {
                // Lightweight debug breadcrumb (no console spam in production)
                r._elapsedMs = nowMs() - t0;
            }
        }

        var winner = pickWinner(results);
        if (winner) {
            applyMappingToReaction(reaction, winner.mapping);
            return {
                mapping: winner.mapping,
                bondChanges: winner.bondChanges,
                score: winner.score,
                strategy: winner.strategy,
                timedOut: timedOut,
                warnings: warnings,
                strategyResults: results,
                bipartiteApplied: !!winner.bipartiteApplied
            };
        }

        return {
            mapping: {},
            bondChanges: [],
            score: 0,
            strategy: null,
            timedOut: timedOut,
            warnings: warnings.concat(['no candidate mapping produced']),
            strategyResults: results
        };
    }

    // -----------------------------------------------------------------------
    // Exports
    // -----------------------------------------------------------------------

    var RDT = {
        mapReaction: mapReaction,
        runMinPairwise: runMinPairwise,
        runMaxPairwise: runMaxPairwise,
        runMixturePairwise: runMixturePairwise,
        runRingPairwise: runRingPairwise,
        annotateBondChanges: annotateBondChanges,
        fitness: fitness,
        STRATEGIES: STRATEGIES.slice(),
        // Internals re-exported for tests
        _splitReactionSides: splitReactionSides,
        _checkBalance: checkBalance,
        _buildSimilarityCache: buildSimilarityCache,
        _pickWinner: pickWinner,
        _applyMappingToReaction: applyMappingToReaction,
        _bipartiteComponentPairing: _bipartiteComponentPairing,
        _runCipOnSides: runCipOnSides,
        version: '1.3.0'
    };

    global.RDT = RDT;

})(typeof window !== 'undefined' ? window : globalThis);
