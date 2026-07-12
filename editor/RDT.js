/**
 * RDT.js - Atom-atom mapping engine for BIME.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Atom-atom mapping (AAM) and bond-change annotation for reactions.
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
 * v1.4.1 additions
 * ----------------
 *   - RDT.deriveSubFragments(result, options) returns the per-MCS-sub-fragment
 *     halo groups for the renderer, where each "sub-fragment" is a maximal
 *     connected sub-graph of mapped reactant atoms whose corresponding mapped
 *     product atoms are also bonded with the same skeleton (preserved-bond
 *     union-find). Bond-change events break sub-fragments, so each one
 *     corresponds to a rigid scaffold piece that survived intact through the
 *     reaction. This matches the classic three-color visualisation in enzyme
 *     mapping diagrams (e.g. blue scaffold + green scaffold + orange
 *     transferred phosphate).
 *
 * v1.4.2 additions
 * ----------------
 *   - Leftover-atom rescue pass after the main bipartite component pairing.
 *     Catches the unbalanced-component case (e.g. benzene + Cl ->
 *     chlorobenzene: 2 reactant components, 1 product component) where the
 *     one-to-one pairing leaves a reactant component unpaired even though its
 *     atoms appear unmapped inside an already-paired product component.
 *     Element + charge + isotope identity required; mapped-neighbour anchor
 *     score breaks ties. Bond-change annotation re-runs after rescue so
 *     newly-mapped atoms participate in `formed` events. componentPairs
 *     descriptor is updated so rescued reactant components show up as
 *     paletteIndex >= 0 with `rescued: true`.
 *
 * v2.0.10 additions (AAM-quality release — all opt-out, defaults additive)
 * -----------------------------------------------------------------------
 *   - options.requireStoichiometricBalance (default true): pre-flight
 *     atom-balance check. Unbalanced reactions return early with
 *     status='unbalanced' + an `elementDelta` map (reactant - product per
 *     element) instead of silently producing a "best-effort" mapping that
 *     hides the imbalance. Set to false to restore the pre-v2.0.10
 *     permissive flow (best-effort mapping + 'atom-balance mismatch'
 *     warning).
 *   - options.mechanisticPlausibility (default true): adds a penalty term
 *     to fitness() that fires when a cleaved or orderChange event sits
 *     ENTIRELY inside a pharmacophore-relevant SMARTS match on the
 *     reactant side (amide, ester, aryl ring, pyranose-like O-heterocycle,
 *     phosphate-like). Each broken match contributes 2.0 to the score.
 *     Pass false to compare scores against pre-v2.0.10 fitness numbers.
 *   - options.topK (default 1): when set > 1, the result includes a
 *     `candidates: [{strategy, mapping, bondChanges, score}, ...]` array
 *     sorted ascending by fitness so reviewers can inspect alternatives.
 *   - Every mapReaction return now carries a `status` field
 *     ('mapped' | 'unbalanced' | 'one-side-empty' | 'empty' | 'no-mapping')
 *     and a `decisiveness` field (1 - clamped gap between winner and
 *     second-best, scale-invariant). Distinct from the existing
 *     coverage-based `confidence` field.
 *
 * Public API
 * ----------
 *   RDT.mapReaction(reaction, options) -> { mapping, bondChanges, score, strategy,
 *                                          status, decisiveness, candidates, ... }
 *   RDT.runMinPairwise(reaction, options)
 *   RDT.runMaxPairwise(reaction, options)
 *   RDT.runMixturePairwise(reaction, options)
 *   RDT.runRingPairwise(reaction, options)
 *   RDT.annotateBondChanges(reactantSide, productSide, atomMap, opts) -> bondChanges[]
 *   RDT.fitness(reaction, atomMap, bondChanges, opts) -> number
 *   RDT._bipartiteComponentPairing(rN, pN, weights) -> [{reactantIdx,productIdx},...]
 *   RDT.deriveSubFragments(result, options) -> [{reactantAtomIds, productAtomIds, paletteIndex}, ...]
 */
(function(global) {
    'use strict';

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    var STRATEGIES = ['MIN', 'MAX', 'MIXTURE', 'RING'];

    // Tie-break order when two strategies post equal fitness:
    //   Lex order chosen to be deterministic (alphabetical).
    var TIEBREAK_ORDER = { MAX: 0, MIN: 1, MIXTURE: 2, RING: 3, MAX_FLEX: 4, MIN_FLEX: 5 };

    var DEFAULT_TIMEOUT_MS = 10000;

    // -----------------------------------------------------------------------
    // v2.0.10 — pharmacophore SMARTS library for mechanistic-plausibility
    // -----------------------------------------------------------------------
    //
    // A small library of substructures whose preservation is a proxy for
    // mechanistic plausibility. A candidate mapping that breaks an amide,
    // ester, aromatic ring, or sugar oxygen heterocycle at the bond-change
    // level is more likely to be wrong than one that leaves these intact —
    // real reactions almost always rearrange a single bond cluster while
    // leaving downstream functional groups untouched.
    //
    // The penalty fires only when a `cleaved` or `orderChange` event sits
    // ENTIRELY INSIDE one of these SMARTS matches on the reactant side
    // (i.e. both bond endpoints participate in the match). It does NOT
    // fire for bonds that merely touch a match at one endpoint — those
    // are normal substituent rearrangements. Each fired event adds
    // PHARMACOPHORE_PENALTY_WEIGHT to the fitness score.
    //
    // Conservative by design: when SmartsMatch is unavailable, when a
    // SMARTS fails to parse, or when bondChanges is empty, the penalty
    // is zero — never negative, never a regression from v2.0.9 fitness.
    var PHARMACOPHORE_SMARTS = [
        '[NX3][CX3](=[OX1])',            // amide        — peptide bonds, anilides
        '[CX3](=[OX1])[OX2H0]',          // ester        — RCOOR'
        'c1ccccc1',                      // aryl ring    — benzene, naphthalene unit
        '[#6]1[OX2][#6][#6][#6][#6]1',   // pyranose-like 6-membered O-heterocycle
        '[PX4](=[OX1])([OX2])[OX2]'      // phosphate-like P(=O)(O)(O)
    ];
    var PHARMACOPHORE_PENALTY_WEIGHT = 2.0;

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

            // Re-normalise chirality to the component's adjacency frame. The
            // component re-adds atoms in sorted-id order, which can permute a
            // stereocentre's neighbour sequence relative to the source. But
            // atom.chirality is stored against the SOURCE getNeighbors() frame
            // (SmilesParser normalises tokens into it, and SmilesWriter reads
            // it back), so without this correction the copied @/@@ token would
            // silently invert on write for any reordered ring/adjacent centre.
            for (j = 0; j < sub.atoms.length; j++) {
                var chAtom = sub.atoms[j];
                if (!chAtom.chirality) { continue; }
                var srcOrder = source.getNeighbors(chAtom.id).filter(function(n) { return idSet[n]; });
                var subOrder = sub.getNeighbors(chAtom.id);
                if (srcOrder.length !== subOrder.length || srcOrder.length < 3) { continue; }
                var cperm = [], cok = true;
                for (var cpi = 0; cpi < srcOrder.length; cpi++) {
                    var cidx = subOrder.indexOf(srcOrder[cpi]);
                    if (cidx < 0) { cok = false; break; }
                    cperm.push(cidx);
                }
                if (!cok) { continue; }
                var cswaps = 0, cseen = {};
                for (var cci = 0; cci < cperm.length; cci++) {
                    if (cseen[cci]) { continue; }
                    var ck = cci, cl = 0;
                    while (!cseen[ck]) { cseen[ck] = true; ck = cperm[ck]; cl++; }
                    cswaps += (cl - 1);
                }
                if (cswaps % 2 === 1) {
                    chAtom.chirality = (chAtom.chirality === '@') ? '@@' : '@';
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

    // v2.0.10: heavy-atom-only counts for the stoichiometric pre-check.
    // Excludes explicit H atoms because the implicit-H count is computed
    // from valence and is not a reliable balance signal across SMILES that
    // mix explicit and implicit representations of the same hydrogen pool
    // (e.g. CC=O.[H][H]>>CCO is heavy-atom-balanced; an H-inclusive count
    // would mark it imbalanced and refuse a perfectly sensible reduction).
    function heavyElementCounts(components) {
        var counts = {};
        for (var i = 0; i < components.length; i++) {
            var atoms = components[i].atoms;
            for (var j = 0; j < atoms.length; j++) {
                var sym = atoms[j].symbol;
                if (sym === 'H') { continue; }
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

    // Per-element difference (a - b) over the union of keys, emitting only
    // the elements whose counts actually differ. Used to report which
    // elements are unbalanced in the stoichiometric pre-check.
    function elementCountDelta(a, b) {
        var delta = {};
        var k;
        for (k in a) {
            if (a.hasOwnProperty(k) && (a[k] || 0) !== (b[k] || 0)) {
                delta[k] = (a[k] || 0) - (b[k] || 0);
            }
        }
        for (k in b) {
            if (b.hasOwnProperty(k) && !delta.hasOwnProperty(k) &&
                (a[k] || 0) !== (b[k] || 0)) {
                delta[k] = (a[k] || 0) - (b[k] || 0);
            }
        }
        return delta;
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

    function buildSimilarityCache(reactants, products, chemOpts, mcsOpts, ringOnly, sharedCache) {
        // v2.0.15: MCS memoization across strategies. The four pairing
        // strategies (MIN/MAX/MIXTURE/RING) differ only in their comparator/
        // filter — the underlying MCS is identical for MIN/MAX/MIXTURE (same
        // cOpts) and a second flavour for RING (ringMatchesRingOnly=true). So
        // there are only TWO distinct caches per mapReaction call, not four.
        // When mapReaction passes a sharedCache object, memoize by flavour and
        // reuse — turning 4x MCS computation into 2x (or 1x when RING is the
        // only or no-op strategy). Result-identical: cache cells hold only
        // numbers + atom-id->atom-id maps (no molecule-object references), and
        // extractComponents preserves atom ids across the per-strategy splits,
        // so a cell computed in one strategy is byte-valid in the next. The
        // single-strategy public entry points (run*Pairwise) pass no
        // sharedCache, so each still builds its own — behaviour unchanged.
        var flavour = (ringOnly ? 'ring' : 'default') +
            ((chemOpts && chemOpts.matchBondOrder === 'any') ? '-flex' : '');
        if (sharedCache && sharedCache[flavour]) { return sharedCache[flavour]; }

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
                    ringSize: countRingMappedAtoms(reactants[i], atomMapping),
                    chainSize: (mcs.size || 0) - countRingMappedAtoms(reactants[i], atomMapping)
                };
            }
        }
        if (sharedCache) { sharedCache[flavour] = cache; }
        return cache;
    }

    // Build the set { atomId: true } of atoms that participate in any ring
    // (size <= 8) across a side. `side` is an array of Molecule components;
    // pass [mol] to scan a single molecule. Centralises the ring-membership
    // scan that fitness / countRingPreserved / countRingMappedAtoms all
    // need (was independently rebuilt in each).
    function buildRingAtomSet(side) {
        var set = {};
        for (var i = 0; i < side.length; i++) {
            var mol = side[i];
            var rings = mol.findRings ? mol.findRings(8) : [];
            for (var r = 0; r < rings.length; r++) {
                var atoms = rings[r].atoms;
                for (var a = 0; a < atoms.length; a++) { set[atoms[a]] = true; }
            }
        }
        return set;
    }

    function countRingMappedAtoms(reactantMol, atomMapping) {
        if (!atomMapping) { return 0; }
        var ringSet = buildRingAtomSet([reactantMol]);
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
            result.mapping, result.bondChanges, options);
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
        // v1.4.1 (bug-fix #7): precompute the reverse mapping pId -> rId once
        // up front so each product bond is checked in O(1) instead of
        // O(|atomMap|). Was O(B_p × M); now O(B_p + M).
        var pBondIds = Object.keys(pBonds);
        pBondIds.sort();
        var pToR = {};
        for (var rkInit in atomMap) {
            if (atomMap.hasOwnProperty(rkInit)) { pToR[atomMap[rkInit]] = rkInit; }
        }
        for (i = 0; i < pBondIds.length; i++) {
            var pb2 = pBonds[pBondIds[i]];
            if (seenProductBonds[pb2.id]) { continue; }
            // Both endpoints must be in the image of atomMap (i.e. mapped
            // from the reactant side) for this to be a formed bond — otherwise
            // it's a bond between unmapped product atoms (unbalanced reaction)
            // and we skip it to avoid spurious events.
            var p1Mapped = (pToR[pb2.atom1] !== undefined);
            var p2Mapped = (pToR[pb2.atom2] !== undefined);
            if (!(p1Mapped && p2Mapped)) { continue; }

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
        // are responsible for CIPStereo.assign() if they want stereo events.
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
    // v2.0.10 — mechanistic-plausibility helpers
    // -----------------------------------------------------------------------

    // For one reactant Molecule, count SMARTS matches per pharmacophore.
    // Returns a map smartsString -> [match, ...] where each match is the
    // SmartsMatch atom-id mapping object (queryAtomId -> targetAtomId).
    // Safe-by-default: returns an empty bucket per SMARTS when SmartsMatch
    // or SmartsParser is unavailable, or when a SMARTS fails to parse.
    function pharmacophoreMatches(mol) {
        var out = {};
        var hasSmarts = (typeof SmartsMatch !== 'undefined' &&
                         typeof SmartsMatch.matchSmarts === 'function');
        for (var i = 0; i < PHARMACOPHORE_SMARTS.length; i++) {
            var smarts = PHARMACOPHORE_SMARTS[i];
            out[smarts] = [];
            if (!hasSmarts) { continue; }
            try {
                var mm = SmartsMatch.matchSmarts(mol, smarts);
                if (Array.isArray(mm)) { out[smarts] = mm; }
            } catch (e) { /* keep empty on parse / runtime error */ }
        }
        return out;
    }

    // Compute the mechanistic-plausibility penalty for a candidate mapping.
    // For every reactant pharmacophore match, check whether any cleaved or
    // orderChange bondChange has BOTH its reactant endpoints inside the
    // match — that is, the bond change happens INSIDE a pharmacophore.
    // Each such match contributes PHARMACOPHORE_PENALTY_WEIGHT once
    // (counting per match, not per event, so a tri-substituted ring
    // doesn't multi-penalise).
    //
    // Returns 0 when bondChanges is empty, when SmartsMatch is absent, or
    // when no event lies inside any match.
    function mechanisticPlausibilityPenalty(rSide, bondChanges) {
        if (!bondChanges || bondChanges.length === 0) { return 0; }
        if (!rSide || rSide.length === 0) { return 0; }
        var broken = 0;
        for (var ri = 0; ri < rSide.length; ri++) {
            var rMol = rSide[ri];
            var hits = pharmacophoreMatches(rMol);
            // Build a reactant-bond lookup once per molecule so we can
            // resolve each event's bond endpoints in O(1).
            var bondById = {};
            if (rMol.bonds) {
                for (var bi = 0; bi < rMol.bonds.length; bi++) {
                    bondById[rMol.bonds[bi].id] = rMol.bonds[bi];
                }
            }
            for (var key in hits) {
                if (!hits.hasOwnProperty(key)) { continue; }
                var matches = hits[key];
                for (var mi = 0; mi < matches.length; mi++) {
                    var atomSet = {};
                    var mapping = matches[mi];
                    for (var qId in mapping) {
                        if (mapping.hasOwnProperty(qId)) {
                            atomSet[mapping[qId]] = true;
                        }
                    }
                    var matchBroken = false;
                    for (var ei = 0; ei < bondChanges.length; ei++) {
                        var ev = bondChanges[ei];
                        if (ev.type !== 'cleaved' && ev.type !== 'orderChange') {
                            continue;
                        }
                        var rb = bondById[ev.reactantBond];
                        if (!rb) { continue; }
                        if (atomSet[rb.atom1] && atomSet[rb.atom2]) {
                            matchBroken = true;
                            break;
                        }
                    }
                    if (matchBroken) { broken += 1; }
                }
            }
        }
        return broken * PHARMACOPHORE_PENALTY_WEIGHT;
    }

    // -----------------------------------------------------------------------
    // Fitness function
    // -----------------------------------------------------------------------

    // Compute the fitness score for a candidate. Lower is better.
    //
    // v2.0.10: accepts an optional `opts` argument. When
    // `opts.mechanisticPlausibility !== false` (default ON), the score
    // includes a penalty term that fires when cleaved/orderChange events
    // occur INSIDE a pharmacophore-relevant SMARTS match on the reactant
    // side (amide, ester, aryl ring, pyranose-like O-heterocycle,
    // phosphate-like). Pass `{ mechanisticPlausibility: false }` to
    // compare scores against the pre-v2.0.10 fitness function.
    function fitness(reactionOrSides, atomMap, bondChanges, opts) {
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
        var rRingSet = buildRingAtomSet(rSide);
        var pRingSet = buildRingAtomSet(pSide);
        var k;
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

        var base = changeCount - 0.5 * ringPreserved + 1.0 * aromaticMismatches;
        var useMP = !opts || opts.mechanisticPlausibility !== false;
        var mp = useMP ? mechanisticPlausibilityPenalty(rSide, bondChanges) : 0;

        // v2.0.46: chemistry-aware fitness terms. Promotes the v2.0.41
        // diagnostic metrics into selection-time penalties:
        //   - 0.3 per WHOLE aromatic ring broken across the mapping
        //     (stronger than the existing per-atom aromaticMismatches
        //     penalty because the per-atom flag misses scatters that
        //     keep individual atoms aromatic but break the ring)
        //   - 0.3 per stereo centre dropped (mapped pair where one
        //     side carries a CIP / chirality descriptor but the
        //     partner does not)
        // Weights deliberately small (smaller than the existing 1.0
        // aromaticMismatches and changeCount terms) so they tie-break
        // rather than flip clear winners.
        //
        // v2.0.57: DEFAULT ON (was opt-in via opts.chemistryAware===true
        // in v2.0.46-v2.0.56). Golden-corpus run with chemistryAware=true
        // passes 138/138, identical to the chemistryAware=false baseline,
        // so this is a safe default-flip. Callers can still opt out with
        // `opts.chemistryAware === false` for byte-identical v2.0.46
        // mapping selection on borderline tie-breaks.
        var chem = 0;
        if (!opts || opts.chemistryAware !== false) {
            chem = chemistryPenalty(rSide, pSide, atomMap, rAtoms, pAtoms);
        }
        return base + mp + chem;
    }

    // v2.0.46: count whole aromatic rings broken + stereo centres dropped
    // across the mapping. Returns a small positive penalty proportional to
    // the chemistry-correctness deficit. Read-only — does NOT mutate any
    // input. Used by `fitness` only when opts.chemistryAware = true.
    function chemistryPenalty(rSide, pSide, atomMap, rAtomIdx, pAtomIdx) {
        var penalty = 0;
        // ----- aromatic-ring conservation -----
        // Use the v2.0.41 helper but inline its core to avoid a
        // result-object allocation in the hot loop.
        var rRings = collectAromaticRings(rSide);
        var pRings = collectAromaticRings(pSide);
        for (var ri = 0; ri < rRings.length; ri++) {
            var ring = rRings[ri];
            var mappedSet = {};
            var fullyMapped = true;
            for (var ai = 0; ai < ring.length; ai++) {
                var pid = atomMap[ring[ai]];
                if (typeof pid === 'undefined') { fullyMapped = false; break; }
                mappedSet[pid] = true;
            }
            if (!fullyMapped) { penalty += 0.3; continue; }
            var foundMatch = false;
            for (var rj = 0; rj < pRings.length; rj++) {
                var pr = pRings[rj];
                if (pr.length !== ring.length) { continue; }
                var match = true;
                for (var pi = 0; pi < pr.length; pi++) {
                    if (!mappedSet[pr[pi]]) { match = false; break; }
                }
                if (match) { foundMatch = true; break; }
            }
            if (!foundMatch) { penalty += 0.3; }
        }
        // ----- stereo-centre conservation -----
        // Per-atom dropped-stereo check. Avoids requiring R == R since
        // valid chemistry (SN2 Walden inversion) flips R↔S.
        for (var aid in atomMap) {
            if (!atomMap.hasOwnProperty(aid)) { continue; }
            var ra = rAtomIdx[+aid];
            if (!ra || !hasStereoDescriptor(ra)) { continue; }
            var pa = pAtomIdx[atomMap[aid]];
            if (!pa || !hasStereoDescriptor(pa)) { penalty += 0.3; }
        }
        return penalty;
    }

    // -----------------------------------------------------------------------
    // Strategy runners
    // -----------------------------------------------------------------------

    // Build a candidate {atomMap, pairs} from a comparator-based strategy.
    // `sharedCache` (optional) memoizes MCS results across strategies within
    // a single mapReaction call; the single-strategy public entry points pass
    // none.
    function runStrategy(reaction, options, strategyName, sharedCache) {
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

        // '<STRATEGY>_FLEX' runs the base strategy but lets the MCS span a
        // reaction-centre bond-ORDER change (single<->double on a bond that is
        // preserved between the same atoms). It is added as an EXTRA candidate
        // and only wins on a strictly better score (fewer bond changes), so it
        // corrects ethanol oxidation written O-first without disturbing
        // reactions where the strict mapping (e.g. a genuinely formed C-C bond
        // in a Claisen condensation) is the better-scoring answer.
        var flexBond = /_FLEX$/.test(strategyName);
        var baseStrategy = flexBond ? strategyName.slice(0, -5) : strategyName;
        var chemOpts = options.chemOpts || {};
        if (flexBond) {
            chemOpts = shallowCopyOpts(chemOpts);
            chemOpts.matchBondType = false;
            chemOpts.matchBondOrder = 'any';
        }
        var mcsOpts = options.mcsOpts || {};
        var ringOnly = (baseStrategy === 'RING');
        var cache;
        try {
            cache = buildSimilarityCache(reactants, products, chemOpts, mcsOpts, ringOnly, sharedCache);
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
        switch (baseStrategy) {
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
        var score = fitness({ reactants: sides.reactants, products: sides.products }, atomMap, bondChanges, options);

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
        var rRingSet = buildRingAtomSet(sides.reactants);
        var pRingSet = buildRingAtomSet(sides.products);
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

    // -----------------------------------------------------------------------
    // v1.4.2: leftover-atom rescue.
    //
    // The bipartite component-pairing in mapReaction is one-to-one: each
    // reactant component pairs with at most one product component. When the
    // reaction is unbalanced in component count (e.g. benzene + Cl ->
    // chlorobenzene: 2 reactant components, 1 product component), one
    // reactant component ends up unpaired and its atoms get no map number,
    // even when the same element appears unmapped inside an already-paired
    // product component.
    //
    // This pass walks the unpaired reactant components, finds their atoms,
    // and matches them to unmapped product atoms by:
    //   1) element symbol identity (required),
    //   2) charge identity (required),
    //   3) isotope identity (required),
    //   4) mapped-neighbour anchor score (preferred): how many already-
    //      mapped reactant neighbours of r have a product partner that is
    //      also a neighbour of p, plus a 0.5-weighted bonus for each
    //      mapped neighbour of p (joined-in semantic — favour atoms that
    //      ended up bonded to the rest of the mapping).
    //
    // The match is greedy: pick the highest-scoring (rId, pId) pair, fix it,
    // remove both from candidate pools, repeat until no compatible pair
    // remains. Element + charge + isotope guards make this safe — a Cl will
    // never be mapped to a non-Cl, etc.
    //
    // Mathematical notes (v1.4.2 panel-review):
    //   - The 0.5 tie-break weight is in (0,1); since the primary
    //     anchor-overlap score is in N, floor(s + 0.5·k) = s, so the
    //     0.5 term acts strictly as a deterministic tie-breaker between
    //     primary-equal candidates.
    //   - Greedy max-weight bipartite matching is a 1/2-approximation in
    //     the worst case (Avis 1983). For chemistry inputs with element
    //     + charge + isotope class partitioning, the matching is usually
    //     vertex-disjoint per class and effectively optimal in practice.
    //   - Monotonicity property: rescue only ADDS keys to `mapping` and
    //     never removes them, so any cost function that is non-increasing
    //     in mapping extension (incl. RDT.fitness) cannot get worse.
    //     `winner.scoreAfterRescue` is therefore <= `winner.score` (lower
    //     is better in this codebase).
    //
    // Returns true when at least one new mapping was added, false otherwise.
    // The caller is expected to re-run annotateBondChanges + re-derive
    // componentPairs / confidence so downstream consumers see the extended
    // mapping consistently.
    // -----------------------------------------------------------------------
    function rescueLeftoverAtoms(winner) {
        if (!winner || !winner.mapping || !winner.sides) { return false; }
        var mapping = winner.mapping;
        var rSide = winner.sides.reactants || [];
        var pSide = winner.sides.products || [];

        var pairedR = {};
        var pairs = winner.pairs || [];
        for (var pi = 0; pi < pairs.length; pi++) { pairedR[pairs[pi].rIdx] = true; }

        // Build adjacency first — we need it for the anchor-eligibility check.
        var rAdj = {};
        for (var i = 0; i < rSide.length; i++) {
            var bonds = rSide[i].bonds || [];
            for (var j = 0; j < bonds.length; j++) {
                var b = bonds[j];
                if (!rAdj[b.atom1]) { rAdj[b.atom1] = []; }
                if (!rAdj[b.atom2]) { rAdj[b.atom2] = []; }
                rAdj[b.atom1].push(b.atom2);
                rAdj[b.atom2].push(b.atom1);
            }
        }
        var pAdj = {};
        for (i = 0; i < pSide.length; i++) {
            var pBonds = pSide[i].bonds || [];
            for (j = 0; j < pBonds.length; j++) {
                var pb = pBonds[j];
                if (!pAdj[pb.atom1]) { pAdj[pb.atom1] = []; }
                if (!pAdj[pb.atom2]) { pAdj[pb.atom2] = []; }
                pAdj[pb.atom1].push(pb.atom2);
                pAdj[pb.atom2].push(pb.atom1);
            }
        }

        // Reactant rescue candidates:
        //   1) atoms in UNPAIRED components (spectator molecules) — always
        //      eligible (the whole component had no MCS partner), OR
        //   2) atoms in PAIRED components that have at least one already-
        //      mapped neighbour ("anchored" partial-component leftover —
        //      e.g. the methyl C of CH3-I when only I got mapped).
        // The anchor rule prevents arbitrary element-matching across two
        // disjoint component pairings (which would otherwise create spurious
        // mappings on, say, CC.OO.NN>>CO.CO).
        var rCands = [];
        for (i = 0; i < rSide.length; i++) {
            var atoms = rSide[i].atoms || [];
            for (var j = 0; j < atoms.length; j++) {
                var a = atoms[j];
                if (mapping[a.id] !== undefined) { continue; }
                if (!pairedR[i]) {
                    // Unpaired component — atom always eligible.
                    rCands.push(a);
                    continue;
                }
                // Paired component — require an already-mapped neighbour.
                var nbrs = rAdj[a.id] || [];
                var anchored = false;
                for (var k = 0; k < nbrs.length; k++) {
                    if (mapping[nbrs[k]] !== undefined) { anchored = true; break; }
                }
                if (anchored) { rCands.push(a); }
            }
        }
        if (rCands.length === 0) { return false; }

        // Product atoms not yet in the image of mapping.
        var imagePids = {};
        for (var kk in mapping) { if (mapping.hasOwnProperty(kk)) { imagePids[mapping[kk]] = true; } }
        var pCands = [];
        for (var pCompIdx = 0; pCompIdx < pSide.length; pCompIdx++) {
            var pAtoms = pSide[pCompIdx].atoms || [];
            for (var pj = 0; pj < pAtoms.length; pj++) {
                if (!imagePids[pAtoms[pj].id]) { pCands.push(pAtoms[pj]); }
            }
        }
        if (pCands.length === 0) { return false; }

        function compatible(rA, pA) {
            if (!rA || !pA) { return false; }
            if (rA.symbol !== pA.symbol) { return false; }
            if ((rA.isotope || 0) !== (pA.isotope || 0)) { return false; }
            // v2.0.13: charge identity is required for multivalent
            // (heavy-skeleton) atoms — relaxing it there could create
            // wrong mappings. But a monatomic ion or a terminal leaving-
            // group atom — degree <= 1 on BOTH sides — is chemically
            // expected to change formal charge across a reaction: e.g.
            // Finkelstein, the SAME chlorine goes from neutral C-Cl to
            // Cl(-) anion, and the SAME iodine goes from I(-) anion to
            // covalent C-I. Allow the charge to differ in that bounded
            // case only; the anchor score below still steers the pairing.
            var rDeg = (rAdj[rA.id] || []).length;
            var pDeg = (pAdj[pA.id] || []).length;
            var bothTerminal = (rDeg <= 1 && pDeg <= 1);
            if (!bothTerminal && (rA.charge || 0) !== (pA.charge || 0)) {
                return false;
            }
            return true;
        }
        function score(rA, pA) {
            // Anchor score: count rA's mapped neighbours whose product partner
            // is a neighbour of pA. Adds half-credit for any mapped neighbour
            // of pA (joined-in bonus) so a leftover atom that bonds to the
            // rest of the mapping wins ties over a free spectator placement.
            var nb1 = rAdj[rA.id] || [];
            var nb2 = pAdj[pA.id] || [];
            var pSet = {};
            for (var k = 0; k < nb2.length; k++) { pSet[nb2[k]] = true; }
            var sc = 0;
            for (k = 0; k < nb1.length; k++) {
                var partner = mapping[nb1[k]];
                if (partner !== undefined && pSet[partner]) { sc++; }
            }
            var pMappedNbr = 0;
            for (k = 0; k < nb2.length; k++) {
                if (imagePids[nb2[k]]) { pMappedNbr++; }
            }
            return sc + 0.5 * pMappedNbr;
        }

        var rescued = false;
        var safety = rCands.length * pCands.length + 1;
        while (safety-- > 0) {
            var bestR = -1, bestP = -1, bestS = -1;
            for (i = 0; i < rCands.length; i++) {
                if (rCands[i] === null) { continue; }
                for (j = 0; j < pCands.length; j++) {
                    if (pCands[j] === null) { continue; }
                    if (!compatible(rCands[i], pCands[j])) { continue; }
                    var s = score(rCands[i], pCands[j]);
                    // Tie-break: prefer lower atom-id pair for determinism.
                    if (s > bestS ||
                        (s === bestS && (bestR < 0 ||
                            rCands[i].id < rCands[bestR].id ||
                            (rCands[i].id === rCands[bestR].id && pCands[j].id < pCands[bestP].id)))) {
                        bestS = s; bestR = i; bestP = j;
                    }
                }
            }
            if (bestR < 0 || bestP < 0) { break; }
            mapping[rCands[bestR].id] = pCands[bestP].id;
            imagePids[pCands[bestP].id] = true;
            rescued = true;
            rCands[bestR] = null;
            pCands[bestP] = null;
        }
        return rescued;
    }

    // -----------------------------------------------------------------------
    // v1.4.2: after a leftover-atom rescue, update the per-component-pair
    // descriptor so each previously-unpaired reactant component shows up as
    // joined to the product component that absorbed its atoms.
    //
    // The descriptor's `paletteIndex` for rescued components is appended at
    // the end of the existing palette indices so renderer halos remain
    // distinguishable from the "main" pairings.
    // -----------------------------------------------------------------------
    function enrichComponentPairsWithRescue(winner, componentPairs) {
        if (!winner || !winner.sides || !winner.mapping) { return componentPairs; }
        var rSide = winner.sides.reactants || [];
        var pSide = winner.sides.products || [];
        var mapping = winner.mapping;

        // Build atom-id -> compIdx maps for quick reverse lookup.
        var rAtomToComp = {};
        for (var i = 0; i < rSide.length; i++) {
            var aR = rSide[i].atoms || [];
            for (var jR = 0; jR < aR.length; jR++) { rAtomToComp[aR[jR].id] = i; }
        }
        var pAtomToComp = {};
        for (i = 0; i < pSide.length; i++) {
            var aP = pSide[i].atoms || [];
            for (var jP = 0; jP < aP.length; jP++) { pAtomToComp[aP[jP].id] = i; }
        }

        // Find next available paletteIndex.
        var maxIdx = -1;
        for (i = 0; i < componentPairs.length; i++) {
            if (componentPairs[i].paletteIndex > maxIdx) { maxIdx = componentPairs[i].paletteIndex; }
        }

        // For each rescued mapping (rId -> pId), figure out the reactant /
        // product components and update the corresponding entry in
        // componentPairs (or insert one if missing).
        var updated = componentPairs.slice();
        var seen = {};
        for (var k in mapping) {
            if (!mapping.hasOwnProperty(k)) { continue; }
            var rCompIdx = rAtomToComp[+k];
            var pCompIdx = pAtomToComp[+mapping[k]];
            if (rCompIdx === undefined || pCompIdx === undefined) { continue; }
            var key = rCompIdx + '|' + pCompIdx;
            if (seen[key]) { continue; }
            seen[key] = true;
            // Is this pair already in the descriptor?
            var found = false;
            for (i = 0; i < updated.length; i++) {
                if (updated[i].reactantCompIdx === rCompIdx &&
                    updated[i].productCompIdx === pCompIdx) {
                    found = true;
                    if (updated[i].paletteIndex < 0) {
                        updated[i].paletteIndex = ++maxIdx;
                    }
                    break;
                }
            }
            if (!found) {
                // Replace the orphan -1 entry for this reactant comp (if any)
                // with a paired entry. Otherwise append a fresh entry.
                var idxOrphanR = -1;
                for (i = 0; i < updated.length; i++) {
                    if (updated[i].reactantCompIdx === rCompIdx &&
                        updated[i].paletteIndex === -1) {
                        idxOrphanR = i;
                        break;
                    }
                }
                if (idxOrphanR >= 0) {
                    updated[idxOrphanR] = {
                        reactantCompIdx: rCompIdx,
                        productCompIdx: pCompIdx,
                        mcsSize: 1,
                        paletteIndex: ++maxIdx,
                        rescued: true
                    };
                } else {
                    updated.push({
                        reactantCompIdx: rCompIdx,
                        productCompIdx: pCompIdx,
                        mcsSize: 1,
                        paletteIndex: ++maxIdx,
                        rescued: true
                    });
                }
            }
        }
        return updated;
    }

    // Run CIP perception on every component of every side. Best-effort: any
    // exception inside CIPStereo.assign is swallowed so a CIP failure on one
    // pathological component never breaks the AAM run as a whole.
    function runCipOnSides(sides) {
        var Cip = (typeof global.CIPStereo !== 'undefined') ? global.CIPStereo : null;
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

    // v2.0.11: shared empty-result shape used by every NON-mapped exit
    // path of mapReaction. Centralises the field surface so v2.0.10's
    // status / decisiveness / candidates / elementDelta additions land
    // in one place instead of five hand-maintained copies — the previous
    // four early-return blocks had drifted by 3-5 line variations and
    // were missing fields each new release. Each caller passes the
    // sides snapshot, the status string, and an optional `extras` bag
    // overriding any of: warnings (defaults to []), strategyResults
    // (defaults to []), elementDelta (set only on 'unbalanced'),
    // score (defaults to 0; the 'unbalanced' path sets Infinity),
    // timedOut (defaults to false).
    function _emptyResult(sides, status, extras) {
        extras = extras || {};
        var out = {
            status: status,
            mapping: {},
            bondChanges: [],
            score: (typeof extras.score === 'number') ? extras.score : 0,
            strategy: null,
            timedOut: !!extras.timedOut,
            warnings: extras.warnings || [],
            strategyResults: extras.strategyResults || [],
            componentPairs: [],
            confidence: 0,
            decisiveness: 0,
            candidates: [],
            sides: sides || { reactants: [], products: [] }
        };
        // Only set elementDelta when extras supplied it (the 'unbalanced'
        // path). Setting `out.elementDelta = undefined` unconditionally
        // would shift `'elementDelta' in result` semantics across the
        // other three non-mapped branches — wire-compatible (JSON-stringify
        // drops undefined) but visible to in-process `in`-checks.
        if (extras.elementDelta) { out.elementDelta = extras.elementDelta; }
        return out;
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
            // v1.4.1 (bug-fix #27): always return a `sides` field so
            // downstream consumers (deriveSubFragments, MolEditor) can call
            // result.sides.reactants without a defensive null check.
            return _emptyResult(null, 'empty', { warnings: ['empty reaction'] });
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
        // v2.0.10: stoichiometric pre-check. When requireStoichiometricBalance
        // (default true) is enabled and the reaction is unbalanced on
        // HEAVY ATOMS, abort with status='unbalanced' and an explicit
        // per-element delta rather than silently producing a best-effort
        // mapping that hides the imbalance. Heavy-atom-only (excluding H)
        // because the implicit-H count is a valence-derived signal that
        // can be unreliable across SMILES strings that mix explicit and
        // implicit H representations (e.g. CC=O.[H][H]>>CCO is heavy-atom-
        // balanced and should map fine). Set requireStoichiometricBalance:
        // false to restore pre-v2.0.10 permissive behaviour (best-effort
        // mapping with an 'atom-balance mismatch' warning).
        var requireBalance = (opts.requireStoichiometricBalance !== false);
        var heavyR = heavyElementCounts(sides0.reactants);
        var heavyP = heavyElementCounts(sides0.products);
        var heavyBalanced = elementCountsEqual(heavyR, heavyP);
        if (!heavyBalanced && requireBalance) {
            var elementDelta = elementCountDelta(heavyR, heavyP);
            var deltaParts = [];
            for (var dk in elementDelta) {
                if (elementDelta.hasOwnProperty(dk)) {
                    var v = elementDelta[dk];
                    deltaParts.push(dk + (v > 0 ? '+' : '') + v);
                }
            }
            // Restore any pre-existing maps the user authored so we don't
            // erase their intent when refusing to map.
            for (var pi0 = 0; pi0 < preMapPairs.length; pi0++) {
                var ra0 = reaction.getAtom(preMapPairs[pi0][0]);
                var pa0 = reaction.getAtom(preMapPairs[pi0][1]);
                if (ra0 && pa0) { ra0.mapNumber = pa0.mapNumber = (pi0 + 1); }
            }
            return _emptyResult(sides0, 'unbalanced', {
                score: Infinity,
                elementDelta: elementDelta,
                warnings: ['unbalanced reaction (element delta reactant-product: ' +
                    deltaParts.join(', ') +
                    '); set options.requireStoichiometricBalance=false to map anyway']
            });
        }
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
            // v1.4.1 (bug-fix #27): always return `sides` so downstream
            // sub-fragment / pair derivation can read it unconditionally.
            return _emptyResult(sides0, 'one-side-empty', { warnings: warnings });
        }

        var results = [];
        var timedOut = false;
        // v2.0.15: per-call MCS memo shared across the strategy loop. NOT
        // global — it must not leak across mapReaction calls (atom ids are
        // only unique within one reaction). MIN/MAX/MIXTURE hit the same
        // 'default' flavour; RING builds 'ring' once. Pass
        // options.shareMcsCache=false to disable the memo and rebuild the
        // MCS per strategy (the pre-v2.0.15 path) — an escape hatch that
        // also lets tests assert cache-on == cache-off result-identity.
        var sharedMcsCache = (opts.shareMcsCache === false) ? null : {};
        for (var i = 0; i < enabledStrategies.length; i++) {
            if ((nowMs() - t0) > timeoutMs) { timedOut = true; break; }
            var s = enabledStrategies[i];
            if (STRATEGIES.indexOf(s) < 0) { continue; }
            var r = runStrategy(reaction, opts, s, sharedMcsCache);

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
                r.score = fitness({ reactants: r.sides.reactants, products: r.sides.products }, r.mapping, r.bondChanges, opts);
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

        // Extra bond-order-flexible candidates. These let the MCS span a
        // reaction-centre bond-ORDER change (a C-O single that becomes a C=O
        // double is the same skeletal bond) and compete purely on score, so
        // they only replace the winner when they produce FEWER bond changes —
        // fixing e.g. ethanol oxidation written O-first (an order change instead
        // of a spurious cleaved+formed pair) while leaving reactions where the
        // strict mapping scores better (a genuinely formed C-C bond) untouched.
        // Skipped when the caller pins an explicit strategy list, so
        // strategy-count assertions are unaffected.
        if (!(opts.strategies && opts.strategies.length > 0)) {
            var flexCache = {};
            var flexList = ['MIN_FLEX', 'MAX_FLEX'];
            for (var fi = 0; fi < flexList.length; fi++) {
                // Run UNGATED — no wall-clock check. The flex candidate is what
                // makes the order-change mapping reachable, so gating it on
                // elapsed time made the RESULT budget-dependent (a cold or slow
                // machine skipped it and returned the wrong ethanol-oxidation
                // mapping). Each flex strategy is still bounded by the MCS's own
                // internal node/time budget, so total work stays bounded.
                var fr = runStrategy(reaction, opts, flexList[fi], flexCache);
                if (includeStereo && fr && fr.sides) {
                    runCipOnSides(fr.sides);
                    fr.bondChanges = annotateBondChanges(fr.sides.reactants, fr.sides.products, fr.mapping, {
                        includeHydrogens: (opts.includeHydrogens !== false),
                        includeStereo: true
                    });
                    fr.score = fitness({ reactants: fr.sides.reactants, products: fr.sides.products }, fr.mapping, fr.bondChanges, opts);
                }
                if (useBipartite && fr && fr.sides) { applyBipartitePostPass(fr, opts); }
                results.push(fr);
            }
        }

        var winner = pickWinner(results);
        if (winner) {
            // v1.4.2: leftover-atom rescue. Catches the unbalanced-component
            // case (e.g. benzene + Cl -> chlorobenzene) where the bipartite
            // pairing leaves a reactant component unpaired even though its
            // atoms appear unmapped inside an already-paired product
            // component. Element + charge + isotope identity is required;
            // mapped-neighbour anchor score breaks ties. Pass
            // `useLeftoverRescue: false` to opt out (full v1.4.1 compat).
            var rescueEnabled = (opts.useLeftoverRescue !== false);
            var rescued = rescueEnabled ? rescueLeftoverAtoms(winner) : false;
            if (rescued && winner.sides) {
                // Re-annotate bond-changes so c4-Cl-style FORMED events pick
                // up the newly-mapped atoms. We preserve `winner.score` as
                // the SELECTION-time fitness (what pickWinner saw); the
                // rescue-extended fitness is exposed as `scoreAfterRescue`
                // so consumers can read either consistently. Strategy
                // selection itself does not re-run.
                var rcOpts = {
                    includeHydrogens: (opts.includeHydrogens !== false),
                    includeStereo: (opts.includeStereo !== false)
                };
                winner.bondChanges = annotateBondChanges(winner.sides.reactants, winner.sides.products, winner.mapping, rcOpts);
                winner.scoreAfterRescue = fitness({ reactants: winner.sides.reactants, products: winner.sides.products }, winner.mapping, winner.bondChanges, opts);
            }
            applyMappingToReaction(reaction, winner.mapping);
            // v1.4.0: precompute the per-component-pair descriptor and a
            // normalised confidence so MolEditor / the renderer halo overlay
            // can read them directly off the result.
            var componentPairs = deriveComponentPairList(winner);
            // v1.4.2: enrich componentPairs so rescued components show up as
            // paired (was previously stuck on paletteIndex = -1). Each rescued
            // reactant component is attributed to whichever product component
            // received the most atoms from it; the entry mirrors the v1.4.0
            // shape so renderer / consumer code is unchanged.
            if (rescued) {
                componentPairs = enrichComponentPairsWithRescue(winner, componentPairs);
            }
            var confidence = deriveConfidence(winner);
            // v2.0.10: top-K candidates + decisiveness. `candidates` is
            // the strategy results sorted ascending by fitness, truncated
            // to opts.topK (default 1 — preserves pre-v2.0.10 contract
            // where only the winner is surfaced). Each entry carries the
            // minimal shape callers actually need to inspect alternatives.
            // `decisiveness` is the gap between the winner and second-
            // best, normalised so a perfect winner (winner=0, second>0)
            // returns 1.0 and a dead tie returns 0. Distinct from the
            // existing coverage-based `confidence`: a high-coverage
            // mapping can still be a low-decisiveness pick when a second
            // strategy ties it. Renderers and reviewers should look at
            // BOTH numbers.
            var topK = (typeof opts.topK === 'number' && opts.topK > 0)
                ? opts.topK : 1;
            var candidates = buildCandidateList(results, topK);
            var decisiveness = deriveDecisiveness(results, winner);
            var mappedResult = {
                status: 'mapped',
                mapping: winner.mapping,
                bondChanges: winner.bondChanges,
                score: winner.score,
                strategy: winner.strategy,
                timedOut: timedOut,
                warnings: warnings,
                strategyResults: results,
                bipartiteApplied: !!winner.bipartiteApplied,
                componentPairs: componentPairs,
                confidence: confidence,
                decisiveness: decisiveness,
                // v2.0.41: AAM quality diagnostic (cheap, ~1ms). Includes
                // categorical grade A/B/C/F, aromatic-ring conservation,
                // stereo-centre conservation, mapped/unmapped counts. Opt
                // out with options.includeDiagnostic = false.
                quality: null,         // populated below if includeDiagnostic
                diagnostic: null,      // populated below if includeDiagnostic
                candidates: candidates,
                // v1.4.0: re-export the winner's per-pair cells and sides so
                // RDT.deriveComponentPairs (atom-id lookup) and tests can
                // operate on the top-level result without rummaging through
                // strategyResults.
                pairs: winner.pairs,
                sides: winner.sides
            };
            // v2.0.41: populate the AAM quality diagnostic by default.
            // Caller can opt out via options.includeDiagnostic = false to
            // save the ~1 ms compute when running mass benchmarks.
            if (opts.includeDiagnostic !== false) {
                mappedResult.diagnostic = deriveDiagnostic(mappedResult);
                mappedResult.quality = mappedResult.diagnostic.quality;
            }
            // v2.0.12: post-classify the reaction into a heuristic family
            // label (substitution, addition, elimination, ...). Default
            // ON; opt out with options.detectReactionClass=false. Purely
            // metadata — the mapping itself is unchanged.
            if (opts.detectReactionClass !== false) {
                var cls = detectReactionClass(mappedResult);
                mappedResult.reactionClass = cls.reactionClass;
                mappedResult.classReason = cls.classReason;
            }
            return mappedResult;
        }

        // v1.4.1 (bug-fix #27): always return `sides` (sides0 captured above).
        return _emptyResult(sides0, 'no-mapping', {
            timedOut: timedOut,
            warnings: warnings.concat(['no candidate mapping produced']),
            strategyResults: results
        });
    }

    // v2.0.10: collect the top-K finite-score results into a compact list,
    // sorted ascending by score (lower is better). When two scores tie,
    // strategy order TIEBREAK_ORDER decides — same convention as
    // pickWinner. Each entry exposes the minimum information a reviewer
    // needs to reason about an alternative mapping without re-running
    // mapReaction.
    function buildCandidateList(results, topK) {
        var sorted = [];
        for (var i = 0; i < results.length; i++) {
            var r = results[i];
            if (!r || !isFinite(r.score)) { continue; }
            sorted.push(r);
        }
        sorted.sort(function (a, b) {
            if (a.score !== b.score) { return a.score - b.score; }
            var ta = (TIEBREAK_ORDER[a.strategy] != null) ?
                TIEBREAK_ORDER[a.strategy] : 999;
            var tb = (TIEBREAK_ORDER[b.strategy] != null) ?
                TIEBREAK_ORDER[b.strategy] : 999;
            return ta - tb;
        });
        var out = [];
        var limit = Math.min(topK, sorted.length);
        for (var k = 0; k < limit; k++) {
            var s = sorted[k];
            out.push({
                strategy: s.strategy,
                mapping: s.mapping,
                bondChanges: s.bondChanges,
                score: s.score
            });
        }
        return out;
    }

    // v2.0.10: decisiveness — how much the winner beat the second-best.
    // Returns 1.0 when (a) no second-best exists (single-strategy result —
    // distinct from "won by a large margin"; the caller may want to look
    // at strategyResults.length to disambiguate) or (b) the winner is a
    // perfect zero-fitness mapping with a positive-fitness second. Returns
    // 0 on dead ties. Otherwise: (second - winner) / (1 + |second|),
    // clamped to [0, 1]. The (1 + |second|) denominator makes the score
    // scale-stable: a 5-vs-10 win gives 5/11 ≈ 0.45, a 50-vs-100 win gives
    // 50/101 ≈ 0.495 — both report "decisive" without exploding when
    // second.score approaches 0. (Not strictly scale-invariant — that
    // would require dropping the +1 — but the +1 is what guards against
    // small-score blow-up, and the trade-off is the right one for the
    // typical fitness magnitudes RDT produces.)
    function deriveDecisiveness(results, winner) {
        if (!winner || !isFinite(winner.score)) { return 0; }
        var second = null;
        for (var i = 0; i < results.length; i++) {
            var r = results[i];
            if (!r || r === winner || !isFinite(r.score)) { continue; }
            if (second == null || r.score < second.score) { second = r; }
        }
        if (!second) { return 1; }
        var gap = second.score - winner.score;
        if (gap <= 0) { return 0; }
        var denom = 1 + Math.abs(second.score);
        var d = gap / denom;
        if (d > 1) { d = 1; }
        if (d < 0) { d = 0; }
        return d;
    }

    // -----------------------------------------------------------------------
    // v2.0.12 — reaction-class detection
    //
    // Heuristic post-classification of a successful mapping. Looks at the
    // bondChange profile + heavy-atom delta + ring delta and assigns a
    // best-guess `reactionClass` string. Purely metadata: the underlying
    // mapping, score, candidates, decisiveness, and confidence fields are
    // unchanged. Reviewers get a one-word label they can scan instead of
    // squinting at bondChange counts to figure out what kind of reaction
    // was just mapped.
    //
    // Detection priority (first match wins; check more specific rules
    // before general ones):
    //
    //   identity           no bond changes of any kind
    //   stereo-only        no formed/cleaved/orderChange events, only
    //                      stereoChange or hydrogenChange events
    //   ring-opening       reactant has more rings than product
    //   ring-closing       product has more rings than reactant
    //   redox              only orderChange events on heavy atoms AND
    //                      hydrogenChange events present (covers
    //                      alcohol↔aldehyde, double-bond hydrogenation,
    //                      etc.; checked BEFORE isomerisation so an
    //                      oxidation isn't misread as a pure
    //                      bond-order change)
    //   addition           >=1 formed, 0 cleaved
    //   elimination        0 formed, >=1 cleaved
    //   substitution       >=1 formed AND >=1 cleaved (covers hydrolysis,
    //                      esterification, SN1, SN2, Friedel-Crafts,
    //                      named rearrangements that the corpus knows
    //                      are mechanistic rearrangements but whose
    //                      bond-change profile is indistinguishable
    //                      from substitution at this heuristic level —
    //                      the v2.0.12 classifier deliberately does not
    //                      try to discriminate Beckmann / Pinacol / etc.
    //                      from a simple substitution; that needs a
    //                      reaction-template library, queued for
    //                      future work)
    //   isomerisation      orderChange only, no formed/cleaved, no
    //                      hydrogenChange (covers cis/trans, tautomers)
    //   unclassified       fell through every rule (rare in practice)
    //
    // The result also carries a `classReason` field — a short human-
    // readable trace of WHICH counts triggered the classification, so
    // the heuristic is transparent rather than opaque. When the
    // `detectReactionClass` option is false (default true), both
    // `reactionClass` and `classReason` are omitted entirely.
    function detectReactionClass(result) {
        if (!result || !result.bondChanges || !result.sides) {
            return { reactionClass: 'unclassified', classReason: 'no result' };
        }
        var rSide = result.sides.reactants || [];
        var pSide = result.sides.products || [];
        var rAtoms = 0, pAtoms = 0, rRings = 0, pRings = 0, i;
        for (i = 0; i < rSide.length; i++) {
            rAtoms += (rSide[i].atoms || []).length;
            if (rSide[i].findRings) { rRings += rSide[i].findRings(8).length; }
        }
        for (i = 0; i < pSide.length; i++) {
            pAtoms += (pSide[i].atoms || []).length;
            if (pSide[i].findRings) { pRings += pSide[i].findRings(8).length; }
        }
        var formed = 0, cleaved = 0, orderChange = 0, stereoChange = 0, hChange = 0;
        var bc = result.bondChanges;
        for (i = 0; i < bc.length; i++) {
            switch (bc[i].type) {
                case 'formed':        formed++; break;
                case 'cleaved':       cleaved++; break;
                case 'orderChange':   orderChange++; break;
                case 'stereoChange':  stereoChange++; break;
                case 'hydrogenChange': hChange++; break;
            }
        }
        var heavyChanges = formed + cleaved + orderChange;
        var counts = 'f=' + formed + ' c=' + cleaved + ' o=' + orderChange +
                     ' s=' + stereoChange + ' h=' + hChange +
                     ' rings=' + rRings + '->' + pRings +
                     ' atoms=' + rAtoms + '->' + pAtoms;

        if (heavyChanges === 0 && stereoChange === 0 && hChange === 0) {
            return { reactionClass: 'identity', classReason: 'no bond changes (' + counts + ')' };
        }
        if (heavyChanges === 0 && (stereoChange > 0 || hChange > 0)) {
            return { reactionClass: 'stereo-only', classReason: counts };
        }
        if (rRings > pRings) {
            return { reactionClass: 'ring-opening', classReason: 'rings lost (' + counts + ')' };
        }
        if (pRings > rRings) {
            return { reactionClass: 'ring-closing', classReason: 'rings gained (' + counts + ')' };
        }
        // Redox MUST come before isomerisation / substitution: an
        // alcohol→aldehyde oxidation produces orderChange=1 + hChange=2,
        // which would otherwise hit the isomerisation rule.
        if (formed === 0 && cleaved === 0 && orderChange > 0 && hChange > 0) {
            return { reactionClass: 'redox', classReason: 'orderChange + hChange (' + counts + ')' };
        }
        if (formed > 0 && cleaved === 0) {
            return { reactionClass: 'addition', classReason: counts };
        }
        if (cleaved > 0 && formed === 0) {
            return { reactionClass: 'elimination', classReason: counts };
        }
        if (formed > 0 && cleaved > 0) {
            return { reactionClass: 'substitution', classReason: counts };
        }
        if (orderChange > 0 && formed === 0 && cleaved === 0 && hChange === 0) {
            return { reactionClass: 'isomerisation', classReason: 'orderChange only (' + counts + ')' };
        }
        return { reactionClass: 'unclassified', classReason: counts };
    }

    // -----------------------------------------------------------------------
    // Exports
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // v1.4.0: derive component pairs from a mapReaction result, for the
    // mol-mol pair-highlight overlay rendered by Renderer / MolEditor.
    //
    // Input  : the object returned by RDT.mapReaction (must include `pairs`,
    //          `sides`, and `mapping`).
    // Output : array of { reactantAtomIds, productAtomIds, color, pairIndex,
    //          reactantComponentIdx, productComponentIdx }.
    //
    // The pairing strategy is the simple bipartite/greedy one already used
    // by the strategy dispatcher: each entry in `result.pairs` carries
    // (rIdx, pIdx) — the component indices into result.sides.reactants /
    // result.sides.products. We translate those back to the atom-id arrays
    // and assign one of eight palette colours by pairIndex.
    // -----------------------------------------------------------------------
    function deriveComponentPairs(result, palette) {
        var defaultPalette = [
            '#0d9488', '#a855f7', '#f97316', '#06b6d4',
            '#84cc16', '#ec4899', '#eab308', '#3b82f6'
        ];
        var pal = (palette && palette.length) ? palette : defaultPalette;
        var out = [];
        if (!result || !result.pairs || !result.sides) { return out; }
        var rSide = result.sides.reactants || [];
        var pSide = result.sides.products || [];
        // Sort by (rIdx, pIdx) for fully deterministic colour assignment.
        var sorted = result.pairs.slice().sort(function(a, b) {
            if (a.rIdx !== b.rIdx) { return a.rIdx - b.rIdx; }
            return a.pIdx - b.pIdx;
        });
        for (var i = 0; i < sorted.length; i++) {
            var cell = sorted[i];
            if (!cell) { continue; }
            var rIdx = cell.rIdx, pIdx = cell.pIdx;
            if (rIdx < 0 || rIdx >= rSide.length) { continue; }
            if (pIdx < 0 || pIdx >= pSide.length) { continue; }
            var rAtoms = rSide[rIdx].atoms.map(function(a) { return a.id; });
            var pAtoms = pSide[pIdx].atoms.map(function(a) { return a.id; });
            out.push({
                reactantAtomIds: rAtoms,
                productAtomIds: pAtoms,
                color: pal[i % pal.length],
                pairIndex: i,
                reactantComponentIdx: rIdx,
                productComponentIdx: pIdx
            });
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // v1.4.0: derive a flat per-component-pair descriptor for the renderer
    // halo overlay. Each entry is:
    //   { reactantCompIdx, productCompIdx, mcsSize, paletteIndex }
    // where paletteIndex is 0..N-1 for paired components and -1 for any
    // component that did not participate in the pairing (e.g. spectator,
    // stoichiometric mismatch). Reactants and products are processed
    // independently so each side's unpaired components still get an entry
    // with paletteIndex = -1.
    // -----------------------------------------------------------------------
    function deriveComponentPairList(result) {
        var out = [];
        if (!result || !result.sides) { return out; }
        var rSide = result.sides.reactants || [];
        var pSide = result.sides.products || [];
        var pairs = result.pairs || [];
        // Sort pairs by (rIdx, pIdx) for deterministic palette assignment.
        var sorted = pairs.slice().sort(function(a, b) {
            if (a.rIdx !== b.rIdx) { return a.rIdx - b.rIdx; }
            return a.pIdx - b.pIdx;
        });
        var pairedR = {};
        var pairedP = {};
        var i;
        for (i = 0; i < sorted.length; i++) {
            var cell = sorted[i];
            if (!cell) { continue; }
            if (cell.rIdx < 0 || cell.rIdx >= rSide.length) { continue; }
            if (cell.pIdx < 0 || cell.pIdx >= pSide.length) { continue; }
            out.push({
                reactantCompIdx: cell.rIdx,
                productCompIdx: cell.pIdx,
                mcsSize: cell.mcsSize || 0,
                paletteIndex: i
            });
            pairedR[cell.rIdx] = true;
            pairedP[cell.pIdx] = true;
        }
        // Emit unpaired reactant components (paletteIndex = -1).
        for (i = 0; i < rSide.length; i++) {
            if (!pairedR[i]) {
                out.push({
                    reactantCompIdx: i,
                    productCompIdx: -1,
                    mcsSize: 0,
                    paletteIndex: -1
                });
            }
        }
        // Emit unpaired product components (paletteIndex = -1).
        for (i = 0; i < pSide.length; i++) {
            if (!pairedP[i]) {
                out.push({
                    reactantCompIdx: -1,
                    productCompIdx: i,
                    mcsSize: 0,
                    paletteIndex: -1
                });
            }
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // Mapping-confidence score in [0, 1] derived from a mapReaction result.
    // Higher = better. Defined as
    //
    //   confidence = coverage × parsimony
    //   coverage   = mappedAtoms / max(rTotal, pTotal)
    //   parsimony  = 1 / (1 + bondChanges / max(rTotal, pTotal))
    //
    // where mappedAtoms is the number of mapped reactant atoms and the
    // denominator is the larger of (reactant heavy-atom count, product
    // heavy-atom count). Used by MolEditor for the status-line
    // "Confidence = N.NNNNN" display and by gradeMapping() for the A/B/C/F
    // letter grade.
    //
    // v2.0.73 RECALIBRATION (closes the standing gap flagged in v2.0.64 /
    // v2.0.72). The pre-v2.0.73 parsimony term was 1/(1 + 0.1·bondChanges)
    // — an ABSOLUTE bond-change penalty. That conflated two different
    // things: "how uncertain is this mapping?" (a property of the mapping)
    // and "how many bonds does this reaction rearrange?" (a property of the
    // reaction). A complex but perfectly-mapped reaction (e.g. adenosine
    // deaminase: coverage 1.00, 6 bond changes) was driven to confidence
    // 0.625 → grade F, purely for being chemically elaborate. The bug was
    // that an absolute count is not scale-stable across reactions of
    // different size.
    //
    // The fix normalises bond changes by reaction size: parsimony now
    // penalises bond-change DENSITY (changes per heavy atom), not the raw
    // count. A few bond changes among many atoms → low density → high
    // confidence; near-total bond turnover (the signature of a scrambled /
    // wrong mapping, OR a genuinely tiny reaction) → high density → low
    // confidence. This is the Principle of Minimal Chemical Distance applied
    // at the right scale.
    //
    // The grade bands (C ≥ 0.50, B ≥ 0.75, A ≥ 0.95) are tuned against an
    // internal benchmark of known-correct mappings. The score is NOT a
    // calibrated probability — only a better-calibrated ordering signal.
    // -----------------------------------------------------------------------
    function deriveConfidence(result) {
        if (!result) { return 0; }
        var sides = result.sides;
        if (!sides || !sides.reactants || !sides.products) { return 0; }
        var i, j, atoms;
        var rTotal = 0, pTotal = 0;
        for (i = 0; i < sides.reactants.length; i++) {
            atoms = sides.reactants[i].atoms || [];
            rTotal += atoms.length;
        }
        for (i = 0; i < sides.products.length; i++) {
            atoms = sides.products[i].atoms || [];
            pTotal += atoms.length;
        }
        var heavy = Math.max(rTotal, pTotal, 1);
        var mappedAtoms = 0;
        if (result.mapping) {
            for (var k in result.mapping) {
                if (result.mapping.hasOwnProperty(k)) { mappedAtoms++; }
            }
        }
        var nChanges = 0;
        var changes = result.bondChanges || [];
        for (j = 0; j < changes.length; j++) {
            var t = changes[j].type;
            if (t === 'formed' || t === 'cleaved' || t === 'orderChange') { nChanges++; }
        }
        var coverage = mappedAtoms / heavy;
        // v2.0.73: size-normalised parsimony (bond-change density), not an
        // absolute count. See the calibration note above.
        var parsimony = 1 / (1 + nChanges / heavy);
        var c = coverage * parsimony;
        if (!isFinite(c) || c < 0) { c = 0; }
        if (c > 1) { c = 1; }
        return c;
    }

    // v2.0.41 — AAM quality diagnostics.
    //
    // These helpers compute READ-ONLY post-hoc metrics on a mapped result;
    // they do NOT participate in selection-time fitness, so adding them
    // cannot regress the golden corpus. They give callers a chemistry-
    // aware view of how good a mapping is beyond raw atom coverage.
    //
    // ARC = aromatic-ring conservation: for each aromatic ring on each
    // reactant fragment, count how many rings on the product side share
    // EXACTLY the same set of mapped atoms. A reaction that preserves all
    // aromatic systems gets ARC = 1.0; a reaction that breaks all of them
    // gets 0.
    function countAromaticRingsPreserved(result) {
        if (!result || !result.sides) { return { preserved: 0, total: 0 }; }
        var sides = result.sides;
        var mapping = result.mapping || {};
        var reactantsRings = collectAromaticRings(sides.reactants);
        var productsRings = collectAromaticRingSets(sides.products);
        var preserved = 0;
        for (var i = 0; i < reactantsRings.length; i++) {
            var ring = reactantsRings[i];
            // Project the reactant ring through the mapping; if any atom
            // is unmapped the ring is broken.
            var mappedSet = {};
            var fullyMapped = true;
            for (var j = 0; j < ring.length; j++) {
                var pid = mapping[ring[j]];
                if (typeof pid === 'undefined') { fullyMapped = false; break; }
                mappedSet[pid] = true;
            }
            if (!fullyMapped) { continue; }
            // Look for a product aromatic ring whose atom-set equals the
            // mappedSet exactly (order-independent).
            for (var k = 0; k < productsRings.length; k++) {
                var pRing = productsRings[k];
                if (pRing.length !== ring.length) { continue; }
                var match = true;
                for (var ki = 0; ki < pRing.length; ki++) {
                    if (!mappedSet[pRing[ki]]) { match = false; break; }
                }
                if (match) { preserved++; break; }
            }
        }
        return { preserved: preserved, total: reactantsRings.length };
    }

    function collectAromaticRings(mols) {
        // Returns an array of ring atom-id arrays for every aromatic ring
        // in every fragment. Molecule.findRings returns objects of shape
        // { atoms: [ids], key, size } — we lift the .atoms array out.
        var out = [];
        if (!mols) { return out; }
        for (var mi = 0; mi < mols.length; mi++) {
            var mol = mols[mi];
            if (!mol || typeof mol.findRings !== 'function') { continue; }
            var rings = mol.findRings();
            for (var ri = 0; ri < rings.length; ri++) {
                var ringObj = rings[ri];
                var ringAtoms = (ringObj && Array.isArray(ringObj.atoms)) ? ringObj.atoms : null;
                if (!ringAtoms || ringAtoms.length < 3) { continue; }
                var allArom = true;
                for (var ai = 0; ai < ringAtoms.length; ai++) {
                    var atom = mol.getAtom ? mol.getAtom(ringAtoms[ai]) : null;
                    if (!atom || !atom.aromatic) { allArom = false; break; }
                }
                if (allArom) { out.push(ringAtoms.slice()); }
            }
        }
        return out;
    }

    function collectAromaticRingSets(mols) {
        // Same as collectAromaticRings but returns rings as atom-id arrays
        // for the product side (so the caller can compare atom-sets).
        return collectAromaticRings(mols);
    }

    // STC = stereo-centre conservation: count atoms with a non-empty
    // CIP / wedge stereo descriptor that are mapped to a partner with
    // a non-empty descriptor. We don't require the descriptor to MATCH
    // (R↔S inversions are valid chemistry in some reactions); only that
    // the partner has SOME stereo descriptor too. This catches the
    // common bug where a mapping silently drops a chiral centre.
    function countStereoPreserved(result) {
        if (!result || !result.sides) { return { preserved: 0, total: 0 }; }
        var sides = result.sides;
        var mapping = result.mapping || {};
        var preserved = 0, total = 0;
        // Reactant stereo atoms.
        var pAtoms = indexAtoms(sides.products);
        for (var mi = 0; mi < sides.reactants.length; mi++) {
            var mol = sides.reactants[mi];
            if (!mol || !mol.atoms) { continue; }
            for (var ai = 0; ai < mol.atoms.length; ai++) {
                var a = mol.atoms[ai];
                if (!a || !hasStereoDescriptor(a)) { continue; }
                total++;
                var partnerId = mapping[a.id];
                if (typeof partnerId === 'undefined') { continue; }
                var partner = pAtoms[partnerId];
                if (partner && hasStereoDescriptor(partner)) { preserved++; }
            }
        }
        return { preserved: preserved, total: total };
    }

    function hasStereoDescriptor(atom) {
        if (!atom) { return false; }
        if (atom.cipLabel === 'R' || atom.cipLabel === 'S') { return true; }
        if (atom.chirality && atom.chirality !== '') { return true; }
        return false;
    }

    function indexAtoms(mols) {
        var byId = {};
        if (!mols) { return byId; }
        for (var mi = 0; mi < mols.length; mi++) {
            var mol = mols[mi];
            if (!mol || !mol.atoms) { continue; }
            for (var ai = 0; ai < mol.atoms.length; ai++) {
                var a = mol.atoms[ai];
                if (a) { byId[a.id] = a; }
            }
        }
        return byId;
    }

    // QUALITY GRADE: a categorical letter summarising a numeric mapping
    // result, layered on top of confidence + (optionally) decisiveness +
    // aromatic-ring + stereo.
    //
    //   A — confidence ≥ 0.95, decisiveness ≥ 0.3, ARC ≥ 0.95, stereo ≥ 0.9
    //   B — confidence ≥ 0.75, ARC ≥ 0.85, stereo ≥ 0.8
    //   C — confidence ≥ 0.5,  ARC ≥ 0.6
    //   F — anything mapped but failing C
    //   'unmapped' — status !== 'mapped'
    //
    // v2.0.72 recalibration. The pre-v2.0.72 thresholds also gated B and C
    // on `decisiveness ≥ 0.2` / `≥ 0.1`, but `deriveDecisiveness` returns
    // 0 whenever the top two strategies TIE at the optimum (a common
    // outcome for everyday reactions — symmetric substrates, simple
    // substitutions, etc.). That meant well-formed mappings with high
    // confidence and full aromatic preservation routinely dropped to F
    // for users running reactions outside the curated golden corpus.
    // Decisiveness is now retained only as the A-grade gate (where ties
    // genuinely should not earn the gold-standard label); B and C are
    // determined by confidence + aromaticity + stereo, which behave
    // smoothly on real-world inputs. The A threshold is unchanged, so
    // the curated golden corpus is unaffected; B/C are now strictly more
    // generous (any pre-v2.0.72 B/C still grades ≥ B/C) so no golden
    // reaction regresses.
    function gradeMapping(result) {
        if (!result || result.status !== 'mapped') { return 'unmapped'; }
        var confidence = (typeof result.confidence === 'number') ? result.confidence : 0;
        var decisiveness = (typeof result.decisiveness === 'number') ? result.decisiveness : 0;
        var arc = countAromaticRingsPreserved(result);
        var arcRate = arc.total > 0 ? arc.preserved / arc.total : 1;
        var stc = countStereoPreserved(result);
        var stcRate = stc.total > 0 ? stc.preserved / stc.total : 1;
        if (confidence >= 0.95 && decisiveness >= 0.3 && arcRate >= 0.95 && stcRate >= 0.9) {
            return 'A';
        }
        if (confidence >= 0.75 && arcRate >= 0.85 && stcRate >= 0.8) {
            return 'B';
        }
        if (confidence >= 0.5 && arcRate >= 0.6) {
            return 'C';
        }
        return 'F';
    }

    // Composite diagnostic returned alongside mapReaction results. Cheap
    // to compute (under 1 ms even for 50-atom reactions) so it ships in
    // every result by default. Callers that need the absolute minimal
    // result can pass `options.includeDiagnostic: false`.
    function deriveDiagnostic(result) {
        if (!result) { return null; }
        var sides = result.sides || {};
        var mapping = result.mapping || {};
        var rTotal = 0, pTotal = 0;
        var i, atoms;
        for (i = 0; i < (sides.reactants || []).length; i++) {
            atoms = sides.reactants[i].atoms || [];
            rTotal += atoms.length;
        }
        for (i = 0; i < (sides.products || []).length; i++) {
            atoms = sides.products[i].atoms || [];
            pTotal += atoms.length;
        }
        var mappedCount = 0;
        for (var k in mapping) {
            if (mapping.hasOwnProperty(k)) { mappedCount++; }
        }
        var arc = countAromaticRingsPreserved(result);
        var stc = countStereoPreserved(result);
        return {
            quality: gradeMapping(result),
            mappedCount: mappedCount,
            reactantAtomCount: rTotal,
            productAtomCount: pTotal,
            unmappedReactantCount: Math.max(0, rTotal - mappedCount),
            unmappedProductCount: Math.max(0, pTotal - mappedCount),
            // v2.0.49: per-atom unmapped-reason map.
            // { reactantAtomId | productAtomId: reasonString } — diagnostic
            // help for "why didn't atom X map?".
            unmappedReasons: deriveUnmappedReasons(result),
            aromaticRingsPreserved: arc.preserved,
            aromaticRingsTotal: arc.total,
            stereoCentresPreserved: stc.preserved,
            stereoCentresTotal: stc.total,
            confidence: (typeof result.confidence === 'number') ? result.confidence : 0,
            decisiveness: (typeof result.decisiveness === 'number') ? result.decisiveness : 0
        };
    }

    // v2.0.49: per-atom unmapped-reason map. For each atom that is NOT
    // present in the mapping (in either direction), classify the most
    // likely reason. The classification is heuristic but useful for
    // debugging AAM failures and powering richer atom-trace UI hints.
    //
    // Reason categories:
    //   'no-partner-element'   — no atom of the matching element exists on the other side
    //   'orphan-after-rescue'  — atom belongs to a component that wasn't paired
    //   'cofactor-side-only'   — atom is on the side whose stoichiometric counterpart is missing
    //   'no-isomorphic-partner' — every same-element candidate failed the structural extension
    //   'timeout'              — mapper hit options.timeoutMs before this atom
    //                            could be reached (v2.0.50)
    //
    // v2.0.50: when result.timedOut is true, atoms that would otherwise
    // have been classified 'no-isomorphic-partner' or 'orphan-after-rescue'
    // are promoted to 'timeout' — the structural-extension verdict isn't
    // trustworthy because the mapper aborted before exhausting candidates.
    // The two stoichiometric verdicts ('no-partner-element' and
    // 'cofactor-side-only') are NOT overridden: they're per-element counts
    // computed up front, independent of the strategy loop, so a timeout
    // doesn't invalidate them. This keeps the timeout category orthogonal
    // to "the reaction is unbalanced" — both can coexist in the same map.
    function deriveUnmappedReasons(result) {
        var reasons = {};
        if (!result || !result.sides) { return reasons; }
        var sides = result.sides;
        var mapping = result.mapping || {};
        // Build mapping reverse-index (productAtomId → reactantAtomId).
        var mappedReactants = {}, mappedProducts = {};
        for (var k in mapping) {
            if (!mapping.hasOwnProperty(k)) { continue; }
            mappedReactants[+k] = true;
            mappedProducts[mapping[k]] = true;
        }
        // Tally element multiplicities per side (for the no-partner-element check).
        var rByEl = elementMultiplicity(sides.reactants);
        var pByEl = elementMultiplicity(sides.products);
        // Scan reactants.
        scanSide(sides.reactants, mappedReactants, true);
        // Scan products.
        scanSide(sides.products, mappedProducts, false);

        function scanSide(mols, mappedSet, isReactant) {
            if (!mols) { return; }
            for (var mi = 0; mi < mols.length; mi++) {
                var mol = mols[mi];
                if (!mol || !mol.atoms) { continue; }
                for (var ai = 0; ai < mol.atoms.length; ai++) {
                    var atom = mol.atoms[ai];
                    if (!atom || mappedSet[atom.id]) { continue; }
                    var sym = atom.symbol || 'C';
                    var thisSide = isReactant ? rByEl : pByEl;
                    var otherSide = isReactant ? pByEl : rByEl;
                    var reason;
                    if ((otherSide[sym] || 0) === 0) {
                        reason = 'no-partner-element';
                    } else if ((otherSide[sym] || 0) < (thisSide[sym] || 0)) {
                        reason = 'cofactor-side-only';
                    } else {
                        // Same-element atoms exist on the other side, but the
                        // structural extension didn't include this atom.
                        // Could be MCS-coverage shortfall (no isomorphic
                        // partner) or a stray after rescue.
                        var hasResult = result.strategyResults &&
                                        result.strategyResults.length > 0;
                        reason = hasResult ? 'no-isomorphic-partner' : 'orphan-after-rescue';
                        // v2.0.50: if the mapper aborted on timeoutMs, the
                        // structural-extension verdict is unreliable (we don't
                        // know whether a longer run would have found a partner).
                        // Promote both no-isomorphic-partner and
                        // orphan-after-rescue to 'timeout'.
                        if (result.timedOut) {
                            reason = 'timeout';
                        }
                    }
                    reasons[atom.id] = reason;
                }
            }
        }
        return reasons;
    }

    function elementMultiplicity(mols) {
        var counts = {};
        if (!mols) { return counts; }
        for (var mi = 0; mi < mols.length; mi++) {
            var mol = mols[mi];
            if (!mol || !mol.atoms) { continue; }
            for (var ai = 0; ai < mol.atoms.length; ai++) {
                var sym = mol.atoms[ai].symbol || 'C';
                counts[sym] = (counts[sym] || 0) + 1;
            }
        }
        return counts;
    }

    // -----------------------------------------------------------------------
    // v1.4.1: derive per-MCS-sub-fragment halo groups for the renderer.
    //
    // A "sub-fragment" is a maximal connected subgraph of mapped reactant
    // atoms whose corresponding mapped product atoms are also connected via
    // the same (preserved) bond skeleton. Two mapped reactant atoms r1 and
    // r2 are placed in the same sub-fragment iff:
    //   1) they are bonded in the reactant molecule, and
    //   2) M[r1] and M[r2] are bonded in the product molecule.
    //
    // Bond-change events (cleavage, formation, order change) break the union,
    // so each sub-fragment corresponds to a "rigid scaffold piece" that
    // survived intact through the reaction. This matches the RDT-style
    // three-colour visualisation in classic enzyme mapping diagrams: e.g.
    // dihydrostreptomycin + ATP -> bisphosphate + ADP yields three groups
    // (streptomycin scaffold preserved, adenosine scaffold preserved, gamma
    // phosphate migrating from ATP to streptomycin).
    //
    // Output: array of
    //   { reactantAtomIds: [..ids], productAtomIds: [..ids],
    //     paletteIndex: N, size: N }
    // ordered by descending size (ties by smallest atom-id) for deterministic
    // colour assignment. Unmapped atoms do not appear in the output (they get
    // no halo). Sub-fragments of size < options.minSize (default 1) are
    // filtered out — useful for less cluttered display when minSize=2.
    //
    // Bond-order changes (e.g. C=C -> C-C reduction): we treat both atom
    // pairs as STILL bonded — the union-find adjacency check uses presence
    // only, not bond order. This is the deliberate "rigid scaffold" reading
    // (atom connectivity preserved = same skeleton, even if pi-system
    // changed). Bond-order deltas are still surfaced separately as
    // `orderChange` events in result.bondChanges, so the user sees both
    // the conserved-skeleton colour AND the bond-order annotation.
    //
    // Complexity: O((R + B_r) log d_max + B_p) where R = mapped atom count,
    //   B_r = reactant bond count, B_p = product bond count, d_max =
    //   maximum atom degree. The log d_max factor comes from sorting each
    //   atom's neighbour list for deterministic union order; for chemistry
    //   (d_max <= 4 typically) it is effectively linear.
    //
    // Determinism: union-find walks atom-id-sorted neighbour lists; bucket
    // ordering uses (size desc, smallest-atom-id asc) as a total order, so
    // re-running with the same input yields byte-identical output.
    // -----------------------------------------------------------------------
    function deriveSubFragments(result, options) {
        var opts = options || {};
        // v1.4.1 (bug-fix #19): accept minSize === 0 ("keep all") and only
        // fall back to the default when undefined / non-numeric. The previous
        // `> 0` check silently coerced 0 to 1.
        var minSize = (typeof opts.minSize === 'number' && opts.minSize >= 0) ? opts.minSize : 1;
        var out = [];
        if (!result || !result.mapping || !result.sides) { return out; }
        var rSide = result.sides.reactants || [];
        var pSide = result.sides.products || [];
        var mapping = result.mapping;

        var rAtomSet = {};
        var pAtomSet = {};
        var rAdj = {};
        var pAdj = {};
        var i, j, mol, a, b, key;

        for (i = 0; i < rSide.length; i++) {
            mol = rSide[i];
            if (!mol || !mol.atoms) { continue; }
            for (j = 0; j < mol.atoms.length; j++) {
                a = mol.atoms[j];
                rAtomSet[a.id] = true;
                if (!rAdj[a.id]) { rAdj[a.id] = []; }
            }
            if (mol.bonds) {
                for (j = 0; j < mol.bonds.length; j++) {
                    b = mol.bonds[j];
                    if (!rAdj[b.atom1]) { rAdj[b.atom1] = []; }
                    if (!rAdj[b.atom2]) { rAdj[b.atom2] = []; }
                    rAdj[b.atom1].push(b.atom2);
                    rAdj[b.atom2].push(b.atom1);
                }
            }
        }
        for (i = 0; i < pSide.length; i++) {
            mol = pSide[i];
            if (!mol || !mol.atoms) { continue; }
            for (j = 0; j < mol.atoms.length; j++) {
                a = mol.atoms[j];
                pAtomSet[a.id] = true;
                if (!pAdj[a.id]) { pAdj[a.id] = []; }
            }
            if (mol.bonds) {
                for (j = 0; j < mol.bonds.length; j++) {
                    b = mol.bonds[j];
                    if (!pAdj[b.atom1]) { pAdj[b.atom1] = []; }
                    if (!pAdj[b.atom2]) { pAdj[b.atom2] = []; }
                    pAdj[b.atom1].push(b.atom2);
                    pAdj[b.atom2].push(b.atom1);
                }
            }
        }

        // Mapped reactant atoms only. rPair: rId -> pId. pAdjSet: O(1) lookup.
        var rPair = {};
        var rIds = [];
        var pAdjSet = {};
        for (key in mapping) {
            if (!mapping.hasOwnProperty(key)) { continue; }
            var rId = +key;
            var pId = +mapping[key];
            if (!rAtomSet[rId] || !pAtomSet[pId]) { continue; }
            rPair[rId] = pId;
            rIds.push(rId);
        }
        // Build product adjacency lookup keyed by 'p1|p2' (smaller id first).
        for (key in pAdj) {
            if (!pAdj.hasOwnProperty(key)) { continue; }
            var src = +key;
            var nbrs = pAdj[key];
            for (j = 0; j < nbrs.length; j++) {
                var dst = nbrs[j];
                var lo = src < dst ? src : dst;
                var hi = src < dst ? dst : src;
                pAdjSet[lo + '|' + hi] = true;
            }
        }
        rIds.sort(function(a, b) { return a - b; });

        // Union-find over mapped reactant atom ids.
        var parent = {};
        var rank = {};
        function find(x) {
            var p = parent[x];
            if (p === x) { return x; }
            // Iterative path-halving (O(alpha) and avoids stack blow-up).
            while (parent[p] !== p) { parent[x] = parent[p]; x = p; p = parent[p]; }
            return p;
        }
        function union(a, b) {
            var ra = find(a), rb = find(b);
            if (ra === rb) { return; }
            if (rank[ra] < rank[rb]) { parent[ra] = rb; }
            else if (rank[ra] > rank[rb]) { parent[rb] = ra; }
            else { parent[rb] = ra; rank[ra]++; }
        }
        for (i = 0; i < rIds.length; i++) {
            parent[rIds[i]] = rIds[i];
            rank[rIds[i]] = 0;
        }

        // Walk every reactant bond r1-r2 where both atoms are mapped. If the
        // corresponding atoms in the product are also bonded, union r1 r2.
        for (i = 0; i < rIds.length; i++) {
            var r1 = rIds[i];
            var p1 = rPair[r1];
            var nbrList = rAdj[r1] || [];
            // Sort neighbour ids ascending for deterministic union order.
            var sortedNbrs = nbrList.slice().sort(function(a, b) { return a - b; });
            for (j = 0; j < sortedNbrs.length; j++) {
                var r2 = sortedNbrs[j];
                if (r2 <= r1) { continue; }
                if (rPair[r2] === undefined) { continue; }
                var p2 = rPair[r2];
                var lo2 = p1 < p2 ? p1 : p2;
                var hi2 = p1 < p2 ? p2 : p1;
                if (pAdjSet[lo2 + '|' + hi2]) { union(r1, r2); }
            }
        }

        // Bucket reactant atoms by union-find representative.
        var buckets = {};
        for (i = 0; i < rIds.length; i++) {
            var rep = find(rIds[i]);
            if (!buckets[rep]) { buckets[rep] = []; }
            buckets[rep].push(rIds[i]);
        }
        var groups = [];
        for (key in buckets) {
            if (buckets.hasOwnProperty(key)) { groups.push(buckets[key]); }
        }
        // Deterministic order: largest first, ties by smallest atom-id.
        // v1.4.1 (panel-review minor): use min over the bucket explicitly so
        // a future refactor that breaks "rIds is ascending so a[0] === min"
        // does not silently break determinism. Each atom-id appears in
        // exactly one bucket (partition property), so min(bucket) is unique
        // across buckets and the sort is total.
        function minOf(arr) {
            var m = arr[0];
            for (var k = 1; k < arr.length; k++) { if (arr[k] < m) { m = arr[k]; } }
            return m;
        }
        groups.sort(function(a, b) {
            if (a.length !== b.length) { return b.length - a.length; }
            return minOf(a) - minOf(b);
        });

        var idx = 0;
        for (i = 0; i < groups.length; i++) {
            var rGroup = groups[i];
            if (rGroup.length < minSize) { continue; }
            // v1.4.1 (bug-fix #8): clone the bucket arrays before sorting so
            // a caller who holds a reference from a previous deriveSubFragments
            // call doesn't see their array re-sorted by a later call. The
            // bucket arrays themselves are owned by the union-find and must
            // not leak into the public output.
            var rGroupOut = rGroup.slice().sort(function(a, b) { return a - b; });
            var pGroupOut = [];
            for (j = 0; j < rGroupOut.length; j++) { pGroupOut.push(rPair[rGroupOut[j]]); }
            pGroupOut.sort(function(a, b) { return a - b; });
            out.push({
                reactantAtomIds: rGroupOut,
                productAtomIds: pGroupOut,
                paletteIndex: idx,
                size: rGroupOut.length
            });
            idx++;
        }
        return out;
    }

    var RDT = {
        mapReaction: mapReaction,
        runMinPairwise: runMinPairwise,
        runMaxPairwise: runMaxPairwise,
        runMixturePairwise: runMixturePairwise,
        runRingPairwise: runRingPairwise,
        annotateBondChanges: annotateBondChanges,
        fitness: fitness,
        deriveComponentPairs: deriveComponentPairs,
        deriveComponentPairList: deriveComponentPairList,
        deriveConfidence: deriveConfidence,
        // v2.0.41: AAM quality diagnostic surface.
        gradeMapping: gradeMapping,
        deriveDiagnostic: deriveDiagnostic,
        _countAromaticRingsPreserved: countAromaticRingsPreserved,
        _countStereoPreserved: countStereoPreserved,
        // v2.0.49: per-atom unmapped-reason classifier.
        deriveUnmappedReasons: deriveUnmappedReasons,
        deriveSubFragments: deriveSubFragments,
        STRATEGIES: STRATEGIES.slice(),
        // Internals re-exported for tests
        _splitReactionSides: splitReactionSides,
        _checkBalance: checkBalance,
        _buildSimilarityCache: buildSimilarityCache,
        _pickWinner: pickWinner,
        _applyMappingToReaction: applyMappingToReaction,
        _bipartiteComponentPairing: _bipartiteComponentPairing,
        _runCipOnSides: runCipOnSides,
        // v1.4.2: leftover-atom rescue helpers (exported for tests).
        _rescueLeftoverAtoms: rescueLeftoverAtoms,
        _enrichComponentPairsWithRescue: enrichComponentPairsWithRescue,
        // v2.0.14: reaction-class detector (exported for tests).
        _detectReactionClass: detectReactionClass,
        version: '3.0.3'
    };

    global.RDT = RDT;

})(typeof window !== 'undefined' ? window : globalThis);
