/**
 * RDT.js - Reaction Decoder Tool: pure-JavaScript port for BIME v1.1.0.
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
 *   4. Annotate bond-change events (formed / cleaved / orderChange / stereoChange).
 *   5. Score each candidate via the fitness function (lower = better).
 *   6. Select winner by fitness, with deterministic tie-break.
 *
 * Public API
 * ----------
 *   RDT.mapReaction(reaction, options) -> { mapping, bondChanges, score, strategy, ... }
 *   RDT.runMinPairwise(reaction, options)
 *   RDT.runMaxPairwise(reaction, options)
 *   RDT.runMixturePairwise(reaction, options)
 *   RDT.runRingPairwise(reaction, options)
 *   RDT.annotateBondChanges(reactantSide, productSide, atomMap) -> bondChanges[]
 *   RDT.fitness(reaction, atomMap, bondChanges) -> number
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
    function annotateBondChanges(reactantSide, productSide, atomMap) {
        var rBonds = buildBondLookupBySide(reactantSide);
        var pBonds = buildBondLookupBySide(productSide);
        var events = [];
        var seenProductBonds = {};

        // Atom lookup for stereo info.
        var rAtoms = {};
        var i;
        for (i = 0; i < reactantSide.length; i++) {
            for (var ji = 0; ji < reactantSide[i].atoms.length; ji++) {
                rAtoms[reactantSide[i].atoms[ji].id] = reactantSide[i].atoms[ji];
            }
        }
        var pAtoms = {};
        for (i = 0; i < productSide.length; i++) {
            for (var ji2 = 0; ji2 < productSide[i].atoms.length; ji2++) {
                pAtoms[productSide[i].atoms[ji2].id] = productSide[i].atoms[ji2];
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

        // Pass 3: stereo changes at mapped atoms
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
                    afterStereo: afterRS
                });
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

        var bondChanges = annotateBondChanges(sides.reactants, sides.products, atomMap);
        var score = fitness({ reactants: sides.reactants, products: sides.products }, atomMap, bondChanges);

        return {
            strategy: strategyName,
            mapping: atomMap,
            bondChanges: bondChanges,
            score: score,
            pairs: pairs,
            sides: sides
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

    function mapReaction(reaction, options) {
        var opts = shallowCopyOpts(options);
        var t0 = nowMs();
        var timeoutMs = opts.timeoutMs || DEFAULT_TIMEOUT_MS;
        var debug = !!opts.debug;
        var enabledStrategies = (opts.strategies && opts.strategies.length > 0) ? opts.strategies : STRATEGIES;

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
                strategyResults: results
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
        version: '1.1.0'
    };

    global.RDT = RDT;

})(typeof window !== 'undefined' ? window : globalThis);
