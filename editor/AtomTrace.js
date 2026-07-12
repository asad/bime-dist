(function(global) {
    'use strict';
    /**
     * AtomTrace.js — chain per-reaction atom–atom mappings into cross-pathway
     * atom traces, so a single atom can be followed from a start metabolite
     * through every downstream metabolite (the "fate of atoms" view used in
     * 13C metabolic-flux tracing).
     *
     * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
     * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
     *
     * How it works
     * ------------
     * Each pathway step is one reaction `from >> to`. We run BIME's RDT atom–
     * atom mapping (RDT.js) on it in balance-tolerant mode (cofactor atoms such
     * as the phosphate from ATP are not in the half-reaction, so the heavy-atom
     * counts need not balance). RDT returns `mapping` = { reactantAtomId ->
     * productAtomId } and `sides` = the reactant/product atom lists. We convert
     * atom ids to *indices within their molecule* (SMILES parse order is
     * deterministic), giving an index->index map per step. Because the shared
     * metabolite is the same SMILES string on both sides of the join (product
     * of step i, reactant of step i+1), its atom indices line up across steps —
     * so a start atom's index can be propagated forward step by step. When a
     * start atom no longer has an image (it forked to a co-product, e.g. the
     * aldolase split, or left as CO2/H2O), its trace ends and we record where.
     *
     * Pure (no DOM); depends on the global Molecule, SmilesParser, and RDT.
     *
     * Caller contract: pass ONE canonical SMILES per metabolite. The chaining
     * relies on a shared metabolite being the identical SMILES string on both
     * sides of a join (so its deterministic parse order — hence atom indices —
     * matches). Supplying the same molecule via two different SMILES (e.g.
     * "CC=O" vs "O=CC") would misalign indices. Within a single tracePathway
     * call this is guaranteed (one list entry feeds both join sides).
     */

    function flattenSideAtoms(side) {
        var out = [];
        if (!side) { return out; }
        for (var i = 0; i < side.length; i++) {
            var atoms = side[i] && side[i].atoms ? side[i].atoms : [];
            for (var j = 0; j < atoms.length; j++) { out.push(atoms[j]); }
        }
        return out;
    }

    // v2.0.56: cofactor SMILES table. Small, hand-curated subset of the
    // metabolic-cofactor universe — enough to wire phosphate / water / proton
    // lineage through glycolysis, TCA, urea, PPP, and β-oxidation. The
    // monster cofactors (ATP/ADP/AMP, NAD+/NADH, NADP+/NADPH, CoA / Acetyl-
    // CoA, FAD / FADH2, GTP / GDP) carry their full canonical SMILES so the
    // adenosine / nicotinamide scaffolds can ride along too. Keys are matched
    // case-insensitively against pathway-step cofactor strings (parsed via
    // parseCofactorString below).
    var COFACTOR_SMILES = {
        // Simple inorganics — most commonly transferred fragments in
        // central metabolism.
        'pi':     'OP(O)(O)=O',
        'ppi':    'OP(=O)(O)OP(O)(O)=O',
        'h2o':    'O',
        'water':  'O',
        'co2':    'O=C=O',
        'h+':     '[H+]',
        'h':      '[H+]',
        'nh3':    'N',
        'nh4+':   '[NH4+]',
        'o2':     'O=O',
        // Adenine nucleotides — minus-water-of-condensation SMILES so the
        // γ-phosphate (ATP) / β-phosphate (ADP) is the rightmost P.
        'atp':    'Nc1ncnc2c1ncn2[C@@H]3O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]3O',
        'adp':    'Nc1ncnc2c1ncn2[C@@H]3O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]3O',
        'amp':    'Nc1ncnc2c1ncn2[C@@H]3O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]3O',
        // Guanine nucleotides.
        'gtp':    'Nc1nc2c(ncn2[C@@H]3O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]3O)c(=O)[nH]1',
        'gdp':    'Nc1nc2c(ncn2[C@@H]3O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]3O)c(=O)[nH]1',
        // NAD+ / NADH redox pair. The C4 of nicotinamide is the formal
        // hydride acceptor.
        'nad+':   'NC(=O)c1ccc[n+](C2O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@@H]3O[C@H](n4cnc5c(N)ncnc54)[C@@H](O)[C@H]3O)[C@@H](O)[C@H]2O)c1',
        'nadh':   'NC(=O)C1=CN(C2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@@H]3O[C@H](n4cnc5c(N)ncnc54)[C@@H](O)[C@H]3O)[C@@H](O)[C@H]2O)CC=C1',
        'nadp+':  'NC(=O)c1ccc[n+](C2O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@@H]3O[C@H](n4cnc5c(N)ncnc54)[C@@H](OP(O)(O)=O)[C@H]3O)[C@@H](O)[C@H]2O)c1',
        'nadph':  'NC(=O)C1=CN(C2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@@H]3O[C@H](n4cnc5c(N)ncnc54)[C@@H](OP(O)(O)=O)[C@H]3O)[C@@H](O)[C@H]2O)CC=C1'
    };

    // Aliases that point at canonical COFACTOR_SMILES keys.
    var COFACTOR_ALIASES = {
        'phosphate':         'pi',
        'orthophosphate':    'pi',
        'pyrophosphate':     'ppi',
        'inorganic phosphate': 'pi',
        'carbon dioxide':    'co2',
        'proton':            'h+',
        'h2 o':              'h2o',
        'ammonia':           'nh3',
        'dioxygen':          'o2',
        'molecular oxygen':  'o2',
        'nicotinamide adenine dinucleotide': 'nad+',
        'nadh + h+':         'nadh'
    };

    function lookupCofactorSmiles(name) {
        if (!name) { return null; }
        var key = String(name).toLowerCase().trim();
        if (COFACTOR_SMILES[key]) { return COFACTOR_SMILES[key]; }
        if (COFACTOR_ALIASES[key]) {
            var canon = COFACTOR_ALIASES[key];
            return COFACTOR_SMILES[canon] || null;
        }
        return null;
    }

    // v2.0.56: parse the MetaboliteLibrary `cofactors` string into a structured
    // {in:[name…], out:[name…]} record. Existing pathway templates carry
    // strings like 'ATP → ADP', '+ Acetyl-CoA + H2O → CoA-SH', 'NAD+ + Pi →
    // NADH', '− H2O'. The grammar is informal but consistent:
    //   - Tokens are separated by ' + ' (in/out same side) and ' → ' (in vs out).
    //   - A leading '+ X' on the IN side means "add X"; e.g. step 1 of TCA
    //     ('+ Acetyl-CoA + H2O → CoA-SH') reads as in={Acetyl-CoA, H2O},
    //     out={CoA-SH}.
    //   - A leading '− X' on a side means "subtract" — when there is no
    //     '→' separator, '− X' is interpreted as out-only (the metabolite
    //     LOSES X). e.g. enolase carries '− H2O' meaning out={H2O}.
    // Unknown cofactor tokens are kept in the structured record but won't
    // contribute SMILES (their atoms simply don't appear in the balanced rxn).
    function parseCofactorString(s) {
        if (!s || typeof s !== 'string') { return { in: [], out: [] }; }
        var raw = s.replace(/→/g, '->').replace(/←/g, '<-').replace(/−/g, '-');
        // Detect a one-side "loses X" pattern: leading '-' with no '->'.
        if (raw.indexOf('->') === -1 && /^\s*-\s*/.test(raw)) {
            return { in: [], out: splitPlusList(raw.replace(/^\s*-\s*/, '')) };
        }
        // Detect a one-side "gains X" pattern: leading '+' with no '->'.
        if (raw.indexOf('->') === -1 && /^\s*\+\s*/.test(raw)) {
            return { in: splitPlusList(raw.replace(/^\s*\+\s*/, '')), out: [] };
        }
        var parts = raw.split('->');
        if (parts.length === 2) {
            return {
                in:  splitPlusList(parts[0].replace(/^\s*\+\s*/, '')),
                out: splitPlusList(parts[1])
            };
        }
        return { in: [], out: [] };
    }
    function splitPlusList(s) {
        if (!s) { return []; }
        // v2.0.56: split STRICTLY on whitespace-plus-whitespace (' + ') so
        // that compound cofactor names ending in '+' or '-' (NAD+, NADP+,
        // H+, etc.) survive the split intact. Previous attempts using
        // `\s*\+\s*` or `\s+\+` ate the '+' in 'NAD+'.
        var tokens = s.split(/\s+\+\s+/);
        var out = [];
        for (var i = 0; i < tokens.length; i++) {
            var t = tokens[i].trim();
            if (t) { out.push(t); }
        }
        return out;
    }

    // v2.0.56: extended mapStep that ALSO parses cofactors into the balanced
    // reaction so AAM can see them. Returns the existing single-side
    // mapping (metabolite→metabolite) PLUS three cross-side maps that
    // surface cofactor lineage:
    //   - cofactorInToProductMetabolite — atoms that entered via a
    //     cofactor IN (e.g. ATP's γ-phosphate) and ended up in the
    //     product metabolite. Keyed by 'cofName:atomIdxInCofactor' → idxInProductMetabolite.
    //   - reactantMetaboliteToCofactorOut — atoms that left the reactant
    //     metabolite into a cofactor OUT (e.g. ADP's β-phosphate carrying
    //     the metabolite-originated atoms). Keyed by reactantMetaboliteIdx →
    //     'cofName:atomIdxInCofactor'.
    //   - cofactorInToCofactorOut — cofactor→cofactor (e.g. ATP→ADP carry-over
    //     atoms). Mostly informational.
    // Backward compat: callers that pass no cofactors get back the same
    // shape mapStep returned in v2.0.20-v2.0.55.
    function mapStepWithCofactors(fromSmiles, toSmiles, cofactorsIn, cofactorsOut, options) {
        options = options || {};
        cofactorsIn = cofactorsIn || [];
        cofactorsOut = cofactorsOut || [];

        // Resolve cofactor SMILES.
        var inSmilesList = [];
        var inResolved = [];
        for (var i = 0; i < cofactorsIn.length; i++) {
            var smi = lookupCofactorSmiles(cofactorsIn[i]);
            if (smi) { inSmilesList.push(smi); inResolved.push(cofactorsIn[i]); }
        }
        var outSmilesList = [];
        var outResolved = [];
        for (var j = 0; j < cofactorsOut.length; j++) {
            var smi2 = lookupCofactorSmiles(cofactorsOut[j]);
            if (smi2) { outSmilesList.push(smi2); outResolved.push(cofactorsOut[j]); }
        }

        var combinedFrom = fromSmiles;
        var combinedTo = toSmiles;
        if (inSmilesList.length) { combinedFrom = fromSmiles + '.' + inSmilesList.join('.'); }
        if (outSmilesList.length) { combinedTo = toSmiles + '.' + outSmilesList.join('.'); }

        var rxn = SmilesParser.parse(combinedFrom + '>>' + combinedTo);
        var res = RDT.mapReaction(rxn, {
            requireStoichiometricBalance: false,
            shareMcsCache: options.shareMcsCache
        });
        var rAtoms = flattenSideAtoms(res.sides && res.sides.reactants);
        var pAtoms = flattenSideAtoms(res.sides && res.sides.products);
        var idToR = {}, idToP = {};
        for (var k = 0; k < rAtoms.length; k++) { idToR[rAtoms[k].id] = k; }
        for (var p = 0; p < pAtoms.length; p++) { idToP[pAtoms[p].id] = p; }

        // Compute partition boundaries: atom index < metaboliteFromCount
        // belongs to the from-metabolite; the rest belong to cofactor-in[i]
        // in their declared order. Same for the product side.
        var metaboliteFromCount = countAtomsInSmiles(fromSmiles);
        var metaboliteToCount = countAtomsInSmiles(toSmiles);
        var inBoundaries = [metaboliteFromCount];
        var inOffset = metaboliteFromCount;
        for (var ic = 0; ic < inSmilesList.length; ic++) {
            inOffset += countAtomsInSmiles(inSmilesList[ic]);
            inBoundaries.push(inOffset);
        }
        var outBoundaries = [metaboliteToCount];
        var outOffset = metaboliteToCount;
        for (var oc = 0; oc < outSmilesList.length; oc++) {
            outOffset += countAtomsInSmiles(outSmilesList[oc]);
            outBoundaries.push(outOffset);
        }

        // Classify each atom-index by side membership.
        function classifyReactant(idx) {
            if (idx < inBoundaries[0]) { return { side: 'metabolite', localIdx: idx }; }
            for (var b = 0; b < inResolved.length; b++) {
                if (idx < inBoundaries[b + 1]) {
                    return { side: 'cofactor', name: inResolved[b], localIdx: idx - inBoundaries[b] };
                }
            }
            return { side: 'unknown', localIdx: idx };
        }
        function classifyProduct(idx) {
            if (idx < outBoundaries[0]) { return { side: 'metabolite', localIdx: idx }; }
            for (var b2 = 0; b2 < outResolved.length; b2++) {
                if (idx < outBoundaries[b2 + 1]) {
                    return { side: 'cofactor', name: outResolved[b2], localIdx: idx - outBoundaries[b2] };
                }
            }
            return { side: 'unknown', localIdx: idx };
        }

        // Walk the AAM mapping and partition into four buckets.
        var f2t = {};   // metabolite -> metabolite (the existing mapStep return)
        var cofactorInToProductMetabolite = {};
        var reactantMetaboliteToCofactorOut = {};
        var cofactorInToCofactorOut = {};
        var m = res.mapping || {};
        for (var rid in m) {
            if (!m.hasOwnProperty(rid)) { continue; }
            var ri = idToR[rid], pi = idToP[m[rid]];
            if (ri === undefined || pi === undefined) { continue; }
            var rClass = classifyReactant(ri);
            var pClass = classifyProduct(pi);
            if (rClass.side === 'metabolite' && pClass.side === 'metabolite') {
                f2t[rClass.localIdx] = pClass.localIdx;
            } else if (rClass.side === 'cofactor' && pClass.side === 'metabolite') {
                cofactorInToProductMetabolite[rClass.name + ':' + rClass.localIdx] = pClass.localIdx;
            } else if (rClass.side === 'metabolite' && pClass.side === 'cofactor') {
                reactantMetaboliteToCofactorOut[rClass.localIdx] = pClass.name + ':' + pClass.localIdx;
            } else if (rClass.side === 'cofactor' && pClass.side === 'cofactor') {
                cofactorInToCofactorOut[rClass.name + ':' + rClass.localIdx] =
                    pClass.name + ':' + pClass.localIdx;
            }
        }

        return {
            fromToIndex: f2t,
            cofactorInToProductMetabolite: cofactorInToProductMetabolite,
            reactantMetaboliteToCofactorOut: reactantMetaboliteToCofactorOut,
            cofactorInToCofactorOut: cofactorInToCofactorOut,
            cofactorsInResolved: inResolved,
            cofactorsOutResolved: outResolved,
            metaboliteFromCount: metaboliteFromCount,
            metaboliteToCount: metaboliteToCount,
            mappedCount: Object.keys(f2t).length,
            cofactorXferIn: Object.keys(cofactorInToProductMetabolite).length,
            cofactorXferOut: Object.keys(reactantMetaboliteToCofactorOut).length,
            status: res.status
        };
    }

    // Helper: parse a single SMILES + return its atom count (for the boundary
    // bookkeeping above). Falls back to 0 on parse failure.
    function countAtomsInSmiles(smi) {
        if (!smi) { return 0; }
        try {
            var mol = new Molecule();
            SmilesParser.parse(smi, mol);
            return mol.atoms.length;
        } catch (e) { return 0; }
    }

    // v2.4.2: SHARED RDT-id -> per-SMILES POSITIONAL-index translation. This is
    // the SINGLE place that runs `from >> to` through RDT.mapReaction and then
    // rewrites RDT's internal atom-ids into the deterministic SMILES-parse-order
    // indices that the trace chainer (tracePathway) consumes — i.e. the
    // `fromToIndex` convention. Both mapStep (the trace's own per-edge step) and
    // workbench's linkImportedStepAAM call THIS function, so a stored
    // edge.aam.mapping is guaranteed to live in EXACTLY the same positional
    // index space tracePathway later reads. NEVER store RDT internal atom-ids:
    // route through here instead. Returns the same shape mapStep's non-cofactor
    // branch always returned, plus the raw `score` (lower-is-better fitness from
    // RDT) so callers that want a confidence proxy don't re-derive it.
    function mapStepIndexMap(fromSmiles, toSmiles, options) {
        options = options || {};
        var rxn = SmilesParser.parse(fromSmiles + '>>' + toSmiles);
        var res = RDT.mapReaction(rxn, {
            requireStoichiometricBalance: false,
            shareMcsCache: options.shareMcsCache
        });
        var rAtoms = flattenSideAtoms(res.sides && res.sides.reactants);
        var pAtoms = flattenSideAtoms(res.sides && res.sides.products);
        var idToR = {}, idToP = {};
        for (var i = 0; i < rAtoms.length; i++) { idToR[rAtoms[i].id] = i; }
        for (var j = 0; j < pAtoms.length; j++) { idToP[pAtoms[j].id] = j; }
        var f2t = {};
        var m = res.mapping || {};
        for (var rid in m) {
            if (m.hasOwnProperty(rid)) {
                var ri = idToR[rid], pi = idToP[m[rid]];
                // Single-image assumption: AAM is 1->1, so each reactant index
                // has one product index. If RDT ever emitted 1->many for an id,
                // this is deterministic last-write-wins (crash-free).
                if (ri !== undefined && pi !== undefined) { f2t[ri] = pi; }
            }
        }
        return {
            fromToIndex: f2t,
            reactantElements: rAtoms.map(function(a) { return a.symbol; }),
            productElements: pAtoms.map(function(a) { return a.symbol; }),
            mappedCount: Object.keys(f2t).length,
            score: (typeof res.score === 'number') ? res.score : null,
            status: res.status
        };
    }

    // Atom-map one reaction step (from >> to) and return the index->index
    // correspondence (reactant atom index -> product atom index), plus per-side
    // element lists and the AAM status.
    //
    // v2.0.56: if options.cofactorsIn / cofactorsOut are present, delegate
    // to mapStepWithCofactors and merge its single-side mapping back into
    // the legacy return shape so existing callers see the same fields.
    // v2.4.2: the bare (no-cofactor) path now delegates to the shared
    // mapStepIndexMap so the id->index translation is defined ONCE and reused
    // by workbench's linkImportedStepAAM (same positional space).
    function mapStep(fromSmiles, toSmiles, options) {
        options = options || {};
        if ((options.cofactorsIn && options.cofactorsIn.length) ||
                (options.cofactorsOut && options.cofactorsOut.length)) {
            var ext = mapStepWithCofactors(fromSmiles, toSmiles,
                options.cofactorsIn || [], options.cofactorsOut || [], options);
            // Element lists are taken from the bare metabolite-side parse so
            // that callers indexing by "metabolite atom index" don't have to
            // think about cofactor padding.
            var rxnPlain = SmilesParser.parse(fromSmiles + '>>' + toSmiles);
            var rPlain = flattenSideAtoms(rxnPlain.sides && rxnPlain.sides.reactants);
            var pPlain = flattenSideAtoms(rxnPlain.sides && rxnPlain.sides.products);
            return {
                fromToIndex: ext.fromToIndex,
                reactantElements: rPlain.map(function(a) { return a.symbol; }),
                productElements: pPlain.map(function(a) { return a.symbol; }),
                mappedCount: ext.mappedCount,
                cofactorXferIn: ext.cofactorXferIn,
                cofactorXferOut: ext.cofactorXferOut,
                cofactorInToProductMetabolite: ext.cofactorInToProductMetabolite,
                reactantMetaboliteToCofactorOut: ext.reactantMetaboliteToCofactorOut,
                status: ext.status
            };
        }
        var indexMap = mapStepIndexMap(fromSmiles, toSmiles, options);
        return {
            fromToIndex: indexMap.fromToIndex,
            reactantElements: indexMap.reactantElements,
            productElements: indexMap.productElements,
            mappedCount: indexMap.mappedCount,
            status: indexMap.status
        };
    }

    // v2.0.53: Kahn's-algorithm topological sort over (n vertices, edges).
    // Returns an array of indices in topological order, or null if a cycle
    // is detected (callers fall back to source-order traversal).
    function topoSort(n, edges) {
        var inDeg = [];
        var adj = [];
        for (var i = 0; i < n; i++) { inDeg.push(0); adj.push([]); }
        for (var e = 0; e < edges.length; e++) {
            adj[edges[e].from].push(edges[e].to);
            inDeg[edges[e].to]++;
        }
        var queue = [];
        for (var k = 0; k < n; k++) { if (inDeg[k] === 0) { queue.push(k); } }
        var out = [];
        while (queue.length) {
            var u = queue.shift();
            out.push(u);
            for (var j = 0; j < adj[u].length; j++) {
                var v = adj[u][j];
                if (--inDeg[v] === 0) { queue.push(v); }
            }
        }
        return (out.length === n) ? out : null;
    }

    // v2.0.53: reachable-set from `src` along forward edges, used to decide
    // whether a downstream metabolite SHOULD have received a label (so we
    // can flag "lost at X" only for reachable metabolites; metabolites
    // outside the connected component are not "losses").
    function reachableFrom(src, n, edges) {
        var reach = [];
        for (var i = 0; i < n; i++) { reach.push(false); }
        reach[src] = true;
        var changed = true;
        while (changed) {
            changed = false;
            for (var e = 0; e < edges.length; e++) {
                if (reach[edges[e].from] && !reach[edges[e].to]) {
                    reach[edges[e].to] = true;
                    changed = true;
                }
            }
        }
        return reach;
    }

    // Trace every atom of the first metabolite forward through a list of
    // metabolite SMILES. Returns labels-per-metabolite + per-start-atom paths.
    //
    // v2.0.20-v2.0.52 default: linear chain — metaboliteSmiles[i] → [i+1].
    // v2.0.53 extension: pass `options.edges` = [{from, to}, …] to walk the
    // actual pathway DAG, so forks (aldolase: F1,6BP → GAP/DHAP) and
    // convergences (TPI: DHAP → GAP) route labels through the real
    // chemistry instead of inventing phantom "GAP → DHAP" / "DHAP → 1,3BPG"
    // links. When `options.edges` is absent, behaviour is byte-identical to
    // v2.0.52.
    //
    // Convergence policy: when atom `b` in metabolite `m` receives a label
    // from two incoming edges, the FIRST writer (in topological order) wins.
    // This is deterministic and matches typical biochemistry: an atom that
    // arrives at GAP via direct F1,6BP → GAP and again via F1,6BP → DHAP →
    // GAP carries the same start-atom label both times (the F1,6BP carbon
    // is the same atom either way), so the policy is invariant for the
    // 1-to-1 fold-back case. For ambiguous merges (different start atoms
    // converge on the same downstream atom), this is a tie-breaker.
    function tracePathway(metaboliteSmiles, options) {
        options = options || {};
        if (!metaboliteSmiles || metaboliteSmiles.length < 2) { return null; }

        var n = metaboliteSmiles.length;
        var mols = [];
        for (var i = 0; i < n; i++) {
            var mol = new Molecule();
            SmilesParser.parse(metaboliteSmiles[i], mol);
            mols.push(mol);
        }
        var n0 = mols[0].atoms.length;

        // v2.0.53: choose between linear chain (default) and explicit DAG.
        // v2.4.2: optional `options.edgeMaps[ei]` carries a PRECOMPUTED
        // fromToIndex map for edge ei (e.g. a workbench edge.aam.mapping that
        // was produced earlier via this module's own mapStepIndexMap). When an
        // entry is present it is used VERBATIM instead of re-running mapStep —
        // a pure cache, never a different index convention. Carried through the
        // same validity filter as edges so the alignment can't drift. Absent or
        // null entries fall back to the live mapStep, so default behaviour is
        // byte-identical to v2.0.53-v2.4.1.
        var edges;
        var edgeMapOverride = [];   // parallel to `edges` after filtering
        var rawEdgeMaps = (options && Array.isArray(options.edgeMaps)) ? options.edgeMaps : null;
        var usingDag = (options && Array.isArray(options.edges) && options.edges.length > 0);
        if (usingDag) {
            edges = [];
            for (var ei = 0; ei < options.edges.length; ei++) {
                var raw = options.edges[ei];
                if (!raw) { continue; }
                var fromIdx = (typeof raw.from === 'number') ? raw.from : raw[0];
                var toIdx   = (typeof raw.to   === 'number') ? raw.to   : raw[1];
                if (typeof fromIdx !== 'number' || typeof toIdx !== 'number') { continue; }
                if (fromIdx < 0 || fromIdx >= n || toIdx < 0 || toIdx >= n) { continue; }
                edges.push({ from: fromIdx, to: toIdx });
                edgeMapOverride.push(rawEdgeMaps ? (rawEdgeMaps[ei] || null) : null);
            }
            if (edges.length === 0) { usingDag = false; }
        }
        if (!usingDag) {
            edges = [];
            edgeMapOverride = [];
            for (var li = 0; li + 1 < n; li++) {
                edges.push({ from: li, to: li + 1 });
                edgeMapOverride.push(null);
            }
        }

        // v2.0.56: optional per-edge cofactor specs. `options.cofactors[e]`
        // = `{ in: [name…], out: [name…] }`. When present, the per-edge AAM
        // includes the cofactor SMILES in the balanced reaction, and we
        // capture the cofactor↔metabolite atom transfers as `cofactorLabels`
        // (parallel to the per-atom `labels` field).
        var cofactorSpecs = options.cofactors;
        var hasCofactors = Array.isArray(cofactorSpecs) && cofactorSpecs.length > 0;

        // Per-edge AAM step.
        var edgeMaps = [];
        var edgeCofactorXfer = [];   // v2.0.56: same length as edges
        var steps = [];
        for (var es = 0; es < edges.length; es++) {
            var edge = edges[es];
            var stepOpts = { shareMcsCache: options.shareMcsCache };
            var hasEdgeCofactors = hasCofactors && cofactorSpecs[es] &&
                (((cofactorSpecs[es].in || []).length) || ((cofactorSpecs[es].out || []).length));
            if (hasEdgeCofactors) {
                stepOpts.cofactorsIn = cofactorSpecs[es].in || [];
                stepOpts.cofactorsOut = cofactorSpecs[es].out || [];
            }
            // v2.4.2: prefer a precomputed map for this edge — but only when no
            // cofactor spec is in play (a cached metabolite->metabolite map
            // can't carry the cofactor↔metabolite transfers cofactorSpecs need).
            var override = (!hasEdgeCofactors && edgeMapOverride[es] &&
                typeof edgeMapOverride[es] === 'object') ? edgeMapOverride[es] : null;
            if (override) {
                edgeMaps.push(override);
                edgeCofactorXfer.push({ in: {}, out: {} });
                steps.push({
                    from: edge.from, to: edge.to,
                    status: 'cached', mapped: Object.keys(override).length,
                    cofactorXferIn: 0, cofactorXferOut: 0
                });
                continue;
            }
            var step = mapStep(metaboliteSmiles[edge.from], metaboliteSmiles[edge.to], stepOpts);
            edgeMaps.push(step.fromToIndex);
            edgeCofactorXfer.push({
                in:  step.cofactorInToProductMetabolite || {},
                out: step.reactantMetaboliteToCofactorOut || {}
            });
            steps.push({
                from: edge.from, to: edge.to,
                status: step.status, mapped: step.mappedCount,
                cofactorXferIn:  step.cofactorXferIn  || 0,
                cofactorXferOut: step.cofactorXferOut || 0
            });
        }

        // labels[metaboliteIndex][atomIndex] = array of originating start-atom
        // indices (empty when not reached). v2.0.53 multi-valued model: under
        // convergence (e.g. F1,6BP → DHAP → GAP AND F1,6BP → GAP direct), a
        // single downstream atom can have ancestors from BOTH incoming edges.
        // Pre-v2.0.53 callers that pattern-matched `labels[m][k] === N`
        // should switch to `labels[m][k].indexOf(N) !== -1`.
        var labels = [];
        for (var lm = 0; lm < n; lm++) {
            var bag = [];
            for (var la = 0; la < mols[lm].atoms.length; la++) { bag.push([]); }
            labels.push(bag);
        }
        // Start metabolite seeds the identity labels.
        for (var seed = 0; seed < n0; seed++) { labels[0][seed].push(seed); }

        // Propagate labels through edges in topological order so each
        // downstream atom is settled before its successors see it. Multi-
        // valued semantics: a label entering an atom via multiple incoming
        // edges is recorded ONCE (de-duplicated), but distinct labels from
        // distinct ancestor paths all accumulate.
        var topo = topoSort(n, edges);
        var order = topo;
        if (!order) {
            // Cycle detected — fall back to source-order traversal (preserves
            // v2.0.20-v2.0.52 behaviour on pathway templates that ever
            // accidentally introduce a back-edge). The trace is still
            // well-formed; only the convergence-tie ordering changes.
            order = [];
            for (var oi = 0; oi < n; oi++) { order.push(oi); }
        }
        // v2.0.56: parallel cofactor-origin labels per atom.
        // cofactorLabels[m][atomIdx] = array of {cofactor, atomIdx, step} records
        // recording every cofactor IN whose atom mapped onto this metabolite atom
        // across the DAG. Empty when no cofactor contribution.
        var cofactorLabels = [];
        for (var clm = 0; clm < n; clm++) {
            var cbag = [];
            for (var cla = 0; cla < mols[clm].atoms.length; cla++) { cbag.push([]); }
            cofactorLabels.push(cbag);
        }
        // Track atoms LEAVING a metabolite into a cofactor OUT — keyed by
        // (metaboliteIdx, atomIdx) -> [{cofactor, atomIdx, step}].
        var cofactorExits = [];
        for (var cem = 0; cem < n; cem++) {
            var ebag = [];
            for (var cea = 0; cea < mols[cem].atoms.length; cea++) { ebag.push([]); }
            cofactorExits.push(ebag);
        }

        for (var po = 0; po < order.length; po++) {
            var m = order[po];
            for (var ec = 0; ec < edges.length; ec++) {
                if (edges[ec].from !== m) { continue; }
                var dst = edges[ec].to;
                var mp = edgeMaps[ec];
                for (var srcId in mp) {
                    if (!mp.hasOwnProperty(srcId)) { continue; }
                    var dstIdx = mp[srcId];
                    var inboundLabels = labels[m][+srcId];
                    if (!inboundLabels || !inboundLabels.length) { continue; }
                    var outboundLabels = labels[dst][dstIdx];
                    for (var lab = 0; lab < inboundLabels.length; lab++) {
                        var L = inboundLabels[lab];
                        if (outboundLabels.indexOf(L) === -1) {
                            outboundLabels.push(L);
                        }
                    }
                }
                // v2.0.56: propagate inbound cofactor labels (so a cofactor
                // atom that landed in metabolite M and survives forward also
                // shows up in downstream metabolites). Atoms NOT in mp are
                // simply not propagated — matches the v2.0.53 semantic.
                for (var srcId2 in mp) {
                    if (!mp.hasOwnProperty(srcId2)) { continue; }
                    var dstIdx2 = mp[srcId2];
                    var inboundCofs = cofactorLabels[m][+srcId2];
                    if (!inboundCofs || !inboundCofs.length) { continue; }
                    var outboundCofs = cofactorLabels[dst][dstIdx2];
                    for (var lci = 0; lci < inboundCofs.length; lci++) {
                        var rec = inboundCofs[lci];
                        var dup = false;
                        for (var oci = 0; oci < outboundCofs.length; oci++) {
                            var ex = outboundCofs[oci];
                            if (ex.cofactor === rec.cofactor && ex.atomIdx === rec.atomIdx) {
                                dup = true; break;
                            }
                        }
                        if (!dup) { outboundCofs.push(rec); }
                    }
                }
                // v2.0.56: seed NEW cofactor labels for atoms that arrived
                // from a cofactor IN at this edge.
                var inXfer = edgeCofactorXfer[ec].in;
                for (var xferKey in inXfer) {
                    if (!inXfer.hasOwnProperty(xferKey)) { continue; }
                    var parts = xferKey.split(':');
                    var cofName = parts[0];
                    var cofAtomIdx = +parts[1];
                    var prodMetIdx = inXfer[xferKey];
                    var bagDst = cofactorLabels[dst][prodMetIdx];
                    if (!bagDst) { continue; }
                    // De-dup on (cofactor, atomIdx) so a cofactor that
                    // touches multiple incoming edges into dst isn't
                    // double-counted.
                    var seenCof = false;
                    for (var sc = 0; sc < bagDst.length; sc++) {
                        if (bagDst[sc].cofactor === cofName &&
                                bagDst[sc].atomIdx === cofAtomIdx) {
                            seenCof = true; break;
                        }
                    }
                    if (!seenCof) {
                        bagDst.push({
                            cofactor: cofName,
                            atomIdx: cofAtomIdx,
                            step: ec
                        });
                    }
                }
                // v2.0.56: record atoms LEAVING the reactant metabolite at
                // this edge (e.g. enolase loses one O of 2PG to H2O).
                var outXfer = edgeCofactorXfer[ec].out;
                for (var srcMetIdxStr in outXfer) {
                    if (!outXfer.hasOwnProperty(srcMetIdxStr)) { continue; }
                    var srcMetIdx = +srcMetIdxStr;
                    var exitTo = outXfer[srcMetIdxStr];
                    var exitParts = exitTo.split(':');
                    var bagExit = cofactorExits[m][srcMetIdx];
                    if (!bagExit) { continue; }
                    bagExit.push({
                        cofactor: exitParts[0],
                        atomIdx: +exitParts[1],
                        step: ec
                    });
                }
            }
        }

        // Per-start-atom trace: enumerate every metabolite where the label
        // appears (in topological order). With forks/convergences a single
        // start label can appear in multiple atoms of the same downstream
        // metabolite — we surface the FIRST hit per metabolite (deterministic
        // by parse-index order), which keeps `traces[t].path` a clean linear
        // sequence the UI can highlight. The full presence map is available
        // via the multi-valued `labels` field for callers that need every hit.
        //
        // v2.0.53 semantics for `reachesEnd` and `lostAtStep`:
        //   - `reachesEnd`  — label appears in the LAST metabolite of the input
        //                     array (metaboliteSmiles[n-1]). This is the
        //                     biology-meaningful question "does this atom
        //                     end up in the final product?". Under forks, an
        //                     atom that takes the OTHER path (skipping a
        //                     side-branch metabolite) still has reachesEnd =
        //                     true as long as it lands in the final node.
        //   - `lostAtStep`  — first reachable metabolite (in topo order) where
        //                     this label is ABSENT. Informational; can be set
        //                     even when reachesEnd is true (the label simply
        //                     bypassed that metabolite via the fork).
        var reach = reachableFrom(0, n, edges);
        var endIdx = n - 1;
        var traces = [];
        for (var t = 0; t < n0; t++) {
            var path = [{ metabolite: 0, atomIndex: t, element: mols[0].atoms[t].symbol }];
            var lostAt = null;
            var hitEnd = false;
            for (var oi2 = 0; oi2 < order.length; oi2++) {
                var mi = order[oi2];
                if (mi === 0) { continue; }
                if (!reach[mi]) { continue; }
                var arr = labels[mi], hit = -1;
                for (var ka = 0; ka < arr.length; ka++) {
                    if (arr[ka].indexOf(t) !== -1) { hit = ka; break; }
                }
                if (hit >= 0) {
                    path.push({ metabolite: mi, atomIndex: hit, element: mols[mi].atoms[hit].symbol });
                    if (mi === endIdx) { hitEnd = true; }
                } else if (lostAt === null) {
                    lostAt = mi;
                }
            }
            traces.push({
                start: t,
                element: mols[0].atoms[t].symbol,
                path: path,
                lostAtStep: lostAt,
                // reachesEnd is the "did this atom land in the final
                // metabolite" question. With v2.0.53 fork semantics this can
                // be true even when lostAt is non-null (the label took the
                // other branch and merged back via TPI / similar).
                reachesEnd: hitEnd
            });
        }

        return {
            metaboliteCount: n,
            startAtoms: n0,
            labels: labels,
            traces: traces,
            steps: steps,
            // v2.0.53: surface the edges actually used so callers can verify
            // the trace ran over the intended DAG (not a phantom linear chain).
            edges: edges,
            // v2.0.56: per-atom cofactor origin labels (atoms that ENTERED
            // via a cofactor IN) and cofactor exits (atoms that LEFT into
            // a cofactor OUT). Empty arrays when no cofactor specs provided.
            cofactorLabels: cofactorLabels,
            cofactorExits: cofactorExits,
            // Per-edge cofactor transfer maps (raw output, for power-callers).
            edgeCofactorXfer: edgeCofactorXfer
        };
    }

    // Compact interactive SVG renderer for the atom-trace strip. Each heavy
    // atom is wrapped in a `<g class="bime-trace-atom" data-atom-index="k">`
    // with a transparent hit-target circle (clickable) and a highlight ring
    // that becomes visible when the parent group carries `.is-traced`. Bonds
    // are drawn as simple lines (offset pairs for double, plus an axial line
    // for triple). Returns an inline SVG STRING ready to be embedded so atoms
    // remain DOM-addressable for click + highlight by the workbench. Independent
    // of ImageExport.toSVG — keeps that mature path untouched.
    function renderInteractiveSvg(mol, opts) {
        opts = opts || {};
        var W = opts.width || 180;
        var H = opts.height || 110;
        var PAD = opts.padding || 12;
        var NS = 'http://www.w3.org/2000/svg';
        if (!mol || !mol.atoms || !mol.atoms.length) {
            return '<svg xmlns="' + NS + '" width="' + W + '" height="' + H + '"></svg>';
        }
        var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            if (a.x < minX) { minX = a.x; }
            if (a.y < minY) { minY = a.y; }
            if (a.x > maxX) { maxX = a.x; }
            if (a.y > maxY) { maxY = a.y; }
        }
        var w = Math.max(maxX - minX, 1);
        var h = Math.max(maxY - minY, 1);
        var s = Math.min((W - 2 * PAD) / w, (H - 2 * PAD) / h);
        var ox = (W - w * s) / 2 - minX * s;
        var oy = (H - h * s) / 2 - minY * s;
        function tx(x) { return x * s + ox; }
        function ty(y) { return y * s + oy; }
        var parts = ['<svg xmlns="' + NS + '" width="' + W + '" height="' + H +
            '" viewBox="0 0 ' + W + ' ' + H + '" class="bime-trace-svg">'];
        // bonds
        for (var b = 0; b < mol.bonds.length; b++) {
            var bond = mol.bonds[b];
            var a1 = mol.getAtom(bond.atom1), a2 = mol.getAtom(bond.atom2);
            if (!a1 || !a2) { continue; }
            var x1 = tx(a1.x), y1 = ty(a1.y), x2 = tx(a2.x), y2 = ty(a2.y);
            var order = bond.order || 1;
            if (order === 2 || order === 3) {
                var dx = x2 - x1, dy = y2 - y1;
                var len = Math.sqrt(dx * dx + dy * dy) || 1;
                var px = -dy / len * 2.5, py = dx / len * 2.5;
                parts.push('<line x1="' + (x1 + px) + '" y1="' + (y1 + py) +
                    '" x2="' + (x2 + px) + '" y2="' + (y2 + py) +
                    '" stroke="#334155" stroke-width="1.4"/>');
                parts.push('<line x1="' + (x1 - px) + '" y1="' + (y1 - py) +
                    '" x2="' + (x2 - px) + '" y2="' + (y2 - py) +
                    '" stroke="#334155" stroke-width="1.4"/>');
                if (order === 3) {
                    parts.push('<line x1="' + x1 + '" y1="' + y1 + '" x2="' + x2 + '" y2="' + y2 +
                        '" stroke="#334155" stroke-width="1.4"/>');
                }
            } else {
                parts.push('<line x1="' + x1 + '" y1="' + y1 + '" x2="' + x2 + '" y2="' + y2 +
                    '" stroke="#334155" stroke-width="1.6"/>');
            }
        }
        // atom palette (matches BIME's existing element colours)
        var palette = { O: '#dc2626', N: '#2563eb', P: '#ea580c', S: '#ca8a04',
            F: '#16a34a', Cl: '#16a34a', Br: '#9333ea', I: '#7c3aed' };
        // v2.0.55 (security hardening): atom.symbol is interpolated into SVG
        // attribute and text content via string concatenation; reject anything
        // that isn't a plausible element symbol so a tampered Molecule object
        // can't smuggle markup. SmilesParser already canonicalises to known
        // element symbols, so this is a belt-and-braces assert — well-formed
        // inputs are always allowed; malformed inputs render as '?' instead
        // of producing raw markup. Documented for future call sites.
        var SYMBOL_RE = /^[A-Za-z][A-Za-z]?$/;
        for (var k = 0; k < mol.atoms.length; k++) {
            var atom = mol.atoms[k];
            var ax = tx(atom.x), ay = ty(atom.y);
            var symbol = (typeof atom.symbol === 'string' && SYMBOL_RE.test(atom.symbol))
                ? atom.symbol : '?';
            var showLabel = symbol !== 'C' || (atom.charge || 0) !== 0 || (atom.isotope || 0) !== 0;
            parts.push('<g class="bime-trace-atom" data-atom-index="' + k +
                '" data-element="' + symbol + '">');
            // transparent hit-target (clickable even when no visible label)
            parts.push('<circle cx="' + ax + '" cy="' + ay + '" r="9" fill="transparent" class="bime-trace-atom-hit"/>');
            // highlight ring (hidden until .is-traced on the parent group)
            parts.push('<circle cx="' + ax + '" cy="' + ay + '" r="8" class="bime-trace-atom-mark" fill="none" pointer-events="none"/>');
            if (showLabel) {
                var col = palette[symbol] || '#333';
                parts.push('<text x="' + ax + '" y="' + (ay + 4) + '" fill="' + col +
                    '" font-size="11" font-weight="700" text-anchor="middle" pointer-events="none">' +
                    symbol + '</text>');
            }
            parts.push('</g>');
        }
        parts.push('</svg>');
        return parts.join('');
    }

    // -----------------------------------------------------------------------
    // v2.0.54: MOIETY TRACING
    //
    // A moiety is a CONNECTED SET of start-metabolite atoms followed through
    // the pathway as a UNIT. The single-atom tracer (above) answers "where
    // does this atom go?". The moiety tracer answers "did this functional
    // group (phosphate, acyl, sugar ring) survive intact, and if not, where
    // did it break?".
    //
    // `deriveMoietyTrace(traceResult, moietyAtoms, options)` is a POST-
    // PROCESSOR over the result of `tracePathway`. It re-uses the trace's
    // multi-valued labels (v2.0.53) AND runs `RDT.deriveSubFragments` once
    // per pathway edge to determine, for every start-bond of the moiety,
    // which downstream metabolites still carry it. The strict aliveness
    // verdict at metabolite M is:
    //
    //   alive(M) =  every moiety atom has an image at M
    //               AND those images form ONE connected component in M
    //               AND every moiety start-bond is preserved at M
    //               (preserved iff bonded along at least one DAG path
    //                from start to M — "OR across paths, AND along an edge")
    //
    // The OR-across-paths join keeps moiety-alive monotone with respect to
    // label presence (which is itself a union under v2.0.53). Without it
    // an atom could claim presence-via-convergence while the connecting
    // bond was reported broken on a different path.
    //
    // Sub-fragment policy: a pair of atoms that lands in the same RDT
    // sub-fragment of one edge IS considered "still bonded" across that
    // edge — even if the bond order changed (C-C -> C=C is preserved).
    // This is the "rigid scaffold" reading and matches deriveSubFragments
    // semantics; bond-order changes are surfaced separately on the
    // result.bondChanges path, never as moiety breakage.
    //
    // Cofactor lineage (ATP / ADP / Pi / NAD+ / NADH) is OUT OF SCOPE for
    // v2.0.54 — a moiety transferred to a cofactor (e.g. PGK moves the
    // 1-phosphate of 1,3-BPG to ADP) registers as broken/absent from
    // that step onward. The v2.0.5x roadmap addresses cofactor following.
    //
    // Returns:
    //   {
    //       moietyAtoms:   [int],           // input, validated + sorted
    //       moietySize:    int,
    //       startBonds:    [{a, b}],        // bonds within the moiety at metabolite 0
    //       criterion:     'strict' | 'present',
    //       perMetabolite: [{               // one entry per metabolite
    //           metabolite:      int,
    //           reachable:       bool,
    //           present:         [int],     // start-atom indices whose label appears at m
    //           missing:         [int],     // start-atom indices NOT at m
    //           imageAtoms:      [int],     // atom indices in m carrying any moiety label (first image per ancestor)
    //           components:      [[int]],   // connected components of imageAtoms in m's bond graph
    //           preservedBonds:  [{a, b}],  // start-bonds still preserved at m
    //           brokenBonds:     [{a, b}],
    //           alive:           bool,
    //           status:          'intact' | 'fragmented' | 'partial-loss' | 'absent' | 'unreachable'
    //       }],
    //       breakAt:    int | null,         // first reachable metabolite where alive=false
    //       survivors:  [int],              // metabolite indices where alive=true
    //       error:      'not-connected' | 'empty-moiety' | undefined
    //   }
    // -----------------------------------------------------------------------

    // Build a Set-like "still-bonded" key map for one edge:
    //   key(a, b) -> true  iff reactant-side atoms at indices a, b land in
    //                       the SAME sub-fragment of the AAM result and are
    //                       therefore "still bonded" at the product side.
    // The keys use the same reactant-side flattened indexing as
    // `mapStep.fromToIndex` (so the moiety bond walk can compose without
    // an id <-> index round trip).
    function buildEdgePreservedSet(fromSmiles, toSmiles, options) {
        options = options || {};
        var rxn = SmilesParser.parse(fromSmiles + '>>' + toSmiles);
        var res = RDT.mapReaction(rxn, {
            requireStoichiometricBalance: false,
            shareMcsCache: options.shareMcsCache
        });
        var preserved = {};
        if (!res || !res.sides) { return preserved; }
        var rAtoms = flattenSideAtoms(res.sides.reactants);
        var idToR = {};
        for (var i = 0; i < rAtoms.length; i++) { idToR[rAtoms[i].id] = i; }
        var frags = (typeof RDT.deriveSubFragments === 'function')
            ? RDT.deriveSubFragments(res, { minSize: 2 }) : [];
        for (var f = 0; f < frags.length; f++) {
            var ids = frags[f].reactantAtomIds || [];
            if (ids.length < 2) { continue; }
            // Mark every unordered pair (a, b) within this fragment as preserved.
            // Chemistry sizes are tiny (a fragment seldom exceeds 20 atoms), so the
            // O(k^2) walk is cheap.
            for (var a = 0; a < ids.length; a++) {
                var ai = idToR[ids[a]];
                if (ai === undefined) { continue; }
                for (var b = a + 1; b < ids.length; b++) {
                    var bi = idToR[ids[b]];
                    if (bi === undefined) { continue; }
                    var lo = Math.min(ai, bi), hi = Math.max(ai, bi);
                    preserved[lo + '|' + hi] = true;
                }
            }
        }
        return preserved;
    }

    // Adjacency utility — returns boolean N×N (sparse) for a molecule's bond list.
    // Used to (a) extract moiety start-bonds and (b) compute connected
    // components of moiety image atoms in each downstream metabolite.
    function adjacencyOf(mol) {
        var adj = {};
        if (!mol || !mol.atoms || !mol.bonds) { return { adj: adj, idToIndex: {} }; }
        var idToIndex = {};
        for (var i = 0; i < mol.atoms.length; i++) { idToIndex[mol.atoms[i].id] = i; }
        for (var b = 0; b < mol.bonds.length; b++) {
            var bond = mol.bonds[b];
            var u = idToIndex[bond.atom1], v = idToIndex[bond.atom2];
            if (u === undefined || v === undefined) { continue; }
            if (!adj[u]) { adj[u] = {}; }
            if (!adj[v]) { adj[v] = {}; }
            adj[u][v] = true; adj[v][u] = true;
        }
        return { adj: adj, idToIndex: idToIndex };
    }

    // Connected components of a SUBSET of atoms under an adjacency map.
    function componentsOfSubset(subsetAtoms, adj) {
        var inSubset = {};
        for (var i = 0; i < subsetAtoms.length; i++) { inSubset[subsetAtoms[i]] = true; }
        var seen = {};
        var out = [];
        for (var s = 0; s < subsetAtoms.length; s++) {
            var start = subsetAtoms[s];
            if (seen[start]) { continue; }
            var stack = [start], comp = [];
            seen[start] = true;
            while (stack.length) {
                var u = stack.pop();
                comp.push(u);
                var n = adj[u] || {};
                for (var v in n) {
                    if (!n.hasOwnProperty(v)) { continue; }
                    var vi = +v;
                    if (inSubset[vi] && !seen[vi]) { seen[vi] = true; stack.push(vi); }
                }
            }
            comp.sort(function(a, b) { return a - b; });
            out.push(comp);
        }
        out.sort(function(a, b) { return a[0] - b[0]; });
        return out;
    }

    function deriveMoietyTrace(traceResult, moietyAtoms, options) {
        options = options || {};
        if (!traceResult || !traceResult.labels || !traceResult.edges) {
            return { error: 'invalid-trace', perMetabolite: [] };
        }
        if (!Array.isArray(moietyAtoms) || moietyAtoms.length === 0) {
            return { error: 'empty-moiety', perMetabolite: [] };
        }
        // De-dup + sort + validate.
        var seenAtom = {};
        var moiety = [];
        for (var ma = 0; ma < moietyAtoms.length; ma++) {
            var idx = moietyAtoms[ma] | 0;
            if (idx < 0) { continue; }
            if (!seenAtom[idx]) { seenAtom[idx] = true; moiety.push(idx); }
        }
        moiety.sort(function(a, b) { return a - b; });

        // metaboliteSmiles is required to derive bond preservation.
        var smiles = options.metaboliteSmiles;
        if (!Array.isArray(smiles) || smiles.length !== traceResult.metaboliteCount) {
            return { error: 'metabolite-smiles-required', perMetabolite: [] };
        }

        // Parse molecules ourselves (cheap) so we don't depend on caller-shared
        // Molecule arrays (preserves the no-DOM invariant for headless tests).
        var n = traceResult.metaboliteCount;
        var mols = [];
        for (var p = 0; p < n; p++) {
            var mol = new Molecule();
            SmilesParser.parse(smiles[p], mol);
            mols.push(mol);
        }

        // Start-metabolite adjacency + connectivity check.
        var startAdjInfo = adjacencyOf(mols[0]);
        var startAdj = startAdjInfo.adj;
        var startComps = componentsOfSubset(moiety, startAdj);
        if (startComps.length !== 1) {
            return {
                error: 'not-connected',
                moietyAtoms: moiety,
                moietySize: moiety.length,
                componentsAtStart: startComps,
                perMetabolite: []
            };
        }
        // Enumerate moiety start-bonds (every pair within the moiety that is
        // bonded in m=0). These are the bonds we'll track across the pathway.
        var startBonds = [];
        for (var ai = 0; ai < moiety.length; ai++) {
            var u = moiety[ai];
            var row = startAdj[u] || {};
            for (var aj = ai + 1; aj < moiety.length; aj++) {
                var v = moiety[aj];
                if (row[v]) { startBonds.push({ a: u, b: v }); }
            }
        }

        var criterion = (options.criterion === 'present') ? 'present' : 'strict';

        // Pre-compute "still-bonded" sets for every pathway edge.
        var edges = traceResult.edges;
        var edgePreserved = [];
        for (var e = 0; e < edges.length; e++) {
            edgePreserved.push(buildEdgePreservedSet(
                smiles[edges[e].from], smiles[edges[e].to], options));
        }

        // DAG walk: bondPreserved[m] = { 'a|b': true } for moiety start-bonds
        // (in start-metabolite index space) still alive at m.
        var bondPreserved = [];
        for (var bm = 0; bm < n; bm++) { bondPreserved.push({}); }
        for (var sb = 0; sb < startBonds.length; sb++) {
            var key0 = startBonds[sb].a + '|' + startBonds[sb].b;
            bondPreserved[0][key0] = true;
        }

        // Topological order over the trace's DAG. Falls back to source order
        // on cycle (matches tracePathway's policy).
        var topo = topoSort(n, edges);
        var order = topo;
        if (!order) {
            order = [];
            for (var oi = 0; oi < n; oi++) { order.push(oi); }
        }

        // For each topo step, propagate preserved bonds across outgoing edges.
        for (var po = 0; po < order.length; po++) {
            var src = order[po];
            // For each outgoing edge (src -> dst), check each currently-preserved
            // bond at src against the edge's still-bonded set.
            for (var ec = 0; ec < edges.length; ec++) {
                if (edges[ec].from !== src) { continue; }
                var dst = edges[ec].to;
                var preserved = edgePreserved[ec];
                for (var startBondKey in bondPreserved[src]) {
                    if (!bondPreserved[src].hasOwnProperty(startBondKey)) { continue; }
                    var parts = startBondKey.split('|');
                    var sa = +parts[0], sb2 = +parts[1];
                    // Find each end's image at the src metabolite via the
                    // existing per-start-atom path.
                    var imgA = imageOf(traceResult, sa, src);
                    var imgB = imageOf(traceResult, sb2, src);
                    if (imgA < 0 || imgB < 0) { continue; }
                    var lo = Math.min(imgA, imgB), hi = Math.max(imgA, imgB);
                    var edgeKey = lo + '|' + hi;
                    if (preserved[edgeKey]) {
                        // OR-union: bond alive at dst if alive along any
                        // incoming edge.
                        bondPreserved[dst][startBondKey] = true;
                    }
                }
            }
        }

        // Build per-metabolite report.
        var reach = reachableFrom(0, n, edges);
        var perMetabolite = [];
        var survivors = [];
        var breakAt = null;

        for (var m = 0; m < n; m++) {
            var entry = {
                metabolite: m,
                reachable: !!reach[m],
                present: [],
                missing: [],
                imageAtoms: [],
                components: [],
                preservedBonds: [],
                brokenBonds: [],
                alive: false,
                status: 'unreachable'
            };
            perMetabolite.push(entry);

            if (!reach[m]) { continue; }

            // Which start atoms have an image at m? Use the multi-valued labels.
            var atomLabels = traceResult.labels[m] || [];
            var imageByStart = {};
            for (var k = 0; k < atomLabels.length; k++) {
                var bag = atomLabels[k];
                if (!bag || !bag.length) { continue; }
                for (var bi = 0; bi < bag.length; bi++) {
                    var startIdx = bag[bi];
                    if (seenAtom[startIdx] && imageByStart[startIdx] === undefined) {
                        imageByStart[startIdx] = k;   // first image wins
                    }
                }
            }
            for (var mai = 0; mai < moiety.length; mai++) {
                var sa2 = moiety[mai];
                if (imageByStart[sa2] !== undefined) {
                    entry.present.push(sa2);
                    entry.imageAtoms.push(imageByStart[sa2]);
                } else {
                    entry.missing.push(sa2);
                }
            }

            // Connected components of the image atoms in m's bond graph.
            var mAdjInfo = adjacencyOf(mols[m]);
            entry.components = componentsOfSubset(entry.imageAtoms.slice(), mAdjInfo.adj);

            // Preserved / broken start-bonds at m.
            for (var bk = 0; bk < startBonds.length; bk++) {
                var sb3 = startBonds[bk];
                var keySB = sb3.a + '|' + sb3.b;
                if (bondPreserved[m][keySB]) {
                    entry.preservedBonds.push(sb3);
                } else {
                    entry.brokenBonds.push(sb3);
                }
            }

            // Aliveness verdict.
            var allPresent = (entry.missing.length === 0);
            var oneComponent = (entry.components.length === 1);
            var allBondsPreserved = (entry.brokenBonds.length === 0);

            if (criterion === 'strict') {
                entry.alive = allPresent && oneComponent && allBondsPreserved;
            } else {
                entry.alive = allPresent;
            }

            if (entry.alive) {
                entry.status = 'intact';
                survivors.push(m);
            } else if (!allPresent && entry.present.length === 0) {
                entry.status = 'absent';
            } else if (!allPresent) {
                entry.status = 'partial-loss';
            } else {
                entry.status = 'fragmented';
            }
            if (m !== 0 && breakAt === null && !entry.alive) {
                breakAt = m;
            }
        }

        return {
            moietyAtoms: moiety,
            moietySize: moiety.length,
            startBonds: startBonds,
            criterion: criterion,
            perMetabolite: perMetabolite,
            breakAt: breakAt,
            survivors: survivors,
            version: '2.0.54'
        };
    }

    // Helper: find the first-image atom index for start atom sa at metabolite m.
    // -1 if no image. Mirrors the per-atom `traces[].path` first-hit policy.
    function imageOf(traceResult, sa, m) {
        if (m === 0) { return sa; }
        var bag = traceResult.labels && traceResult.labels[m];
        if (!bag) { return -1; }
        for (var k = 0; k < bag.length; k++) {
            var labs = bag[k];
            if (!labs) { continue; }
            for (var li = 0; li < labs.length; li++) {
                if (labs[li] === sa) { return k; }
            }
        }
        return -1;
    }

    global.AtomTrace = {
        tracePathway: tracePathway,
        mapStep: mapStep,
        // v2.4.2: shared RDT-id -> positional-index translation. Exposed so
        // workbench's linkImportedStepAAM stores edge.aam.mapping in EXACTLY
        // the same fromToIndex positional space tracePathway consumes (never
        // RDT internal atom-ids).
        mapStepIndexMap: mapStepIndexMap,
        renderInteractiveSvg: renderInteractiveSvg,
        // v2.0.54: connected-subgraph moiety tracing.
        deriveMoietyTrace: deriveMoietyTrace,
        // v2.0.53: exposed for tests / external callers that want to verify
        // the topological ordering used by tracePathway.
        _topoSort: topoSort,
        _reachableFrom: reachableFrom,
        // v2.0.54: exposed for tests.
        _buildEdgePreservedSet: buildEdgePreservedSet,
        _adjacencyOf: adjacencyOf,
        _componentsOfSubset: componentsOfSubset,
        // v2.0.56: cofactor lineage primitives.
        mapStepWithCofactors: mapStepWithCofactors,
        parseCofactorString: parseCofactorString,
        lookupCofactorSmiles: lookupCofactorSmiles,
        COFACTOR_SMILES: COFACTOR_SMILES,
        // Release-tracker stamp: surfaced by `bime version` alongside
        // RDT.version / ToolbarPrefs.version, so it tracks the BIME release
        // (pinned in tests/test_v2_0_70_stamp_version_lockstep.js). NOTE:
        // this is distinct from the moiety-trace RESULT-format version
        // returned by deriveMoietyTrace() (frozen at 2.0.54 — that records
        // the data-format revision, not the release).
        version: '3.0.3'
    };
})(typeof window !== 'undefined' ? window : this);
