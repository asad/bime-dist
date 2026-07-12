/**
 * tests/test_v2_0_54_moiety_trace.js — connected-subgraph moiety tracing
 * (v2.0.54).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.54 adds `AtomTrace.deriveMoietyTrace(traceResult, moietyAtoms,
 * options)` — a post-processor over a v2.0.53 trace result. The user
 * supplies a CONNECTED set of start-metabolite atom indices (a "moiety"
 * — phosphate group, acyl, sugar ring) and the post-processor reports
 * for each downstream metabolite whether the moiety survived intact
 * (atoms present + connected + every start-bond preserved).
 *
 * Bond preservation is computed per-edge via `RDT.deriveSubFragments`
 * and joined across DAG paths with OR semantics: a bond is preserved
 * at M if it survives along at least ONE 0 → M path. This is the
 * monotone join consistent with v2.0.53's multi-valued labels.
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');
require('../editor/AtomTrace.js');
require('../editor/MetaboliteLibrary.js');

var runner = shim.makeRunner('Moiety trace (v2.0.54)');
var test = runner.test;

console.log('Moiety trace (v2.0.54)');

// Glycolysis setup re-used across tests ---------------------------------

function glycolysisSetup() {
    var pw = MetaboliteLibrary.pathways.glycolysis;
    var order = [], seen = {};
    pw.steps.forEach(function (s) {
        [s.from, s.to].forEach(function (n) {
            if (!seen[n]) { seen[n] = true; order.push(n); }
        });
    });
    var nameToIdx = {};
    order.forEach(function (n, i) { nameToIdx[n.toLowerCase()] = i; });
    var smiles = order.map(function (n) {
        return (MetaboliteLibrary.find(n) || {}).smiles || '';
    });
    var edges = pw.steps.map(function (s) {
        return {
            from: nameToIdx[s.from.toLowerCase()],
            to:   nameToIdx[s.to.toLowerCase()]
        };
    });
    return { order: order, nameToIdx: nameToIdx, smiles: smiles, edges: edges };
}

// Identify the phosphate group (P + adjacent Os) inside a parsed molecule.
function findPhosphateAtomIndices(smiles) {
    var mol = new Molecule();
    SmilesParser.parse(smiles, mol);
    var idToIdx = {};
    for (var i = 0; i < mol.atoms.length; i++) { idToIdx[mol.atoms[i].id] = i; }
    var pIdxs = [];
    for (var k = 0; k < mol.atoms.length; k++) {
        if (mol.atoms[k].symbol === 'P') { pIdxs.push(k); }
    }
    if (pIdxs.length !== 1) { return null; }
    var P = pIdxs[0];
    var pAtomId = mol.atoms[P].id;
    var atoms = [P];
    for (var b = 0; b < mol.bonds.length; b++) {
        var bond = mol.bonds[b];
        if (bond.atom1 === pAtomId && mol.atoms[idToIdx[bond.atom2]].symbol === 'O') {
            atoms.push(idToIdx[bond.atom2]);
        } else if (bond.atom2 === pAtomId && mol.atoms[idToIdx[bond.atom1]].symbol === 'O') {
            atoms.push(idToIdx[bond.atom1]);
        }
    }
    atoms.sort(function (a, b) { return a - b; });
    return atoms;
}

// ---------------------------------------------------------------------------
// A. API contract
// ---------------------------------------------------------------------------
test('A1 deriveMoietyTrace returns version 2.0.54 + mirrors moietyAtoms + emits startBonds', function () {
    var smi = ['CCO', 'CC=O'];
    var tr = AtomTrace.tracePathway(smi);
    var r = AtomTrace.deriveMoietyTrace(tr, [0, 1], { metaboliteSmiles: smi });
    assert.strictEqual(r.version, '2.0.54');
    assert.deepStrictEqual(r.moietyAtoms, [0, 1]);
    assert.strictEqual(r.moietySize, 2);
    assert.deepStrictEqual(r.startBonds, [{ a: 0, b: 1 }]);
    assert.strictEqual(r.criterion, 'strict');
    assert.strictEqual(r.perMetabolite.length, 2);
});

// ---------------------------------------------------------------------------
// B. Validation
// ---------------------------------------------------------------------------
test('B1 deriveMoietyTrace rejects an empty moiety with error="empty-moiety"', function () {
    var tr = AtomTrace.tracePathway(['CCO', 'CC=O']);
    var r = AtomTrace.deriveMoietyTrace(tr, [], { metaboliteSmiles: ['CCO', 'CC=O'] });
    assert.strictEqual(r.error, 'empty-moiety');
    assert.deepStrictEqual(r.perMetabolite, []);
});

test('B2 deriveMoietyTrace rejects a disconnected moiety with error="not-connected"', function () {
    // CCO atoms: 0=C, 1=C, 2=O. Atoms 0 and 2 are NOT directly bonded
    // (the C-O bond is 1-2, not 0-2), so {0, 2} is disconnected.
    var smi = ['CCO', 'CC=O'];
    var tr = AtomTrace.tracePathway(smi);
    var r = AtomTrace.deriveMoietyTrace(tr, [0, 2], { metaboliteSmiles: smi });
    assert.strictEqual(r.error, 'not-connected');
    // Reports the components at the start to help the caller debug.
    assert.ok(Array.isArray(r.componentsAtStart));
    assert.strictEqual(r.componentsAtStart.length, 2);
});

test('B3 deriveMoietyTrace requires metaboliteSmiles in options', function () {
    var tr = AtomTrace.tracePathway(['CCO', 'CC=O']);
    var r = AtomTrace.deriveMoietyTrace(tr, [0, 1], {});
    assert.strictEqual(r.error, 'metabolite-smiles-required');
});

// ---------------------------------------------------------------------------
// C. Linear-chain behaviour — bond preservation across simple oxidations
// ---------------------------------------------------------------------------
test('C1 C-C bond of ethanol → acetaldehyde → acetic acid stays preserved end-to-end', function () {
    // The C-C scaffold is rigid across all three; only the C-O changes.
    var smi = ['CCO', 'CC=O', 'CC(=O)O'];
    var tr = AtomTrace.tracePathway(smi);
    var r = AtomTrace.deriveMoietyTrace(tr, [0, 1], { metaboliteSmiles: smi });
    assert.strictEqual(r.breakAt, null, 'C-C scaffold never breaks');
    assert.strictEqual(r.survivors.length, 3, 'alive at all three metabolites');
    r.perMetabolite.forEach(function (e, i) {
        assert.strictEqual(e.alive, true, 'metabolite ' + i + ' alive');
        assert.strictEqual(e.status, 'intact');
    });
});

test('C2 bond-order change (C-O -> C=O) is NOT counted as moiety breakage', function () {
    // Ethanol C-O atoms = [1, 2]. They stay bonded across CC=O (order 2)
    // and CC(=O)O (order 2). RDT's "rigid scaffold" reading preserves
    // them — this matches deriveSubFragments semantics.
    var smi = ['CCO', 'CC=O', 'CC(=O)O'];
    var tr = AtomTrace.tracePathway(smi);
    var r = AtomTrace.deriveMoietyTrace(tr, [1, 2], { metaboliteSmiles: smi });
    assert.strictEqual(r.breakAt, null, 'order change is not breakage');
    assert.strictEqual(r.survivors.length, 3);
});

// ---------------------------------------------------------------------------
// D. Multi-valued labels — a 1-atom moiety reduces to single-atom tracing
// ---------------------------------------------------------------------------
test('D1 single-atom moiety reduces to single-atom trace presence', function () {
    var smi = ['CCO', 'CC=O', 'CC(=O)O'];
    var tr = AtomTrace.tracePathway(smi);
    var r = AtomTrace.deriveMoietyTrace(tr, [0], { metaboliteSmiles: smi });
    // No start-bonds (only one atom in the moiety).
    assert.deepStrictEqual(r.startBonds, []);
    // alive iff atom is present at m. tr.traces[0] is the single-atom view.
    for (var m = 0; m < 3; m++) {
        var presentInTrace = tr.traces[0].path.some(function (p) { return p.metabolite === m; });
        assert.strictEqual(r.perMetabolite[m].alive, presentInTrace,
            'moiety alive matches single-atom presence at m=' + m);
    }
});

// ---------------------------------------------------------------------------
// E. Criterion: strict vs present
// ---------------------------------------------------------------------------
test('E1 both criteria are honoured; phosphate is intact through PEP, lost at pyruvate', function () {
    // v2.4.17: the AAM now maps glycolysis correctly (the flex bond-order
    // candidate lets the MCS span reaction-centre order changes), so the
    // phosphate ester is preserved ALL the way to PEP — enolase (2-PG -> PEP)
    // dehydrates without touching the phosphate — and only leaves at the
    // PEP -> pyruvate step. Before the fix the mapper spuriously reported a
    // broken phosphate bond at PEP; that was a mapping artefact, not chemistry.
    // With a correct mapping the group is either fully intact or fully gone, so
    // strict and present agree on glycolysis (the criterion field is still
    // applied and reported).
    var setup = glycolysisSetup();
    // Start from G6P (index 1) so the phosphate exists at start.
    var sliceSmi = setup.smiles.slice(1);
    var sliceEdges = [];
    setup.edges.forEach(function (e) {
        if (e.from >= 1 && e.to >= 1) {
            sliceEdges.push({ from: e.from - 1, to: e.to - 1 });
        }
    });
    var tr = AtomTrace.tracePathway(sliceSmi, { edges: sliceEdges });
    var moiety = findPhosphateAtomIndices(sliceSmi[0]);
    var strict = AtomTrace.deriveMoietyTrace(tr, moiety, {
        metaboliteSmiles: sliceSmi, criterion: 'strict'
    });
    var present = AtomTrace.deriveMoietyTrace(tr, moiety, {
        metaboliteSmiles: sliceSmi, criterion: 'present'
    });
    assert.strictEqual(strict.criterion, 'strict');
    assert.strictEqual(present.criterion, 'present');
    // PEP (relative index 8): all 5 phosphate atoms present, group intact.
    var strictPep = strict.perMetabolite[8];
    assert.strictEqual(strictPep.present.length, 5,
        'all 5 phosphate atoms present at PEP');
    assert.strictEqual(strictPep.alive, true,
        'phosphate is strict-intact at PEP (enolase leaves the ester untouched)');
    // Pyruvate (relative index 9): the phosphate has left (transferred to ADP);
    // both criteria agree it is no longer alive because the atoms are gone.
    assert.strictEqual(strict.perMetabolite[9].alive, false,
        'strict: phosphate gone at pyruvate');
    assert.strictEqual(present.perMetabolite[9].alive, false,
        'present: phosphate atoms no longer present at pyruvate');
});

// ---------------------------------------------------------------------------
// F. Glycolysis: phosphate moiety survives ~8 steps then fragments
// ---------------------------------------------------------------------------
test('F1 glycolysis: G6P phosphate (PO3) survives strict-intact through ~8 metabolites', function () {
    var setup = glycolysisSetup();
    var sliceSmi = setup.smiles.slice(1);
    var sliceEdges = [];
    setup.edges.forEach(function (e) {
        if (e.from >= 1 && e.to >= 1) {
            sliceEdges.push({ from: e.from - 1, to: e.to - 1 });
        }
    });
    var tr = AtomTrace.tracePathway(sliceSmi, { edges: sliceEdges });
    var moiety = findPhosphateAtomIndices(sliceSmi[0]);
    assert.ok(moiety && moiety.length === 5, 'phosphate moiety has 5 atoms (P + 4 O)');
    var r = AtomTrace.deriveMoietyTrace(tr, moiety, { metaboliteSmiles: sliceSmi });
    // We expect: intact through G6P, F6P, F1,6BP, GAP, DHAP, 1,3BPG, 3PG, 2PG.
    // Breakage at PEP / Pyr (enolase dehydration + PK transfer).
    assert.ok(r.survivors.length >= 7,
        'phosphate moiety survives at least 7 metabolites (got ' + r.survivors.length + ')');
    assert.ok(r.breakAt !== null, 'eventually breaks');
});

test('F2 glycolysis: at pyruvate the phosphate moiety is "partial-loss" (most atoms gone)', function () {
    var setup = glycolysisSetup();
    var sliceSmi = setup.smiles.slice(1);
    var sliceEdges = [];
    setup.edges.forEach(function (e) {
        if (e.from >= 1 && e.to >= 1) {
            sliceEdges.push({ from: e.from - 1, to: e.to - 1 });
        }
    });
    var tr = AtomTrace.tracePathway(sliceSmi, { edges: sliceEdges });
    var moiety = findPhosphateAtomIndices(sliceSmi[0]);
    var r = AtomTrace.deriveMoietyTrace(tr, moiety, { metaboliteSmiles: sliceSmi });
    // Pyruvate is the LAST metabolite in the slice.
    var lastEntry = r.perMetabolite[r.perMetabolite.length - 1];
    assert.strictEqual(lastEntry.alive, false);
    assert.ok(lastEntry.missing.length >= 1, 'some phosphate atoms have left by pyruvate');
});

// ---------------------------------------------------------------------------
// G. Per-metabolite report shape
// ---------------------------------------------------------------------------
test('G1 perMetabolite entries carry the documented fields', function () {
    var smi = ['CCO', 'CC=O'];
    var tr = AtomTrace.tracePathway(smi);
    var r = AtomTrace.deriveMoietyTrace(tr, [0, 1], { metaboliteSmiles: smi });
    var e = r.perMetabolite[0];
    var keys = Object.keys(e).sort();
    var expected = ['alive', 'brokenBonds', 'components', 'imageAtoms',
                    'metabolite', 'missing', 'preservedBonds', 'present',
                    'reachable', 'status'].sort();
    assert.deepStrictEqual(keys, expected);
});

// ---------------------------------------------------------------------------
// H. Backward-compat with v2.0.20 / v2.0.53
// ---------------------------------------------------------------------------
test('H1 tracePathway output unchanged for non-moiety callers (v2.0.53 backward-compat)', function () {
    // The v2.0.53 baseline test for the linear chain should still produce
    // identical traces — adding deriveMoietyTrace doesn't perturb tracePathway.
    var r = AtomTrace.tracePathway(['CCO', 'CC=O', 'CC(=O)O']);
    assert.strictEqual(r.metaboliteCount, 3);
    assert.strictEqual(r.startAtoms, 3);
    assert.deepStrictEqual(r.edges, [{ from: 0, to: 1 }, { from: 1, to: 2 }]);
});

test('H2 AtomTrace.version is at least 2.0.54 (forward-compatible check)', function () {
    // The v2.0.54 feature surface is what we're testing; later patches may
    // bump AtomTrace.version while keeping the feature intact.
    assert.ok(typeof AtomTrace.version === 'string');
    var parts = AtomTrace.version.split('.').map(Number);
    var atLeast = (parts[0] > 2) || (parts[0] === 2 && (parts[1] > 0 || parts[2] >= 54));
    assert.ok(atLeast, 'AtomTrace.version >= 2.0.54, got ' + AtomTrace.version);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
