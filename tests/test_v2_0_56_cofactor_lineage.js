/**
 * tests/test_v2_0_56_cofactor_lineage.js — cofactor lineage scaffolding
 * (v2.0.56).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Prior to v2.0.56, `AtomTrace.mapStep` parsed only `from + ">>" + to`,
 * so the phosphate that lands on glucose-6-phosphate at step 1 of
 * glycolysis (donated by ATP's γ-phosphate) was invisible to the
 * tracer. v2.0.56 adds:
 *
 *   - A small `COFACTOR_SMILES` table (Pi, H2O, CO2, H+, ATP/ADP/AMP,
 *     GTP/GDP, NAD+/NADH, NADP+/NADPH, …).
 *   - `parseCofactorString('ATP → ADP')` → `{in:['ATP'], out:['ADP']}`.
 *   - `mapStepWithCofactors(from, to, cofactorsIn, cofactorsOut, opts)`
 *     that builds a balanced reaction including the cofactors and
 *     returns four cross-side maps:
 *       fromToIndex                      — metabolite → metabolite (legacy)
 *       cofactorInToProductMetabolite    — cofactor IN → product metabolite
 *       reactantMetaboliteToCofactorOut  — reactant metabolite → cofactor OUT
 *       cofactorInToCofactorOut          — cofactor IN → cofactor OUT
 *   - `tracePathway(.., {cofactors:[{in:..,out:..}]})` propagates the
 *     cofactor lineage in parallel with the v2.0.53 multi-valued labels,
 *     producing `result.cofactorLabels[m][k] = [{cofactor, atomIdx,
 *     step}, ...]` and `result.cofactorExits[m][k] = [{cofactor,
 *     atomIdx, step}, ...]`.
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');
require('../editor/AtomTrace.js');

var runner = shim.makeRunner('Cofactor lineage (v2.0.56)');
var test = runner.test;

console.log('Cofactor lineage (v2.0.56)');

// ---------------------------------------------------------------------------
// A. cofactor SMILES + lookup
// ---------------------------------------------------------------------------
test('A1 COFACTOR_SMILES has the eight common metabolic cofactors', function () {
    var must = ['atp', 'adp', 'amp', 'pi', 'h2o', 'co2', 'h+', 'nad+'];
    for (var i = 0; i < must.length; i++) {
        assert.ok(typeof AtomTrace.COFACTOR_SMILES[must[i]] === 'string' &&
                  AtomTrace.COFACTOR_SMILES[must[i]].length > 0,
                  'COFACTOR_SMILES[' + must[i] + '] is non-empty');
    }
});

test('A2 lookupCofactorSmiles handles canonical names + aliases case-insensitively', function () {
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('Pi'), 'OP(O)(O)=O');
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('phosphate'), 'OP(O)(O)=O');
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('Orthophosphate'), 'OP(O)(O)=O');
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('H2O'), 'O');
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('water'), 'O');
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('co2'), 'O=C=O');
    // ATP SMILES contains 3 phosphate residues, so its string contains 'P'.
    assert.ok(AtomTrace.lookupCofactorSmiles('ATP').indexOf('P') >= 0,
        'ATP SMILES contains phosphate atoms');
    assert.strictEqual(AtomTrace.lookupCofactorSmiles('unknown-cofactor'), null);
});

// ---------------------------------------------------------------------------
// B. parseCofactorString
// ---------------------------------------------------------------------------
test('B1 parseCofactorString handles "ATP → ADP" (the common arrow form)', function () {
    var r = AtomTrace.parseCofactorString('ATP → ADP');
    assert.deepStrictEqual(r, { in: ['ATP'], out: ['ADP'] });
});

test('B2 parseCofactorString handles multi-token "NAD+ + Pi → NADH"', function () {
    var r = AtomTrace.parseCofactorString('NAD+ + Pi → NADH');
    assert.deepStrictEqual(r, { in: ['NAD+', 'Pi'], out: ['NADH'] });
});

test('B3 parseCofactorString handles leading "+ X + Y → Z" (TCA step 1 style)', function () {
    var r = AtomTrace.parseCofactorString('+ Acetyl-CoA + H2O → CoA-SH');
    assert.deepStrictEqual(r, { in: ['Acetyl-CoA', 'H2O'], out: ['CoA-SH'] });
});

test('B4 parseCofactorString handles one-sided "− H2O" loss', function () {
    var r = AtomTrace.parseCofactorString('− H2O');
    assert.deepStrictEqual(r, { in: [], out: ['H2O'] });
});

test('B5 parseCofactorString handles one-sided "+ H2O" gain', function () {
    var r = AtomTrace.parseCofactorString('+ H2O');
    assert.deepStrictEqual(r, { in: ['H2O'], out: [] });
});

test('B6 parseCofactorString returns empty on empty/non-string', function () {
    assert.deepStrictEqual(AtomTrace.parseCofactorString(''),    { in: [], out: [] });
    assert.deepStrictEqual(AtomTrace.parseCofactorString(null),  { in: [], out: [] });
    assert.deepStrictEqual(AtomTrace.parseCofactorString(undefined), { in: [], out: [] });
});

// ---------------------------------------------------------------------------
// C. mapStepWithCofactors — Pi → 1,3BPG (GAPDH step)
// ---------------------------------------------------------------------------
test('C1 GAPDH (GAP + Pi → 1,3BPG): Pi atoms appear in the product metabolite', function () {
    // Glyceraldehyde 3-phosphate >> 1,3-Bisphosphoglycerate. The new
    // phosphate on 1,3BPG comes from Pi (free phosphate). Under the
    // cofactor-aware AAM we expect ≥ 1 atom transferred from Pi to 1,3BPG.
    var GAP = 'O=C[C@@H](O)COP(O)(O)=O';
    var BPG = 'O[C@@H](COP(O)(O)=O)C(=O)OP(O)(O)=O';
    var step = AtomTrace.mapStepWithCofactors(GAP, BPG, ['Pi'], [], { shareMcsCache: false });
    assert.ok(step.cofactorXferIn >= 1,
        'at least 1 Pi atom maps into 1,3BPG, got ' + step.cofactorXferIn);
    // Verify the transferred atom names mention Pi.
    var sawPi = false;
    for (var key in step.cofactorInToProductMetabolite) {
        if (key.indexOf('Pi:') === 0) { sawPi = true; break; }
    }
    assert.ok(sawPi, 'cofactorInToProductMetabolite has a Pi: key');
});

// ---------------------------------------------------------------------------
// D. mapStepWithCofactors — enolase loses H2O
// ---------------------------------------------------------------------------
test('D1 enolase (2PG → PEP + H2O): a reactant 2PG atom leaves into H2O', function () {
    var PG2 = 'OC[C@@H](OP(O)(O)=O)C(O)=O';
    var PEP = 'OC(=O)C(=C)OP(O)(O)=O';
    var step = AtomTrace.mapStepWithCofactors(PG2, PEP, [], ['H2O'], { shareMcsCache: false });
    assert.ok(step.cofactorXferOut >= 1,
        'at least 1 reactant atom leaves into H2O, got ' + step.cofactorXferOut);
    var sawH2O = false;
    for (var srcIdx in step.reactantMetaboliteToCofactorOut) {
        if (step.reactantMetaboliteToCofactorOut[srcIdx].indexOf('H2O:') === 0) {
            sawH2O = true; break;
        }
    }
    assert.ok(sawH2O, 'reactantMetaboliteToCofactorOut has an H2O: target');
});

// ---------------------------------------------------------------------------
// E. mapStep backward-compat — no cofactor opts ≡ v2.0.55 behaviour
// ---------------------------------------------------------------------------
test('E1 mapStep without cofactor opts returns the v2.0.55 shape', function () {
    var step = AtomTrace.mapStep('CCO', 'CC=O');
    assert.ok(step.fromToIndex);
    assert.ok(Array.isArray(step.reactantElements));
    assert.ok(Array.isArray(step.productElements));
    assert.strictEqual(typeof step.mappedCount, 'number');
    assert.strictEqual(typeof step.status, 'string');
    // No cofactor fields when no cofactor opts.
    assert.strictEqual(step.cofactorXferIn, undefined);
});

test('E2 mapStep WITH cofactor opts forwards to mapStepWithCofactors', function () {
    var step = AtomTrace.mapStep('O=C[C@@H](O)COP(O)(O)=O',
                                  'O[C@@H](COP(O)(O)=O)C(=O)OP(O)(O)=O',
                                  { cofactorsIn: ['Pi'] });
    assert.ok(step.cofactorXferIn >= 1, 'cofactor transfer surfaced via legacy mapStep');
});

// ---------------------------------------------------------------------------
// F. tracePathway end-to-end with cofactors
// ---------------------------------------------------------------------------
test('F1 tracePathway with cofactors fills cofactorLabels and cofactorExits', function () {
    // Three-step micro-pathway: GAP -> 1,3BPG (+Pi) -> 3PG (-ATP-like phosphate)
    // Using the same SMILES as the glycolysis fragment.
    var GAP = 'O=C[C@@H](O)COP(O)(O)=O';
    var BPG = 'O[C@@H](COP(O)(O)=O)C(=O)OP(O)(O)=O';
    var PG3 = 'O[C@@H](COP(O)(O)=O)C(O)=O';
    var smi = [GAP, BPG, PG3];
    var edges = [{ from: 0, to: 1 }, { from: 1, to: 2 }];
    var cofactors = [
        { in: ['Pi'], out: [] },   // GAP + Pi -> 1,3BPG
        { in: [],     out: [] }    // 1,3BPG -> 3PG (no cofactor used here)
    ];
    var r = AtomTrace.tracePathway(smi, { edges: edges, cofactors: cofactors });
    assert.ok(r && Array.isArray(r.cofactorLabels), 'cofactorLabels populated');
    assert.strictEqual(r.cofactorLabels.length, 3, 'one per metabolite');
    // 1,3BPG should have at least one atom whose cofactorLabels contains a 'Pi' record.
    var bpgAtoms = r.cofactorLabels[1];
    var sawPi = false;
    for (var k = 0; k < bpgAtoms.length; k++) {
        for (var li = 0; li < bpgAtoms[k].length; li++) {
            if (bpgAtoms[k][li].cofactor === 'Pi') { sawPi = true; break; }
        }
        if (sawPi) { break; }
    }
    assert.ok(sawPi, 'at least one 1,3BPG atom carries a Pi cofactorLabel');
});

test('F2 cofactor labels propagate downstream (Pi-origin atoms in BPG persist into 3PG)', function () {
    var GAP = 'O=C[C@@H](O)COP(O)(O)=O';
    var BPG = 'O[C@@H](COP(O)(O)=O)C(=O)OP(O)(O)=O';
    var PG3 = 'O[C@@H](COP(O)(O)=O)C(O)=O';
    var smi = [GAP, BPG, PG3];
    var edges = [{ from: 0, to: 1 }, { from: 1, to: 2 }];
    var cofactors = [
        { in: ['Pi'], out: [] },
        { in: [],     out: [] }
    ];
    var r = AtomTrace.tracePathway(smi, { edges: edges, cofactors: cofactors });
    // 3PG should ALSO have a Pi-origin atom (forward-propagated from 1,3BPG).
    var pg3Atoms = r.cofactorLabels[2];
    var sawPi = false;
    for (var k = 0; k < pg3Atoms.length; k++) {
        for (var li = 0; li < pg3Atoms[k].length; li++) {
            if (pg3Atoms[k][li].cofactor === 'Pi') { sawPi = true; break; }
        }
        if (sawPi) { break; }
    }
    assert.ok(sawPi, '3PG inherits a Pi-origin cofactorLabel via forward propagation');
});

test('F3 tracePathway without cofactors leaves cofactorLabels EMPTY (backward compat)', function () {
    var r = AtomTrace.tracePathway(['CCO', 'CC=O', 'CC(=O)O']);
    assert.ok(Array.isArray(r.cofactorLabels));
    assert.strictEqual(r.cofactorLabels.length, 3);
    for (var m = 0; m < r.cofactorLabels.length; m++) {
        for (var k = 0; k < r.cofactorLabels[m].length; k++) {
            assert.strictEqual(r.cofactorLabels[m][k].length, 0,
                'm=' + m + ' k=' + k + ' has no cofactor labels');
        }
    }
});

// ---------------------------------------------------------------------------
// G. version stamp
// ---------------------------------------------------------------------------
test('G1 AtomTrace.version is at least 2.0.56 (forward-compatible check)', function () {
    // The v2.0.56 cofactor-lineage feature surface is what we're testing;
    // later releases may bump AtomTrace.version (v2.0.59 wired the UI etc.).
    assert.ok(typeof AtomTrace.version === 'string');
    var parts = AtomTrace.version.split('.').map(Number);
    var atLeast = (parts[0] > 2) || (parts[0] === 2 && (parts[1] > 0 || parts[2] >= 56));
    assert.ok(atLeast, 'AtomTrace.version >= 2.0.56, got ' + AtomTrace.version);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
