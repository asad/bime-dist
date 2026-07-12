/**
 * tests/test_cip_ring_stereocentres.js — CIP R/S correctness on ring-embedded
 * stereocentres (regression armour).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists
 * ---------------
 * A review once claimed BIME emitted inverted R/S on 1,2-disubstituted
 * alicyclic rings. That claim was traced to mislabelled structure/identifier
 * pairs, NOT an engine fault: the engine's descriptors match the modern
 * IUPAC-2013 reference (cross-checked against an independent CIP labeller and
 * PubChem isomeric SMILES) on every genuine-R/S ring stereocentre tested here.
 *
 * This suite PINS that correctness so a future well-meaning "fix" to the
 * ring-closure digraph / parity path cannot silently invert these labels.
 * Assertions compare the sorted multiset of emitted descriptors against the
 * reference, which is robust to atom-input ordering. Expected values are the
 * reference ground truth (IUPAC 2013 P-91), independent of engine output.
 *
 * Scope note: pseudo-asymmetric ring centres (lowercase r/s) and the
 * abstention behaviour on symmetric ring polyols are covered separately by
 * tests/test_v2_4_12_cyclitol_cip.js; this suite covers the uppercase-R/S
 * cases the engine resolves.
 *
 * Pure Node (real engine via shim.loadAll), no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var SmilesParser = global.SmilesParser;
var CIPStereo = global.CIPStereo;

var runner = shim.makeRunner('CIP ring stereocentre correctness');
var test = runner.test;
console.log('CIP ring stereocentre correctness');

// Sorted multiset of uppercase R/S descriptors emitted for a SMILES.
function descriptors(smi) {
    var m = SmilesParser.parse(smi);
    assert.ok(m && m.atoms.length && !(m.parseErrors && m.parseErrors.length), 'parsed ' + smi);
    CIPStereo.assign(m);
    var out = [];
    for (var i = 0; i < m.atoms.length; i++) {
        var a = m.atoms[i];
        if (a.chirality && a.cipLabel) out.push(a.cipLabel);
    }
    return out.sort().join('');
}

// Reference ground truth (IUPAC 2013; validated against an independent CIP
// labeller + PubChem). [label, SMILES, expected sorted R/S multiset].
var CASES = [
    ['trans-(1R,2R)-cyclohexane-1,2-diol', 'O[C@@H]1CCCC[C@H]1O', 'RR'],
    ['trans-(1S,2S)-cyclohexane-1,2-diol', 'O[C@H]1CCCC[C@@H]1O', 'SS'],
    ['cis-cyclohexane-1,2-diol (meso)',    'O[C@@H]1CCCC[C@@H]1O', 'RS'],
    ['(-)-menthol (1R,2S,5R)',             'C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C', 'RRS'],
    ['cis-2-methylcyclohexan-1-ol',        'C[C@H]1CCCC[C@H]1O', 'RS'],
    ['trans-2-methylcyclohexan-1-ol',      'C[C@@H]1CCCC[C@H]1O', 'RR'],
    ['alpha-D-glucopyranose',              'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', 'RRSSS'],
    ['(R)-3-methylcyclohexan-1-one',       'C[C@@H]1CCCC(=O)C1', 'R'],
    ['(S)-3-methylcyclohexan-1-one',       'C[C@H]1CCCC(=O)C1', 'S'],
    ['trans-2-decalol (fused ring)',       'O[C@H]1CC[C@@H]2CCCC[C@@H]2C1', 'RSS']
];

CASES.forEach(function (c, idx) {
    test('R' + (idx + 1) + ' ' + c[0] + ' -> ' + c[2], function () {
        var got = descriptors(c[1]);
        var exp = c[2].split('').sort().join('');
        assert.strictEqual(got, exp,
            c[0] + ': engine emitted [' + got + '] but the IUPAC reference is [' + exp + ']');
    });
});

// Adjacency guard: the original (false) report claimed adjacent 1,2-ring
// centres were inverted. Pin that BOTH centres of the trans-1,2-diol resolve
// (no abstention) and neither is dropped.
test('R11 adjacent 1,2-ring centres both resolve (no spurious abstention)', function () {
    var m = SmilesParser.parse('O[C@@H]1CCCC[C@H]1O');
    CIPStereo.assign(m);
    var labelled = 0, undet = 0;
    for (var i = 0; i < m.atoms.length; i++) {
        if (!m.atoms[i].chirality) continue;
        if (m.atoms[i].cipLabel) labelled++;
        else if (m.atoms[i].cipUndetermined) undet++;
    }
    assert.strictEqual(labelled, 2, 'both ring stereocentres should be labelled');
    assert.strictEqual(undet, 0, 'no ring stereocentre should abstain here');
});

// Acyclic anchors must stay correct (calibration sentinels) — proves the
// sign map is sound and shared by ring and acyclic centres alike.
test('R12 acyclic anchors unchanged (L-alanine S, D-glyceraldehyde R)', function () {
    assert.strictEqual(descriptors('N[C@@H](C)C(=O)O'), 'S', 'L-alanine = S');
    assert.strictEqual(descriptors('OC[C@@H](O)C=O'), 'R', 'D-glyceraldehyde = R');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
