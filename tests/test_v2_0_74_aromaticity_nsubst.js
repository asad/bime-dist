/**
 * tests/test_v2_0_74_aromaticity_nsubst.js — N-substituted pyrrole-type
 * aromatic nitrogen perception (v2.0.74 fix).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * A neutral aromatic nitrogen donates its lone pair to the pi system
 * (pyrrole-type, 2 e-) when it is sigma-saturated with 3 bonds — whether
 * that third bond is an H ([nH]) OR a C/N substituent (N-alkyl / N-aryl).
 * Before v2.0.74, Molecule.perceiveAromaticity credited the 2-electron
 * pyrrole contribution only when an explicit H was present, so every
 * N-substituted azole (N-methylpyrrole, N-methylimidazole, caffeine,
 * theobromine, 1-substituted indoles / benzimidazoles) failed Hueckel and
 * perceived as NON-aromatic. That corrupted rendering (no aromatic circle),
 * aromatic SMARTS matching, the AAM aromatic-ring-conservation metric, and
 * the canonical SMILES writer (it de-aromatised the ring on output).
 *
 * Pyridine-type N (degree-2, in-plane lone pair) stays 1 e-, and a
 * positively-charged N (pyridinium [nH+], lone pair protonated) stays
 * pyridine-type — both verified below as regression guards.
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Aromaticity N-substituted (v2.0.74)');
var test = runner.test;

console.log('Aromaticity N-substituted (v2.0.74)');

function perceiveCount(smi) {
    var m = SmilesParser.parse(smi);
    if (m.parseErrors.length) { throw new Error('parse: ' + m.parseErrors.join('; ')); }
    var per = Molecule.perceiveAromaticity(m);
    var c = 0;
    for (var k in per) { if (per.hasOwnProperty(k) && per[k]) { c++; } }
    return c;
}

// ---------------------------------------------------------------------------
// A. The fix — N-substituted pyrrole-type N now donates its lone pair.
// ---------------------------------------------------------------------------
test('A1 N-methylpyrrole (aromatic input) → 5 aromatic atoms', function () {
    assert.strictEqual(perceiveCount('Cn1cccc1'), 5);
});

test('A2 N-methylimidazole → 5 aromatic atoms', function () {
    assert.strictEqual(perceiveCount('Cn1ccnc1'), 5);
});

test('A3 1-methylindole → 9 aromatic atoms (both fused rings)', function () {
    assert.strictEqual(perceiveCount('Cn1ccc2ccccc12'), 9);
});

test('A4 Kekule N-methylpyrrole CN1C=CC=C1 → 5 aromatic atoms', function () {
    // Exercises the non-aromatic (Kekule) branch of perceiveAromaticity.
    assert.strictEqual(perceiveCount('CN1C=CC=C1'), 5);
});

test('A5 caffeine imidazole ring is perceived aromatic (>= 5 atoms)', function () {
    // Caffeine: imidazole ring aromatic; the pyrimidinedione ring is NOT
    // (cyclic diamide). At least the 5-membered imidazole must be aromatic.
    assert.ok(perceiveCount('Cn1cnc2c1c(=O)n(C)c(=O)n2C') >= 5,
        'caffeine must perceive its imidazole ring as aromatic');
});

// ---------------------------------------------------------------------------
// B. Regression guards — pyridine-type N and neutral/charged cases unchanged.
// ---------------------------------------------------------------------------
test('B1 benzene → 6', function () { assert.strictEqual(perceiveCount('c1ccccc1'), 6); });
test('B2 pyridine → 6 (degree-2 N stays pyridine-type)', function () {
    assert.strictEqual(perceiveCount('c1ccncc1'), 6);
});
test('B3 pyrrole [nH] → 5 (unchanged)', function () {
    assert.strictEqual(perceiveCount('c1cc[nH]c1'), 5);
});
test('B4 pyridinium [nH+] → 6 (protonated N stays pyridine-type, ring aromatic)',
    function () {
        assert.strictEqual(perceiveCount('c1cc[nH+]cc1'), 6);
    });
test('B5 naphthalene → 10 (unchanged)', function () {
    assert.strictEqual(perceiveCount('c1ccc2ccccc2c1'), 10);
});
test('B6 piperidine → 0 (saturated N-ring NOT aromatised)', function () {
    assert.strictEqual(perceiveCount('C1CCNCC1'), 0);
});
test('B7 N-methylpiperidine → 0 (saturated, sp3 C guard holds)', function () {
    assert.strictEqual(perceiveCount('CN1CCCCC1'), 0);
});

// ---------------------------------------------------------------------------
// C. Writer cascade — the canonical writer now emits aromatic lowercase for
// the newly-aromatic ring instead of de-aromatising it.
// ---------------------------------------------------------------------------
test('C1 writer emits aromatic lowercase for N-methylimidazole', function () {
    var out = SmilesWriter.write(SmilesParser.parse('Cn1ccnc1'));
    // Must contain a lowercase aromatic ring atom (n or c), not a fully
    // upper-cased Kekule de-aromatisation.
    assert.ok(/[a-z]/.test(out.replace(/Cl|Br/g, '')),
        'expected aromatic lowercase in writer output, got: ' + out);
    var rt = SmilesParser.parse(out);
    assert.strictEqual(rt.parseErrors.length, 0, 'writer output must re-parse');
    assert.strictEqual(rt.atoms.length, 6, 'atom count preserved through write');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
