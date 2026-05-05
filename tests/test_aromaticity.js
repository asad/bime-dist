/**
 * tests/test_aromaticity.js — Hückel aromaticity regressions through the
 * SmilesParser → Molecule pipeline. Locks in BIME's current policy.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Aromaticity');
var test = runner.test;

console.log('Aromaticity');

function P(smi) { return SmilesParser.parse(smi); }
function aromaticCount(smi) {
    var m = P(smi);
    return m.atoms.filter(function(a) { return a.aromatic; }).length;
}
function totalAtoms(smi) { return P(smi).atoms.length; }

test('benzene → 6 aromatic atoms', function() {
    assert.strictEqual(aromaticCount('c1ccccc1'), 6);
});

test('pyridine → 6 aromatic atoms', function() {
    assert.strictEqual(aromaticCount('c1ccncc1'), 6);
});

test('pyrrole → 5 aromatic atoms (NH lone-pair donates)', function() {
    assert.strictEqual(aromaticCount('c1cc[nH]c1'), 5);
});

test('furan → 5 aromatic atoms', function() {
    assert.strictEqual(aromaticCount('c1ccoc1'), 5);
});

test('thiophene → 5 aromatic atoms', function() {
    assert.strictEqual(aromaticCount('c1ccsc1'), 5);
});

test('tetrazole c1nnn[nH]1 → all 5 aromatic', function() {
    assert.strictEqual(aromaticCount('c1nnn[nH]1'), 5);
});

test('imidazole c1nc[nH]c1 → all 5 aromatic', function() {
    assert.strictEqual(aromaticCount('c1nc[nH]c1'), 5);
});

test('cyclopentadiene C1=CCC=C1 → NOT aromatic (0 aromatic atoms)', function() {
    assert.strictEqual(aromaticCount('C1=CCC=C1'), 0);
});

test('cyclooctatetraene C1=CC=CC=CC=C1 → NOT aromatic (8 atoms, 0 aromatic)', function() {
    var m = P('C1=CC=CC=CC=C1');
    assert.strictEqual(m.atoms.length, 8);
    assert.strictEqual(aromaticCount('C1=CC=CC=CC=C1'), 0);
});

test('naphthalene → 10 atoms, all aromatic', function() {
    assert.strictEqual(totalAtoms('c1ccc2ccccc2c1'), 10);
    assert.strictEqual(aromaticCount('c1ccc2ccccc2c1'), 10);
});

test('indole → 9 atoms, all aromatic', function() {
    assert.strictEqual(totalAtoms('c1ccc2[nH]ccc2c1'), 9);
    assert.strictEqual(aromaticCount('c1ccc2[nH]ccc2c1'), 9);
});

test('cyclohexane C1CCCCC1 → NOT aromatic (0 aromatic)', function() {
    assert.strictEqual(aromaticCount('C1CCCCC1'), 0);
});

// --- BIME policy on 2-pyridone ---
// Document, don't fix. The TEST-REPORT notes BIME excludes pyridone aromaticity
// while RDKit includes it. This test pins the current behaviour:
//   - When parsed in lowercase aromatic form, parser preserves aromatic flags.
//   - When parsed in Kekule form (O=C1C=CC=CN1), no aromatic flags are set.
//   - When written, perceiveAromaticity drops the aromatic ring (writer outputs uppercase).
test('Kekule 2-pyridone (O=C1C=CC=CN1) → 0 aromatic atoms (BIME policy)', function() {
    assert.strictEqual(aromaticCount('O=C1C=CC=CN1'), 0);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
