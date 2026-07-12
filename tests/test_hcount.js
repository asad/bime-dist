/**
 * tests/test_hcount.js — Implicit hydrogen calculation regressions for
 * Molecule.calcHydrogens().
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('H Count');
var test = runner.test;

console.log('H Count');

function P(smi) { return SmilesParser.parse(smi); }

function calcFor(smi, predicate) {
    var m = P(smi);
    var atom = m.atoms.find(predicate);
    assert.ok(atom, 'no atom matched predicate in ' + smi);
    return m.calcHydrogens(atom.id);
}
function calcFirst(smi) { return calcFor(smi, function() { return true; }); }

test('C → 4 H', function() {
    assert.strictEqual(calcFirst('C'), 4);
});

test('O → 2 H', function() {
    assert.strictEqual(calcFirst('O'), 2);
});

test('N → 3 H', function() {
    assert.strictEqual(calcFirst('N'), 3);
});

test('[NH4+] → 4 H (explicit, calcHydrogens returns explicit)', function() {
    assert.strictEqual(calcFirst('[NH4+]'), 4);
});

// Document the current policy: parser stores hydrogens=0 for [O-]
// (no H token in bracket → hydrogens stays at 0). calcHydrogens returns
// the stored count when ≥0 — so [O-] reports 0 H. This pins that behavior.
test('[O-] → 0 H (BIME stores hydrogens=0 for bracket atoms without H token)', function() {
    assert.strictEqual(calcFirst('[O-]'), 0);
});

test('aromatic c (in benzene) → 1 H', function() {
    var m = P('c1ccccc1');
    assert.strictEqual(m.calcHydrogens(m.atoms[0].id), 1);
});

test('[nH] → 1 H (explicit)', function() {
    assert.strictEqual(calcFor('c1cc[nH]c1', function(a) { return a.symbol === 'N'; }), 1);
});

test('pyridine n → 0 H (lone-pair donor, no NH)', function() {
    assert.strictEqual(calcFor('c1ccncc1', function(a) { return a.symbol === 'N'; }), 0);
});

// Documented behavior: methylsulfonyl-H S(=O)(=O)C has bondOrderSum=5 on S,
// adjBondSum=5, smallest valence ≥5 is 6 → 6-5 = 1 H. Locks this in.
test('S(=O)(=O)C → S has 1 H (hypervalent sulfur, valence-table fallback)', function() {
    assert.strictEqual(calcFor('S(=O)(=O)C', function(a) { return a.symbol === 'S'; }), 1);
});

test('round-trip methane: parse(write(parse("C"))) gives calcH=4', function() {
    var m1 = P('C');
    var s = SmilesWriter.write(m1);
    var m2 = P(s);
    assert.strictEqual(m2.calcHydrogens(m2.atoms[0].id), 4);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
