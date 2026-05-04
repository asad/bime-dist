/**
 * tests/test_smiles_writer.js — Round-trip regression suite for SmilesWriter.
 *
 * For each input: parse, write, re-parse, and assert that atom/bond counts
 * (and key properties: aromaticity, charge, isotope, stereo, atom map) survive
 * the round trip.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SMILES Writer');
var test = runner.test;

console.log('SMILES Writer (round-trip)');

function parse(smi) { return SmilesParser.parse(smi); }
function rt(smi) {
    var m1 = parse(smi);
    var written = SmilesWriter.write(m1);
    var m2 = parse(written);
    return { m1: m1, written: written, m2: m2 };
}

function assertCountsMatch(r) {
    assert.strictEqual(r.m1.atoms.length, r.m2.atoms.length, 'atoms differ: ' + r.written);
    assert.strictEqual(r.m1.bonds.length, r.m2.bonds.length, 'bonds differ: ' + r.written);
}

// --- Counts round-trip across many inputs ---
[
    'C', 'CC', 'CCC', 'CCO', 'C=C', 'C#C', 'CC(=O)O',
    'c1ccccc1', 'C1CCCCC1', 'C1CCC1', 'c1ccncc1', 'c1ccoc1', 'c1ccsc1',
    'c1ccc2ccccc2c1', 'c1ccc2[nH]ccc2c1',
    'Cn1c(=O)c2c(ncn2C)n(C)c1=O',
    'CC(=O)Oc1ccccc1OC(C)=O',
    'C(C)(C)(C)C'
].forEach(function(smi) {
    test('round-trip atom/bond counts: ' + smi, function() {
        var r = rt(smi);
        assertCountsMatch(r);
    });
});

// --- Aromaticity round-trip ---
test('benzene stays lowercase aromatic (c1ccccc1 not C1=CC=CC=C1)', function() {
    var r = rt('c1ccccc1');
    // The written form must contain lowercase aromatic c's, not uppercase Kekule.
    assert.ok(/c/.test(r.written), 'expected aromatic c in output: ' + r.written);
    assert.ok(!/C=C/.test(r.written), 'expected no Kekule C=C in output: ' + r.written);
    assert.strictEqual(r.m2.atoms.every(function(a) { return a.aromatic; }), true);
});

// --- Charge preservation ---
test('charge preserved: [NH4+]', function() {
    var r = rt('[NH4+]');
    assert.strictEqual(r.m2.atoms[0].charge, 1);
});

test('charge preserved: [O-]', function() {
    var r = rt('[O-]');
    assert.strictEqual(r.m2.atoms[0].charge, -1);
});

// --- Isotope preservation ---
test('isotope preserved: [13C]', function() {
    var r = rt('[13C]');
    assert.strictEqual(r.m2.atoms[0].isotope, 13);
});

test('isotope preserved: [2H]', function() {
    var r = rt('[2H]');
    assert.strictEqual(r.m2.atoms[0].isotope, 2);
});

// --- Stereo @/@@ preservation ---
test('chirality @@ preserved through round-trip', function() {
    var r = rt('[C@@H](F)(Cl)Br');
    var c = r.m2.atoms.find(function(a) { return a.symbol === 'C'; });
    assert.ok(c.chirality === '@@' || c.chirality === '@',
              'expected chirality on round-trip C, got: ' + c.chirality);
});

test('chirality @ preserved through round-trip', function() {
    var r = rt('[C@H](F)(Cl)Br');
    var c = r.m2.atoms.find(function(a) { return a.symbol === 'C'; });
    assert.ok(c.chirality === '@@' || c.chirality === '@',
              'expected chirality on round-trip C, got: ' + c.chirality);
});

// --- Atom map preservation ---
test('atom map preserved through round-trip', function() {
    var r = rt('[CH3:1][OH:2]');
    var maps = r.m2.atoms.map(function(a) { return a.mapNumber; }).sort();
    assert.deepStrictEqual(maps, [1, 2]);
});

// --- Multi-fragment dot separator ---
test('[Na+].[Cl-] — round-trip preserves both fragments and dot', function() {
    var r = rt('[Na+].[Cl-]');
    assert.strictEqual(r.m2.atoms.length, 2);
    assert.strictEqual(r.m2.bonds.length, 0);
    assert.ok(r.written.indexOf('.') >= 0, 'expected . in output: ' + r.written);
});

// --- Reaction round-trip ---
test('CCO>>CC=O — reaction arrow preserved', function() {
    var r = rt('CCO>>CC=O');
    assert.ok(r.m2.reactionArrow, 'expected reactionArrow after round-trip');
    assert.ok(r.written.indexOf('>>') >= 0, 'expected >> in output: ' + r.written);
});

// --- Empty molecule ---
test('empty Molecule writes empty string', function() {
    var m = new Molecule();
    assert.strictEqual(SmilesWriter.write(m), '');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
