/**
 * tests/test_smiles_parser.js — Regression suite for SmilesParser.
 *
 * Locks in the behaviours catalogued in TEST-REPORT.md §1 (BIME v1.0.1+).
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SMILES Parser');
var test = runner.test;

console.log('SMILES Parser');

function parse(smi) {
    return SmilesParser.parse(smi);
}

// --- Basic atoms ---
test('parses single C as one atom, no bonds', function() {
    var m = parse('C');
    assert.strictEqual(m.atoms.length, 1);
    assert.strictEqual(m.bonds.length, 0);
    assert.strictEqual(m.atoms[0].symbol, 'C');
});

test('parses single N as one atom', function() {
    var m = parse('N');
    assert.strictEqual(m.atoms.length, 1);
    assert.strictEqual(m.atoms[0].symbol, 'N');
});

test('parses single O as one atom', function() {
    var m = parse('O');
    assert.strictEqual(m.atoms.length, 1);
    assert.strictEqual(m.atoms[0].symbol, 'O');
});

// --- Bonds ---
test('parses C=C as 2 atoms, double bond', function() {
    var m = parse('C=C');
    assert.strictEqual(m.atoms.length, 2);
    assert.strictEqual(m.bonds.length, 1);
    assert.strictEqual(m.bonds[0].type, Molecule.BOND_DOUBLE);
});

test('parses C#C as 2 atoms, triple bond', function() {
    var m = parse('C#C');
    assert.strictEqual(m.atoms.length, 2);
    assert.strictEqual(m.bonds.length, 1);
    assert.strictEqual(m.bonds[0].type, Molecule.BOND_TRIPLE);
});

// --- Rings ---
test('parses benzene c1ccccc1 as 6 aromatic atoms, 6 bonds', function() {
    var m = parse('c1ccccc1');
    assert.strictEqual(m.atoms.length, 6);
    assert.strictEqual(m.bonds.length, 6);
    assert.strictEqual(m.atoms.every(function(a) { return a.aromatic; }), true);
});

test('parses cyclohexane C1CCCCC1 as 6 atoms, 6 bonds', function() {
    var m = parse('C1CCCCC1');
    assert.strictEqual(m.atoms.length, 6);
    assert.strictEqual(m.bonds.length, 6);
    assert.strictEqual(m.atoms.some(function(a) { return a.aromatic; }), false);
});

test('parses cyclobutane C1CCC1 as 4 atoms, 4 bonds', function() {
    var m = parse('C1CCC1');
    assert.strictEqual(m.atoms.length, 4);
    assert.strictEqual(m.bonds.length, 4);
});

// --- Heterocycles ---
test('parses pyridine c1ccncc1 as 5C + 1N, all aromatic', function() {
    var m = parse('c1ccncc1');
    assert.strictEqual(m.atoms.length, 6);
    assert.strictEqual(m.atoms.filter(function(a) { return a.symbol === 'N'; }).length, 1);
    assert.strictEqual(m.atoms.every(function(a) { return a.aromatic; }), true);
});

test('parses pyrrole c1cc[nH]c1 — 5 atoms, [nH] H=1', function() {
    var m = parse('c1cc[nH]c1');
    assert.strictEqual(m.atoms.length, 5);
    var n = m.atoms.find(function(a) { return a.symbol === 'N'; });
    assert.strictEqual(n.hydrogens, 1);
    assert.strictEqual(n.aromatic, true);
});

test('parses furan c1ccoc1 — aromatic O', function() {
    var m = parse('c1ccoc1');
    assert.strictEqual(m.atoms.length, 5);
    var o = m.atoms.find(function(a) { return a.symbol === 'O'; });
    assert.strictEqual(o.aromatic, true);
});

test('parses thiophene c1ccsc1 — aromatic S', function() {
    var m = parse('c1ccsc1');
    assert.strictEqual(m.atoms.length, 5);
    var s = m.atoms.find(function(a) { return a.symbol === 'S'; });
    assert.strictEqual(s.aromatic, true);
});

// --- Fused rings ---
test('parses naphthalene as 10 aromatic atoms, 11 bonds', function() {
    var m = parse('c1ccc2ccccc2c1');
    assert.strictEqual(m.atoms.length, 10);
    assert.strictEqual(m.bonds.length, 11);
    assert.strictEqual(m.atoms.filter(function(a) { return a.aromatic; }).length, 10);
});

test('parses indole as 9 atoms', function() {
    var m = parse('c1ccc2[nH]ccc2c1');
    assert.strictEqual(m.atoms.length, 9);
});

test('parses caffeine — 14 atoms, 15 bonds', function() {
    var m = parse('Cn1c(=O)c2c(ncn2C)n(C)c1=O');
    assert.strictEqual(m.atoms.length, 14);
    assert.strictEqual(m.bonds.length, 15);
    assert.strictEqual(m.parseErrors.length, 0);
});

// --- Charged ---
test('parses [NH4+] — charge=+1, hydrogens=4', function() {
    var m = parse('[NH4+]');
    assert.strictEqual(m.atoms.length, 1);
    assert.strictEqual(m.atoms[0].charge, 1);
    assert.strictEqual(m.atoms[0].hydrogens, 4);
});

test('parses [O-] — charge=-1', function() {
    var m = parse('[O-]');
    assert.strictEqual(m.atoms[0].charge, -1);
});

test('parses [Na+].[Cl-] — 2 atoms, 0 bonds, dot fragment', function() {
    var m = parse('[Na+].[Cl-]');
    assert.strictEqual(m.atoms.length, 2);
    assert.strictEqual(m.bonds.length, 0);
});

// --- Isotopes ---
test('parses [2H] as deuterium with isotope=2', function() {
    var m = parse('[2H]');
    assert.strictEqual(m.atoms[0].symbol, 'H');
    assert.strictEqual(m.atoms[0].isotope, 2);
});

test('parses [13C] with isotope=13', function() {
    var m = parse('[13C]');
    assert.strictEqual(m.atoms[0].symbol, 'C');
    assert.strictEqual(m.atoms[0].isotope, 13);
});

// --- Stereo ---
test('parses [C@@H](F)(Cl)Br — chirality stored on C, H=1', function() {
    var m = parse('[C@@H](F)(Cl)Br');
    assert.strictEqual(m.atoms.length, 4);
    var c = m.atoms.find(function(a) { return a.symbol === 'C'; });
    assert.strictEqual(c.chirality, '@@');
    assert.strictEqual(c.hydrogens, 1);
});

test('parses /C=C/ — 2 atoms, 1 double bond', function() {
    var m = parse('/C=C/');
    assert.strictEqual(m.atoms.length, 2);
    assert.strictEqual(m.bonds.length, 1);
    assert.strictEqual(m.bonds[0].type, Molecule.BOND_DOUBLE);
});

// --- Reactions ---
test('parses CCO>>CC=O — reactant + product + arrow', function() {
    var m = parse('CCO>>CC=O');
    assert.strictEqual(m.atoms.length, 6);
    assert.ok(m.reactionArrow);
});

test('parses [CH3:1][OH:2]>>[CH2:1]=[O:2] — atom maps preserved', function() {
    var m = parse('[CH3:1][OH:2]>>[CH2:1]=[O:2]');
    assert.strictEqual(m.atoms.length, 4);
    var maps = m.atoms.map(function(a) { return a.mapNumber; }).sort();
    assert.deepStrictEqual(maps, [1, 1, 2, 2]);
});

// --- Edge cases ---
test('parses empty string — 0 atoms, no errors', function() {
    var m = parse('');
    assert.strictEqual(m.atoms.length, 0);
});

test('parses 50-atom chain (CCC...C) — 50 atoms, 49 bonds', function() {
    var m = parse('C'.repeat(50));
    assert.strictEqual(m.atoms.length, 50);
    assert.strictEqual(m.bonds.length, 49);
    assert.strictEqual(m.parseErrors.length, 0);
});

test('parses unclosed ring C1CCC — captures error', function() {
    var m = parse('C1CCC');
    assert.ok(m.parseErrors.length > 0);
    assert.ok(m.parseErrors.some(function(e) { return /[Uu]nclosed.*ring/.test(e); }));
});

test('parses [unclosed — captures bracket error', function() {
    var m = parse('[unclosed');
    assert.ok(m.parseErrors.length > 0);
});

// --- Atom map preservation through parse ---
test('atom map :7 preserved on bracket atom', function() {
    var m = parse('[CH4:7]');
    assert.strictEqual(m.atoms[0].mapNumber, 7);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
