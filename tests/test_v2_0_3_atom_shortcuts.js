/**
 * tests/test_v2_0_3_atom_shortcuts.js — v2.0.3 multi-letter atom shortcuts.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Exercises the SmilesParser abbreviation table introduced in v2.0.3.
 * Two flavours share one dispatch:
 *
 *   Organic substructure aliases  Ph, Me, Et, iPr, tBu, Bn, Bz
 *     -> expand inline into full atom/bond graph; Ph IS a benzene ring.
 *   Amino-acid residue placeholders  Gly, Ala, Ser, Val, Leu
 *     -> single labelled atom (R-group form, same family as [R]/[R1]).
 *        GlyAla parses as two residue nodes joined by one single bond;
 *        the inner peptide-bond atoms are not drawn — the label IS the
 *        residue. SmilesWriter emits bracket form [Gly] on round-trip
 *        and the bracket-atom parser accepts 3-char element symbols so
 *        every residue round-trips cleanly.
 *
 * Ring digits inside an organic-alias expansion live in their own
 * namespace (the expansion functions build atoms+bonds directly rather
 * than re-entering the parser), so the aromatic ring inside Ph cannot
 * collide with any digit the outer SMILES is using.
 *
 * Regression: every shorthand must NOT shadow an existing SMILES token.
 * The dispatch is case-sensitive and longest-first, so Br / Cl / Pt
 * still parse as today.
 */
'use strict';

var assert = require('assert');
var path   = require('path');
var shim   = require(path.join(__dirname, 'shim.js'));
shim.loadAll();

var runner = shim.makeRunner('Atom shortcuts v2.0.3');
var test   = runner.test;

console.log('Atom shortcuts v2.0.3');

// Composition helper: returns a sorted "C:n N:n O:n C(a):n" descriptor of
// the atom multiset. Bond-count check is separate (asserted directly).
function comp(smi) {
    var m = SmilesParser.parse(smi);
    var c = {};
    for (var i = 0; i < m.atoms.length; i++) {
        var a = m.atoms[i];
        var k = a.symbol + (a.aromatic ? '(a)' : '');
        c[k] = (c[k] || 0) + 1;
    }
    var keys = Object.keys(c).sort();
    return {
        mol: m,
        signature: keys.map(function (k) { return k + ':' + c[k]; }).join(' '),
        errors: m.parseErrors.slice()
    };
}

// =========================================================================
// (A) Organic aliases
// =========================================================================

test('A1. Ph expands to a 6-atom aromatic ring (C6H5-)', function () {
    var r = comp('Ph');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C(a):6');
    assert.strictEqual(r.mol.bonds.length, 6);
});

test('A2. Me expands to a single carbon', function () {
    var r = comp('Me');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:1');
    assert.strictEqual(r.mol.bonds.length, 0);
});

test('A3. Et expands to C-C (two-carbon chain)', function () {
    var r = comp('Et');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:2');
    assert.strictEqual(r.mol.bonds.length, 1);
});

test('A4. iPr expands to a Y-shaped C-C-C (central CH with two CH3)',
function () {
    var r = comp('iPr');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:3');
    assert.strictEqual(r.mol.bonds.length, 2);
    // Atom 0 is the central CH; degree must be 2.
    assert.strictEqual(r.mol.getNeighbors(r.mol.atoms[0].id).length, 2);
});

test('A5. tBu expands to a 4-carbon star (central C with three CH3)',
function () {
    var r = comp('tBu');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:4');
    assert.strictEqual(r.mol.bonds.length, 3);
    assert.strictEqual(r.mol.getNeighbors(r.mol.atoms[0].id).length, 3);
});

test('A6. Bn expands to CH2-phenyl (methylene + 6-aromatic ring)',
function () {
    var r = comp('Bn');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:1 C(a):6');
    assert.strictEqual(r.mol.bonds.length, 7);   // 1 CH2-Caryl + 6 ring
});

test('A7. Bz expands to C(=O)-phenyl (benzoyl group)', function () {
    var r = comp('Bz');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:1 C(a):6 O:1');
    assert.strictEqual(r.mol.bonds.length, 8);   // 1 C=O + 1 C-Caryl + 6 ring
    var doubles = r.mol.bonds.filter(function (b) { return b.type === 2; });
    assert.strictEqual(doubles.length, 1, 'exactly one double bond (the C=O)');
});

// =========================================================================
// (B) Amino-acid residues — R-group / labelled-atom form
//
// Each residue collapses to ONE atom whose symbol is the three-letter
// code (Gly / Ala / Ser / Val / Leu) and whose hydrogens are pinned at
// 0 (it's a residue label, not an element with a valence). Round-trip
// through SmilesWriter produces the bracket form ([Gly] etc.) which
// the bracket-atom parser accepts because the codes are listed in
// ALL_ELEMENTS alongside generic [R] / [R1].
// =========================================================================

test('B1. Gly is a single labelled atom (R-group form)', function () {
    var r = comp('Gly');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.mol.atoms.length, 1);
    assert.strictEqual(r.mol.atoms[0].symbol, 'Gly');
    assert.strictEqual(r.mol.atoms[0].hydrogens, 0);
    assert.strictEqual(r.mol.bonds.length, 0);
});

test('B2. Ala / Ser / Val / Leu each collapse to one labelled atom',
function () {
    ['Ala', 'Ser', 'Val', 'Leu'].forEach(function (code) {
        var r = comp(code);
        assert.deepStrictEqual(r.errors, [], code + ' should parse');
        assert.strictEqual(r.mol.atoms.length, 1, code + ' is one atom');
        assert.strictEqual(r.mol.atoms[0].symbol, code,
            code + ' atom symbol matches the code');
        assert.strictEqual(r.mol.atoms[0].hydrogens, 0,
            code + ' has hydrogens pinned at 0');
        assert.strictEqual(r.mol.bonds.length, 0, code + ' has no bonds');
    });
});

test('B3. GlyAla is two residue nodes joined by a single bond',
function () {
    var r = comp('GlyAla');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.mol.atoms.length, 2);
    assert.strictEqual(r.mol.atoms[0].symbol, 'Gly');
    assert.strictEqual(r.mol.atoms[1].symbol, 'Ala');
    assert.strictEqual(r.mol.bonds.length, 1);
    assert.strictEqual(r.mol.bonds[0].type, 1);
    // The bond connects the two residues directly.
    var nbrs0 = r.mol.getNeighbors(r.mol.atoms[0].id);
    assert.ok(nbrs0.indexOf(r.mol.atoms[1].id) >= 0,
        'Gly node must bond to Ala node');
});

test('B4. GlyAlaSer is a 3-node residue chain', function () {
    var r = comp('GlyAlaSer');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.mol.atoms.length, 3);
    assert.strictEqual(r.mol.atoms.map(function (a) { return a.symbol; }).join('-'),
        'Gly-Ala-Ser');
    assert.strictEqual(r.mol.bonds.length, 2);
    // Chain topology: Gly-Ala and Ala-Ser, NOT Gly-Ser.
    assert.strictEqual(r.mol.getNeighbors(r.mol.atoms[1].id).length, 2,
        'middle Ala node has two neighbours');
});

test('B5. Residue placeholders round-trip through SmilesWriter', function () {
    ['Gly', 'Ala', 'Ser', 'Val', 'Leu', 'GlyAla', 'GlyAlaSer'].forEach(
    function (smi) {
        var m1 = SmilesParser.parse(smi);
        assert.strictEqual(m1.parseErrors.length, 0, smi + ' parses');
        var written = SmilesWriter.write(m1);
        var m2 = SmilesParser.parse(written);
        assert.strictEqual(m2.parseErrors.length, 0,
            smi + ' round-trips: written form ' + written + ' must re-parse');
        assert.strictEqual(m2.atoms.length, m1.atoms.length,
            smi + ' round-trip preserves atom count');
        assert.strictEqual(m2.bonds.length, m1.bonds.length,
            smi + ' round-trip preserves bond count');
    });
});

test('B6. [Gly] bracket form parses to the same residue atom as bare Gly',
function () {
    var bare    = SmilesParser.parse('Gly');
    var bracket = SmilesParser.parse('[Gly]');
    assert.strictEqual(bare.parseErrors.length, 0);
    assert.strictEqual(bracket.parseErrors.length, 0);
    assert.strictEqual(bare.atoms[0].symbol, 'Gly');
    assert.strictEqual(bracket.atoms[0].symbol, 'Gly');
});

test('B7. Residues attach cleanly to ordinary chemistry (CGly, PhGly)',
function () {
    var cGly = comp('CGly');
    assert.deepStrictEqual(cGly.errors, []);
    assert.strictEqual(cGly.signature, 'C:1 Gly:1');
    assert.strictEqual(cGly.mol.bonds.length, 1);

    var phGly = comp('PhGly');
    assert.deepStrictEqual(phGly.errors, []);
    // 6 aromatic C from Ph + 1 Gly node = 7 atoms; 6 ring bonds + 1 attach.
    assert.strictEqual(phGly.signature, 'C(a):6 Gly:1');
    assert.strictEqual(phGly.mol.bonds.length, 7);
});

// =========================================================================
// (C) Integration with bonds, branches, and ring digits
// =========================================================================

test('C1. Phenyl as a branch (CC(Ph)C) bonds the parent C to ring atom 0',
function () {
    var r = comp('CC(Ph)C');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:3 C(a):6');   // 3 aliphatic + 6 aromatic
    assert.strictEqual(r.mol.bonds.length, 9);       // 2 chain + 6 ring + 1 Ph-attach
});

test('C2. PhCC continues from Ph (treating Ph as one node)', function () {
    var r = comp('PhCC');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:2 C(a):6');
    assert.strictEqual(r.mol.bonds.length, 8);   // 6 ring + 1 Ph-C + 1 C-C
});

test('C3. Ring digits inside an abbreviation do NOT collide with outer rings',
function () {
    // Ph's internal aromatic ring is built directly (no parser ring digit
    // consumed). Outer ring 1 still closes correctly.
    var r = comp('C1CCCCC1Ph');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:6 C(a):6');
    // 6 aliphatic ring bonds + 6 aromatic ring bonds + 1 C-Ph linker = 13.
    assert.strictEqual(r.mol.bonds.length, 13);
});

test('C4. Screenshot SMILES: \'C(O)(Ph)Ph\' now parses with two phenyls',
function () {
    // The screenshot used [N+]3(CCCCc1ccccc1)(CCOc2ccccc2)CCC(CC3)OC(=O)C(O)(Ph)Ph
    // (with one user typo C0 -> CCO corrected here). Only the abbreviation
    // expansion is on test; the layout is exercised by the diagnostic.
    var smi = '[N+]3(CCCCc1ccccc1)(CCOc2ccccc2)CCC(CC3)OC(=O)C(O)(Ph)Ph';
    var r = comp(smi);
    assert.deepStrictEqual(r.errors, []);
    // Three benzene rings (c1, c2, Ph, Ph) = 24 aromatic C; 13 aliphatic C;
    // 1 N; 4 O. 4 rings each contribute 6 ring bonds + 1 attachment each
    // = 24 + 4 = 28 ring-system bonds. Plus the rest of the backbone.
    var hasN = r.mol.atoms.filter(function (a) { return a.symbol === 'N'; });
    assert.strictEqual(hasN.length, 1);
    var aromaticC = r.mol.atoms.filter(function (a) {
        return a.symbol === 'C' && a.aromatic;
    });
    assert.strictEqual(aromaticC.length, 24,
        'four benzene rings = 24 aromatic carbons');
});

// =========================================================================
// (D) Regression — abbreviations MUST NOT shadow existing SMILES tokens
// =========================================================================

test('D1. Br (bromine, organic two-char) is unaffected', function () {
    var r = comp('BrCCBr');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'Br:2 C:2');
});

test('D2. Cl (chlorine) is unaffected', function () {
    var r = comp('ClCCl');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C:1 Cl:2');
});

test('D3. Bracket atoms (e.g. [Pt]) are unaffected — abbreviation lookup '
+ 'is BEFORE the case-based dispatch but Pt is bracketed', function () {
    var r = comp('[Pt](Cl)(Cl)(N)N');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'Cl:2 N:2 Pt:1');
});

test('D4. Single-letter organic atoms still work (P phosphorus etc.)',
function () {
    var r = comp('OP(=O)(O)O');           // phosphate
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'O:4 P:1');
});

test('D5. Aromatic ring c1ccccc1 (lowercase atom literals) still works',
function () {
    var r = comp('c1ccccc1');
    assert.deepStrictEqual(r.errors, []);
    assert.strictEqual(r.signature, 'C(a):6');
    assert.strictEqual(r.mol.bonds.length, 6);
});

test('D6. Alanine canonical SMILES N[C@@H](C)C(=O)O still parses with chirality',
function () {
    var r = comp('N[C@@H](C)C(=O)O');
    assert.deepStrictEqual(r.errors, []);
    var chiralAtoms = r.mol.atoms.filter(function (a) { return a.chirality; });
    assert.strictEqual(chiralAtoms.length, 1, 'alpha-C is chiral');
});

test('D7. Morphine round-trip (v1.8.34 canary) is unaffected — no atom '
+ 'in the molecule starts with an abbreviation prefix', function () {
    var r = comp('CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5');
    assert.deepStrictEqual(r.errors, []);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
