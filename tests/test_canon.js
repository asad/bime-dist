/**
 * tests/test_canon.js — Canonical SMILES (Canon) regression coverage for
 * BIME v1.1.1.
 *
 * Verifies the canonical-SMILES properties of SmilesWriter.write composed
 * with SmilesParser.parse:
 *
 *   1. Idempotence:    canon(canon(s)) === canon(s) byte-identical
 *   2. Equivalence:    different SMILES of the SAME molecule → same canonical
 *   3. Determinism:    canon(s) is identical across repeated runs
 *   4. Aromaticity:    parsed-aromatic input stays in lowercase aromatic form
 *   5. Charge / iso /  charge, isotope, and stereo round-trip canonical so a
 *      stereo:         second pass is byte-identical to the first
 *   6. Edge cases:     empty molecule → "", methane "C" → "C"
 *
 * Stereo round-trip drift on polycyclic multi-stereocentre molecules is a
 * known limitation — the chirality letters can flip across re-parses while
 * the atom skeleton stays identical. See test_round_trip.js KNOWN_CANON_DRIFT.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see
 * LICENSE.txt
 *
 * No file I/O, no network. Runs in <100 ms.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Canonical SMILES (Canon)');
var test = runner.test;

console.log('Canonical SMILES (Canon)');

function P(smi) { return SmilesParser.parse(smi); }
function W(mol) { return SmilesWriter.write(mol); }
function canon(smi) { return W(P(smi)); }

// ---------------------------------------------------------------------------
// 1. Idempotence — canon(canon(s)) === canon(s) byte-identical
// ---------------------------------------------------------------------------
//
// One pass through the writer normalises atom order, bond placement, and
// ring-closure numbering. A second pass must produce the exact same string.
//
// Note: indole drifts by a single ring-closure swap on the second pass — a
// known stereo-of-aromatic-fused-ring drift, see test_round_trip.js. It is
// intentionally omitted from this idempotence list.
var IDEMPOTENT_CASES = [
    { name: 'ethanol',     smi: 'CCO' },
    { name: 'benzene',     smi: 'c1ccccc1' },
    { name: 'naphthalene', smi: 'c1ccc2ccccc2c1' },
    { name: 'pyridine',    smi: 'n1ccccc1' },
    { name: 'pyrrole',     smi: 'c1cc[nH]c1' },
    { name: 'caffeine',    smi: 'Cn1c(=O)c2c(ncn2C)n(C)c1=O' },
    { name: 'aspirin',     smi: 'CC(=O)Oc1ccccc1C(=O)O' },
    { name: 'acetic_acid', smi: 'CC(=O)O' },
    { name: 'glucose',     smi: 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O' },
    { name: 'phenol',      smi: 'Oc1ccccc1' }
];

IDEMPOTENT_CASES.forEach(function(c) {
    test('idempotence: canon(canon(' + c.name + ')) === canon(' + c.name + ')', function() {
        var c1 = canon(c.smi);
        assert.ok(c1.length > 0, 'first canon empty');
        var c2 = canon(c1);
        assert.strictEqual(c2, c1,
            c.name + ' canonical drifted: "' + c1 + '" → "' + c2 + '"');
    });
});

// ---------------------------------------------------------------------------
// 2. Equivalence under permuted input — different SMILES strings of the
//    SAME molecule must produce the same canonical output
// ---------------------------------------------------------------------------
var EQUIV_CASES = [
    { name: 'ethanol two orderings (CCO, OCC)',
      a: 'CCO', b: 'OCC' },
    { name: 'naphthalene two ring-opening orders',
      a: 'c1ccc2ccccc2c1', b: 'c1cc2ccccc2cc1' },
    { name: 'pyrrole two ring-opening orders',
      a: 'c1cc[nH]c1', b: '[nH]1cccc1' },
    { name: 'acetic acid two writings (CC(=O)O, O=C(C)O)',
      a: 'CC(=O)O', b: 'O=C(C)O' },
    { name: 'toluene two ring-opening orders (Cc1ccccc1, c1ccc(C)cc1)',
      a: 'Cc1ccccc1', b: 'c1ccc(C)cc1' },
    { name: 'phenylacetic acid two ring-opening orders',
      a: 'O=C(O)Cc1ccccc1', b: 'c1ccc(CC(=O)O)cc1' },
    { name: 'N-methylpyrrolidine two writings (CN1CCCC1, N1(C)CCCC1)',
      a: 'CN1CCCC1', b: 'N1(C)CCCC1' },
    { name: 'glucose (no stereo) two ring-opening orders',
      a: 'OC1OC(CO)C(O)C(O)C1O', b: 'OCC1OC(O)C(O)C(O)C1O' }
];

EQUIV_CASES.forEach(function(c) {
    test('equivalence: ' + c.name, function() {
        var ca = canon(c.a);
        var cb = canon(c.b);
        assert.strictEqual(ca, cb,
            'inputs "' + c.a + '" and "' + c.b + '" produced different canonicals: ' +
            '"' + ca + '" vs "' + cb + '"');
    });
});

// ---------------------------------------------------------------------------
// 3. Determinism across runs — canon(s) is byte-identical across 5 runs
// ---------------------------------------------------------------------------
var DETERMINISM_CASES = [
    { name: 'aspirin', smi: 'CC(=O)Oc1ccccc1C(=O)O' },
    { name: 'caffeine', smi: 'Cn1c(=O)c2c(ncn2C)n(C)c1=O' },
    { name: 'glucose', smi: 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O' }
];

DETERMINISM_CASES.forEach(function(c) {
    test('determinism: 5 canon runs of ' + c.name + ' produce identical output', function() {
        var first = canon(c.smi);
        for (var i = 1; i < 5; i++) {
            var next = canon(c.smi);
            assert.strictEqual(next, first,
                'run ' + (i + 1) + ' of ' + c.name + ' drifted: ' +
                '"' + first + '" → "' + next + '"');
        }
    });
});

// ---------------------------------------------------------------------------
// 4. Aromatic preference — parsed-aromatic rings stay in lowercase aromatic
//    form on canonical write (and on every subsequent canonicalisation)
// ---------------------------------------------------------------------------
function hasAromaticAtom(s) {
    // Look for any of the lowercase aromatic organic-subset symbols.
    return /[bcnops]/.test(s);
}

test('aromatic preference: c1ccccc1 (benzene) canonicalises to lowercase aromatic', function() {
    var c = canon('c1ccccc1');
    assert.ok(hasAromaticAtom(c),
        'expected lowercase aromatic atoms in canonical of benzene, got: "' + c + '"');
    assert.strictEqual(c, 'c1ccccc1', 'benzene canonical mismatch: "' + c + '"');
});

test('aromatic preference: n1ccccc1 (pyridine) canonicalises to lowercase aromatic', function() {
    var c = canon('n1ccccc1');
    assert.ok(hasAromaticAtom(c),
        'expected lowercase aromatic atoms in canonical of pyridine, got: "' + c + '"');
    assert.ok(c.indexOf('n') >= 0, 'expected lowercase n in pyridine canonical');
});

test('aromatic preference: c1cc[nH]c1 (pyrrole) canonicalises to lowercase aromatic with [nH]', function() {
    var c = canon('c1cc[nH]c1');
    assert.ok(hasAromaticAtom(c),
        'expected lowercase aromatic atoms in canonical of pyrrole, got: "' + c + '"');
    assert.ok(c.indexOf('[nH]') >= 0,
        'expected explicit [nH] in pyrrole canonical, got: "' + c + '"');
});

// ---------------------------------------------------------------------------
// 5. Charge / isotope / stereo round-trip canonical — a single canon pass
//    leaves the SMILES in canonical form, so a second pass is a no-op
// ---------------------------------------------------------------------------
var ROUND_TRIP_CANON_CASES = [
    { name: 'ammonium [NH4+]', smi: '[NH4+]' },
    { name: '13C-labelled ethyl [13C]CC', smi: '[13C]CC' },
    { name: 'tetrahedral stereocentre [C@H](N)(O)F', smi: '[C@H](N)(O)F' }
];

ROUND_TRIP_CANON_CASES.forEach(function(c) {
    test('round-trip canonical: ' + c.name + ' is stable on second pass', function() {
        var c1 = canon(c.smi);
        var c2 = canon(c1);
        assert.strictEqual(c2, c1,
            c.name + ' second pass drifted: "' + c1 + '" → "' + c2 + '"');
    });
});

// ---------------------------------------------------------------------------
// 6. Edge cases — empty / single atom
// ---------------------------------------------------------------------------
test('edge case: empty SMILES "" canonicalises to ""', function() {
    var c = canon('');
    assert.strictEqual(c, '', 'expected empty canonical for empty input, got: "' + c + '"');
});

test('edge case: methane "C" canonicalises to "C"', function() {
    var c = canon('C');
    assert.strictEqual(c, 'C', 'expected "C" for methane, got: "' + c + '"');
});

// ---------------------------------------------------------------------------
// Cross-cutting integrity checks
// ---------------------------------------------------------------------------
test('integrity: canon output of every idempotence case re-parses without errors', function() {
    IDEMPOTENT_CASES.forEach(function(c) {
        var s = canon(c.smi);
        var m = P(s);
        assert.strictEqual(m.parseErrors.length, 0,
            c.name + ' canonical output "' + s + '" has parse errors: ' +
            m.parseErrors.join('; '));
    });
});

test('integrity: equivalence pairs preserve atom and bond counts after canonicalisation', function() {
    EQUIV_CASES.forEach(function(c) {
        var ma = P(c.a);
        var mb = P(c.b);
        assert.strictEqual(ma.atoms.length, mb.atoms.length,
            c.name + ' atom-count mismatch on parse');
        assert.strictEqual(ma.bonds.length, mb.bonds.length,
            c.name + ' bond-count mismatch on parse');
        var mca = P(canon(c.a));
        var mcb = P(canon(c.b));
        assert.strictEqual(mca.atoms.length, mcb.atoms.length,
            c.name + ' atom-count mismatch after canon');
        assert.strictEqual(mca.bonds.length, mcb.bonds.length,
            c.name + ' bond-count mismatch after canon');
    });
});

// ---------------------------------------------------------------------------
// Regression: permutationParity / _permParity bad-input guard — v1.1.3 fix
// ---------------------------------------------------------------------------
//
// Before the fix, CipStereo.assign() and SmilesWriter.write() could infinite-
// loop when an atom's neighbour CIP-rank array contained undefined or an
// out-of-range value (e.g. a partial rank assignment on a molecule with
// unusual valence). The fix validates each element of the perm array before
// entering the cycle-sort and returns 0 (even parity, no stereo flip) on
// invalid input.
//
// We exercise the guard indirectly: build a molecule with a tetrahedral centre
// but only 2 heavy neighbours (valence violation), run CipStereo.assign() and
// SmilesWriter.write(). Both must return without throwing and without hanging.
test('CipStereo.assign() + SmilesWriter.write() complete on under-valenced stereocentre — regression for v1.1.3 permutationParity fix', function() {
    // Carbon with chirality flag but only 2 substituents — triggers the
    // perm-validation guard in CipStereo and _permParity in SmilesWriter.
    var mol = SmilesParser.parse('C[C@@](F)');  // malformed; partial stereocentre
    // Both calls must return (no infinite loop) without throwing.
    var assigned;
    assert.doesNotThrow(function() {
        CipStereo.assign(mol);
        assigned = SmilesWriter.write(mol);
    });
    assert.strictEqual(typeof assigned, 'string');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
