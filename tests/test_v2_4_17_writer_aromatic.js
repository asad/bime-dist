/**
 * tests/test_v2_4_17_writer_aromatic.js — the SMILES writer must not DE-aromatise
 * (saturate) an aromatic system on round-trip (v2.4.17, partial).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background: SmilesParser stores aromatic ring bonds as single (deferring
 * Kekulisation) and never assigns Kekulé orders; the writer re-perceives
 * aromaticity with a conservative Hückel pass that cannot re-derive every
 * aromatic system BIME can PARSE (tetrazoles, azulene, some fused azines, the
 * nucleobase lactams). Those atoms were emitted as saturated uppercase-single,
 * so parse -> write -> parse silently changed the molecular formula: uracil ->
 * dihydrouracil, azulene C10H8 -> C10H18, the sartan tetrazoles gaining 5 H.
 * Across the shipped 1181-molecule library this affected roughly 100 entries.
 *
 * Fix (SmilesWriter.perceiveAromaticity): for an atom the parser flagged
 * aromatic whose bonds are still stored single, trust that flag rather than the
 * incomplete re-perception. calcHydrogens already treats these atoms via
 * AROMATIC_VALENCE, so the formula is recovered exactly. (A residual class —
 * Kekulé-input aromatic N-H rings re-aromatised without an [nH] — is tracked
 * separately.)
 *
 * These tests pin formula conservation for the affected teaching-relevant
 * aromatics. Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();

var SmilesParser = globalThis.SmilesParser;
var SmilesWriter = globalThis.SmilesWriter;
var Molecule = globalThis.Molecule;

var runner = shim.makeRunner('Writer aromatic conservation (v2.4.17)');
var test = runner.test;
console.log('Writer aromatic conservation (v2.4.17)');

function formulaOf(m) {
    var counts = (typeof m.elementCounts === 'function') ? m.elementCounts() : null;
    return counts ? Molecule.formulaString(counts) : '?';
}
function roundtripFormula(smi) {
    return formulaOf(SmilesParser.parse(SmilesWriter.write(SmilesParser.parse(smi))));
}

// name -> [smiles, expected Hill formula]. Each of these USED to saturate on
// write before the v2.4.17 writer fix.
var CASES = {
    'azulene':          ['c1ccc2cccc2cc1', 'C10H8'],
    'uracil':           ['O=c1cc[nH]c(=O)[nH]1', 'C4H4N2O2'],
    'cytosine':         ['Nc1cc[nH]c(=O)n1', 'C4H5N3O'],
    'thymine':          ['Cc1c[nH]c(=O)[nH]c1=O', 'C5H6N2O2'],
    'tetrazole (sartan motif)': ['c1ccc(-c2ccccc2c2nnn[nH]2)cc1', 'C13H10N4'],
    'imidazole':        ['c1c[nH]cn1', 'C3H4N2'],
    'benzimidazole':    ['c1ccc2[nH]cnc2c1', 'C7H6N2'],
    'pyrrole':          ['c1cc[nH]c1', 'C4H5N'],
    'indole':           ['c1ccc2[nH]ccc2c1', 'C8H7N']
};
// Benzenoids / standard heterocycles must remain conserved too (no regression).
var CONTROLS = {
    'benzene':      ['c1ccccc1', 'C6H6'],
    'pyridine':     ['c1ccncc1', 'C5H5N'],
    'furan':        ['c1ccoc1', 'C4H4O'],
    'thiophene':    ['c1ccsc1', 'C4H4S'],
    'naphthalene':  ['c1ccc2ccccc2c1', 'C10H8'],
    'caffeine':     ['Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C8H10N4O2']
};

test('A previously-saturated aromatics keep their formula on parse->write->parse', function () {
    Object.keys(CASES).forEach(function (name) {
        var smi = CASES[name][0], want = CASES[name][1];
        assert.strictEqual(formulaOf(SmilesParser.parse(smi)), want, name + ' parses to ' + want);
        assert.strictEqual(roundtripFormula(smi), want, name + ' must not saturate on round-trip');
    });
});

test('B benzenoids + standard heterocycles remain conserved (no regression)', function () {
    Object.keys(CONTROLS).forEach(function (name) {
        assert.strictEqual(roundtripFormula(CONTROLS[name][0]), CONTROLS[name][1],
            name + ' stays formula-stable');
    });
});

test('C azulene is not saturated to C10H18', function () {
    var w = SmilesWriter.write(SmilesParser.parse('c1ccc2cccc2cc1'));
    assert.strictEqual(roundtripFormula('c1ccc2cccc2cc1'), 'C10H8',
        'azulene round-trips as C10H8, not C10H18; written=' + w);
});

test('D Kekulé-input aromatic N-H rings keep the NH (indole/serotonin -> [nH])', function () {
    // Parsed as Kekulé (uppercase N, explicit C=C); the writer re-aromatises and
    // must emit [nH], not bare n, or the reader (bare n = pyridine, 0 H) drops it.
    var cases = {
        'serotonin':       ['c12cc(O)ccc2NC=C1CCN', 'C10H12N2O'],
        'indole (Kekulé)': ['C1=CC2=CC=CC=C2N1', 'C8H7N'],
        'tryptamine':      ['c12ccccc2NC=C1CCN', 'C10H12N2']
    };
    Object.keys(cases).forEach(function (name) {
        assert.strictEqual(roundtripFormula(cases[name][0]), cases[name][1],
            name + ' keeps its NH on round-trip');
    });
});

test('E the previously-malformed library purines parse to the correct formula', function () {
    // These 4 used the non-standard c1=2ncnc=2 encoding and parsed to a wrong H
    // count (e.g. theophylline C7H7N4O2 instead of C7H8N4O2); corrected in v2.4.17.
    var lib = {
        'Cn1cnc2c1c(=O)n(C)c(=O)[nH]2': 'C7H8N4O2',          // Theophylline
        'NC1=Nc2[nH]cnc2C(=O)N1': 'C5H5N5O',                 // Guanine (Kekulé)
        'O=c1[nH]cnc2[nH]ncc12': 'C5H4N4O',                  // Allopurinol
        'CCc1cccc2c1[nH]c1c2CCOC1(CC)CC(=O)O': 'C17H21NO3'   // Etodolac
    };
    Object.keys(lib).forEach(function (smi) {
        var m = SmilesParser.parse(smi);
        assert.strictEqual((m.parseErrors || []).length, 0, smi + ' parses cleanly');
        assert.strictEqual(formulaOf(m), lib[smi], smi + ' parses to ' + lib[smi]);
        assert.strictEqual(roundtripFormula(smi), lib[smi], smi + ' round-trip stable');
    });
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
