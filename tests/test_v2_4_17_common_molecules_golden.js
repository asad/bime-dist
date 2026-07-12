/**
 * tests/test_v2_4_17_common_molecules_golden.js — golden fingerprints for the
 * six common-molecules SMILES corrected in v2.4.17.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Six entries (Theophylline, Guanine, Guanine-Kekulé, Allopurinol, Etodolac,
 * Sitagliptin) used a non-standard aromatic encoding and parsed to the wrong
 * molecular formula.
 * The v2.4.17 corrections were first verified by FORMULA only — which is blind
 * to enantiomers (Sitagliptin: (S) distomer vs (R)) and to constitutional
 * isomers (dimethylxanthines all share C7H8N4O2, so "Theophylline" was actually
 * paraxanthine). This test pins the connectivity + absolute configuration, not
 * just the formula, so those error classes cannot regress unnoticed.
 *
 * Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'CIPStereo.js'));

var SmilesParser = globalThis.SmilesParser;
var CIPStereo = globalThis.CIPStereo;
var Molecule = globalThis.Molecule;

/* eslint-disable no-eval */
// strict-mode eval scopes the file's `var COMMON_MOLECULES` internally, so
// capture it via a trailing expression.
var COMMON_MOLECULES = eval(fs.readFileSync(path.join(__dirname, '..', 'common-molecules.js'), 'utf8') + '\n;COMMON_MOLECULES;');

var runner = shim.makeRunner('Common-molecules golden (v2.4.17)');
var test = runner.test;
console.log('Common-molecules golden (v2.4.17)');

function entry(name) {
    for (var i = 0; i < COMMON_MOLECULES.length; i++) {
        var e = COMMON_MOLECULES[i];
        var n = Array.isArray(e) ? e[0] : (e.name || e.n);
        if (n === name) { return { name: n, smiles: Array.isArray(e) ? e[1] : (e.smiles || e.smi) }; }
    }
    return null;
}
function formulaOf(m) { return Molecule.formulaString(m.elementCounts()); }
function cipMultiset(m) { CIPStereo.assignRS(m); return m.atoms.map(function (a) { return a.cipLabel; }).filter(Boolean).sort().join(','); }

// name -> { formula, cip } — formula + sorted CIP R/S multiset, both verified
// against the authoritative reference compound records.
var GOLD = {
    'Theophylline':     { formula: 'C7H8N4O2', cip: '' },
    'Guanine':          { formula: 'C5H5N5O',  cip: '' },
    'Guanine (Kekule)': { formula: 'C5H5N5O',  cip: '' },
    'Allopurinol':      { formula: 'C5H4N4O',  cip: '' },
    'Etodolac':         { formula: 'C17H21NO3', cip: '' },
    'Sitagliptin':      { formula: 'C16H15F6N5O', cip: 'R' }
};

Object.keys(GOLD).forEach(function (name) {
    var gold = GOLD[name];
    test(name + ' has the correct formula + absolute configuration', function () {
        var e = entry(name);
        assert.ok(e && e.smiles, name + ' is in common-molecules.js');
        var m = SmilesParser.parse(e.smiles);
        assert.strictEqual((m.parseErrors || []).length, 0, name + ' parses cleanly');
        assert.strictEqual(formulaOf(m), gold.formula, name + ' molecular formula');
        assert.strictEqual(cipMultiset(m), gold.cip, name + ' CIP R/S multiset (absolute config)');
    });
});

test('Sitagliptin is the (R) eutomer, not the (S) distomer', function () {
    var m = SmilesParser.parse(entry('Sitagliptin').smiles);
    CIPStereo.assignRS(m);
    var labels = m.atoms.map(function (a) { return a.cipLabel; }).filter(Boolean);
    assert.deepStrictEqual(labels, ['R'], 'the single stereocentre must be R (the active enantiomer)');
});

test('Theophylline is 1,3-dimethyl (imidazole N-H), NOT paraxanthine (1,7-dimethyl)', function () {
    // Theophylline: both methyls on the pyrimidine ring N, the imidazole N
    // carries an H. Paraxanthine puts a methyl on the imidazole N. Count the
    // N-methyl and N-H nitrogens to distinguish them despite the shared formula.
    var m = SmilesParser.parse(entry('Theophylline').smiles);
    var nMethyl = 0, nH = 0;
    m.atoms.forEach(function (a) {
        if (a.symbol !== 'N') { return; }
        var nbrs = m.getBondsOfAtom(a.id).map(function (b) { return m.getAtom(b.otherAtom(a.id)); });
        if (nbrs.some(function (x) { return x.symbol === 'C' && m.getBondsOfAtom(x.id).length === 1; })) { nMethyl++; }
        if (m.calcHydrogens(a.id) > 0) { nH++; }
    });
    assert.strictEqual(nMethyl, 2, 'theophylline has exactly two N-methyl groups');
    assert.strictEqual(nH, 1, 'theophylline has exactly one N-H (the imidazole N7-H)');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
