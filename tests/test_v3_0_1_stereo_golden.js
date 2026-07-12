/**
 * tests/test_v3_0_1_stereo_golden.js — pin the authoritative stereochemistry of
 * the 8 library molecules that gained isomeric stereo in v3.0.1.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v3.0.0 corrected 38 malformed-aromatic library molecules; 8 of them were left
 * flat (constitution correct, stereo unassigned). v3.0.1 added each molecule's
 * full authoritative stereochemistry — every one independently confirmed to match
 * its authoritative InChIKey (all three blocks, including the stereo layer).
 *
 * The library-wide valence test only asserts "chemically possible", and the
 * common-molecules golden test covers a different set — so nothing pinned the
 * ACTUAL stereo of these 8. This test locks it: for each, the sorted CIP R/S
 * multiset (tetrahedral centres) AND the sorted E/Z double-bond multiset must
 * match the authoritative configuration, so a future edit that flips a centre or
 * a double-bond geometry (e.g. an all-trans retinol turned 13-cis, or a
 * corticosteroid epimerised at C17) fails immediately. Dependency-free: uses
 * BIME's own CIP engine (no external toolkit), mirroring the v2.4.17
 * Sitagliptin "eutomer not distomer" golden pin.
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
var COMMON_MOLECULES = eval(
    fs.readFileSync(path.join(__dirname, '..', 'common-molecules.js'), 'utf8') +
    '\n;COMMON_MOLECULES;');

var runner = shim.makeRunner('Stereo golden (v3.0.1)');
var test = runner.test;
console.log('Stereo golden (v3.0.1)');

// name -> { formula, cip (sorted R/S multiset), ez (sorted E/Z multiset) }.
// Each value is the authoritative configuration, cross-checked against the
// compound's full InChIKey stereo layer when the SMILES was assigned.
var GOLD = {
    'Digoxin':            { formula: 'C41H64O14',   cip: 'R,R,R,R,R,R,R,R,R,S,S,S,S,S,S,S,S,S,S,S,S', ez: '' },
    'Beclomethasone':     { formula: 'C22H29ClO5',  cip: 'R,R,S,S,S,S,S,S', ez: '' },
    'Doxycycline':        { formula: 'C22H24N2O8',  cip: 'R,R,R,S,S,S',     ez: '' },
    'Prednisolone':       { formula: 'C21H28O5',    cip: 'R,R,S,S,S,S,S',   ez: '' },
    'Dexamethasone':      { formula: 'C22H29FO5',   cip: 'R,R,R,S,S,S,S,S', ez: '' },
    'Triamcinolone':      { formula: 'C21H27FO6',   cip: 'R,R,S,S,S,S,S,S', ez: '' },
    'Ethinyl Estradiol':  { formula: 'C20H24O2',    cip: 'R,R,S,S,S',       ez: '' },
    'Vitamin A (Retinol)':{ formula: 'C20H30O',     cip: '',                ez: 'E,E,E,E' }
};

function entry(name) {
    for (var i = 0; i < COMMON_MOLECULES.length; i++) {
        var e = COMMON_MOLECULES[i];
        var n = Array.isArray(e) ? e[0] : (e.name || e.n);
        if (n === name) { return Array.isArray(e) ? e[1] : (e.smiles || e.smi); }
    }
    return null;
}
function cipMultiset(m) {
    CIPStereo.assignRS(m);
    return m.atoms.map(function (a) { return a.cipLabel; }).filter(Boolean).sort().join(',');
}
function ezMultiset(m) {
    var out = [];
    for (var i = 0; i < m.bonds.length; i++) {
        var b = m.bonds[i];
        if (b.type === Molecule.BOND_DOUBLE && typeof CIPStereo.doubleBondEZ === 'function') {
            var v = CIPStereo.doubleBondEZ(m, b);
            if (v === 'E' || v === 'Z') { out.push(v); }
        }
    }
    return out.sort().join(',');
}

Object.keys(GOLD).forEach(function (name) {
    var gold = GOLD[name];
    test(name + ' has the authoritative stereochemistry (formula + CIP + E/Z)', function () {
        var smi = entry(name);
        assert.ok(smi, name + ' is present in common-molecules.js');
        var m = SmilesParser.parse(smi);
        assert.strictEqual((m.parseErrors || []).length, 0, name + ' parses cleanly');
        assert.strictEqual(Molecule.formulaString(m.elementCounts()), gold.formula,
            name + ' molecular formula (constitution)');
        assert.strictEqual(cipMultiset(m), gold.cip,
            name + ' CIP R/S multiset (tetrahedral absolute configuration)');
        assert.strictEqual(ezMultiset(m), gold.ez,
            name + ' E/Z double-bond multiset (geometry)');
    });
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
