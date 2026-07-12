/**
 * tests/test_v2_4_16_ez_structure_key.js — geometry-free canonical E/Z + the
 * cis/trans mislabel fix for compound labelling (v2.4.16).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background: MetaboliteLibrary.findBySmiles keyed on (constitutional canonical
 * SMILES + atom CIP R/S multiset). That omitted double-bond E/Z, so cis-maleate
 * was mislabelled as its trans partner Fumarate — a wrong name on a figure.
 *
 * Fix:
 *   A. CIPStereo.ezDescriptors(mol) — a canonical, order-invariant, GEOMETRY-FREE
 *      E/Z descriptor per double bond, derived from the parsed directional bonds
 *      (bond.stereo 1='/' 6='\\') + the existing CIP substituent ranking. (The
 *      pre-existing assignEZ reads 2D coordinates, which a re-layout corrupts.)
 *   B. findBySmiles folds the sorted E/Z multiset into the structure key, so
 *      cis and trans never collide; the Fumarate library entry carries its real
 *      (E) configuration so genuine fumarate still resolves.
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
require(path.join(__dirname, '..', 'editor', 'MetaboliteLibrary.js'));

var SmilesParser = globalThis.SmilesParser;
var CIPStereo = globalThis.CIPStereo;
var MetaboliteLibrary = globalThis.MetaboliteLibrary;

var runner = shim.makeRunner('E/Z structure key (v2.4.16)');
var test = runner.test;
console.log('E/Z structure key (v2.4.16)');

var CIP_SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'CIPStereo.js'), 'utf8');
var ML_SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'MetaboliteLibrary.js'), 'utf8');

function ez(smi) { return CIPStereo.ezDescriptors(SmilesParser.parse(smi)).join(','); }
function name(smi) { var e = MetaboliteLibrary.findBySmiles(smi); return e ? e.name : null; }

// ---------------------------------------------------------------------------
// A. CIPStereo.ezDescriptors — correct + geometry-free + order-invariant.
// ---------------------------------------------------------------------------
test('A1 E/Z is read correctly from the directional bonds (no layout needed)', function () {
    assert.strictEqual(ez('OC(=O)/C=C/C(=O)O'), 'E', 'fumarate is (E)');
    assert.strictEqual(ez('OC(=O)/C=C\\C(=O)O'), 'Z', 'maleate is (Z)');
    assert.strictEqual(ez('C/C=C/C'), 'E', '(E)-2-butene');
    assert.strictEqual(ez('C/C=C\\C'), 'Z', '(Z)-2-butene');
});

test('A2 an unspecified double bond yields no descriptor (no guessing)', function () {
    assert.strictEqual(ez('CC=CC'), '', 'plain C=C has no E/Z');
    assert.strictEqual(ez('c1ccccc1'), '', 'aromatic ring contributes no E/Z');
    assert.strictEqual(ez('CC(=O)C(O)=O'), '', 'carbonyls are not E/Z double bonds');
});

test('A3 the descriptor is order-invariant across equivalent writings', function () {
    assert.strictEqual(ez('OC(=O)/C=C/C(=O)O'), ez('O=C(O)/C=C/C(=O)O'), 'fumarate two ways agree');
    // trans written with backslashes is the same molecule -> same descriptor
    assert.strictEqual(ez('OC(=O)\\C=C\\C(=O)O'), 'E', 'backslash-trans is still (E)');
});

// ---------------------------------------------------------------------------
// B. The mislabel is fixed; genuine fumarate still resolves.
// ---------------------------------------------------------------------------
test('B1 cis-maleate is NO LONGER mislabelled as trans-Fumarate', function () {
    assert.strictEqual(name('OC(=O)/C=C\\C(=O)O'), null,
        'maleate (Z) must not resolve to Fumarate (the bug)');
});

test('B2 genuine (E)-fumarate still resolves (incl. its @/@@-style drift form)', function () {
    assert.strictEqual(name('OC(=O)/C=C/C(=O)O'), 'Fumarate', 'fumarate resolves');
    assert.strictEqual(name('O=C(O)/C=C/C(=O)O'), 'Fumarate', 'a re-written fumarate still resolves');
});

test('B3 the Fumarate library entry carries its real (E) configuration', function () {
    var fum = MetaboliteLibrary.find('Fumarate');
    assert.ok(fum && fum.smiles, 'Fumarate is in the library');
    assert.strictEqual(ez(fum.smiles), 'E', 'the stored Fumarate SMILES is (E)');
});

// ---------------------------------------------------------------------------
// C. No regression: every structure entry still self-resolves (no new
//    collisions from the added E/Z key component).
// ---------------------------------------------------------------------------
test('C1 every library metabolite still resolves to itself', function () {
    var mets = MetaboliteLibrary.metabolites.filter(function (m) { return m.smiles; });
    var fail = [];
    mets.forEach(function (m) {
        var e = MetaboliteLibrary.findBySmiles(m.smiles);
        if (!e || e.name !== m.name) { fail.push(m.name + ' -> ' + (e ? e.name : 'null')); }
    });
    assert.strictEqual(fail.length, 0, 'all self-resolve; mismatches: ' + fail.join(' ; '));
    assert.ok(mets.length >= 40, 'a meaningful number of structure entries were checked (' + mets.length + ')');
});

test('C2 achiral + chiral metabolites without C=C stereo are unaffected', function () {
    assert.strictEqual(name('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'), 'D-Glucose');
    assert.strictEqual(name('CC(=O)C(O)=O'), 'Pyruvate');
});

// ---------------------------------------------------------------------------
// D. Source-shape pins.
// ---------------------------------------------------------------------------
test('D1 CIPStereo exports the geometry-free ezDescriptors', function () {
    assert.ok(/ezDescriptors: ezDescriptors/.test(CIP_SRC), 'ezDescriptors is exported');
    assert.ok(/function _directionalSide\(/.test(CIP_SRC), 'directional-side helper present');
    assert.ok(/bond\.stereo|b\.stereo/.test(CIP_SRC), 'E/Z is read from bond.stereo, not coordinates');
});

test('D2 findBySmiles folds the E/Z multiset into the structure key', function () {
    assert.ok(/ezDescriptors\(m\)/.test(ML_SRC), 'the structure key calls ezDescriptors');
    assert.ok(/constitutional \+ '\|' \+ stereo \+ '\|' \+ ez/.test(ML_SRC), 'E/Z is part of the key');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
