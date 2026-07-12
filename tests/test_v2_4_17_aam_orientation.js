/**
 * tests/test_v2_4_17_aam_orientation.js — the atom-atom mapper must be correct
 * and machine-independent regardless of how the reactant/product SMILES are
 * written (v2.4.17).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Two defects:
 *   1. bondsCompat() keys off matchBondOrder, not matchBondType, so RDT's
 *      matchBondType:false intent was dropped and single<->double was rejected.
 *      The MCS then could not span the reaction-centre bond-ORDER change of an
 *      oxidation (C-O -> C=O); it stopped at the change and RDT extended the
 *      leftover atoms by an order-dependent heuristic. So ethanol oxidation was
 *      mapped correctly as an order change when written C-first (CCO>>CC=O) but
 *      as a spurious cleaved+formed pair when written O-first (OCC>>CC=O).
 *   2. The winner could depend on a wall-clock timeout budget.
 *
 * Fix: RDT now runs extra bond-order-flexible candidate strategies (MIN_FLEX /
 * MAX_FLEX) that can span the order change and compete purely on score, so the
 * fewer-change (correct) mapping wins for an oxidation, while reactions that
 * genuinely form/cleave a bond (Claisen) keep their strict mapping.
 *
 * Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'RDT.js'));

var SmilesParser = globalThis.SmilesParser;
var RDT = globalThis.RDT;

var runner = shim.makeRunner('AAM orientation independence (v2.4.17)');
var test = runner.test;
console.log('AAM orientation independence (v2.4.17)');

function changeTypes(rxn, opts) {
    var res = RDT.mapReaction(SmilesParser.parse(rxn), opts || {});
    return (res.bondChanges || []).map(function (b) { return b.type; });
}

test('A ethanol oxidation is an order change in ALL four SMILES orientations', function () {
    ['CCO>>CC=O', 'OCC>>CC=O', 'CCO>>O=CC', 'OCC>>O=CC'].forEach(function (rxn) {
        var t = changeTypes(rxn);
        assert.ok(t.indexOf('orderChange') !== -1, rxn + ' must report the C-O -> C=O order change');
        assert.strictEqual(t.indexOf('cleaved'), -1, rxn + ' must NOT report a cleaved bond');
        assert.strictEqual(t.indexOf('formed'), -1, rxn + ' must NOT report a formed bond');
    });
});

test('B the mapping is budget-independent (same result at 1ms and 60000ms)', function () {
    ['CCO>>CC=O', 'OCC>>CC=O'].forEach(function (rxn) {
        var fast = changeTypes(rxn, { timeoutMs: 1 }).slice().sort().join(',');
        var slow = changeTypes(rxn, { timeoutMs: 60000 }).slice().sort().join(',');
        assert.strictEqual(fast, slow, rxn + ' bond-change set must not depend on the timeout budget');
    });
});

test('C a genuine bond FORMATION is still detected (Claisen keeps its strict map)', function () {
    // Two ethyl acetates -> ethyl acetoacetate + ethanol: a new C-C bond forms.
    // The flex candidate over-spans this (more changes, worse score) and loses,
    // so the strict mapping that reports the formed bond must still win.
    var t = changeTypes('CCOC(=O)C.CCOC(=O)C>>CCOC(=O)CC(=O)C.OCC');
    assert.ok(t.indexOf('formed') !== -1, 'the new C-C bond in the Claisen condensation must be reported as formed');
});

test('D a genuine bond CLEAVAGE is still detected (ring opening)', function () {
    // Cyclopropane -> propene breaks a ring C-C bond.
    var t = changeTypes('C1CC1>>C=CC');
    assert.ok(t.indexOf('cleaved') !== -1, 'the ring bond opening must still be reported as cleaved');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
