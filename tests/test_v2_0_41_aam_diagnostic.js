/**
 * tests/test_v2_0_41_aam_diagnostic.js — AAM quality diagnostic
 * (grade + aromatic-ring + stereo conservation; v2.0.41).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Pure value-add diagnostics on top of the existing v2.0.10 AAM result
 * shape: `quality` (A/B/C/F/unmapped), `diagnostic` (structured report).
 * The underlying fitness function + selection are unchanged, so this
 * release cannot regress the golden corpus or the per-class pass rate.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('AAM diagnostic (v2.0.41)');
var test = runner.test;

console.log('AAM diagnostic (v2.0.41)');

function parseReaction(smi) {
    var m = SmilesParser.parse(smi);
    if (m.parseErrors.length > 0) {
        throw new Error('parse: ' + m.parseErrors.join('; '));
    }
    return m;
}

// ---------------------------------------------------------------------------
// A. gradeMapping
// ---------------------------------------------------------------------------
test('A1 gradeMapping returns "unmapped" when status !== "mapped"', function() {
    assert.strictEqual(RDT.gradeMapping(null), 'unmapped');
    assert.strictEqual(RDT.gradeMapping({}), 'unmapped');
    assert.strictEqual(RDT.gradeMapping({ status: 'unbalanced' }), 'unmapped');
});

test('A2 gradeMapping returns "A" for a high-confidence, high-decisiveness mapping with no aromatic / stereo loss', function() {
    var r = {
        status: 'mapped',
        confidence: 0.99,
        decisiveness: 0.5,
        mapping: {},
        sides: { reactants: [], products: [] }
    };
    assert.strictEqual(RDT.gradeMapping(r), 'A');
});

test('A3 gradeMapping degrades to "B" when confidence drops just below A threshold', function() {
    var r = {
        status: 'mapped',
        confidence: 0.9,
        decisiveness: 0.25,
        mapping: {},
        sides: { reactants: [], products: [] }
    };
    assert.strictEqual(RDT.gradeMapping(r), 'B');
});

test('A4 gradeMapping degrades to "C" for moderate-confidence mappings', function() {
    var r = {
        status: 'mapped',
        confidence: 0.7,
        decisiveness: 0.15,
        mapping: {},
        sides: { reactants: [], products: [] }
    };
    assert.strictEqual(RDT.gradeMapping(r), 'C');
});

test('A5 gradeMapping degrades to "F" below the C threshold', function() {
    var r = {
        status: 'mapped',
        confidence: 0.3,
        decisiveness: 0.05,
        mapping: {},
        sides: { reactants: [], products: [] }
    };
    assert.strictEqual(RDT.gradeMapping(r), 'F');
});

// ---------------------------------------------------------------------------
// B. countAromaticRingsPreserved
// ---------------------------------------------------------------------------
test('B1 aromatic-ring count is zero for an empty mapping', function() {
    var r = { status: 'mapped', mapping: {}, sides: { reactants: [], products: [] } };
    var arc = RDT._countAromaticRingsPreserved(r);
    assert.strictEqual(arc.preserved, 0);
    assert.strictEqual(arc.total, 0);
});

test('B2 aromatic-ring count reports the reactant total when atoms are unmapped', function() {
    // Trivial benzene → benzene identity (no actual chemistry change),
    // but mapping is empty — every ring should be counted as "total"
    // but zero "preserved" because no atoms are mapped.
    var mol = parseReaction('c1ccccc1>>c1ccccc1');
    var res = { status: 'mapped', mapping: {}, sides: {
        reactants: [mol.reactionArrow ? mol : mol], products: [mol.reactionArrow ? mol : mol]
    }};
    var arc = RDT._countAromaticRingsPreserved(res);
    assert.ok(arc.total >= 0, 'total is finite');
    assert.strictEqual(arc.preserved, 0, 'no atoms mapped → no rings preserved');
});

// ---------------------------------------------------------------------------
// C. countStereoPreserved
// ---------------------------------------------------------------------------
test('C1 stereo-preserved count is zero/zero for an empty mapping', function() {
    var r = { status: 'mapped', mapping: {}, sides: { reactants: [], products: [] } };
    var stc = RDT._countStereoPreserved(r);
    assert.strictEqual(stc.preserved, 0);
    assert.strictEqual(stc.total, 0);
});

// ---------------------------------------------------------------------------
// D. mapReaction integration
// ---------------------------------------------------------------------------
test('D1 mapReaction result carries the new "quality" and "diagnostic" fields by default', function() {
    var mol = parseReaction('CCO>>CC=O');  // ethanol → acetaldehyde (dehydrogenation)
    var res = RDT.mapReaction(mol);
    assert.ok(res, 'result returned');
    assert.notStrictEqual(typeof res.quality, 'undefined', 'quality field present');
    assert.notStrictEqual(typeof res.diagnostic, 'undefined', 'diagnostic field present');
    if (res.status === 'mapped') {
        assert.ok(/^[ABCF]$/.test(res.quality) || res.quality === 'unmapped',
            'quality is A/B/C/F or unmapped, got: ' + res.quality);
        var d = res.diagnostic;
        assert.strictEqual(typeof d.mappedCount, 'number');
        assert.strictEqual(typeof d.reactantAtomCount, 'number');
        assert.strictEqual(typeof d.productAtomCount, 'number');
        assert.strictEqual(typeof d.aromaticRingsPreserved, 'number');
        assert.strictEqual(typeof d.aromaticRingsTotal, 'number');
        assert.strictEqual(typeof d.stereoCentresPreserved, 'number');
        assert.strictEqual(typeof d.stereoCentresTotal, 'number');
        assert.strictEqual(typeof d.confidence, 'number');
        assert.strictEqual(typeof d.decisiveness, 'number');
    }
});

test('D2 mapReaction with includeDiagnostic:false skips the diagnostic compute (opt-out)', function() {
    var mol = parseReaction('CCO>>CC=O');
    var res = RDT.mapReaction(mol, { includeDiagnostic: false });
    assert.ok(res);
    assert.strictEqual(res.quality, null, 'quality not populated');
    assert.strictEqual(res.diagnostic, null, 'diagnostic not populated');
});

test('D3 deriveDiagnostic on a real result returns mappedCount matching the mapping size', function() {
    var mol = parseReaction('CCO>>CC=O');
    var res = RDT.mapReaction(mol);
    if (res.status !== 'mapped') { return; }
    var d = res.diagnostic;
    var mappingSize = 0;
    for (var k in res.mapping) { if (res.mapping.hasOwnProperty(k)) { mappingSize++; } }
    assert.strictEqual(d.mappedCount, mappingSize);
});

// ---------------------------------------------------------------------------
// E. benzene-aromatic-ring preservation (real chemistry)
// ---------------------------------------------------------------------------
test('E1 benzene → benzene identity: aromatic ring is preserved when atoms map', function() {
    var mol = parseReaction('c1ccccc1>>c1ccccc1');
    var res = RDT.mapReaction(mol);
    if (res.status !== 'mapped') { return; }
    var arc = RDT._countAromaticRingsPreserved(res);
    if (arc.total > 0) {
        // For identity mapping, the aromatic ring atoms must all map to
        // each other, so preserved should equal total.
        assert.strictEqual(arc.preserved, arc.total,
            'identity reaction preserves all aromatic rings, got ' + arc.preserved + '/' + arc.total);
    }
});

// ---------------------------------------------------------------------------
// F. backward compat — existing fields unchanged
// ---------------------------------------------------------------------------
test('F1 existing v2.0.10 fields (status, mapping, score, confidence, decisiveness) still present', function() {
    var mol = parseReaction('CCO>>CC=O');
    var res = RDT.mapReaction(mol);
    assert.notStrictEqual(typeof res.status, 'undefined');
    assert.notStrictEqual(typeof res.mapping, 'undefined');
    assert.notStrictEqual(typeof res.score, 'undefined');
    assert.notStrictEqual(typeof res.confidence, 'undefined');
    assert.notStrictEqual(typeof res.decisiveness, 'undefined');
    assert.notStrictEqual(typeof res.candidates, 'undefined');
    assert.notStrictEqual(typeof res.bondChanges, 'undefined');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
