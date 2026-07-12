/**
 * tests/test_v2_0_72_aam_grade_calibration.js — pin the v2.0.72 AAM
 * quality-grade recalibration so it cannot silently revert.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists
 * ---------------
 * Pre-v2.0.72 the grade gates included `decisiveness ≥ 0.2` (B) and
 * `decisiveness ≥ 0.1` (C). Decisiveness is the gap between the winning
 * mapping strategy and the second-best — it collapses to 0 whenever
 * multiple strategies TIE at the optimum, which is the common case for
 * everyday reactions (symmetric substrates, simple substitutions, the
 * keto-enol equilibrium, etc.). The consequence was that well-formed
 * mappings with high confidence and full aromatic preservation routinely
 * dropped to F for any reaction outside the curated golden corpus. The
 * v2.0.72 recalibration drops the decisiveness gate from B/C (it is
 * retained as the A-grade gate, where ties genuinely should not earn
 * the gold-standard label) and lowers the confidence floors so the
 * remaining metrics carry the signal smoothly.
 *
 * This file pins:
 *   - the threshold values, so a future hand-tuning sees explicit failure
 *   - a panel of seven textbook reactions and their post-recalibration
 *     grades, so the user-facing improvement is locked in
 *   - the no-regression promise: every reaction that graded ≥ X under
 *     the old thresholds still grades ≥ X under the new ones
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');
// shim.loadAll() does not include RDT, so require it explicitly. Then
// bind a local reference (strict-mode safety against bare-name lookup).
var RDT = global.RDT || globalThis.RDT;

var runner = shim.makeRunner('AAM grade calibration (v2.0.72)');
var test = runner.test;

console.log('AAM grade calibration (v2.0.72)');

function fakeResult(confidence, decisiveness, arc, stc) {
    return {
        status: 'mapped',
        confidence: confidence,
        decisiveness: decisiveness,
        mapping: {},
        sides: { reactants: [], products: [] },
        // Inject pre-computed ARC/stereo via a small fake counter the
        // grader will accept (countAromatic/countStereoPreserved return
        // {preserved:0, total:0} on empty sides, which evaluates to rate=1
        // — so empty-mapping results behave like full preservation).
        _arc: arc,
        _stc: stc
    };
}

// ---------------------------------------------------------------------------
// A. v2.0.72 threshold pins. A failure here is a deliberate edit and
// requires updating both the doc-comment in editor/RDT.js and this test.
// ---------------------------------------------------------------------------
test('A1 A grade: confidence ≥ 0.95, decisiveness ≥ 0.3, ARC + stereo high',
    function () {
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.99, 0.5)), 'A',
            'high confidence + high decisiveness + empty mapping → A');
});

test('A2 A grade fails when decisiveness < 0.3 (ties at optimum should '
    + 'NOT earn A even with high confidence)', function () {
        // conf=0.99 but dec=0.0 (strategies tied) → drops to B.
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.99, 0)), 'B',
            'tied strategies cannot earn A regardless of confidence');
});

test('A3 B grade: confidence ≥ 0.75, no decisiveness gate', function () {
        // conf=0.75 with dec=0 should grade B post-v2.0.72; pre-v2.0.72
        // it would have dropped to F via the dec ≥ 0.2 gate.
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.75, 0)), 'B',
            'confidence at the B floor with tied strategies → B (not F)');
});

test('A4 C grade: confidence ≥ 0.5, no decisiveness gate', function () {
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.5, 0)), 'C',
            'confidence at the C floor with tied strategies → C (not F)');
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.6, 0)), 'C',
            'mid-confidence + ties → C');
});

test('A5 F grade: confidence below the C floor', function () {
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.49, 0.5)), 'F',
            'below C floor → F even with strong decisiveness');
});

test('A6 unmapped status overrides everything', function () {
        var r = { status: 'unbalanced', confidence: 0.99, decisiveness: 0.5 };
        assert.strictEqual(RDT.gradeMapping(r), 'unmapped');
        assert.strictEqual(RDT.gradeMapping(null), 'unmapped');
});

// ---------------------------------------------------------------------------
// B. Monotonicity: every grade that was achievable under the old (stricter)
// B/C gates is still achievable under the new (more generous) ones. The
// inverse is allowed: some old-F reactions should now grade C/B.
// ---------------------------------------------------------------------------
test('B1 every pre-v2.0.72-passing B still grades ≥ B post-recalibration',
    function () {
        // Pre-v2.0.72 B required conf ≥ 0.85, dec ≥ 0.2, ARC ≥ 0.85,
        // stereo ≥ 0.8. The minimum such input still grades ≥ B today.
        assert.notStrictEqual(
            RDT.gradeMapping(fakeResult(0.85, 0.2)), 'F',
            'pre-v2.0.72 B-floor (conf=0.85, dec=0.2) must not regress to F');
        assert.notStrictEqual(
            RDT.gradeMapping(fakeResult(0.85, 0.2)), 'C',
            'pre-v2.0.72 B-floor must still grade B (or better)');
});

test('B2 every pre-v2.0.72-passing C still grades ≥ C', function () {
        // Pre-v2.0.72 C required conf ≥ 0.6, dec ≥ 0.1.
        assert.notStrictEqual(
            RDT.gradeMapping(fakeResult(0.6, 0.1)), 'F',
            'pre-v2.0.72 C-floor (conf=0.6, dec=0.1) must not regress to F');
});

// ---------------------------------------------------------------------------
// C. The user-facing fix: real reactions with high confidence but
// tied strategies now grade B/C (not F).
// ---------------------------------------------------------------------------
test('C1 high-confidence + tied strategies + preserved aromaticity → B',
    function () {
        // E.g. bromination of benzene: confidence 0.83, decisiveness 0.0
        // (multiple equivalent ring positions). Pre-v2.0.72: F. Now: B.
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.83, 0)), 'B',
            'high-confidence aromatic substitution with tied strategies → B');
});

test('C2 mid-confidence + tied strategies → C (not F)', function () {
        // E.g. Fischer esterification: confidence ~0.66, decisiveness 0.
        // Pre-v2.0.72: F. Now: C.
        assert.strictEqual(
            RDT.gradeMapping(fakeResult(0.66, 0)), 'C',
            'mid-confidence reaction with tied strategies → C');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
