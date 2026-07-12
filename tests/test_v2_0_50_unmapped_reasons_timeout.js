/**
 * tests/test_v2_0_50_unmapped_reasons_timeout.js — 'timeout' reason
 * category in unmappedReasons (v2.0.50).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.50 extends v2.0.49's per-atom unmapped-reason classifier with a
 * fifth category, 'timeout'. When `result.timedOut === true` (the mapper
 * aborted on its options.timeoutMs budget), atoms that would otherwise
 * have been classified 'no-isomorphic-partner' or 'orphan-after-rescue'
 * are promoted to 'timeout' — the structural-extension verdict isn't
 * trustworthy because the loop bailed early.
 *
 * The two stoichiometric verdicts ('no-partner-element' and
 * 'cofactor-side-only') are NOT overridden by a timeout: they're
 * computed from per-element multiplicities up front, independent of
 * the strategy loop, so a timeout doesn't invalidate them.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('Unmapped reasons timeout (v2.0.50)');
var test = runner.test;

console.log('Unmapped reasons timeout (v2.0.50)');

// Synthetic fake-result builder so tests don't depend on the mapper
// actually timing out (which would be slow and timing-flaky).
function fakeMol(atomSpecs, idBase) {
    return {
        atoms: atomSpecs.map(function (s, i) {
            return { id: idBase + i, symbol: s };
        })
    };
}

// ---------------------------------------------------------------------------
// D. timeout promotion
// ---------------------------------------------------------------------------
test('D1 timedOut=true promotes "no-isomorphic-partner" to "timeout"', function() {
    // 3 reactant C, 3 product C, nothing mapped, strategyResults present.
    // Without timedOut: classification is 'no-isomorphic-partner'.
    // With    timedOut: classification is 'timeout'.
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C', 'C'], 100)],
            products:  [fakeMol(['C', 'C', 'C'], 200)]
        },
        mapping: {},
        strategyResults: [{}],
        timedOut: true
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    assert.strictEqual(reasons[100], 'timeout');
    assert.strictEqual(reasons[101], 'timeout');
    assert.strictEqual(reasons[102], 'timeout');
    assert.strictEqual(reasons[200], 'timeout');
});

test('D2 timedOut=true promotes "orphan-after-rescue" to "timeout" (no strategyResults)', function() {
    // Same shape as D1 but strategyResults is empty/missing — the original
    // classifier would produce 'orphan-after-rescue'. timedOut still promotes
    // it to 'timeout'.
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 300)],
            products:  [fakeMol(['C', 'C'], 400)]
        },
        mapping: {},
        strategyResults: [],          // no results -> orphan branch
        timedOut: true
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    assert.strictEqual(reasons[300], 'timeout');
    assert.strictEqual(reasons[301], 'timeout');
});

test('D3 timedOut=true does NOT override "no-partner-element"', function() {
    // Reactant has Cl; product side has none. Stoichiometric reason, not
    // a timeout symptom — must survive the promotion.
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'Cl'], 500)],
            products:  [fakeMol(['C'], 600)]
        },
        mapping: {},
        strategyResults: [{}],
        timedOut: true
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    // The Cl atom is a stoichiometric mismatch even if the mapper had
    // infinite time: keep the original verdict.
    assert.strictEqual(reasons[501], 'no-partner-element');
});

test('D4 timedOut=true does NOT override "cofactor-side-only"', function() {
    // 2 reactant C, 1 product C, nothing mapped. The two reactant carbons
    // both get 'cofactor-side-only' (stoichiometric). The lone product C
    // gets the structural-extension reason — which, under timeout, promotes
    // to 'timeout'. So this case checks BOTH:
    //   - stoich verdicts survive timeout, and
    //   - the structural verdict on the other side still promotes.
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 700)],
            products:  [fakeMol(['C'], 800)]
        },
        mapping: {},
        strategyResults: [{}],
        timedOut: true
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    assert.strictEqual(reasons[700], 'cofactor-side-only');
    assert.strictEqual(reasons[701], 'cofactor-side-only');
    // Product side: same-element pool exists on the reactant side; mapper
    // timed out before reaching it.
    assert.strictEqual(reasons[800], 'timeout');
});

test('D5 timedOut=false leaves classifications unchanged (back-compat)', function() {
    // The v2.0.49 baseline: when the mapper completed (timedOut absent or
    // false), no atom should be classified 'timeout'. This guards against
    // a regression that silently classifies completed-mapper failures as
    // timeouts.
    var resultA = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 900)],
            products:  [fakeMol(['C', 'C'], 1000)]
        },
        mapping: {},
        strategyResults: [{}]
        // no timedOut field
    };
    var reasonsA = RDT.deriveUnmappedReasons(resultA);
    assert.strictEqual(reasonsA[900], 'no-isomorphic-partner');
    assert.strictEqual(reasonsA[901], 'no-isomorphic-partner');

    var resultB = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 1100)],
            products:  [fakeMol(['C', 'C'], 1200)]
        },
        mapping: {},
        strategyResults: [{}],
        timedOut: false
    };
    var reasonsB = RDT.deriveUnmappedReasons(resultB);
    assert.strictEqual(reasonsB[1100], 'no-isomorphic-partner');
    assert.strictEqual(reasonsB[1101], 'no-isomorphic-partner');
});

test('D6 mapped atoms remain absent under timedOut=true', function() {
    // Sanity: timeout promotion must not introduce 'timeout' entries for
    // atoms that did get mapped. The mapping table still wins.
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C', 'C'], 1300)],
            products:  [fakeMol(['C', 'C', 'C'], 1400)]
        },
        mapping: { 1300: 1400, 1301: 1401 },   // 1302 / 1402 unmapped
        strategyResults: [{}],
        timedOut: true
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    // Mapped atoms absent from reasons.
    assert.strictEqual(typeof reasons[1300], 'undefined');
    assert.strictEqual(typeof reasons[1400], 'undefined');
    // Unmapped atoms get 'timeout'.
    assert.strictEqual(reasons[1302], 'timeout');
    assert.strictEqual(reasons[1402], 'timeout');
});

// ---------------------------------------------------------------------------
// E. integration with the diagnostic
// ---------------------------------------------------------------------------
test('E1 deriveDiagnostic propagates timeout reasons via unmappedReasons', function() {
    // Build a fake mapReaction-shaped result and run deriveDiagnostic.
    // deriveDiagnostic reads result.sides/mapping/timedOut, so a hand-rolled
    // result is sufficient — no actual mapping needed.
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 1500)],
            products:  [fakeMol(['C', 'C'], 1600)]
        },
        mapping: {},
        strategyResults: [{}],
        timedOut: true
    };
    var diag = RDT.deriveDiagnostic(result);
    assert.ok(diag && diag.unmappedReasons, 'diagnostic exposes unmappedReasons');
    assert.strictEqual(diag.unmappedReasons[1500], 'timeout');
    assert.strictEqual(diag.unmappedReasons[1501], 'timeout');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
