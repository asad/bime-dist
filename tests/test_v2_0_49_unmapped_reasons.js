/**
 * tests/test_v2_0_49_unmapped_reasons.js — per-atom unmapped-reason
 * reporting in AAM (v2.0.49).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * For each atom not present in `result.mapping`, classify the most
 * likely reason. Reason categories: 'no-partner-element',
 * 'cofactor-side-only', 'no-isomorphic-partner', 'orphan-after-rescue'.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('Unmapped reasons (v2.0.49)');
var test = runner.test;

console.log('Unmapped reasons (v2.0.49)');

function parse(smi) {
    var m = SmilesParser.parse(smi);
    if (m.parseErrors.length > 0) {
        throw new Error('parse: ' + m.parseErrors.join('; '));
    }
    return m;
}

// ---------------------------------------------------------------------------
// A. unmappedReasons exists on the diagnostic
// ---------------------------------------------------------------------------
test('A1 result.diagnostic.unmappedReasons is an object (empty when fully mapped)', function() {
    var rxn = parse('CCO>>CC=O');
    var res = RDT.mapReaction(rxn);
    assert.ok(res.diagnostic, 'diagnostic populated');
    assert.strictEqual(typeof res.diagnostic.unmappedReasons, 'object');
    // For a clean dehydrogenation the mapping is full — no unmapped atoms.
    // unmappedReasons should be empty.
    var count = 0;
    for (var k in res.diagnostic.unmappedReasons) {
        if (res.diagnostic.unmappedReasons.hasOwnProperty(k)) { count++; }
    }
    if (res.status === 'mapped' && res.diagnostic.unmappedReactantCount === 0 &&
            res.diagnostic.unmappedProductCount === 0) {
        assert.strictEqual(count, 0, 'no unmapped atoms → empty reasons map');
    }
});

// ---------------------------------------------------------------------------
// B. deriveUnmappedReasons classifier
// ---------------------------------------------------------------------------
test('B1 deriveUnmappedReasons returns an empty map when given an empty result', function() {
    var reasons = RDT.deriveUnmappedReasons({});
    assert.deepStrictEqual(reasons, {});
});

test('B2 deriveUnmappedReasons classifies "no-partner-element" when the other side lacks the element', function() {
    // Synthetic result: reactant has a Cl atom, product side has none.
    // No mapping at all → every reactant atom should be classified.
    var fakeMol = function (atomSpecs) {
        return {
            atoms: atomSpecs.map(function (s, i) {
                return { id: 100 + i, symbol: s };
            })
        };
    };
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'Cl'])],
            products:  [fakeMol(['C'])]
        },
        mapping: {},   // nothing mapped
        strategyResults: [{}]
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    // The Cl atom (id 101) has no partner element on the product side
    // → 'no-partner-element'.
    assert.strictEqual(reasons[101], 'no-partner-element');
});

test('B3 deriveUnmappedReasons classifies "cofactor-side-only" when this side has MORE of the element than the other', function() {
    // Reactant: 2 carbons. Product: 1 carbon. Both unmapped.
    // The reactant carbons get 'cofactor-side-only' (this side has more
    // of the element). Use disjoint id ranges so reactant/product
    // entries don't collide in the reasons map.
    var fakeMol = function (atomSpecs, idBase) {
        return {
            atoms: atomSpecs.map(function (s, i) {
                return { id: idBase + i, symbol: s };
            })
        };
    };
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 200)],
            products:  [fakeMol(['C'], 300)]
        },
        mapping: {},
        strategyResults: [{}]
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    assert.strictEqual(reasons[200], 'cofactor-side-only');
    assert.strictEqual(reasons[201], 'cofactor-side-only');
    // The lone product C has FEWER of the element than the reactant side,
    // so its classification is no-isomorphic-partner (other side has
    // matching element, just couldn't extend a mapping to it).
    assert.strictEqual(reasons[300], 'no-isomorphic-partner');
});

test('B4 deriveUnmappedReasons classifies "no-isomorphic-partner" when same-element atoms exist but mapping skipped', function() {
    // Reactant: 3 carbons. Product: 3 carbons. Mapping empty. Each side has
    // matching multiplicities, but no mapping was made — strategyResults
    // present (mapper TRIED).
    var fakeMol = function (atomSpecs, idBase) {
        return {
            atoms: atomSpecs.map(function (s, i) {
                return { id: idBase + i, symbol: s };
            })
        };
    };
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C', 'C'], 300)],
            products:  [fakeMol(['C', 'C', 'C'], 400)]
        },
        mapping: {},
        strategyResults: [{}]
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    assert.strictEqual(reasons[300], 'no-isomorphic-partner');
    assert.strictEqual(reasons[301], 'no-isomorphic-partner');
});

test('B5 mapped atoms do NOT appear in unmappedReasons', function() {
    var fakeMol = function (atomSpecs, idBase) {
        return {
            atoms: atomSpecs.map(function (s, i) {
                return { id: idBase + i, symbol: s };
            })
        };
    };
    var result = {
        status: 'mapped',
        sides: {
            reactants: [fakeMol(['C', 'C'], 500)],
            products:  [fakeMol(['C', 'C'], 600)]
        },
        mapping: { 500: 600, 501: 601 },
        strategyResults: [{}]
    };
    var reasons = RDT.deriveUnmappedReasons(result);
    assert.strictEqual(typeof reasons[500], 'undefined', 'mapped reactant absent');
    assert.strictEqual(typeof reasons[600], 'undefined', 'mapped product absent');
});

// ---------------------------------------------------------------------------
// C. opt-out via includeDiagnostic: false
// ---------------------------------------------------------------------------
test('C1 includeDiagnostic:false skips unmappedReasons computation', function() {
    var rxn = parse('CCO>>CC=O');
    var res = RDT.mapReaction(rxn, { includeDiagnostic: false });
    assert.strictEqual(res.diagnostic, null);
    // Caller can still invoke deriveUnmappedReasons manually if they want.
    var reasons = RDT.deriveUnmappedReasons(res);
    assert.strictEqual(typeof reasons, 'object');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
