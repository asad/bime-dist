/**
 * tests/test_mcs_extended.js — extended MCS regressions for BIME v1.1.1.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Adds drug-pair, disconnected, timeout, bond-consistency, threshold, Tanimoto,
 * symmetry, identity, and single-atom test cases for SMSDMCS.findMCS /
 * SMSDMCS.findMCSFromSmiles.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SMSD MCS (extended)');
var test = runner.test;

console.log('SMSD MCS (extended)');

function mcsSmi(a, b, mcsOpts) {
    return SMSDMCS.findMCSFromSmiles(a, b, null, mcsOpts || {});
}
function size(a, b) { return mcsSmi(a, b).size; }

function buildGraph(smi) {
    var SG = SMSDGraph.SMSDGraph;
    return new SG(SmilesParser.parse(smi));
}

// Count MCS bonds: for every (i,k) pair in g1 with i<k that are both mapped,
// is there a bond between mapping[i] and mapping[k] in g2?
function countMappedBonds(g1, g2, mapping) {
    var n = 0;
    for (var i = 0; i < g1.n; i++) {
        if (!(i in mapping)) continue;
        var nbrs = g1.neighbors[i];
        for (var j = 0; j < nbrs.length; j++) {
            var k = nbrs[j];
            if (k <= i) continue;
            if (!(k in mapping)) continue;
            if (g2.neighbors[mapping[i]].indexOf(mapping[k]) >= 0) n++;
        }
    }
    return n;
}

// =========================================================================
// (1) Real-drug molecules
// =========================================================================
test('aspirin vs salicylic acid → MCS = 9-10 atoms (the salicylate scaffold)', function() {
    var r = mcsSmi('CC(=O)Oc1ccccc1C(=O)O', 'Oc1ccccc1C(=O)O');
    assert.ok(r.size >= 9, 'expected >=9 atoms, got ' + r.size);
    assert.ok(r.size <= 11, 'expected <=11 atoms, got ' + r.size);
});

test('caffeine vs theobromine → MCS >= 12 atoms (purine core)', function() {
    var r = mcsSmi('Cn1c(=O)c2c(ncn2C)n(C)c1=O', 'Cn1c(=O)c2c(nc[nH]2)n(C)c1=O');
    assert.ok(r.size >= 12, 'expected ≥12 atoms, got ' + r.size);
});

test('ibuprofen vs naproxen → MCS >= 10 atoms', function() {
    var r = mcsSmi('CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'COc1ccc2cc(C(C)C(=O)O)ccc2c1');
    assert.ok(r.size >= 10, 'expected ≥10 atoms, got ' + r.size);
});

// =========================================================================
// (2) Disconnected query molecule
// =========================================================================
test('mcs(CC.NN, CCNN) returns >= 2 atoms (some shared fragment)', function() {
    var r = mcsSmi('CC.NN', 'CCNN');
    assert.ok(r.size >= 2, 'expected ≥2 atoms, got ' + r.size);
});

// =========================================================================
// (3) MCS with timeout — finishes regardless on tiny inputs
// =========================================================================
test('mcs(CCO, CCN) with timeoutMs=1 still returns a result', function() {
    var r = mcsSmi('CCO', 'CCN', { timeoutMs: 1 });
    assert.ok(r.size >= 1);
});

// =========================================================================
// (4) Bond-mapping consistency: mapped atoms imply mapped bonds
// =========================================================================
test('aspirin-salicylic-acid MCS has matching bonds in g2', function() {
    var g1 = buildGraph('CC(=O)Oc1ccccc1C(=O)O');
    var g2 = buildGraph('Oc1ccccc1C(=O)O');
    var r = SMSDMCS.findMCS(g1, g2);
    var nb = countMappedBonds(g1, g2, r.mapping);
    // For a connected MCS of N atoms we expect at least N-1 mapped bonds.
    assert.ok(nb >= r.size - 1, 'expected >= ' + (r.size - 1) + ' bonds, got ' + nb);
});

// =========================================================================
// (5) MCS minimum size threshold (mcsOpts.minSize is advisory)
// =========================================================================
// minSize is *not* a hard truncation in the BIME v1.x policy: the engine
// still returns its best discovered size; the option informs early termination
// only when the best-found exceeds the threshold. This test documents that.
test('mcsOpts.minSize is advisory: result.size may be smaller than minSize', function() {
    // Two tiny mols where the true MCS is 2 atoms; even minSize=4 still
    // returns size 2 (the policy is not to lie about what was found).
    var r = mcsSmi('CCO', 'CCN', { minSize: 4 });
    assert.ok(r.size >= 1);
    assert.ok(r.size <= 2);
});

// =========================================================================
// (6) Tanimoto from MCS size
// =========================================================================
test('aspirin vs salicylic acid Tanimoto in (0.5, 1.0)', function() {
    var r = mcsSmi('CC(=O)Oc1ccccc1C(=O)O', 'Oc1ccccc1C(=O)O');
    assert.ok(r.tanimoto > 0.5);
    assert.ok(r.tanimoto < 1.0);
});

test('mcs(c1ccccc1, c1ccccc1) Tanimoto = 1.0 (identical)', function() {
    var r = mcsSmi('c1ccccc1', 'c1ccccc1');
    assert.ok(Math.abs(r.tanimoto - 1.0) < 1e-9);
});

// =========================================================================
// (7) Symmetric MCS: mcs(A,B).size === mcs(B,A).size on 5 pairs
// =========================================================================
test('Symmetric MCS on 5 pairs', function() {
    var pairs = [
        ['c1ccccc1', 'Cc1ccccc1'],
        ['CCO', 'CCN'],
        ['CCC', 'OCC'],
        ['c1ccccc1', 'c1ccncc1'],
        ['C1CCCCC1', 'C1CCCCC1Cl']
    ];
    for (var i = 0; i < pairs.length; i++) {
        var ab = size(pairs[i][0], pairs[i][1]);
        var ba = size(pairs[i][1], pairs[i][0]);
        assert.strictEqual(ab, ba, 'asymmetry on ' + pairs[i][0] + ' / ' + pairs[i][1]);
    }
});

// =========================================================================
// (8) Stereochemistry preservation (v1.x policy: atom-type match is the core
// invariant; CIP labels are a downstream concern)
// =========================================================================
test('mcs([C@H](N)(C)O, [C@@H](N)(C)O) maps all 4 heavy atoms (stereo not enforced at v1.x)', function() {
    // BIME v1.x SMSD does not enforce CIP equivalence by default; this
    // test pins the policy: the inverted chiral center still yields a full
    // 4-atom MCS so downstream callers know what to expect.
    var r = mcsSmi('[C@H](N)(C)O', '[C@@H](N)(C)O');
    assert.strictEqual(r.size, 4);
});

// =========================================================================
// (9) Identical molecules → MCS is the full graph
// =========================================================================
test('mcs(aspirin, aspirin).size = 13 (all heavy atoms)', function() {
    var r = mcsSmi('CC(=O)Oc1ccccc1C(=O)O', 'CC(=O)Oc1ccccc1C(=O)O');
    assert.strictEqual(r.size, 13);
});

// =========================================================================
// (10) Single-atom edge cases
// =========================================================================
test('mcs(C, CC) = 1 (single C atom)', function() {
    assert.strictEqual(size('C', 'CC'), 1);
});

test('mcs(C, C) = 1 (degenerate single-atom identity)', function() {
    assert.strictEqual(size('C', 'C'), 1);
});

test('mcs(N, C) = 0 (different element, no overlap)', function() {
    assert.strictEqual(size('N', 'C'), 0);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
