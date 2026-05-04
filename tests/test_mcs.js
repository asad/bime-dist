/**
 * tests/test_mcs.js — Maximum Common Substructure (SMSDMCS) regressions.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SMSD MCS');
var test = runner.test;

console.log('SMSD MCS');

function mcs(a, b) {
    return SMSDMCS.findMCSFromSmiles(a, b);
}
function size(a, b) { return mcs(a, b).size; }

// --- Identity ---
test('mcs(c1ccccc1, c1ccccc1) = full benzene (6 atoms)', function() {
    assert.strictEqual(size('c1ccccc1', 'c1ccccc1'), 6);
});

test('mcs(CCO, CCO) = 3 atoms (full mol)', function() {
    assert.strictEqual(size('CCO', 'CCO'), 3);
});

// --- Disjoint atoms ---
test('mcs(C, N) = 0 atoms (different elements)', function() {
    assert.strictEqual(size('C', 'N'), 0);
});

// --- Substituted benzene ---
test('mcs(c1ccccc1, c1ccccc1Cl) = 6 atoms (the benzene ring)', function() {
    assert.strictEqual(size('c1ccccc1', 'c1ccccc1Cl'), 6);
});

// --- Symmetry ---
test('MCS is symmetric: |mcs(A,B)| = |mcs(B,A)| for CCO and CCN', function() {
    assert.strictEqual(size('CCO', 'CCN'), size('CCN', 'CCO'));
});

test('MCS is symmetric: |mcs(benzene,toluene)| = |mcs(toluene,benzene)|', function() {
    assert.strictEqual(size('c1ccccc1', 'Cc1ccccc1'), size('Cc1ccccc1', 'c1ccccc1'));
});

// --- Determinism ---
test('MCS is deterministic: 3 runs of mcs(aspirin, salicylic acid) all equal', function() {
    var a = 'CC(=O)Oc1ccccc1OC(C)=O';
    var b = 'OC(=O)c1ccccc1O';
    var s1 = size(a, b);
    var s2 = size(a, b);
    var s3 = size(a, b);
    assert.strictEqual(s1, s2);
    assert.strictEqual(s2, s3);
});

// --- Lower bound check via mcs_size (uses upper bounds + greedy probe) ---
test('mcs_size on identical graphs returns full atom count', function() {
    var bz = SmilesParser.parse('c1ccccc1');
    var g1 = new SMSDGraph.SMSDGraph(bz);
    var g2 = new SMSDGraph.SMSDGraph(SmilesParser.parse('c1ccccc1'));
    assert.strictEqual(SMSDMCS.mcs_size(g1, g2), 6);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
