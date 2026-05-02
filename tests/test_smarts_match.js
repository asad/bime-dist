/**
 * tests/test_smarts_match.js — SMARTS / VF2 substructure matching regressions.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SMARTS Match');
var test = runner.test;

console.log('SMARTS Match');

function P(smi) { return SmilesParser.parse(smi); }
function Q(sma) { return SmartsParser.parse(sma); }
function nMatches(targetSmi, smarts) {
    return SmartsMatch.match(P(targetSmi), Q(smarts)).length;
}

// --- Element / wildcard ---
test('[#6] on benzene → 6 matches', function() {
    assert.strictEqual(nMatches('c1ccccc1', '[#6]'), 6);
});

test('wildcard [*] on benzene → 6 matches', function() {
    assert.strictEqual(nMatches('c1ccccc1', '[*]'), 6);
});

// --- H count / connectivity ---
test('[OH] on phenol → 1 match (the hydroxyl O)', function() {
    assert.strictEqual(nMatches('c1ccc(O)cc1', '[OH]'), 1);
});

test('[NX3;H2] on ethylamine CCN → 1 match', function() {
    assert.strictEqual(nMatches('CCN', '[NX3;H2]'), 1);
});

// --- Aromaticity ---
test('benzene SMARTS c1ccccc1 on benzene → at least 1 match', function() {
    var n = nMatches('c1ccccc1', 'c1ccccc1');
    assert.ok(n >= 1, 'expected ≥1 mapping, got ' + n);
});

test('[a] on benzene → 6 matches (every C is aromatic)', function() {
    assert.strictEqual(nMatches('c1ccccc1', '[a]'), 6);
});

test('[A] on benzene → 0 matches (no aliphatic atom)', function() {
    assert.strictEqual(nMatches('c1ccccc1', '[A]'), 0);
});

// --- Ring membership ---
test('[R] on cyclohexane → 6 matches', function() {
    assert.strictEqual(nMatches('C1CCCCC1', '[R]'), 6);
});

test('[R] on ethane CC → 0 matches (no ring atoms)', function() {
    assert.strictEqual(nMatches('CC', '[R]'), 0);
});

// --- Recursive SMARTS ---
test('[$(*~[#7])] on pyridine → at least 1 match (C adjacent to N)', function() {
    var n = nMatches('c1ccncc1', '[$(*~[#7])]');
    assert.ok(n >= 1, 'expected ≥1 atom adjacent to N, got ' + n);
});

// --- More coverage ---
test('[+1] on [NH4+] → 1 match', function() {
    assert.strictEqual(nMatches('[NH4+]', '[+1]'), 1);
});

test('[r6] on naphthalene → 10 matches', function() {
    assert.strictEqual(nMatches('c1ccc2ccccc2c1', '[r6]'), 10);
});

test('[D2] on benzene → 6 matches (every C has degree 2)', function() {
    assert.strictEqual(nMatches('c1ccccc1', '[D2]'), 6);
});

test('hasMatch returns true for [#6] on benzene', function() {
    assert.strictEqual(SmartsMatch.hasMatch(P('c1ccccc1'), Q('[#6]')), true);
});

test('hasMatch returns false for N in benzene', function() {
    assert.strictEqual(SmartsMatch.hasMatch(P('c1ccccc1'), Q('[#7]')), false);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
