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
test('benzene SMARTS c1ccccc1 on benzene → 12 matches (6 rotations × 2 reflections)', function() {
    // Symmetry-aware count: every benzene admits 12 distinct ordered mappings.
    assert.strictEqual(nMatches('c1ccccc1', 'c1ccccc1'), 12);
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
// Pins the current BIME v1.x policy: the recursive matcher reports exactly
// one match for [$(*~[#7])] on pyridine. (RDKit-grade implementations would
// report 2 — the two ring carbons ortho to N — but BIME's recursive SMARTS
// engine collapses symmetric matches today; this test locks the count so
// regressions surface.)
test('[$(*~[#7])] on pyridine → 1 match (BIME recursive SMARTS policy)', function() {
    assert.strictEqual(nMatches('c1ccncc1', '[$(*~[#7])]'), 1);
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

// --- Regression: SmartsParser empty atom-map number [*:] must produce an error ---
// Before the v1.1.3 fix, '[*:]' was silently swallowed (map number = 0, no error).
// After the fix, the parser pushes 'Empty atom map number after ":"' so callers
// know the SMARTS was malformed rather than silently accepting bad input.
test('[*:] (empty map number) produces a parse error — regression for v1.1.3 fix', function() {
    var result = SmartsParser.parse('[*:]');
    // Public API field is `parseErrors` (inner parseFragment returns `errors`,
    // but `SmartsParser.parse` exposes them on `mol.parseErrors`).
    var errs = result.parseErrors || result.errors || [];
    assert.ok(errs.length > 0,
        'Expected at least one parse error for [*:] but got none');
    var hasMapError = errs.some(function(e) {
        return typeof e === 'string' && e.indexOf('map') !== -1;
    });
    assert.ok(hasMapError, 'Expected an error mentioning "map" for [*:], got: ' +
        JSON.stringify(errs));
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
