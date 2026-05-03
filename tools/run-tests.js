/**
 * tools/run-tests.js — BIME v1.0.2 regression test runner.
 *
 * Loads each tests/test_*.js in order, prints results, and exits 1 on any
 * failure. Plain Node, no external dependencies.
 *
 * Usage:  node tools/run-tests.js
 */
'use strict';

var path = require('path');
var fs = require('fs');

var TESTS_DIR = path.join(__dirname, '..', 'tests');
var FILES = [
    'test_smiles_parser.js',
    'test_smiles_writer.js',
    'test_smarts_match.js',
    'test_substructure_vf2.js',
    'test_mcs.js',
    'test_aromaticity.js',
    'test_hcount.js',
    'test_round_trip.js',
    'test_history.js',
    'test_chem_torture.js',
    'test_aam_round_trip.js',
    'test_charge_iso_round_trip.js',
    'test_rdt.js',
    'test_aam.js',
    'test_mcs_extended.js',
    'test_substructure_extended.js',
    'test_bondchange.js',
    'test_canon.js',
    'test_layout.js',
    'test_v1_2_0_features.js',
    'test_v1_4_0_features.js',
    'test_v1_4_1_features.js',
    'test_v1_4_2_features.js',
    'test_v1_4_3_features.js'
];

var startedAt = Date.now();
var totalPassed = 0;
var totalFailed = 0;
var perFile = [];

console.log('=== BIME v1.0.2 Regression Test Suite ===');
console.log('');

for (var i = 0; i < FILES.length; i++) {
    var f = FILES[i];
    var fp = path.join(TESTS_DIR, f);
    if (!fs.existsSync(fp)) {
        console.log('[skip] ' + f + ' (not found)');
        continue;
    }
    // Each test file exports a `summary` function that returns
    // { label, passed, failed, failures }.
    delete require.cache[require.resolve(fp)];
    var summaryFn = require(fp);
    var s = (typeof summaryFn === 'function') ? summaryFn() : summaryFn;
    perFile.push({ file: f, summary: s });
    totalPassed += s.passed;
    totalFailed += s.failed;
    console.log('  → ' + s.label + ': ' + s.passed + ' passed, ' + s.failed + ' failed');
    console.log('');
}

var elapsed = Date.now() - startedAt;

console.log('-------------------------------------------');
console.log(totalPassed + ' passed, ' + totalFailed + ' failed, total ' + elapsed + ' ms');

if (totalFailed > 0) {
    console.log('');
    console.log('Failures:');
    perFile.forEach(function(pf) {
        pf.summary.failures.forEach(function(fail) {
            console.log('  ' + pf.file + ' :: ' + fail.name);
            console.log('    ' + fail.reason);
        });
    });
    process.exit(1);
}

process.exit(0);
