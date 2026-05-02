/**
 * tests/run-against-bundle.js - run the BIME regression suite against the
 * built bundle in dist/bime.js (and optionally dist/bime.min.js) instead
 * of the editor/*.js source files.
 *
 * The shim normally loads each editor module separately. This runner
 * monkey-patches Module._resolveFilename so that any require('../editor/<X>.js')
 * (or require_editor('X')) resolves to dist/bime.js, which on Node populates
 * globalThis with the same global symbols (Molecule, SmilesParser, ...) by
 * way of the bottom shim emitted by tools/build.js.
 *
 * Usage:    node tests/run-against-bundle.js [--min]
 * Exits 0 on full pass, 1 on any failure.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var path = require('path');
var fs = require('fs');
var Module = require('module');

var ROOT = path.join(__dirname, '..');
var DIST_DIR = path.join(ROOT, 'dist');
var EDITOR_DIR = path.join(ROOT, 'editor');

var useMin = process.argv.indexOf('--min') !== -1;
var bundlePath = path.join(DIST_DIR, useMin ? 'bime.min.js' : 'bime.js');

if (!fs.existsSync(bundlePath)) {
    console.error('bundle not built: ' + bundlePath);
    console.error('run `node tools/build.js` first.');
    process.exit(1);
}

// Build a set of the editor files the shim and tests require directly.
var editorFiles = {};
fs.readdirSync(EDITOR_DIR).forEach(function(f) {
    if (/\.js$/.test(f)) {
        editorFiles[path.join(EDITOR_DIR, f)] = true;
    }
});

// Monkey-patch resolution: any request that resolves to editor/*.js gets
// redirected to the bundle. The bundle is loaded once; subsequent requires
// are no-ops because Node caches it.
var origResolve = Module._resolveFilename;
Module._resolveFilename = function(request, parent, isMain, options) {
    var resolved;
    try {
        resolved = origResolve.call(this, request, parent, isMain, options);
    } catch (e) {
        throw e;
    }
    if (editorFiles[resolved]) {
        return bundlePath;
    }
    return resolved;
};

// Load the shim first (sets up window/document stubs and globalThis.window).
// Then load the bundle once so its global side effects (window.Molecule = ..,
// window.SmilesParser = .., ...) populate globalThis. Subsequent
// require_editor calls from tests/shim.js are short-circuited via the
// monkey-patched _resolveFilename above.
require('./shim.js');
require(bundlePath);

// Run the standard regression suite. tools/run-tests.js shells through the
// shim, which now sees the bundle in place of editor/*.js.
var TESTS_DIR = path.join(ROOT, 'tests');
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
    'test_layout.js'
];

var startedAt = Date.now();
var totalPassed = 0;
var totalFailed = 0;
var perFile = [];

console.log('=== BIME v1.1.2 Regression Suite (against ' + path.relative(ROOT, bundlePath) + ') ===');
console.log('');

for (var i = 0; i < FILES.length; i++) {
    var f = FILES[i];
    var fp = path.join(TESTS_DIR, f);
    if (!fs.existsSync(fp)) {
        console.log('[skip] ' + f + ' (not found)');
        continue;
    }
    delete require.cache[require.resolve(fp)];
    var summaryFn = require(fp);
    var s = (typeof summaryFn === 'function') ? summaryFn() : summaryFn;
    perFile.push({ file: f, summary: s });
    totalPassed += s.passed;
    totalFailed += s.failed;
    console.log('  -> ' + s.label + ': ' + s.passed + ' passed, ' + s.failed + ' failed');
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
