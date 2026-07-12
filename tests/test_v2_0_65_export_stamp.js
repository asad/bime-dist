/**
 * tests/test_v2_0_65_export_stamp.js — provenance + tamper-evidence
 * stamps on every BIME export (v2.0.65).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 */
'use strict';

var assert = require('assert');
var path = require('path');
var childProcess = require('child_process');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/ExportStamp.js');

var runner = shim.makeRunner('Export stamp (v2.0.65)');
var test = runner.test;

console.log('Export stamp (v2.0.65)');

var ROOT = path.join(__dirname, '..');
var BIME = path.join(ROOT, 'bin', 'bime');

function runCli(args) {
    var r = childProcess.spawnSync(BIME, args, { encoding: 'utf8' });
    return { stdout: r.stdout || '', stderr: r.stderr || '', code: r.status };
}

// ---------------------------------------------------------------------------
// A. stampRecord — pure data
// ---------------------------------------------------------------------------
test('A1 stampRecord on empty content emits a sha384 (or fallback) hex digest', function () {
    var rec = ExportStamp.stampRecord('', { version: '2.0.65', now: '2026-05-25T00:00:00Z' });
    assert.strictEqual(rec.project, 'BIME');
    assert.strictEqual(rec.version, '2.0.65');
    assert.strictEqual(rec.license, 'Apache-2.0');
    assert.strictEqual(rec.generatedAt, '2026-05-25T00:00:00Z');
    assert.ok(rec.algo === 'sha384' || rec.algo === 'fnv1a64',
        'algo is sha384 (Node) or fnv1a64 (browser fallback)');
    assert.ok(/^[0-9a-f]+$/i.test(rec.hash), 'hash is hex');
    assert.strictEqual(rec.url, 'https://github.com/asad/bime-dist');
});

test('A2 stampRecord is deterministic for the same content + opts', function () {
    var a = ExportStamp.stampRecord('hello world', { version: '2.0.65', now: '2026-05-25T00:00:00Z' });
    var b = ExportStamp.stampRecord('hello world', { version: '2.0.65', now: '2026-05-25T00:00:00Z' });
    assert.deepStrictEqual(a, b);
});

// ---------------------------------------------------------------------------
// B. stampSvg — SVG metadata + watermark
// ---------------------------------------------------------------------------
test('B1 stampSvg injects <metadata id="bime-stamp"> with BIME version + hash', function () {
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"><circle cx="50" cy="30" r="10"/></svg>';
    var stamped = ExportStamp.stampSvg(svg, { version: '2.0.65', now: '2026-05-25T00:00:00Z' });
    assert.ok(stamped.indexOf('<metadata id="bime-stamp">') !== -1);
    assert.ok(stamped.indexOf('<bime:version>2.0.65</bime:version>') !== -1);
    assert.ok(/<bime:fingerprint algo="[^"]+">[0-9a-f]+<\/bime:fingerprint>/.test(stamped));
    assert.ok(stamped.indexOf('Apache-2.0') !== -1);
    assert.ok(stamped.indexOf('https://github.com/asad/bime-dist') !== -1);
});

test('B2 stampSvg injects a bottom-right BIME wordmark + 6-vertex glyph at low opacity', function () {
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"></svg>';
    var stamped = ExportStamp.stampSvg(svg);
    assert.ok(stamped.indexOf('class="bime-watermark"') !== -1, 'watermark group present');
    assert.ok(stamped.indexOf('opacity="0.55"') !== -1, 'low opacity');
    assert.ok(stamped.indexOf('>BIME</text>') !== -1, 'BIME wordmark text');
});

test('B3 stampSvg is idempotent — restamping replaces the existing block', function () {
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"></svg>';
    var once = ExportStamp.stampSvg(svg, { version: '2.0.65' });
    var twice = ExportStamp.stampSvg(once, { version: '2.0.65' });
    // Both should still have exactly one stamp + one watermark.
    var metaCount = (twice.match(/<metadata id="bime-stamp">/g) || []).length;
    var wmCount = (twice.match(/class="bime-watermark"/g) || []).length;
    assert.strictEqual(metaCount, 1, 'exactly one stamp metadata block');
    assert.strictEqual(wmCount, 1, 'exactly one watermark group');
});

test('B4 stampSvg preserves a pre-existing <title> (user-authored)', function () {
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"><title>Aspirin</title><circle cx="50" cy="30" r="10"/></svg>';
    var stamped = ExportStamp.stampSvg(svg);
    assert.ok(stamped.indexOf('<title>Aspirin</title>') !== -1,
        'pre-existing <title> survives stamping');
});

test('B5 stampSvg with watermark:false skips the wordmark but keeps the metadata', function () {
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"></svg>';
    var stamped = ExportStamp.stampSvg(svg, { watermark: false });
    assert.ok(stamped.indexOf('<metadata id="bime-stamp">') !== -1, 'metadata still present');
    assert.strictEqual(stamped.indexOf('class="bime-watermark"'), -1, 'no watermark');
});

// ---------------------------------------------------------------------------
// C. stampText — MOL / SDF / RXN / SMILES / generic
// ---------------------------------------------------------------------------
test('C1 stampText on MOL prepends a BIME-STAMP-BEGIN/END comment block', function () {
    var mol = '  BIME    05252611:30  2D\n\n  3  2  0\n...';
    var stamped = ExportStamp.stampText('mol', mol);
    assert.ok(/^# BIME-STAMP-BEGIN/.test(stamped), 'stamp at the top');
    assert.ok(stamped.indexOf('# BIME-STAMP-END') !== -1);
    assert.ok(stamped.indexOf('# Version:') !== -1);
    assert.ok(stamped.indexOf('# Copyright:') !== -1);
    assert.ok(stamped.indexOf(mol) !== -1, 'original MOL content preserved');
});

test('C2 stampText on SMILES appends the stamp BELOW the SMILES (first line stays clean)', function () {
    var smi = 'c1ccccc1';
    var stamped = ExportStamp.stampText('smiles', smi);
    var firstLine = stamped.split('\n')[0];
    assert.strictEqual(firstLine, 'c1ccccc1', 'first line is the bare SMILES');
    assert.ok(stamped.indexOf('# BIME-STAMP-BEGIN') !== -1, 'stamp present below');
});

test('C3 stampText is idempotent — restamping replaces the block', function () {
    var rxn = '$RXN\n\n  BIME    1.0\n\n  1  1\n';
    var once  = ExportStamp.stampText('rxn', rxn);
    var twice = ExportStamp.stampText('rxn', once);
    var beginCount = (twice.match(/# BIME-STAMP-BEGIN/g) || []).length;
    var endCount   = (twice.match(/# BIME-STAMP-END/g)   || []).length;
    assert.strictEqual(beginCount, 1, 'exactly one BEGIN marker');
    assert.strictEqual(endCount,   1, 'exactly one END marker');
});

test('C4 stampText hash is computed over the UN-stamped content (verifiable by stripping)', function () {
    var content = 'CCO';
    var stamped = ExportStamp.stampText('smiles', content);
    // Strip the stamp.
    var stripped = stamped.replace(/# BIME-STAMP-BEGIN[\s\S]*?# BIME-STAMP-END\s*/, '');
    // Trailing newline tolerance.
    stripped = stripped.replace(/\n+$/, '');
    assert.strictEqual(stripped, content,
        'stripping the stamp recovers the original content');
});

// ---------------------------------------------------------------------------
// D. extractStamp round-trip
// ---------------------------------------------------------------------------
test('D1 extractStamp recovers project + version + hash from a stamped SVG', function () {
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"></svg>';
    var stamped = ExportStamp.stampSvg(svg, { version: '2.0.65' });
    var rec = ExportStamp.extractStamp(stamped);
    assert.ok(rec, 'stamp extracted');
    assert.strictEqual(rec.medium, 'svg');
    assert.strictEqual(rec.version, '2.0.65');
    assert.ok(/^[0-9a-f]+$/i.test(rec.hash));
});

test('D2 extractStamp recovers project + version + hash from a stamped MOL', function () {
    var stamped = ExportStamp.stampText('mol', 'mol content here');
    var rec = ExportStamp.extractStamp(stamped);
    assert.ok(rec, 'stamp extracted');
    assert.strictEqual(rec.medium, 'text');
    assert.ok(/^[0-9a-f]+$/i.test(rec.hash));
});

test('D3 extractStamp returns null on unstamped content', function () {
    assert.strictEqual(ExportStamp.extractStamp('<svg></svg>'), null);
    assert.strictEqual(ExportStamp.extractStamp('hello world'), null);
    assert.strictEqual(ExportStamp.extractStamp(''), null);
    assert.strictEqual(ExportStamp.extractStamp(null), null);
});

// ---------------------------------------------------------------------------
// E. CLI integration — default-on stamping + --no-stamp opt-out
// ---------------------------------------------------------------------------
test('E1 bime export --format svg emits a stamped SVG by default', function () {
    var r = runCli(['export', 'c1ccccc1', '--format', 'svg']);
    assert.strictEqual(r.code, 0);
    assert.ok(r.stdout.indexOf('<metadata id="bime-stamp">') !== -1,
        'svg output carries the BIME provenance metadata');
    assert.ok(r.stdout.indexOf('class="bime-watermark"') !== -1,
        'svg output carries the BIME watermark');
});

test('E2 bime export --format mol emits a stamped MOL by default', function () {
    var r = runCli(['export', 'Cc1ccccc1', '--format', 'mol']);
    assert.strictEqual(r.code, 0);
    assert.ok(/^# BIME-STAMP-BEGIN/.test(r.stdout),
        'mol output starts with a BIME stamp');
});

test('E3 bime export --no-stamp opt-out drops the stamp', function () {
    var r = runCli(['export', 'c1ccccc1', '--format', 'svg', '--no-stamp']);
    assert.strictEqual(r.code, 0);
    assert.strictEqual(r.stdout.indexOf('<metadata id="bime-stamp">'), -1,
        'svg has no stamp metadata under --no-stamp');
    assert.strictEqual(r.stdout.indexOf('class="bime-watermark"'), -1,
        'svg has no watermark under --no-stamp');
});

test('E4 bime export --format smiles --no-stamp returns a bare SMILES', function () {
    var r = runCli(['export', 'c1ccccc1', '--format', 'smiles', '--no-stamp']);
    assert.strictEqual(r.code, 0);
    var trimmed = r.stdout.trim();
    assert.strictEqual(trimmed.indexOf('# BIME-STAMP'), -1,
        'no stamp comments in bare-SMILES mode');
});

// ---------------------------------------------------------------------------
// F. version stamp
// ---------------------------------------------------------------------------
test('F1 ExportStamp.version exists', function () {
    assert.strictEqual(typeof ExportStamp.version, 'string');
    assert.ok(/^\d+\.\d+\.\d+/.test(ExportStamp.version));
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
