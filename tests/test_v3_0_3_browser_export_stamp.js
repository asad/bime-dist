/**
 * tests/test_v3_0_3_browser_export_stamp.js — browser SVG exports carry the
 * BIME provenance stamp (v3.0.3).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * The CLI (`bime export`) already runs its output through ExportStamp. The
 * browser publication-SVG download (js/workbench.js runExportPreset) did NOT, so
 * a figure exported from the workbench shipped without the bime:version metadata,
 * SHA fingerprint, Apache-2.0 attribution, and wordmark that the CLI export
 * carries. v3.0.3 stamps the browser download at the WORKBENCH layer
 * (stampExportSvg → ExportStamp.stampSvg), deliberately NOT inside ImageExport —
 * so ImageExport, the on-screen preview, and the regression tests stay unstamped
 * and the CLI is unchanged.
 *
 * This test pins (1) the source-shape wiring and (2) the separation: the raw
 * ImageExport SVG is unstamped, ExportStamp adds the stamp, and it is idempotent.
 *
 * Real engine, no DOM (the download itself needs a browser; the stamp logic does not).
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'ImageExport.js'));
require(path.join(__dirname, '..', 'editor', 'ExportStamp.js'));

var SmilesParser = globalThis.SmilesParser;
var Layout = globalThis.Layout;
var ImageExport = globalThis.ImageExport;
var ExportStamp = globalThis.ExportStamp;

var runner = shim.makeRunner('Browser export stamp (v3.0.3)');
var test = runner.test;
console.log('Browser export stamp (v3.0.3)');

var WB = fs.readFileSync(path.join(__dirname, '..', 'js', 'workbench.js'), 'utf8');

// ---------------------------------------------------------------------------
// A. source-shape: the browser SVG download is stamped at the workbench layer.
// ---------------------------------------------------------------------------
test('A1 stampExportSvg is defined and delegates to ExportStamp.stampSvg', function () {
    assert.ok(/function\s+stampExportSvg\s*\(/.test(WB), 'stampExportSvg() defined in js/workbench.js');
    var body = WB.slice(WB.indexOf('function stampExportSvg'));
    body = body.slice(0, body.indexOf('function runExportPreset'));
    assert.ok(body.indexOf('ExportStamp.stampSvg') !== -1,
        'stampExportSvg delegates to ExportStamp.stampSvg (workbench-layer stamp)');
});

test('A2 the publication-SVG download stamps before ImageExport.downloadSVG', function () {
    assert.ok(WB.indexOf('ImageExport.downloadSVG(stampExportSvg(svg)') !== -1,
        'runExportPreset pub-svg wraps the SVG in stampExportSvg() before download');
});

test('A3 ImageExport is NOT modified to stamp (CLI + preview stay unstamped)', function () {
    var IE = fs.readFileSync(path.join(__dirname, '..', 'editor', 'ImageExport.js'), 'utf8');
    assert.strictEqual(IE.indexOf('ExportStamp'), -1,
        'ImageExport.js must not reference ExportStamp (that would restamp the CLI)');
});

// ---------------------------------------------------------------------------
// B. functional: ImageExport is unstamped, ExportStamp adds provenance.
// ---------------------------------------------------------------------------
function benzeneSvg() {
    var m = SmilesParser.parse('c1ccccc1');
    if (Layout && Layout.layout) { Layout.layout(m); }
    return ImageExport.toSVG(m, { width: 320, height: 240, padding: 16, background: 'transparent' });
}

test('B1 raw ImageExport SVG carries no bime provenance stamp', function () {
    var svg = benzeneSvg();
    assert.ok(svg && svg.indexOf('<svg') !== -1, 'ImageExport produced an SVG');
    assert.strictEqual(svg.indexOf('id="bime-stamp"'), -1, 'the raw ImageExport SVG is unstamped');
});

test('B2 ExportStamp.stampSvg adds the provenance stamp (version metadata)', function () {
    var canonical = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'versions.json'), 'utf8')).bime;
    var stamped = ExportStamp.stampSvg(benzeneSvg());
    assert.ok(stamped.indexOf('id="bime-stamp"') !== -1, 'stamped SVG carries the bime-stamp metadata block');
    assert.ok(stamped.indexOf(canonical) !== -1, 'stamped SVG mentions the current BIME version');
});

test('B3 the stamp is idempotent (re-export never double-stamps)', function () {
    var once = ExportStamp.stampSvg(benzeneSvg());
    var twice = ExportStamp.stampSvg(once);
    var count = twice.split('id="bime-stamp"').length - 1;
    assert.strictEqual(count, 1, 're-stamping replaces, not appends, the stamp block');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
