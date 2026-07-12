/**
 * tests/test_v2_4_7_publication_export.js — publication-quality export (v2.4.7):
 * a journal-standard PUBLICATION preset + a transparent-background option.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Two additive capabilities on the mature ImageExport engine:
 *   1. TRANSPARENT BACKGROUND. Passing background:'transparent' now produces a
 *      truly transparent SVG: the full-canvas background <rect> is OMITTED, and
 *      the atom-label "knockout" halos (which paint the background colour to mask
 *      bonds behind a label) switch to no-fill (the bold glyph drawn on top stays
 *      legible). toPNG already left the canvas unfilled for transparent.
 *   2. PUBLICATION PRESET (journal-standard). opts.publication / toPublicationSVG
 *      folds in journal-standard Helvetica/Arial labels + CMYK-safe element
 *      colours (via printMode), Kekule double bonds (showAromaticCircles:false),
 *      heavier bonds + larger labels, and no watermark. The caller's explicit
 *      options still win over the preset (e.g. publication + transparent).
 *
 * This suite drives the REAL engine (SmilesParser -> Layout -> ImageExport, the
 * same stack test_layout.js loads), plus a source-shape pin on ImageExport.js so
 * the preset definition + the transparent/merge wiring cannot silently regress.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();

// Layout + its SDG dependencies are not part of shim.loadAll(); require them
// explicitly (mirrors test_layout.js) so ImageExport gets laid-out coordinates.
[
    'Templates.js', 'SDGLayout.js', 'sdg/MacroCycleLayout.js', 'sdg/OverlapResolver.js',
    'sdg/HydrogenPlacer.js', 'sdg/CorrectGeometricConfiguration.js', 'sdg/LayoutRefiner.js',
    'sdg/RingPlacer.js', 'sdg/TemplateHandler.js', 'sdg/IdentityTemplateLibrary.js',
    'sdg/NonplanarBonds.js', 'Layout.js', 'ImageExport.js'
].forEach(function (f) { require(path.join(__dirname, '..', 'editor', f)); });

var SmilesParser = globalThis.SmilesParser;
var Layout = globalThis.Layout;
var ImageExport = globalThis.ImageExport;

var runner = shim.makeRunner('Publication-quality export (v2.4.7)');
var test = runner.test;
console.log('Publication-quality export (v2.4.7)');

var SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'ImageExport.js'), 'utf8');

function P(smi) {
    var m = SmilesParser.parse(smi);
    assert.ok(m && m.atoms && m.atoms.length, 'parsed ' + smi);
    Layout.layout(m);
    return m;
}
// Full-canvas background rect: <rect width="W" height="H" fill="..."/>
var BG_RECT = /<rect width="[0-9.]+" height="[0-9.]+" fill="[^"]+"\/>/;
function circles(svg) { return (svg.match(/<circle/g) || []).length; }
function whiteFills(svg) { return (svg.match(/fill="#ffffff"/g) || []).length; }

var PHENOL = P('Oc1ccccc1');   // has an O label -> exercises the label halo
var BENZENE = P('c1ccccc1');   // aromatic ring -> exercises Kekule vs circles

// ---------------------------------------------------------------------------
// A. Transparent background.
// ---------------------------------------------------------------------------
test('A1 background:transparent OMITS the full-canvas background rect', function () {
    var white = ImageExport.toSVG(PHENOL, { background: '#ffffff' });
    var clear = ImageExport.toSVG(PHENOL, { background: 'transparent' });
    assert.ok(BG_RECT.test(white), 'an opaque background still paints the full-canvas rect');
    assert.strictEqual(BG_RECT.test(clear), false, 'transparent drops the full-canvas background rect');
});

test('A2 transparent leaves NO opaque white fills (background rect + label halos all clear)', function () {
    var white = ImageExport.toSVG(PHENOL, { background: '#ffffff' });
    var clear = ImageExport.toSVG(PHENOL, { background: 'transparent' });
    assert.ok(whiteFills(white) >= 1, 'opaque export paints white background + label knockout boxes');
    assert.strictEqual(whiteFills(clear), 0,
        'transparent export paints zero #ffffff fills (halos switch to no-fill)');
    // Still a real depiction: the molecule (bonds/atoms) is present.
    assert.ok(clear.indexOf('<text') !== -1 || clear.indexOf('<line') !== -1 || clear.indexOf('<path') !== -1,
        'the transparent SVG still contains the drawn structure');
});

// ---------------------------------------------------------------------------
// B. Publication preset (journal-standard).
// ---------------------------------------------------------------------------
test('B1 publication labels use the journal Helvetica/Arial font; the default does not', function () {
    var screen = ImageExport.toSVG(PHENOL, {});
    var pub = ImageExport.toPublicationSVG(PHENOL, {});
    assert.strictEqual(/Helvetica/.test((screen.match(/font-family="[^"]+"/) || [''])[0]), false,
        'the default screen export keeps the UI sans-serif stack (no Helvetica)');
    assert.ok(/Helvetica/.test(pub), 'the publication export uses the journal Helvetica/Arial stack');
});

test('B2 publication draws Kekule double bonds (no aromatic circles)', function () {
    var aromOn = ImageExport.toSVG(BENZENE, { showAromaticCircles: true });
    var aromOff = ImageExport.toSVG(BENZENE, { showAromaticCircles: false });
    var pub = ImageExport.toPublicationSVG(BENZENE, {});
    assert.ok(circles(aromOn) > circles(aromOff), 'the aromatic-circle option actually changes the depiction');
    assert.strictEqual(circles(aromOff), 0, 'benzene with circles off has no aromatic ring circle');
    assert.strictEqual(circles(pub), circles(aromOff),
        'publication matches circles-off (Kekule): no aromatic circle');
});

test('B3 publication uses heavier bonds than the screen default', function () {
    var screen = ImageExport.toSVG(PHENOL, {});
    var pub = ImageExport.toPublicationSVG(PHENOL, {});
    assert.ok(screen.indexOf('stroke-width="1.8"') !== -1, 'screen default bond width is 1.8');
    assert.ok(pub.indexOf('stroke-width="2.2"') !== -1, 'publication bond width is the heavier 2.2');
});

// ---------------------------------------------------------------------------
// C. Merge order: caller options win over the preset, which wins over DEFAULTS.
// ---------------------------------------------------------------------------
test('C1 publication composes with transparent (Helvetica AND no background rect)', function () {
    var svg = ImageExport.toSVG(PHENOL, { publication: true, background: 'transparent' });
    assert.ok(/Helvetica/.test(svg), 'the publication preset still applies the journal font');
    assert.strictEqual(BG_RECT.test(svg), false, "the caller's transparent background still wins");
});

test('C2 an explicit caller option overrides the preset value', function () {
    var svg = ImageExport.toSVG(PHENOL, { publication: true, bondWidth: 9 });
    assert.ok(svg.indexOf('stroke-width="9"') !== -1,
        'caller bondWidth:9 overrides the preset 2.2 (caller > preset > DEFAULTS)');
});

// ---------------------------------------------------------------------------
// D. No regression: the default export path is byte-for-byte unchanged in shape.
// ---------------------------------------------------------------------------
test('D1 the default (no-options) export is unchanged: white rect, screen font, 1.8 bonds, aromatic circle', function () {
    var screen = ImageExport.toSVG(PHENOL, {});
    var benz = ImageExport.toSVG(BENZENE, {});
    assert.ok(BG_RECT.test(screen), 'default still paints a background rect');
    assert.ok(whiteFills(screen) >= 1, 'default background is white');
    assert.strictEqual(/Helvetica/.test((screen.match(/font-family="[^"]+"/) || [''])[0]), false,
        'default keeps the screen font (publication is opt-in)');
    assert.ok(screen.indexOf('stroke-width="1.8"') !== -1, 'default bond width unchanged (1.8)');
    assert.ok(circles(benz) >= 1, 'default benzene still shows the aromatic circle (opt-in Kekule only)');
});

// ---------------------------------------------------------------------------
// E. Source-shape pins on editor/ImageExport.js.
// ---------------------------------------------------------------------------
test('E1 a PUBLICATION preset is defined with printMode + Kekule + no watermark', function () {
    var i = SRC.indexOf('var PUBLICATION = {');
    assert.ok(i !== -1, 'PUBLICATION preset object exists');
    var block = SRC.slice(i, SRC.indexOf('};', i) + 2);
    assert.ok(/printMode:\s*true/.test(block), 'PUBLICATION sets printMode:true (journal font + CMYK colours)');
    assert.ok(/showAromaticCircles:\s*false/.test(block), 'PUBLICATION sets showAromaticCircles:false (Kekule)');
    assert.ok(/watermark:\s*false/.test(block), 'PUBLICATION sets watermark:false (clean figure)');
    // The preset MUST NOT hardcode a background — the caller picks white/transparent.
    assert.strictEqual(/background:/.test(block), false, 'PUBLICATION leaves background to the caller');
});

test('E2 toSVG and toPNG both resolve options through _resolveOpts (preset folds in)', function () {
    assert.ok(/function _resolveOpts\(options\)/.test(SRC), '_resolveOpts helper exists');
    assert.ok(/toSVG = function\(mol, options\) \{\s*return _buildSVG\(mol, _resolveOpts\(options\)\)/.test(SRC),
        'toSVG routes through _resolveOpts');
    assert.ok(/toPNG = function\(mol, options\) \{\s*var opts = _resolveOpts\(options\)/.test(SRC),
        'toPNG routes through _resolveOpts');
    assert.ok(SRC.indexOf('ImageExport.toPublicationSVG = function') !== -1,
        'a toPublicationSVG convenience entry point exists');
});

test('E3 _buildSVGImpl omits the bg rect and neutralises label halos for transparent', function () {
    assert.ok(/var haloFill = \(background === 'transparent'\) \? 'none' : background;/.test(SRC),
        'haloFill is no-fill when the background is transparent');
    assert.ok(/var bgRect = \(background === 'transparent'\)\s*\?\s*''/.test(SRC),
        'the full-canvas bgRect is omitted (empty string) for transparent');
    // The halos actually consume haloFill (not the raw background) at both sites.
    assert.ok((SRC.match(/fill="' \+ haloFill \+ '"/g) || []).length >= 2,
        'both label-knockout rects fill with haloFill');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
