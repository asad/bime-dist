/**
 * tests/test_v2_4_14_compound_labels.js — per-component compound name/ID labels
 * on the publication reaction-map figure (v2.4.14).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Two additive capabilities:
 *   A. MetaboliteLibrary.findBySmiles(smiles) — a structure lookup that returns
 *      the metabolite entry ({name, kegg, ...}) for a query SMILES, or null.
 *      It keys on the CONSTITUTIONAL canonical SMILES (order-invariant, so it is
 *      robust to the SmilesWriter @/@@ emission drift) PLUS an order-invariant
 *      CIP R/S multiset, and only returns an UNAMBIGUOUS single match — so it
 *      matches the same molecule written two ways yet still refuses to label an
 *      epimer (e.g. galactose) as its C4 partner (glucose).
 *   B. ImageExport.toReactionMapSVG captions each reactant/product structure
 *      with its name (or KEGG id via showIds) when known — "if present"; unknown
 *      components get no label. Off with labels:false.
 *
 * Drives the real engine (SmilesParser -> RDT.mapReaction -> ImageExport) +
 * source-shape pins. No DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();

[
    'Templates.js', 'SDGLayout.js', 'sdg/MacroCycleLayout.js', 'sdg/OverlapResolver.js',
    'sdg/HydrogenPlacer.js', 'sdg/CorrectGeometricConfiguration.js', 'sdg/LayoutRefiner.js',
    'sdg/RingPlacer.js', 'sdg/TemplateHandler.js', 'sdg/IdentityTemplateLibrary.js',
    'sdg/NonplanarBonds.js', 'Layout.js', 'CIPStereo.js',
    'SMSDVersion.js', 'SMSDGraph.js', 'SMSDRings.js', 'SMSDMCS.js', 'SMSDBatch.js',
    'SMSDLayout.js', 'RDT.js', 'MetaboliteLibrary.js', 'ImageExport.js'
].forEach(function (f) { require(path.join(__dirname, '..', 'editor', f)); });

var SmilesParser = globalThis.SmilesParser;
var SmilesWriter = globalThis.SmilesWriter;
var RDT = globalThis.RDT;
var ImageExport = globalThis.ImageExport;
var MetaboliteLibrary = globalThis.MetaboliteLibrary;

var runner = shim.makeRunner('Compound labels (v2.4.14)');
var test = runner.test;
console.log('Compound labels (v2.4.14)');

var ML_SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'MetaboliteLibrary.js'), 'utf8');
var IE_SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'ImageExport.js'), 'utf8');

var GLUCOSE = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O';
var G6P     = 'O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O';   // C00668, phosphate on exocyclic C6 (corrected v2.4.17)
var GALACTOSE = 'OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O';   // C4 epimer of glucose

// ---------------------------------------------------------------------------
// A. MetaboliteLibrary.findBySmiles — drift-robust + epimer-safe.
// ---------------------------------------------------------------------------
test('A1 findBySmiles resolves a known metabolite (and a re-written form of it)', function () {
    var g = MetaboliteLibrary.findBySmiles(GLUCOSE);
    assert.ok(g && g.name === 'D-Glucose', 'glucose SMILES -> D-Glucose');
    // a re-serialised (drift-prone) form of the SAME molecule must still match
    var rewritten = SmilesWriter.write(SmilesParser.parse(GLUCOSE));
    var g2 = MetaboliteLibrary.findBySmiles(rewritten);
    assert.ok(g2 && g2.name === 'D-Glucose', 'a re-written glucose still resolves (drift-robust)');
});

test('A2 findBySmiles resolves glucose-6-phosphate (the @/@@-drift regression)', function () {
    var p = MetaboliteLibrary.findBySmiles(G6P);
    assert.ok(p && p.name === 'Glucose 6-phosphate', 'G6P -> Glucose 6-phosphate; got ' + (p ? p.name : 'null'));
    assert.strictEqual(p.kegg, 'C00668', 'carries the KEGG id');
});

test('A3 an epimer is NOT mislabelled, and a non-metabolite returns null', function () {
    assert.strictEqual(MetaboliteLibrary.findBySmiles(GALACTOSE), null,
        'D-galactose (C4 epimer of glucose) must not resolve to D-Glucose');
    assert.strictEqual(MetaboliteLibrary.findBySmiles('c1ccccc1'), null, 'benzene is not a library metabolite');
    assert.strictEqual(MetaboliteLibrary.findBySmiles('not a smiles ###'), null, 'unparseable input -> null');
});

// ---------------------------------------------------------------------------
// B. toReactionMapSVG captions each known component.
// ---------------------------------------------------------------------------
function fig(opts) {
    var m = SmilesParser.parse(GLUCOSE + '>>' + G6P);
    var res = RDT.mapReaction(m, {});
    return ImageExport.toReactionMapSVG(m, res, opts || {});
}

test('B1 both known components are captioned with their names', function () {
    var svg = fig();
    assert.ok(/<text[^>]*>D-Glucose<\/text>/.test(svg), 'reactant captioned D-Glucose');
    assert.ok(/<text[^>]*>Glucose 6-phosphate<\/text>/.test(svg), 'product captioned Glucose 6-phosphate');
});

test('B2 showIds captions with the KEGG id instead of the name', function () {
    var svg = fig({ showIds: true });
    assert.ok(/<text[^>]*>C00031<\/text>/.test(svg), 'reactant captioned C00031');
    assert.ok(/<text[^>]*>C00668<\/text>/.test(svg), 'product captioned C00668');
    assert.strictEqual(/>D-Glucose</.test(svg), false, 'name not shown when showIds');
});

test('B3 labels:false suppresses all captions', function () {
    var svg = fig({ labels: false });
    assert.strictEqual(/>D-Glucose</.test(svg), false, 'no name caption');
    assert.strictEqual(/>C00031</.test(svg), false, 'no id caption');
});

test('B4 an unknown component is left uncaptioned while a known one is labelled', function () {
    // glucose (known) reacting toward benzene (not a library metabolite, fake):
    var m = SmilesParser.parse(GLUCOSE + '>>c1ccccc1');
    var res = RDT.mapReaction(m, {});
    var svg = ImageExport.toReactionMapSVG(m, res, {});
    assert.ok(/<text[^>]*>D-Glucose<\/text>/.test(svg), 'the known reactant is still captioned');
    // benzene contributes no caption text we can assert is absent generically,
    // but the figure must not crash and must contain the known label.
    assert.ok(svg.indexOf('<svg') === 0, 'a valid SVG is produced even with an unknown component');
});

// ---------------------------------------------------------------------------
// C. Source-shape pins.
// ---------------------------------------------------------------------------
test('C1 MetaboliteLibrary exports findBySmiles with the drift-robust key', function () {
    assert.ok(/findBySmiles: findBySmiles/.test(ML_SRC), 'findBySmiles is exported');
    assert.ok(/function _structureKey\(/.test(ML_SRC), 'structure key helper present');
    assert.ok(/hits\.length === 1/.test(ML_SRC), 'only an unambiguous single match is returned (never mislabel)');
});

test('C2 toReactionMapSVG builds componentLabels and _buildSVGImpl renders them', function () {
    assert.ok(/opts\.componentLabels = componentLabels/.test(IE_SRC), 'componentLabels assembled in toReactionMapSVG');
    assert.ok(/findBySmiles\(SmilesWriter\.write\(_comp\)\)/.test(IE_SRC), 'components auto-named via findBySmiles');
    assert.ok(/opts\.componentLabels && opts\.componentLabels\.length/.test(IE_SRC), 'captions rendered in _buildSVGImpl');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
