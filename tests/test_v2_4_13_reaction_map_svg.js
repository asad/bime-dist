/**
 * tests/test_v2_4_13_reaction_map_svg.js — publication-quality mapped-reaction
 * figure export (v2.4.13).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * ImageExport.toReactionMapSVG renders a MAPPED reaction as a static,
 * publication-quality figure — the batchable equivalent of the workbench's
 * on-screen reaction-mapping view. Four overlays, all reusing the real engine:
 *   1. COLOUR MOL  — per-atom soft halos tint mapped sub-fragments (the trace
 *      colouring) WITHOUT touching the CPK atom-label colours.
 *   2. ATOM-ATOM MAP NUMBERS — a reactant atom and its product partner share
 *      one number (the AAM read-out + the trace pairing).
 *   3. BOND CHANGES (BC) — changed bonds are recoloured: red cleaved, green
 *      formed, amber order change.
 *   4. TRACE — the sub-fragment halos carry the same tint across the arrow, so
 *      an atom's fate reads left-to-right.
 * Confidence is intentionally NOT drawn (per the feature brief).
 *
 * This suite drives the REAL engine (SmilesParser -> RDT.mapReaction -> Layout
 * -> ImageExport), plus source-shape pins on the overlay API + the CLI wiring.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();

// Layout + SDG deps and the SMSD/RDT mapping stack are not part of loadAll();
// require them explicitly (mirrors test_v2_4_7 + the AAM suites).
[
    'Templates.js', 'SDGLayout.js', 'sdg/MacroCycleLayout.js', 'sdg/OverlapResolver.js',
    'sdg/HydrogenPlacer.js', 'sdg/CorrectGeometricConfiguration.js', 'sdg/LayoutRefiner.js',
    'sdg/RingPlacer.js', 'sdg/TemplateHandler.js', 'sdg/IdentityTemplateLibrary.js',
    'sdg/NonplanarBonds.js', 'Layout.js',
    'SMSDVersion.js', 'SMSDGraph.js', 'SMSDRings.js', 'SMSDMCS.js', 'SMSDBatch.js',
    'SMSDLayout.js', 'RDT.js', 'ImageExport.js'
].forEach(function (f) { require(path.join(__dirname, '..', 'editor', f)); });

var SmilesParser = globalThis.SmilesParser;
var Layout = globalThis.Layout;
var RDT = globalThis.RDT;
var ImageExport = globalThis.ImageExport;

var runner = shim.makeRunner('Reaction-map SVG export (v2.4.13)');
var test = runner.test;
console.log('Reaction-map SVG export (v2.4.13)');

var IE_SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'ImageExport.js'), 'utf8');
var CLI_SRC = fs.readFileSync(path.join(__dirname, '..', 'tools', 'bime-cli.js'), 'utf8');

// Esterification: benzoic acid + ethanol -> ethyl benzoate + water. A clean
// mapped reaction with a ring, a cleaved bond and a formed bond.
var RXN = 'c1ccccc1C(=O)O.OCC>>c1ccccc1C(=O)OCC.O';

function mapAndRender(smi, opts) {
    var mol = SmilesParser.parse(smi);
    assert.ok(mol && mol.atoms.length, 'parsed ' + smi);
    var res = RDT.mapReaction(mol, {});
    // toReactionMapSVG owns the layout (refine + correctly place the arrow);
    // do NOT pre-run Layout.layout — it would scramble the reaction scheme.
    return { mol: mol, res: res, svg: ImageExport.toReactionMapSVG(mol, res, opts || {}) };
}

// ---------------------------------------------------------------------------
// A. The figure carries all four overlays.
// ---------------------------------------------------------------------------
test('A1 colour-mol trace halos are painted (pale soft circles behind atoms)', function () {
    var out = mapAndRender(RXN);
    assert.ok(out.res.status === 'mapped', 'the esterification maps');
    var halos = (out.svg.match(/<circle[^>]*opacity="0.85"/g) || []).length;
    assert.ok(halos >= 4, 'several sub-fragment trace halos are drawn; got ' + halos);
});

test('A2 atom-atom map numbers are drawn, and a reactant + its product partner share one', function () {
    var out = mapAndRender(RXN);
    assert.ok(/font-weight="bold" text-anchor="middle">[0-9]+<\/text>/.test(out.svg),
        'map-number labels are rendered');
    // toReactionMapSVG sets atom.mapNumber from the mapping; verify a pair shares it.
    var keys = Object.keys(out.res.mapping || {});
    assert.ok(keys.length > 0, 'the reaction produced a mapping');
    var rId = parseInt(keys[0], 10);
    var pId = out.res.mapping[keys[0]];
    var ra = out.mol.getAtom(rId), pa = out.mol.getAtom(pId);
    assert.ok(ra && pa, 'mapped pair atoms exist in the laid-out reaction');
    assert.ok(ra.mapNumber > 0 && ra.mapNumber === pa.mapNumber,
        'a reactant atom and its product partner share one map number (' +
            ra.mapNumber + ' vs ' + pa.mapNumber + ')');
});

test('A3 bond changes recolour the changed bonds (cleaved + formed)', function () {
    var out = mapAndRender(RXN);
    // publication preset => printMode => CMYK-safe BC variants.
    assert.ok(/stroke="#cc0000"/.test(out.svg), 'a cleaved bond is drawn red');
    assert.ok(/stroke="#008800"/.test(out.svg), 'a formed bond is drawn green');
});

test('A4 the reaction is laid out as a scheme (arrow present, components separated)', function () {
    var out = mapAndRender(RXN);
    assert.ok(out.mol.reactionArrow, 'the reaction has an arrow');
    assert.ok(out.mol.getComponents().length >= 3, 'multiple reactant/product components');
    // The figure is auto-sized to a tight aspect (wide scheme, not a tiny square).
    var vb = out.svg.match(/viewBox="0 0 ([0-9.]+) ([0-9.]+)"/);
    assert.ok(vb, 'svg has a viewBox');
    assert.ok(parseFloat(vb[1]) > parseFloat(vb[2]), 'a reaction scheme is wider than it is tall');
});

test('A5 publication quality by default: Helvetica labels + white background', function () {
    var out = mapAndRender(RXN);
    assert.ok(/Helvetica/.test(out.svg), 'publication font stack');
    assert.ok(/<rect width="[0-9.]+" height="[0-9.]+" fill="#ffffff"\/>/.test(out.svg),
        'opaque white figure background by default');
});

// ---------------------------------------------------------------------------
// B. Confidence is intentionally omitted (per the brief).
// ---------------------------------------------------------------------------
test('B1 NO confidence caption is drawn in the figure', function () {
    var out = mapAndRender(RXN);
    assert.strictEqual(/[Cc]onfidence/.test(out.svg), false,
        'the reaction-map figure must not render a confidence caption');
});

// ---------------------------------------------------------------------------
// C. Overlays are additive — the single-molecule export path is unchanged.
// ---------------------------------------------------------------------------
test('C1 a plain single-molecule toSVG is overlay-free (no halos / BC / no regression)', function () {
    var m = SmilesParser.parse('c1ccccc1C(=O)O'); Layout.layout(m);
    var plain = ImageExport.toSVG(m, {});
    assert.strictEqual(/<circle[^>]*opacity="0.85"/.test(plain), false, 'no trace halos on a plain export');
    assert.strictEqual(/stroke="#cc0000"|stroke="#008800"/.test(plain), false, 'no bond-change recolour on a plain export');
});

test('C2 caller can opt into a transparent figure and the reaction-centre rings', function () {
    var out = mapAndRender(RXN, { background: 'transparent', showReactionCenter: true });
    assert.strictEqual(/fill="#ffffff"\/>/.test(out.svg), false, 'transparent: no white background rect');
    assert.ok(/stroke-dasharray="3,2"/.test(out.svg), 'reaction-centre dashed rings when opted in');
});

// ---------------------------------------------------------------------------
// D. Source-shape pins — the overlay API + palettes + CLI wiring.
// ---------------------------------------------------------------------------
test('D1 ImageExport exposes toReactionMapSVG + the overlay palettes', function () {
    assert.ok(IE_SRC.indexOf('ImageExport.toReactionMapSVG = function') !== -1, 'toReactionMapSVG entry point');
    assert.ok(/var COMPONENT_PAIR_PALETTE = \[/.test(IE_SRC), 'sub-fragment trace palette present');
    assert.ok(/var BOND_CHANGE_COLORS = \{/.test(IE_SRC), 'bond-change colour map present');
    assert.ok(/function _bcColorKey\(/.test(IE_SRC), 'bondChange.type -> colour normaliser present');
    assert.ok(/opts\.componentPairs/.test(IE_SRC) && /opts\.bondChanges/.test(IE_SRC),
        '_buildSVGImpl reads the componentPairs + bondChanges overlay opts');
});

test('D2 the CLI wires `aam --format svg` through toReactionMapSVG', function () {
    assert.ok(/format === 'svg'/.test(CLI_SRC), 'cmdAam has an svg format branch');
    assert.ok(/ImageExport\.toReactionMapSVG\(rxn, res/.test(CLI_SRC),
        'the svg branch renders via toReactionMapSVG');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
