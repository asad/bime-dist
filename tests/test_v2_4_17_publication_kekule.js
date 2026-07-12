/**
 * tests/test_v2_4_17_publication_kekule.js — the circles-off / publication export
 * must draw Kekulé double bonds, not a bare (saturated-looking) ring (v2.4.17).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background: aromatic ring bonds are stored single. The PUBLICATION preset sets
 * showAromaticCircles:false but nothing kekulised, so benzene exported as a bare
 * hexagon (6 single lines, no circle) — chemically it read as cyclohexane. The
 * prior test only asserted `circles === 0`, so it passed over the false picture.
 *
 * Fix: a perfect-matching (backtracking) kekuliser assigns alternating double bonds
 * over the aromatic subgraph when circles are off; those bonds render double.
 *
 * These tests assert that the publication export shows unsaturation (extra lines
 * for the double bonds) and no aromatic circle, while the default export keeps
 * the circle. Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'CIPStereo.js'));
require(path.join(__dirname, '..', 'editor', 'Layout.js'));
require(path.join(__dirname, '..', 'editor', 'ImageExport.js'));

var SmilesParser = globalThis.SmilesParser;
var Layout = globalThis.Layout;
var ImageExport = globalThis.ImageExport;

var runner = shim.makeRunner('publication Kekulé export (v2.4.17)');
var test = runner.test;
console.log('publication Kekulé export (v2.4.17)');

function laidOut(smi) { var m = SmilesParser.parse(smi); Layout.layout(m); return m; }
function pubSVG(m) {
    return ImageExport.toPublicationSVG ? ImageExport.toPublicationSVG(m)
        : ImageExport.toSVG(m, { showAromaticCircles: false, printMode: true });
}
function count(svg, re) { return (svg.match(re) || []).length; }

// molecule -> minimum number of <line> elements once kekulised (ring bonds +
// one extra line per Kekulé double bond). Lower bounds, robust to substituents.
var CASES = { 'benzene': ['c1ccccc1', 9], 'pyridine': ['c1ccncc1', 9], 'naphthalene': ['c1ccc2ccccc2c1', 16], 'pyrrole': ['c1cc[nH]c1', 7] };

test('A publication export draws Kekulé double bonds (no bare aromatic ring)', function () {
    Object.keys(CASES).forEach(function (name) {
        var svg = pubSVG(laidOut(CASES[name][0]));
        assert.strictEqual(count(svg, /<circle /g), 0, name + ': publication has no aromatic circle');
        assert.ok(count(svg, /<line /g) >= CASES[name][1],
            name + ': publication must draw double bonds (>= ' + CASES[name][1] + ' lines, got ' + count(svg, /<line /g) + ')');
    });
});

test('B benzene: publication (Kekulé) has more lines than a bare hexagon would', function () {
    var svg = pubSVG(laidOut('c1ccccc1'));
    // 6 ring bonds + 3 Kekulé doubles = 9; a bare hexagon would be exactly 6.
    assert.ok(count(svg, /<line /g) >= 9, 'benzene publication shows 3 double bonds');
});

test('C pyrrole gets exactly 2 ring double bonds (N is a lone-pair donor, not =N)', function () {
    var svg = pubSVG(laidOut('c1cc[nH]c1'));
    // 5 ring bonds + 2 doubles = 7
    assert.strictEqual(count(svg, /<line /g), 7, 'pyrrole = 5 ring bonds + 2 doubles');
});

test('D default export still uses the aromatic circle (unchanged)', function () {
    var svg = ImageExport.toSVG(laidOut('c1ccccc1'), {});
    assert.strictEqual(count(svg, /<circle /g), 1, 'default benzene keeps its aromatic circle');
    assert.strictEqual(count(svg, /<line /g), 6, 'default benzene draws 6 single ring bonds (circle carries aromaticity)');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
