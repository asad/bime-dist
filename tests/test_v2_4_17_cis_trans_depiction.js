/**
 * tests/test_v2_4_17_cis_trans_depiction.js — the 2D layout must draw cis/trans
 * (Z/E) double bonds with the correct geometry (v2.4.17).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background: Layout.layout wipes all coordinates and rebuilds a trans zigzag,
 * and CorrectGeometricConfiguration.correct() only fired when bond.cipLabel was
 * already set — which the layout never did (assignEZ reads coordinates, so it is
 * circular after a re-layout). Result: EVERY C=C was drawn trans — maleate was
 * depicted identically to fumarate, cis-2-butene like trans, oleic like elaidic.
 * A student would be shown the wrong geometric isomer.
 *
 * Fix: before the geometry-correction pass, populate bond.cipLabel from the
 * geometry-free directional markers (SMILES / \) via CIPStereo.doubleBondEZ, so
 * the corrector reflects cis bonds to the correct side.
 *
 * This test lays out real molecules and MEASURES the drawn side of the two
 * priority substituents across the double-bond axis. Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'CIPStereo.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'CorrectGeometricConfiguration.js'));
require(path.join(__dirname, '..', 'editor', 'Layout.js'));

var SmilesParser = globalThis.SmilesParser;
var Layout = globalThis.Layout;
var CIPStereo = globalThis.CIPStereo;

var runner = shim.makeRunner('cis/trans depiction (v2.4.17)');
var test = runner.test;
console.log('cis/trans depiction (v2.4.17)');

// Return 'cis' | 'trans' | null for the first C=C in the laid-out molecule,
// measured from the drawn coordinates of the CIP-priority substituent on each end.
function drawnGeometry(smi) {
    var m = SmilesParser.parse(smi);
    Layout.layout(m);
    var db = null;
    for (var i = 0; i < m.bonds.length; i++) {
        var b = m.bonds[i];
        var a1 = m.getAtom(b.atom1), a2 = m.getAtom(b.atom2);
        if (b.type === 2 && a1 && a2 && a1.symbol === 'C' && a2.symbol === 'C') { db = b; break; }
    }
    if (!db) return null;
    var e1 = m.getAtom(db.atom1), e2 = m.getAtom(db.atom2);
    function pri(aid, exc) {
        var es = CIPStereo.cipPriorities(m, aid) || [];
        for (var k = 0; k < es.length; k++) { if (es[k].neighborId !== exc && es[k].neighborId > 0) return es[k].neighborId; }
        return null;
    }
    var p1 = m.getAtom(pri(e1.id, e2.id)), p2 = m.getAtom(pri(e2.id, e1.id));
    if (!p1 || !p2) return null;
    var bdx = e2.x - e1.x, bdy = e2.y - e1.y, bl = Math.hypot(bdx, bdy);
    if (bl < 1e-9) return null;
    var nx = -bdy / bl, ny = bdx / bl;
    var s1 = (p1.x - e1.x) * nx + (p1.y - e1.y) * ny;
    var s2 = (p2.x - e2.x) * nx + (p2.y - e2.y) * ny;
    return (s1 * s2 > 0) ? 'cis' : 'trans';
}

test('A cis (Z) double bonds are drawn cis, trans (E) drawn trans', function () {
    var Z = { 'maleate': 'OC(=O)/C=C\\C(=O)O', 'cis-2-butene': 'C/C=C\\C', 'oleate': 'CCCCCCCC/C=C\\CCCCCCCC(=O)O' };
    var E = { 'fumarate': 'OC(=O)/C=C/C(=O)O', 'trans-2-butene': 'C/C=C/C', 'elaidate': 'CCCCCCCC/C=C/CCCCCCCC(=O)O' };
    Object.keys(Z).forEach(function (n) { assert.strictEqual(drawnGeometry(Z[n]), 'cis', n + ' must be drawn cis'); });
    Object.keys(E).forEach(function (n) { assert.strictEqual(drawnGeometry(E[n]), 'trans', n + ' must be drawn trans'); });
});

test('B maleate and fumarate are drawn DIFFERENTLY (the headline bug)', function () {
    assert.notStrictEqual(drawnGeometry('OC(=O)/C=C\\C(=O)O'), drawnGeometry('OC(=O)/C=C/C(=O)O'),
        'cis-maleate and trans-fumarate must not depict identically');
});

test('C an unspecified double bond still lays out (no marker => no forced flip)', function () {
    // plain butene, no directional markers — must not throw and must produce a C=C
    assert.ok(drawnGeometry('CC=CC') !== null, 'unmarked C=C still lays out');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
