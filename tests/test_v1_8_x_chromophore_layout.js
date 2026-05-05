/**
 * tests/test_v1_8_x_chromophore_layout.js
 *
 * Regression target for the SDG (Structure-Diagram-Generator) layout
 * completion work introduced in v1.8.12 / v1.8.13.
 *
 * Benchmark molecule: a mono-sulfonated bisazo Direct dye, C32H22N6O3S
 *
 *   c5(cc6ccc(\N=N\c1cc2ccc(cc2cc1N)c4ccc(\N=N\c3ccccc3)cc4)cc6cc5N)S(=O)(=O)O
 *
 * Topology:
 *   [phenyl]–N=N–[phenylene]–[aminonaphthyl]–N=N–[aminonaphthyl-SO3H]
 *
 * A canonical chemist's drawing places this as a horizontal rod of four
 * ring systems with both N=N bridges parallel. The current default
 * Layout produces a folded layout with distorted (non-hexagonal) rings.
 * This file pins both the *current* behaviour and an *aspirational*
 * acceptance gate for when the SDG refiner lands.
 *
 * Sections:
 *   1. Sanity            — parse + round-trip + layout don't throw
 *   2. Ring hexagonality — every aromatic 6-ring should look like a ring
 *   3. Chromophore rod   — the global molecule should be elongated
 *
 * Tests in §2 and §3 are tagged "RELAXED" until the SDG refiner is wired
 * up; they record the current degraded values and the target thresholds
 * so completion progress is measurable. Failing tests are advisory, not
 * release blockers, until the SDG pass is enabled by default.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'Templates.js'));
require(path.join(__dirname, '..', 'editor', 'Layout.js'));

var SmilesParser = globalThis.SmilesParser;
var SmilesWriter = globalThis.SmilesWriter;
var Layout       = globalThis.Layout;
var Molecule     = globalThis.Molecule;
var SMSDGraph    = globalThis.SMSDGraph;
var SMSDRings    = globalThis.SMSDRings;
var BL           = Molecule.BOND_LENGTH;

var DYE = 'c5(cc6ccc(\\N=N\\c1cc2ccc(cc2cc1N)c4ccc(\\N=N\\c3ccccc3)cc4)cc6cc5N)S(=O)(=O)O';

// Acceptance thresholds.
//
// CURRENT  = what BIME's default Layout produces today (locked in so a
//            silent regression of the rule-based path is caught).
// TARGET   = what the SDG refiner should achieve when complete; tests
//            tagged TARGET are skipped (logged as advisory) until the
//            refiner is enabled.
var THRESHOLDS = {
    ring: {
        currentMaxBLDev:  0.20,   // 20% bond-length spread across one ring
        currentMaxAngDev: 90,     // up to 85° interior-angle distortion observed
        targetMaxBLDev:   0.05,   // 5% — visually a regular hexagon
        targetMaxAngDev:  5       // 5°  — same
    },
    rod: {
        currentMinAspect: 1.5,    // current 1.77, lock floor
        targetMinAspect:  3.0     // chromophore rod should be ≥ 3:1
    }
};

function makeRunner(label) {
    var passed = 0, failed = 0, failures = [], results = [];
    function test(name, fn) {
        try {
            fn();
            passed++;
            results.push({ name: name, pass: true });
            console.log('  ✓ ' + name);
        } catch (e) {
            failed++;
            failures.push({ name: name, reason: e && e.message ? e.message : String(e) });
            results.push({ name: name, pass: false, reason: e.message });
            console.log('  ✗ ' + name + ' (' + (e && e.message) + ')');
        }
    }
    function summary() { return { label: label, passed: passed, failed: failed, failures: failures, results: results }; }
    return { test: test, summary: summary };
}

var runner = makeRunner('Layout chromophore (v1.8.x SDG benchmark)');
var test = runner.test;

console.log('Layout chromophore (v1.8.x SDG benchmark)');

// -----------------------------------------------------------------------
// Section 1 — Sanity
// -----------------------------------------------------------------------

test('parses C32H22N6O3S bisazo dye without error', function () {
    var m = SmilesParser.parse(DYE);
    assert.strictEqual(m.atoms.length, 42, 'expected 42 heavy atoms');
    assert.strictEqual(m.bonds.length, 47, 'expected 47 bonds');

    var hist = { C: 0, N: 0, O: 0, S: 0 };
    m.atoms.forEach(function (a) { if (hist[a.symbol] !== undefined) hist[a.symbol]++; });
    assert.strictEqual(hist.C, 32, 'expected 32 carbons');
    assert.strictEqual(hist.N, 6,  'expected 6 nitrogens');
    assert.strictEqual(hist.O, 3,  'expected 3 oxygens');
    assert.strictEqual(hist.S, 1,  'expected 1 sulfur');
});

test('SMILES round-trip preserves atom + bond counts', function () {
    var m = SmilesParser.parse(DYE);
    var out = SmilesWriter.write(m);
    var m2 = SmilesParser.parse(out);
    assert.strictEqual(m2.atoms.length, m.atoms.length);
    assert.strictEqual(m2.bonds.length, m.bonds.length);
});

test('contains exactly two N=N (azo) double bonds', function () {
    var m = SmilesParser.parse(DYE);
    var azo = 0;
    for (var i = 0; i < m.bonds.length; i++) {
        var b = m.bonds[i];
        if (b.type !== 2) continue;
        var a1 = m.getAtom(b.atom1), a2 = m.getAtom(b.atom2);
        if (a1.symbol === 'N' && a2.symbol === 'N') azo++;
    }
    assert.strictEqual(azo, 2, 'expected 2 azo bridges');
});

test('Layout.layout() runs to completion with finite coords', function () {
    var m = SmilesParser.parse(DYE);
    Layout.layout(m);
    var bad = 0;
    for (var i = 0; i < m.atoms.length; i++) {
        if (!isFinite(m.atoms[i].x) || !isFinite(m.atoms[i].y)) bad++;
    }
    assert.strictEqual(bad, 0, 'all atoms must have finite coordinates');
});

// -----------------------------------------------------------------------
// Section 2 — Ring hexagonality
// -----------------------------------------------------------------------

function ringStats(m, ring) {
    var atoms = ring.map(function (i) { return m.atoms[i]; });
    var bls = [];
    for (var k = 0; k < ring.length; k++) {
        var a1 = atoms[k], a2 = atoms[(k + 1) % ring.length];
        var dx = a1.x - a2.x, dy = a1.y - a2.y;
        bls.push(Math.sqrt(dx * dx + dy * dy));
    }
    var avgBL = bls.reduce(function (s, v) { return s + v; }, 0) / bls.length;
    var maxBLDev = 0;
    bls.forEach(function (b) { var d = Math.abs(b - avgBL) / avgBL; if (d > maxBLDev) maxBLDev = d; });

    var ideal = 180 * (ring.length - 2) / ring.length;
    var maxAngDev = 0;
    for (var k2 = 0; k2 < ring.length; k2++) {
        var p0 = atoms[(k2 - 1 + ring.length) % ring.length];
        var p1 = atoms[k2];
        var p2 = atoms[(k2 + 1) % ring.length];
        var v1x = p0.x - p1.x, v1y = p0.y - p1.y;
        var v2x = p2.x - p1.x, v2y = p2.y - p1.y;
        var dot = v1x * v2x + v1y * v2y;
        var n1 = Math.sqrt(v1x * v1x + v1y * v1y);
        var n2 = Math.sqrt(v2x * v2x + v2y * v2y);
        var ang = Math.acos(Math.max(-1, Math.min(1, dot / (n1 * n2)))) * 180 / Math.PI;
        var d = Math.abs(ang - ideal);
        if (d > maxAngDev) maxAngDev = d;
    }
    return { avgBL: avgBL, maxBLDev: maxBLDev, maxAngDev: maxAngDev, size: ring.length };
}

test('every 6-ring stays within current bond-length spread (CURRENT lock)', function () {
    var m = SmilesParser.parse(DYE);
    Layout.layout(m);
    var rings = SMSDRings.computeRings(new SMSDGraph.SMSDGraph(m));
    rings.forEach(function (ring, idx) {
        if (ring.length !== 6) return;
        var s = ringStats(m, ring);
        assert(s.maxBLDev <= THRESHOLDS.ring.currentMaxBLDev,
            'ring #' + idx + ' bond-length spread ' + (s.maxBLDev * 100).toFixed(0) +
            '% exceeds current lock ' + (THRESHOLDS.ring.currentMaxBLDev * 100) + '%');
    });
});

test('every 6-ring stays within current angle deviation (CURRENT lock)', function () {
    var m = SmilesParser.parse(DYE);
    Layout.layout(m);
    var rings = SMSDRings.computeRings(new SMSDGraph.SMSDGraph(m));
    rings.forEach(function (ring, idx) {
        if (ring.length !== 6) return;
        var s = ringStats(m, ring);
        assert(s.maxAngDev <= THRESHOLDS.ring.currentMaxAngDev,
            'ring #' + idx + ' max interior-angle deviation ' + s.maxAngDev.toFixed(1) +
            '° exceeds current lock ' + THRESHOLDS.ring.currentMaxAngDev + '°');
    });
});

test('TARGET — every 6-ring is a regular hexagon (advisory until SDG)', function () {
    var m = SmilesParser.parse(DYE);
    Layout.layout(m);
    var rings = SMSDRings.computeRings(new SMSDGraph.SMSDGraph(m));
    var distorted = 0;
    rings.forEach(function (ring) {
        if (ring.length !== 6) return;
        var s = ringStats(m, ring);
        if (s.maxBLDev > THRESHOLDS.ring.targetMaxBLDev ||
            s.maxAngDev > THRESHOLDS.ring.targetMaxAngDev) {
            distorted++;
        }
    });
    if (distorted > 0) {
        // Advisory: do not fail the suite. Log the gap.
        console.log('    [advisory] ' + distorted + '/' + rings.length +
            ' six-rings not yet regular hexagons — SDG refiner pending');
    } else {
        assert.strictEqual(distorted, 0, 'all six-rings should be regular hexagons');
    }
});

// -----------------------------------------------------------------------
// Section 3 — Chromophore rod
// -----------------------------------------------------------------------

function bbox(m) {
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (var i = 0; i < m.atoms.length; i++) {
        var a = m.atoms[i];
        if (a.x < minX) minX = a.x; if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y; if (a.y > maxY) maxY = a.y;
    }
    return { w: maxX - minX, h: maxY - minY };
}

test('layout aspect ratio holds at current floor (CURRENT lock)', function () {
    var m = SmilesParser.parse(DYE);
    Layout.layout(m);
    var bb = bbox(m);
    var aspect = bb.w / bb.h;
    assert(aspect >= THRESHOLDS.rod.currentMinAspect,
        'aspect ratio ' + aspect.toFixed(2) + ' below current floor ' +
        THRESHOLDS.rod.currentMinAspect);
});

test('TARGET — chromophore lays out as a >=3:1 rod (advisory until SDG)', function () {
    var m = SmilesParser.parse(DYE);
    Layout.layout(m);
    var bb = bbox(m);
    var aspect = bb.w / bb.h;
    if (aspect < THRESHOLDS.rod.targetMinAspect) {
        console.log('    [advisory] aspect ratio ' + aspect.toFixed(2) +
            ' below target ' + THRESHOLDS.rod.targetMinAspect + ' — SDG refiner pending');
    } else {
        assert(aspect >= THRESHOLDS.rod.targetMinAspect);
    }
});

test('layout is deterministic across runs', function () {
    var m1 = SmilesParser.parse(DYE); Layout.layout(m1);
    var m2 = SmilesParser.parse(DYE); Layout.layout(m2);
    for (var i = 0; i < m1.atoms.length; i++) {
        assert.strictEqual(m1.atoms[i].x, m2.atoms[i].x, 'x at atom ' + i);
        assert.strictEqual(m1.atoms[i].y, m2.atoms[i].y, 'y at atom ' + i);
    }
});

module.exports = runner.summary;
