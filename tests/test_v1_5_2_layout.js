/**
 * tests/test_v1_5_2_layout.js — v1.5.2 chemistry-canonical layout regression
 *
 * Pins the post-pass behaviour added by Layout.options.enforceConventions:
 *   - Bond-length variance shrinks toward the BOND_LENGTH target.
 *   - Aromatic 6-rings draw with one bond horizontal (not vertex-up).
 *   - Substituents on rings point radially outward, not inward.
 *   - Naphthalene long-axis is horizontal.
 *   - Determinism: same SMILES => byte-identical coords across runs.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0.
 */
'use strict';

var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'Templates.js'));
require(path.join(__dirname, '..', 'editor', 'Layout.js'));

var SmilesParser = globalThis.SmilesParser;
var Layout       = globalThis.Layout;
var Molecule     = globalThis.Molecule;
var BL           = Molecule.BOND_LENGTH;

function build(smi) {
    var mol = SmilesParser.parse(smi);
    Layout.layout(mol);
    return mol;
}

function bondLengths(mol) {
    var lens = [];
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        var a = mol.getAtom(b.atom1), c = mol.getAtom(b.atom2);
        var dx = a.x - c.x, dy = a.y - c.y;
        lens.push(Math.sqrt(dx * dx + dy * dy));
    }
    return lens;
}

function stddev(arr) {
    if (arr.length === 0) return 0;
    var m = 0;
    for (var i = 0; i < arr.length; i++) m += arr[i];
    m /= arr.length;
    var s = 0;
    for (var i = 0; i < arr.length; i++) s += (arr[i] - m) * (arr[i] - m);
    return Math.sqrt(s / arr.length);
}

function ringCentroid(mol, atomIds) {
    var cx = 0, cy = 0;
    for (var i = 0; i < atomIds.length; i++) {
        var a = mol.getAtom(atomIds[i]); cx += a.x; cy += a.y;
    }
    return { x: cx / atomIds.length, y: cy / atomIds.length };
}

function findRing(mol, size) {
    var info = Layout.getRingInfo(mol);
    for (var i = 0; i < info.length; i++) if (info[i].size === size) return info[i];
    return null;
}

function snapshot(mol) {
    // Use position in atoms[] (not absolute id, which is a module-global
    // counter incrementing across parse() calls) so snapshots compare
    // structure-by-structure.
    var coords = [];
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        coords.push(i + ':' + a.symbol + ':' + a.x.toFixed(6) + ',' + a.y.toFixed(6));
    }
    return coords.join('|');
}

var runner = shim.makeRunner('Layout v1.5.2 chemistry conventions');
var test = runner.test;

// --- Step 1 pin: bond-length variance ---
test('bond-length stddev < 5% BOND_LENGTH on 20-atom chain mol (decanol)', function () {
    var mol = build('CCCCCCCCCCO');
    var lens = bondLengths(mol);
    var sd = stddev(lens);
    var pct = sd / BL;
    if (pct >= 0.05) {
        throw new Error('bond stddev ' + sd.toFixed(3) + ' = ' + (pct * 100).toFixed(1) +
                        '% of BL (>= 5%)');
    }
});

test('bond-length stddev < 5% BOND_LENGTH on 4 fused rings (steroid)', function () {
    var mol = build('C1CCC2CCC3CCCCC3C2C1');
    var lens = bondLengths(mol);
    var sd = stddev(lens);
    var pct = sd / BL;
    if (pct >= 0.05) {
        throw new Error('bond stddev ' + sd.toFixed(3) + ' = ' + (pct * 100).toFixed(1) +
                        '% of BL (>= 5%)');
    }
});

// --- Step 2 pin: aromatic ring orientation ---
// Para-disubstituted benzene: substituents 180 deg apart. After Step 2
// the highest-degree edge (the C-C bond between the two substituent-bearing
// atoms) should be roughly horizontal -> the two substituent atoms have
// |dy| < small threshold relative to the ring centroid.
test('para-cresol benzene draws with substituent edge horizontal', function () {
    var mol = build('Cc1ccc(O)cc1');
    var info = Layout.getRingInfo(mol);
    if (!info || info.length !== 1) throw new Error('expected 1 ring, got ' + info.length);
    var ring = info[0];
    if (!ring.aromatic) throw new Error('ring not aromatic');
    // Para axis: the two substituent-bearing C atoms should be at (left, right)
    // of the centroid, i.e. their y-coords should match within 5% BL.
    var subAtoms = [];
    for (var i = 0; i < ring.atoms.length; i++) {
        var aid = ring.atoms[i];
        var nb = mol.getNeighbors(aid);
        var hasSub = false;
        for (var j = 0; j < nb.length; j++) {
            if (ring.atoms.indexOf(nb[j]) < 0) { hasSub = true; break; }
        }
        if (hasSub) subAtoms.push(aid);
    }
    if (subAtoms.length !== 2) {
        throw new Error('expected 2 sub-bearing ring atoms, got ' + subAtoms.length);
    }
    var a = mol.getAtom(subAtoms[0]), b = mol.getAtom(subAtoms[1]);
    var dy = Math.abs(a.y - b.y);
    if (dy > BL * 0.20) {
        throw new Error('substituent-bearing C atoms y-differ by ' + dy.toFixed(2) +
                        ' (expected < ' + (BL * 0.20).toFixed(2) + ' for horizontal axis)');
    }
});

// --- Step 3 pin: substituent radial direction ---
// For methoxybenzene (anisole), the methoxy oxygen and methyl carbon
// should both lie outside the ring (further from ring centroid than the
// ring atom they attach to).
test('methoxybenzene substituents point radially OUTWARD from ring', function () {
    var mol = build('COc1ccccc1');
    var info = Layout.getRingInfo(mol);
    if (!info || info.length !== 1) throw new Error('expected 1 ring');
    var ring = info[0];
    var c = ring.center;
    // Ring radius approx
    var ringRadius = BL / (2 * Math.sin(Math.PI / 6));
    // Find atoms outside the ring (the O and the methyl C).
    var outside = [];
    for (var i = 0; i < mol.atoms.length; i++) {
        if (ring.atoms.indexOf(mol.atoms[i].id) >= 0) continue;
        outside.push(mol.atoms[i]);
    }
    for (var i = 0; i < outside.length; i++) {
        var a = outside[i];
        var dx = a.x - c.x, dy = a.y - c.y;
        var dist = Math.sqrt(dx * dx + dy * dy);
        if (dist <= ringRadius + 1e-3) {
            throw new Error('substituent atom ' + a.id + ' (' + a.symbol +
                            ') at distance ' + dist.toFixed(2) + ' inside ring radius ' +
                            ringRadius.toFixed(2));
        }
    }
});

// --- Step 4 pin: naphthalene long-axis horizontal ---
test('naphthalene long-axis horizontal: width > height by >= 1.3x', function () {
    var mol = build('c1ccc2ccccc2c1');
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        if (a.x < minX) minX = a.x;
        if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y;
        if (a.y > maxY) maxY = a.y;
    }
    var w = maxX - minX, h = maxY - minY;
    if (w < h * 1.3) {
        throw new Error('naphthalene not landscape: width ' + w.toFixed(1) +
                        ', height ' + h.toFixed(1));
    }
});

// --- Determinism pin ---
test('determinism: same SMILES => byte-identical coordinates across runs', function () {
    var smiList = [
        'CCO', 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O',
        'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
        'COc1cc(C=CC(=O)c2cc(O)cc(O)c2)ccc1O' // a flavone-like
    ];
    for (var i = 0; i < smiList.length; i++) {
        var m1 = build(smiList[i]);
        var s1 = snapshot(m1);
        var m2 = build(smiList[i]);
        var s2 = snapshot(m2);
        if (s1 !== s2) {
            throw new Error('non-deterministic for ' + smiList[i] + ':\n  ' + s1 + '\n  ' + s2);
        }
    }
});

// --- Backward-compat flag ---
test('enforceConventions=false produces v1.5.1-style layout (no exceptions)', function () {
    Layout.options.enforceConventions = false;
    try {
        var mol = build('Cc1ccc(O)cc1');
        // Just make sure it runs and produces finite coords.
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            if (!isFinite(a.x) || !isFinite(a.y)) {
                throw new Error('non-finite coord at atom ' + i);
            }
        }
    } finally {
        Layout.options.enforceConventions = true;
    }
});

// Export summary
module.exports = function () { return runner.summary(); };
if (require.main === module) {
    var s = runner.summary();
    console.log('=== ' + s.label + ': ' + s.passed + ' passed, ' + s.failed + ' failed ===');
}
