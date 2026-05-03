/**
 * tests/test_layout.js — 2D coordinate generation regression tests
 *
 * Validates that Layout.layout produces clean, finite coordinates for a
 * representative set of polycyclic and acyclic molecules:
 *   - All atoms have finite numeric (x,y) — no NaN, no Infinity, no
 *     absurdly large values (a known bug used to send acyclic sugars to
 *     ±8.7M after collision resolution divided by EPSILON).
 *   - Median bond length is within 10% of BOND_LENGTH (no global scale
 *     distortion).
 *   - No two non-bonded atoms are closer than BOND_LENGTH * 0.4
 *     (no atom collisions).
 *   - For specific scaffolds (sugars, steroids, simple aromatics) bond
 *     crossings stay below documented thresholds.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0.
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

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function buildAndLayout(smi) {
    var mol = SmilesParser.parse(smi);
    if (!mol || !mol.atoms || mol.atoms.length === 0) {
        throw new Error('parse failed for ' + smi);
    }
    Layout.layout(mol);
    return mol;
}

function bondLengths(mol) {
    var lens = [];
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        var a = mol.getAtom(b.atom1);
        var c = mol.getAtom(b.atom2);
        lens.push(Math.sqrt((a.x - c.x) * (a.x - c.x) + (a.y - c.y) * (a.y - c.y)));
    }
    lens.sort(function (a, b) { return a - b; });
    return lens;
}

function median(arr) {
    if (arr.length === 0) return 0;
    return arr[Math.floor(arr.length / 2)];
}

function countCollisions(mol, threshold) {
    var bonded = {};
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        bonded[Math.min(b.atom1, b.atom2) + ':' + Math.max(b.atom1, b.atom2)] = true;
    }
    var count = 0;
    for (var i = 0; i < mol.atoms.length; i++) {
        for (var j = i + 1; j < mol.atoms.length; j++) {
            var a = mol.atoms[i], c = mol.atoms[j];
            var key = Math.min(a.id, c.id) + ':' + Math.max(a.id, c.id);
            if (bonded[key]) continue;
            var dx = a.x - c.x, dy = a.y - c.y;
            if (Math.sqrt(dx * dx + dy * dy) < threshold) count++;
        }
    }
    return count;
}

function countBondCrossings(mol) {
    function ccw(A, B, C) {
        return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x);
    }
    var n = 0;
    for (var i = 0; i < mol.bonds.length; i++) {
        for (var j = i + 1; j < mol.bonds.length; j++) {
            var b1 = mol.bonds[i], b2 = mol.bonds[j];
            if (b1.atom1 === b2.atom1 || b1.atom1 === b2.atom2 ||
                b1.atom2 === b2.atom1 || b1.atom2 === b2.atom2) continue;
            var p1 = mol.getAtom(b1.atom1), p2 = mol.getAtom(b1.atom2);
            var p3 = mol.getAtom(b2.atom1), p4 = mol.getAtom(b2.atom2);
            if (ccw(p1, p3, p4) !== ccw(p2, p3, p4) &&
                ccw(p1, p2, p3) !== ccw(p1, p2, p4)) n++;
        }
    }
    return n;
}

function assertFiniteCoords(mol) {
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        if (typeof a.x !== 'number' || !isFinite(a.x)) {
            throw new Error('atom ' + i + ' x is non-finite: ' + a.x);
        }
        if (typeof a.y !== 'number' || !isFinite(a.y)) {
            throw new Error('atom ' + i + ' y is non-finite: ' + a.y);
        }
        // Sanity: no atom should be more than a few thousand pixels from origin
        if (Math.abs(a.x) > 5000 || Math.abs(a.y) > 5000) {
            throw new Error('atom ' + i + ' coordinates absurd: (' +
                            a.x.toFixed(2) + ', ' + a.y.toFixed(2) + ')');
        }
    }
}

function assertMedianBondLength(mol, tolerance) {
    var lens = bondLengths(mol);
    if (lens.length === 0) return;
    var med = median(lens);
    var ratio = med / BL;
    if (ratio < 1 - tolerance || ratio > 1 + tolerance) {
        throw new Error('median bond length ' + med.toFixed(2) +
                        ' deviates more than ' + (tolerance * 100).toFixed(0) +
                        '% from ' + BL);
    }
}

function assertNoCollisions(mol) {
    var collisions = countCollisions(mol, BL * 0.4);
    if (collisions > 0) {
        throw new Error(collisions + ' non-bonded atom collisions detected');
    }
}

function assertCrossingsAtMost(mol, max) {
    var crossings = countBondCrossings(mol);
    if (crossings > max) {
        throw new Error(crossings + ' bond crossings (max allowed ' + max + ')');
    }
}

// ---------------------------------------------------------------------------
// Test runner
// ---------------------------------------------------------------------------

var runner = shim.makeRunner('Layout (2D coordinate generation)');
var test = runner.test;

// 1. Simple aromatics
test('benzene: 6-ring places cleanly', function () {
    var mol = buildAndLayout('c1ccccc1');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.05);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

test('naphthalene: fused 6-6 places cleanly', function () {
    var mol = buildAndLayout('c1ccc2ccccc2c1');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.05);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

test('biphenyl: two benzenes connected by a single bond', function () {
    var mol = buildAndLayout('c1ccc(-c2ccccc2)cc1');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

// 2. Heterocycles
test('aspirin: substituted benzene with carboxylic acid + ester', function () {
    var mol = buildAndLayout('CC(=O)Oc1ccccc1C(=O)O');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.1);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

test('caffeine: fused 5-6 N-rich heterocycle', function () {
    var mol = buildAndLayout('Cn1cnc2n(C)c(=O)n(C)c(=O)c12');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 1);
});

// 3. Sugars (Haworth-like single ring placement)
test('glucose: pyranose with ring O at sensible position', function () {
    var mol = buildAndLayout('OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.1);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

test('fructose: open-chain ketohexose (no ring)', function () {
    // Regression: this used to send atoms to ±8.7M because of a
    // divide-by-EPSILON in collision resolution.
    var mol = buildAndLayout('OCC(=O)[C@@H](O)[C@H](O)[C@H](O)CO');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.15);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

test('ribose: furanose 5-ring with O', function () {
    var mol = buildAndLayout('OCC1OC(O)C(O)C1O');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.1);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

// 4. Steroids (4 fused rings, classic textbook A-B-C-D scaffold)
test('cholesterol: steroid skeleton with side chain', function () {
    var mol = buildAndLayout('C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 1);
});

test('testosterone: 4 fused rings + ketone', function () {
    var mol = buildAndLayout('CC12CCC3C(C1CCC2O)CCC1=CC(=O)CCC31C');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 2);
});

test('estradiol: aromatic A-ring steroid', function () {
    var mol = buildAndLayout('CC12CCC3c4ccc(O)cc4CCC3C1CCC2O');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 0);
});

// 5. Bridged bicyclic (penicillin, atropine)
test('penicillin G: beta-lactam thiazolidine with acyl side chain', function () {
    var mol = buildAndLayout('CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 1);
});

test('atropine: tropane bicyclic + ester side chain', function () {
    var mol = buildAndLayout('CN1[C@@H]2CC[C@H]1C[C@@H](OC(=O)C(CO)c1ccccc1)C2');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    // Bridged bicyclic + branched chain — allow some crossings
    assertCrossingsAtMost(mol, 4);
});

// 6. Morphine / morphinan family — 5 fused rings, bridged + aromatic
test('morphine: 5 fused rings, finite coords, no collisions', function () {
    var mol = buildAndLayout('CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    // Crossings bounded by 4 (bridged morphinan inherently produces 1-2)
    assertCrossingsAtMost(mol, 4);
});

test('codeine: morphine methyl ether, same skeleton', function () {
    var mol = buildAndLayout('CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5');
    assertFiniteCoords(mol);
    assertNoCollisions(mol);
    assertCrossingsAtMost(mol, 4);
});

// 7. Acyclic and disconnected
test('ethanol: simple acyclic', function () {
    var mol = buildAndLayout('CCO');
    assertFiniteCoords(mol);
    assertMedianBondLength(mol, 0.1);
    assertNoCollisions(mol);
});

test('disconnected: two methanes separated', function () {
    var mol = buildAndLayout('C.C');
    assertFiniteCoords(mol);
});

// Export summary
module.exports = function () { return runner.summary(); };

// If run directly, print results
if (require.main === module) {
    console.log('=== ' + runner.summary().label + ' ===');
}
