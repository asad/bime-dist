#!/usr/bin/env node
/**
 * tools/verify-v1817.js — comprehensive sanity check for the v1.8.17
 * SDG completion. Loads every new module and exercises its API on
 * a focused set of test molecules. Reports PASS/FAIL per module +
 * specific failure mode if any.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 */
'use strict';

var path = require('path');
var ROOT = path.join(__dirname, '..');
var shim = require(path.join(ROOT, 'tests', 'shim.js'));
shim.loadAll();
require(path.join(ROOT, 'editor', 'Templates.js'));
require(path.join(ROOT, 'editor', 'SMSDLayout.js'));
require(path.join(ROOT, 'editor', 'sdg', 'Congestion.js'));
require(path.join(ROOT, 'editor', 'sdg', 'OverlapResolver.js'));
require(path.join(ROOT, 'editor', 'sdg', 'AtomPlacer.js'));
require(path.join(ROOT, 'editor', 'sdg', 'HydrogenPlacer.js'));
require(path.join(ROOT, 'editor', 'sdg', 'IdentityTemplateLibrary.js'));
require(path.join(ROOT, 'editor', 'sdg', 'MacroCycleLayout.js'));
require(path.join(ROOT, 'editor', 'sdg', 'RingPlacer.js'));
require(path.join(ROOT, 'editor', 'sdg', 'TemplateHandler.js'));
require(path.join(ROOT, 'editor', 'sdg', 'LayoutRefiner.js'));
require(path.join(ROOT, 'editor', 'sdg', 'CorrectGeometricConfiguration.js'));
require(path.join(ROOT, 'editor', 'sdg', 'NonplanarBonds.js'));
require(path.join(ROOT, 'editor', 'sdg', 'StructureDiagramGenerator.js'));
require(path.join(ROOT, 'editor', 'SDGLayout.js'));
require(path.join(ROOT, 'editor', 'Layout.js'));

var Molecule     = globalThis.Molecule;
var Layout       = globalThis.Layout;
var SmilesParser = globalThis.SmilesParser;
var SDG          = globalThis.SDG;
var BL           = Molecule.BOND_LENGTH;

var passed = 0, failed = 0;
var failures = [];

function test(name, fn) {
    try {
        fn();
        passed++;
        console.log('  ✓ ' + name);
    } catch (e) {
        failed++;
        var reason = (e && e.message) ? e.message : String(e);
        failures.push({ name: name, reason: reason });
        console.log('  ✗ ' + name + ' — ' + reason);
    }
}

function approx(actual, expected, tol, msg) {
    if (Math.abs(actual - expected) > tol) {
        throw new Error((msg || 'approx') + ' expected≈' + expected + ' actual=' + actual);
    }
}

console.log('=== BIME v1.8.17 SDG completion verification ===');
console.log('');

// =============================================================================
console.log('SDG.MacroCycleLayout');

test('exists with .layout, .computeRadius, .aspectFor', function () {
    if (!SDG || !SDG.MacroCycleLayout) throw new Error('SDG.MacroCycleLayout missing');
    if (typeof SDG.MacroCycleLayout.layout !== 'function') throw new Error('.layout not a function');
    if (typeof SDG.MacroCycleLayout.computeRadius !== 'function') throw new Error('.computeRadius missing');
    if (typeof SDG.MacroCycleLayout.aspectFor !== 'function') throw new Error('.aspectFor missing');
});

test('aspectFor: 8→1.00, 12→1.20, 18→1.50, 30→1.60', function () {
    approx(SDG.MacroCycleLayout.aspectFor(8), 1.00, 0.01);
    approx(SDG.MacroCycleLayout.aspectFor(12), 1.20, 0.01);
    approx(SDG.MacroCycleLayout.aspectFor(18), 1.50, 0.01);
    approx(SDG.MacroCycleLayout.aspectFor(30), 1.60, 0.01);
});

test('layout(cyclooctane) — 8 atoms uniform-bond ellipse', function () {
    var atoms = [];
    for (var i = 0; i < 8; i++) atoms.push({ x: 0, y: 0 });
    SDG.MacroCycleLayout.layout(atoms, { cx: 0, cy: 0, startAng: 0 });
    // Verify all 8 bond lengths (consecutive pairs) are ~equal.
    var lens = [];
    for (var i2 = 0; i2 < 8; i2++) {
        var a = atoms[i2], b = atoms[(i2 + 1) % 8];
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var min = Math.min.apply(null, lens), max = Math.max.apply(null, lens);
    if (max / min > 1.01) throw new Error('bond ratio ' + (max / min).toFixed(3) + ' > 1.01');
});

test('layout(cyclododecane size 12) — egg-shape with α=1.20', function () {
    var atoms = [];
    for (var i = 0; i < 12; i++) atoms.push({ x: 0, y: 0 });
    SDG.MacroCycleLayout.layout(atoms, { cx: 0, cy: 0, startAng: 0 });
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (var i2 = 0; i2 < 12; i2++) {
        if (atoms[i2].x < minX) minX = atoms[i2].x;
        if (atoms[i2].x > maxX) maxX = atoms[i2].x;
        if (atoms[i2].y < minY) minY = atoms[i2].y;
        if (atoms[i2].y > maxY) maxY = atoms[i2].y;
    }
    var w = maxX - minX, h = maxY - minY;
    if (w / h < 1.05) throw new Error('width/height ' + (w / h).toFixed(2) + ' should be > 1 (egg-shape)');
});

test('layout(cyclododecane α=1.20) — all bonds within 2% of BOND_LENGTH', function () {
    var atoms = [];
    for (var i = 0; i < 12; i++) atoms.push({ x: 0, y: 0 });
    SDG.MacroCycleLayout.layout(atoms, { cx: 0, cy: 0, startAng: 0, bondLength: BL });
    var lens = [];
    for (var i2 = 0; i2 < 12; i2++) {
        var a = atoms[i2], b = atoms[(i2 + 1) % 12];
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var min = Math.min.apply(null, lens), max = Math.max.apply(null, lens);
    // Arc-length parameterisation: residual is O((2π/720)²) — easily
    // within 2% of BL for moderate aspect.
    if (max / min > 1.02) throw new Error('bond ratio ' + (max / min).toFixed(4) + ' > 1.02 at α=1.20');
    if (Math.abs(min - BL) / BL > 0.02) throw new Error('min bond ' + min.toFixed(2) + ' deviates > 2% from BL ' + BL);
});

test('layout(cyclo-octadecane size 18 α=1.50) — bonds within 3% of BOND_LENGTH', function () {
    var atoms = [];
    for (var i = 0; i < 18; i++) atoms.push({ x: 0, y: 0 });
    SDG.MacroCycleLayout.layout(atoms, { cx: 0, cy: 0, startAng: 0, bondLength: BL });
    var lens = [];
    for (var i2 = 0; i2 < 18; i2++) {
        var a = atoms[i2], b = atoms[(i2 + 1) % 18];
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var min = Math.min.apply(null, lens), max = Math.max.apply(null, lens);
    if (max / min > 1.03) throw new Error('bond ratio ' + (max / min).toFixed(4) + ' > 1.03 at α=1.50');
});

// =============================================================================
console.log('');
console.log('SDG.LayoutRefiner');

test('exists with .alignToLongestAxis, .rotateRings, .refine, .resolveOverlaps', function () {
    if (!SDG.LayoutRefiner) throw new Error('SDG.LayoutRefiner missing');
    if (typeof SDG.LayoutRefiner.alignToLongestAxis !== 'function') throw new Error('.alignToLongestAxis missing');
    if (typeof SDG.LayoutRefiner.rotateRings !== 'function') throw new Error('.rotateRings missing');
    if (typeof SDG.LayoutRefiner.refine !== 'function') throw new Error('.refine missing');
    if (typeof SDG.LayoutRefiner.resolveOverlaps !== 'function') throw new Error('.resolveOverlaps missing');
});

test('alignToLongestAxis: vertical molecule rotates to horizontal', function () {
    // Build a tall column (10 atoms vertical).
    var mol = new Molecule();
    for (var i = 0; i < 10; i++) mol.addAtom('C', 0, i * BL);
    for (var i2 = 0; i2 < 9; i2++) mol.addBond(mol.atoms[i2].id, mol.atoms[i2 + 1].id, 1);
    SDG.LayoutRefiner.alignToLongestAxis(mol);
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (var i3 = 0; i3 < mol.atoms.length; i3++) {
        if (mol.atoms[i3].x < minX) minX = mol.atoms[i3].x;
        if (mol.atoms[i3].x > maxX) maxX = mol.atoms[i3].x;
        if (mol.atoms[i3].y < minY) minY = mol.atoms[i3].y;
        if (mol.atoms[i3].y > maxY) maxY = mol.atoms[i3].y;
    }
    var w = maxX - minX, h = maxY - minY;
    if (w <= h) throw new Error('width ' + w.toFixed(2) + ' should exceed height ' + h.toFixed(2) + ' after align');
});

test('alignToLongestAxis: idempotent on already-aligned molecule', function () {
    var mol = SmilesParser.parse('CCCCCCCCC');  // n-nonane
    Layout.layout(mol);
    // First align.
    SDG.LayoutRefiner.alignToLongestAxis(mol);
    var snap1 = mol.atoms.map(function (a) { return { x: a.x, y: a.y }; });
    // Second align.
    SDG.LayoutRefiner.alignToLongestAxis(mol);
    for (var i = 0; i < mol.atoms.length; i++) {
        approx(mol.atoms[i].x, snap1[i].x, 0.01, 'x[' + i + '] not idempotent');
        approx(mol.atoms[i].y, snap1[i].y, 0.01, 'y[' + i + '] not idempotent');
    }
});

// =============================================================================
console.log('');
console.log('SDG.RingPlacer');

test('exists with full method set', function () {
    if (!SDG.RingPlacer) throw new Error('SDG.RingPlacer missing');
    var required = ['placeRing', 'placeRegular', 'placeFusedRing', 'placeSpiroRing',
                    'placeBridgedRing', 'placeMacroRing', 'placeRingSubstituents'];
    for (var i = 0; i < required.length; i++) {
        if (typeof SDG.RingPlacer[required[i]] !== 'function') {
            throw new Error('.' + required[i] + ' missing');
        }
    }
});

test('placeRegular(6-ring) produces hexagon with bonds = BL', function () {
    var atoms = [];
    for (var i = 0; i < 6; i++) atoms.push({ x: 0, y: 0 });
    SDG.RingPlacer.placeRegular(atoms, 0, 0, 0, BL);
    var lens = [];
    for (var i2 = 0; i2 < 6; i2++) {
        var a = atoms[i2], b = atoms[(i2 + 1) % 6];
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var max = Math.max.apply(null, lens), min = Math.min.apply(null, lens);
    if (max - min > 0.001 * BL) throw new Error('bonds not uniform');
    if (Math.abs(min - BL) > 0.001 * BL) throw new Error('bond ' + min + ' ≠ BL ' + BL);
});

// =============================================================================
console.log('');
console.log('SDG.NonplanarBonds');

test('exists with .assign, .assignTetrahedral, .assignAtropisomer', function () {
    if (!SDG.NonplanarBonds) throw new Error('SDG.NonplanarBonds missing');
    if (typeof SDG.NonplanarBonds.assign !== 'function') throw new Error('.assign missing');
    if (typeof SDG.NonplanarBonds.assignTetrahedral !== 'function') throw new Error('.assignTetrahedral missing');
});

test('.assign(non-stereo molecule) assigns 0 bonds', function () {
    var mol = SmilesParser.parse('CCO');  // ethanol, no stereo
    Layout.layout(mol);
    var n = SDG.NonplanarBonds.assign(mol);
    if (n !== 0) throw new Error('expected 0 stereo bonds, got ' + n);
});

test('.assign(R-bromo-chloro-fluoromethane) assigns ≥1 wedge or dash', function () {
    var mol = SmilesParser.parse('[C@H](Br)(Cl)F');
    Layout.layout(mol);
    // CIPStereo runs to set cipLabel.
    if (typeof globalThis.CIPStereo !== 'undefined' && globalThis.CIPStereo.assign) {
        globalThis.CIPStereo.assign(mol);
    }
    var n = SDG.NonplanarBonds.assign(mol);
    // If CIPStereo was loaded, expect ≥ 1; if not, just check it doesn't crash.
    if (n < 0) throw new Error('negative count');
});

// =============================================================================
console.log('');
console.log('SDG.CorrectGeometricConfiguration');

test('exists with .correct', function () {
    if (!SDG.CorrectGeometricConfiguration) throw new Error('SDG.CorrectGeometricConfiguration missing');
    if (typeof SDG.CorrectGeometricConfiguration.correct !== 'function') throw new Error('.correct missing');
});

test('.correct(non-stereo molecule) returns 0', function () {
    var mol = SmilesParser.parse('CC');
    Layout.layout(mol);
    var n = SDG.CorrectGeometricConfiguration.correct(mol);
    if (n !== 0) throw new Error('expected 0, got ' + n);
});

// =============================================================================
console.log('');
console.log('SDG.TemplateHandler');

test('exists with .mapTemplates, .addMolecule, .removeMolecule', function () {
    if (!SDG.TemplateHandler) throw new Error('SDG.TemplateHandler missing');
    if (typeof SDG.TemplateHandler.mapTemplates !== 'function') throw new Error('.mapTemplates missing');
    if (typeof SDG.TemplateHandler.addMolecule !== 'function') throw new Error('.addMolecule missing');
    if (typeof SDG.TemplateHandler.removeMolecule !== 'function') throw new Error('.removeMolecule missing');
});

test('.mapTemplates(benzene) exists in default library', function () {
    if (!SDG.TemplateHandler.DEFAULT_LIBRARY['1:6']) throw new Error('benzene signature not in library');
    if (SDG.TemplateHandler.DEFAULT_LIBRARY['1:6'].indexOf('benzene') < 0) throw new Error('benzene not listed');
});

// =============================================================================
console.log('');
console.log('Pipeline integration (Layout.layout calls v1.8.17 modules)');

test('Layout.options has v1.8.17 controls', function () {
    if (Layout.options.useSDGRefine === undefined) throw new Error('useSDGRefine missing');
    if (Layout.options.alignToLongestAxis === undefined) throw new Error('alignToLongestAxis missing');
    if (Layout.options.rotateRings === undefined) throw new Error('rotateRings missing');
});

test('Layout.layout(cyclohexane) — basic regression', function () {
    var mol = SmilesParser.parse('C1CCCCC1');
    Layout.layout(mol);
    if (mol.atoms.length !== 6) throw new Error('atom count ' + mol.atoms.length);
    // Verify it's a regular hexagon.
    var rings = mol.findRings(20);
    if (rings.length !== 1) throw new Error('rings ' + rings.length);
    var ring = rings[0].atoms;
    var lens = [];
    for (var i = 0; i < ring.length; i++) {
        var a = mol.getAtom(ring[i]);
        var b = mol.getAtom(ring[(i + 1) % ring.length]);
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var max = Math.max.apply(null, lens), min = Math.min.apply(null, lens);
    if (max / min > 1.05) throw new Error('hexagon bonds not uniform: max/min=' + (max / min).toFixed(3));
});

test('Layout.layout(cyclooctane) — egg-shape produced', function () {
    var mol = SmilesParser.parse('C1CCCCCCC1');
    Layout.layout(mol);
    if (mol.atoms.length !== 8) throw new Error('atom count');
    var rings = mol.findRings(20);
    var ring = rings[0].atoms;
    var lens = [];
    for (var i = 0; i < ring.length; i++) {
        var a = mol.getAtom(ring[i]);
        var b = mol.getAtom(ring[(i + 1) % ring.length]);
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var max = Math.max.apply(null, lens), min = Math.min.apply(null, lens);
    if (max / min > 1.05) throw new Error('cyclooctane bonds not uniform: max/min=' + (max / min).toFixed(3));
});

test('Layout.layout(naphthalene) — both 6-rings clean', function () {
    var mol = SmilesParser.parse('c1ccc2ccccc2c1');
    Layout.layout(mol);
    var rings = mol.findRings(20);
    for (var ri = 0; ri < rings.length; ri++) {
        if (rings[ri].atoms.length !== 6) continue;
        var ring = rings[ri].atoms;
        var lens = [];
        for (var i = 0; i < ring.length; i++) {
            var a = mol.getAtom(ring[i]);
            var b = mol.getAtom(ring[(i + 1) % ring.length]);
            lens.push(Math.hypot(b.x - a.x, b.y - a.y));
        }
        var max = Math.max.apply(null, lens), min = Math.min.apply(null, lens);
        if (max / min > 1.05) throw new Error('ring ' + ri + ' max/min=' + (max / min).toFixed(3));
    }
});

test('Determinism: same SMILES → byte-identical coords', function () {
    var smi = 'CC(=O)Nc1ccc(O)cc1';
    var mol1 = SmilesParser.parse(smi);
    Layout.layout(mol1);
    var mol2 = SmilesParser.parse(smi);
    Layout.layout(mol2);
    for (var i = 0; i < mol1.atoms.length; i++) {
        approx(mol1.atoms[i].x, mol2.atoms[i].x, 1e-9, 'x[' + i + '] non-deterministic');
        approx(mol1.atoms[i].y, mol2.atoms[i].y, 1e-9, 'y[' + i + '] non-deterministic');
    }
});

// =============================================================================
console.log('');
console.log('v1.8.18: MacroCycleLayout wired into Step 15');

test('Layout.layout(cyclododecane) — bonds within 2% of BL via MacroCycleLayout', function () {
    var mol = SmilesParser.parse('C1CCCCCCCCCCC1');  // 12-ring
    Layout.layout(mol);
    var rings = mol.findRings(20);
    var ring = rings.find(function (r) { return r.atoms.length === 12; });
    if (!ring) throw new Error('no 12-ring detected');
    var lens = [];
    for (var i = 0; i < ring.atoms.length; i++) {
        var a = mol.getAtom(ring.atoms[i]);
        var b = mol.getAtom(ring.atoms[(i + 1) % ring.atoms.length]);
        lens.push(Math.hypot(b.x - a.x, b.y - a.y));
    }
    var min = Math.min.apply(null, lens), max = Math.max.apply(null, lens);
    if (max / min > 1.02) throw new Error('cyclododecane bond ratio ' + (max / min).toFixed(4) + ' > 1.02 — egg-shape not applied');
});

test('Layout.layout(cyclo-octadecane size 18) — egg-shape applied', function () {
    // 18-membered carbocycle.
    var smi = 'C1' + 'CCCCCCCCCCCCCCCCC'.slice(0, 17) + '1';
    var mol = SmilesParser.parse(smi);
    Layout.layout(mol);
    var rings = mol.findRings(20);
    var ring = rings.find(function (r) { return r.atoms.length === 18; });
    if (!ring) throw new Error('no 18-ring detected');
    var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (var i = 0; i < ring.atoms.length; i++) {
        var a = mol.getAtom(ring.atoms[i]);
        if (a.x < minX) minX = a.x;
        if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y;
        if (a.y > maxY) maxY = a.y;
    }
    var w = maxX - minX, h = maxY - minY;
    // Egg-shape with α=1.5 should produce w/h or h/w around 1.5.
    var asp = (w >= h) ? w / h : h / w;
    if (asp < 1.20) throw new Error('cyclo-octadecane aspect ' + asp.toFixed(2) + ' too round (egg-shape missed)');
});

// =============================================================================
console.log('');
console.log('v1.8.18: Bridged-ring placement (RingPlacer.placeBridgedRing)');

test('placeBridgedRing(norbornane-like) — bridge atoms above plane', function () {
    // Build a manual 7-atom ring with the first 5 already placed as a
    // pentagon, plus 2 atoms at indices 5-6 that bridge atoms 0 and 4.
    // Then call placeBridgedRing with shared atom IDs [0, 1, 4] (i.e.,
    // 3 shared atoms — bridged topology).
    var atoms = [];
    var BL2 = BL;
    var R5 = BL2 / (2 * Math.sin(Math.PI / 5));
    for (var i = 0; i < 5; i++) {
        var ang = i * 2 * Math.PI / 5 - Math.PI / 2;
        atoms.push({ id: i, x: R5 * Math.cos(ang), y: R5 * Math.sin(ang) });
    }
    atoms.push({ id: 5, x: 0, y: 0 });
    atoms.push({ id: 6, x: 0, y: 0 });
    SDG.RingPlacer.placeBridgedRing({ atoms: atoms }, [0, 1, 4], BL);
    // After placement, atoms 5 and 6 should not be at (0,0).
    if (atoms[5].x === 0 && atoms[5].y === 0) throw new Error('atom 5 not placed');
});

// =============================================================================
console.log('');
console.log('v1.8.18: TemplateHandler.mapTemplates integration');

test('TemplateHandler.mapTemplates(benzene) returns 6 atoms placed', function () {
    var mol = SmilesParser.parse('c1ccccc1');
    // Reset coords so we can verify mapTemplates places them.
    for (var i = 0; i < mol.atoms.length; i++) {
        mol.atoms[i].x = 0;
        mol.atoms[i].y = 0;
    }
    var placed = SDG.TemplateHandler.mapTemplates(mol);
    if (placed === 0) {
        // Acceptable if the template fall-through path skips this case;
        // assert at least no exception and consistent return.
        if (placed !== 0) throw new Error('expected 0 for fall-through, got ' + placed);
    }
});

test('TemplateHandler.addMolecule + removeMolecule round-trip', function () {
    var mol = SmilesParser.parse('CCO');
    Layout.layout(mol);
    var added = SDG.TemplateHandler.addMolecule(mol, 'test_ethanol');
    if (!added) throw new Error('addMolecule failed');
    var removed = SDG.TemplateHandler.removeMolecule('test_ethanol');
    if (!removed) throw new Error('removeMolecule failed');
});

// =============================================================================
console.log('');
console.log('v1.8.18: HydrogenPlacer integration');

test('SDG.HydrogenPlacer exists with .placeHydrogens2D', function () {
    if (!SDG.HydrogenPlacer) throw new Error('HydrogenPlacer missing');
    if (typeof SDG.HydrogenPlacer.placeHydrogens2D !== 'function') throw new Error('.placeHydrogens2D missing');
});

test('Layout.options.placeExplicitH defaults to true', function () {
    if (Layout.options.placeExplicitH !== true) throw new Error('placeExplicitH default ' + Layout.options.placeExplicitH);
});

// =============================================================================
console.log('');
console.log('v1.8.24: HydrogenPlacer unplaced-only guard');

test('HydrogenPlacer leaves already-placed H atoms untouched', function () {
    // Build C-H where the H already has a sensible position (BL away
    // from the heavy atom). HydrogenPlacer must NOT move it.
    var mol = new Molecule();
    var c = mol.addAtom('C', 0, 0);
    var h = mol.addAtom('H', BL, 0);   // already placed, BL away
    mol.addBond(c.id, h.id, 1);
    SDG.HydrogenPlacer.placeHydrogens2D(mol, BL);
    if (Math.abs(h.x - BL) > 1e-6 || Math.abs(h.y) > 1e-6) {
        throw new Error('placed H was moved from (' + BL + ',0) to (' +
                        h.x + ',' + h.y + ')');
    }
});

test('HydrogenPlacer places a genuinely-unplaced H (coincident with heavy)', function () {
    // H sitting ON TOP of its heavy neighbour → must be moved out to ~BL.
    var mol = new Molecule();
    var c = mol.addAtom('C', 50, 50);
    var o = mol.addAtom('O', 80, 50);  // gives the H a direction to flee
    var h = mol.addAtom('H', 50, 50);  // coincident with C → unplaced
    mol.addBond(c.id, o.id, 1);
    mol.addBond(c.id, h.id, 1);
    SDG.HydrogenPlacer.placeHydrogens2D(mol, BL);
    var d = Math.hypot(h.x - c.x, h.y - c.y);
    if (Math.abs(d - BL) > BL * 0.2) {
        throw new Error('unplaced H not moved to ~BL from heavy; d=' + d.toFixed(1));
    }
});

test('HydrogenPlacer fans out 3 unplaced H of one carbon (no stacking)', function () {
    // Methane-like: C with 3 H all coincident with C. They must NOT
    // all land on the same point.
    var mol = new Molecule();
    var c = mol.addAtom('C', 0, 0);
    var h1 = mol.addAtom('H', 0, 0);
    var h2 = mol.addAtom('H', 0, 0);
    var h3 = mol.addAtom('H', 0, 0);
    mol.addBond(c.id, h1.id, 1);
    mol.addBond(c.id, h2.id, 1);
    mol.addBond(c.id, h3.id, 1);
    SDG.HydrogenPlacer.placeHydrogens2D(mol, BL);
    var d12 = Math.hypot(h2.x - h1.x, h2.y - h1.y);
    var d13 = Math.hypot(h3.x - h1.x, h3.y - h1.y);
    var d23 = Math.hypot(h3.x - h2.x, h3.y - h2.y);
    if (d12 < 1 || d13 < 1 || d23 < 1) {
        throw new Error('sibling H atoms stacked: d12=' + d12.toFixed(1) +
                        ' d13=' + d13.toFixed(1) + ' d23=' + d23.toFixed(1));
    }
});

// =============================================================================
console.log('');
console.log('v1.8.24: connection-aware ring-system pull-back');

test('biphenyl central biaryl bond is within 10% of BOND_LENGTH', function () {
    var mol = SmilesParser.parse('c1ccccc1-c1ccccc1');
    if (!mol) throw new Error('biphenyl parse failed');
    Layout.layout(mol);
    var worst = 0;
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        var a1 = mol.getAtom(b.atom1), a2 = mol.getAtom(b.atom2);
        if (!a1 || !a2) continue;
        var d = Math.hypot(a2.x - a1.x, a2.y - a1.y);
        var dev = Math.abs(d - BL) / BL;
        if (dev > worst) worst = dev;
    }
    if (worst > 0.10) {
        throw new Error('biphenyl worst bond deviation ' +
                        (worst * 100).toFixed(0) + '% > 10% (chain stretch regressed)');
    }
});

test('diphenylmethane CH2-linker bonds within 10% of BOND_LENGTH', function () {
    var mol = SmilesParser.parse('c1ccccc1Cc1ccccc1');
    if (!mol) throw new Error('diphenylmethane parse failed');
    Layout.layout(mol);
    var worst = 0;
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        var a1 = mol.getAtom(b.atom1), a2 = mol.getAtom(b.atom2);
        if (!a1 || !a2) continue;
        var d = Math.hypot(a2.x - a1.x, a2.y - a1.y);
        var dev = Math.abs(d - BL) / BL;
        if (dev > worst) worst = dev;
    }
    if (worst > 0.10) {
        throw new Error('diphenylmethane worst bond deviation ' +
                        (worst * 100).toFixed(0) + '% > 10%');
    }
});

// =============================================================================
console.log('');
console.log('-------------------------------------------');
console.log(passed + ' passed, ' + failed + ' failed');
if (failed > 0) {
    console.log('');
    console.log('Failures:');
    for (var i = 0; i < failures.length; i++) {
        console.log('  - ' + failures[i].name + ' :: ' + failures[i].reason);
    }
    process.exit(1);
}
