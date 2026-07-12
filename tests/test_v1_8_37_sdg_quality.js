/**
 * tests/test_v1_8_37_sdg_quality.js
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Guards SDG layout scoring and template best-fit transforms.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));

shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'Templates.js'));
require(path.join(__dirname, '..', 'editor', 'SDGLayout.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'MacroCycleLayout.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'OverlapResolver.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'HydrogenPlacer.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'CorrectGeometricConfiguration.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'LayoutRefiner.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'RingPlacer.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'TemplateHandler.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'IdentityTemplateLibrary.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'NonplanarBonds.js'));
require(path.join(__dirname, '..', 'editor', 'sdg', 'LayoutQuality.js'));
require(path.join(__dirname, '..', 'editor', 'Layout.js'));

var runner = shim.makeRunner('SDG layout quality and template fitting v1.8.37');
var test = runner.test;

var Molecule = globalThis.Molecule;
var SmilesParser = globalThis.SmilesParser;
var Layout = globalThis.Layout;
var Templates = globalThis.Templates;
var LayoutQuality = globalThis.SDG.LayoutQuality;
var TemplateHandler = globalThis.SDG.TemplateHandler;
var BL = Molecule.BOND_LENGTH;

function near(a, b, eps) {
    assert.ok(Math.abs(a - b) <= eps, a + ' not within ' + eps + ' of ' + b);
}

test('LayoutQuality flags impossible bonds and severe atom clashes', function () {
    var mol = new Molecule();
    var a = mol.addAtom('C', 0, 0);
    var b = mol.addAtom('C', BL * 2.2, 0);
    var c = mol.addAtom('O', 1, 1);
    mol.addBond(a.id, b.id, 1);

    var q = LayoutQuality.evaluate(mol);
    assert.ok(q.longBonds >= 1, 'expected long bond to be counted');
    assert.ok(q.severeClashes >= 1 || q.duplicateAtoms >= 1, 'expected close nonbonded atoms to be counted');
    assert.ok(q.hardFailures >= 1, 'expected hard failure');
    assert.ok(q.needsAdaptiveRescue, 'expected rescue flag');
    assert.ok(q.penalty > 100, 'expected meaningful penalty');
});

test('LayoutQuality counts non-shared bond crossings', function () {
    var mol = new Molecule();
    var a = mol.addAtom('C', 0, 0);
    var b = mol.addAtom('C', BL, BL);
    var c = mol.addAtom('C', 0, BL);
    var d = mol.addAtom('C', BL, 0);
    mol.addBond(a.id, b.id, 1);
    mol.addBond(c.id, d.id, 1);

    var q = LayoutQuality.evaluate(mol);
    assert.strictEqual(q.crossings, 1);
    assert.ok(q.penalty >= 6);
});

test('LayoutQuality accepts saturated cage projection crossings', function () {
    var mol = SmilesParser.parse('CC12CC3CC(C1)(CC(C3)(C2)N)C');
    Layout.layout(mol);
    var q = Layout.quality(mol);
    assert.strictEqual(q.crossings, 0);
    assert.strictEqual(q.cageCrossings, 1);
    assert.strictEqual(q.penalty, 0);
});

test('Layout.quality reports clean hard status after normal scaffold layout', function () {
    var mol = SmilesParser.parse('c1ccc2ccccc2c1');
    Layout.layout(mol);
    var q = Layout.quality(mol);

    assert.strictEqual(q.hardFailures, 0);
    assert.strictEqual(q.needsAdaptiveRescue, false);
    assert.strictEqual(q.longBonds, 0);
    assert.strictEqual(q.shortBonds, 0);
});

test('TemplateHandler best-fit rotates and translates a mapped template', function () {
    var mol = new Molecule();
    var a = mol.addAtom('C', 10, 20);
    var b = mol.addAtom('C', 10, 50);
    mol.addBond(a.id, b.id, 1);

    var tmpl = {
        atoms: [
            { symbol: 'C', x: 0, y: 0 },
            { symbol: 'C', x: BL, y: 0 }
        ],
        bonds: [{ a1: 0, a2: 1, type: 1 }]
    };
    var placed = {};
    var transform = TemplateHandler.applyTemplateBestFit(mol, { 0: a.id, 1: b.id }, tmpl, placed);

    assert.ok(transform, 'missing transform');
    assert.strictEqual(placed[a.id], true);
    assert.strictEqual(placed[b.id], true);
    near(mol.getAtom(a.id).x, 10, 1e-6);
    near(mol.getAtom(a.id).y, 20, 1e-6);
    near(mol.getAtom(b.id).x, 10, 1e-6);
    near(mol.getAtom(b.id).y, 50, 1e-6);
});

test('TemplateHandler degenerate target fallback keeps template expanded and centered', function () {
    var mol = new Molecule();
    var tmpl = Templates.benzene();
    var mapping = {};
    for (var i = 0; i < tmpl.atoms.length; i++) {
        var atom = mol.addAtom(tmpl.atoms[i].symbol, 0, 0);
        mapping[i] = atom.id;
    }

    var transform = TemplateHandler.applyTemplateBestFit(mol, mapping, tmpl, {});
    assert.ok(transform && transform.degenerate, 'expected degenerate fallback');

    var cx = 0, cy = 0, maxR = 0;
    for (var ai = 0; ai < mol.atoms.length; ai++) {
        cx += mol.atoms[ai].x;
        cy += mol.atoms[ai].y;
    }
    cx /= mol.atoms.length;
    cy /= mol.atoms.length;
    for (var aj = 0; aj < mol.atoms.length; aj++) {
        var dx = mol.atoms[aj].x - cx;
        var dy = mol.atoms[aj].y - cy;
        var r = Math.sqrt(dx * dx + dy * dy);
        if (r > maxR) maxR = r;
    }

    near(cx, 0, 1e-6);
    near(cy, 0, 1e-6);
    assert.ok(maxR > BL * 0.8, 'template collapsed during degenerate fallback');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
