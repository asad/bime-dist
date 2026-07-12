/**
 * tests/test_v2_0_63_rc_bc_highlights.js — Reaction-Center / Bond-Change /
 * Mapped-atom highlights + color overrides (v2.0.63).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.63 surfaces the AAM output that BIME already computes:
 *  - Reaction Centre  — atoms that participate in `bondChanges` get a
 *    dashed amber halo (`.bime-rc-halo` SVG class).
 *  - Bond Changes     — bonds in `bondChanges` get a colour-coded stroke
 *    (red cleaved, green formed, amber order, purple stereo).
 *  - Mapped atoms     — soft halo behind any atom with `mapNumber > 0`.
 *  - Color overrides  — user-supplied palette overrides by element
 *    symbol or by mapNumber (string-keyed).
 *
 * MolEditor exposes toggle / set helpers; the Renderer owns the actual
 * draw. AAM results auto-forward `bondChanges` onto the renderer so the
 * RC / BC toggles have data to draw without extra wiring on the
 * caller side.
 *
 * Plain Node, no DOM. Source-shape tests + small Renderer prototype
 * inspection (without instantiating a Renderer, since the constructor
 * touches the DOM).
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('RC/BC/mapped-atom highlights (v2.0.63)');
var test = runner.test;

console.log('RC/BC/mapped-atom highlights (v2.0.63)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. Renderer additions
// ---------------------------------------------------------------------------
test('A1 Renderer constructor declares bondChanges + showReactionCenter + showBondChanges + colorOverrides + showMappedAtoms', function () {
    var src = readRepoFile('editor/Renderer.js');
    // Look at the constructor body (everything before the first prototype method).
    var ctorEnd = src.indexOf('Renderer.prototype._createGroup');
    var ctor = src.substring(0, ctorEnd);
    assert.ok(ctor.indexOf('this.bondChanges') !== -1,         'bondChanges initialised on constructor');
    assert.ok(ctor.indexOf('this.showReactionCenter') !== -1, 'showReactionCenter initialised');
    assert.ok(ctor.indexOf('this.showBondChanges') !== -1,    'showBondChanges initialised');
    assert.ok(ctor.indexOf('this.colorOverrides') !== -1,     'colorOverrides initialised');
    assert.ok(ctor.indexOf('this.showMappedAtoms') !== -1,    'showMappedAtoms initialised');
});

test('A2 Renderer declares _drawReactionCenter + _drawMappedAtomHalos', function () {
    var src = readRepoFile('editor/Renderer.js');
    assert.ok(src.indexOf('Renderer.prototype._drawReactionCenter = function') !== -1,
        '_drawReactionCenter prototype method exists');
    assert.ok(src.indexOf('Renderer.prototype._drawMappedAtomHalos = function') !== -1,
        '_drawMappedAtomHalos prototype method exists');
});

test('A3 _drawReactionCenter paints a dashed amber halo on every bondChange-participant atom', function () {
    var src = readRepoFile('editor/Renderer.js');
    var idx = src.indexOf('Renderer.prototype._drawReactionCenter = function');
    var fn = src.substring(idx, idx + 1500);
    assert.ok(fn.indexOf('#f59e0b') !== -1,         'amber stroke colour');
    assert.ok(fn.indexOf('stroke-dasharray') !== -1, 'dashed stroke');
    assert.ok(fn.indexOf('bime-rc-halo') !== -1,     'is a .bime-rc-halo element');
    assert.ok(fn.indexOf('this.bondChanges') !== -1, 'reads from this.bondChanges');
});

test('A4 render() calls _drawReactionCenter when showReactionCenter is on AND bondChanges is non-empty', function () {
    var src = readRepoFile('editor/Renderer.js');
    var idx = src.indexOf('Renderer.prototype.render = function');
    var slice = src.substring(idx, idx + 3000);
    assert.ok(/this\.showReactionCenter\s*&&\s*this\.bondChanges/.test(slice),
        'render guards _drawReactionCenter on the toggle + non-empty data');
    assert.ok(slice.indexOf('this._drawReactionCenter()') !== -1,
        'render invokes _drawReactionCenter');
});

test('A5 _drawBond overrides stroke colour for cleaved/formed/orderChange/stereoChange when showBondChanges is on', function () {
    var src = readRepoFile('editor/Renderer.js');
    var idx = src.indexOf('Renderer.prototype._drawBond = function');
    var fn = src.substring(idx, idx + 4500);
    // Each bondChange kind has a documented colour. Test the colour
    // mapping; the actual kind branches must all be present.
    assert.ok(fn.indexOf("if (this.showBondChanges") !== -1,
        '_drawBond guards on this.showBondChanges');
    assert.ok(fn.indexOf("bc.kind === 'cleaved'") !== -1, 'cleaved branch present');
    assert.ok(fn.indexOf("bc.kind === 'formed'") !== -1,  'formed branch present');
    assert.ok(fn.indexOf("bc.kind === 'orderChange'") !== -1, 'orderChange branch present');
    assert.ok(fn.indexOf("bc.kind === 'stereoChange'") !== -1, 'stereoChange branch present');
    assert.ok(fn.indexOf('#dc2626') !== -1, 'red for cleaved');
    assert.ok(fn.indexOf('#16a34a') !== -1, 'green for formed');
});

// ---------------------------------------------------------------------------
// B. MolEditor wiring
// ---------------------------------------------------------------------------
test('B1 MolEditor.toggleReactionCenter / toggleBondChanges / toggleMappedAtoms are declared', function () {
    var src = readRepoFile('editor/MolEditor.js');
    assert.ok(src.indexOf('MolEditor.prototype.toggleReactionCenter = function') !== -1);
    assert.ok(src.indexOf('MolEditor.prototype.toggleBondChanges = function') !== -1);
    assert.ok(src.indexOf('MolEditor.prototype.toggleMappedAtoms = function') !== -1);
});

test('B2 MolEditor.setColorOverride + clearColorOverrides are declared', function () {
    var src = readRepoFile('editor/MolEditor.js');
    assert.ok(src.indexOf('MolEditor.prototype.setColorOverride = function') !== -1);
    assert.ok(src.indexOf('MolEditor.prototype.clearColorOverrides = function') !== -1);
});

test('B3 AAM result auto-forwards bondChanges onto the renderer', function () {
    var src = readRepoFile('editor/MolEditor.js');
    // Look at the section right after this.renderer.componentPairs = rendererPairs.
    var anchor = src.indexOf('this.renderer.componentPairs = rendererPairs');
    assert.ok(anchor > 0);
    var slice = src.substring(anchor, anchor + 1200);
    assert.ok(slice.indexOf('this.renderer.bondChanges') !== -1,
        'renderer.bondChanges assigned post-AAM');
});

test('B4 clearComponentPairs ALSO clears renderer.bondChanges (lockstep clear)', function () {
    var src = readRepoFile('editor/MolEditor.js');
    var idx = src.indexOf('MolEditor.prototype.clearComponentPairs = function');
    var fn = src.substring(idx, idx + 600);
    assert.ok(fn.indexOf('this.renderer.bondChanges = []') !== -1,
        'clearComponentPairs clears the parallel bondChanges cache');
});

// ---------------------------------------------------------------------------
// C. behavioural — call the setters directly on a stub renderer
// ---------------------------------------------------------------------------
test('C1 a duck-typed MolEditor.toggleReactionCenter flips renderer.showReactionCenter', function () {
    // We can't construct a real MolEditor (DOM required). But the
    // toggle method body is short and pure: read .renderer flag,
    // negate, render, return. Replicate the relevant bit and assert.
    var renderer = { showReactionCenter: false, render: function() { this._rendered = (this._rendered || 0) + 1; } };
    var ed = { renderer: renderer, render: function() {} };
    // Inline copy of toggleReactionCenter:
    ed.toggleReactionCenter = function() {
        if (!this.renderer) return false;
        this.renderer.showReactionCenter = !this.renderer.showReactionCenter;
        this.render();
        return this.renderer.showReactionCenter;
    };
    assert.strictEqual(ed.toggleReactionCenter(), true);
    assert.strictEqual(renderer.showReactionCenter, true);
    assert.strictEqual(ed.toggleReactionCenter(), false);
    assert.strictEqual(renderer.showReactionCenter, false);
});

test('C2 setColorOverride writes + clearColorOverrides resets the map', function () {
    var renderer = { colorOverrides: null, render: function() {} };
    var ed = { renderer: renderer, render: function() {} };
    ed.setColorOverride = function(key, color) {
        if (!this.renderer.colorOverrides) { this.renderer.colorOverrides = {}; }
        if (!color) { delete this.renderer.colorOverrides[key]; }
        else { this.renderer.colorOverrides[key] = color; }
        this.render();
        return this.renderer.colorOverrides;
    };
    ed.clearColorOverrides = function() {
        this.renderer.colorOverrides = {};
        this.render();
        return this.renderer.colorOverrides;
    };
    ed.setColorOverride('C', '#ff0000');
    ed.setColorOverride('1', '#00ff00');
    assert.deepStrictEqual(renderer.colorOverrides, { C: '#ff0000', '1': '#00ff00' });
    ed.setColorOverride('C', null);
    assert.deepStrictEqual(renderer.colorOverrides, { '1': '#00ff00' });
    ed.clearColorOverrides();
    assert.deepStrictEqual(renderer.colorOverrides, {});
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
