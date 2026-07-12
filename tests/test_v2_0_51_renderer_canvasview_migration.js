/**
 * tests/test_v2_0_51_renderer_canvasview_migration.js — Renderer.js
 * call-site migration to the v2.0.48 CanvasView adapter layer (v2.0.51).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.51 migrates `Renderer.prototype.screenToMol` and
 * `Renderer.prototype.centerMolecule` to consume the shared
 * `CanvasView` module via the v2.0.48 adapter layer. Behaviour must be
 * byte-identical to v2.0.50: no golden-corpus drift, no pixel motion,
 * no rounding-error introduction.
 *
 * These tests pin the equivalences:
 *
 *   _tx(x)            ≡ CanvasView.worldToScreen(canvasView(), x, _).x
 *   screenToMol math  ≡ CanvasView.screenToWorld(canvasView(), vx, vy)
 *   centerMolecule    ≡ v2.0.50 inline arithmetic (single + reaction cases)
 *
 * Plain Node, no DOM. We duck-type the Renderer (its constructor needs
 * DOM, but the four migrated methods are pure arithmetic on `this`).
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/CanvasView.js');

var runner = shim.makeRunner('Renderer ↔ CanvasView call-site migration (v2.0.51)');
var test = runner.test;

console.log('Renderer ↔ CanvasView call-site migration (v2.0.51)');

// ---------------------------------------------------------------------------
// Build a duck-typed Renderer with just the four migrated methods. We
// can't `new Renderer()` in Node because the constructor touches the
// DOM, but `_tx`/`_ty`/`screenToMol`/`centerMolecule` only read
// numeric fields — so an object with those fields + the v2.0.48
// `canvasView` / `applyCanvasView` methods is sufficient.
// ---------------------------------------------------------------------------

function makeRenderer(opts) {
    opts = opts || {};
    var r = {
        scale:   (typeof opts.scale   === 'number') ? opts.scale   : 1,
        offsetX: (typeof opts.offsetX === 'number') ? opts.offsetX : 0,
        offsetY: (typeof opts.offsetY === 'number') ? opts.offsetY : 0,
        width:   (typeof opts.width   === 'number') ? opts.width   : 800,
        height:  (typeof opts.height  === 'number') ? opts.height  : 600,
        molecule: opts.molecule || null,
        svg: { getBoundingClientRect: function() {
            return { left: 0, top: 0, width: r.width, height: r.height };
        } }
    };
    // Inline copies of the four migrated methods (lifted verbatim from
    // editor/Renderer.js for v2.0.51).
    r._tx = function(x) { return (x + this.offsetX) * this.scale; };
    r._ty = function(y) { return (y + this.offsetY) * this.scale; };
    r.canvasView = function() {
        return { scale: this.scale, offsetX: this.offsetX * this.scale, offsetY: this.offsetY * this.scale };
    };
    r.applyCanvasView = function(view) {
        if (!view || !(view.scale > 0)) { return; }
        this.scale = view.scale;
        this.offsetX = view.offsetX / view.scale;
        this.offsetY = view.offsetY / view.scale;
    };
    r.screenToMol = function(sx, sy) {
        var rect = this.svg.getBoundingClientRect();
        var rw = rect.width || this.width || 1;
        var rh = rect.height || this.height || 1;
        var vx = (sx - rect.left) / rw * this.width;
        var vy = (sy - rect.top) / rh * this.height;
        var w = CanvasView.screenToWorld(this.canvasView(), vx, vy);
        if (w) { return w; }
        return { x: vx / this.scale - this.offsetX, y: vy / this.scale - this.offsetY };
    };
    r.centerMolecule = function() {
        if (!this.molecule || (this.molecule.isEmpty && this.molecule.isEmpty())) return;
        var bounds = this.molecule.getBounds();
        if (this.molecule.reactionArrow) {
            var BL = (typeof opts.bondLength === 'number') ? opts.bondLength : 30;
            var pad = BL * 2;
            var scaleX = (this.width  - pad) / (bounds.w > 0 ? bounds.w : 1);
            var scaleY = (this.height - pad) / (bounds.h > 0 ? bounds.h : 1);
            var fitScale = Math.min(scaleX, scaleY, 1.6);
            fitScale = Math.max(fitScale, 0.3);
            this.scale = fitScale;
        }
        var s = this.scale;
        var offsetScreenX = (this.width  - bounds.w * s) / 2 - bounds.x * s;
        var offsetScreenY = (this.height - bounds.h * s) / 2 - bounds.y * s;
        this.applyCanvasView({ scale: s, offsetX: offsetScreenX, offsetY: offsetScreenY });
    };
    return r;
}

// Fake molecule shape that exposes getBounds() / isEmpty() / reactionArrow.
function makeMol(bounds, hasReactionArrow) {
    return {
        getBounds: function() { return bounds; },
        isEmpty:   function() { return false; },
        reactionArrow: hasReactionArrow || null
    };
}

// ---------------------------------------------------------------------------
// A. _tx / _ty ≡ CanvasView.worldToScreen
// ---------------------------------------------------------------------------
test('A1 _tx(x) equals CanvasView.worldToScreen(canvasView(), x, _).x at identity', function() {
    var r = makeRenderer({ scale: 1, offsetX: 0, offsetY: 0 });
    for (var x = -100; x <= 100; x += 10) {
        var direct = r._tx(x);
        var via    = CanvasView.worldToScreen(r.canvasView(), x, 0).x;
        assert.strictEqual(direct, via,
            '_tx vs CanvasView at x=' + x + ': ' + direct + ' vs ' + via);
    }
});

test('A2 _tx / _ty equal CanvasView.worldToScreen at non-identity view', function() {
    var r = makeRenderer({ scale: 2.5, offsetX: 12, offsetY: -8 });
    for (var x = -50; x <= 50; x += 5) {
        for (var y = -50; y <= 50; y += 5) {
            var sx = r._tx(x);
            var sy = r._ty(y);
            var w  = CanvasView.worldToScreen(r.canvasView(), x, y);
            assert.ok(Math.abs(sx - w.x) < 1e-9, '_tx mismatch at (' + x + ',' + y + ')');
            assert.ok(Math.abs(sy - w.y) < 1e-9, '_ty mismatch at (' + x + ',' + y + ')');
        }
    }
});

// ---------------------------------------------------------------------------
// B. screenToMol ≡ CanvasView.screenToWorld (after viewBox compensation)
// ---------------------------------------------------------------------------
test('B1 screenToMol returns the same world coords as the v2.0.50 inline formula', function() {
    // Build TWO renderers with identical state. One uses the v2.0.50
    // inline math, one uses the v2.0.51 CanvasView path. Their outputs
    // must match for every (sx, sy) probe.
    var r51 = makeRenderer({ scale: 2, offsetX: 10, offsetY: 5, width: 800, height: 600 });
    var r50 = function(sx, sy) {
        var rect = { left: 0, top: 0, width: 800, height: 600 };
        var vx = (sx - rect.left) / 800 * 800;
        var vy = (sy - rect.top) / 600 * 600;
        return { x: vx / 2 - 10, y: vy / 2 - 5 };
    };
    for (var sx = 0; sx <= 800; sx += 50) {
        for (var sy = 0; sy <= 600; sy += 50) {
            var a = r51.screenToMol(sx, sy);
            var b = r50(sx, sy);
            assert.ok(Math.abs(a.x - b.x) < 1e-9, 'x mismatch at (' + sx + ',' + sy + ')');
            assert.ok(Math.abs(a.y - b.y) < 1e-9, 'y mismatch at (' + sx + ',' + sy + ')');
        }
    }
});

test('B2 screenToMol round-trips: screen -> mol -> _tx/_ty -> screen', function() {
    var r = makeRenderer({ scale: 1.75, offsetX: 20, offsetY: -15, width: 800, height: 600 });
    for (var sx = 100; sx <= 700; sx += 100) {
        for (var sy = 100; sy <= 500; sy += 100) {
            var w = r.screenToMol(sx, sy);
            var back = { x: r._tx(w.x), y: r._ty(w.y) };
            assert.ok(Math.abs(back.x - sx) < 1e-9, 'round-trip x at (' + sx + ',' + sy + ')');
            assert.ok(Math.abs(back.y - sy) < 1e-9, 'round-trip y at (' + sx + ',' + sy + ')');
        }
    }
});

// ---------------------------------------------------------------------------
// C. centerMolecule ≡ v2.0.50 inline formula (single molecule)
// ---------------------------------------------------------------------------
test('C1 centerMolecule on a non-reaction molecule matches the v2.0.50 inline result', function() {
    // A single molecule (no reactionArrow) must NOT change scale; only
    // offsetX/offsetY shift to centre the content.
    var bounds = { x: -50, y: -30, w: 100, h: 60 };
    var r = makeRenderer({
        scale: 1, offsetX: 0, offsetY: 0, width: 800, height: 600,
        molecule: makeMol(bounds, null)
    });
    r.centerMolecule();
    var expectedOffsetX = (800 / 1 - 100) / 2 - (-50);   // = 350 + 50 = 400
    var expectedOffsetY = (600 / 1 -  60) / 2 - (-30);   // = 270 + 30 = 300
    assert.strictEqual(r.scale, 1, 'scale unchanged for single molecule');
    assert.ok(Math.abs(r.offsetX - expectedOffsetX) < 1e-9, 'offsetX matches');
    assert.ok(Math.abs(r.offsetY - expectedOffsetY) < 1e-9, 'offsetY matches');
});

test('C2 centerMolecule on a reaction auto-scales and centres', function() {
    // A wide reaction (R + arrow + P) should auto-scale to fit viewport,
    // capped at [0.3, 1.6].
    var bounds = { x: 0, y: 0, w: 500, h: 100 };
    var r = makeRenderer({
        scale: 1, offsetX: 0, offsetY: 0, width: 800, height: 600,
        molecule: makeMol(bounds, { x1: 0, y1: 50, x2: 500, y2: 50 }),
        bondLength: 30
    });
    r.centerMolecule();
    // Expected scale: min((800-60)/500, (600-60)/100, 1.6) = min(1.48, 5.4, 1.6) = 1.48
    var expectedScale = Math.max(Math.min((800 - 60) / 500, (600 - 60) / 100, 1.6), 0.3);
    assert.ok(Math.abs(r.scale - expectedScale) < 1e-9, 'auto-scale correct');
    // Centred at the new scale.
    var s = r.scale;
    var expectedOffsetX = (800 / s - 500) / 2 - 0;
    var expectedOffsetY = (600 / s - 100) / 2 - 0;
    assert.ok(Math.abs(r.offsetX - expectedOffsetX) < 1e-9, 'reaction offsetX');
    assert.ok(Math.abs(r.offsetY - expectedOffsetY) < 1e-9, 'reaction offsetY');
});

test('C3 centerMolecule reaction-scale is clamped to the [0.3, 1.6] band', function() {
    // Tiny content → would naturally scale 100x, must clamp to 1.6.
    var rTiny = makeRenderer({
        scale: 1, width: 800, height: 600,
        molecule: makeMol({ x: 0, y: 0, w: 1, h: 1 }, { x1: 0, y1: 0, x2: 1, y2: 1 }),
        bondLength: 30
    });
    rTiny.centerMolecule();
    assert.strictEqual(rTiny.scale, 1.6, 'tiny reaction clamped to 1.6');

    // Huge content → would naturally scale to ~0.01, must clamp to 0.3.
    var rHuge = makeRenderer({
        scale: 1, width: 800, height: 600,
        molecule: makeMol({ x: 0, y: 0, w: 100000, h: 100000 },
                          { x1: 0, y1: 0, x2: 100000, y2: 100000 }),
        bondLength: 30
    });
    rHuge.centerMolecule();
    assert.strictEqual(rHuge.scale, 0.3, 'huge reaction clamped to 0.3');
});

test('C4 centerMolecule is a no-op on null / empty molecule', function() {
    var r = makeRenderer({ scale: 2, offsetX: 17, offsetY: -3, width: 800, height: 600 });
    r.molecule = null;
    r.centerMolecule();
    assert.strictEqual(r.scale, 2);
    assert.strictEqual(r.offsetX, 17);
    assert.strictEqual(r.offsetY, -3);

    r.molecule = { isEmpty: function() { return true; },
                   getBounds: function() { return { x: 0, y: 0, w: 0, h: 0 }; } };
    r.centerMolecule();
    assert.strictEqual(r.scale, 2);
    assert.strictEqual(r.offsetX, 17);
    assert.strictEqual(r.offsetY, -3);
});

// ---------------------------------------------------------------------------
// D. equivalence to the v2.0.50 inline formulas (regression guard)
// ---------------------------------------------------------------------------
test('D1 centerMolecule offset formula equals v2.0.50 inline arithmetic (single molecule)', function() {
    var probes = [
        { bounds: { x:    0, y:    0, w:  50, h:  50 }, scale: 1 },
        { bounds: { x:  -20, y:  -10, w: 100, h:  60 }, scale: 1 },
        { bounds: { x:  500, y:  300, w: 200, h: 120 }, scale: 1 }
    ];
    for (var i = 0; i < probes.length; i++) {
        var p = probes[i];
        var r = makeRenderer({
            scale: p.scale, width: 800, height: 600,
            molecule: makeMol(p.bounds, null)
        });
        r.centerMolecule();
        var oldX = (800 / p.scale - p.bounds.w) / 2 - p.bounds.x;
        var oldY = (600 / p.scale - p.bounds.h) / 2 - p.bounds.y;
        assert.ok(Math.abs(r.offsetX - oldX) < 1e-9, 'probe ' + i + ' offsetX');
        assert.ok(Math.abs(r.offsetY - oldY) < 1e-9, 'probe ' + i + ' offsetY');
    }
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
