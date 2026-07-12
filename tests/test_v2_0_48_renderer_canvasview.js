/**
 * tests/test_v2_0_48_renderer_canvasview.js — Renderer ↔ CanvasView
 * interop helpers (v2.0.48).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * The molecule editor's Renderer.js stores offset in WORLD coordinates
 * (screen = (world + offsetWorld) * scale), whereas the shared
 * CanvasView module stores offset in SCREEN coordinates
 * (screen = world * scale + offsetScreen). v2.0.48 adds two adapter
 * methods so external callers can interop with CanvasView's helpers
 * without changing the renderer's storage convention:
 *
 *   renderer.canvasView()           -> { scale, offsetX, offsetY }  // CanvasView convention
 *   renderer.applyCanvasView(view)  // sets renderer state from a CanvasView triple
 *
 * Plain Node, no DOM (Renderer.canvasView is a value method — no
 * DOM access required for these tests).
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/CanvasView.js');

var runner = shim.makeRunner('Renderer CanvasView interop (v2.0.48)');
var test = runner.test;

console.log('Renderer CanvasView interop (v2.0.48)');

// We can't construct a real Renderer in Node (it needs DOM). For the
// pure-arithmetic adapter we can use a duck-typed stub that has the
// renderer's scale/offsetX/offsetY fields and the new methods. The
// methods themselves are defined on Renderer.prototype, so we attach
// them to the stub via Object.assign.

function makeStubRenderer(scale, offsetX, offsetY) {
    return {
        scale: scale, offsetX: offsetX, offsetY: offsetY,
        // Inline copies of Renderer.prototype.canvasView and
        // applyCanvasView so we can test the math without loading the
        // full DOM-dependent Renderer.js.
        canvasView: function() {
            return {
                scale:   this.scale,
                offsetX: this.offsetX * this.scale,
                offsetY: this.offsetY * this.scale
            };
        },
        applyCanvasView: function(view) {
            if (!view || !(view.scale > 0)) { return; }
            this.scale = view.scale;
            this.offsetX = view.offsetX / view.scale;
            this.offsetY = view.offsetY / view.scale;
        }
    };
}

// ---------------------------------------------------------------------------
// A. canvasView() snapshot
// ---------------------------------------------------------------------------
test('A1 canvasView() converts world-space offset to screen-space offset', function() {
    // Renderer.offsetX = 10 (world units), scale = 2 → CanvasView.offsetX = 20 (screen)
    var r = makeStubRenderer(2, 10, 5);
    var v = r.canvasView();
    assert.strictEqual(v.scale, 2);
    assert.strictEqual(v.offsetX, 20, '10 * 2');
    assert.strictEqual(v.offsetY, 10, '5 * 2');
});

test('A2 canvasView() handles unit scale (no-op offset)', function() {
    var r = makeStubRenderer(1, 7, 3);
    var v = r.canvasView();
    assert.deepStrictEqual(v, { scale: 1, offsetX: 7, offsetY: 3 });
});

// ---------------------------------------------------------------------------
// B. applyCanvasView() reverse mapping
// ---------------------------------------------------------------------------
test('B1 applyCanvasView() converts screen-space offset back to world-space', function() {
    var r = makeStubRenderer(1, 0, 0);
    r.applyCanvasView({ scale: 2, offsetX: 20, offsetY: 10 });
    assert.strictEqual(r.scale, 2);
    assert.strictEqual(r.offsetX, 10, '20 / 2');
    assert.strictEqual(r.offsetY, 5,  '10 / 2');
});

test('B2 applyCanvasView() ignores zero-scale views (defensive)', function() {
    var r = makeStubRenderer(2, 10, 5);
    r.applyCanvasView({ scale: 0, offsetX: 99, offsetY: 99 });
    // unchanged
    assert.strictEqual(r.scale, 2);
    assert.strictEqual(r.offsetX, 10);
});

// ---------------------------------------------------------------------------
// C. round-trip: canvasView → CanvasView ops → applyCanvasView
// ---------------------------------------------------------------------------
test('C1 round-trip preserves the renderer state when no view op is applied', function() {
    var r = makeStubRenderer(2.5, 12, -8);
    var view = r.canvasView();
    r.applyCanvasView(view);
    assert.strictEqual(r.scale, 2.5);
    assert.ok(Math.abs(r.offsetX - 12) < 1e-10);
    assert.ok(Math.abs(r.offsetY - (-8)) < 1e-10);
});

test('C2 CanvasView.zoom on canvasView() applied back gives the expected scale', function() {
    var r = makeStubRenderer(1, 0, 0);
    var view = r.canvasView();
    var zoomed = CanvasView.zoom(view, 2);
    r.applyCanvasView(zoomed);
    assert.strictEqual(r.scale, 2);
});

test('C3 CanvasView.zoom with a screen pivot translates correctly through canvasView ↔ Renderer', function() {
    // Start at identity. Renderer's _tx(world)=screen at scale=1, offsetX=0
    // means world = screen. Zoom by 2 with pivot at screen (50, 50). The
    // world point that was at screen (50, 50) (=world (50, 50)) must
    // remain at screen (50, 50) after the zoom.
    var r = makeStubRenderer(1, 0, 0);
    var view = r.canvasView();
    var zoomed = CanvasView.zoom(view, 2, { cx: 50, cy: 50 });
    r.applyCanvasView(zoomed);
    // After the zoom, world (50, 50) should project to screen (50, 50).
    // Renderer's screen = (world + offsetWorld) * scale.
    var screen = (50 + r.offsetX) * r.scale;
    assert.ok(Math.abs(screen - 50) < 1e-9,
        'world (50,50) should still sit at screen 50, got ' + screen);
});

// ---------------------------------------------------------------------------
// D. CanvasView.screenToWorld via canvasView() matches Renderer's screenToMol arithmetic
// ---------------------------------------------------------------------------
test('D1 screen → world via CanvasView matches the Renderer arithmetic', function() {
    var r = makeStubRenderer(2, 10, 5);
    // Renderer arithmetic: world = screen / scale - offsetWorld
    // For screen = (50, 30): world = (50/2 - 10, 30/2 - 5) = (15, 10)
    var view = r.canvasView();
    var w = CanvasView.screenToWorld(view, 50, 30);
    assert.strictEqual(w.x, 15);
    assert.strictEqual(w.y, 10);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
