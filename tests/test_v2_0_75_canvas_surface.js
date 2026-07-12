/**
 * tests/test_v2_0_75_canvas_surface.js — CanvasSurface shared interaction
 * primitives (Stage 1 of the editor/pathway merge; v2.0.75).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * CanvasSurface is an ADDITIVE module: it re-exports editor/CanvasView.js
 * verbatim (the pan/zoom math stays in one place) and adds two primitives the
 * unified focus-lens surface needs — a viewport transform-emit helper and a
 * world-coordinate hit-test registry. These tests pin:
 *   A. re-export identity (CanvasSurface.zoom === CanvasView.zoom, ...), so the
 *      math is delegated, never forked.
 *   B. transformString byte-identity vs the legacy inline literal, so the
 *      pathway viewport transform is unchanged.
 *   C. applyTransform sets that exact string on a target element.
 *   D. HitRegistry geometry (inside / edge / miss / topmost-wins).
 *   E. the tier-2 frozen version stamp (NOT release-tracked).
 *
 * Plain Node, no DOM (applyTransform is tested with a tiny inline recording
 * stub — no shared shim change needed).
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
// Require for side effects only — these set globalThis.CanvasView /
// .CanvasSurface. Under the bundle test runner, require('../editor/<X>.js')
// resolves to the whole bundle (not the individual object), so we read the
// modules off the global exactly like tests/test_v2_0_40_canvas_history.js.
require('../editor/CanvasView.js');
require('../editor/CanvasSurface.js');
var CanvasView = globalThis.CanvasView;
var CanvasSurface = globalThis.CanvasSurface;

var runner = shim.makeRunner('CanvasSurface (v2.0.75)');
var test = runner.test;

console.log('CanvasSurface (v2.0.75)');

// The legacy literal emitted at js/workbench.js renderPathwayCanvas, reproduced
// here so the test fails loudly if transformString ever diverges from it.
function legacyTransform(v) {
    return 'translate(' + v.offsetX + ' ' + v.offsetY + ') scale(' + v.scale + ')';
}

// ---------------------------------------------------------------------------
// A. Re-export identity — CanvasView math is delegated, never forked.
// ---------------------------------------------------------------------------
test('A1 CanvasSurface re-exports every CanvasView function by reference',
    function () {
        var fns = ['worldToScreen', 'screenToWorld', 'zoom', 'fit',
            'visibleWorldRect', 'makeView'];
        for (var i = 0; i < fns.length; i++) {
            assert.strictEqual(typeof CanvasView[fns[i]], 'function',
                'CanvasView.' + fns[i] + ' should exist');
            assert.strictEqual(CanvasSurface[fns[i]], CanvasView[fns[i]],
                'CanvasSurface.' + fns[i] + ' must be the SAME reference as '
                    + 'CanvasView.' + fns[i] + ' (delegate, not fork)');
        }
        assert.strictEqual(CanvasSurface.DEFAULT_VIEW, CanvasView.DEFAULT_VIEW,
            'DEFAULT_VIEW re-exported by reference');
    });

test('A2 delegated math behaves identically through CanvasSurface', function () {
    var v = CanvasSurface.makeView(0.835, 12.5, -7);
    var s = CanvasSurface.worldToScreen(v, 10, 20);
    var s2 = CanvasView.worldToScreen(v, 10, 20);
    assert.deepStrictEqual(s, s2);
    var w = CanvasSurface.screenToWorld(v, s.x, s.y);
    assert.ok(Math.abs(w.x - 10) < 1e-9 && Math.abs(w.y - 20) < 1e-9,
        'screen->world inverts world->screen');
});

// ---------------------------------------------------------------------------
// B. transformString byte-identity vs the legacy inline literal.
// ---------------------------------------------------------------------------
test('B1 transformString is byte-identical to the legacy pathway literal',
    function () {
        var probes = [
            { offsetX: 0, offsetY: 0, scale: 1 },
            { offsetX: 12.5, offsetY: -7, scale: 0.835 },
            { offsetX: -300.25, offsetY: 618, scale: 2.2 },
            { offsetX: 1e-7, offsetY: 0, scale: 1.0000001 },
            { offsetX: 1000000, offsetY: -999999.5, scale: 0.35 }
        ];
        for (var i = 0; i < probes.length; i++) {
            assert.strictEqual(
                CanvasSurface.transformString(probes[i]),
                legacyTransform(probes[i]),
                'transformString must equal the legacy literal for '
                    + JSON.stringify(probes[i]));
        }
    });

test('B2 transformString tolerates a missing/partial view', function () {
    assert.strictEqual(CanvasSurface.transformString(null),
        'translate(0 0) scale(1)');
    assert.strictEqual(CanvasSurface.transformString({ scale: 2 }),
        'translate(0 0) scale(2)');
});

// ---------------------------------------------------------------------------
// C. applyTransform sets the exact string on a target element.
// ---------------------------------------------------------------------------
test('C1 applyTransform writes transformString to the element', function () {
    var recorded = {};
    var el = {
        setAttribute: function (name, value) { recorded[name] = value; }
    };
    var v = { offsetX: 5, offsetY: 6, scale: 2 };
    var ret = CanvasSurface.applyTransform(el, v);
    assert.strictEqual(recorded.transform, CanvasSurface.transformString(v));
    assert.strictEqual(recorded.transform, 'translate(5 6) scale(2)');
    assert.strictEqual(ret, el, 'applyTransform returns the element');
});

test('C2 applyTransform is a no-op on a null/setAttribute-less target',
    function () {
        assert.doesNotThrow(function () { CanvasSurface.applyTransform(null, {}); });
        assert.doesNotThrow(function () { CanvasSurface.applyTransform({}, {}); });
    });

// ---------------------------------------------------------------------------
// D. HitRegistry — world-coordinate bbox pick, topmost wins.
// ---------------------------------------------------------------------------
test('D1 pick returns the item whose bbox contains the point', function () {
    var r = CanvasSurface.HitRegistry();
    r.add('a', 'node', 0, 0, 100, 60);
    assert.deepStrictEqual(r.pick(10, 10), { id: 'a', type: 'node' });
    assert.deepStrictEqual(r.pick(0, 0), { id: 'a', type: 'node' }, 'top-left edge inclusive');
    assert.deepStrictEqual(r.pick(100, 60), { id: 'a', type: 'node' }, 'bottom-right edge inclusive');
    assert.strictEqual(r.pick(101, 10), null, 'just outside → null');
    assert.strictEqual(r.pick(-1, -1), null, 'before origin → null');
});

test('D2 pick returns the TOPMOST (last-added) item on overlap', function () {
    var r = CanvasSurface.HitRegistry();
    r.add('back', 'node', 0, 0, 100, 100);
    r.add('front', 'node', 40, 40, 100, 100);
    assert.deepStrictEqual(r.pick(50, 50), { id: 'front', type: 'node' },
        'overlap → later-added wins (matches paint order)');
    assert.deepStrictEqual(r.pick(10, 10), { id: 'back', type: 'node' },
        'point only in the back item → back');
});

test('D3 clear empties the registry; count/all reflect contents', function () {
    var r = CanvasSurface.HitRegistry();
    r.add('a', 'node', 0, 0, 10, 10).add('b', 'edge', 5, 5, 10, 10);
    assert.strictEqual(r.count(), 2);
    assert.strictEqual(r.all().length, 2);
    r.clear();
    assert.strictEqual(r.count(), 0);
    assert.strictEqual(r.pick(1, 1), null);
});

test('D4 mutating all() does not corrupt the registry', function () {
    var r = CanvasSurface.HitRegistry();
    r.add('a', 'node', 0, 0, 10, 10);
    var snapshot = r.all();
    snapshot.push({ id: 'x' });
    assert.strictEqual(r.count(), 1, 'all() returns a copy');
});

// ---------------------------------------------------------------------------
// E. Version stamp is tier-2 frozen (NOT release-tracked).
// ---------------------------------------------------------------------------
test('E1 CanvasSurface.version is a frozen tier-2 stamp', function () {
    assert.strictEqual(typeof CanvasSurface.version, 'string');
    assert.ok(/^\d+\.\d+\.\d+/.test(CanvasSurface.version));
    // It must NOT inherit CanvasView's version (own stamp preserved).
    assert.notStrictEqual(CanvasSurface.version, CanvasView.version,
        'CanvasSurface keeps its own creation stamp, not CanvasView.version');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
