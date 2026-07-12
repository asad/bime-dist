/**
 * editor/CanvasView.js — shared pan / zoom / fit viewport math
 * (R3-phase-2; v2.0.40).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Both the molecule editor (`editor/Renderer.js`) and the pathway canvas
 * (`js/workbench.js`) maintain a `{scale, offsetX, offsetY}` viewport
 * triple and apply the same pan/zoom transforms. v2.0.40 extracts the
 * arithmetic into a single module so the same numerical edge cases are
 * handled identically in both surfaces and so any new tool that needs
 * "world → screen" or "screen → world" can ask one helper.
 *
 * Pure data model — no DOM. Each function is side-effect-free; callers
 * pass a `view` triple and receive a new triple (or a scalar). Trivially
 * testable in Node.
 *
 * Convention: `screen = world * scale + offset`.
 * "world" = molecule / pathway intrinsic coords (atoms, pathway nodes).
 * "screen" = SVG / DOM coords (post-transform).
 *
 * No external dependencies; pure ES5; no DOM access.
 */
(function (global) {
    'use strict';

    var DEFAULT_VIEW = { scale: 1, offsetX: 0, offsetY: 0 };

    function makeView(scale, offsetX, offsetY) {
        return {
            scale:   (typeof scale   === 'number' && isFinite(scale)   && scale > 0) ? scale   : 1,
            offsetX: (typeof offsetX === 'number' && isFinite(offsetX)) ? offsetX : 0,
            offsetY: (typeof offsetY === 'number' && isFinite(offsetY)) ? offsetY : 0
        };
    }

    function worldToScreen(view, x, y) {
        var v = view || DEFAULT_VIEW;
        return { x: x * v.scale + v.offsetX, y: y * v.scale + v.offsetY };
    }

    function screenToWorld(view, x, y) {
        var v = view || DEFAULT_VIEW;
        if (!(v.scale > 0)) { return null; }
        return { x: (x - v.offsetX) / v.scale, y: (y - v.offsetY) / v.scale };
    }

    function zoom(view, factor, opts) {
        var v = view || DEFAULT_VIEW;
        opts = opts || {};
        if (!(factor > 0) || !isFinite(factor)) { return makeView(v.scale, v.offsetX, v.offsetY); }
        var newScale = v.scale * factor;
        if (typeof opts.minScale === 'number' && newScale < opts.minScale) { newScale = opts.minScale; }
        if (typeof opts.maxScale === 'number' && newScale > opts.maxScale) { newScale = opts.maxScale; }
        if (newScale <= 0) { newScale = v.scale; }
        var newOffsetX = v.offsetX;
        var newOffsetY = v.offsetY;
        if (typeof opts.cx === 'number' && typeof opts.cy === 'number' && v.scale > 0) {
            // Keep the world point under (cx, cy) stationary across the zoom.
            var wx = (opts.cx - v.offsetX) / v.scale;
            var wy = (opts.cy - v.offsetY) / v.scale;
            newOffsetX = opts.cx - wx * newScale;
            newOffsetY = opts.cy - wy * newScale;
        }
        return makeView(newScale, newOffsetX, newOffsetY);
    }

    function fit(worldBounds, viewport, opts) {
        opts = opts || {};
        var b = worldBounds;
        var vw = viewport && viewport.w;
        var vh = viewport && viewport.h;
        if (!b || !(b.w > 0) || !(b.h > 0) || !(vw > 0) || !(vh > 0)) {
            return makeView(1, 0, 0);
        }
        var padding = (typeof opts.padding === 'number') ? opts.padding : 0.05;
        if (padding < 0) { padding = 0; } else if (padding > 0.45) { padding = 0.45; }
        var usableW = vw * (1 - 2 * padding);
        var usableH = vh * (1 - 2 * padding);
        var sx = usableW / b.w;
        var sy = usableH / b.h;
        var scale = (sx < sy) ? sx : sy;
        if (!(scale > 0) || !isFinite(scale)) { return makeView(1, 0, 0); }
        var contentScreenW = b.w * scale;
        var contentScreenH = b.h * scale;
        var offsetX = (vw - contentScreenW) / 2 - b.x * scale;
        var offsetY = (vh - contentScreenH) / 2 - b.y * scale;
        return makeView(scale, offsetX, offsetY);
    }

    function visibleWorldRect(view, viewport) {
        var v = view || DEFAULT_VIEW;
        var vw = (viewport && viewport.w) || 0;
        var vh = (viewport && viewport.h) || 0;
        if (!(v.scale > 0) || !(vw > 0) || !(vh > 0)) { return { x: 0, y: 0, w: 0, h: 0 }; }
        return {
            x: -v.offsetX / v.scale,
            y: -v.offsetY / v.scale,
            w: vw / v.scale,
            h: vh / v.scale
        };
    }

    var CanvasView = {
        DEFAULT_VIEW: DEFAULT_VIEW,
        makeView: makeView,
        worldToScreen: worldToScreen,
        screenToWorld: screenToWorld,
        zoom: zoom,
        fit: fit,
        visibleWorldRect: visibleWorldRect,
        version: '2.0.40'
    };

    global.CanvasView = CanvasView;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = CanvasView;
    }

})(typeof window !== 'undefined' ? window : this);
