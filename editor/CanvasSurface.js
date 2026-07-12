/**
 * editor/CanvasSurface.js — shared interaction-surface primitives for the
 * unified editor canvas (Stage 1 of the molecule-editor / pathway-canvas
 * merge; v2.0.75).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists
 * ---------------
 * The molecule editor (editor/Renderer.js) and the reaction/pathway canvas
 * (js/workbench.js) are being unified into one surface (no mode toggle). The
 * end-state UX is a "focus lens": the pathway is the document, and selecting a
 * structure node expands it in place into an atom-editable inset rendered into
 * that node's <g>. That focus lens needs two primitives the two surfaces do not
 * yet share:
 *
 *   1. transformString / applyTransform — emit the viewport <g> transform
 *      ("translate(ox oy) scale(s)") in ONE place, byte-identical to the
 *      literal both surfaces hand-roll today. Stage 2 reuses this to position
 *      a molecule render inside a node's <g> at node-local world coords.
 *
 *   2. HitRegistry — a geometric (world-coordinate) pick over registered
 *      bounding boxes. The pathway today hit-tests via DOM
 *      `target.closest('[data-pathway-id]')`, which the focus lens cannot use
 *      (it needs to pick a node from a world point, e.g. to decide which node a
 *      drawing gesture lands in). This is the registry for that.
 *
 * This module is intentionally ADDITIVE. It re-exports every editor/CanvasView.js
 * member verbatim (delegating, never forking — the pan/zoom math stays in one
 * place and its symmetry tests are unaffected) and adds the two primitives
 * above. It moves NO existing code: callers consume it behind a
 * `typeof CanvasSurface !== 'undefined'` guard with the legacy inline path kept
 * as a fallback, exactly like the existing CanvasView adoption.
 *
 * Pure ES5; the math/registry are DOM-free and Node-testable. `applyTransform`
 * is the only DOM-touching helper (it calls gEl.setAttribute).
 *
 * Version: tier-2 frozen provenance stamp.
 * NOT surfaced by `bime version`, NOT enrolled in the stamp-lockstep test —
 * frozen at its creation release (2.0.75) until its own logic changes.
 */
(function (global) {
    'use strict';

    // Resolve CanvasView: the browser bundle loads it immediately before this
    // file (tools/editor-files.js), so `global.CanvasView` is set; in Node the
    // global may be unset, so fall back to a direct require. Never reimplement.
    var CV = (global && global.CanvasView) ? global.CanvasView : null;
    if (!CV && typeof require !== 'undefined') {
        try { CV = require('./CanvasView.js'); } catch (e) { CV = null; }
    }

    // --- viewport transform emit ------------------------------------------
    // Byte-identical to the literal at js/workbench.js renderPathwayCanvas:
    //   'translate(' + offsetX + ' ' + offsetY + ') scale(' + scale + ')'
    // No toFixed / rounding — the raw numbers reproduce the exact same string.
    function transformString(view) {
        var v = view || { offsetX: 0, offsetY: 0, scale: 1 };
        var ox = (typeof v.offsetX === 'number') ? v.offsetX : 0;
        var oy = (typeof v.offsetY === 'number') ? v.offsetY : 0;
        var s = (typeof v.scale === 'number') ? v.scale : 1;
        return 'translate(' + ox + ' ' + oy + ') scale(' + s + ')';
    }

    function applyTransform(gEl, view) {
        if (gEl && typeof gEl.setAttribute === 'function') {
            gEl.setAttribute('transform', transformString(view));
        }
        return gEl;
    }

    // --- geometric hit-test registry --------------------------------------
    // Axis-aligned bounding-box pick in WORLD coordinates. pick() returns the
    // topmost hit: items added later sit on top, matching SVG paint order
    // (renderPathwayCanvas draws edges before nodes, nodes before overlays),
    // so the last-registered item covering the point wins.
    function HitRegistry() {
        var items = [];
        return {
            clear: function () { items.length = 0; return this; },
            add: function (id, type, x, y, w, h) {
                items.push({
                    id: id, type: type,
                    x: x, y: y,
                    w: (w >= 0) ? w : 0,
                    h: (h >= 0) ? h : 0
                });
                return this;
            },
            pick: function (wx, wy) {
                for (var i = items.length - 1; i >= 0; i--) {
                    var it = items[i];
                    if (wx >= it.x && wx <= it.x + it.w &&
                        wy >= it.y && wy <= it.y + it.h) {
                        return { id: it.id, type: it.type };
                    }
                }
                return null;
            },
            count: function () { return items.length; },
            all: function () { return items.slice(); }
        };
    }

    var CanvasSurface = {
        // New shared primitives.
        transformString: transformString,
        applyTransform: applyTransform,
        HitRegistry: HitRegistry,
        // Tier-2 frozen provenance stamp (see header).
        version: '2.0.75'
    };

    // Re-export every CanvasView member verbatim (delegate, never fork). Skips
    // keys CanvasSurface already defines, so its own `version` is preserved
    // while CanvasView's math (worldToScreen / screenToWorld / zoom / fit /
    // visibleWorldRect / makeView / DEFAULT_VIEW) is exposed by reference —
    // CanvasSurface.zoom === CanvasView.zoom, etc.
    if (CV) {
        for (var k in CV) {
            if (CV.hasOwnProperty(k) && !CanvasSurface.hasOwnProperty(k)) {
                CanvasSurface[k] = CV[k];
            }
        }
    }

    global.CanvasSurface = CanvasSurface;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = CanvasSurface;
    }

})(typeof window !== 'undefined' ? window : this);
