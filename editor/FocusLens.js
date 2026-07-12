/**
 * editor/FocusLens.js — focus-lens state machine + geometry for the unified
 * editor canvas (Stage 2 of the molecule-editor / pathway-canvas merge; v2.3.0).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists
 * ---------------
 * The molecule editor and the pathway canvas are being unified into ONE surface
 * (the "focus-lens continuum"): the pathway is the document, and selecting a
 * structure node expands it IN PLACE into an atom-editable inset — the molecule
 * editor positioned as an overlay over that node — then re-rasterizes to the
 * node image and collapses on escape/click-away. No mode toggle.
 *
 * This module is the DOM-free heart of that interaction, kept separate (and
 * Node-testable) from the DOM glue in js/workbench.js, exactly as CanvasView /
 * CanvasSurface keep their math DOM-free:
 *
 *   1. FocusLens()        — a tiny state machine owning ONLY which node is
 *                           focused. The glue reads it to (a) position the editor
 *                           overlay and (b) arbitrate event routing: while a lens
 *                           is open, pathway gestures defer to the editor. Atoms
 *                           and pathway nodes never share a hit-test plane; focus
 *                           arbitrates. open()/close() are idempotent.
 *   2. lensRectForNode()  — the on-canvas rectangle a node occupies under the
 *                           current pan/zoom view, in the pathway SVG's viewBox
 *                           coordinate space (worldToScreen of the node's
 *                           top-left corner, size scaled by the view). Pure. The
 *                           glue maps viewBox-space → CSS px via the SVG's
 *                           getBoundingClientRect / viewBox ratio.
 *
 * The Renderer positions every atom by per-element arithmetic (no single
 * viewport <g> to reparent), so the lens is realised as an OVERLAY — the
 * always-mounted editor container repositioned over the focused node — rather
 * than by reparenting the editor's SVG into the node's <g>. lensRectForNode is
 * the geometry that overlay tracks.
 *
 * Pure ES5; no DOM. Depends on CanvasView's world→screen convention by value
 * (re-derived here to stay DOM/require-free), matching renderPathwayNode's use
 * of node.{x,y} as the TOP-LEFT corner in world coords.
 *
 * Version: tier-2 frozen provenance stamp.
 * NOT surfaced by `bime version`, NOT enrolled in the stamp-lockstep test —
 * frozen at its creation release (2.3.0) until its own logic changes.
 */
(function (global) {
    'use strict';

    // ---------------------------------------------------------------------
    // Pure geometry — the viewBox-space rectangle a node occupies under a view.
    //
    // node.{x,y} is the node rect's TOP-LEFT corner in world coords and
    // node.{w,h} its world-space size (see js/workbench.js renderPathwayNode,
    // which draws <rect x=node.x y=node.y width=node.w height=node.h> inside the
    // viewport <g transform="translate(ox oy) scale(s)">). So the on-canvas
    // position is the CanvasView world→screen map of the top-left, and the size
    // scales by the view. Returns {x,y,w,h} in the pathway SVG's viewBox space;
    // the glue converts that to CSS pixels.
    // ---------------------------------------------------------------------
    function lensRectForNode(view, node) {
        if (!node) { return null; }
        var scale = (view && typeof view.scale === 'number') ? view.scale : 1;
        var ox = (view && typeof view.offsetX === 'number') ? view.offsetX : 0;
        var oy = (view && typeof view.offsetY === 'number') ? view.offsetY : 0;
        var nx = (typeof node.x === 'number') ? node.x : 0;
        var ny = (typeof node.y === 'number') ? node.y : 0;
        var nw = (typeof node.w === 'number') ? node.w : 0;
        var nh = (typeof node.h === 'number') ? node.h : 0;
        return {
            x: nx * scale + ox,
            y: ny * scale + oy,
            w: (nw > 0 ? nw : 0) * scale,
            h: (nh > 0 ? nh : 0) * scale
        };
    }

    // ---------------------------------------------------------------------
    // Focus-lens state machine — owns ONLY the focused node id.
    //
    // The glue calls open(id) when a structure node is focused (double-click /
    // Enter), close() on collapse. isOpen()/focused() drive overlay positioning
    // and event arbitration. Deliberately minimal: no DOM, no editor handle, no
    // history — those live in the glue and in EditorHistory/PathwayHistory. A
    // later stage (history unification) adds compose hooks here.
    // ---------------------------------------------------------------------
    function FocusLens() {
        var focusedId = null;
        return {
            open: function (nodeId) {
                if (nodeId === null || typeof nodeId === 'undefined' || nodeId === '') {
                    return this;
                }
                focusedId = nodeId;
                return this;
            },
            close: function () { focusedId = null; return this; },
            focused: function () { return focusedId; },
            isOpen: function () { return focusedId !== null; },
            // Convenience: the lens rect for the focused node under a view.
            rectForNode: function (view, node) { return lensRectForNode(view, node); }
        };
    }

    var FocusLensAPI = {
        // Factory — the glue does `FocusLens.create()` to make its _lens.
        create: FocusLens,
        // Pure geometry, also exposed standalone for testing.
        lensRectForNode: lensRectForNode,
        // Tier-2 frozen provenance stamp (see header).
        version: '2.3.0'
    };

    global.FocusLens = FocusLensAPI;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = FocusLensAPI;
    }

})(typeof window !== 'undefined' ? window : this);
