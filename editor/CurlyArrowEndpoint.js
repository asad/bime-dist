/**
 * editor/CurlyArrowEndpoint.js — shared endpoint-ref helpers for curly
 * arrows (v2.0.33).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Endpoint refs come in two shapes:
 *   { kind: 'atom', atomId: <id> }
 *   { kind: 'bond', bondId: <id> }
 *
 * Both editor/CurlyArrows.js (creating + cloning arrows) and
 * editor/Molecule.js (cloning whole molecules with atom-id + bond-id
 * remapping) had a kind-dispatch sequence each — `if ep.kind === 'bond'
 * { … } else { … atom … }`. v2.0.31's panel flagged the duplication; v2.0.33
 * extracts the helpers here so future endpoint kinds (lone-pair, conformer,
 * etc.) can be added in one place.
 *
 * Four pure functions:
 *
 *   cloneEndpoint(ep)             — shallow copy with same kind + same id
 *   remapEndpoint(ep, atomMap, bondMap)
 *                                  — clone through an id remap (atom or bond)
 *   endpointKey(ep)               — deterministic integer key for parity /
 *                                   hash use, dispatching on kind
 *   normaliseEndpoint(ep, label)  — validate + return a fresh ref; throws
 *                                   a descriptive Error if the shape is bad
 *
 * No DOM access; pure ES5; no external dependencies.
 */
(function (global) {
    'use strict';

    function cloneEndpoint(ep) {
        if (!ep) { return ep; }
        if (ep.kind === 'bond') {
            var c = { kind: 'bond', bondId: ep.bondId };
            // v2.0.36: preserve optional t for biased-along-bond placement.
            if (typeof ep.t === 'number') { c.t = ep.t; }
            return c;
        }
        return { kind: 'atom', atomId: ep.atomId };
    }

    function remapEndpoint(ep, atomMap, bondMap) {
        if (!ep) { return ep; }
        if (ep.kind === 'bond') {
            var r = { kind: 'bond', bondId: (ep.bondId in bondMap) ? bondMap[ep.bondId] : ep.bondId };
            // v2.0.36: t is intrinsic to the arrow (not an id), pass through unchanged.
            if (typeof ep.t === 'number') { r.t = ep.t; }
            return r;
        }
        return { kind: 'atom', atomId: (ep.atomId in atomMap) ? atomMap[ep.atomId] : ep.atomId };
    }

    function endpointKey(ep) {
        if (!ep) { return 0; }
        if (ep.kind === 'bond') { return (ep.bondId | 0); }
        return (ep.atomId | 0);
    }

    function normaliseEndpoint(ep, label) {
        if (!ep || !ep.kind) {
            throw new Error('CurlyArrowEndpoint.normaliseEndpoint: ' + label + '.kind required');
        }
        if (ep.kind === 'atom') {
            if (typeof ep.atomId === 'undefined') {
                throw new Error('CurlyArrowEndpoint.normaliseEndpoint: ' + label + '.atomId required');
            }
            return { kind: 'atom', atomId: ep.atomId };
        }
        if (ep.kind === 'bond') {
            if (typeof ep.bondId === 'undefined') {
                throw new Error('CurlyArrowEndpoint.normaliseEndpoint: ' + label + '.bondId required');
            }
            var out = { kind: 'bond', bondId: ep.bondId };
            // v2.0.36: optional t. Allow caller to pass it through; clamp at
            // resolution time (endpointXY) rather than here so user input
            // round-trips byte-equivalent through serialisation.
            if (typeof ep.t === 'number' && isFinite(ep.t)) { out.t = ep.t; }
            return out;
        }
        throw new Error('CurlyArrowEndpoint.normaliseEndpoint: ' + label + '.kind must be "atom" or "bond"');
    }

    var CurlyArrowEndpoint = {
        cloneEndpoint: cloneEndpoint,
        remapEndpoint: remapEndpoint,
        endpointKey: endpointKey,
        normaliseEndpoint: normaliseEndpoint,
        version: '2.0.36'
    };

    global.CurlyArrowEndpoint = CurlyArrowEndpoint;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = CurlyArrowEndpoint;
    }

})(typeof window !== 'undefined' ? window : this);
