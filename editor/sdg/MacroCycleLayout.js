/**
 * editor/sdg/MacroCycleLayout.js — CDK MacroCycleLayout port (scaffold).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Reference (Java original, ~14 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/MacroCycleLayout.java
 *
 * Specialised layout for ≥ 8-membered rings. Regular polygons get visually
 * cramped at large ring sizes, so CDK uses an "egg-shape" / oval that
 * makes substituent attachment visible. RingPlacer delegates to this for
 * macrocycles.
 *
 * Status: SCAFFOLD. v1.8.x patches will fill in the CDK-faithful
 * macrocycle ovalisation. BIME's existing Layout.js falls back to a
 * regular polygon for any ring; the result is acceptable up to ~10
 * atoms but visually crowded for cyclodecane and beyond.
 */
(function (global) {
    'use strict';

    var TWO_PI = 2 * Math.PI;

    var MacroCycleLayout = {};

    function _BL() {
        return (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
    }

    // ---------------------------------------------------------------------
    // layout(ring) — place a macrocycle. Stub: regular polygon for now.
    // TODO(v1.8.x): port CDK's egg-shape interpolation.
    // ---------------------------------------------------------------------
    MacroCycleLayout.layout = function (ring, bondLength) {
        if (!ring) return;
        bondLength = bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 8) return; // regular polygon path is fine
        var n = atoms.length;
        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        for (var i = 0; i < n; i++) {
            var ang = i * TWO_PI / n;
            atoms[i].x = radius * Math.cos(ang);
            atoms[i].y = radius * Math.sin(ang);
        }
    };

    global.SDG = global.SDG || {};
    global.SDG.MacroCycleLayout = MacroCycleLayout;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MacroCycleLayout;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
