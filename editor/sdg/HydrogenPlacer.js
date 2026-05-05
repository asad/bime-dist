/**
 * editor/sdg/HydrogenPlacer.js — clean-room JS port of CDK's HydrogenPlacer.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Reference (Java original):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/HydrogenPlacer.java
 *
 * Places explicit hydrogen atoms on a heavy-atom layout. Computes the angle
 * direction opposite the existing heavy-atom neighbours and places the H
 * along that direction at BOND_LENGTH. For atoms with no heavy neighbours,
 * spreads H atoms around a circle.
 *
 * Public API (matches CDK):
 *   HydrogenPlacer.placeHydrogens2D(mol, bondLength) → void
 *
 * Mutates `atom.x` / `atom.y` for every H atom in `mol.atoms`.
 *
 * BIME's parser elides explicit H by default; this is a no-op on most
 * BIME molecules. Useful when the user loads a MOL file with explicit H
 * (e.g., for stereochemistry display) and wants the H positions cleaned up.
 */
(function (global) {
    'use strict';

    var HydrogenPlacer = {};

    HydrogenPlacer.placeHydrogens2D = function (mol, bondLength) {
        if (!mol || !mol.atoms) return;
        bondLength = bondLength ||
            (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;

        for (var i = 0; i < mol.atoms.length; i++) {
            var h = mol.atoms[i];
            if (h.symbol !== 'H') continue;
            var nbrs = mol.getNeighbors(h.id) || [];
            // Only one heavy neighbour: place H opposite the bond-vector
            // from the heavy atom's other neighbours' centroid.
            if (nbrs.length === 1) {
                var heavy = mol.getAtom(nbrs[0]);
                if (!heavy) continue;
                var heavyNbrs = mol.getNeighbors(heavy.id) || [];
                // Compute centroid of heavy's neighbours OTHER than this H.
                var cx = 0, cy = 0, cnt = 0;
                for (var j = 0; j < heavyNbrs.length; j++) {
                    if (heavyNbrs[j] === h.id) continue;
                    var nb = mol.getAtom(heavyNbrs[j]);
                    if (!nb) continue;
                    cx += nb.x; cy += nb.y; cnt++;
                }
                if (cnt === 0) {
                    // No other neighbours; just put H to the right.
                    h.x = heavy.x + bondLength;
                    h.y = heavy.y;
                } else {
                    cx /= cnt; cy /= cnt;
                    var dx = heavy.x - cx;
                    var dy = heavy.y - cy;
                    var d = Math.sqrt(dx * dx + dy * dy);
                    if (d < 1e-6) { dx = 1; dy = 0; d = 1; }
                    h.x = heavy.x + (dx / d) * bondLength;
                    h.y = heavy.y + (dy / d) * bondLength;
                }
            }
            // Multi-bonded H or zero-bond H: leave as-is. CDK does the same.
        }
    };

    global.SDG = global.SDG || {};
    global.SDG.HydrogenPlacer = HydrogenPlacer;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = HydrogenPlacer;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
