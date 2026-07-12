/**
 * editor/sdg/HydrogenPlacer.js — native explicit-hydrogen placer.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Places explicit hydrogen atoms on a heavy-atom layout. Computes the angle
 * direction opposite the existing heavy-atom neighbours and places the H
 * along that direction at BOND_LENGTH. For atoms with no heavy neighbours,
 * spreads H atoms around a circle.
 *
 * Public API:
 *   HydrogenPlacer.placeHydrogens2D(mol, bondLength) → void
 *
 * Mutates `atom.x` / `atom.y` ONLY for H atoms that the heavy-atom layout
 * left UNPLACED — i.e. an H still coincident with its heavy neighbour.
 *
 * v1.8.24 fix: previously this function repositioned EVERY H atom
 * unconditionally. BIME's main layout pipeline (Steps 1-15) already
 * lays out hydrogens as ordinary graph atoms, and imported MOL files
 * can arrive with valid explicit-H coordinates — so the unconditional
 * re-placement was throwing away good coordinates AND stacking every
 * H of a CH2/CH3 group along the same direction vector, producing
 * H-H overlaps in ~65 % of the curated benchmark. The guard below
 * restores the documented intent: a no-op on any H that already has
 * a sensible position, active only on genuinely-unplaced hydrogens.
 */
(function (global) {
    'use strict';

    var HydrogenPlacer = {};

    HydrogenPlacer.placeHydrogens2D = function (mol, bondLength) {
        if (!mol || !mol.atoms) return;
        bondLength = bondLength ||
            (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;

        // Track how many H of each heavy atom we have already placed in
        // THIS call, so multiple unplaced H on one heavy atom fan out to
        // distinct angles instead of stacking on the same vector.
        var placedPerHeavy = {};

        for (var i = 0; i < mol.atoms.length; i++) {
            var h = mol.atoms[i];
            if (h.symbol !== 'H') continue;
            var nbrs = mol.getNeighbors(h.id) || [];
            // Only one heavy neighbour: place H opposite the bond-vector
            // from the heavy atom's other neighbours' centroid.
            if (nbrs.length === 1) {
                var heavy = mol.getAtom(nbrs[0]);
                if (!heavy) continue;
                // GUARD: skip H atoms that are already placed. An H that
                // sits a sensible distance from its heavy neighbour was
                // positioned by the main layout pipeline (or imported
                // with valid coordinates); re-placing it discards good
                // data and risks stacking sibling H atoms.
                var hHeavyDist = Math.sqrt(
                    (h.x - heavy.x) * (h.x - heavy.x) +
                    (h.y - heavy.y) * (h.y - heavy.y));
                if (hHeavyDist > bondLength * 0.1) continue;
                var heavyNbrs = mol.getNeighbors(heavy.id) || [];
                // Compute centroid of heavy's neighbours OTHER than this H.
                var cx = 0, cy = 0, cnt = 0;
                for (var j = 0; j < heavyNbrs.length; j++) {
                    if (heavyNbrs[j] === h.id) continue;
                    var nb = mol.getAtom(heavyNbrs[j]);
                    if (!nb) continue;
                    cx += nb.x; cy += nb.y; cnt++;
                }
                // Sibling-index for fan-out: 0 for the first unplaced H
                // of this heavy atom, 1 for the second, etc.
                var sibIdx = placedPerHeavy[heavy.id] || 0;
                placedPerHeavy[heavy.id] = sibIdx + 1;
                // Small alternating angular offset (±25°, ±50° …) so two
                // or three H on the same heavy atom do not stack.
                var fan = 0;
                if (sibIdx > 0) {
                    var step = (Math.PI / 180) * 25;
                    fan = ((sibIdx % 2 === 1) ? 1 : -1) *
                          Math.ceil(sibIdx / 2) * step;
                }
                if (cnt === 0) {
                    // No other neighbours; put H to the right, fanned.
                    h.x = heavy.x + bondLength * Math.cos(fan);
                    h.y = heavy.y + bondLength * Math.sin(fan);
                } else {
                    cx /= cnt; cy /= cnt;
                    var dx = heavy.x - cx;
                    var dy = heavy.y - cy;
                    var d = Math.sqrt(dx * dx + dy * dy);
                    if (d < 1e-6) { dx = 1; dy = 0; d = 1; }
                    var baseAng = Math.atan2(dy / d, dx / d) + fan;
                    h.x = heavy.x + Math.cos(baseAng) * bondLength;
                    h.y = heavy.y + Math.sin(baseAng) * bondLength;
                }
            }
            // Multi-bonded H or zero-bond H: leave as-is.
        }
    };

    global.SDG = global.SDG || {};
    global.SDG.HydrogenPlacer = HydrogenPlacer;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = HydrogenPlacer;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
