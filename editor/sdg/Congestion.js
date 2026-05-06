/**
 * editor/sdg/Congestion.js — clean-room JS port of CDK's Congestion class.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Reference (Java original):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/Congestion.java
 *
 * Computes a "congestion score" — a measure of how crowded a layout is.
 * Lower is better. Defined as the sum over every pair of non-bonded atoms
 * of (1 / d²) where d is the Euclidean distance between them. Used by
 * LayoutRefiner to compare two candidate layouts (e.g. a 180° rotation)
 * and pick the less-congested one.
 *
 * Public API (matches CDK):
 *   Congestion.score(mol, atomIds, bondLength) → number
 *
 * Determinism: O(N²) pairwise sum, iteration order is by atom-id ascending.
 */
(function (global) {
    'use strict';

    var Congestion = {};

    Congestion.score = function (mol, atomIds, bondLength) {
        if (!mol || !mol.atoms || !atomIds || atomIds.length < 2) return 0;
        bondLength = bondLength ||
            (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
        var minD2 = (bondLength * 0.05) * (bondLength * 0.05); // floor

        // Build bonded set so we exclude bonded pairs (CDK convention).
        var bonded = {};
        for (var i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            bonded[b.atom1 + ',' + b.atom2] = true;
            bonded[b.atom2 + ',' + b.atom1] = true;
        }

        var sum = 0;
        for (var p = 0; p < atomIds.length; p++) {
            var ap = mol.getAtom(atomIds[p]);
            if (!ap) continue;
            for (var q = p + 1; q < atomIds.length; q++) {
                var aq = mol.getAtom(atomIds[q]);
                if (!aq) continue;
                if (bonded[atomIds[p] + ',' + atomIds[q]]) continue;
                var dx = ap.x - aq.x;
                var dy = ap.y - aq.y;
                var d2 = dx * dx + dy * dy;
                if (d2 < minD2) d2 = minD2;
                sum += 1 / d2;
            }
        }
        return sum;
    };

    global.SDG = global.SDG || {};
    global.SDG.Congestion = Congestion;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = Congestion;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
