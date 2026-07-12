/**
 * SDGLayout.js — v1.8.12 Structure Diagram Generator phases
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Native BIME diagram-generation phases:
 *
 * Phase mapping:
 *
 *   AtomPlacer.placeLinearChain           → SDGLayout.placeLinearChain
 *   AtomPlacer.distributePartners         → SDGLayout.distributePartners
 *   AtomPlacer.populatePolygonCorners     → SDGLayout.populatePolygonCorners
 *   AtomPlacer.markNotPlaced / markPlaced → tracked via a per-atom 'placed' flag
 *   OverlapResolver.resolveOverlap         → SDGLayout.resolveOverlap
 *       (v1.8.25: now a thin wrapper that delegates to the single
 *        canonical implementation in editor/sdg/OverlapResolver.js)
 *   StructureDiagramGenerator.generateCoordinates(IAtomContainer, boolean firstBondVector, boolean isConnected)
 *                                          → SDGLayout.layout(mol)
 *
 * BIME's existing editor/Layout.js implements similar phases via a different
 * code path. SDGLayout is a corrective pass for cases where Layout.js
 * produces collapsed rings, overlapping ring systems, or chain-tangle
 * (the "broken depict" reports in v1.8.x).
 *
 * Activate via Layout.options.useSDG = true. When active, after the
 * existing Layout.js pipeline completes its rough placement, SDGLayout
 * re-walks the placed atoms applying:
 *   1. distributePartners on every atom with un-placed neighbours
 *   2. placeLinearChain for any chain that wandered off-direction
 *   3. resolveOverlap with summed forces until stable
 *
 * Wired into editor/Layout.js step 14 (gated by Layout.options.useSDGRefine,
 * default ON since the regression suite baselined favourably).
 *
 * Determinism: every iteration order is by atom-id (stable, no reliance
 * on hash-table iteration order), so same input → byte-identical output.
 */
(function (global) {
    'use strict';

    var BL = (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
    var TWO_PI = 2 * Math.PI;
    var DEG60  = Math.PI / 3;
    var DEG120 = 2 * Math.PI / 3;

    var SDGLayout = {};

    // ---------------------------------------------------------------------
    // distributePartners
    //
    // Given an atom A whose placed neighbours occupy a known angular range,
    // distribute its un-placed neighbours evenly within the *open arc* on
    // the opposite side. This is the canonical algorithm that makes
    // substituents on a ring point radially outward.
    //
    // Inputs:
    //   atom            : the focal atom (already placed)
    //   placedNeighbours: array of Atom objects (already placed, give us
    //                     the angular range to AVOID)
    //   unplacedNeighbours: array of Atom objects to place around `atom`
    //   bondLength      : target bond length (defaults to BOND_LENGTH)
    //
    // Output: each `unplacedNeighbours[i]` gets its (x, y) set; placed flag
    // is set on the SDGLayout-side ledger.
    // ---------------------------------------------------------------------
    SDGLayout.distributePartners = function (atom, placedNeighbours,
                                             unplacedNeighbours, bondLength) {
        if (!unplacedNeighbours || unplacedNeighbours.length === 0) return;
        bondLength = bondLength || BL;

        var n_placed   = placedNeighbours ? placedNeighbours.length : 0;
        var n_unplaced = unplacedNeighbours.length;

        if (n_placed === 0) {
            // No existing neighbours — distribute all around full circle.
            var step = TWO_PI / n_unplaced;
            var startAng = -Math.PI / 2; // start at top
            for (var i = 0; i < n_unplaced; i++) {
                var ang = startAng + i * step;
                unplacedNeighbours[i].x = atom.x + bondLength * Math.cos(ang);
                unplacedNeighbours[i].y = atom.y + bondLength * Math.sin(ang);
            }
            return;
        }

        // Compute angle of each placed neighbour relative to the focal atom.
        // Sort ascending in [0, 2π).
        var placedAngles = new Array(n_placed);
        for (var p = 0; p < n_placed; p++) {
            var dx = placedNeighbours[p].x - atom.x;
            var dy = placedNeighbours[p].y - atom.y;
            var ang2 = Math.atan2(dy, dx);
            if (ang2 < 0) ang2 += TWO_PI;
            placedAngles[p] = ang2;
        }
        placedAngles.sort(function (a, b) { return a - b; });

        // Find the largest angular gap between consecutive placed angles
        // (wrapping around). The unplaced neighbours go in this gap.
        var bestGapStart = 0;
        var bestGapWidth = 0;
        for (var k = 0; k < n_placed; k++) {
            var a1 = placedAngles[k];
            var a2 = placedAngles[(k + 1) % n_placed];
            var gap = a2 - a1;
            if (gap < 0) gap += TWO_PI;          // wrap
            if (gap > bestGapWidth) {
                bestGapWidth = gap;
                bestGapStart = a1;
            }
        }

        // Distribute n_unplaced atoms inside the gap (excluding endpoints
        // — those are the placed neighbours' angles).
        // Spacing = gap / (n_unplaced + 1) for even distribution.
        var spacing = bestGapWidth / (n_unplaced + 1);
        for (var u = 0; u < n_unplaced; u++) {
            var theta = bestGapStart + spacing * (u + 1);
            unplacedNeighbours[u].x = atom.x + bondLength * Math.cos(theta);
            unplacedNeighbours[u].y = atom.y + bondLength * Math.sin(theta);
        }
    };

    // ---------------------------------------------------------------------
    // placeLinearChain
    //
    // Place a linear chain of atoms in a 120° zigzag. The first atom is
    // already placed (`startAtom`); we walk along `chainAtoms` placing
    // each next atom with bond-length BL and an angle alternating ±60°
    // from the previous bond direction.
    //
    // The first bond direction is given by `initialAngle` (radians).
    // Subsequent bonds rotate the direction by ±60° to produce the
    // 120° interior angle that's standard for sp³ chain depictions.
    //
    // The alternation sign is chosen deterministically — atom-id parity
    // — so re-runs on the same input produce byte-identical coordinates.
    // ---------------------------------------------------------------------
    SDGLayout.placeLinearChain = function (chainAtoms, startAtom, initialAngle,
                                            bondLength) {
        if (!chainAtoms || chainAtoms.length === 0) return;
        bondLength = bondLength || BL;

        var prevX = startAtom.x;
        var prevY = startAtom.y;
        var dir = initialAngle;
        // Sign alternates ±60° to produce 120° interior angles.
        // Initial sign chosen by parity of startAtom.id for determinism.
        var sign = (startAtom.id & 1) ? +1 : -1;

        for (var i = 0; i < chainAtoms.length; i++) {
            var a = chainAtoms[i];
            a.x = prevX + bondLength * Math.cos(dir);
            a.y = prevY + bondLength * Math.sin(dir);
            prevX = a.x;
            prevY = a.y;
            dir += sign * DEG60;
            sign = -sign;
        }
    };

    // ---------------------------------------------------------------------
    // populatePolygonCorners
    //
    // Place atoms at the corners of a regular polygon centred on (cx, cy)
    // with given radius and starting angle. Used to lay out a ring as a
    // regular n-gon. BIME's Layout.js has placeRingAsPolygon which does
    // the same; this is the SDG corrective-pass version.
    // ---------------------------------------------------------------------
    SDGLayout.populatePolygonCorners = function (ringAtoms, cx, cy, radius,
                                                  startAngle) {
        if (!ringAtoms || ringAtoms.length === 0) return;
        startAngle = startAngle || 0;
        var n = ringAtoms.length;
        var step = TWO_PI / n;
        for (var i = 0; i < n; i++) {
            var ang = startAngle + i * step;
            ringAtoms[i].x = cx + radius * Math.cos(ang);
            ringAtoms[i].y = cy + radius * Math.sin(ang);
        }
    };

    // ---------------------------------------------------------------------
    // resolveOverlap — delegates to editor/sdg/OverlapResolver.js, the
    // single source of truth for the summed-force algorithm (consolidated
    // in v1.8.25). Falls through silently if the SDG module is absent
    // (an undocumented bare-file load).
    // ---------------------------------------------------------------------
    SDGLayout.resolveOverlap = function (mol, atomIds, opts) {
        if (global.SDG && global.SDG.OverlapResolver &&
            global.SDG.OverlapResolver.resolveOverlap) {
            return global.SDG.OverlapResolver.resolveOverlap(mol, atomIds, opts);
        }
        return 0;
    };

    // ---------------------------------------------------------------------
    // refineLayout — entry point for the v1.8.12 corrective pass.
    //
    // Run AFTER editor/Layout.js's main pipeline. Applies, in order:
    //   1. resolveOverlap with summed forces
    //   2. distributePartners on atoms whose substituent fan looks bad
    //
    // Idempotent: re-running converges to the same fixed point given the
    // same input.
    // ---------------------------------------------------------------------
    SDGLayout.refineLayout = function (mol, atomIds, opts) {
        if (!mol || !mol.atoms || mol.atoms.length < 2) return;
        opts = opts || {};
        atomIds = atomIds || mol.atoms.map(function (a) { return a.id; });

        // Phase 1: summed-force overlap resolution.
        SDGLayout.resolveOverlap(mol, atomIds, opts);

        // Phase 2: substituent fan repair. For every atom with degree ≥ 3,
        // check if any pair of bonded neighbours has angular separation
        // < 30°. If so, redistribute via distributePartners.
        // (Conservative: only fix atoms where the fan is visibly cramped.)
        for (var i = 0; i < atomIds.length; i++) {
            var a = mol.getAtom(atomIds[i]);
            if (!a) continue;
            var nbrs = mol.getNeighbors(a.id);
            if (nbrs.length < 3) continue;

            // Compute angles of each neighbour relative to a.
            var angs = nbrs.map(function (id) {
                var nb = mol.getAtom(id);
                return { id: id, angle: Math.atan2(nb.y - a.y, nb.x - a.x) };
            });
            angs.sort(function (x, y) { return x.angle - y.angle; });

            // Find min angular gap.
            var minGap = TWO_PI;
            for (var g = 0; g < angs.length; g++) {
                var gap = angs[(g + 1) % angs.length].angle - angs[g].angle;
                if (gap < 0) gap += TWO_PI;
                if (gap < minGap) minGap = gap;
            }

            // If any pair is < 30° apart, the fan is cramped. Redistribute
            // — but only the neighbours not in a ring (so we don't shift
            // ring atoms which have positional constraints).
            if (minGap < Math.PI / 6) {
                var ringSet = SDGLayout._buildInRingSet(mol);
                var inRing = [], outRing = [];
                for (var nn = 0; nn < nbrs.length; nn++) {
                    var nbAtom = mol.getAtom(nbrs[nn]);
                    if (ringSet[nbrs[nn]]) inRing.push(nbAtom);
                    else outRing.push(nbAtom);
                }
                // Redistribute the non-ring neighbours in the open arc.
                if (outRing.length > 0 && inRing.length > 0) {
                    SDGLayout.distributePartners(a, inRing, outRing, BL);
                }
            }
        }

        // Phase 3: final overlap pass to clean up after redistribution.
        SDGLayout.resolveOverlap(mol, atomIds,
            Object.assign({}, opts, { maxIters: 30 }));
    };

    // ---------------------------------------------------------------------
    // Internal: build a set of atom-ids that are in any ring up to size 8.
    // ---------------------------------------------------------------------
    SDGLayout._buildInRingSet = function (mol) {
        var set = {};
        if (!mol.findRings) return set;
        var rings = mol.findRings(8) || [];
        for (var i = 0; i < rings.length; i++) {
            var atoms = rings[i].atoms || rings[i];
            for (var j = 0; j < atoms.length; j++) set[atoms[j]] = true;
        }
        return set;
    };

    // ---------------------------------------------------------------------
    // Export.
    // ---------------------------------------------------------------------
    global.SDGLayout = SDGLayout;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = SDGLayout;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
