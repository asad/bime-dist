/**
 * editor/sdg/AtomPlacer.js — clean-room JS port of CDK's AtomPlacer.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Reference (Java original, ~30 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/AtomPlacer.java
 *
 * The "AtomPlacer" in CDK is a utility class with static methods that
 * place atoms one-at-a-time relative to already-placed neighbours. It's
 * the primary geometry primitive used by StructureDiagramGenerator,
 * RingPlacer, and LayoutRefiner.
 *
 * Public API (matches CDK):
 *   AtomPlacer.distributePartners(atom, placedNeighbours,
 *                                 sharedAtomsCenter, unplacedNeighbours,
 *                                 bondLength) → void
 *   AtomPlacer.populatePolygonCorners(atoms, center, startAngle, addAngle,
 *                                      radius) → void
 *   AtomPlacer.placeLinearChain(atomContainer, initialBondVector,
 *                                bondLength) → void
 *   AtomPlacer.markNotPlaced(atomContainer) → void
 *   AtomPlacer.markPlaced(atomContainer) → void
 *   AtomPlacer.getPlacedAtoms(atomContainer) → array
 *   AtomPlacer.getUnplacedAtoms(atomContainer) → array
 *   AtomPlacer.getNextAtomWithUnplacedNeighbours(atomContainer) → atom
 *   AtomPlacer.shouldBeLinear(atom, atomContainer) → boolean
 *   AtomPlacer.getNextBondVector(atom, previousAtom, distanceMeasure,
 *                                trans) → Vector2d
 *
 * Implementation status:
 *   ✅ distributePartners          — implemented (also exposed in SDGLayout)
 *   ✅ populatePolygonCorners      — implemented
 *   ✅ placeLinearChain            — implemented (deterministic ±60° alternation)
 *   ✅ markPlaced / markNotPlaced  — implemented (uses _placed boolean on atom)
 *   ✅ getPlacedAtoms / getUnplaced — implemented
 *   ✅ getNextAtomWithUnplacedNeighbours — implemented
 *   ⚠ shouldBeLinear              — heuristic stub (CDK uses formal-charge +
 *                                    bond-order + valence; we use degree only)
 *   ⚠ getNextBondVector           — minimal implementation (perpendicular
 *                                    to the previous bond, alternating sign)
 *
 * v1.8.x patches will replace the stub methods with full CDK ports.
 */
(function (global) {
    'use strict';

    var TWO_PI = 2 * Math.PI;
    var DEG60  = Math.PI / 3;
    var EPS = 1e-6;

    function _BL() {
        return (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
    }

    var AtomPlacer = {};

    // ---------------------------------------------------------------------
    // distributePartners — distribute unplaced neighbours in the open
    // angular arc opposite the placed neighbours' centroid.
    // CDK signature includes a `sharedAtomsCenter` (Point2d) for cases
    // where the placed neighbours are part of a fused ring system; we
    // accept it as the 3rd arg but compute the angular range from the
    // placed atoms directly when null.
    // ---------------------------------------------------------------------
    AtomPlacer.distributePartners = function (atom, placedNeighbours,
                                              sharedAtomsCenter,
                                              unplacedNeighbours, bondLength) {
        if (!unplacedNeighbours || unplacedNeighbours.length === 0) return;
        bondLength = bondLength || _BL();

        var n_placed   = placedNeighbours ? placedNeighbours.length : 0;
        var n_unplaced = unplacedNeighbours.length;

        if (n_placed === 0) {
            var step = TWO_PI / n_unplaced;
            var startAng = -Math.PI / 2;
            for (var i = 0; i < n_unplaced; i++) {
                var ang = startAng + i * step;
                unplacedNeighbours[i].x = atom.x + bondLength * Math.cos(ang);
                unplacedNeighbours[i].y = atom.y + bondLength * Math.sin(ang);
            }
            return;
        }

        var placedAngles = new Array(n_placed);
        for (var p = 0; p < n_placed; p++) {
            var dx = placedNeighbours[p].x - atom.x;
            var dy = placedNeighbours[p].y - atom.y;
            var ang2 = Math.atan2(dy, dx);
            if (ang2 < 0) ang2 += TWO_PI;
            placedAngles[p] = ang2;
        }
        placedAngles.sort(function (a, b) { return a - b; });

        var bestGapStart = 0;
        var bestGapWidth = 0;
        for (var k = 0; k < n_placed; k++) {
            var a1 = placedAngles[k];
            var a2 = placedAngles[(k + 1) % n_placed];
            var gap = a2 - a1;
            if (gap < 0) gap += TWO_PI;
            if (gap > bestGapWidth) {
                bestGapWidth = gap;
                bestGapStart = a1;
            }
        }

        var spacing = bestGapWidth / (n_unplaced + 1);
        for (var u = 0; u < n_unplaced; u++) {
            var theta = bestGapStart + spacing * (u + 1);
            unplacedNeighbours[u].x = atom.x + bondLength * Math.cos(theta);
            unplacedNeighbours[u].y = atom.y + bondLength * Math.sin(theta);
        }
    };

    // ---------------------------------------------------------------------
    // populatePolygonCorners — place atoms at the corners of a regular
    // n-gon centred at (cx, cy) with given radius and start angle. CDK's
    // `addAngle` lets you flip the rotation direction; we honour it.
    // ---------------------------------------------------------------------
    AtomPlacer.populatePolygonCorners = function (atoms, center, startAngle,
                                                   addAngle, radius) {
        if (!atoms || atoms.length === 0) return;
        startAngle = startAngle || 0;
        addAngle = addAngle || (TWO_PI / atoms.length);
        if (radius === undefined) radius = _BL();
        for (var i = 0; i < atoms.length; i++) {
            var ang = startAngle + i * addAngle;
            atoms[i].x = center.x + radius * Math.cos(ang);
            atoms[i].y = center.y + radius * Math.sin(ang);
        }
    };

    // ---------------------------------------------------------------------
    // placeLinearChain — 120° zigzag chain placement with deterministic
    // ±60° alternation seeded by the first atom's id parity.
    // ---------------------------------------------------------------------
    AtomPlacer.placeLinearChain = function (chainAtoms, startAtom,
                                             initialAngle, bondLength) {
        if (!chainAtoms || chainAtoms.length === 0) return;
        bondLength = bondLength || _BL();

        var prevX = startAtom.x, prevY = startAtom.y;
        var dir = initialAngle;
        var sign = (startAtom.id & 1) ? +1 : -1;

        for (var i = 0; i < chainAtoms.length; i++) {
            var a = chainAtoms[i];
            a.x = prevX + bondLength * Math.cos(dir);
            a.y = prevY + bondLength * Math.sin(dir);
            prevX = a.x; prevY = a.y;
            dir += sign * DEG60;
            sign = -sign;
        }
    };

    // ---------------------------------------------------------------------
    // markPlaced / markNotPlaced / getPlacedAtoms / getUnplacedAtoms
    // ---------------------------------------------------------------------
    AtomPlacer.markPlaced = function (mol, atomIds) {
        var ids = atomIds || mol.atoms.map(function (a) { return a.id; });
        for (var i = 0; i < ids.length; i++) {
            var a = mol.getAtom(ids[i]);
            if (a) a._placed = true;
        }
    };

    AtomPlacer.markNotPlaced = function (mol, atomIds) {
        var ids = atomIds || mol.atoms.map(function (a) { return a.id; });
        for (var i = 0; i < ids.length; i++) {
            var a = mol.getAtom(ids[i]);
            if (a) a._placed = false;
        }
    };

    AtomPlacer.getPlacedAtoms = function (mol) {
        return mol.atoms.filter(function (a) { return !!a._placed; });
    };

    AtomPlacer.getUnplacedAtoms = function (mol) {
        return mol.atoms.filter(function (a) { return !a._placed; });
    };

    AtomPlacer.getNextAtomWithUnplacedNeighbours = function (mol) {
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            if (!a._placed) continue;
            var nbrs = mol.getNeighbors(a.id) || [];
            for (var j = 0; j < nbrs.length; j++) {
                var n = mol.getAtom(nbrs[j]);
                if (n && !n._placed) return a;
            }
        }
        return null;
    };

    // ---------------------------------------------------------------------
    // shouldBeLinear — heuristic stub. CDK uses formal-charge + bond-order
    // + valence to decide if an atom should sit on a 180° (sp) chain.
    // BIME stub: returns true if every neighbour bond is triple-bonded
    // (acetylene-like). v1.8.x will port the full CDK heuristic.
    // ---------------------------------------------------------------------
    AtomPlacer.shouldBeLinear = function (mol, atom) {
        var nbrs = mol.getNeighbors(atom.id) || [];
        if (nbrs.length !== 2) return false;
        // Check both bonds; if both are triple, force linear.
        for (var i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            if (b.atom1 !== atom.id && b.atom2 !== atom.id) continue;
            if (b.type !== Molecule.BOND_TRIPLE) return false;
        }
        return true;
    };

    // ---------------------------------------------------------------------
    // getNextBondVector — given an atom and its previously-placed neighbour,
    // return a unit vector for the next bond. Minimal implementation:
    // alternate ±60° from the previous direction. Full CDK port pending.
    // ---------------------------------------------------------------------
    AtomPlacer.getNextBondVector = function (atom, previousAtom, trans) {
        var dx = atom.x - previousAtom.x;
        var dy = atom.y - previousAtom.y;
        var len = Math.sqrt(dx * dx + dy * dy) || 1;
        var ux = dx / len, uy = dy / len;
        // Rotate by ±60°.
        var rot = trans ? +DEG60 : -DEG60;
        var c = Math.cos(rot), s = Math.sin(rot);
        return { x: ux * c - uy * s, y: ux * s + uy * c };
    };

    global.SDG = global.SDG || {};
    global.SDG.AtomPlacer = AtomPlacer;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = AtomPlacer;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
