/**
 * editor/sdg/RingPlacer.js — CDK RingPlacer port (v1.8.17 full).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * IP NOTICE — clean-room re-implementation, no source copy:
 *   This file (and the v1.8.15 fused-ring snap helpers in
 *   editor/Layout.js — `snapRingAtCentroid15`, `snapFusedRingByReflection15`,
 *   `snapSpiroRingFromAtom15`) is an independent JavaScript implementation
 *   of the algorithm described in CDK's RingPlacer (LGPL-2.1) and
 *   Helson, T.; "Structure Diagram Generation"; Reviews in Computational
 *   Chemistry vol. 13 (1999). NO CDK source code has been copied,
 *   transliterated, or machine-translated into BIME. The CDK URL in this
 *   file is a reference to the upstream public API + algorithm
 *   description, not a code source. Apache-2.0 licensing is preserved.
 *
 * Reference (Java original, ~42 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/RingPlacer.java
 *
 * RingPlacer handles ring-system placement: regular polygons, fused
 * rings (sharing 2 atoms), spiro junctions (sharing 1 atom), bridged
 * ring templates (≥ 3 shared atoms), and macrocycles. v1.8.17 fills
 * in the bridged-ring + macrocycle paths; the fused / spiro paths are
 * implemented in editor/Layout.js Step 15 helpers (see above).
 */
(function (global) {
    'use strict';

    function _BL() {
        return (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
    }

    var TWO_PI = 2 * Math.PI;
    var EPS = 1e-9;

    var RingPlacer = {};

    /**
     * placeRing(ring, sharedAtoms, sharedAtomsCenter, ringCenterVector, bondLength)
     * — top-level ring placement entry. Dispatches based on the
     * sharing topology:
     *   sharedAtoms.length === 0 → regular polygon at origin
     *   sharedAtoms.length === 1 → spiro (delegates to placeSpiroRing)
     *   sharedAtoms.length === 2 → edge-fused (delegates to placeFusedRing)
     *   sharedAtoms.length >= 3  → bridged (delegates to placeBridgedRing)
     */
    RingPlacer.placeRing = function (ring, sharedAtoms, sharedAtomsCenter,
                                      ringCenterVector, bondLength) {
        if (!ring) return;
        bondLength = bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 3) return;

        if (!sharedAtoms || sharedAtoms.length === 0) {
            return RingPlacer.placeRegular(ring, 0, 0, -Math.PI / 2, bondLength);
        }
        if (sharedAtoms.length === 1) {
            return RingPlacer.placeSpiroRing(ring, sharedAtoms[0], bondLength);
        }
        if (sharedAtoms.length === 2) {
            return RingPlacer.placeFusedRing(ring, sharedAtoms, bondLength);
        }
        return RingPlacer.placeBridgedRing(ring, sharedAtoms, bondLength);
    };

    /**
     * placeRegular(ring, cx, cy, startAng, bondLength) — place ring as
     * a regular n-gon centred at (cx, cy) with vertex 0 at angle
     * startAng from the centre.
     */
    RingPlacer.placeRegular = function (ring, cx, cy, startAng, bondLength) {
        bondLength = bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 3) return;
        var n = atoms.length;
        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        for (var i = 0; i < n; i++) {
            var ang = startAng + i * TWO_PI / n;
            atoms[i].x = cx + radius * Math.cos(ang);
            atoms[i].y = cy + radius * Math.sin(ang);
        }
    };

    /**
     * placeFusedRing(ring, sharedAtomIds, bondLength) — place a ring
     * that shares EXACTLY 2 contiguous atoms with an already-placed
     * neighbour. Uses edge reflection: the new ring's centre is on
     * the opposite side of the shared edge from the parent, at apothem
     * distance.
     *
     * The two shared atoms must already have positions; the rest of
     * the ring is placed by regular-polygon vertex assignment around
     * the reflected centre.
     *
     * NOTE: this is a STANDALONE invocation — for the in-pipeline
     * Step 15 fused-ring snap, see editor/Layout.js
     * snapFusedRingByReflection15.
     */
    RingPlacer.placeFusedRing = function (ring, sharedAtomIds, bondLength) {
        bondLength = bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 3) return;
        if (!sharedAtomIds || sharedAtomIds.length < 2) return;
        var n = atoms.length;

        // Find the shared atoms' indices in the ring traversal.
        var idx1 = -1, idx2 = -1;
        for (var i = 0; i < n; i++) {
            var aid = atoms[i].id !== undefined ? atoms[i].id : atoms[i];
            if (aid === sharedAtomIds[0]) idx1 = i;
            else if (aid === sharedAtomIds[1]) idx2 = i;
        }
        if (idx1 < 0 || idx2 < 0) return;
        var a1 = atoms[idx1], a2 = atoms[idx2];

        // Compute reflected centre.
        var mx = (a1.x + a2.x) / 2, my = (a1.y + a2.y) / 2;
        var edx = a2.x - a1.x, edy = a2.y - a1.y;
        var nrm = Math.hypot(edx, edy);
        if (nrm < EPS) return;
        var nx = -edy / nrm, ny = edx / nrm;
        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        var apothem = radius * Math.cos(Math.PI / n);
        var cx = mx + nx * apothem, cy = my + ny * apothem;

        // Place each non-shared atom at its polygon vertex.
        var ang1 = Math.atan2(a1.y - cy, a1.x - cx);
        var ang2 = Math.atan2(a2.y - cy, a2.x - cx);
        var step = TWO_PI / n;
        var stepsForward = (idx2 - idx1 + n) % n;
        var angleDiff = ang2 - ang1;
        while (angleDiff > Math.PI) angleDiff -= TWO_PI;
        while (angleDiff < -Math.PI) angleDiff += TWO_PI;
        var clockwise = (Math.abs(angleDiff - stepsForward * step) <
                         Math.abs(angleDiff + stepsForward * step));
        for (var i2 = 0; i2 < n; i2++) {
            if (i2 === idx1 || i2 === idx2) continue;  // shared atoms keep position
            var steps = (i2 - idx1 + n) % n;
            var ang = clockwise ? ang1 + steps * step : ang1 - steps * step;
            atoms[i2].x = cx + radius * Math.cos(ang);
            atoms[i2].y = cy + radius * Math.sin(ang);
        }
    };

    /**
     * placeSpiroRing(ring, sharedAtomId, bondLength) — place a ring
     * sharing a single atom (spiro junction). The new ring is rotated
     * 90° relative to the existing ring's orientation at the spiro
     * atom (axis perpendicular to the bisector of the existing bonds).
     *
     * NOTE: standalone variant. For the in-pipeline version, see
     * editor/Layout.js snapSpiroRingFromAtom15.
     */
    RingPlacer.placeSpiroRing = function (ring, sharedAtomId, bondLength) {
        bondLength = bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 3) return;
        var n = atoms.length;
        var spiroIdx = -1;
        for (var i = 0; i < n; i++) {
            var aid = atoms[i].id !== undefined ? atoms[i].id : atoms[i];
            if (aid === sharedAtomId) { spiroIdx = i; break; }
        }
        if (spiroIdx < 0) return;
        var spiroAtom = atoms[spiroIdx];

        var radius = bondLength / (2 * Math.sin(Math.PI / n));
        // Default: axis at 0° (no parent reference available in this
        // standalone API). Caller is responsible for orienting if
        // needed; the in-pipeline variant uses the spiro atom's
        // existing bonds to compute the bisector.
        var awayAng = 0;
        var cx = spiroAtom.x + radius * Math.cos(awayAng);
        var cy = spiroAtom.y + radius * Math.sin(awayAng);
        var startAng = Math.atan2(spiroAtom.y - cy, spiroAtom.x - cx);
        var step = TWO_PI / n;
        for (var i2 = 0; i2 < n; i2++) {
            if (i2 === spiroIdx) continue;
            var steps = (i2 - spiroIdx + n) % n;
            var ang = startAng + steps * step;
            atoms[i2].x = cx + radius * Math.cos(ang);
            atoms[i2].y = cy + radius * Math.sin(ang);
        }
    };

    /**
     * placeBridgedRing(ring, sharedAtomIds, bondLength) — handle a ring
     * that shares ≥ 3 atoms with an already-placed neighbour. Uses the
     * "two bridgeheads + bridge axis" geometry that BIME's existing
     * `fuseBridgedRing` (Layout.js) implements but extracted here so
     * external SDG callers can reach it directly.
     *
     * Algorithm:
     *   1. Find the two BRIDGEHEAD atoms = the shared-atom pair with
     *      the maximum ring-walk distance.
     *   2. Compute the perpendicular axis to the bridgehead axis.
     *   3. Place each non-shared atom along the perpendicular axis on
     *      the side opposite the existing bridge.
     *
     * For norbornane (5+5 sharing 3 atoms), this places the bridge
     * methylene above the cyclopentane plane.
     */
    RingPlacer.placeBridgedRing = function (ring, sharedAtomIds, bondLength) {
        bondLength = bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 4) return;
        if (!sharedAtomIds || sharedAtomIds.length < 3) return;
        var n = atoms.length;
        var sharedSet = {};
        for (var s = 0; s < sharedAtomIds.length; s++) sharedSet[sharedAtomIds[s]] = true;
        // Map shared atoms to ring indices.
        var sharedIdx = [];
        for (var i = 0; i < n; i++) {
            var aid = atoms[i].id !== undefined ? atoms[i].id : atoms[i];
            if (sharedSet[aid]) sharedIdx.push(i);
        }
        if (sharedIdx.length < 3) return;

        // Find the two bridgeheads = pair of shared indices with max
        // ring-walk distance (i.e., farthest apart in the cycle).
        var bh1 = sharedIdx[0], bh2 = sharedIdx[sharedIdx.length - 1];
        var maxDist = 0;
        for (var i2 = 0; i2 < sharedIdx.length; i2++) {
            for (var j2 = i2 + 1; j2 < sharedIdx.length; j2++) {
                var d = Math.min((sharedIdx[j2] - sharedIdx[i2] + n) % n,
                                 (sharedIdx[i2] - sharedIdx[j2] + n) % n);
                if (d > maxDist) {
                    maxDist = d;
                    bh1 = sharedIdx[i2];
                    bh2 = sharedIdx[j2];
                }
            }
        }
        var bha = atoms[bh1], bhb = atoms[bh2];

        // Bridge axis (midpoint, perpendicular).
        var mx = (bha.x + bhb.x) / 2, my = (bha.y + bhb.y) / 2;
        var edx = bhb.x - bha.x, edy = bhb.y - bha.y;
        var edLen = Math.hypot(edx, edy);
        if (edLen < EPS) return;
        var nx = -edy / edLen, ny = edx / edLen;

        // Determine which side of the bridge axis is less crowded.
        // "Crowding" on each side = sum of distances to placed atoms
        // dotted with the normal.
        var sumPos = 0, sumNeg = 0;
        for (var k = 0; k < n; k++) {
            if (k === bh1 || k === bh2) continue;
            if (sharedSet[atoms[k].id !== undefined ? atoms[k].id : atoms[k]]) {
                var dot = (atoms[k].x - mx) * nx + (atoms[k].y - my) * ny;
                if (dot > 0) sumPos++;
                else sumNeg++;
            }
        }
        // Less-crowded side = where to put NEW atoms.
        var sign = (sumPos <= sumNeg) ? 1 : -1;

        // Walk forward from bh1 to bh2 placing non-shared atoms along
        // an arc on the less-crowded side.
        var stepsForward = (bh2 - bh1 + n) % n;
        var stepsBackward = (bh1 - bh2 + n) % n;
        // Pick the shorter "bridge path" (fewer atoms to place along
        // the perpendicular axis).
        var bridgePath = [];
        var pathStart = bh1, pathLen = stepsForward;
        if (stepsBackward < stepsForward) {
            pathStart = bh2;
            pathLen = stepsBackward;
        }
        for (var p = 1; p < pathLen; p++) {
            var idx = (pathStart + p) % n;
            if (sharedSet[atoms[idx].id !== undefined ? atoms[idx].id : atoms[idx]]) continue;
            bridgePath.push(idx);
        }
        if (bridgePath.length === 0) return;

        // Place the bridge path along an arc perpendicular to the
        // bridgehead axis on the chosen side.
        var arcRadius = (bondLength * (bridgePath.length + 1)) / 2;
        for (var bp = 0; bp < bridgePath.length; bp++) {
            var t = (bp + 1) / (bridgePath.length + 1);
            // Position along the bridgehead axis, then offset by arc.
            var alongX = bha.x + (bhb.x - bha.x) * t;
            var alongY = bha.y + (bhb.y - bha.y) * t;
            // Arc-height proportional to position (bell shape).
            var arcH = arcRadius * Math.sin(Math.PI * t);
            atoms[bridgePath[bp]].x = alongX + sign * nx * arcH;
            atoms[bridgePath[bp]].y = alongY + sign * ny * arcH;
        }
    };

    /**
     * placeRingSubstituents(mol, ringAtoms, bondLength) — distribute
     * non-ring substituents on each ring atom evenly into the open
     * angular arc opposite the ring interior.
     *
     * Delegates to SDG.AtomPlacer.distributePartners if available.
     */
    RingPlacer.placeRingSubstituents = function (mol, ringAtoms, bondLength) {
        bondLength = bondLength || _BL();
        var AtomPlacer = global.SDG && global.SDG.AtomPlacer;
        if (!AtomPlacer || !AtomPlacer.distributePartners) return;
        // For each ring atom, find non-ring neighbours and distribute.
        var ringSet = {};
        for (var i = 0; i < ringAtoms.length; i++) ringSet[ringAtoms[i]] = true;
        for (var i2 = 0; i2 < ringAtoms.length; i2++) {
            var nbrs = mol.getNeighbors(ringAtoms[i2]) || [];
            var subs = [];
            for (var j = 0; j < nbrs.length; j++) {
                if (!ringSet[nbrs[j]]) subs.push(nbrs[j]);
            }
            if (subs.length === 0) continue;
            AtomPlacer.distributePartners(mol, ringAtoms[i2], subs, bondLength);
        }
    };

    /**
     * placeMacroRing(ring, bondLength) — delegate to MacroCycleLayout
     * for ≥ 8-membered rings (egg-shape ovalisation).
     */
    RingPlacer.placeMacroRing = function (ring, bondLength) {
        var MacroCycleLayout = global.SDG && global.SDG.MacroCycleLayout;
        if (!MacroCycleLayout || !MacroCycleLayout.layout) {
            // Fallback to regular polygon.
            return RingPlacer.placeRegular(ring, 0, 0, -Math.PI / 2, bondLength);
        }
        return MacroCycleLayout.layout(ring, { bondLength: bondLength });
    };

    global.SDG = global.SDG || {};
    global.SDG.RingPlacer = RingPlacer;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = RingPlacer;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
