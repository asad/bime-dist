/**
 * editor/sdg/LayoutRefiner.js — CDK LayoutRefiner port (v1.8.17 full).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * IP NOTICE — clean-room re-implementation, no source copy:
 *   Independent JavaScript implementation of the algorithm described
 *   in CDK's LayoutRefiner (LGPL-2.1) and Helson '99 SDG. NO CDK
 *   source code copied. The CDK URL is an algorithm reference only.
 *
 * Reference (Java original, ~42 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/LayoutRefiner.java
 *
 * Post-pass refinement that improves a "rough" layout by:
 *   1. alignToLongestAxis(mol)  — PCA rotation so the principal axis
 *      is horizontal (publication convention).
 *   2. rotateRings(mol, rings)  — for each ring, try a 180° flip and
 *      pick whichever orientation has lower Congestion.score.
 *   3. resolveOverlaps(mol)     — delegates to SDG.OverlapResolver.
 *   4. refine(mol, rings, options) — orchestrator that runs 1+2+3.
 *
 * Status (v1.8.17): ALL FOUR FUNCTIONS FULL. The previous v1.8.13
 * scaffold simply returned without acting; this version implements
 * each phase as production-grade JS.
 */
(function (global) {
    'use strict';

    var TWO_PI = 2 * Math.PI;
    var EPS = 1e-9;

    var LayoutRefiner = {};

    /**
     * refine(mol, rings, options) — orchestrator.
     *
     * Runs alignToLongestAxis → rotateRings → resolveOverlaps in that
     * order. Each phase is a no-op if the corresponding option is set
     * to false. Default: all phases run.
     */
    LayoutRefiner.refine = function (mol, rings, options) {
        if (!mol || !mol.atoms) return mol;
        options = options || {};
        if (options.alignToLongestAxis !== false) {
            LayoutRefiner.alignToLongestAxis(mol);
        }
        if (options.rotateRings !== false && rings && rings.length > 0) {
            LayoutRefiner.rotateRings(mol, rings, options);
        }
        if (options.resolveOverlaps !== false) {
            LayoutRefiner.resolveOverlaps(mol, options);
        }
        return mol;
    };

    /**
     * alignToLongestAxis(mol) — rotate the entire molecule about its
     * centroid so the longest principal axis (computed via 2D
     * covariance / PCA) is horizontal. Convention: wider-than-tall.
     *
     * Algorithm:
     *   1. Compute centroid (cx, cy).
     *   2. Build 2D covariance matrix: Sxx = Σ(xᵢ-cx)², Syy = Σ(yᵢ-cy)²,
     *      Sxy = Σ(xᵢ-cx)(yᵢ-cy).
     *   3. Principal angle θ = ½·atan2(2·Sxy, Sxx-Syy).
     *   4. Rotate every atom by -θ about (cx, cy) so the principal
     *      axis aligns with the x-axis.
     *   5. Verify width ≥ height after rotation; if not, rotate by π/2.
     *
     * The rotation preserves all bond lengths, ring shapes, and stereo
     * orientations exactly (it's a rigid transformation).
     *
     * Idempotent: a molecule already aligned produces θ ≈ 0 and the
     * rotation is a no-op.
     */
    LayoutRefiner.alignToLongestAxis = function (mol) {
        if (!mol || !mol.atoms || mol.atoms.length < 2) return;
        var atoms = mol.atoms;
        var n = atoms.length;

        var cx = 0, cy = 0;
        for (var i = 0; i < n; i++) { cx += atoms[i].x; cy += atoms[i].y; }
        cx /= n; cy /= n;

        var sxx = 0, syy = 0, sxy = 0;
        for (var i2 = 0; i2 < n; i2++) {
            var dx = atoms[i2].x - cx, dy = atoms[i2].y - cy;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }
        // Skip if degenerate (all atoms co-located).
        if (sxx + syy < EPS) return;
        var theta = 0.5 * Math.atan2(2 * sxy, sxx - syy);
        // Rotate all atoms by -theta about (cx, cy).
        var ct = Math.cos(-theta), st = Math.sin(-theta);
        for (var i3 = 0; i3 < n; i3++) {
            var ddx = atoms[i3].x - cx, ddy = atoms[i3].y - cy;
            atoms[i3].x = cx + ddx * ct - ddy * st;
            atoms[i3].y = cy + ddx * st + ddy * ct;
        }
        // Verify width ≥ height; if not (i.e., principal axis was
        // aligned to vertical), rotate by π/2.
        var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        for (var i4 = 0; i4 < n; i4++) {
            if (atoms[i4].x < minX) minX = atoms[i4].x;
            if (atoms[i4].x > maxX) maxX = atoms[i4].x;
            if (atoms[i4].y < minY) minY = atoms[i4].y;
            if (atoms[i4].y > maxY) maxY = atoms[i4].y;
        }
        var w = maxX - minX, h = maxY - minY;
        if (h > w) {
            var ct2 = Math.cos(-Math.PI / 2), st2 = Math.sin(-Math.PI / 2);
            for (var i5 = 0; i5 < n; i5++) {
                var ddx2 = atoms[i5].x - cx, ddy2 = atoms[i5].y - cy;
                atoms[i5].x = cx + ddx2 * ct2 - ddy2 * st2;
                atoms[i5].y = cy + ddx2 * st2 + ddy2 * ct2;
            }
        }
    };

    /**
     * rotateRings(mol, rings, options) — for each ring, try rotating
     * it 180° about its centroid and keep whichever orientation has
     * lower Congestion.score (computed locally over the ring's atoms
     * + their non-ring neighbours).
     *
     * Use case: a ring with bulky substituents on top might depict
     * better with substituents on the bottom. CDK uses this for
     * publication-style depictions.
     *
     * Conservative: only flips a ring if Congestion improves by ≥ 5%.
     * Skips rings shared between systems (fused or spiro) since
     * flipping one would deform the rest.
     */
    LayoutRefiner.rotateRings = function (mol, rings, options) {
        if (!mol || !rings || rings.length === 0) return;
        options = options || {};
        var IMPROVEMENT_THRESHOLD = options.improvementThreshold || 0.05;
        var Congestion = global.SDG && global.SDG.Congestion;
        if (!Congestion || !Congestion.score) return;

        // Build atom→ring-count map; only flip rings whose every atom
        // is in exactly THIS ring (i.e., not part of a fused system).
        var atomRingCount = {};
        for (var ri = 0; ri < rings.length; ri++) {
            var r = rings[ri].atoms || rings[ri];
            for (var ai = 0; ai < r.length; ai++) {
                atomRingCount[r[ai]] = (atomRingCount[r[ai]] || 0) + 1;
            }
        }

        for (var ri2 = 0; ri2 < rings.length; ri2++) {
            var ring = rings[ri2].atoms || rings[ri2];
            // Only flip if no atom is shared with another ring.
            var anyShared = false;
            for (var ai2 = 0; ai2 < ring.length; ai2++) {
                if (atomRingCount[ring[ai2]] > 1) { anyShared = true; break; }
            }
            if (anyShared) continue;

            // Compute ring centroid.
            var cx = 0, cy = 0, cn = 0;
            for (var ai3 = 0; ai3 < ring.length; ai3++) {
                var a = mol.getAtom(ring[ai3]);
                if (a) { cx += a.x; cy += a.y; cn++; }
            }
            if (cn === 0) continue;
            cx /= cn; cy /= cn;

            // Collect ring atoms + their non-ring neighbours (substituent
            // atoms) — Congestion.score over this set.
            var localIds = ring.slice();
            var localSet = {};
            for (var ai4 = 0; ai4 < ring.length; ai4++) localSet[ring[ai4]] = true;
            for (var ai5 = 0; ai5 < ring.length; ai5++) {
                var nbrs = mol.getNeighbors(ring[ai5]) || [];
                for (var nb = 0; nb < nbrs.length; nb++) {
                    if (!localSet[nbrs[nb]]) {
                        localSet[nbrs[nb]] = true;
                        localIds.push(nbrs[nb]);
                    }
                }
            }

            // Score the current orientation.
            var scoreA = Congestion.score(mol, localIds);

            // Snapshot ring atom positions, flip 180° about (cx, cy),
            // score, and either keep the flip or restore.
            var snap = {};
            for (var ai6 = 0; ai6 < ring.length; ai6++) {
                var aa = mol.getAtom(ring[ai6]);
                if (aa) snap[ring[ai6]] = { x: aa.x, y: aa.y };
            }
            // Also flip substituents that hang off ring atoms (so the
            // local geometry stays connected).
            var subSnap = {};
            var subIds = [];
            for (var li = 0; li < localIds.length; li++) {
                if (!{}.hasOwnProperty.call(snap, localIds[li])) {
                    var sa = mol.getAtom(localIds[li]);
                    if (sa) {
                        subSnap[localIds[li]] = { x: sa.x, y: sa.y };
                        subIds.push(localIds[li]);
                    }
                }
            }
            for (var ai7 = 0; ai7 < ring.length; ai7++) {
                var ab = mol.getAtom(ring[ai7]);
                if (!ab) continue;
                ab.x = 2 * cx - ab.x;
                ab.y = 2 * cy - ab.y;
            }
            for (var si = 0; si < subIds.length; si++) {
                var sb = mol.getAtom(subIds[si]);
                if (!sb) continue;
                sb.x = 2 * cx - sb.x;
                sb.y = 2 * cy - sb.y;
            }

            var scoreB = Congestion.score(mol, localIds);

            // Keep the flip iff it improves congestion by the threshold.
            if (scoreB < scoreA * (1 - IMPROVEMENT_THRESHOLD)) {
                // Keep flip — leave coordinates as-is.
            } else {
                // Restore.
                for (var k in snap) {
                    var ka = mol.getAtom(parseInt(k, 10));
                    if (ka) { ka.x = snap[k].x; ka.y = snap[k].y; }
                }
                for (var k2 in subSnap) {
                    var kb = mol.getAtom(parseInt(k2, 10));
                    if (kb) { kb.x = subSnap[k2].x; kb.y = subSnap[k2].y; }
                }
            }
        }
    };

    /**
     * resolveOverlaps(mol, options) — delegates to SDG.OverlapResolver.
     */
    LayoutRefiner.resolveOverlaps = function (mol, options) {
        if (!global.SDG || !global.SDG.OverlapResolver) return 0;
        if (!mol || !mol.atoms) return 0;
        var ids = mol.atoms.map(function (a) { return a.id; });
        return global.SDG.OverlapResolver.resolveOverlap(mol, ids, options || {});
    };

    global.SDG = global.SDG || {};
    global.SDG.LayoutRefiner = LayoutRefiner;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = LayoutRefiner;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
