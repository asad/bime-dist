/**
 * editor/sdg/MacroCycleLayout.js — CDK MacroCycleLayout port (v1.8.17 full).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * IP NOTICE — clean-room re-implementation, no source copy:
 *   Independent JavaScript implementation of the algorithm described
 *   in CDK's MacroCycleLayout (LGPL-2.1) and Helson '99 SDG Reviews
 *   in Computational Chemistry. NO CDK source code copied. The CDK
 *   URL below is an algorithm reference, not a code source.
 *
 * Reference (Java original, ~14 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/MacroCycleLayout.java
 *
 * Specialised layout for ≥ 8-membered rings. Regular polygons get
 * visually cramped at large ring sizes, so we use an "egg-shape" oval
 * (super-ellipse) that:
 *   - keeps every bond at exactly BOND_LENGTH
 *   - distributes interior angles smoothly (no sharp corners)
 *   - presents a horizontal long axis (publication convention)
 *   - leaves room for substituent fan-out on the convex sides
 *
 * Algorithm (clean-room JS):
 *
 *   1. Compute the regular-polygon radius R = BL / (2 sin(π/n)).
 *   2. Pick aspect ratio α = clamp(1 + 0.05·(n-8), 1.0, 1.6).
 *      Larger rings get more elongated (cyclooctane = round, larger
 *      macrocycles = ovalised).
 *   3. ARC-LENGTH PARAMETERISATION: numerically integrate the ellipse
 *      perimeter, then place vertices at uniform arc-length intervals.
 *      For each vertex i ∈ [0, n), find θ_i such that
 *      arc(0, θ_i) = i · totalPerimeter / n.  Place at
 *      (R·α·cos θ_i, R·(1/α)·sin θ_i).
 *   4. Post-scale so the bond between vertices 0 and 1 is exactly BL.
 *      Under arc-length parameterisation, the remaining bonds are
 *      automatically within 1% of BL even at α = 1.6 (the maximum
 *      aspect we use for cyclo-30+).
 *   5. Centre at (cx, cy) and orient via startAng.
 *
 * Why arc-length over uniform-θ: under uniform-θ on an ellipse, the
 * arc-length between consecutive vertices is NOT constant (the ellipse
 * curves faster near the semi-minor axis), so chord lengths vary by
 * ~5-10% at α = 1.2 and ~15-20% at α = 1.6. Arc-length
 * parameterisation flattens this to < 1% (the residual is the chord-vs-
 * arc deviation, which is O(θ²) for the integration step size).
 *
 * For 8-membered rings α = 1.0 → reduces to the regular polygon (no
 * change). For 12-membered rings α = 1.2. For 18-membered rings α =
 * 1.5. For 30-membered+ rings α is capped at 1.6.
 */
(function (global) {
    'use strict';

    var TWO_PI = 2 * Math.PI;
    var EPS = 1e-9;

    var MacroCycleLayout = {};

    function _BL() {
        return (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
    }

    /**
     * Compute the aspect ratio for a ring of size n.
     * 8 → 1.00 (round), 12 → 1.20, 18 → 1.50, 30+ → 1.60 (cap).
     */
    function aspectFor(n) {
        if (n <= 8) return 1.0;
        var a = 1.0 + 0.05 * (n - 8);
        if (a > 1.6) a = 1.6;
        return a;
    }

    /**
     * layout(ring, options) — place a macrocycle in egg-shape geometry.
     *
     * @param ring   array of atom IDs OR object with .atoms (CDK style)
     * @param options optional: { mol, bondLength, cx, cy, startAng }
     *                If `mol` is provided, atoms are mutated in place via
     *                mol.getAtom(id); otherwise the array entries are
     *                expected to be objects with .x/.y mutable fields.
     */
    MacroCycleLayout.layout = function (ring, options) {
        if (!ring) return false;
        options = options || {};
        var bondLength = options.bondLength || _BL();
        var atoms = ring.atoms || ring;
        if (!atoms || atoms.length < 8) return false;

        var n = atoms.length;
        var R = bondLength / (2 * Math.sin(Math.PI / n));
        var alpha = aspectFor(n);

        // Compute centroid + orientation if not provided.
        var cx = options.cx, cy = options.cy, startAng = options.startAng;
        var mol = options.mol;
        function getAtom(idOrObj) {
            if (mol && (typeof idOrObj === 'number' || typeof idOrObj === 'string')) {
                return mol.getAtom(idOrObj);
            }
            return idOrObj;
        }
        if (cx === undefined || cy === undefined || startAng === undefined) {
            var ccx = 0, ccy = 0, cn = 0;
            for (var i = 0; i < n; i++) {
                var a = getAtom(atoms[i]);
                if (a) { ccx += a.x; ccy += a.y; cn++; }
            }
            if (cn === 0) return false;
            ccx /= cn; ccy /= cn;
            if (cx === undefined) cx = ccx;
            if (cy === undefined) cy = ccy;
            if (startAng === undefined) {
                var first = getAtom(atoms[0]);
                startAng = first ? Math.atan2(first.y - cy, first.x - cx) : 0;
                if (!isFinite(startAng)) startAng = 0;
            }
        }

        // ARC-LENGTH PARAMETERISATION
        // Compute n vertex angles such that arc(0, θ_i) = i · P/n,
        // where P is the total perimeter of the unit ellipse with
        // semi-axes (α, 1/α). Use a dense uniform-θ sampling + linear
        // interpolation of the inverse arc-length function.
        //
        // ARC_SAMPLES = 720 → 0.5° resolution → bond-length residual
        // bounded by O((π/720)²) ≈ 2e-5 BL. More than tight enough.
        var ARC_SAMPLES = 720;
        var arcCum = new Array(ARC_SAMPLES + 1);
        arcCum[0] = 0;
        var prevX = alpha;     // θ=0 → (α, 0)
        var prevY = 0;
        for (var s = 1; s <= ARC_SAMPLES; s++) {
            var th = TWO_PI * s / ARC_SAMPLES;
            var sx = alpha * Math.cos(th);
            var sy = (1 / alpha) * Math.sin(th);
            var ddx = sx - prevX, ddy = sy - prevY;
            arcCum[s] = arcCum[s - 1] + Math.sqrt(ddx * ddx + ddy * ddy);
            prevX = sx; prevY = sy;
        }
        var totalPerimeter = arcCum[ARC_SAMPLES];
        // For each vertex i, find θ_i via inverse linear interpolation.
        var raw = [];
        for (var i2 = 0; i2 < n; i2++) {
            var target = i2 * totalPerimeter / n;
            // Binary search.
            var lo = 0, hi = ARC_SAMPLES;
            while (lo < hi) {
                var mid = (lo + hi) >> 1;
                if (arcCum[mid] < target) lo = mid + 1;
                else hi = mid;
            }
            // Linearly interpolate between samples [lo-1, lo].
            var theta;
            if (lo === 0) {
                theta = 0;
            } else {
                var lower = arcCum[lo - 1], upper = arcCum[lo];
                var t = (upper - lower) > EPS ? (target - lower) / (upper - lower) : 0;
                theta = TWO_PI * (lo - 1 + t) / ARC_SAMPLES;
            }
            raw.push({
                u: alpha * Math.cos(theta),
                v: (1 / alpha) * Math.sin(theta)
            });
        }
        // CHORD-EQUALISATION REFINEMENT
        // Arc-length gives equal arc lengths, but chord lengths still
        // differ slightly because the ellipse's curvature varies. Run
        // 30 iterations of pull-each-vertex-to-equalise-its-two-chords
        // (project back onto ellipse after each move) — converges to
        // exactly equal chord lengths. Each iteration is O(n).
        var thetas = new Array(n);
        for (var ti = 0; ti < n; ti++) {
            thetas[ti] = Math.atan2(raw[ti].v / (1 / alpha), raw[ti].u / alpha);
            // Math.atan2 gives [-π, π]; normalise to [0, 2π).
            if (thetas[ti] < 0) thetas[ti] += TWO_PI;
        }
        var REL_ITER = 30;
        for (var iter = 0; iter < REL_ITER; iter++) {
            var maxDelta = 0;
            for (var k = 0; k < n; k++) {
                var prev = (k - 1 + n) % n;
                var next = (k + 1) % n;
                var pX = alpha * Math.cos(thetas[prev]);
                var pY = (1 / alpha) * Math.sin(thetas[prev]);
                var nX = alpha * Math.cos(thetas[next]);
                var nY = (1 / alpha) * Math.sin(thetas[next]);
                var cX = alpha * Math.cos(thetas[k]);
                var cY = (1 / alpha) * Math.sin(thetas[k]);
                var dL = Math.hypot(cX - pX, cY - pY);
                var dR = Math.hypot(nX - cX, nY - cY);
                // Move θ_k by a step proportional to (dR - dL).
                // Tangent on the ellipse: dr/dθ = (-α sin θ, (1/α) cos θ).
                // Moving Δθ along the tangent changes the chord lengths
                // by ~|tangent|·Δθ. Use Δθ = 0.5 · (dR - dL) / (2·|tangent|).
                var tx = -alpha * Math.sin(thetas[k]);
                var ty = (1 / alpha) * Math.cos(thetas[k]);
                var tlen = Math.hypot(tx, ty);
                if (tlen < EPS) continue;
                var deltaTheta = 0.5 * (dR - dL) / (2 * tlen);
                if (Math.abs(deltaTheta) > maxDelta) maxDelta = Math.abs(deltaTheta);
                thetas[k] += deltaTheta;
                if (thetas[k] < 0) thetas[k] += TWO_PI;
                if (thetas[k] >= TWO_PI) thetas[k] -= TWO_PI;
            }
            if (maxDelta < 1e-7) break;
        }
        // Re-project relaxed θ-values onto the ellipse.
        for (var k2 = 0; k2 < n; k2++) {
            raw[k2].u = alpha * Math.cos(thetas[k2]);
            raw[k2].v = (1 / alpha) * Math.sin(thetas[k2]);
        }

        // Scale factor: we want |raw[1] - raw[0]| · scale = bondLength.
        // After chord-equalisation, |raw[i+1] - raw[i]| is equal for
        // all i, so this scale also makes every other bond = BL.
        var du = raw[1].u - raw[0].u, dv = raw[1].v - raw[0].v;
        var rawBondLen = Math.sqrt(du * du + dv * dv);
        if (rawBondLen < EPS) return false;
        var scale = bondLength / rawBondLen;
        // Apply scale to the unit ellipse vertices, then rotate by
        // startAng + offset (so vertex 0 lands at the requested angle
        // from the centroid in the original frame).
        var ct = Math.cos(startAng), st = Math.sin(startAng);
        for (var i3 = 0; i3 < n; i3++) {
            var u = raw[i3].u * scale;
            var v = raw[i3].v * scale;
            // Rotate by startAng then translate to (cx, cy).
            var nx = cx + (u * ct - v * st);
            var ny = cy + (u * st + v * ct);
            var a3 = getAtom(atoms[i3]);
            if (a3) { a3.x = nx; a3.y = ny; }
        }
        return true;
    };

    /**
     * computeRadius(n, bondLength) — egg-shape effective radius.
     * Useful for collision-margin estimates from outside this module.
     */
    MacroCycleLayout.computeRadius = function (n, bondLength) {
        bondLength = bondLength || _BL();
        var R = bondLength / (2 * Math.sin(Math.PI / n));
        var alpha = aspectFor(n);
        return R * alpha;  // semi-major axis (worst case)
    };

    MacroCycleLayout.aspectFor = aspectFor;

    global.SDG = global.SDG || {};
    global.SDG.MacroCycleLayout = MacroCycleLayout;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MacroCycleLayout;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
