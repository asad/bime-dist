/**
 * editor/sdg/OverlapResolver.js — clean-room JS port of CDK's OverlapResolver.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Reference (Java original):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/OverlapResolver.java
 *
 * Iterative summed-force overlap resolver. Unlike pairwise sequential
 * pushes (which create secondary overlaps), this computes ALL repulsion
 * + spring forces SIMULTANEOUSLY then applies them in a single
 * integration step per iteration.
 *
 * This module supersedes the inline copy in editor/SDGLayout.js. SDGLayout
 * delegates to SDG.OverlapResolver when both are loaded (full v1.8.x port).
 *
 * Public API (matches CDK):
 *   OverlapResolver.resolveOverlap(mol, atomIds, options) → number iters used
 *
 * Options:
 *   bondLength    : target bond length (default Molecule.BOND_LENGTH)
 *   minDist       : repulsion threshold (default bondLength * 0.6)
 *   maxIters      : iteration cap (default 60)
 *   stepSize      : displacement scale per iteration (default 0.5)
 *   convergeTol   : max-displacement convergence threshold (default 0.05)
 *   springCoeff   : spring-force scaling for bonded pairs (default 1.0)
 *
 * v1.8.12 noted that springCoeff=1.0 was too strong for long polycyclic
 * chains (created secondary overlaps). v1.8.13+ tuning targets:
 *   springCoeff: 0.3–0.5
 *   stepSize: cosine-decay 0.5 → 0.1
 */
(function (global) {
    'use strict';

    var OverlapResolver = {};

    OverlapResolver.resolveOverlap = function (mol, atomIds, opts) {
        opts = opts || {};
        var BL = (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
        var bondLength = opts.bondLength || BL;
        var minDist    = opts.minDist || (bondLength * 0.6);
        var maxIters   = opts.maxIters || 60;
        var stepSize   = opts.stepSize || 0.5;
        var springC    = (opts.springCoeff !== undefined) ? opts.springCoeff : 1.0;
        var convergeTol = opts.convergeTol || 0.05;

        if (!mol || !atomIds || atomIds.length < 2) return 0;

        var bondedSet = {};
        for (var bi = 0; bi < mol.bonds.length; bi++) {
            var b = mol.bonds[bi];
            if (atomIds.indexOf(b.atom1) < 0 || atomIds.indexOf(b.atom2) < 0) continue;
            bondedSet[b.atom1 + ',' + b.atom2] = true;
            bondedSet[b.atom2 + ',' + b.atom1] = true;
        }

        var n = atomIds.length;
        var fx = new Float64Array(n);
        var fy = new Float64Array(n);
        var EPS = 1e-6;

        for (var iter = 0; iter < maxIters; iter++) {
            for (var z = 0; z < n; z++) { fx[z] = 0; fy[z] = 0; }

            for (var i = 0; i < n; i++) {
                var ai = mol.getAtom(atomIds[i]);
                if (!ai) continue;
                for (var j = i + 1; j < n; j++) {
                    var aj = mol.getAtom(atomIds[j]);
                    if (!aj) continue;
                    var dx = aj.x - ai.x;
                    var dy = aj.y - ai.y;
                    var d2 = dx * dx + dy * dy;
                    var d = Math.sqrt(d2);
                    if (d < EPS) {
                        var phi = (atomIds[i] + atomIds[j]) * 0.5;
                        dx = Math.cos(phi);
                        dy = Math.sin(phi);
                        d = 1;
                    }
                    var ux = dx / d, uy = dy / d;
                    if (bondedSet[atomIds[i] + ',' + atomIds[j]]) {
                        var spring = springC * (d - bondLength);
                        fx[i] += spring * ux;
                        fy[i] += spring * uy;
                        fx[j] -= spring * ux;
                        fy[j] -= spring * uy;
                    } else if (d < minDist) {
                        var rep = (minDist - d);
                        fx[i] -= rep * ux;
                        fy[i] -= rep * uy;
                        fx[j] += rep * ux;
                        fy[j] += rep * uy;
                    }
                }
            }

            var maxDisp = 0;
            for (var k = 0; k < n; k++) {
                var ak = mol.getAtom(atomIds[k]);
                if (!ak) continue;
                var ddx = fx[k] * stepSize;
                var ddy = fy[k] * stepSize;
                ak.x += ddx;
                ak.y += ddy;
                var disp = Math.sqrt(ddx * ddx + ddy * ddy);
                if (disp > maxDisp) maxDisp = disp;
            }
            if (maxDisp < convergeTol) return iter + 1;
        }
        return maxIters;
    };

    global.SDG = global.SDG || {};
    global.SDG.OverlapResolver = OverlapResolver;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = OverlapResolver;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
