/**
 * MLDepict.js — v1.6.0 tiny ML 2D-coordinate residual predictor.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Pure-JS forward pass for the v1.6.0 residual coord network. Reads the
 * shipped weights from editor/ml-depict-weights.json (or a globally injected
 * MLDepictWeights object), and exposes:
 *
 *   MLDepict.ready()                        -> bool, are weights loaded
 *   MLDepict.featuresFor(mol, atomId, ringSet) -> 20-ch feature vector
 *   MLDepict.predictResidual(focal, nbrPool)   -> [dx, dy] in BOND_LENGTH units
 *   MLDepict.refineLayout(mol, alpha)          -> blends ML coords into atoms
 *
 * The model:
 *   Input  : 40 ch  (20 focal-atom features + 20 neighbour-pooled features)
 *   FC1    : 40 -> 64   (Linear + ReLU)
 *   FC2    : 64 -> 32   (Linear + ReLU)
 *   Out    : 32 -> 2    (Linear)        -> (dx, dy) residual in BL units
 *
 * Final coord = neighbour_centroid + residual * BOND_LENGTH.
 *
 * Default Layout.options.useMLDepict = false; refineLayout is OFF unless
 * the caller flips that flag (and supplies an alpha > 0 blend weight).
 *
 * Determinism: forward pass is pure float32 arithmetic, no Math.random,
 * same input -> byte-identical output.
 */
(function (global) {
    'use strict';

    var MLDepict = {};
    var W = null;
    var FEATURE_DIM = 20;
    var ELEM_ONE_HOT = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'H'];

    // ---------------------------------------------------------------------
    // Weight loader. Tries (in order):
    //   1. global.MLDepictWeights (injected by host page)
    //   2. fetch('editor/ml-depict-weights.json') in browser
    //   3. require('./ml-depict-weights.json') in Node
    // ---------------------------------------------------------------------
    function loadWeightsSync() {
        if (W) return W;
        if (global.MLDepictWeights) { W = global.MLDepictWeights; return W; }
        // Node fallback for tests.
        if (typeof require === 'function' && typeof module !== 'undefined') {
            try {
                W = require('./ml-depict-weights.json');
                return W;
            } catch (e) { /* not found */ }
        }
        return null;
    }

    MLDepict.loadWeights = function (weights) {
        W = weights;
    };

    MLDepict.ready = function () {
        return loadWeightsSync() !== null;
    };

    MLDepict.weights = function () { return loadWeightsSync(); };

    // ---------------------------------------------------------------------
    // Feature encoder (must match tools/extract-ml-depict-dataset.js).
    // ---------------------------------------------------------------------
    MLDepict.buildInRingSet = function (mol) {
        var set = {};
        if (!mol || !mol.atoms || mol.atoms.length >= 200) return set;
        var rings = mol.findRings ? mol.findRings(8) : [];
        for (var i = 0; i < rings.length; i++) {
            for (var j = 0; j < rings[i].atoms.length; j++) {
                set[rings[i].atoms[j]] = true;
            }
        }
        return set;
    };

    MLDepict.featuresFor = function (mol, atomId, inRingSet) {
        var atom = mol.getAtom(atomId);
        var nbrs = mol.getNeighbors(atomId);
        var f = new Array(FEATURE_DIM);
        for (var z = 0; z < FEATURE_DIM; z++) f[z] = 0;
        var elemIdx = ELEM_ONE_HOT.indexOf(atom.symbol);
        if (elemIdx >= 0) f[elemIdx] = 1;
        f[10] = Math.min(nbrs.length, 4) / 4;
        f[11] = atom.aromatic ? 1 : 0;
        f[12] = (inRingSet && inRingSet[atomId]) ? 1 : 0;
        var deg = nbrs.length;
        if (atom.aromatic || deg === 3) f[14] = 1;
        else if (deg <= 2) f[13] = 1;
        else f[15] = 1;
        f[16] = (atom.charge || 0) / 2;
        f[17] = (atom.symbol !== 'C' && atom.symbol !== 'H') ? 1 : 0;
        f[18] = (atom.symbol === 'F' || atom.symbol === 'Cl' ||
                 atom.symbol === 'Br' || atom.symbol === 'I') ? 1 : 0;
        f[19] = Math.min(atom.hydrogens >= 0 ? atom.hydrogens : 0, 4) / 4;
        return f;
    };

    // ---------------------------------------------------------------------
    // Forward pass. Architecture is described by the bundled weights blob.
    // The forward pass auto-detects the arch via w.arch.in (40 vs 84).
    // Older weights still work; new weights need new inputs.
    // ---------------------------------------------------------------------
    function relu(x) { return x > 0 ? x : 0; }

    // predictResidual is variadic on extra inputs. Callers can pass:
    //   (focalFeat, nbrPool)                                   v1.6.x / v1.7.x
    //   (focalFeat, nbrPool1, nbrPool2, nbrPool3, sourceTag)   v1.8.0
    //
    // sourceTag is one of arch.sources (default 'unknown'). Pass null /
    // undefined to use the unknown channel.
    MLDepict.predictResidual = function (focalFeat, nbrPool1, nbrPool2,
                                         nbrPool3, sourceTag) {
        var w = loadWeightsSync();
        if (!w) return [0, 0];
        var IN = w.arch.in, H1 = w.arch.h1, H2 = w.arch.h2, OUT = w.arch.out;
        var W1 = w.W1, b1 = w.b1, W2 = w.W2, b2 = w.b2, W3 = w.W3, b3 = w.b3;
        var FD = w.arch.feature_dim || FEATURE_DIM;

        // Build input vector. The shape is recorded by w.arch.in:
        //   IN == 40  -> [focal(20) | 1-hop(20)]                 v1.6.x/v1.7.x
        //   IN == 84  -> [focal(20) | 1-hop(20) | 2-hop(20) |
        //                 3-hop(20) | source-one-hot(4)]         v1.8.0
        var x = new Array(IN);
        for (var z0 = 0; z0 < IN; z0++) x[z0] = 0;

        // Channel 0..FD-1: focal features (always present).
        for (var i = 0; i < FD; i++) x[i] = focalFeat[i] || 0;

        // Channel FD..2FD-1: 1-hop pool (always present).
        if (nbrPool1) {
            for (var i2 = 0; i2 < FD; i2++) x[FD + i2] = nbrPool1[i2] || 0;
        }

        if (IN >= 4 * FD) {
            // 2-hop, 3-hop pools.
            if (nbrPool2) {
                for (var i3 = 0; i3 < FD; i3++) x[2*FD + i3] = nbrPool2[i3] || 0;
            }
            if (nbrPool3) {
                for (var i4 = 0; i4 < FD; i4++) x[3*FD + i4] = nbrPool3[i4] || 0;
            }
            // Source one-hot at the tail.
            var sources = (w.arch && w.arch.sources) || ['bime','external','curated','unknown'];
            var tag = sourceTag;
            if (tag == null || sources.indexOf(tag) < 0) tag = 'unknown';
            var srcIdx = sources.indexOf(tag);
            if (srcIdx < 0) srcIdx = sources.length - 1;
            x[4*FD + srcIdx] = 1;
        }

        // FC1
        var h1 = new Array(H1);
        for (var j = 0; j < H1; j++) {
            var zz = b1[j];
            for (var k = 0; k < IN; k++) zz += x[k] * W1[k * H1 + j];
            h1[j] = relu(zz);
        }
        // FC2
        var h2 = new Array(H2);
        for (var j2 = 0; j2 < H2; j2++) {
            var z2 = b2[j2];
            for (var k2 = 0; k2 < H1; k2++) z2 += h1[k2] * W2[k2 * H2 + j2];
            h2[j2] = relu(z2);
        }
        // OUT
        var out = new Array(OUT);
        for (var jo = 0; jo < OUT; jo++) {
            var zo = b3[jo];
            for (var ko = 0; ko < H2; ko++) zo += h2[ko] * W3[ko * OUT + jo];
            out[jo] = zo;
        }
        return out;
    };

    // ---------------------------------------------------------------------
    // refineLayout(mol, alpha) — blend ML residual into existing 2D coords.
    //
    //   alpha = 0    -> coords unchanged (default for safety)
    //   alpha = 1    -> atoms snap exactly to (centroid + ML residual)
    //   alpha 0..1   -> linear blend
    //
    // Caller is responsible for passing alpha. Layout.js looks up
    // Layout.options.mlDepictWeight (default 0).
    //
    // Returns count of atoms actually refined.
    // ---------------------------------------------------------------------
    MLDepict.refineLayout = function (mol, alpha, sourceTag) {
        if (!mol || !mol.atoms || alpha <= 0) return 0;
        var w = loadWeightsSync();
        if (!w) return 0;
        var BL = (typeof Molecule !== 'undefined' && Molecule.BOND_LENGTH) || 30;
        var threeHop = (w.arch.in >= 4 * (w.arch.feature_dim || FEATURE_DIM));

        // Pre-compute focal features for every atom.
        var inRingSet = MLDepict.buildInRingSet(mol);
        var feats = new Array(mol.atoms.length);
        var idToIdx = {};
        for (var i = 0; i < mol.atoms.length; i++) {
            feats[i] = MLDepict.featuresFor(mol, mol.atoms[i].id, inRingSet);
            idToIdx[mol.atoms[i].id] = i;
        }

        // BFS-pool helper for 2-hop / 3-hop. Identical maths to the
        // training-time extractor in tools/extract-ml-depict-*.js.
        function khopPool(rootId, hop) {
            var pool = new Array(FEATURE_DIM);
            for (var z = 0; z < FEATURE_DIM; z++) pool[z] = 0;
            var seen = {}; seen[rootId] = 0;
            var frontier = [rootId];
            for (var step = 1; step <= hop; step++) {
                var next = [];
                for (var fi = 0; fi < frontier.length; fi++) {
                    var ns = mol.getNeighbors(frontier[fi]);
                    for (var ni = 0; ni < ns.length; ni++) {
                        if (seen[ns[ni]] !== undefined) continue;
                        seen[ns[ni]] = step;
                        next.push(ns[ni]);
                    }
                }
                frontier = next;
                if (!frontier.length) break;
            }
            if (!frontier.length) return pool;
            for (var pi = 0; pi < frontier.length; pi++) {
                var pf = feats[idToIdx[frontier[pi]]];
                for (var c = 0; c < FEATURE_DIM; c++) pool[c] += pf[c];
            }
            for (var c2 = 0; c2 < FEATURE_DIM; c2++) pool[c2] /= frontier.length;
            return pool;
        }

        var refined = 0;
        for (var ai = 0; ai < mol.atoms.length; ai++) {
            var atom = mol.atoms[ai];
            var nbrs = mol.getNeighbors(atom.id);
            if (nbrs.length === 0) continue;

            // Neighbour centroid + 1-hop pool.
            var cx = 0, cy = 0;
            var nf1 = new Array(FEATURE_DIM);
            for (var nz = 0; nz < FEATURE_DIM; nz++) nf1[nz] = 0;
            for (var k = 0; k < nbrs.length; k++) {
                var nbAtom = mol.getAtom(nbrs[k]);
                cx += nbAtom.x; cy += nbAtom.y;
                var nbf = feats[idToIdx[nbrs[k]]];
                for (var c = 0; c < FEATURE_DIM; c++) nf1[c] += nbf[c];
            }
            cx /= nbrs.length; cy /= nbrs.length;
            for (var c2 = 0; c2 < FEATURE_DIM; c2++) nf1[c2] /= nbrs.length;

            var nf2 = null, nf3 = null;
            if (threeHop) {
                nf2 = khopPool(atom.id, 2);
                nf3 = khopPool(atom.id, 3);
            }

            var r = MLDepict.predictResidual(feats[ai], nf1, nf2, nf3,
                                              sourceTag || 'unknown');
            // Sanity: clamp to [-2, 2] BL — anything bigger is a model
            // glitch and we'd rather keep the rule-based coord.
            var rx = Math.max(-2, Math.min(2, r[0]));
            var ry = Math.max(-2, Math.min(2, r[1]));
            var mlx = cx + rx * BL;
            var mly = cy + ry * BL;

            atom.x = (1 - alpha) * atom.x + alpha * mlx;
            atom.y = (1 - alpha) * atom.y + alpha * mly;
            refined++;
        }
        return refined;
    };

    // ---------------------------------------------------------------------
    // Diagnostics: report best-effort summary about loaded weights.
    // ---------------------------------------------------------------------
    MLDepict.summary = function () {
        var w = loadWeightsSync();
        if (!w) return { ready: false };
        return {
            ready: true,
            version: w.version,
            arch: w.arch,
            training: w.training,
            param_count:
                w.W1.length + w.b1.length +
                w.W2.length + w.b2.length +
                w.W3.length + w.b3.length
        };
    };

    // ---------------------------------------------------------------------
    // Export.
    // ---------------------------------------------------------------------
    global.MLDepict = MLDepict;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = MLDepict;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
