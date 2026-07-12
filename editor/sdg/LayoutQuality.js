/**
 * editor/sdg/LayoutQuality.js — shared 2D depiction quality metrics.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Scores a drawn molecule using the same failure modes chemists notice on
 * screen: invalid coordinates, impossible bond lengths, atom crowding, acute
 * bond angles, and bond crossings. Lower penalty is better.
 */
(function (global) {
    'use strict';

    var LayoutQuality = {};
    var EPS = 1e-9;

    function bondLengthDefault() {
        return (global.Molecule && global.Molecule.BOND_LENGTH) || 30;
    }

    function dist(a, b) {
        var dx = a.x - b.x;
        var dy = a.y - b.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

    function key(a, b) {
        return a < b ? a + ':' + b : b + ':' + a;
    }

    function finiteAtom(a) {
        return !!a && isFinite(a.x) && isFinite(a.y);
    }

    function orient(a, b, c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    function segmentsCross(a, b, c, d) {
        var o1 = orient(a, b, c);
        var o2 = orient(a, b, d);
        var o3 = orient(c, d, a);
        var o4 = orient(c, d, b);
        if (Math.abs(o1) < EPS || Math.abs(o2) < EPS ||
            Math.abs(o3) < EPS || Math.abs(o4) < EPS) {
            return false;
        }
        return (o1 > 0) !== (o2 > 0) && (o3 > 0) !== (o4 > 0);
    }

    function atomIdsFor(mol, atomIds) {
        if (atomIds && atomIds.length) return atomIds.slice();
        if (!mol || !mol.atoms) return [];
        var ids = [];
        for (var i = 0; i < mol.atoms.length; i++) ids.push(mol.atoms[i].id);
        return ids;
    }

    function countComponents(ids, selectedBonds) {
        var adj = {};
        var seen = {};
        var count = 0;
        for (var i = 0; i < ids.length; i++) adj[ids[i]] = [];
        for (var bi = 0; bi < selectedBonds.length; bi++) {
            var bond = selectedBonds[bi].bond;
            if (!adj[bond.atom1] || !adj[bond.atom2]) continue;
            adj[bond.atom1].push(bond.atom2);
            adj[bond.atom2].push(bond.atom1);
        }
        for (var ii = 0; ii < ids.length; ii++) {
            var start = ids[ii];
            if (seen[start]) continue;
            count++;
            var stack = [start];
            seen[start] = true;
            while (stack.length) {
                var current = stack.pop();
                var nbrs = adj[current] || [];
                for (var ni = 0; ni < nbrs.length; ni++) {
                    if (seen[nbrs[ni]]) continue;
                    seen[nbrs[ni]] = true;
                    stack.push(nbrs[ni]);
                }
            }
        }
        return count;
    }

    function isSaturatedCageCrossing(mol, b1, b2, degree, ringRank) {
        if (ringRank < 3 || b1.type !== 1 || b2.type !== 1) return false;
        var atomIds = [b1.atom1, b1.atom2, b2.atom1, b2.atom2];
        var branchAtoms = 0;
        for (var i = 0; i < atomIds.length; i++) {
            var atom = mol.getAtom ? mol.getAtom(atomIds[i]) : null;
            if (!atom || atom.symbol !== 'C' || atom.aromatic) return false;
            var deg = degree[atomIds[i]] || 0;
            if (deg < 2 || deg > 4) return false;
            if (deg >= 3) branchAtoms++;
        }
        return branchAtoms >= 2;
    }

    /**
     * evaluate(mol, atomIds, options)
     *
     * Returns:
     *   {
     *     penalty,
     *     hardFailures,
     *     needsAdaptiveRescue,
     *     duplicateAtoms, severeClashes, closePairs,
     *     invalidBonds, nonFiniteAtoms,
     *     shortBonds, longBonds, crossings, cageCrossings, acuteAngles,
     *     bounds
     *   }
     */
    LayoutQuality.evaluate = function (mol, atomIds, options) {
        options = options || {};
        var ids = atomIdsFor(mol, atomIds);
        var idSet = {};
        var metrics = {
            penalty: 0,
            hardFailures: 0,
            needsAdaptiveRescue: false,
            duplicateAtoms: 0,
            severeClashes: 0,
            closePairs: 0,
            invalidBonds: 0,
            nonFiniteAtoms: 0,
            shortBonds: 0,
            longBonds: 0,
            bondWarnings: 0,
            crossings: 0,
            cageCrossings: 0,
            acuteAngles: 0,
            atoms: ids.length,
            bonds: 0,
            bounds: null
        };

        if (!mol || !mol.atoms) {
            metrics.hardFailures++;
            metrics.penalty += 1000;
            return metrics;
        }

        var bondLength = options.bondLength || bondLengthDefault();
        var duplicateDist = options.duplicateDist || 0.25;
        var severeClashDist = options.severeClashDist || bondLength * 0.35;
        var closePairDist = options.closePairDist || bondLength * 0.45;
        var rescueShort = options.rescueShort || bondLength * 0.55;
        var rescueLong = options.rescueLong || bondLength * 1.45;
        var warnShort = options.warnShort || bondLength * 0.65;
        var warnLong = options.warnLong || bondLength * 1.25;
        var acuteRadians = options.acuteRadians || (Math.PI / 18); // 10 deg

        for (var i = 0; i < ids.length; i++) idSet[ids[i]] = true;

        var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        for (var ai = 0; ai < ids.length; ai++) {
            var atom = mol.getAtom ? mol.getAtom(ids[ai]) : null;
            if (!finiteAtom(atom)) {
                metrics.nonFiniteAtoms++;
                metrics.hardFailures++;
                metrics.penalty += 1000;
                continue;
            }
            if (atom.x < minX) minX = atom.x;
            if (atom.x > maxX) maxX = atom.x;
            if (atom.y < minY) minY = atom.y;
            if (atom.y > maxY) maxY = atom.y;
        }
        if (minX < Infinity) {
            metrics.bounds = {
                minX: minX,
                minY: minY,
                maxX: maxX,
                maxY: maxY,
                width: maxX - minX,
                height: maxY - minY
            };
        }

        var bonded = {};
        var degree = {};
        var selectedBonds = [];
        for (var bi = 0; bi < mol.bonds.length; bi++) {
            var bond = mol.bonds[bi];
            if (!bond || bond.atom1 === bond.atom2 ||
                !mol.getAtom(bond.atom1) || !mol.getAtom(bond.atom2)) {
                metrics.invalidBonds++;
                metrics.hardFailures++;
                metrics.penalty += 1000;
                continue;
            }
            bonded[key(bond.atom1, bond.atom2)] = true;
            degree[bond.atom1] = (degree[bond.atom1] || 0) + 1;
            degree[bond.atom2] = (degree[bond.atom2] || 0) + 1;
            if (!idSet[bond.atom1] || !idSet[bond.atom2]) continue;

            var ba = mol.getAtom(bond.atom1);
            var bb = mol.getAtom(bond.atom2);
            if (!finiteAtom(ba) || !finiteAtom(bb)) continue;
            var len = dist(ba, bb);
            metrics.bonds++;
            selectedBonds.push({ bond: bond, a: ba, b: bb });

            if (len < rescueShort) {
                metrics.shortBonds++;
                metrics.needsAdaptiveRescue = true;
            }
            if (len > rescueLong) {
                metrics.longBonds++;
                metrics.needsAdaptiveRescue = true;
            }
            if (len < warnShort) {
                metrics.bondWarnings++;
                metrics.penalty += (warnShort - len) / bondLength * 10;
            } else if (len > warnLong) {
                metrics.bondWarnings++;
                metrics.penalty += (len - warnLong) / bondLength * 10;
            }
        }

        for (var pi = 0; pi < ids.length; pi++) {
            var pa = mol.getAtom(ids[pi]);
            if (!finiteAtom(pa)) continue;
            for (var pj = pi + 1; pj < ids.length; pj++) {
                if (bonded[key(ids[pi], ids[pj])]) continue;
                var pb = mol.getAtom(ids[pj]);
                if (!finiteAtom(pb)) continue;
                var d = dist(pa, pb);
                if (d < duplicateDist) {
                    metrics.duplicateAtoms++;
                    metrics.hardFailures++;
                    metrics.needsAdaptiveRescue = true;
                    metrics.penalty += 500;
                } else if (d < severeClashDist) {
                    metrics.severeClashes++;
                    metrics.hardFailures++;
                    metrics.needsAdaptiveRescue = true;
                    metrics.penalty += 200 + (severeClashDist - d) / bondLength * 50;
                } else if (d < closePairDist) {
                    metrics.closePairs++;
                    metrics.penalty += 2 + (closePairDist - d) / bondLength * 20;
                }
            }
        }

        var componentCount = countComponents(ids, selectedBonds);
        var ringRank = Math.max(0, selectedBonds.length - ids.length + componentCount);

        for (var ci = 0; ci < selectedBonds.length; ci++) {
            var b1 = selectedBonds[ci].bond;
            for (var cj = ci + 1; cj < selectedBonds.length; cj++) {
                var b2 = selectedBonds[cj].bond;
                if (b1.atom1 === b2.atom1 || b1.atom1 === b2.atom2 ||
                    b1.atom2 === b2.atom1 || b1.atom2 === b2.atom2) {
                    continue;
                }
                if (segmentsCross(selectedBonds[ci].a, selectedBonds[ci].b,
                                  selectedBonds[cj].a, selectedBonds[cj].b)) {
                    if (isSaturatedCageCrossing(mol, b1, b2, degree, ringRank)) {
                        metrics.cageCrossings++;
                        continue;
                    }
                    metrics.crossings++;
                    metrics.penalty += 6;
                }
            }
        }

        for (var angi = 0; angi < ids.length; angi++) {
            var center = mol.getAtom(ids[angi]);
            if (!finiteAtom(center) || !mol.getNeighbors) continue;
            var nbrs = mol.getNeighbors(ids[angi]) || [];
            var local = [];
            for (var ni = 0; ni < nbrs.length; ni++) {
                if (idSet[nbrs[ni]]) local.push(nbrs[ni]);
            }
            for (var n1 = 0; n1 < local.length; n1++) {
                var a1 = mol.getAtom(local[n1]);
                if (!finiteAtom(a1)) continue;
                var v1x = a1.x - center.x, v1y = a1.y - center.y;
                var l1 = Math.sqrt(v1x * v1x + v1y * v1y);
                if (l1 < EPS) continue;
                for (var n2 = n1 + 1; n2 < local.length; n2++) {
                    var a2 = mol.getAtom(local[n2]);
                    if (!finiteAtom(a2)) continue;
                    var v2x = a2.x - center.x, v2y = a2.y - center.y;
                    var l2 = Math.sqrt(v2x * v2x + v2y * v2y);
                    if (l2 < EPS) continue;
                    var cosang = (v1x * v2x + v1y * v2y) / (l1 * l2);
                    if (cosang > 1) cosang = 1;
                    if (cosang < -1) cosang = -1;
                    var angle = Math.acos(cosang);
                    if (angle < acuteRadians) {
                        metrics.acuteAngles++;
                        metrics.penalty += 4 + (acuteRadians - angle) * 10;
                    }
                }
            }
        }

        return metrics;
    };

    LayoutQuality.penalty = function (mol, atomIds, options) {
        return LayoutQuality.evaluate(mol, atomIds, options).penalty;
    };

    LayoutQuality.needsAdaptiveRescue = function (mol, atomIds, options) {
        return LayoutQuality.evaluate(mol, atomIds, options).needsAdaptiveRescue;
    };

    LayoutQuality.isBetter = function (candidate, incumbent, tolerance) {
        tolerance = tolerance === undefined ? 1e-6 : tolerance;
        if (!incumbent) return true;
        if (!candidate) return false;
        if (candidate.hardFailures !== incumbent.hardFailures) {
            return candidate.hardFailures < incumbent.hardFailures;
        }
        return candidate.penalty < incumbent.penalty - tolerance;
    };

    global.SDG = global.SDG || {};
    global.SDG.LayoutQuality = LayoutQuality;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = LayoutQuality;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
