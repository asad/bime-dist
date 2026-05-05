/**
 * editor/sdg/CorrectGeometricConfiguration.js — CDK port (v1.8.17 full).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * IP NOTICE — clean-room re-implementation, no source copy:
 *   Independent JavaScript implementation of the algorithm described
 *   in CDK's CorrectGeometricConfiguration (LGPL-2.1) and Helson '99
 *   SDG. NO CDK source code copied. The CDK URL is an algorithm
 *   reference only.
 *
 * Reference (Java original, ~13 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/CorrectGeometricConfiguration.java
 *
 * For each E/Z (cis/trans) double bond, checks whether the depicted 2D
 * coordinates match the chemical configuration encoded in
 * bond.cipLabel (set by editor/CipStereo.js). If the depicted
 * configuration is inverted, reflects the smaller-half substituents
 * across the double-bond axis.
 *
 * Status (v1.8.17): FULL implementation. Replaces the v1.8.13 stub.
 */
(function (global) {
    'use strict';

    var EPS = 1e-9;

    var CorrectGeometricConfiguration = {};

    /**
     * correct(mol) — top-level entry. Returns the count of corrections
     * applied (number of double bonds reflected).
     */
    CorrectGeometricConfiguration.correct = function (mol) {
        if (!mol || !mol.atoms || !mol.bonds) return 0;
        var corrected = 0;
        for (var b = 0; b < mol.bonds.length; b++) {
            var bnd = mol.bonds[b];
            if (bnd.type !== 2) continue;
            if (bnd.cipLabel !== 'E' && bnd.cipLabel !== 'Z') continue;
            if (CorrectGeometricConfiguration._correctOne(mol, bnd)) corrected++;
        }
        return corrected;
    };

    /**
     * _correctOne(mol, bond) — correct one E/Z double bond.
     */
    CorrectGeometricConfiguration._correctOne = function (mol, bond) {
        var a1 = mol.getAtom(bond.atom1);
        var a2 = mol.getAtom(bond.atom2);
        if (!a1 || !a2) return false;

        // Find the priority-1 substituent on each side.
        var sub1 = CorrectGeometricConfiguration._priorityNeighbour(mol, a1.id, a2.id);
        var sub2 = CorrectGeometricConfiguration._priorityNeighbour(mol, a2.id, a1.id);
        if (!sub1 || !sub2) return false;
        var s1 = mol.getAtom(sub1);
        var s2 = mol.getAtom(sub2);
        if (!s1 || !s2) return false;

        // Bond axis direction.
        var bdx = a2.x - a1.x, bdy = a2.y - a1.y;
        var blen = Math.hypot(bdx, bdy);
        if (blen < EPS) return false;
        var nx = -bdy / blen, ny = bdx / blen;  // perpendicular

        // Side of axis (sign of perpendicular projection).
        var s1Side = (s1.x - a1.x) * nx + (s1.y - a1.y) * ny;
        var s2Side = (s2.x - a2.x) * nx + (s2.y - a2.y) * ny;

        // Same sign → cis (Z); opposite → trans (E).
        var measuredZ = (s1Side * s2Side > 0);
        var wantsZ = bond.cipLabel === 'Z';
        if (measuredZ === wantsZ) return false;  // already correct

        // Need to flip ONE side. Pick the smaller-subtree side.
        var subtree1 = CorrectGeometricConfiguration._collectSubtree(mol, sub1, a1.id);
        var subtree2 = CorrectGeometricConfiguration._collectSubtree(mol, sub2, a2.id);
        var flipSide, pivotAtom;
        if (subtree1.length <= subtree2.length) {
            flipSide = subtree1; pivotAtom = a1;
        } else {
            flipSide = subtree2; pivotAtom = a2;
        }

        // Reflect each atom across the line through the bond axis.
        var dx = bdx / blen, dy = bdy / blen;
        for (var i = 0; i < flipSide.length; i++) {
            var pa = mol.getAtom(flipSide[i]);
            if (!pa) continue;
            var px = pa.x - pivotAtom.x, py = pa.y - pivotAtom.y;
            var t = px * dx + py * dy;
            var prx = pivotAtom.x + t * dx;
            var pry = pivotAtom.y + t * dy;
            pa.x = 2 * prx - pa.x;
            pa.y = 2 * pry - pa.y;
        }
        return true;
    };

    /**
     * _priorityNeighbour(mol, atomId, excludeId) — highest-priority
     * neighbour by atomic number (CIP proxy).
     */
    CorrectGeometricConfiguration._priorityNeighbour = function (mol, atomId, excludeId) {
        var nbrs = mol.getNeighbors(atomId) || [];
        var best = null, bestPri = -1;
        var ATOMIC_NUMBERS = {
            H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, P: 15, S: 16,
            Cl: 17, Br: 35, I: 53
        };
        for (var i = 0; i < nbrs.length; i++) {
            if (nbrs[i] === excludeId) continue;
            var na = mol.getAtom(nbrs[i]);
            if (!na) continue;
            var pri = ATOMIC_NUMBERS[na.symbol] || 0;
            if (pri > bestPri || best === null) {
                best = nbrs[i];
                bestPri = pri;
            }
        }
        return best;
    };

    /**
     * _collectSubtree(mol, startId, blockedId) — BFS reachability set.
     */
    CorrectGeometricConfiguration._collectSubtree = function (mol, startId, blockedId) {
        var visited = {};
        var queue = [startId];
        var result = [];
        visited[blockedId] = true;
        visited[startId] = true;
        while (queue.length > 0) {
            var cur = queue.shift();
            result.push(cur);
            var nbrs = mol.getNeighbors(cur) || [];
            for (var i = 0; i < nbrs.length; i++) {
                if (visited[nbrs[i]]) continue;
                visited[nbrs[i]] = true;
                queue.push(nbrs[i]);
            }
        }
        return result;
    };

    global.SDG = global.SDG || {};
    global.SDG.CorrectGeometricConfiguration = CorrectGeometricConfiguration;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = CorrectGeometricConfiguration;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
