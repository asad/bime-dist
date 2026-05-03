/**
 * Templates.js — Pre-computed coordinate templates for common molecular scaffolds
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Provides ready-to-use 2D coordinate sets for frequently drawn structures.
 * Each template is an object with atoms [{symbol, x, y}] and bonds [{a1, a2, type}]
 * where a1/a2 are zero-based indices into the atoms array.
 * All coordinates use BOND_LENGTH = 30px as the standard bond length.
 */
(function(global) {
    'use strict';

    var BL = Molecule.BOND_LENGTH; // 30
    var SINGLE = Molecule.BOND_SINGLE;
    var DOUBLE = Molecule.BOND_DOUBLE;

    // Helper: generate regular polygon coordinates
    function polygon(n, cx, cy, startAngle) {
        var radius = BL / (2 * Math.sin(Math.PI / n));
        var angle = startAngle || -Math.PI / 2;
        if (n % 2 === 0) angle += Math.PI / n;
        var coords = [];
        for (var i = 0; i < n; i++) {
            coords.push({
                x: cx + radius * Math.cos(angle + i * 2 * Math.PI / n),
                y: cy + radius * Math.sin(angle + i * 2 * Math.PI / n)
            });
        }
        return coords;
    }

    // Helper: ring bonds (alternating double for aromatic)
    function ringBonds(n, startIdx, alternateDouble) {
        var bonds = [];
        for (var i = 0; i < n; i++) {
            var type = SINGLE;
            if (alternateDouble && i % 2 === 0) type = DOUBLE;
            bonds.push({ a1: startIdx + i, a2: startIdx + (i + 1) % n, type: type });
        }
        return bonds;
    }

    var Templates = {};

    // =========================================================================
    // Benzene (6-membered aromatic ring)
    // =========================================================================
    Templates.benzene = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true); // alternating double bonds
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Cyclohexane (6-membered all single bonds)
    // =========================================================================
    Templates.cyclohexane = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Cyclopentane (5-membered ring)
    // =========================================================================
    Templates.cyclopentane = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(5, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Naphthalene (two fused benzene rings)
    // =========================================================================
    Templates.naphthalene = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        // First ring bonds
        var bonds = ringBonds(6, 0, true);

        // Second ring shares edge between atoms 1-2 (right side)
        // Compute positions for the 4 new atoms
        var shared1 = atoms[1];
        var shared2 = atoms[2];
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        // Normal pointing outward (right)
        var nx = -dy;
        var ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        nx /= len; ny /= len;

        var radius2 = BL / (2 * Math.sin(Math.PI / 6));
        var apothem = radius2 * Math.cos(Math.PI / 6);
        var cx2 = mx + nx * apothem;
        var cy2 = my + ny * apothem;

        var r6b = polygon(6, cx2, cy2);
        // Find the 4 atoms of ring 2 that don't overlap with ring 1
        var newAtomStart = atoms.length;
        // The shared atoms from ring2 closest to atoms[1] and atoms[2]
        // We need to find the correct ring2 vertices
        var ring2Order = reorderToMatch(r6b, shared1, shared2);
        // ring2Order[0] matches shared1, ring2Order[1] matches shared2
        // Add atoms 2,3,4,5 from ring2
        for (var i = 2; i < 6; i++) {
            atoms.push({ symbol: 'C', x: ring2Order[i].x, y: ring2Order[i].y });
        }
        // Ring2 bonds: shared1 -> ring2[2] ... ring2[5] -> shared2
        bonds.push({ a1: 1, a2: newAtomStart, type: DOUBLE });
        for (var i = 0; i < 3; i++) {
            bonds.push({ a1: newAtomStart + i, a2: newAtomStart + i + 1, type: i % 2 === 0 ? SINGLE : DOUBLE });
        }
        bonds.push({ a1: newAtomStart + 3, a2: 2, type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Pyridine (6-membered ring with one nitrogen)
    // =========================================================================
    Templates.pyridine = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Pyrimidine (6-membered ring with N at positions 1 and 3)
    // =========================================================================
    Templates.pyrimidine = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i === 0 || i === 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Indole (benzene fused to pyrrole = 6+5 fused ring)
    // =========================================================================
    Templates.indole = function() {
        // Start with benzene ring
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        // Fuse 5-membered ring sharing edge 2-3
        var shared1 = atoms[2];
        var shared2 = atoms[3];
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        var nx = -dy, ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        nx /= len; ny /= len;

        var radius5 = BL / (2 * Math.sin(Math.PI / 5));
        var apothem5 = radius5 * Math.cos(Math.PI / 5);
        var cx5 = mx + nx * apothem5;
        var cy5 = my + ny * apothem5;

        var r5 = polygon(5, cx5, cy5);
        var ring5Ordered = reorderToMatch(r5, shared1, shared2);

        var newStart = atoms.length;
        // Add 3 new atoms from the 5-ring (index 2 is NH)
        atoms.push({ symbol: 'C', x: ring5Ordered[2].x, y: ring5Ordered[2].y });
        atoms.push({ symbol: 'N', x: ring5Ordered[3].x, y: ring5Ordered[3].y }); // NH
        atoms.push({ symbol: 'C', x: ring5Ordered[4].x, y: ring5Ordered[4].y });

        // Bonds for 5-ring: 2 -> new0 -> new1(N) -> new2 -> 3
        bonds.push({ a1: 2, a2: newStart, type: DOUBLE });
        bonds.push({ a1: newStart, a2: newStart + 1, type: SINGLE });
        bonds.push({ a1: newStart + 1, a2: newStart + 2, type: SINGLE });
        bonds.push({ a1: newStart + 2, a2: 3, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Steroid skeleton (4 fused rings: A,B,C = 6-membered; D = 5-membered)
    // Gonane: cyclopenta[a]phenanthrene backbone
    // =========================================================================
    Templates.steroid = function() {
        // Cyclopenta[a]phenanthrene scaffold (gonane): A,B,C are 6-rings, D is
        // 5-ring, all fused linearly so the system extends left -> right.
        // Standard angular drawing convention used in textbooks.
        // Ring A (6): leftmost
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, false);

        // Ring B (6): shares edge 0-1 of ring A (bottom-right edge)
        var shared = [atoms[0], atoms[1]];
        var ring = fusedHexagon(shared[0], shared[1], atoms);
        var bStart = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ring[i].x, y: ring[i].y });
        }
        bonds.push({ a1: 0, a2: bStart, type: SINGLE });
        bonds.push({ a1: bStart, a2: bStart + 1, type: SINGLE });
        bonds.push({ a1: bStart + 1, a2: bStart + 2, type: SINGLE });
        bonds.push({ a1: bStart + 2, a2: bStart + 3, type: SINGLE });
        bonds.push({ a1: bStart + 3, a2: 1, type: SINGLE });

        // Ring C (6): shares edge bStart - (bStart+1) of ring B
        shared = [atoms[bStart], atoms[bStart + 1]];
        ring = fusedHexagon(shared[0], shared[1], atoms);
        var cStart = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ring[i].x, y: ring[i].y });
        }
        bonds.push({ a1: bStart, a2: cStart, type: SINGLE });
        bonds.push({ a1: cStart, a2: cStart + 1, type: SINGLE });
        bonds.push({ a1: cStart + 1, a2: cStart + 2, type: SINGLE });
        bonds.push({ a1: cStart + 2, a2: cStart + 3, type: SINGLE });
        bonds.push({ a1: cStart + 3, a2: bStart + 1, type: SINGLE });

        // Ring D (5): shares edge cStart - (cStart+1) of ring C
        shared = [atoms[cStart], atoms[cStart + 1]];
        ring = fusedPentagon(shared[0], shared[1], atoms);
        var dStart = atoms.length;
        for (var i = 0; i < 3; i++) {
            atoms.push({ symbol: 'C', x: ring[i].x, y: ring[i].y });
        }
        bonds.push({ a1: cStart, a2: dStart, type: SINGLE });
        bonds.push({ a1: dStart, a2: dStart + 1, type: SINGLE });
        bonds.push({ a1: dStart + 1, a2: dStart + 2, type: SINGLE });
        bonds.push({ a1: dStart + 2, a2: cStart + 1, type: SINGLE });

        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Amino acid backbone: N-Ca-C(=O) with R group
    // =========================================================================
    Templates.aminoAcid = function() {
        // Standard backbone layout: zigzag
        var sin60 = Math.sin(Math.PI / 3);
        var cos60 = Math.cos(Math.PI / 3);

        var atoms = [
            { symbol: 'N', x: 0, y: 0 },                                          // 0: amino N
            { symbol: 'C', x: BL, y: 0 },                                          // 1: C-alpha
            { symbol: 'C', x: BL + BL * cos60, y: BL * sin60 },                    // 2: carbonyl C
            { symbol: 'O', x: BL + BL * cos60 + BL, y: BL * sin60 },              // 3: carbonyl O
            { symbol: 'O', x: BL + BL * cos60, y: BL * sin60 + BL },              // 4: hydroxyl O (COOH)
            { symbol: 'C', x: BL + BL * cos60, y: -BL * sin60 }                   // 5: R group (side chain)
        ];

        var bonds = [
            { a1: 0, a2: 1, type: SINGLE },   // N - Ca
            { a1: 1, a2: 2, type: SINGLE },   // Ca - C
            { a1: 2, a2: 3, type: DOUBLE },   // C = O
            { a1: 2, a2: 4, type: SINGLE },   // C - OH
            { a1: 1, a2: 5, type: SINGLE }    // Ca - R
        ];

        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Purine (fused pyrimidine + imidazole, 9 atoms)
    // Numbering: 6-ring [0..5], 5-ring shares edge 3-4, adds [6,7,8]
    // =========================================================================
    Templates.purine = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p, i) {
            // Pyrimidine ring: N at positions 0 and 2
            return { symbol: (i === 0 || i === 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);

        // Fuse 5-ring sharing edge 3-4
        var shared1 = atoms[3];
        var shared2 = atoms[4];
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        var nx = -dy, ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        nx /= len; ny /= len;

        var radius5 = BL / (2 * Math.sin(Math.PI / 5));
        var apothem5 = radius5 * Math.cos(Math.PI / 5);
        var cx5 = mx + nx * apothem5;
        var cy5 = my + ny * apothem5;

        var r5 = polygon(5, cx5, cy5);
        var ring5Ordered = reorderToMatch(r5, shared1, shared2);

        var ns = atoms.length;
        // Imidazole portion: C, N, N
        atoms.push({ symbol: 'C', x: ring5Ordered[2].x, y: ring5Ordered[2].y }); // 6
        atoms.push({ symbol: 'N', x: ring5Ordered[3].x, y: ring5Ordered[3].y }); // 7 (N7)
        atoms.push({ symbol: 'C', x: ring5Ordered[4].x, y: ring5Ordered[4].y }); // 8 (C8)

        bonds.push({ a1: 3, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 4, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Quinoline (fused benzene + pyridine, 10 atoms)
    // =========================================================================
    Templates.quinoline = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        // Fuse second ring sharing edge 1-2
        var shared1 = atoms[1];
        var shared2 = atoms[2];
        var ring = fusedHexagon(shared1, shared2);
        var ns = atoms.length;
        // Pyridine: N at position 1 of second ring (index ns+1)
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'N', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });

        bonds.push({ a1: 1, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: DOUBLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 2, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Isoquinoline (fused benzene + pyridine, N at different position)
    // =========================================================================
    Templates.isoquinoline = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        var shared1 = atoms[1];
        var shared2 = atoms[2];
        var ring = fusedHexagon(shared1, shared2);
        var ns = atoms.length;
        // N at position 2 of second ring (index ns+2)
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'N', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });

        bonds.push({ a1: 1, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: DOUBLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 2, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Quinazoline (6+6 fused, N at 1 and 3 of second ring, 10 atoms)
    // =========================================================================
    Templates.quinazoline = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        var shared1 = atoms[1];
        var shared2 = atoms[2];
        var ring = fusedHexagon(shared1, shared2);
        var ns = atoms.length;
        // N at positions 1 and 3 of second ring
        atoms.push({ symbol: 'N', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'N', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });

        bonds.push({ a1: 1, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: DOUBLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: DOUBLE });
        bonds.push({ a1: ns + 3, a2: 2, type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Benzimidazole (benzene fused to imidazole, 9 atoms)
    // =========================================================================
    Templates.benzimidazole = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        // Fuse 5-ring sharing edge 2-3
        var shared1 = atoms[2];
        var shared2 = atoms[3];
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        var nx = -dy, ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        nx /= len; ny /= len;

        var radius5 = BL / (2 * Math.sin(Math.PI / 5));
        var apothem5 = radius5 * Math.cos(Math.PI / 5);
        var cx5 = mx + nx * apothem5;
        var cy5 = my + ny * apothem5;

        var r5 = polygon(5, cx5, cy5);
        var ring5Ordered = reorderToMatch(r5, shared1, shared2);

        var ns = atoms.length;
        atoms.push({ symbol: 'N', x: ring5Ordered[2].x, y: ring5Ordered[2].y }); // NH
        atoms.push({ symbol: 'C', x: ring5Ordered[3].x, y: ring5Ordered[3].y });
        atoms.push({ symbol: 'N', x: ring5Ordered[4].x, y: ring5Ordered[4].y });

        bonds.push({ a1: 2, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: DOUBLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 3, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Benzofuran (benzene fused to furan, 9 atoms)
    // =========================================================================
    Templates.benzofuran = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        var shared1 = atoms[2];
        var shared2 = atoms[3];
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        var nx = -dy, ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        nx /= len; ny /= len;

        var radius5 = BL / (2 * Math.sin(Math.PI / 5));
        var apothem5 = radius5 * Math.cos(Math.PI / 5);
        var cx5 = mx + nx * apothem5;
        var cy5 = my + ny * apothem5;

        var r5 = polygon(5, cx5, cy5);
        var ring5Ordered = reorderToMatch(r5, shared1, shared2);

        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring5Ordered[2].x, y: ring5Ordered[2].y });
        atoms.push({ symbol: 'O', x: ring5Ordered[3].x, y: ring5Ordered[3].y });
        atoms.push({ symbol: 'C', x: ring5Ordered[4].x, y: ring5Ordered[4].y });

        bonds.push({ a1: 2, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 3, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Benzothiazole (benzene fused to thiazole, 9 atoms)
    // =========================================================================
    Templates.benzothiazole = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        var shared1 = atoms[2];
        var shared2 = atoms[3];
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        var nx = -dy, ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        nx /= len; ny /= len;

        var radius5 = BL / (2 * Math.sin(Math.PI / 5));
        var apothem5 = radius5 * Math.cos(Math.PI / 5);
        var cx5 = mx + nx * apothem5;
        var cy5 = my + ny * apothem5;

        var r5 = polygon(5, cx5, cy5);
        var ring5Ordered = reorderToMatch(r5, shared1, shared2);

        var ns = atoms.length;
        atoms.push({ symbol: 'S', x: ring5Ordered[2].x, y: ring5Ordered[2].y });
        atoms.push({ symbol: 'C', x: ring5Ordered[3].x, y: ring5Ordered[3].y });
        atoms.push({ symbol: 'N', x: ring5Ordered[4].x, y: ring5Ordered[4].y });

        bonds.push({ a1: 2, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: DOUBLE });
        bonds.push({ a1: ns + 2, a2: 3, type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Pteridine (6+6 fused, 4 nitrogens: pyrazine + pyrimidine, 10 atoms)
    // =========================================================================
    Templates.pteridine = function() {
        var r6 = polygon(6, 0, 0);
        // Pyrimidine ring: N at 0 and 2
        var atoms = r6.map(function(p, i) {
            return { symbol: (i === 0 || i === 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);

        // Fuse pyrazine ring sharing edge 3-4, N at positions 1 and 2
        var shared1 = atoms[3];
        var shared2 = atoms[4];
        var ring = fusedHexagon(shared1, shared2);
        var ns = atoms.length;
        atoms.push({ symbol: 'N', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'N', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });

        bonds.push({ a1: 3, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: DOUBLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: DOUBLE });
        bonds.push({ a1: ns + 3, a2: 4, type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // =========================================================================
    // Phenothiazine (tricyclic: 6+6+6, S and N bridges, 13 atoms)
    // =========================================================================
    Templates.phenothiazine = function() {
        // Central 6-ring with S at 0, N at 3
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p, i) {
            if (i === 0) return { symbol: 'S', x: p.x, y: p.y };
            if (i === 3) return { symbol: 'N', x: p.x, y: p.y };
            return { symbol: 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);

        // Left ring fused at edge 4-5
        var ringL = fusedHexagon(atoms[4], atoms[5]);
        var ls = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringL[i].x, y: ringL[i].y });
        }
        bonds.push({ a1: 4, a2: ls, type: DOUBLE });
        bonds.push({ a1: ls, a2: ls + 1, type: SINGLE });
        bonds.push({ a1: ls + 1, a2: ls + 2, type: DOUBLE });
        bonds.push({ a1: ls + 2, a2: ls + 3, type: SINGLE });
        bonds.push({ a1: ls + 3, a2: 5, type: DOUBLE });

        // Right ring fused at edge 1-2
        var ringR = fusedHexagon(atoms[1], atoms[2]);
        var rs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringR[i].x, y: ringR[i].y });
        }
        bonds.push({ a1: 1, a2: rs, type: DOUBLE });
        bonds.push({ a1: rs, a2: rs + 1, type: SINGLE });
        bonds.push({ a1: rs + 1, a2: rs + 2, type: DOUBLE });
        bonds.push({ a1: rs + 2, a2: rs + 3, type: SINGLE });
        bonds.push({ a1: rs + 3, a2: 2, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // =========================================================================
    // Acridine (tricyclic 6+6+6, one N in central ring, 13 atoms)
    // =========================================================================
    Templates.acridine = function() {
        // Central ring with N at position 0
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);

        // Left ring fused at edge 4-5
        var ringL = fusedHexagon(atoms[4], atoms[5]);
        var ls = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringL[i].x, y: ringL[i].y });
        }
        bonds.push({ a1: 4, a2: ls, type: DOUBLE });
        bonds.push({ a1: ls, a2: ls + 1, type: SINGLE });
        bonds.push({ a1: ls + 1, a2: ls + 2, type: DOUBLE });
        bonds.push({ a1: ls + 2, a2: ls + 3, type: SINGLE });
        bonds.push({ a1: ls + 3, a2: 5, type: DOUBLE });

        // Right ring fused at edge 1-2
        var ringR = fusedHexagon(atoms[1], atoms[2]);
        var rs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringR[i].x, y: ringR[i].y });
        }
        bonds.push({ a1: 1, a2: rs, type: DOUBLE });
        bonds.push({ a1: rs, a2: rs + 1, type: SINGLE });
        bonds.push({ a1: rs + 1, a2: rs + 2, type: DOUBLE });
        bonds.push({ a1: rs + 2, a2: rs + 3, type: SINGLE });
        bonds.push({ a1: rs + 3, a2: 2, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // =========================================================================
    // Phenanthrene (tricyclic 6+6+6 angular fusion, 14 atoms)
    // =========================================================================
    Templates.phenanthrene = function() {
        // Ring A
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);

        // Ring B fused at edge 1-2 (linear)
        var ringB = fusedHexagon(atoms[1], atoms[2]);
        var bs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringB[i].x, y: ringB[i].y });
        }
        bonds.push({ a1: 1, a2: bs, type: SINGLE });
        bonds.push({ a1: bs, a2: bs + 1, type: DOUBLE });
        bonds.push({ a1: bs + 1, a2: bs + 2, type: SINGLE });
        bonds.push({ a1: bs + 2, a2: bs + 3, type: DOUBLE });
        bonds.push({ a1: bs + 3, a2: 2, type: SINGLE });

        // Ring C fused at edge bs+1 - bs+2 (angular)
        var ringC = fusedHexagon(atoms[bs + 1], atoms[bs + 2]);
        var cs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringC[i].x, y: ringC[i].y });
        }
        bonds.push({ a1: bs + 1, a2: cs, type: DOUBLE });
        bonds.push({ a1: cs, a2: cs + 1, type: SINGLE });
        bonds.push({ a1: cs + 1, a2: cs + 2, type: DOUBLE });
        bonds.push({ a1: cs + 2, a2: cs + 3, type: SINGLE });
        bonds.push({ a1: cs + 3, a2: bs + 2, type: DOUBLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // =========================================================================
    // Carbazole (tricyclic: benzene + pyrrole + benzene, 13 atoms)
    // =========================================================================
    Templates.carbazole = function() {
        // Central 5-ring (pyrrole) with N at top
        var r5 = polygon(5, 0, 0);
        var atoms = r5.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE; // N-C single
        bonds[1].type = DOUBLE;
        bonds[2].type = SINGLE;
        bonds[3].type = DOUBLE;
        bonds[4].type = SINGLE; // C-N single

        // Left benzene fused at edge 3-4
        var ringL = fusedHexagon(atoms[3], atoms[4]);
        var ls = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringL[i].x, y: ringL[i].y });
        }
        bonds.push({ a1: 3, a2: ls, type: SINGLE });
        bonds.push({ a1: ls, a2: ls + 1, type: DOUBLE });
        bonds.push({ a1: ls + 1, a2: ls + 2, type: SINGLE });
        bonds.push({ a1: ls + 2, a2: ls + 3, type: DOUBLE });
        bonds.push({ a1: ls + 3, a2: 4, type: SINGLE });

        // Right benzene fused at edge 1-2
        var ringR = fusedHexagon(atoms[1], atoms[2]);
        var rs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringR[i].x, y: ringR[i].y });
        }
        bonds.push({ a1: 1, a2: rs, type: SINGLE });
        bonds.push({ a1: rs, a2: rs + 1, type: DOUBLE });
        bonds.push({ a1: rs + 1, a2: rs + 2, type: SINGLE });
        bonds.push({ a1: rs + 2, a2: rs + 3, type: DOUBLE });
        bonds.push({ a1: rs + 3, a2: 2, type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // =========================================================================
    // Pyrrole (5-membered, 1 N)
    // =========================================================================
    Templates.pyrrole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, true);
        bonds[0].type = SINGLE; // N-C single
        bonds[4].type = SINGLE; // C-N single
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Imidazole (5-membered, N at 0 and 2)
    // =========================================================================
    Templates.imidazole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i === 0 || i === 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE;
        bonds[1].type = DOUBLE;
        bonds[2].type = SINGLE;
        bonds[3].type = SINGLE;
        bonds[4].type = DOUBLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Pyrazole (5-membered, N at 0 and 1)
    // =========================================================================
    Templates.pyrazole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i === 0 || i === 1) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE; // N-N
        bonds[1].type = DOUBLE;
        bonds[2].type = SINGLE;
        bonds[3].type = DOUBLE;
        bonds[4].type = SINGLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Thiazole (5-membered, S at 0, N at 2)
    // =========================================================================
    Templates.thiazole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            if (i === 0) return { symbol: 'S', x: p.x, y: p.y };
            if (i === 2) return { symbol: 'N', x: p.x, y: p.y };
            return { symbol: 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE;
        bonds[1].type = DOUBLE;
        bonds[2].type = SINGLE;
        bonds[3].type = DOUBLE;
        bonds[4].type = SINGLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Thiophene (5-membered, S at 0)
    // =========================================================================
    Templates.thiophene = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'S' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, true);
        bonds[0].type = SINGLE; // S-C single
        bonds[4].type = SINGLE; // C-S single
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Furan (5-membered, O at 0)
    // =========================================================================
    Templates.furan = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'O' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, true);
        bonds[0].type = SINGLE; // O-C single
        bonds[4].type = SINGLE; // C-O single
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // 1,2,3-Triazole (5-membered, N at 0, 1, 2)
    // =========================================================================
    Templates.triazole_1_2_3 = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i <= 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE; // N-N
        bonds[1].type = DOUBLE; // N=N
        bonds[2].type = SINGLE;
        bonds[3].type = DOUBLE;
        bonds[4].type = SINGLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // 1,2,4-Triazole (5-membered, N at 0, 1, 3)
    // =========================================================================
    Templates.triazole_1_2_4 = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i === 0 || i === 1 || i === 3) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE; // N-N
        bonds[1].type = DOUBLE;
        bonds[2].type = SINGLE;
        bonds[3].type = SINGLE;
        bonds[4].type = DOUBLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Tetrazole (5-membered, N at 0, 1, 2, 3)
    // =========================================================================
    Templates.tetrazole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i <= 3 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[0].type = SINGLE;
        bonds[1].type = DOUBLE;
        bonds[2].type = SINGLE;
        bonds[3].type = DOUBLE;
        bonds[4].type = SINGLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Piperidine (saturated 6-membered, 1 N)
    // =========================================================================
    Templates.piperidine = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Piperazine (saturated 6-membered, N at 0 and 3)
    // =========================================================================
    Templates.piperazine = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i === 0 || i === 3) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Morpholine (saturated 6-membered, N at 0, O at 3)
    // =========================================================================
    Templates.morpholine = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p, i) {
            if (i === 0) return { symbol: 'N', x: p.x, y: p.y };
            if (i === 3) return { symbol: 'O', x: p.x, y: p.y };
            return { symbol: 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Pyranose (6-membered sugar ring with one O, Haworth-style)
    // O placed at the upper-right vertex; C1 (anomeric) at the right;
    // remaining C2-C5 walking around the ring. For glucose/galactose/etc.
    // =========================================================================
    Templates.pyranose = function() {
        var pts = polygon(6, 0, 0);
        // Standard polygon(6) gives vertex 0 at top-right area for even-size
        // (after offset += step/2). We pick index 0 as the ring O
        // (top-right Haworth convention) and then C1, C2, C3, C4, C5 around.
        var atoms = pts.map(function(p, i) {
            return {
                symbol: i === 0 ? 'O' : 'C',
                x: p.x,
                y: p.y
            };
        });
        var bonds = ringBonds(6, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Furanose (5-membered sugar ring with one O, Haworth-style)
    // For ribose/deoxyribose. O at upper-right vertex.
    // =========================================================================
    Templates.furanose = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return {
                symbol: i === 0 ? 'O' : 'C',
                x: p.x,
                y: p.y
            };
        });
        var bonds = ringBonds(5, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Tropane (8-aza-bicyclo[3.2.1]octane) — bridged bicyclic
    // Used for atropine/cocaine/scopolamine. 7-ring fused with 5-ring
    // sharing 3 atoms (two bridgeheads + bridge N).
    // 8 atoms total: 2 bridgeheads + 5 ring atoms + 1 bridging atom
    // We model it as: cyclohexane (6-ring) with N-bridge across C1-C5.
    // =========================================================================
    Templates.tropane = function() {
        // Place a 6-ring (cyclohexane chair-like flat) plus a single
        // bridging atom on the upper side connecting two non-adjacent atoms.
        var pts = polygon(7, 0, 0);
        // Use a 7-ring for layout simplicity: tropane is conventionally
        // drawn as a cyclohexane with the N bridging atoms 1 and 5; we
        // fake the bridge with a single atom placed above the ring centre.
        // Actually use 8 atoms: 6-ring + bridging N + the bridge C-N.
        var hexPts = polygon(6, 0, 0);
        var atoms = hexPts.map(function(p) {
            return { symbol: 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);

        // Bridge: N atom positioned above the midpoint of the ring,
        // bonded to atoms 1 and 4 (across the ring). Bridge atoms get a
        // perpendicular bow above the ring plane.
        var bridgeN = atoms.length;
        var a1 = atoms[1], a4 = atoms[4];
        var midX = (a1.x + a4.x) / 2;
        var midY = (a1.y + a4.y) / 2;
        atoms.push({ symbol: 'N', x: midX, y: midY - BL * 0.7 });
        bonds.push({ a1: 1, a2: bridgeN, type: SINGLE });
        bonds.push({ a1: 4, a2: bridgeN, type: SINGLE });

        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Beta-lactam (4-membered ring, azetidinone core)
    // =========================================================================
    Templates.beta_lactam = function() {
        var pts = polygon(4, 0, 0);
        var atoms = [
            { symbol: 'N', x: pts[0].x, y: pts[0].y },
            { symbol: 'C', x: pts[1].x, y: pts[1].y },
            { symbol: 'C', x: pts[2].x, y: pts[2].y },
            { symbol: 'C', x: pts[3].x, y: pts[3].y }
        ];
        var bonds = ringBonds(4, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Cyclopropane (3-membered ring)
    // =========================================================================
    Templates.cyclopropane = function() {
        var pts = polygon(3, 0, 0);
        var atoms = pts.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(3, 0, false);
        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Adamantane (cage structure, 10 atoms, 4 fused chair-form 6-rings)
    // Pre-computed 2D projection of the diamond-lattice cage
    // =========================================================================
    Templates.adamantane = function() {
        // Classic 2D adamantane layout: three radial spokes from a central triangle
        var s3 = Math.sqrt(3);
        var h = BL * s3 / 2;
        var atoms = [
            // Outer top bridgehead
            { symbol: 'C', x: 0,          y: -BL - h },          // 0
            // Upper triangle
            { symbol: 'C', x: -BL / 2,    y: -h },               // 1
            { symbol: 'C', x: BL / 2,     y: -h },               // 2
            // Middle ring bridgeheads
            { symbol: 'C', x: -BL,        y: 0 },                // 3
            { symbol: 'C', x: BL,         y: 0 },                // 4
            // Lower triangle
            { symbol: 'C', x: -BL / 2,    y: h },                // 5
            { symbol: 'C', x: BL / 2,     y: h },                // 6
            // Outer bottom bridgehead
            { symbol: 'C', x: 0,          y: BL + h },           // 7
            // Central cross atoms
            { symbol: 'C', x: 0,          y: 0 },                // 8
            { symbol: 'C', x: 0,          y: -BL * 0.4 }         // 9
        ];

        var bonds = [
            { a1: 0, a2: 1, type: SINGLE },
            { a1: 0, a2: 2, type: SINGLE },
            { a1: 1, a2: 3, type: SINGLE },
            { a1: 2, a2: 4, type: SINGLE },
            { a1: 1, a2: 8, type: SINGLE },
            { a1: 2, a2: 8, type: SINGLE },
            { a1: 3, a2: 5, type: SINGLE },
            { a1: 4, a2: 6, type: SINGLE },
            { a1: 5, a2: 7, type: SINGLE },
            { a1: 6, a2: 7, type: SINGLE },
            { a1: 5, a2: 8, type: SINGLE },
            { a1: 6, a2: 8, type: SINGLE },
            { a1: 3, a2: 9, type: SINGLE },
            { a1: 4, a2: 9, type: SINGLE },
            { a1: 0, a2: 9, type: SINGLE }
        ];

        return { atoms: atoms, bonds: bonds };
    };

    // =========================================================================
    // Chromone (4H-chromen-4-one core; flavone / flavanone / chroman scaffold)
    // 10 ring atoms: benzene fused to a pyran (one ring O).
    // Indices:
    //   0: C8a (top of fusion edge, shared with ring B)
    //   1..4: C8, C7, C6, C5 (benzene ring back to bottom of fusion edge)
    //   5: C4a (bottom of fusion edge, shared with ring B)
    //   6: O1 (ring oxygen, top-right of pyran)
    //   7: C2 (the C2-aryl attachment in flavanones / isoflavone C3-aryl)
    //   8: C3
    //   9: C4 (the carbonyl carbon in flavones; saturated in chroman)
    // Drawn with benzene on the LEFT, pyran on the RIGHT, ring O at top —
    // standard ChemDraw / Wikipedia depiction of flavanones / flavones.
    // The exocyclic C4=O ketone is *not* part of the template; the chain
    // layout adds it as a substituent so the same coordinates serve
    // chroman, chromone, flavanone and flavone scaffolds.
    // =========================================================================
    Templates.chromone = function() {
        // Ring A: aromatic benzene (left ring), default polygon orientation.
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true); // aromatic benzene

        // Ring B (pyran/pyranone): fuse on edge 1-2 (right side of benzene),
        // matching the naphthalene-style horizontal fusion. New atoms placed
        // to the right. The shared edge ends are atoms[1] (top) and atoms[2]
        // (bottom). After fusedHexagon + reorderToMatch:
        //   ring2[0] == atoms[1] (shared, C8a)
        //   ring2[1] == atoms[2] (shared, C4a)
        //   ring2[2..5] are the 4 new vertices traversing the right ring
        //   from C8a back around to C4a. Standard convention: O1 sits
        //   adjacent to C8a (i.e. ring2[2] = O1, ring2[3] = C2,
        //   ring2[4] = C3, ring2[5] = C4).
        var ringB = fusedHexagon(atoms[1], atoms[2], atoms);
        // fusedHexagon returns the 4 new vertices walking from the vertex
        // adjacent to atoms[2] (C4a side) around to the vertex adjacent to
        // atoms[1] (C8a side). So:
        //   ringB[0] = vertex adjacent to C4a → C4
        //   ringB[1] = next around → C3
        //   ringB[2] = next around → C2
        //   ringB[3] = vertex adjacent to C8a → O1
        var nb = atoms.length;
        atoms.push({ symbol: 'C', x: ringB[0].x, y: ringB[0].y }); // 6: C4
        atoms.push({ symbol: 'C', x: ringB[1].x, y: ringB[1].y }); // 7: C3
        atoms.push({ symbol: 'C', x: ringB[2].x, y: ringB[2].y }); // 8: C2
        atoms.push({ symbol: 'O', x: ringB[3].x, y: ringB[3].y }); // 9: O1

        // Pyran/pyranone bonds (single bonds — bond orders are taken from the
        // parsed molecule, the template only seeds coordinates):
        //   C4a (2) — C4 (6) — C3 (7) — C2 (8) — O1 (9) — C8a (1)
        bonds.push({ a1: 2,      a2: nb,     type: SINGLE });
        bonds.push({ a1: nb,     a2: nb + 1, type: SINGLE });
        bonds.push({ a1: nb + 1, a2: nb + 2, type: SINGLE });
        bonds.push({ a1: nb + 2, a2: nb + 3, type: SINGLE });
        bonds.push({ a1: nb + 3, a2: 1,      type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // =========================================================================
    // Coumarin (2H-chromen-2-one): benzene fused to alpha-pyrone, C2 is the
    // carbonyl carbon (different position vs chromone). Same 2:6,6 signature
    // and 9C+1O ring composition. We expose this as a separate template name
    // so the scoring can pick the better atom-position fit when needed; in
    // practice the coordinates are equivalent — the carbonyl is exocyclic
    // and added by the chain layout.
    // =========================================================================
    Templates.coumarin = function() {
        return Templates.chromone();
    };

    // =========================================================================
    // Morphinan (5-ring system: phenanthrene + piperidine + oxide bridge)
    // Classic morphine/codeine scaffold, 17 atoms
    // Rings: A(aromatic 6) + B(6) + C(6) + D(piperidine 6) + E(oxide bridge 5)
    // =========================================================================
    Templates.morphinan = function() {
        // Phenanthrene-derived 4-ring skeleton (A,B,C are 6, A aromatic),
        // plus piperidine ring D fused on top of ring C, plus the ether
        // E ring (5-membered) bridging A and C through an oxygen. This gives
        // the classic morphine/codeine display: A flat at left, B/C
        // descending stair-step, D rising back up, E spanning the front.
        // 17 ring atoms + 1 N-methyl substituent are common (we leave the
        // N-methyl for the chain layout to add since substituents vary).

        // Ring A: aromatic benzene (leftmost)
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true); // aromatic

        // Ring B (cyclohexene): fused at edge 0-1 of A (lower right edge)
        var ringB = fusedHexagon(atoms[0], atoms[1], atoms);
        var bStart = atoms.length; // 6
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringB[i].x, y: ringB[i].y });
        }
        // Bonds: A0 -- B0 -- B1 -- B2 = B3 -- A1
        bonds.push({ a1: 0, a2: bStart, type: SINGLE });
        bonds.push({ a1: bStart, a2: bStart + 1, type: SINGLE });
        bonds.push({ a1: bStart + 1, a2: bStart + 2, type: SINGLE });
        bonds.push({ a1: bStart + 2, a2: bStart + 3, type: DOUBLE });
        bonds.push({ a1: bStart + 3, a2: 1, type: SINGLE });

        // Ring C (cyclohexane): fused at edge bStart-(bStart+1) of B
        var ringC = fusedHexagon(atoms[bStart], atoms[bStart + 1], atoms);
        var cStart = atoms.length; // 10
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringC[i].x, y: ringC[i].y });
        }
        bonds.push({ a1: bStart, a2: cStart, type: SINGLE });
        bonds.push({ a1: cStart, a2: cStart + 1, type: SINGLE });
        bonds.push({ a1: cStart + 1, a2: cStart + 2, type: SINGLE });
        bonds.push({ a1: cStart + 2, a2: cStart + 3, type: SINGLE });
        bonds.push({ a1: cStart + 3, a2: bStart + 1, type: SINGLE });

        // Ring D (piperidine): fused at edge cStart-(cStart+1) of C
        var ringD = fusedHexagon(atoms[cStart], atoms[cStart + 1], atoms);
        var dStart = atoms.length; // 14
        atoms.push({ symbol: 'C', x: ringD[0].x, y: ringD[0].y });
        atoms.push({ symbol: 'C', x: ringD[1].x, y: ringD[1].y });
        atoms.push({ symbol: 'N', x: ringD[2].x, y: ringD[2].y }); // piperidine N
        atoms.push({ symbol: 'C', x: ringD[3].x, y: ringD[3].y });

        bonds.push({ a1: cStart, a2: dStart, type: SINGLE });
        bonds.push({ a1: dStart, a2: dStart + 1, type: SINGLE });
        bonds.push({ a1: dStart + 1, a2: dStart + 2, type: SINGLE });
        bonds.push({ a1: dStart + 2, a2: dStart + 3, type: SINGLE });
        bonds.push({ a1: dStart + 3, a2: cStart + 1, type: SINGLE });

        // Ring E (ether bridge): 5-membered O-bearing ring connecting
        // ring A position 0 (lower-right) to ring C atom cStart+3 (front)
        // through an oxygen. This is the dihydrofuran fused under rings
        // A and C in morphine/codeine. Position O so two bonds (~30 px) form.
        var oIdx = atoms.length; // 18
        var aA = atoms[0];           // ring A attachment
        var aC = atoms[cStart + 3];  // ring C attachment closest to A
        var midOx = (aA.x + aC.x) / 2;
        var midOy = (aA.y + aC.y) / 2;
        // Pull O outward (away from the molecule centroid) by a small
        // offset so the E ring forms a clean envelope below rings A-B-C.
        var cxAll = 0, cyAll = 0;
        for (var i = 0; i < atoms.length; i++) { cxAll += atoms[i].x; cyAll += atoms[i].y; }
        cxAll /= atoms.length; cyAll /= atoms.length;
        var ox = midOx - cxAll, oy = midOy - cyAll;
        var oLen = Math.sqrt(ox * ox + oy * oy) || 1;
        atoms.push({
            symbol: 'O',
            x: midOx + (ox / oLen) * BL * 0.4,
            y: midOy + (oy / oLen) * BL * 0.4
        });
        bonds.push({ a1: 0, a2: oIdx, type: SINGLE });
        bonds.push({ a1: oIdx, a2: cStart + 3, type: SINGLE });

        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };


    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 1: simple gap-fills) =====
    Templates.pyrrolidine = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        return { atoms: atoms, bonds: ringBonds(5, 0, false) };
    };

    Templates.cyclobutane = function() {
        var pts = polygon(4, 0, 0);
        var atoms = pts.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        return { atoms: atoms, bonds: ringBonds(4, 0, false) };
    };

    Templates.cycloheptane = function() {
        var pts = polygon(7, 0, 0);
        var atoms = pts.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        return { atoms: atoms, bonds: ringBonds(7, 0, false) };
    };

    Templates.tetrahydropyran = function() {
        var pts = polygon(6, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i === 0 ? 'O' : 'C', x: p.x, y: p.y };
        });
        return { atoms: atoms, bonds: ringBonds(6, 0, false) };
    };

    Templates.oxazole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            if (i === 0) return { symbol: 'O', x: p.x, y: p.y };
            if (i === 2) return { symbol: 'N', x: p.x, y: p.y };
            return { symbol: 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[1].type = DOUBLE; bonds[3].type = DOUBLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    Templates.isoxazole = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            if (i === 0) return { symbol: 'O', x: p.x, y: p.y };
            if (i === 1) return { symbol: 'N', x: p.x, y: p.y };
            return { symbol: 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[1].type = DOUBLE; bonds[3].type = DOUBLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    Templates.triazole123 = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: i <= 2 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[1].type = DOUBLE; bonds[3].type = DOUBLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    Templates.triazole124 = function() {
        var pts = polygon(5, 0, 0);
        var atoms = pts.map(function(p, i) {
            return { symbol: (i === 0 || i === 1 || i === 3) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        bonds[1].type = DOUBLE; bonds[4].type = DOUBLE;
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    Templates.betalactam = function() {
        var pts = polygon(4, 0, 0);
        var atoms = [
            { symbol: 'N', x: pts[0].x, y: pts[0].y },
            { symbol: 'C', x: pts[1].x, y: pts[1].y },
            { symbol: 'C', x: pts[2].x, y: pts[2].y },
            { symbol: 'C', x: pts[3].x, y: pts[3].y }
        ];
        return { atoms: atoms, bonds: ringBonds(4, 0, false) };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 2: fused heterocycles) =====

    // Xanthine — 2,6-dioxopurine; same ring topology as purine.
    Templates.xanthine = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p, i) {
            return { symbol: (i === 0 || i === 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedPentagon(atoms[3], atoms[4]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'N', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 3, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 4, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // 7-Deazapurine — pyrimidine fused to pyrrole.
    Templates.pyrrolopyrimidine = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p, i) {
            return { symbol: (i === 0 || i === 2) ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedPentagon(atoms[3], atoms[4]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'N', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 3, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 4, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // 1,4-Benzodiazepine — benzene + 7-ring with two N (diazepam scaffold).
    // We compute the 7-ring centre manually because the existing helpers
    // fusedHexagon / fusedPentagon only handle 6- and 5-rings; the geometry
    // follows the same midpoint+normal construction generalised to 7.
    Templates.benzodiazepine = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var s1 = atoms[1], s2 = atoms[2];
        var mx = (s1.x + s2.x) / 2;
        var my = (s1.y + s2.y) / 2;
        var dx = s2.x - s1.x, dy = s2.y - s1.y;
        var nx = -dy, ny = dx;
        var nlen = Math.sqrt(nx * nx + ny * ny) || 1;
        nx /= nlen; ny /= nlen;
        var radius7 = BL / (2 * Math.sin(Math.PI / 7));
        var apothem7 = radius7 * Math.cos(Math.PI / 7);
        var sign = (nx * mx + ny * my) >= 0 ? 1 : -1;
        var cx7 = mx + sign * nx * apothem7;
        var cy7 = my + sign * ny * apothem7;
        var pts7 = polygon(7, cx7, cy7);
        var ordered = reorderToMatch(pts7, s1, s2);
        var ns = atoms.length;
        atoms.push({ symbol: 'N', x: ordered[2].x, y: ordered[2].y });
        atoms.push({ symbol: 'C', x: ordered[3].x, y: ordered[3].y });
        atoms.push({ symbol: 'C', x: ordered[4].x, y: ordered[4].y });
        atoms.push({ symbol: 'N', x: ordered[5].x, y: ordered[5].y });
        atoms.push({ symbol: 'C', x: ordered[6].x, y: ordered[6].y });
        bonds.push({ a1: 2, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: DOUBLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: ns + 4, type: SINGLE });
        bonds.push({ a1: ns + 4, a2: 1, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    // Benzothiophene — benzene fused to thiophene.
    Templates.benzothiophene = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedPentagon(atoms[2], atoms[3]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'S', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 2, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 3, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 3: quinoline family) =====

    Templates.quinoxaline = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedHexagon(atoms[1], atoms[2]);
        var ns = atoms.length;
        atoms.push({ symbol: 'N', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'N', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: 1, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: DOUBLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: DOUBLE });
        bonds.push({ a1: ns + 3, a2: 2, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    Templates.cinnoline = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedHexagon(atoms[1], atoms[2]);
        var ns = atoms.length;
        atoms.push({ symbol: 'N', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'N', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: 1, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: DOUBLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 2, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    Templates.phthalazine = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedHexagon(atoms[1], atoms[2]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'N', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'N', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: 1, a2: ns, type: DOUBLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: DOUBLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 2, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    // β-Carboline — pyrido[3,4-b]indole — fused 6+5+6 with N in indole and pyridine.
    Templates.beta_carboline = function() {
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ringB = fusedPentagon(atoms[2], atoms[3]);
        var bs = atoms.length;
        atoms.push({ symbol: 'C', x: ringB[0].x, y: ringB[0].y });
        atoms.push({ symbol: 'N', x: ringB[1].x, y: ringB[1].y });
        atoms.push({ symbol: 'C', x: ringB[2].x, y: ringB[2].y });
        bonds.push({ a1: 2, a2: bs, type: DOUBLE });
        bonds.push({ a1: bs, a2: bs + 1, type: SINGLE });
        bonds.push({ a1: bs + 1, a2: bs + 2, type: SINGLE });
        bonds.push({ a1: bs + 2, a2: 3, type: DOUBLE });
        var ring = fusedHexagon(atoms[bs], atoms[bs + 2]);
        var cs = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'N', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: bs, a2: cs, type: SINGLE });
        bonds.push({ a1: cs, a2: cs + 1, type: DOUBLE });
        bonds.push({ a1: cs + 1, a2: cs + 2, type: SINGLE });
        bonds.push({ a1: cs + 2, a2: cs + 3, type: DOUBLE });
        bonds.push({ a1: cs + 3, a2: bs + 2, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 4: polycyclic aromatics) =====

    Templates.anthracene = function() {
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ringB = fusedHexagon(atoms[1], atoms[2]);
        var bs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringB[i].x, y: ringB[i].y });
        }
        bonds.push({ a1: 1, a2: bs, type: SINGLE });
        bonds.push({ a1: bs, a2: bs + 1, type: DOUBLE });
        bonds.push({ a1: bs + 1, a2: bs + 2, type: SINGLE });
        bonds.push({ a1: bs + 2, a2: bs + 3, type: DOUBLE });
        bonds.push({ a1: bs + 3, a2: 2, type: SINGLE });
        var ringC = fusedHexagon(atoms[bs + 1], atoms[bs + 2]);
        var cs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringC[i].x, y: ringC[i].y });
        }
        bonds.push({ a1: bs + 1, a2: cs, type: SINGLE });
        bonds.push({ a1: cs, a2: cs + 1, type: DOUBLE });
        bonds.push({ a1: cs + 1, a2: cs + 2, type: SINGLE });
        bonds.push({ a1: cs + 2, a2: cs + 3, type: DOUBLE });
        bonds.push({ a1: cs + 3, a2: bs + 2, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    Templates.pyrene = function() {
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ringB = fusedHexagon(atoms[1], atoms[2]);
        var bs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringB[i].x, y: ringB[i].y });
        }
        bonds.push({ a1: 1, a2: bs, type: SINGLE });
        bonds.push({ a1: bs, a2: bs + 1, type: DOUBLE });
        bonds.push({ a1: bs + 1, a2: bs + 2, type: SINGLE });
        bonds.push({ a1: bs + 2, a2: bs + 3, type: DOUBLE });
        bonds.push({ a1: bs + 3, a2: 2, type: SINGLE });
        var ringC = fusedHexagon(atoms[bs + 2], atoms[bs + 3]);
        var cs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringC[i].x, y: ringC[i].y });
        }
        bonds.push({ a1: bs + 2, a2: cs, type: SINGLE });
        bonds.push({ a1: cs, a2: cs + 1, type: DOUBLE });
        bonds.push({ a1: cs + 1, a2: cs + 2, type: SINGLE });
        bonds.push({ a1: cs + 2, a2: cs + 3, type: DOUBLE });
        bonds.push({ a1: cs + 3, a2: bs + 3, type: SINGLE });
        var ringD = fusedHexagon(atoms[cs + 2], atoms[cs + 3]);
        var ds = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringD[i].x, y: ringD[i].y });
        }
        bonds.push({ a1: cs + 2, a2: ds, type: SINGLE });
        bonds.push({ a1: ds, a2: ds + 1, type: DOUBLE });
        bonds.push({ a1: ds + 1, a2: ds + 2, type: SINGLE });
        bonds.push({ a1: ds + 2, a2: ds + 3, type: DOUBLE });
        bonds.push({ a1: ds + 3, a2: cs + 3, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2, 3] };
    };

    Templates.fluorene = function() {
        var r5 = polygon(5, 0, 0);
        var atoms = r5.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(5, 0, false);
        var ringL = fusedHexagon(atoms[3], atoms[4]);
        var ls = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringL[i].x, y: ringL[i].y });
        }
        bonds.push({ a1: 3, a2: ls, type: DOUBLE });
        bonds.push({ a1: ls, a2: ls + 1, type: SINGLE });
        bonds.push({ a1: ls + 1, a2: ls + 2, type: DOUBLE });
        bonds.push({ a1: ls + 2, a2: ls + 3, type: SINGLE });
        bonds.push({ a1: ls + 3, a2: 4, type: DOUBLE });
        var ringR = fusedHexagon(atoms[1], atoms[2]);
        var rs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringR[i].x, y: ringR[i].y });
        }
        bonds.push({ a1: 1, a2: rs, type: DOUBLE });
        bonds.push({ a1: rs, a2: rs + 1, type: SINGLE });
        bonds.push({ a1: rs + 1, a2: rs + 2, type: DOUBLE });
        bonds.push({ a1: rs + 2, a2: rs + 3, type: SINGLE });
        bonds.push({ a1: rs + 3, a2: 2, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [1, 2] };
    };

    Templates.biphenyl = function() {
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        // Hex radius is BL/(2 sin(π/6)) = BL. Two hex rings touch at radius+BL+radius
        // = 3 BL between centres if we want a SINGLE bond bridging the two
        // closest carbons.
        var rB = polygon(6, BL * 3, 0);
        var bs = atoms.length;
        for (var i = 0; i < 6; i++) {
            atoms.push({ symbol: 'C', x: rB[i].x, y: rB[i].y });
        }
        for (var i = 0; i < 6; i++) {
            bonds.push({ a1: bs + i, a2: bs + (i + 1) % 6, type: i % 2 === 0 ? DOUBLE : SINGLE });
        }
        bonds.push({ a1: 1, a2: bs + 4, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1] };
    };

    Templates.terphenyl = function() {
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var rB = polygon(6, BL * 3, 0);
        var bs = atoms.length;
        for (var i = 0; i < 6; i++) {
            atoms.push({ symbol: 'C', x: rB[i].x, y: rB[i].y });
        }
        for (var i = 0; i < 6; i++) {
            bonds.push({ a1: bs + i, a2: bs + (i + 1) % 6, type: i % 2 === 0 ? DOUBLE : SINGLE });
        }
        bonds.push({ a1: 1, a2: bs + 4, type: SINGLE });
        var rC = polygon(6, BL * 6, 0);
        var cs = atoms.length;
        for (var i = 0; i < 6; i++) {
            atoms.push({ symbol: 'C', x: rC[i].x, y: rC[i].y });
        }
        for (var i = 0; i < 6; i++) {
            bonds.push({ a1: cs + i, a2: cs + (i + 1) % 6, type: i % 2 === 0 ? DOUBLE : SINGLE });
        }
        bonds.push({ a1: bs + 1, a2: cs + 4, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 5: β-lactam cores) =====

    Templates.penam = function() {
        var p4 = polygon(4, 0, 0);
        var atoms = [
            { symbol: 'N', x: p4[0].x, y: p4[0].y },
            { symbol: 'C', x: p4[1].x, y: p4[1].y },
            { symbol: 'C', x: p4[2].x, y: p4[2].y },
            { symbol: 'C', x: p4[3].x, y: p4[3].y }
        ];
        var bonds = ringBonds(4, 0, false);
        var ring = fusedPentagon(atoms[0], atoms[3]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'S', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 0, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 3, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    Templates.cepham = function() {
        var p4 = polygon(4, 0, 0);
        var atoms = [
            { symbol: 'N', x: p4[0].x, y: p4[0].y },
            { symbol: 'C', x: p4[1].x, y: p4[1].y },
            { symbol: 'C', x: p4[2].x, y: p4[2].y },
            { symbol: 'C', x: p4[3].x, y: p4[3].y }
        ];
        var bonds = ringBonds(4, 0, false);
        var ring = fusedHexagon(atoms[0], atoms[3]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'S', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: 0, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: DOUBLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 3, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    Templates.carbapenem = function() {
        var p4 = polygon(4, 0, 0);
        var atoms = [
            { symbol: 'N', x: p4[0].x, y: p4[0].y },
            { symbol: 'C', x: p4[1].x, y: p4[1].y },
            { symbol: 'C', x: p4[2].x, y: p4[2].y },
            { symbol: 'C', x: p4[3].x, y: p4[3].y }
        ];
        var bonds = ringBonds(4, 0, false);
        var ring = fusedPentagon(atoms[0], atoms[3]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 0, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: DOUBLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 3, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 6: tetracycline) =====

    Templates.tetracycline = function() {
        var rA = polygon(6, 0, 0);
        var atoms = rA.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, false);
        var ringB = fusedHexagon(atoms[1], atoms[2]);
        var bs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringB[i].x, y: ringB[i].y });
        }
        bonds.push({ a1: 1, a2: bs, type: SINGLE });
        bonds.push({ a1: bs, a2: bs + 1, type: SINGLE });
        bonds.push({ a1: bs + 1, a2: bs + 2, type: SINGLE });
        bonds.push({ a1: bs + 2, a2: bs + 3, type: SINGLE });
        bonds.push({ a1: bs + 3, a2: 2, type: SINGLE });
        var ringC = fusedHexagon(atoms[bs + 1], atoms[bs + 2]);
        var cs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringC[i].x, y: ringC[i].y });
        }
        bonds.push({ a1: bs + 1, a2: cs, type: SINGLE });
        bonds.push({ a1: cs, a2: cs + 1, type: SINGLE });
        bonds.push({ a1: cs + 1, a2: cs + 2, type: SINGLE });
        bonds.push({ a1: cs + 2, a2: cs + 3, type: SINGLE });
        bonds.push({ a1: cs + 3, a2: bs + 2, type: SINGLE });
        var ringD = fusedHexagon(atoms[cs + 1], atoms[cs + 2]);
        var ds = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringD[i].x, y: ringD[i].y });
        }
        bonds.push({ a1: cs + 1, a2: ds, type: DOUBLE });
        bonds.push({ a1: ds, a2: ds + 1, type: SINGLE });
        bonds.push({ a1: ds + 1, a2: ds + 2, type: DOUBLE });
        bonds.push({ a1: ds + 2, a2: ds + 3, type: SINGLE });
        bonds.push({ a1: ds + 3, a2: cs + 2, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [3] };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 7: alkaloid scaffolds) =====

    Templates.quinuclidine = function() {
        var hex = polygon(6, 0, 0);
        var atoms = hex.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        atoms[0] = { symbol: 'N', x: atoms[0].x, y: atoms[0].y };
        var bonds = ringBonds(6, 0, false);
        var c0 = atoms[0], c3 = atoms[3];
        var midX = (c0.x + c3.x) / 2;
        var midY = (c0.y + c3.y) / 2;
        atoms.push({ symbol: 'C', x: midX + BL * 0.5, y: midY });
        atoms.push({ symbol: 'C', x: midX - BL * 0.5, y: midY });
        bonds.push({ a1: 0, a2: 6, type: SINGLE });
        bonds.push({ a1: 6, a2: 7, type: SINGLE });
        bonds.push({ a1: 7, a2: 3, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    Templates.pyrrolizidine = function() {
        var p5 = polygon(5, 0, 0);
        var atoms = p5.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(5, 0, false);
        var ring = fusedPentagon(atoms[0], atoms[1]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 0, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 1, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    Templates.indolizidine = function() {
        var p6 = polygon(6, 0, 0);
        var atoms = p6.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);
        var ring = fusedPentagon(atoms[0], atoms[1]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        bonds.push({ a1: 0, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: 1, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    Templates.quinolizidine = function() {
        var p6 = polygon(6, 0, 0);
        var atoms = p6.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(6, 0, false);
        var ring = fusedHexagon(atoms[0], atoms[1]);
        var ns = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ring[i].x, y: ring[i].y });
        }
        bonds.push({ a1: 0, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 1, type: SINGLE });
        return { atoms: atoms, bonds: bonds };
    };

    // Aporphine — phenanthrene + saturated 6-ring with N (apomorphine class).
    Templates.aporphine = function() {
        var ph = Templates.phenanthrene();
        var atoms = ph.atoms.map(function(a) { return { symbol: a.symbol, x: a.x, y: a.y }; });
        var bonds = ph.bonds.map(function(b) { return { a1: b.a1, a2: b.a2, type: b.type }; });
        var ring = fusedHexagon(atoms[11], atoms[12]);
        var ns = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'N', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'C', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: 11, a2: ns, type: SINGLE });
        bonds.push({ a1: ns, a2: ns + 1, type: SINGLE });
        bonds.push({ a1: ns + 1, a2: ns + 2, type: SINGLE });
        bonds.push({ a1: ns + 2, a2: ns + 3, type: SINGLE });
        bonds.push({ a1: ns + 3, a2: 12, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0, 1, 2] };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 8: dibenzazepine) =====

    Templates.dibenzazepine = function() {
        var p7 = polygon(7, 0, 0);
        var atoms = p7.map(function(p, i) {
            return { symbol: i === 0 ? 'N' : 'C', x: p.x, y: p.y };
        });
        var bonds = ringBonds(7, 0, false);
        var ringL = fusedHexagon(atoms[4], atoms[5]);
        var ls = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringL[i].x, y: ringL[i].y });
        }
        bonds.push({ a1: 4, a2: ls, type: DOUBLE });
        bonds.push({ a1: ls, a2: ls + 1, type: SINGLE });
        bonds.push({ a1: ls + 1, a2: ls + 2, type: DOUBLE });
        bonds.push({ a1: ls + 2, a2: ls + 3, type: SINGLE });
        bonds.push({ a1: ls + 3, a2: 5, type: DOUBLE });
        var ringR = fusedHexagon(atoms[1], atoms[2]);
        var rs = atoms.length;
        for (var i = 0; i < 4; i++) {
            atoms.push({ symbol: 'C', x: ringR[i].x, y: ringR[i].y });
        }
        bonds.push({ a1: 1, a2: rs, type: DOUBLE });
        bonds.push({ a1: rs, a2: rs + 1, type: SINGLE });
        bonds.push({ a1: rs + 1, a2: rs + 2, type: DOUBLE });
        bonds.push({ a1: rs + 2, a2: rs + 3, type: SINGLE });
        bonds.push({ a1: rs + 3, a2: 2, type: DOUBLE });
        return { atoms: atoms, bonds: bonds, aromatic: [1, 2] };
    };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 9: flavonoids backstop) =====

    // Chroman / flavanone / flavone share a chromene-class topology
    // (benzene fused to a saturated 6-ring with one O). Agent C is
    // concurrently adding chromone / coumarin to the parent branch under
    // different names — these chroman/flavanone/flavone aliases are kept
    // as a backstop so this worktree builds standalone.
    Templates.chroman = function() {
        var r6 = polygon(6, 0, 0);
        var atoms = r6.map(function(p) { return { symbol: 'C', x: p.x, y: p.y }; });
        var bonds = ringBonds(6, 0, true);
        var ring = fusedHexagon(atoms[1], atoms[2]);
        var nb = atoms.length;
        atoms.push({ symbol: 'C', x: ring[0].x, y: ring[0].y });
        atoms.push({ symbol: 'C', x: ring[1].x, y: ring[1].y });
        atoms.push({ symbol: 'C', x: ring[2].x, y: ring[2].y });
        atoms.push({ symbol: 'O', x: ring[3].x, y: ring[3].y });
        bonds.push({ a1: 2, a2: nb, type: SINGLE });
        bonds.push({ a1: nb, a2: nb + 1, type: SINGLE });
        bonds.push({ a1: nb + 1, a2: nb + 2, type: SINGLE });
        bonds.push({ a1: nb + 2, a2: nb + 3, type: SINGLE });
        bonds.push({ a1: nb + 3, a2: 1, type: SINGLE });
        return { atoms: atoms, bonds: bonds, aromatic: [0] };
    };

    Templates.flavanone = function() { return Templates.chroman(); };
    Templates.flavone = function() { return Templates.chroman(); };

    // ===== v1.5.2 SCAFFOLD-COVERAGE EXPANSION (Part 10: open-chain sugars) =====

    Templates.linear_hexose = function() {
        var atoms = [
            { symbol: 'C', x: 0, y: 0 },
            { symbol: 'C', x: 0, y: BL },
            { symbol: 'C', x: 0, y: BL * 2 },
            { symbol: 'C', x: 0, y: BL * 3 },
            { symbol: 'C', x: 0, y: BL * 4 },
            { symbol: 'C', x: 0, y: BL * 5 }
        ];
        var bonds = [
            { a1: 0, a2: 1, type: SINGLE },
            { a1: 1, a2: 2, type: SINGLE },
            { a1: 2, a2: 3, type: SINGLE },
            { a1: 3, a2: 4, type: SINGLE },
            { a1: 4, a2: 5, type: SINGLE }
        ];
        return { atoms: atoms, bonds: bonds };
    };

    Templates.linear_pentose = function() {
        var atoms = [
            { symbol: 'C', x: 0, y: 0 },
            { symbol: 'C', x: 0, y: BL },
            { symbol: 'C', x: 0, y: BL * 2 },
            { symbol: 'C', x: 0, y: BL * 3 },
            { symbol: 'C', x: 0, y: BL * 4 }
        ];
        var bonds = [
            { a1: 0, a2: 1, type: SINGLE },
            { a1: 1, a2: 2, type: SINGLE },
            { a1: 2, a2: 3, type: SINGLE },
            { a1: 3, a2: 4, type: SINGLE }
        ];
        return { atoms: atoms, bonds: bonds };
    };


    // ============== END OF v1.5.2 SCAFFOLD-COVERAGE EXPANSION ===============

    // =========================================================================
    // Apply a template to a Molecule object at a given position
    // =========================================================================

    /**
     * Stamp a template into the molecule at position (cx, cy).
     * Returns an array of the newly created atom IDs.
     */
    Templates.apply = function(mol, templateName, cx, cy) {
        var tmpl;
        if (typeof templateName === 'string') {
            tmpl = Templates[templateName];
            if (typeof tmpl === 'function') tmpl = tmpl();
            if (!tmpl) return [];
        } else {
            tmpl = templateName; // direct template object
        }

        cx = cx || 0;
        cy = cy || 0;

        // Center the template
        var sumX = 0, sumY = 0;
        for (var i = 0; i < tmpl.atoms.length; i++) {
            sumX += tmpl.atoms[i].x;
            sumY += tmpl.atoms[i].y;
        }
        var avgX = sumX / tmpl.atoms.length;
        var avgY = sumY / tmpl.atoms.length;

        var atomIds = [];
        for (var i = 0; i < tmpl.atoms.length; i++) {
            var ta = tmpl.atoms[i];
            var atom = mol.addAtom(ta.symbol, cx + ta.x - avgX, cy + ta.y - avgY);
            atomIds.push(atom.id);
        }

        for (var i = 0; i < tmpl.bonds.length; i++) {
            var tb = tmpl.bonds[i];
            mol.addBond(atomIds[tb.a1], atomIds[tb.a2], tb.type);
        }

        return atomIds;
    };

    /**
     * List all available template names.
     */
    Templates.list = function() {
        var names = [];
        for (var key in Templates) {
            if (typeof Templates[key] === 'function' && key !== 'apply' && key !== 'list') {
                names.push(key);
            }
        }
        return names;
    };

    // =========================================================================
    // Internal helpers for fused ring generation
    // =========================================================================

    /**
     * Generate 4 new vertices for a hexagon fused to edge shared1-shared2.
     * `existingAtoms` (optional) is a list of {x,y} atoms already placed; the
     * function picks the side of the edge farthest from them so the fused ring
     * doesn't overlap. Without `existingAtoms` the right-hand normal is used.
     */
    function fusedHexagon(shared1, shared2, existingAtoms) {
        return fusedRingVertices(shared1, shared2, 6, existingAtoms);
    }

    function fusedPentagon(shared1, shared2, existingAtoms) {
        return fusedRingVertices(shared1, shared2, 5, existingAtoms);
    }

    function fusedRingVertices(shared1, shared2, ringSize, existingAtoms) {
        var mx = (shared1.x + shared2.x) / 2;
        var my = (shared1.y + shared2.y) / 2;
        var dx = shared2.x - shared1.x;
        var dy = shared2.y - shared1.y;
        var nx = -dy, ny = dx;
        var len = Math.sqrt(nx * nx + ny * ny);
        if (len > 0) { nx /= len; ny /= len; }

        var radius = BL / (2 * Math.sin(Math.PI / ringSize));
        var apothem = radius * Math.cos(Math.PI / ringSize);

        // Two candidate centres (either side of the shared edge)
        var cxA = mx + nx * apothem, cyA = my + ny * apothem;
        var cxB = mx - nx * apothem, cyB = my - ny * apothem;

        // Pick the centre farthest from existing atoms (other than the
        // shared edge endpoints which are by definition close).
        var cx = cxA, cy = cyA;
        if (existingAtoms && existingAtoms.length > 0) {
            var minA = Infinity, minB = Infinity;
            for (var i = 0; i < existingAtoms.length; i++) {
                var a = existingAtoms[i];
                if (a === shared1 || a === shared2) continue;
                if (Math.abs(a.x - shared1.x) < 1e-3 && Math.abs(a.y - shared1.y) < 1e-3) continue;
                if (Math.abs(a.x - shared2.x) < 1e-3 && Math.abs(a.y - shared2.y) < 1e-3) continue;
                var dA = (a.x - cxA) * (a.x - cxA) + (a.y - cyA) * (a.y - cyA);
                var dB = (a.x - cxB) * (a.x - cxB) + (a.y - cyB) * (a.y - cyB);
                if (dA < minA) minA = dA;
                if (dB < minB) minB = dB;
            }
            if (minB > minA) { cx = cxB; cy = cyB; }
        }

        var pts = polygon(ringSize, cx, cy);
        var ordered = reorderToMatch(pts, shared1, shared2);
        return ordered.slice(2);
    }

    /**
     * Reorder polygon points so the closest to target1 is first,
     * closest to target2 is second, and remaining follow in order.
     */
    function reorderToMatch(pts, target1, target2) {
        // Find vertex closest to target1
        var best1 = 0, bestDist1 = Infinity;
        for (var i = 0; i < pts.length; i++) {
            var d = distPt(pts[i], target1);
            if (d < bestDist1) { bestDist1 = d; best1 = i; }
        }

        // Find vertex closest to target2 (not the same as best1)
        var best2 = 0, bestDist2 = Infinity;
        for (var i = 0; i < pts.length; i++) {
            if (i === best1) continue;
            var d = distPt(pts[i], target2);
            if (d < bestDist2) { bestDist2 = d; best2 = i; }
        }

        // Determine traversal direction
        var n = pts.length;
        var fwd = (best1 + 1) % n === best2;
        var result = [];
        if (fwd) {
            for (var i = 0; i < n; i++) {
                result.push(pts[(best1 + i) % n]);
            }
        } else {
            for (var i = 0; i < n; i++) {
                result.push(pts[(best1 - i + n) % n]);
            }
        }
        return result;
    }

    function distPt(p1, p2) {
        var dx = p1.x - p2.x, dy = p1.y - p2.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

    // =========================================================================
    // Export
    // =========================================================================

    global.Templates = Templates;

})(typeof window !== 'undefined' ? window : this);
