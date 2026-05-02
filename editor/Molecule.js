/**
 * Molecule.js — Molecular graph data model
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Represents a chemical structure as a graph of atoms and bonds with 2D coordinates.
 * Supports reactions (multiple components separated by a reaction arrow).
 */
(function(global) {
    'use strict';

    var _nextAtomId = 1;
    var _nextBondId = 1;

    // Element properties: valence, colour, mass
    var ELEMENTS = {
        'C':  { valence: 4, color: '#333333', mass: 12 },
        'N':  { valence: 3, color: '#2563eb', mass: 14 },
        'O':  { valence: 2, color: '#dc2626', mass: 16 },
        'S':  { valence: 2, color: '#ca8a04', mass: 32 },
        'P':  { valence: 3, color: '#ea580c', mass: 31 },
        'F':  { valence: 1, color: '#16a34a', mass: 19 },
        'Cl': { valence: 1, color: '#16a34a', mass: 35 },
        'Br': { valence: 1, color: '#9333ea', mass: 80 },
        'I':  { valence: 1, color: '#7c3aed', mass: 127 },
        'H':  { valence: 1, color: '#333333', mass: 1 },
        'Si': { valence: 4, color: '#64748b', mass: 28 },
        'B':  { valence: 3, color: '#f59e0b', mass: 11 },
        'Se': { valence: 2, color: '#ca8a04', mass: 79 },
        'As': { valence: 3, color: '#bd93f9', mass: 75 },
        'Te': { valence: 2, color: '#ca8a04', mass: 128 },
        'Na': { valence: 1, color: '#ab5cf2', mass: 23 },
        'K':  { valence: 1, color: '#8f40d4', mass: 39 },
        'Ca': { valence: 2, color: '#3dff00', mass: 40 },
        'Mg': { valence: 2, color: '#8aff00', mass: 24 },
        'Fe': { valence: 3, color: '#e06633', mass: 56 },
        'Zn': { valence: 2, color: '#7d80b0', mass: 65 },
        'Cu': { valence: 2, color: '#c88033', mass: 64 },
        'Pt': { valence: 4, color: '#d0d0e0', mass: 195 },
        'Li': { valence: 1, color: '#cc80ff', mass: 7 },
        'Al': { valence: 3, color: '#bfa6a6', mass: 27 },
        'Sn': { valence: 4, color: '#668080', mass: 119 },
        'R':  { valence: 1, color: '#64748b', mass: 0 }
    };

    var BOND_SINGLE = 1;
    var BOND_DOUBLE = 2;
    var BOND_TRIPLE = 3;

    var STEREO_NONE = 0;
    var STEREO_WEDGE = 1;
    var STEREO_DASH = 6;

    function Atom(symbol, x, y) {
        this.id = _nextAtomId++;
        this.symbol = symbol || 'C';
        // FIX: use explicit undefined check instead of || to avoid discarding valid 0 coordinates
        this.x = (x !== undefined && x !== null) ? x : 0;
        this.y = (y !== undefined && y !== null) ? y : 0;
        this.charge = 0;
        this.isotope = 0;
        this.mapNumber = 0;          // atom-atom mapping for reactions
        this.hydrogens = -1;         // -1 = auto-calculate
        this.aromatic = false;
        this.chirality = '';
        this.radical = 0;
        this.cipLabel = '';          // CIP label: 'R' or 'S' (assigned by CipStereo)
        this.selected = false;
        this.highlighted = false;
        this.bgColor = null;         // background highlight colour
        this.data = null;            // user-attached data
    }

    Atom.prototype.getElement = function() {
        return ELEMENTS[this.symbol] || ELEMENTS['C'];
    };

    Atom.prototype.clone = function() {
        var a = new Atom(this.symbol, this.x, this.y);
        a.charge = this.charge;
        a.isotope = this.isotope;
        a.mapNumber = this.mapNumber;
        a.hydrogens = this.hydrogens;
        a.aromatic = this.aromatic;
        a.chirality = this.chirality;
        a.radical = this.radical;
        a.cipLabel = this.cipLabel;
        a.bgColor = this.bgColor;
        return a;
    };

    function Bond(atom1Id, atom2Id, type) {
        this.id = _nextBondId++;
        this.atom1 = atom1Id;
        this.atom2 = atom2Id;
        this.type = type || BOND_SINGLE;  // 1, 2, 3
        this.stereo = STEREO_NONE;        // 0=none, 1=wedge, 6=dash
        this.cipLabel = '';          // CIP label: 'E' or 'Z' (assigned by CipStereo)
        this.selected = false;
        this.highlighted = false;
        this.bgColor = null;
        this.data = null;
    }

    Bond.prototype.clone = function() {
        var b = new Bond(this.atom1, this.atom2, this.type);
        b.stereo = this.stereo;
        b.cipLabel = this.cipLabel;
        b.bgColor = this.bgColor;
        return b;
    };

    Bond.prototype.otherAtom = function(atomId) {
        return this.atom1 === atomId ? this.atom2 : this.atom1;
    };

    // =========================================================================
    // Molecule
    // =========================================================================
    function Molecule() {
        this.atoms = [];
        this.bonds = [];
        this.name = '';               // molecule name (line 1 of MOL header / SDF > <NAME>)
        this.program = '';             // line 2 of MOL header (program/timestamp stamp)
        this.comment = '';             // line 3 of MOL header (free-text comment) / SDF > <COMMENT>
        this.reactionArrow = null;    // { x1, y1, x2, y2 }
        this.reactionPlusSigns = null; // FIX: initialize property read by getBounds() to avoid undefined access
        this._atomMap = {};           // id → atom
        this._bondMap = {};           // id → bond
        this._adjacency = {};         // atomId → [bondId, ...]
    }

    // --- Atom operations ---

    Molecule.prototype.addAtom = function(symbol, x, y) {
        var atom = new Atom(symbol, x, y);
        this.atoms.push(atom);
        this._atomMap[atom.id] = atom;
        this._adjacency[atom.id] = [];
        return atom;
    };

    Molecule.prototype.removeAtom = function(atomId) {
        // Remove all bonds connected to this atom
        var bondsToRemove = (this._adjacency[atomId] || []).slice();
        for (var i = 0; i < bondsToRemove.length; i++) {
            this.removeBond(bondsToRemove[i]);
        }
        // Remove atom
        this.atoms = this.atoms.filter(function(a) { return a.id !== atomId; });
        delete this._atomMap[atomId];
        delete this._adjacency[atomId];
    };

    Molecule.prototype.getAtom = function(atomId) {
        return this._atomMap[atomId] || null;
    };

    Molecule.prototype.getAtomAt = function(x, y, radius) {
        radius = radius || 15;
        var rSq = radius * radius;
        for (var i = this.atoms.length - 1; i >= 0; i--) {
            var a = this.atoms[i];
            var dx = a.x - x, dy = a.y - y;
            if (dx * dx + dy * dy <= rSq) return a;
        }
        return null;
    };

    // --- Bond operations ---

    Molecule.prototype.addBond = function(atom1Id, atom2Id, type) {
        // Check if bond already exists
        var existing = this.getBondBetween(atom1Id, atom2Id);
        if (existing) return existing;

        var bond = new Bond(atom1Id, atom2Id, type);
        this.bonds.push(bond);
        this._bondMap[bond.id] = bond;
        if (!this._adjacency[atom1Id]) this._adjacency[atom1Id] = [];
        if (!this._adjacency[atom2Id]) this._adjacency[atom2Id] = [];
        this._adjacency[atom1Id].push(bond.id);
        this._adjacency[atom2Id].push(bond.id);
        return bond;
    };

    Molecule.prototype.removeBond = function(bondId) {
        var bond = this._bondMap[bondId];
        if (!bond) return;
        this.bonds = this.bonds.filter(function(b) { return b.id !== bondId; });
        delete this._bondMap[bondId];
        // Remove from adjacency
        var adj1 = this._adjacency[bond.atom1];
        if (adj1) this._adjacency[bond.atom1] = adj1.filter(function(id) { return id !== bondId; });
        var adj2 = this._adjacency[bond.atom2];
        if (adj2) this._adjacency[bond.atom2] = adj2.filter(function(id) { return id !== bondId; });
    };

    Molecule.prototype.getBond = function(bondId) {
        return this._bondMap[bondId] || null;
    };

    Molecule.prototype.getBondBetween = function(atom1Id, atom2Id) {
        var adj = this._adjacency[atom1Id] || [];
        for (var i = 0; i < adj.length; i++) {
            var bond = this._bondMap[adj[i]];
            if (bond && (bond.atom1 === atom2Id || bond.atom2 === atom2Id)) return bond;
        }
        return null;
    };

    Molecule.prototype.getBondAt = function(x, y, radius) {
        radius = radius || 8;
        var best = null, bestDist = radius;
        for (var i = 0; i < this.bonds.length; i++) {
            var b = this.bonds[i];
            var a1 = this._atomMap[b.atom1];
            var a2 = this._atomMap[b.atom2];
            if (!a1 || !a2) continue;
            var dist = pointToSegmentDist(x, y, a1.x, a1.y, a2.x, a2.y);
            if (dist < bestDist) {
                bestDist = dist;
                best = b;
            }
        }
        return best;
    };

    Molecule.prototype.getBondsOfAtom = function(atomId) {
        var bondIds = this._adjacency[atomId] || [];
        var self = this;
        return bondIds.map(function(id) { return self._bondMap[id]; }).filter(Boolean);
    };

    Molecule.prototype.getNeighbors = function(atomId) {
        var bonds = this.getBondsOfAtom(atomId);
        return bonds.map(function(b) { return b.otherAtom(atomId); });
    };

    Molecule.prototype.degree = function(atomId) {
        return (this._adjacency[atomId] || []).length;
    };

    Molecule.prototype.bondOrderSum = function(atomId) {
        var bonds = this.getBondsOfAtom(atomId);
        var sum = 0;
        for (var i = 0; i < bonds.length; i++) sum += bonds[i].type;
        return sum;
    };

    // --- Hydrogen calculation ---

    // Aromatic valences: when an atom is flagged as aromatic and its bonds
    // are stored as single (from aromatic SMILES input), one valence is
    // reserved for the pi-system contribution.  This table gives the
    // effective valence for implicit-H calculation on aromatic atoms.
    // Multi-valence table: lists allowed valences (ascending) for elements
    // that can be hypervalent.  calcHydrogens picks the smallest valence
    // that accommodates the current bond order sum + charge.
    var MULTI_VALENCE = {
        'N': [3, 5], 'P': [3, 5], 'As': [3, 5],
        'S': [2, 4, 6], 'Se': [2, 4, 6], 'Te': [2, 4, 6],
        'B': [3], 'C': [4], 'O': [2], 'F': [1], 'Cl': [1, 3, 5, 7],
        'Br': [1, 3, 5, 7], 'I': [1, 3, 5, 7]
    };

    var AROMATIC_VALENCE = {
        'C': 3,   // 4 - 1 pi electron => 3 sigma bonds
        'N': 2,   // pyridine-type N: donates 1 pi electron, 2 sigma bonds
                  // (pyrrole-type [nH] uses explicit H path, not this table)
        'O': 2,   // furan-type O: contributes lone pair, 2 ring bonds
        'S': 2,   // thiophene-type S: same as O
        'B': 2,   // boron aromatic: 3 - 1 pi = 2 sigma bonds
        'P': 2    // phosphorus aromatic: same pattern as N
    };

    Molecule.prototype.calcHydrogens = function(atomId) {
        var atom = this._atomMap[atomId];
        if (!atom) return 0;
        if (atom.hydrogens >= 0) return atom.hydrogens;
        var elem = ELEMENTS[atom.symbol];
        if (!elem) return 0;
        var bondSum = this.bondOrderSum(atomId);

        // For aromatic atoms whose bonds are stored as single (from aromatic
        // SMILES input like c1ccccc1), use the aromatic valence which accounts
        // for the pi electron donated to the aromatic system.
        if (atom.aromatic && AROMATIC_VALENCE[atom.symbol] !== undefined) {
            // Only apply aromatic valence when the atom has no explicit
            // double/triple bonds (i.e. aromatic SMILES input, not Kekule).
            var bonds = this.getBondsOfAtom(atomId);
            var hasMultipleBond = false;
            for (var i = 0; i < bonds.length; i++) {
                if (bonds[i].type > BOND_SINGLE) { hasMultipleBond = true; break; }
            }
            if (!hasMultipleBond) {
                var arVal = AROMATIC_VALENCE[atom.symbol];
                return Math.max(0, arVal - bondSum - Math.abs(atom.charge));
            }
        }

        // Use multi-valence table: pick the smallest standard valence that
        // accommodates the current bond order sum + absolute charge.
        var valences = MULTI_VALENCE[atom.symbol] || [elem.valence];
        var hCount = 0;
        var adjBondSum = bondSum + Math.abs(atom.charge);
        for (var vi = 0; vi < valences.length; vi++) {
            if (valences[vi] >= adjBondSum) {
                hCount = valences[vi] - adjBondSum;
                break;
            }
        }
        return hCount;
    };

    // --- Molecule operations ---

    Molecule.prototype.clear = function() {
        this.atoms = [];
        this.bonds = [];
        this.name = '';
        this.program = '';      // FIX: also clear MOL header program line on reset
        this.comment = '';      // FIX: also clear MOL header comment line on reset
        this._atomMap = {};
        this._bondMap = {};
        this._adjacency = {};
        this.reactionArrow = null;
        this.reactionPlusSigns = null; // FIX: clear plus signs to avoid stale data after reset
    };

    Molecule.prototype.isEmpty = function() {
        return this.atoms.length === 0;
    };

    Molecule.prototype.clone = function() {
        var mol = new Molecule();
        mol.name = this.name;
        mol.program = this.program; // FIX: clone preserves MOL header program/comment
        mol.comment = this.comment;
        var idMap = {};
        for (var i = 0; i < this.atoms.length; i++) {
            var a = this.atoms[i].clone();
            var oldId = this.atoms[i].id;
            mol.atoms.push(a);
            mol._atomMap[a.id] = a;
            mol._adjacency[a.id] = [];
            idMap[oldId] = a.id;
        }
        for (var i = 0; i < this.bonds.length; i++) {
            var b = this.bonds[i].clone();
            // FIX: use 'in' check instead of || to handle idMap value of 0 (theoretically possible)
            b.atom1 = (this.bonds[i].atom1 in idMap) ? idMap[this.bonds[i].atom1] : b.atom1;
            b.atom2 = (this.bonds[i].atom2 in idMap) ? idMap[this.bonds[i].atom2] : b.atom2;
            mol.bonds.push(b);
            mol._bondMap[b.id] = b;
            if (mol._adjacency[b.atom1]) mol._adjacency[b.atom1].push(b.id);
            if (mol._adjacency[b.atom2]) mol._adjacency[b.atom2].push(b.id);
        }
        if (this.reactionArrow) {
            mol.reactionArrow = { x1: this.reactionArrow.x1, y1: this.reactionArrow.y1,
                                  x2: this.reactionArrow.x2, y2: this.reactionArrow.y2,
                                  type: this.reactionArrow.type,
                                  conditions: this.reactionArrow.conditions };
        }
        // FIX: also copy plus signs; previously clone() dropped them silently
        if (this.reactionPlusSigns) {
            mol.reactionPlusSigns = this.reactionPlusSigns.map(function(p) {
                return { x: p.x, y: p.y };
            });
        }
        return mol;
    };

    Molecule.prototype.getBounds = function() {
        if (this.atoms.length === 0) return { x: 0, y: 0, w: 0, h: 0 };
        var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        for (var i = 0; i < this.atoms.length; i++) {
            var a = this.atoms[i];
            if (a.x < minX) minX = a.x;
            if (a.y < minY) minY = a.y;
            if (a.x > maxX) maxX = a.x;
            if (a.y > maxY) maxY = a.y;
        }
        // Include reaction arrow endpoints in bounds so zoom-to-fit
        // keeps the arrow visible in the viewport
        if (this.reactionArrow) {
            var arr = this.reactionArrow;
            if (arr.x1 < minX) minX = arr.x1;
            if (arr.x2 < minX) minX = arr.x2;
            if (arr.y1 < minY) minY = arr.y1;
            if (arr.y2 < minY) minY = arr.y2;
            if (arr.x1 > maxX) maxX = arr.x1;
            if (arr.x2 > maxX) maxX = arr.x2;
            if (arr.y1 > maxY) maxY = arr.y1;
            if (arr.y2 > maxY) maxY = arr.y2;
        }
        // Include plus sign positions in bounds
        if (this.reactionPlusSigns) {
            for (var pi = 0; pi < this.reactionPlusSigns.length; pi++) {
                var ps = this.reactionPlusSigns[pi];
                if (ps.x < minX) minX = ps.x;
                if (ps.y < minY) minY = ps.y;
                if (ps.x > maxX) maxX = ps.x;
                if (ps.y > maxY) maxY = ps.y;
            }
        }
        return { x: minX, y: minY, w: maxX - minX, h: maxY - minY };
    };

    // Get connected components (for reaction support)
    Molecule.prototype.getComponents = function() {
        var visited = {};
        var components = [];
        var self = this;

        function dfs(atomId, comp) {
            if (visited[atomId]) return;
            visited[atomId] = true;
            comp.push(atomId);
            var neighbors = self.getNeighbors(atomId);
            for (var i = 0; i < neighbors.length; i++) {
                dfs(neighbors[i], comp);
            }
        }

        for (var i = 0; i < this.atoms.length; i++) {
            if (!visited[this.atoms[i].id]) {
                var comp = [];
                dfs(this.atoms[i].id, comp);
                components.push(comp);
            }
        }
        return components;
    };

    // Find all rings using Smallest Set of Smallest Rings (SSSR)
    Molecule.prototype.findRings = function(maxSize) {
        maxSize = maxSize || 8;
        var rings = [];
        var seenKeys = {};
        var self = this;

        // DFS-based ring detection: for each atom, do a depth-limited DFS
        // looking for paths that return to the start atom.
        var visitCount = 0;
        var MAX_VISITS = 50000;

        for (var i = 0; i < this.atoms.length; i++) {
            var startId = this.atoms[i].id;
            // DFS stack: [currentAtomId, path, visitedSet]
            var stack = [[startId, [startId], {}]];
            stack[0][2][startId] = true;

            while (stack.length > 0) {
                if (++visitCount > MAX_VISITS) break;
                var frame = stack.pop();
                var current = frame[0];
                var path = frame[1];
                var visited = frame[2];

                if (path.length > maxSize + 1) continue;

                var neighbors = self.getNeighbors(current);
                for (var j = 0; j < neighbors.length; j++) {
                    var n = neighbors[j];

                    // Ring closure: reached start atom via 3+ atoms
                    if (n === startId && path.length >= 3) {
                        var ring = path.slice();
                        // Normalise: start from smallest id, consistent winding
                        var minIdx = 0;
                        for (var k = 1; k < ring.length; k++) {
                            if (ring[k] < ring[minIdx]) minIdx = k;
                        }
                        var norm = ring.slice(minIdx).concat(ring.slice(0, minIdx));
                        // Also check reverse direction
                        var rev = [norm[0]].concat(norm.slice(1).reverse());
                        var key1 = norm.join(',');
                        var key2 = rev.join(',');
                        var key = key1 < key2 ? key1 : key2;
                        if (!seenKeys[key]) {
                            seenKeys[key] = true;
                            rings.push({ atoms: norm, key: key, size: ring.length });
                        }
                        continue;
                    }

                    // Don't revisit atoms already in this path (except startId handled above)
                    if (visited[n]) continue;
                    if (path.length >= maxSize) continue;

                    // Extend path
                    var newVisited = {};
                    for (var vk in visited) newVisited[vk] = true;
                    newVisited[n] = true;
                    stack.push([n, path.concat(n), newVisited]);
                }
            }
        }
        return rings;
    };

    // --- Ring templates ---

    var BOND_LENGTH = 30; // default bond length in pixels

    Molecule.prototype.addRing = function(size, cx, cy, fuseAtomId) {
        var angle = 2 * Math.PI / size;
        var radius = BOND_LENGTH / (2 * Math.sin(Math.PI / size));
        var atoms = [];
        var startAngle = -Math.PI / 2; // start at top

        if (fuseAtomId) {
            // Fuse ring onto existing atom
            var fuseAtom = this._atomMap[fuseAtomId];
            if (!fuseAtom) return;
            cx = fuseAtom.x;
            cy = fuseAtom.y - radius;
            atoms.push(fuseAtomId);
            startAngle = Math.PI / 2 - angle; // start from fused atom position
        }

        for (var i = (fuseAtomId ? 1 : 0); i < size; i++) {
            var a = startAngle + i * angle;
            var x = cx + radius * Math.cos(a);
            var y = cy + radius * Math.sin(a);
            var atom = this.addAtom('C', x, y);
            atoms.push(atom.id);
        }

        // Add bonds in ring
        for (var i = 0; i < atoms.length; i++) {
            var next = (i + 1) % atoms.length;
            this.addBond(atoms[i], atoms[next], BOND_SINGLE);
        }

        return atoms;
    };

    // --- Serialisation helpers ---

    Molecule.prototype.toJSON = function() {
        return {
            name: this.name,
            // FIX: include MOL header program/comment so undo/redo doesn't
            // strip header metadata after a MOL load.
            program: this.program || '',
            comment: this.comment || '',
            atoms: this.atoms.map(function(a) {
                // FIX: persist `aromatic` so undo/redo doesn't drop SMILES-derived
                // aromatic flags (which the renderer uses for the inner-circle style).
                // FIX: also persist `chirality` (otherwise stereo `@` / `@@` is silently
                // dropped on undo/redo), `mapHighlighted` + `bgColor` (so AAM
                // cross-highlight pins survive undo/redo), and `cipLabel` (so
                // R/S labels stick after undo).
                // FIX: also persist `radical` so undo/redo preserves radical
                // multiplicity. Radicals are read from MOL `M  RAD` / V3000 `RAD=`
                // and would otherwise be silently dropped on the first history op.
                return { id: a.id, symbol: a.symbol, x: a.x, y: a.y, charge: a.charge,
                         isotope: a.isotope, mapNumber: a.mapNumber, hydrogens: a.hydrogens,
                         aromatic: !!a.aromatic,
                         chirality: a.chirality || '',
                         cipLabel: a.cipLabel || '',
                         radical: a.radical || 0,
                         mapHighlighted: !!a.mapHighlighted,
                         bgColor: a.bgColor || null };
            }),
            bonds: this.bonds.map(function(b) {
                // FIX: persist `bgColor` (SMARTS / map highlight halo) and
                // `cipLabel` (E/Z) so undo/redo doesn't strip them.
                return { id: b.id, atom1: b.atom1, atom2: b.atom2, type: b.type,
                         stereo: b.stereo,
                         cipLabel: b.cipLabel || '',
                         bgColor: b.bgColor || null };
            }),
            reactionArrow: this.reactionArrow,
            // FIX: include plus signs so undo/redo round-trips reactions correctly
            reactionPlusSigns: this.reactionPlusSigns
                ? this.reactionPlusSigns.map(function(p) { return { x: p.x, y: p.y }; })
                : null
        };
    };

    // --- Utility ---

    function pointToSegmentDist(px, py, x1, y1, x2, y2) {
        var dx = x2 - x1, dy = y2 - y1;
        var lenSq = dx * dx + dy * dy;
        if (lenSq === 0) return Math.sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1));
        var t = Math.max(0, Math.min(1, ((px - x1) * dx + (py - y1) * dy) / lenSq));
        var projX = x1 + t * dx, projY = y1 + t * dy;
        return Math.sqrt((px - projX) * (px - projX) + (py - projY) * (py - projY));
    }

    // --- Exports ---

    Molecule.Atom = Atom;
    Molecule.Bond = Bond;
    Molecule.ELEMENTS = ELEMENTS;
    Molecule.BOND_SINGLE = BOND_SINGLE;
    Molecule.BOND_DOUBLE = BOND_DOUBLE;
    Molecule.BOND_TRIPLE = BOND_TRIPLE;
    Molecule.STEREO_NONE = STEREO_NONE;
    Molecule.STEREO_WEDGE = STEREO_WEDGE;
    Molecule.STEREO_DASH = STEREO_DASH;
    Molecule.BOND_LENGTH = BOND_LENGTH;

    global.Molecule = Molecule;

})(typeof window !== 'undefined' ? window : this);
