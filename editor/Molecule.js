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
        this.cipLabel = '';          // CIP label: 'R' or 'S' (assigned by CIPStereo)
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
        // FIX: previously clone() dropped mapHighlighted and data, so
        // History.toJSON()/fromJSON drifted from the live atom on undo/redo.
        a.mapHighlighted = !!this.mapHighlighted;
        a.data = this.data;
        return a;
    };

    function Bond(atom1Id, atom2Id, type) {
        this.id = _nextBondId++;
        this.atom1 = atom1Id;
        this.atom2 = atom2Id;
        this.type = type || BOND_SINGLE;  // 1, 2, 3
        this.stereo = STEREO_NONE;        // 0=none, 1=wedge, 6=dash
        this.cipLabel = '';          // CIP label: 'E' or 'Z' (assigned by CIPStereo)
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
        // FIX: also copy `data` so user-attached metadata survives clone().
        b.data = this.data;
        return b;
    };

    Bond.prototype.otherAtom = function(atomId) {
        return this.atom1 === atomId ? this.atom2 : this.atom1;
    };

    // v2.0.31: curly-arrow endpoint remap helper. Used by Molecule.clone()
    // to chase atom-endpoint refs through the per-clone idMap and
    // bond-endpoint refs through the parallel bondIdMap.
    // v2.0.33: delegate to editor/CurlyArrowEndpoint.js when available so the
    // kind-dispatch lives in one place. Fallback closure keeps Molecule
    // self-sufficient if loaded outside the bundle.
    function remapEndpoint(ep, idMap, bondIdMap) {
        if (typeof CurlyArrowEndpoint !== 'undefined' && CurlyArrowEndpoint.remapEndpoint) {
            return CurlyArrowEndpoint.remapEndpoint(ep, idMap, bondIdMap);
        }
        if (!ep) { return ep; }
        if (ep.kind === 'bond') {
            return { kind: 'bond', bondId: (ep.bondId in bondIdMap) ? bondIdMap[ep.bondId] : ep.bondId };
        }
        return { kind: 'atom', atomId: (ep.atomId in idMap) ? idMap[ep.atomId] : ep.atomId };
    }

    // =========================================================================
    // Molecule
    // =========================================================================
    function Molecule() {
        this.atoms = [];
        this.bonds = [];
        this.name = '';               // molecule name (line 1 of MOL header / SDF > <NAME>)
        this.program = '';             // line 2 of MOL header (program/timestamp stamp)
        this.comment = '';             // line 3 of MOL header (free-text comment) / SDF > <COMMENT>
        this.data = null;              // optional molecule-level metadata
        this.reactionArrow = null;    // { x1, y1, x2, y2 }
        this.reactionPlusSigns = null; // FIX: initialize property read by getBounds() to avoid undefined access
        // v2.0.29: Mechanism curly-arrow annotations — curly electron-flow
        // arrows. Each entry: { id, style: 'pair'|'single',
        //                       from: { kind, atomId }, to: { kind, atomId },
        //                       controlOffset: <number|null> }
        // See editor/CurlyArrows.js for the geometry + SVG factory.
        this.curlyArrows = [];
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
        // v2.0.29: also drop any curly arrows that reference this atom so
        // we never end up rendering an arrow to a deleted endpoint.
        // v2.0.31: cascade covers BOTH endpoint kinds — atom endpoints by
        // direct id match, bond endpoints by checking if the to-be-removed
        // bonds list contains the arrow's bondId.
        if (this.curlyArrows && this.curlyArrows.length) {
            var doomedBondIds = {};
            for (var bi = 0; bi < bondsToRemove.length; bi++) { doomedBondIds[bondsToRemove[bi]] = true; }
            this.curlyArrows = this.curlyArrows.filter(function(ca) {
                if (ca.from.kind === 'atom' && ca.from.atomId === atomId) { return false; }
                if (ca.to.kind   === 'atom' && ca.to.atomId   === atomId) { return false; }
                if (ca.from.kind === 'bond' && doomedBondIds[ca.from.bondId]) { return false; }
                if (ca.to.kind   === 'bond' && doomedBondIds[ca.to.bondId])   { return false; }
                return true;
            });
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
        // v2.0.31: cascade into bond-endpoint curly arrows so a removed bond
        // doesn't leave an arrow pointing into thin air.
        if (this.curlyArrows && this.curlyArrows.length) {
            this.curlyArrows = this.curlyArrows.filter(function(ca) {
                if (ca.from.kind === 'bond' && ca.from.bondId === bondId) { return false; }
                if (ca.to.kind   === 'bond' && ca.to.bondId   === bondId) { return false; }
                return true;
            });
        }
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

    // v1.8.26: STANDARD_VALENCES + getDefaultValence were each defined
    // byte-identically (modulo whitespace) in BOTH SmilesParser.js and
    // SmilesWriter.js. They are the implicit-H valence model used by the
    // SMILES round-trip layer; this is the single canonical home.
    // (Distinct from MULTI_VALENCE above — e.g. Cl is [1] here, [1,3,5,7]
    // there — because the SMILES layer uses the conservative neutral
    // valence set, not the hypervalent set.)
    var STANDARD_VALENCES = {
        'B':  [3],        'C':  [4],        'N':  [3, 5],
        'O':  [2],        'P':  [3, 5],     'S':  [2, 4, 6],
        'F':  [1],        'Cl': [1],        'Br': [1],
        'I':  [1],        'H':  [1],        'Si': [4],
        'Se': [2, 4, 6],  'As': [3, 5],     'Te': [2, 4, 6]
    };

    // getDefaultValence(symbol, bondSum, charge) → the smallest standard
    // valence that covers the atom's current bond order + |charge|, or
    // the largest standard valence if none is large enough. Returns 0
    // for elements not in STANDARD_VALENCES.
    function getDefaultValence(symbol, bondSum, charge) {
        var vList = STANDARD_VALENCES[symbol];
        if (!vList) return 0;
        var target = bondSum + Math.abs(charge || 0);
        for (var i = 0; i < vList.length; i++) {
            if (vList[i] >= target) return vList[i];
        }
        return vList[vList.length - 1];
    }

    // v1.8.26: symbol → atomic-number lookup. Previously duplicated four
    // times (CIPStereo.js as ATOMIC_NUMBER, SMSDGraph.js / SmilesWriter.js
    // / SmartsParser.js as ATOMIC_NUMBERS) with three different coverage
    // levels. This is the canonical superset: H(1)–U(92) plus the
    // R-group placeholder 'R': 0. All four modules now alias this.
    var ATOMIC_NUMBERS = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
        'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
        'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
        'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
        'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
        'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
        'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
        'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
        'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
        'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
        'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
        'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
        'R': 0
    };

    // Canonical Hückel aromaticity perception. This is the single shared
    // implementation for writer and matcher paths; it inspects both ring
    // bonds at each atom, so the result is independent of traversal order.
    //
    // Returns a plain object { atomId: true, ... } of every atom that
    // belongs to a perceived-aromatic ring (rings up to size 8, Hückel
    // 4n+2 pi-electron count, with charge- and heteroatom-aware
    // pi-electron contributions).
    function perceiveAromaticity(mol) {
        var aromaticAtoms = {};
        var aromaticElements = { 'C': true, 'N': true, 'O': true, 'S': true };
        var rings = mol.findRings(8);

        for (var r = 0; r < rings.length; r++) {
            var ring = rings[r].atoms;
            if (ring.length < 3) continue;

            // Every ring atom must be a potential aromatic element.
            var allAromatic = true;
            for (var i = 0; i < ring.length; i++) {
                var atom = mol.getAtom(ring[i]);
                if (!atom || !aromaticElements[atom.symbol]) {
                    allAromatic = false;
                    break;
                }
            }
            if (!allAromatic) continue;

            // Count pi electrons, then apply the Hückel 4n+2 rule.
            var piElectrons = 0;
            var valid = true;
            for (var i = 0; i < ring.length; i++) {
                var atomId = ring[i];
                var atom = mol.getAtom(atomId);
                var sym = atom.symbol;

                // Inspect BOTH ring bonds at this atom for a double bond,
                // not just the bond to (i+1) — a Kekulé ring places the
                // double bond on alternating sides, so checking only one
                // direction under-counts. This order-independence is the
                // fix that the divergent copies lacked.
                var prevId = ring[(i - 1 + ring.length) % ring.length];
                var nextId = ring[(i + 1) % ring.length];
                var bondPrev = mol.getBondBetween(atomId, prevId);
                var bondNext = mol.getBondBetween(atomId, nextId);
                var hasDoubleBondInRing =
                    (bondPrev && bondPrev.type === BOND_DOUBLE) ||
                    (bondNext && bondNext.type === BOND_DOUBLE);

                if (sym === 'C') {
                    if (hasDoubleBondInRing) {
                        piElectrons += 1;            // one electron from the ring double bond
                    } else {
                        // Exocyclic double bond → this C contributes no pi
                        // electron to the ring system (e.g. a ring ketone C).
                        var exoBonds = mol.getBondsOfAtom(atomId);
                        var hasExoDouble = false;
                        for (var j = 0; j < exoBonds.length; j++) {
                            var other = exoBonds[j].otherAtom(atomId);
                            if (ring.indexOf(other) < 0 &&
                                exoBonds[j].type === BOND_DOUBLE) {
                                hasExoDouble = true;
                                break;
                            }
                        }
                        if (hasExoDouble) { valid = false; break; }
                        // Charged carbon with no in-ring double bond: the
                        // cyclopentadienyl anion [cH-] donates its lone pair
                        // (2 e-); tropylium [c+] has an empty p-orbital (0 e-).
                        var cCharge = atom.charge || 0;
                        if (cCharge === -1) {
                            piElectrons += 2;
                        } else if (cCharge === 1) {
                            piElectrons += 0;
                        } else if (atom.aromatic) {
                            piElectrons += 1;        // aromatic-flagged 'c'
                        } else {
                            // Neutral sp3 C, no in-ring and no exocyclic
                            // double bond → not part of any pi system.
                            valid = false; break;
                        }
                    }
                } else if (sym === 'N') {
                    if (hasDoubleBondInRing) {
                        piElectrons += 1;            // pyridine-type imine N
                    } else {
                        // Distinguish pyrrole-N (lone pair, 2 e-) from
                        // pyridine-N (1 e-) when ring bonds are all single
                        // (aromatic SMILES input or fully Kekulé-free ring).
                        var explH = atom.hydrogens;  // >=0 for bracket atoms, -1 for organic
                        if (atom.aromatic) {
                            // A neutral aromatic N donates its lone pair to the
                            // pi system (pyrrole-type, 2 e-) when it is sigma-
                            // saturated with 3 bonds (2 ring + H OR a C/N
                            // substituent) and has no in-ring double bond. This
                            // covers pyrrole '[nH]' AND N-alkyl / N-aryl
                            // pyrroles & imidazoles (organic 'n' with a
                            // substituent and no H, e.g. N-methylimidazole,
                            // caffeine, 1-methylindole). A degree-2 aromatic N
                            // keeps its lone pair in-plane (pyridine-type, 1 e-).
                            // A positively-charged N (pyridinium '[nH+]') is
                            // pyridine-type regardless — the lone pair is
                            // protonated / consumed.
                            var nCharge = atom.charge || 0;
                            var nSigma = mol.getNeighbors(atomId).length +
                                (explH > 0 ? explH : 0);
                            piElectrons += (nCharge <= 0 && nSigma >= 3) ? 2 : 1;
                        } else {
                            // Non-aromatic (Kekule) input with no in-ring double
                            // bond. Pyrrole-type (2 e-) when: bracket '[nH]';
                            // OR a neutral N that is sigma-saturated with 3
                            // all-single bonds (e.g. Kekule N-methylpyrrole
                            // 'CN1C=CC=C1'); OR the organic 2-neighbour all-
                            // single heuristic (implicit-H pyrrole). Otherwise
                            // pyridine-type (1 e-).
                            var bondSum = mol.bondOrderSum(atomId);
                            var totalNeighbors = mol.getNeighbors(atomId).length;
                            var nCharge2 = atom.charge || 0;
                            var nSigma2 = totalNeighbors + (explH > 0 ? explH : 0);
                            if (explH > 0 ||
                                (nCharge2 <= 0 && nSigma2 >= 3 && bondSum === totalNeighbors) ||
                                (explH < 0 && totalNeighbors === 2 && bondSum === 2)) {
                                piElectrons += 2;
                            } else {
                                piElectrons += 1;
                            }
                        }
                    }
                } else if (sym === 'O' || sym === 'S') {
                    piElectrons += 2;                // furan / thiophene lone pair
                }
            }
            if (!valid) continue;

            var isHueckel = false;
            for (var n = 0; n <= 5; n++) {
                if (piElectrons === 4 * n + 2) { isHueckel = true; break; }
            }
            if (isHueckel) {
                for (var i = 0; i < ring.length; i++) {
                    aromaticAtoms[ring[i]] = true;
                }
            }
        }

        return aromaticAtoms;
    }

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
        this.data = null;
        this._atomMap = {};
        this._bondMap = {};
        this._adjacency = {};
        this.reactionArrow = null;
        this.reactionPlusSigns = null; // FIX: clear plus signs to avoid stale data after reset
        this.curlyArrows = [];         // v2.0.29: also clear mechanism arrows on reset
    };

    Molecule.prototype.isEmpty = function() {
        return this.atoms.length === 0;
    };

    Molecule.prototype.clone = function() {
        var mol = new Molecule();
        mol.name = this.name;
        mol.program = this.program; // FIX: clone preserves MOL header program/comment
        mol.comment = this.comment;
        mol.data = this.data;
        var idMap = {};
        for (var i = 0; i < this.atoms.length; i++) {
            var a = this.atoms[i].clone();
            var oldId = this.atoms[i].id;
            mol.atoms.push(a);
            mol._atomMap[a.id] = a;
            mol._adjacency[a.id] = [];
            idMap[oldId] = a.id;
        }
        // v2.0.31: keep a bond-id remap so bond-endpoint curly arrows
        // (cloned below) can chase their new bond ids in the new molecule.
        var bondIdMap = {};
        for (var i = 0; i < this.bonds.length; i++) {
            var b = this.bonds[i].clone();
            // FIX: use 'in' check instead of || to handle idMap value of 0 (theoretically possible)
            b.atom1 = (this.bonds[i].atom1 in idMap) ? idMap[this.bonds[i].atom1] : b.atom1;
            b.atom2 = (this.bonds[i].atom2 in idMap) ? idMap[this.bonds[i].atom2] : b.atom2;
            bondIdMap[this.bonds[i].id] = b.id;
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
            if (this.reactionArrow.labelParts) {
                mol.reactionArrow.labelParts = {
                    conditions: this.reactionArrow.labelParts.conditions || '',
                    yield: this.reactionArrow.labelParts.yield || '',
                    note: this.reactionArrow.labelParts.note || ''
                };
            }
        }
        // FIX: also copy plus signs; previously clone() dropped them silently
        if (this.reactionPlusSigns) {
            mol.reactionPlusSigns = this.reactionPlusSigns.map(function(p) {
                return { x: p.x, y: p.y };
            });
        }
        // v2.0.29 + v2.0.31: clone mechanism curly arrows, remapping endpoint
        // refs through the idMap (atom endpoints) or bondIdMap (bond
        // endpoints) so they stay valid in the cloned molecule.
        if (this.curlyArrows && this.curlyArrows.length) {
            mol.curlyArrows = this.curlyArrows.map(function(ca) {
                return {
                    id: ca.id,
                    style: ca.style,
                    from: remapEndpoint(ca.from, idMap, bondIdMap),
                    to:   remapEndpoint(ca.to,   idMap, bondIdMap),
                    controlOffset: ca.controlOffset
                };
            });
        }
        return mol;
    };

    // --- v2.0.29: Curly-arrow (mechanism) operations ---
    // These are wrappers around CurlyArrows.createArrow that enforce
    // referential integrity (the addressed atoms must exist on this
    // molecule) and keep removeAtom from leaving orphaned arrows behind.

    // v2.0.34: soft cap on the per-molecule curly-arrow count. Tied to the
    // v2.0.31 security panel deferral (#3) — a malicious JSON import that
    // tried to push hundreds of thousands of arrows would otherwise consume
    // browser memory linearly. 10k is well above any realistic mechanism
    // diagram (e.g. a full enzyme reaction usually has 4-8 arrows).
    Molecule.MAX_CURLY_ARROWS = 10000;

    Molecule.prototype.addCurlyArrow = function(spec) {
        if (typeof CurlyArrows === 'undefined' || !CurlyArrows.createArrow) {
            throw new Error('CurlyArrows module not loaded; cannot add curly arrow');
        }
        if (this.curlyArrows.length >= Molecule.MAX_CURLY_ARROWS) {
            throw new Error('addCurlyArrow: per-molecule cap of ' + Molecule.MAX_CURLY_ARROWS + ' arrows reached');
        }
        this._assertEndpointResolves(spec.from, 'from');
        this._assertEndpointResolves(spec.to,   'to');
        var arrow = CurlyArrows.createArrow(spec);
        this.curlyArrows.push(arrow);
        return arrow;
    };

    // v2.0.31: endpoint validator. Accepts {kind:'atom', atomId} or
    // {kind:'bond', bondId}; for either kind, the referenced object must
    // already exist on this molecule.
    Molecule.prototype._assertEndpointResolves = function(ref, label) {
        if (!ref || !ref.kind) {
            throw new Error('addCurlyArrow: ' + label + '.kind required');
        }
        if (ref.kind === 'atom') {
            if (!this.getAtom(ref.atomId)) {
                throw new Error('addCurlyArrow: ' + label + ' atom not found: ' + ref.atomId);
            }
        } else if (ref.kind === 'bond') {
            if (!this.getBond(ref.bondId)) {
                throw new Error('addCurlyArrow: ' + label + ' bond not found: ' + ref.bondId);
            }
        } else {
            throw new Error('addCurlyArrow: ' + label + '.kind must be "atom" or "bond"');
        }
    };

    Molecule.prototype.removeCurlyArrow = function(arrowId) {
        this.curlyArrows = this.curlyArrows.filter(function(a) { return a.id !== arrowId; });
    };

    Molecule.prototype.clearCurlyArrows = function() {
        this.curlyArrows = [];
    };

    Molecule.prototype.getCurlyArrow = function(arrowId) {
        for (var i = 0; i < this.curlyArrows.length; i++) {
            if (this.curlyArrows[i].id === arrowId) { return this.curlyArrows[i]; }
        }
        return null;
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
    // v1.4.1 (bug-fix #5): replaced recursive DFS with an iterative stack so
    // long-chain SMILES (>10k atoms — polymers, peptides) no longer overflow
    // the JS call stack. RDT.splitReactionSides, Layout.layout and Renderer
    // all sit on this hot path; a stack overflow here crashes Auto-map.
    Molecule.prototype.getComponents = function() {
        var visited = {};
        var components = [];
        for (var i = 0; i < this.atoms.length; i++) {
            var startId = this.atoms[i].id;
            if (visited[startId]) { continue; }
            var comp = [];
            var stack = [startId];
            while (stack.length > 0) {
                var id = stack.pop();
                if (visited[id]) { continue; }
                visited[id] = true;
                comp.push(id);
                var neighbors = this.getNeighbors(id);
                for (var k = 0; k < neighbors.length; k++) {
                    if (!visited[neighbors[k]]) { stack.push(neighbors[k]); }
                }
            }
            components.push(comp);
        }
        return components;
    };

    // =====================================================================
    // v2.0.16: reaction-equation helpers
    //
    // These are pure (no DOM, no mutation of `this`) so the bundle and the
    // headless test-suite can exercise them directly. The workbench page
    // controller calls them to render the "+ between fragments" glyphs and
    // the coefficient-collapsed equation in the Editor Output.
    // =====================================================================

    // Hill-system element counts for this molecule, including implicit H
    // (via calcHydrogens). Single home for formula counting in the editor
    // layer; the workbench delegates to this.
    Molecule.prototype.elementCounts = function() {
        var counts = {};
        for (var i = 0; i < this.atoms.length; i++) {
            var atom = this.atoms[i];
            counts[atom.symbol] = (counts[atom.symbol] || 0) + 1;
            var h = 0;
            try { h = this.calcHydrogens(atom.id); } catch (e) { h = 0; }
            if (h > 0) { counts.H = (counts.H || 0) + h; }
        }
        return counts;
    };

    // Hill-notation formula string from an element-count map: carbon first,
    // hydrogen second, then the remaining elements alphabetically. Counts of
    // 1 are written without a trailing digit (H2O, not H2O1).
    Molecule.formulaString = function(counts) {
        counts = counts || {};
        var keys = Object.keys(counts).sort();
        if (counts.C) {
            keys = ['C'].concat(counts.H ? ['H'] : [],
                keys.filter(function(k) { return k !== 'C' && k !== 'H'; }));
        }
        return keys.map(function(k) {
            return k + (counts[k] === 1 ? '' : counts[k]);
        }).join('');
    };

    // Compute "+" glyph positions between same-side disconnected fragments of
    // a drawn reaction. A fragment is a reactant when its centroid lies left
    // of the arrow midpoint, otherwise a product; each side is ordered left to
    // right and a "+" is placed in the horizontal gap between consecutive
    // same-side fragments. Returns [] when there is no arrow or fewer than two
    // fragments (a lone reactant -> lone product needs no "+").
    Molecule.prototype.computeReactionPlusSigns = function(arrow) {
        arrow = arrow || this.reactionArrow;
        if (!arrow || this.atoms.length === 0) { return []; }
        var comps = this.getComponents();
        if (comps.length < 2) { return []; }
        var axMid = (arrow.x1 + arrow.x2) / 2;
        var self = this;
        var frags = comps.map(function(ids) {
            var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
            for (var i = 0; i < ids.length; i++) {
                var a = self.getAtom(ids[i]);
                if (!a) { continue; }
                if (a.x < minX) { minX = a.x; }
                if (a.y < minY) { minY = a.y; }
                if (a.x > maxX) { maxX = a.x; }
                if (a.y > maxY) { maxY = a.y; }
            }
            return { minX: minX, maxX: maxX, cx: (minX + maxX) / 2, cy: (minY + maxY) / 2 };
        });
        var reactants = [], products = [];
        for (var f = 0; f < frags.length; f++) {
            (frags[f].cx < axMid ? reactants : products).push(frags[f]);
        }
        var byX = function(a, b) { return a.cx - b.cx; };
        reactants.sort(byX);
        products.sort(byX);
        var plus = [];
        var addBetween = function(side) {
            for (var i = 0; i + 1 < side.length; i++) {
                var a = side[i], b = side[i + 1];
                plus.push({ x: (a.maxX + b.minX) / 2, y: (a.cy + b.cy) / 2 });
            }
        };
        addBetween(reactants);
        addBetween(products);
        return plus;
    };

    // Assemble a human-readable, coefficient-collapsed reaction equation from
    // arrays of reactant/product Molecule objects. Identical species (same
    // keyFn(mol)) collapse into one term with a stoichiometric coefficient,
    // e.g. "2 H2 + O2 -> 2 H2O" (rendered with a Unicode arrow). `balanced`
    // compares summed element counts across both sides. keyFn defaults to the
    // molecular formula (so structural isomers merge); callers that can
    // canonicalise (the workbench) pass canonical SMILES for exact grouping.
    Molecule.reactionEquation = function(reactantMols, productMols, keyFn) {
        reactantMols = reactantMols || [];
        productMols = productMols || [];
        keyFn = keyFn || function(m) { return Molecule.formulaString(m.elementCounts()); };
        var ARROW = '→';

        function groupSide(mols, totalCounts) {
            var order = [];
            var byKey = {};
            for (var i = 0; i < mols.length; i++) {
                var m = mols[i];
                var counts = m.elementCounts();
                Object.keys(counts).forEach(function(k) {
                    totalCounts[k] = (totalCounts[k] || 0) + counts[k];
                });
                var key;
                try { key = keyFn(m); } catch (e) { key = ''; }
                if (!key) { key = Molecule.formulaString(counts); }
                if (byKey[key]) {
                    byKey[key].coeff += 1;
                } else {
                    byKey[key] = { coeff: 1, formula: Molecule.formulaString(counts), key: key };
                    order.push(key);
                }
            }
            return order.map(function(k) { return byKey[k]; });
        }

        function sideText(groups) {
            return groups.map(function(g) {
                return g.coeff > 1 ? (g.coeff + ' ' + g.formula) : g.formula;
            }).join(' + ');
        }

        function sideCharge(mols) {
            var total = 0;
            for (var i = 0; i < mols.length; i++) {
                var atoms = mols[i].atoms;
                for (var j = 0; j < atoms.length; j++) { total += (atoms[j].charge || 0); }
            }
            return total;
        }

        var rCounts = {}, pCounts = {};
        var rGroups = groupSide(reactantMols, rCounts);
        var pGroups = groupSide(productMols, pCounts);

        var allKeys = {};
        Object.keys(rCounts).forEach(function(k) { allKeys[k] = true; });
        Object.keys(pCounts).forEach(function(k) { allKeys[k] = true; });
        var elementsBalanced = Object.keys(allKeys).every(function(k) {
            return (rCounts[k] || 0) === (pCounts[k] || 0);
        });
        // A complete chemical balance conserves both mass (element counts) AND
        // formal charge, so a redox half-reaction (e.g. [Fe+2] -> [Fe+3]) reads
        // as unbalanced rather than mass-balanced-only "Balanced". Scope note:
        // this intentionally checks mass + charge only, not radical/spin count
        // (homolysis like Cl2 -> 2 Cl* is not flagged) — the SMILES parser does
        // not auto-perceive radicals, so that conservation law is unmodelled.
        var rCharge = sideCharge(reactantMols), pCharge = sideCharge(productMols);
        var chargeBalanced = rCharge === pCharge;
        var balanced = elementsBalanced && chargeBalanced;

        var rText = sideText(rGroups);
        var pText = sideText(pGroups);
        return {
            reactants: rGroups,
            products: pGroups,
            reactantText: rText,
            productText: pText,
            equation: (rText || '?') + ' ' + ARROW + ' ' + (pText || '?'),
            balanced: balanced,
            elementsBalanced: elementsBalanced,
            chargeBalanced: chargeBalanced,
            netCharge: { reactants: rCharge, products: pCharge }
        };
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
    Molecule.STANDARD_VALENCES = STANDARD_VALENCES;
    Molecule.getDefaultValence = getDefaultValence;
    Molecule.ATOMIC_NUMBERS = ATOMIC_NUMBERS;
    Molecule.perceiveAromaticity = perceiveAromaticity;
    Molecule.BOND_SINGLE = BOND_SINGLE;
    Molecule.BOND_DOUBLE = BOND_DOUBLE;
    Molecule.BOND_TRIPLE = BOND_TRIPLE;
    Molecule.STEREO_NONE = STEREO_NONE;
    Molecule.STEREO_WEDGE = STEREO_WEDGE;
    Molecule.STEREO_DASH = STEREO_DASH;
    Molecule.BOND_LENGTH = BOND_LENGTH;

    global.Molecule = Molecule;

})(typeof window !== 'undefined' ? window : this);
