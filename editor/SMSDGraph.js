/**
 * SMSDGraph.js — Molecular graph representation for substructure search
 *
 * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * Ported from SMSD upstream, mol_graph.hpp
 *
 * Builds flat arrays from a BIME Molecule for use by SMSDVF2 substructure matcher.
 */
(function() {
    'use strict';

    // ========================================================================
    // ATOMIC_NUMBERS — symbol to atomic number lookup
    // ========================================================================

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

    // ========================================================================
    // ChemOptions — default configuration for chemical matching
    // ========================================================================

    function ChemOptions(overrides) {
        this.matchAtomType = true;
        this.matchFormalCharge = true;
        this.useChirality = false;
        this.ringMatchesRingOnly = true;
        this.matchBondOrder = 'strict';   // 'strict' | 'loose' | 'any'
        this.aromaticityMode = 'flexible'; // 'strict' | 'flexible'
        this.tautomerAware = false;
        this.pH = 7.4;
        this.useTwoHopNLF = true;
        this.useThreeHopNLF = false;

        if (overrides) {
            for (var key in overrides) {
                if (overrides.hasOwnProperty(key) && this.hasOwnProperty(key)) {
                    this[key] = overrides[key];
                }
            }
        }
    }

    /** Strict matching profile. */
    ChemOptions.strict = function() {
        return new ChemOptions({
            matchBondOrder: 'strict',
            aromaticityMode: 'strict'
        });
    };

    /** Loose/compat profile. */
    ChemOptions.loose = function() {
        return new ChemOptions({
            matchBondOrder: 'loose',
            aromaticityMode: 'flexible',
            ringMatchesRingOnly: false
        });
    };

    /** Tautomer-aware profile. */
    ChemOptions.tautomer = function() {
        return new ChemOptions({
            tautomerAware: true,
            matchBondOrder: 'loose',
            aromaticityMode: 'flexible',
            matchFormalCharge: false,
            ringMatchesRingOnly: false
        });
    };

    // ========================================================================
    // Tautomer weights — pKa-informed relevance at pH 7.4 (Sitzmann 2010)
    // ========================================================================

    var TW_KETO_ENOL         = 0.98;
    var TW_AMIDE_IMIDIC      = 0.97;
    var TW_LACTAM_LACTIM     = 0.97;
    var TW_UREA              = 0.96;
    var TW_PYRIDINONE        = 0.95;
    var TW_THIOAMIDE         = 0.94;
    var TW_THIONE_THIOL      = 0.94;
    var TW_ENAMINONE          = 0.85;
    var TW_HYDROXAMIC        = 0.85;
    var TW_PHENOL_QUINONE    = 0.88;
    var TW_HYDROXYPYRIMIDINE = 0.90;
    var TW_PURINE_NH         = 0.78;
    var TW_DIKETONE_ENOL     = 0.72;
    var TW_NITROSO_OXIME     = 0.70;
    var TW_IMIDE             = 0.92;
    var TW_CYANAMIDE         = 0.80;
    var TW_AMIDINE           = 0.50;
    var TW_GUANIDINE         = 0.50;
    var TW_IMIDAZOLE_NH      = 0.50;
    var TW_TRIAZOLE_NH       = 0.50;
    var TW_TETRAZOLE_NH      = 0.50;
    var TW_NITRO_ACI         = 0.25;

    // Tier 2 tautomer weights (T16-T30)
    var TW_RING_CHAIN        = 0.65;
    var TW_ALLYL_SHIFT       = 0.60;
    var TW_SELENOL           = 0.90;
    var TW_SULFOXIDE         = 0.30; // disabled — not a tautomeric equilibrium
    var TW_PHOSPHORYL        = 0.70;
    var TW_NITRILE           = 0.20; // disabled — isomerism, not tautomerism
    var TW_NITROSAMINE       = 0.55;
    var TW_VINYL_AMINE       = 0.75;
    var TW_OXIME_NITROSO     = 0.68;
    var TW_BETA_KETO_ACID    = 0.70;
    var TW_ENOLATE           = 0.80;
    var TW_CYANURIC          = 0.45;
    var TW_DIIMINE           = 0.55;
    var TW_THIOKETONE        = 0.85;
    var TW_EXTENDED_ENOL     = 0.62;

    // ========================================================================
    // SMSDGraph constructor — builds flat arrays from a BIME Molecule
    // ========================================================================

    function SMSDGraph(mol) {
        var atoms = mol.atoms;
        var bonds = mol.bonds;
        var n = atoms.length;

        this.n = n;
        this.atomicNum = new Array(n);
        this.formalCharge = new Array(n);
        this.aromatic = new Array(n);
        this.ring = new Array(n);
        this.minRingSize = new Array(n);
        this.degree = new Array(n);
        this.label = new Array(n);
        this.neighbors = new Array(n);

        // Sparse bond property lookups
        this.bondOrd = {};
        this.bondRing = {};
        this.bondArom = {};

        // ID <-> index maps
        this.idToIdx = {};
        this.idxToId = new Array(n);
        this.atomId = new Array(n);  // original Molecule atom IDs for stable translation

        // Build index maps
        var i, j;
        for (i = 0; i < n; i++) {
            this.idToIdx[atoms[i].id] = i;
            this.idxToId[i] = atoms[i].id;
            this.atomId[i] = atoms[i].id;
            this.neighbors[i] = [];
        }

        // Ring detection via Molecule.findRings()
        var ringAtomSet = {};
        var ringBondSet = {};
        var rings = mol.findRings ? mol.findRings(20) : [];
        for (i = 0; i < rings.length; i++) {
            var ringAtoms = rings[i].atoms;
            for (j = 0; j < ringAtoms.length; j++) {
                ringAtomSet[ringAtoms[j]] = true;
                var nextJ = (j + 1) % ringAtoms.length;
                var a1 = ringAtoms[j], a2 = ringAtoms[nextJ];
                var bk1 = Math.min(a1, a2) + ',' + Math.max(a1, a2);
                ringBondSet[bk1] = true;
            }
        }

        // Compute minimum ring size per atom
        var minRingSizeMap = {};
        for (i = 0; i < rings.length; i++) {
            var rSize = rings[i].atoms.length;
            for (j = 0; j < rings[i].atoms.length; j++) {
                var rAtomId = rings[i].atoms[j];
                if (minRingSizeMap[rAtomId] === undefined || rSize < minRingSizeMap[rAtomId]) {
                    minRingSizeMap[rAtomId] = rSize;
                }
            }
        }

        // Aromatic bond detection via bond aromaticity or atom aromaticity
        var aromBondSet = {};

        // Fill atom arrays
        for (i = 0; i < n; i++) {
            var atom = atoms[i];
            var anum = ATOMIC_NUMBERS[atom.symbol] || 0;
            this.atomicNum[i] = anum;
            this.formalCharge[i] = atom.charge || 0;
            this.aromatic[i] = !!atom.aromatic;
            this.ring[i] = !!ringAtomSet[atom.id];
            this.minRingSize[i] = minRingSizeMap[atom.id] || 0;
            this.neighbors[i] = [];
        }

        // Build adjacency and bond properties
        for (i = 0; i < bonds.length; i++) {
            var bond = bonds[i];
            var idx1 = this.idToIdx[bond.atom1];
            var idx2 = this.idToIdx[bond.atom2];
            if (idx1 === undefined || idx2 === undefined) continue;

            this.neighbors[idx1].push(idx2);
            this.neighbors[idx2].push(idx1);

            var bkey = idx1 + ',' + idx2;
            var bkeyR = idx2 + ',' + idx1;
            var order = bond.type || 1;
            this.bondOrd[bkey] = order;
            this.bondOrd[bkeyR] = order;

            // Bond ring membership: check if both atoms are in a ring AND
            // the bond connects adjacent ring atoms
            var atomId1 = bond.atom1, atomId2 = bond.atom2;
            var idLo = Math.min(atomId1, atomId2);
            var idHi = Math.max(atomId1, atomId2);
            var inRing = !!ringBondSet[idLo + ',' + idHi];
            this.bondRing[bkey] = inRing;
            this.bondRing[bkeyR] = inRing;

            // Bond aromaticity: both atoms aromatic and bond in ring
            var isArom = this.aromatic[idx1] && this.aromatic[idx2] && inRing;
            this.bondArom[bkey] = isArom;
            this.bondArom[bkeyR] = isArom;
            if (isArom) {
                aromBondSet[bkey] = true;
                aromBondSet[bkeyR] = true;
            }
        }

        // Compute degrees and labels
        for (i = 0; i < n; i++) {
            this.degree[i] = this.neighbors[i].length;
            // label = (atomicNum << 2) | (aromatic ? 2 : 0) | (ring ? 1 : 0)
            this.label[i] = (this.atomicNum[i] << 2) |
                            (this.aromatic[i] ? 2 : 0) |
                            (this.ring[i] ? 1 : 0);
        }

        // Tautomer class and weight arrays (computed lazily)
        this.tautomerClass = null;
        this.tautomerWeight = null;

        // NLF caches (computed lazily)
        this._nlf1 = null;
        this._nlf2 = null;
    }

    // ========================================================================
    // Bond accessors
    // ========================================================================

    SMSDGraph.prototype.bondOrder = function(i, j) {
        return this.bondOrd[i + ',' + j] || 0;
    };

    SMSDGraph.prototype.bondInRing = function(i, j) {
        return !!this.bondRing[i + ',' + j];
    };

    SMSDGraph.prototype.bondAromatic = function(i, j) {
        return !!this.bondArom[i + ',' + j];
    };

    SMSDGraph.prototype.hasBond = function(i, j) {
        return this.bondOrder(i, j) !== 0;
    };

    // ========================================================================
    // NLF label: atomicNum + aromaticity (no ring bit — ring is one-directional)
    // ========================================================================

    function nlfLabel(graph, idx) {
        return (graph.atomicNum[idx] << 1) | (graph.aromatic[idx] ? 1 : 0);
    }

    // ========================================================================
    // buildNLF1 — 1-hop neighbor label frequencies
    // ========================================================================

    function buildNLF1(graph) {
        var n = graph.n;
        var result = new Array(n);
        var i, k, p, pairs, kl, kf, jj;

        for (i = 0; i < n; i++) {
            var freq = {};
            var nb = graph.neighbors[i];
            for (k = 0; k < nb.length; k++) {
                var lbl = nlfLabel(graph, nb[k]);
                freq[lbl] = (freq[lbl] || 0) + 1;
            }
            // Pack into sorted (label, freq) pair array
            var keys = [];
            for (var fk in freq) {
                if (freq.hasOwnProperty(fk)) keys.push(parseInt(fk, 10));
            }
            var arr = new Array(keys.length * 2);
            p = 0;
            for (k = 0; k < keys.length; k++) {
                arr[p++] = keys[k];
                arr[p++] = freq[keys[k]];
            }
            // Insertion sort by label
            pairs = keys.length;
            for (k = 1; k < pairs; k++) {
                kl = arr[k * 2];
                kf = arr[k * 2 + 1];
                jj = k - 1;
                while (jj >= 0 && arr[jj * 2] > kl) {
                    arr[(jj + 1) * 2] = arr[jj * 2];
                    arr[(jj + 1) * 2 + 1] = arr[jj * 2 + 1];
                    jj--;
                }
                arr[(jj + 1) * 2] = kl;
                arr[(jj + 1) * 2 + 1] = kf;
            }
            result[i] = arr;
        }
        return result;
    }

    // ========================================================================
    // buildNLF2 — 2-hop neighbor label frequencies
    // ========================================================================

    function buildNLF2(graph) {
        var n = graph.n;
        var result = new Array(n);
        var i, k, m, p, pairs, kl, kf, jj;

        for (i = 0; i < n; i++) {
            var seen = new Array(n);
            for (k = 0; k < n; k++) seen[k] = false;
            // Mark direct neighbors and self
            seen[i] = true;
            var nb = graph.neighbors[i];
            for (k = 0; k < nb.length; k++) seen[nb[k]] = true;
            // Walk 2-hop
            var freq = {};
            for (k = 0; k < nb.length; k++) {
                var nbk = graph.neighbors[nb[k]];
                for (m = 0; m < nbk.length; m++) {
                    var j = nbk[m];
                    if (!seen[j]) {
                        seen[j] = true;
                        var lbl = nlfLabel(graph, j);
                        freq[lbl] = (freq[lbl] || 0) + 1;
                    }
                }
            }
            // Pack into sorted array
            var keys = [];
            for (var fk in freq) {
                if (freq.hasOwnProperty(fk)) keys.push(parseInt(fk, 10));
            }
            var arr = new Array(keys.length * 2);
            p = 0;
            for (k = 0; k < keys.length; k++) {
                arr[p++] = keys[k];
                arr[p++] = freq[keys[k]];
            }
            pairs = keys.length;
            for (k = 1; k < pairs; k++) {
                kl = arr[k * 2];
                kf = arr[k * 2 + 1];
                jj = k - 1;
                while (jj >= 0 && arr[jj * 2] > kl) {
                    arr[(jj + 1) * 2] = arr[jj * 2];
                    arr[(jj + 1) * 2 + 1] = arr[jj * 2 + 1];
                    jj--;
                }
                arr[(jj + 1) * 2] = kl;
                arr[(jj + 1) * 2 + 1] = kf;
            }
            result[i] = arr;
        }
        return result;
    }

    // ========================================================================
    // getNLF1 / getNLF2 — lazy-cached NLF accessors
    // ========================================================================

    SMSDGraph.prototype.getNLF1 = function() {
        if (!this._nlf1) this._nlf1 = buildNLF1(this);
        return this._nlf1;
    };

    SMSDGraph.prototype.getNLF2 = function() {
        if (!this._nlf2) this._nlf2 = buildNLF2(this);
        return this._nlf2;
    };

    // ========================================================================
    // nlfOk — sorted-array merge check: every (label, freq) in fq must
    //         have a matching entry in ft with freq >= fq's freq
    // ========================================================================

    function nlfOk(fq, ft) {
        var fi = 0, ti = 0;
        var fqLen = fq.length, ftLen = ft.length;
        while (fi < fqLen) {
            var fLabel = fq[fi], fFreq = fq[fi + 1];
            while (ti < ftLen && ft[ti] < fLabel) ti += 2;
            if (ti >= ftLen || ft[ti] !== fLabel || ft[ti + 1] < fFreq) return false;
            fi += 2;
        }
        return true;
    }

    // ========================================================================
    // atomsCompat — check atom compatibility between query[i] and target[j]
    // ========================================================================

    function atomsCompat(g1, i, g2, j, opts) {
        // Tautomer relaxation: C/N/O can interchange within tautomeric regions
        var tautRelax = opts.tautomerAware &&
            g1.tautomerClass && g2.tautomerClass &&
            g1.tautomerClass[i] !== -1 && g2.tautomerClass[j] !== -1;

        if (tautRelax) {
            var aq = g1.atomicNum[i], at = g2.atomicNum[j];
            if (!((aq === 6 || aq === 7 || aq === 8) &&
                  (at === 6 || at === 7 || at === 8))) {
                tautRelax = false;
            }
            // Degree guard: genuine tautomers shift at most 1 proton
            if (tautRelax && aq === at &&
                Math.abs(g1.degree[i] - g2.degree[j]) > 1) {
                tautRelax = false;
            }
        }

        if (!tautRelax && opts.matchAtomType &&
            g1.atomicNum[i] !== g2.atomicNum[j]) return false;

        if (opts.matchFormalCharge &&
            g1.formalCharge[i] !== g2.formalCharge[j]) return false;

        if (opts.aromaticityMode === 'strict' &&
            g1.aromatic[i] !== g2.aromatic[j]) return false;

        if (opts.ringMatchesRingOnly &&
            g1.ring[i] && !g2.ring[j]) return false;

        return true;
    }

    // ========================================================================
    // bondsCompat — check bond compatibility between
    //   query bond (i1-j1) and target bond (i2-j2)
    // ========================================================================

    function bondsCompat(g1, i1, j1, g2, i2, j2, opts) {
        var qOrd = g1.bondOrder(i1, j1);
        var tOrd = g2.bondOrder(i2, j2);
        if (qOrd === 0 || tOrd === 0) return false;

        // Tautomer-aware: relax bond ORDER for tautomeric bonds
        var tautBondRelax = opts.tautomerAware &&
            g1.tautomerClass && g2.tautomerClass &&
            g1.tautomerClass[i1] !== -1 && g1.tautomerClass[j1] !== -1 &&
            g2.tautomerClass[i2] !== -1 && g2.tautomerClass[j2] !== -1;

        // Strict aromaticity
        if (opts.aromaticityMode === 'strict' &&
            g1.bondAromatic(i1, j1) !== g2.bondAromatic(i2, j2)) return false;

        // Ring-only bond check
        if (opts.ringMatchesRingOnly &&
            g1.bondInRing(i1, j1) && !g2.bondInRing(i2, j2)) return false;

        // Tautomeric bonds: order relaxed
        if (tautBondRelax) return true;

        if (opts.matchBondOrder === 'any') return true;
        if (qOrd === tOrd) return true;
        if (opts.matchBondOrder === 'loose') return true;

        // Flexible aromaticity: aromatic bonds match single/double
        if (opts.aromaticityMode === 'flexible') {
            var qa = g1.bondAromatic(i1, j1), ta = g2.bondAromatic(i2, j2);
            if ((qa && ta) ||
                (qa && (tOrd === 1 || tOrd === 2)) ||
                (ta && (qOrd === 1 || qOrd === 2))) {
                return true;
            }
        }

        return false;
    }

    // ========================================================================
    // computeTautomerClasses — 30 tautomeric transforms (Tier 1 + Tier 2)
    // ========================================================================

    function computeTautomerClasses(graph) {
        var n = graph.n;
        graph.tautomerClass = new Array(n);
        graph.tautomerWeight = new Array(n);
        var i;
        for (i = 0; i < n; i++) {
            graph.tautomerClass[i] = -1;
            graph.tautomerWeight[i] = 1.0;
        }
        var classId = 0;

        // Helper: assign class + weight to a set of atoms, merging existing classes
        function tc(atoms, cls, w) {
            var a;
            for (a = 0; a < atoms.length; a++) {
                if (graph.tautomerClass[atoms[a]] !== -1) {
                    cls = graph.tautomerClass[atoms[a]];
                    break;
                }
            }
            for (a = 0; a < atoms.length; a++) {
                graph.tautomerClass[atoms[a]] = cls;
                if (graph.tautomerWeight[atoms[a]] === 1.0) {
                    graph.tautomerWeight[atoms[a]] = w;
                }
            }
        }

        for (i = 0; i < n; i++) {
            if (graph.tautomerClass[i] !== -1) continue;

            var nb_i = graph.neighbors[i];
            var j, k, m, o, cls, w, nb_j, nb_k, nb_m;

            // ---- T1: Keto/enol  C=O <-> C-OH ----
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i]) {
                for (j = 0; j < nb_i.length; j++) {
                    var j1 = nb_i[j];
                    if (graph.atomicNum[j1] === 6 && graph.bondOrder(i, j1) === 2) {
                        var hasN = false;
                        nb_j = graph.neighbors[j1];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 7) { hasN = true; break; }
                        }
                        w = hasN ? TW_AMIDE_IMIDIC : TW_KETO_ENOL;
                        cls = -1;
                        for (k = 0; k < nb_j.length; k++) {
                            var k1 = nb_j[k];
                            if (k1 !== i && graph.atomicNum[k1] === 6 && !graph.aromatic[k1]) {
                                if (cls === -1) { cls = classId++; graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = w; }
                                if (graph.tautomerClass[k1] === -1) { graph.tautomerClass[k1] = cls; if (graph.tautomerWeight[k1] === 1.0) graph.tautomerWeight[k1] = w; }
                            }
                        }
                        for (k = 0; k < nb_j.length; k++) {
                            var k2 = nb_j[k];
                            if (k2 !== i && graph.atomicNum[k2] === 7) {
                                if (cls === -1) { cls = classId++; graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = w; }
                                if (graph.tautomerClass[k2] === -1) { graph.tautomerClass[k2] = cls; if (graph.tautomerWeight[k2] === 1.0) graph.tautomerWeight[k2] = w; }
                            }
                        }
                        if (cls !== -1) break;
                    }
                }
            }

            // ---- T2: Thione/thiol  C=S <-> C-SH ----
            if (graph.atomicNum[i] === 16 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j2 = nb_i[j];
                    if (graph.atomicNum[j2] === 6 && graph.bondOrder(i, j2) === 2) {
                        var hasN2 = false;
                        nb_j = graph.neighbors[j2];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 7) { hasN2 = true; break; }
                        }
                        w = hasN2 ? TW_THIOAMIDE : TW_THIONE_THIOL;
                        cls = -1;
                        for (k = 0; k < nb_j.length; k++) {
                            var k3 = nb_j[k];
                            if (k3 !== i && (graph.atomicNum[k3] === 6 || graph.atomicNum[k3] === 7)) {
                                if (cls === -1) { cls = classId++; graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = w; }
                                if (graph.tautomerClass[k3] === -1) { graph.tautomerClass[k3] = cls; if (graph.tautomerWeight[k3] === 1.0) graph.tautomerWeight[k3] = w; }
                            }
                        }
                        if (cls !== -1) break;
                    }
                }
            }

            // ---- T3: Nitroso/oxime  C=N-OH <-> C(=O)-NH ----
            if (graph.atomicNum[i] === 7 && !graph.ring[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j3 = nb_i[j];
                    if (graph.atomicNum[j3] === 6 && graph.bondOrder(i, j3) === 2) {
                        nb_j = graph.neighbors[j3];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8) {
                                cls = classId++;
                                tc([i, nb_j[k]], cls, TW_NITROSO_OXIME);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T4: Phenol/quinone  arom-C-OH <-> para-C=O ----
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j4 = nb_i[j];
                    if (graph.atomicNum[j4] === 6 && graph.aromatic[j4] && graph.ring[j4]) {
                        nb_j = graph.neighbors[j4];
                        for (var a = 0; a < nb_j.length; a++) {
                            var na = nb_j[a];
                            if (na !== i && graph.atomicNum[na] === 6 && graph.aromatic[na] && graph.ring[na]) {
                                var nb_a = graph.neighbors[na];
                                for (var b = 0; b < nb_a.length; b++) {
                                    var nb2 = nb_a[b];
                                    if (nb2 !== j4 && graph.atomicNum[nb2] === 6 && graph.aromatic[nb2] && graph.ring[nb2]) {
                                        var nb_b = graph.neighbors[nb2];
                                        for (var c = 0; c < nb_b.length; c++) {
                                            var nc = nb_b[c];
                                            if (nc !== na && nc !== j4 && graph.atomicNum[nc] === 6 && graph.aromatic[nc] && graph.ring[nc]) {
                                                cls = classId++;
                                                tc([i, j4, nc], cls, TW_PHENOL_QUINONE);
                                                var nb_c = graph.neighbors[nc];
                                                for (var d = 0; d < nb_c.length; d++) {
                                                    if (graph.atomicNum[nb_c[d]] === 8 && !graph.aromatic[nb_c[d]] && graph.tautomerClass[nb_c[d]] === -1) {
                                                        graph.tautomerClass[nb_c[d]] = cls;
                                                        graph.tautomerWeight[nb_c[d]] = TW_PHENOL_QUINONE;
                                                    }
                                                }
                                                break;
                                            }
                                        }
                                        if (graph.tautomerClass[i] !== -1) break;
                                    }
                                }
                                if (graph.tautomerClass[i] !== -1) break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T5: 1,3-diketone  C(=O)-C-C(=O) enolisation ----
            // Both oxygens must have C=O double bonds.
            // 1,3-path-length check: exclude 1,2-diketones (j5 and m5 must NOT
            // be directly bonded — ensures exactly one bridging carbon).
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j5 = nb_i[j];
                    if (graph.atomicNum[j5] === 6 && graph.bondOrder(i, j5) === 2) {
                        nb_j = graph.neighbors[j5];
                        for (k = 0; k < nb_j.length; k++) {
                            var bridge = nb_j[k];
                            if (bridge !== i && graph.atomicNum[bridge] === 6 && !graph.aromatic[bridge]) {
                                nb_k = graph.neighbors[bridge];
                                for (m = 0; m < nb_k.length; m++) {
                                    var m5 = nb_k[m];
                                    if (m5 !== j5 && graph.atomicNum[m5] === 6) {
                                        // 1,3-path check: j5 and m5 must not be directly bonded
                                        // (that would be a 1,2-diketone O=C-C=O, not 1,3)
                                        if (graph.bondOrder(j5, m5) !== 0) { continue; }
                                        nb_m = graph.neighbors[m5];
                                        for (o = 0; o < nb_m.length; o++) {
                                            var o2 = nb_m[o];
                                            if (o2 !== bridge && graph.atomicNum[o2] === 8 && !graph.aromatic[o2] && graph.bondOrder(m5, o2) === 2) {
                                                cls = (graph.tautomerClass[o2] !== -1) ? graph.tautomerClass[o2] : classId++;
                                                w = TW_DIKETONE_ENOL;
                                                if (graph.tautomerClass[i] === -1)      { graph.tautomerClass[i] = cls;      graph.tautomerWeight[i] = w; }
                                                if (graph.tautomerClass[bridge] === -1) { graph.tautomerClass[bridge] = cls; graph.tautomerWeight[bridge] = w; }
                                                if (graph.tautomerClass[o2] === -1)     { graph.tautomerClass[o2] = cls;     graph.tautomerWeight[o2] = w; }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ---- T6: Pyridinone / lactam ----
            if (graph.atomicNum[i] === 7 && graph.aromatic[i] && graph.ring[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j6 = nb_i[j];
                    if (graph.atomicNum[j6] === 6 && graph.aromatic[j6]) {
                        nb_j = graph.neighbors[j6];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && !graph.aromatic[nb_j[k]] && graph.tautomerClass[nb_j[k]] === -1) {
                                cls = classId++;
                                tc([i, nb_j[k]], cls, TW_PYRIDINONE);
                            }
                        }
                    }
                }
            }

            // ---- T7: Imidazole / pyrazole (requires 5-membered ring) ----
            if (graph.atomicNum[i] === 7 && graph.aromatic[i] && graph.ring[i] && graph.minRingSize[i] === 5 && graph.tautomerClass[i] === -1) {
                // Direct N-N bond (pyrazole)
                for (j = 0; j < nb_i.length; j++) {
                    var j7 = nb_i[j];
                    if (graph.atomicNum[j7] === 7 && graph.aromatic[j7] && graph.ring[j7] && graph.minRingSize[j7] === 5 && graph.tautomerClass[j7] === -1) {
                        cls = classId++;
                        tc([i, j7], cls, TW_IMIDAZOLE_NH);
                        break;
                    }
                }
                // N-C-N through aromatic C (imidazole)
                if (graph.tautomerClass[i] === -1) {
                    for (j = 0; j < nb_i.length; j++) {
                        var j7b = nb_i[j];
                        if (graph.atomicNum[j7b] === 6 && graph.aromatic[j7b] && graph.ring[j7b] && graph.minRingSize[j7b] === 5) {
                            nb_j = graph.neighbors[j7b];
                            for (k = 0; k < nb_j.length; k++) {
                                if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 7 && graph.aromatic[nb_j[k]] && graph.ring[nb_j[k]] && graph.minRingSize[nb_j[k]] === 5 && graph.tautomerClass[nb_j[k]] === -1) {
                                    cls = classId++;
                                    tc([i, nb_j[k]], cls, TW_IMIDAZOLE_NH);
                                    break;
                                }
                            }
                            if (graph.tautomerClass[i] !== -1) break;
                        }
                    }
                }
            }

            // ---- T8: Amidine  N=C-N <-> N-C=N ----
            if (graph.atomicNum[i] === 7 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j8 = nb_i[j];
                    if (graph.atomicNum[j8] === 6 && graph.bondOrder(i, j8) === 2) {
                        nb_j = graph.neighbors[j8];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 7 && !graph.aromatic[nb_j[k]]) {
                                cls = classId++;
                                tc([i, j8, nb_j[k]], cls, TW_AMIDINE);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T9: Guanidine  C connected to >= 3 N ----
            if (graph.atomicNum[i] === 6 && graph.tautomerClass[i] === -1) {
                var nNbrs = [];
                for (j = 0; j < nb_i.length; j++) {
                    if (graph.atomicNum[nb_i[j]] === 7 && !graph.aromatic[nb_i[j]]) nNbrs.push(nb_i[j]);
                }
                if (nNbrs.length >= 3) {
                    cls = classId++;
                    graph.tautomerClass[i] = cls;
                    graph.tautomerWeight[i] = TW_GUANIDINE;
                    for (j = 0; j < nNbrs.length; j++) {
                        if (graph.tautomerClass[nNbrs[j]] === -1) {
                            graph.tautomerClass[nNbrs[j]] = cls;
                            graph.tautomerWeight[nNbrs[j]] = TW_GUANIDINE;
                        }
                    }
                }
            }

            // ---- T10: Urea  O=C(-N)(-N) ----
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j10 = nb_i[j];
                    if (graph.atomicNum[j10] === 6 && graph.bondOrder(i, j10) === 2) {
                        var nn = [];
                        nb_j = graph.neighbors[j10];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 7) nn.push(nb_j[k]);
                        }
                        if (nn.length >= 2) {
                            cls = classId++;
                            graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = TW_UREA;
                            graph.tautomerClass[j10] = cls; graph.tautomerWeight[j10] = TW_UREA;
                            for (k = 0; k < nn.length; k++) {
                                if (graph.tautomerClass[nn[k]] === -1) {
                                    graph.tautomerClass[nn[k]] = cls;
                                    graph.tautomerWeight[nn[k]] = TW_UREA;
                                }
                            }
                        }
                        break;
                    }
                }
            }

            // ---- T11: Enaminone  N-C=C-C=O ----
            if (graph.atomicNum[i] === 7 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                var found11 = false;
                for (j = 0; j < nb_i.length && !found11; j++) {
                    var j11 = nb_i[j];
                    if (graph.atomicNum[j11] === 6 && graph.bondOrder(i, j11) === 1) {
                        nb_j = graph.neighbors[j11];
                        for (k = 0; k < nb_j.length && !found11; k++) {
                            var k11 = nb_j[k];
                            if (k11 !== i && graph.atomicNum[k11] === 6 && graph.bondOrder(j11, k11) === 2) {
                                nb_k = graph.neighbors[k11];
                                for (m = 0; m < nb_k.length && !found11; m++) {
                                    var m11 = nb_k[m];
                                    if (m11 !== j11 && graph.atomicNum[m11] === 6) {
                                        nb_m = graph.neighbors[m11];
                                        for (o = 0; o < nb_m.length && !found11; o++) {
                                            if (graph.atomicNum[nb_m[o]] === 8 && graph.bondOrder(m11, nb_m[o]) === 2) {
                                                cls = classId++;
                                                tc([i, j11, k11, m11, nb_m[o]], cls, TW_ENAMINONE);
                                                found11 = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ---- T12: Imide  C(=O)-N-C(=O) ----
            if (graph.atomicNum[i] === 7 && graph.tautomerClass[i] === -1) {
                var carbonyls = [];
                for (j = 0; j < nb_i.length; j++) {
                    var j12 = nb_i[j];
                    if (graph.atomicNum[j12] === 6) {
                        nb_j = graph.neighbors[j12];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && graph.bondOrder(j12, nb_j[k]) === 2) {
                                carbonyls.push(j12);
                                break;
                            }
                        }
                    }
                }
                if (carbonyls.length >= 2) {
                    cls = classId++;
                    graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = TW_IMIDE;
                    for (j = 0; j < carbonyls.length; j++) {
                        if (graph.tautomerClass[carbonyls[j]] === -1) {
                            graph.tautomerClass[carbonyls[j]] = cls;
                            graph.tautomerWeight[carbonyls[j]] = TW_IMIDE;
                        }
                    }
                }
            }

            // ---- T13: Hydroxamic acid  N(-OH)-C=O <-> NH-C(=O)-OH ----
            if (graph.atomicNum[i] === 7 && graph.tautomerClass[i] === -1) {
                var oh = -1;
                for (j = 0; j < nb_i.length; j++) {
                    if (graph.atomicNum[nb_i[j]] === 8 && !graph.aromatic[nb_i[j]] && graph.bondOrder(i, nb_i[j]) === 1) {
                        oh = nb_i[j]; break;
                    }
                }
                if (oh !== -1) {
                    for (j = 0; j < nb_i.length; j++) {
                        var j13 = nb_i[j];
                        if (graph.atomicNum[j13] === 6) {
                            nb_j = graph.neighbors[j13];
                            for (k = 0; k < nb_j.length; k++) {
                                if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && graph.bondOrder(j13, nb_j[k]) === 2) {
                                    cls = classId++;
                                    tc([i, oh, j13, nb_j[k]], cls, TW_HYDROXAMIC);
                                    break;
                                }
                            }
                            if (graph.tautomerClass[i] !== -1) break;
                        }
                    }
                }
            }

            // ---- T14: Hydroxypyrimidine / purine ----
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j14 = nb_i[j];
                    if (graph.atomicNum[j14] === 6 && graph.aromatic[j14] && graph.ring[j14] && graph.bondOrder(i, j14) === 1) {
                        nb_j = graph.neighbors[j14];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 7 && graph.aromatic[nb_j[k]] && graph.ring[nb_j[k]]) {
                                cls = classId++;
                                tc([i, j14, nb_j[k]], cls, TW_HYDROXYPYRIMIDINE);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T15: Tetrazole / triazole N-H ----
            if (graph.atomicNum[i] === 7 && graph.aromatic[i] && graph.ring[i] && graph.tautomerClass[i] === -1) {
                var adjN = [];
                for (j = 0; j < nb_i.length; j++) {
                    if (graph.atomicNum[nb_i[j]] === 7 && graph.aromatic[nb_i[j]] && graph.ring[nb_i[j]]) {
                        adjN.push(nb_i[j]);
                    }
                }
                if (adjN.length >= 2) {
                    for (j = 0; j < adjN.length; j++) {
                        if (graph.tautomerClass[adjN[j]] === -1) {
                            cls = classId++;
                            tc([i, adjN[j]], cls, TW_TETRAZOLE_NH);
                            break;
                        }
                    }
                }
            }

            // ================================================================
            // Tier 2 tautomers (T16-T30)
            // ================================================================

            // ---- T16: Ring-chain tautomer  cyclic hemiketal/hemiacetal ----
            // O in ring bonded to C also bonded to exocyclic OH
            if (graph.atomicNum[i] === 8 && graph.ring[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j16 = nb_i[j];
                    if (graph.atomicNum[j16] === 6 && graph.ring[j16]) {
                        nb_j = graph.neighbors[j16];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && !graph.ring[nb_j[k]] && graph.bondOrder(j16, nb_j[k]) === 1) {
                                cls = classId++;
                                tc([i, j16, nb_j[k]], cls, TW_RING_CHAIN);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T17: Allyl 1,3-shift  C=C-XH <-> XH-C=C ----
            // X = N or O, with alternating single/double through carbon chain.
            // Bond-connectivity validation: ensure continuous path i-j17-k17-m17.
            // H-count validation: both heteroatom endpoints must carry at least
            // one implicit H (proton source for the tautomeric shift).
            if ((graph.atomicNum[i] === 7 || graph.atomicNum[i] === 8) && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                // H-count check on atom i
                var iZ17 = graph.atomicNum[i];
                var iVal17 = (iZ17 === 7) ? 3 : (iZ17 === 8) ? 2 : 0;
                var iDeg17 = graph.neighbors[i].length;
                var iChg17 = graph.formalCharge[i] || 0;
                var iImplH17 = Math.max(0, iVal17 - iDeg17 + iChg17);
                if (iImplH17 > 0) {
                    for (j = 0; j < nb_i.length; j++) {
                        var j17 = nb_i[j];
                        if (graph.atomicNum[j17] === 6 && graph.bondOrder(i, j17) === 1) {
                            nb_j = graph.neighbors[j17];
                            for (k = 0; k < nb_j.length; k++) {
                                var k17 = nb_j[k];
                                if (k17 !== i && graph.atomicNum[k17] === 6 && graph.bondOrder(j17, k17) === 2) {
                                    // Bond-connectivity: j17 and k17 must be bonded (already checked above)
                                    nb_k = graph.neighbors[k17];
                                    for (m = 0; m < nb_k.length; m++) {
                                        var m17 = nb_k[m];
                                        if (m17 !== j17 && (graph.atomicNum[m17] === 7 || graph.atomicNum[m17] === 8) && !graph.aromatic[m17]) {
                                            // Bond-connectivity: k17 and m17 must be bonded
                                            if (graph.bondOrder(k17, m17) === 0) { continue; }
                                            // H-count check on the other heteroatom endpoint
                                            var mZ17 = graph.atomicNum[m17];
                                            var mVal17 = (mZ17 === 7) ? 3 : (mZ17 === 8) ? 2 : 0;
                                            var mDeg17 = graph.neighbors[m17].length;
                                            var mChg17 = graph.formalCharge[m17] || 0;
                                            var mImplH17 = Math.max(0, mVal17 - mDeg17 + mChg17);
                                            if (mImplH17 <= 0) { continue; }
                                            cls = classId++;
                                            tc([i, j17, k17, m17], cls, TW_ALLYL_SHIFT);
                                            break;
                                        }
                                    }
                                    if (graph.tautomerClass[i] !== -1) break;
                                }
                            }
                            if (graph.tautomerClass[i] !== -1) break;
                        }
                    }
                }
            }

            // ---- T18: Sulfoxide (DISABLED - not a tautomeric equilibrium) ----
            // S(=O) <-> S-OH — kept for completeness but not applied
            // if (graph.atomicNum[i] === 16 && ...) { ... }

            // ---- T19: Selenol/selenone  C=Se <-> C-SeH ----
            if (graph.atomicNum[i] === 34 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j19 = nb_i[j];
                    if (graph.atomicNum[j19] === 6 && graph.bondOrder(i, j19) === 2) {
                        nb_j = graph.neighbors[j19];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && (graph.atomicNum[nb_j[k]] === 6 || graph.atomicNum[nb_j[k]] === 7)) {
                                if (graph.tautomerClass[nb_j[k]] === -1) {
                                    cls = classId++;
                                    tc([i, j19, nb_j[k]], cls, TW_SELENOL);
                                    break;
                                }
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T20: Nitrile/isonitrile (DISABLED - isomerism, not tautomerism) ----
            // C#N <-> C=N — functional group isomerism, not proton transfer
            // if (graph.atomicNum[i] === 7 && ...) { ... }

            // ---- T21: Phosphoryl  P(=O)(OH) <-> P(OH)3 ----
            if (graph.atomicNum[i] === 15 && graph.tautomerClass[i] === -1) {
                var oP = [];
                for (j = 0; j < nb_i.length; j++) {
                    if (graph.atomicNum[nb_i[j]] === 8) oP.push(nb_i[j]);
                }
                if (oP.length >= 2) {
                    var hasDouble = false;
                    for (j = 0; j < oP.length; j++) {
                        if (graph.bondOrder(i, oP[j]) === 2) { hasDouble = true; break; }
                    }
                    if (hasDouble) {
                        cls = classId++;
                        graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = TW_PHOSPHORYL;
                        for (j = 0; j < oP.length; j++) {
                            if (graph.tautomerClass[oP[j]] === -1) {
                                graph.tautomerClass[oP[j]] = cls;
                                graph.tautomerWeight[oP[j]] = TW_PHOSPHORYL;
                            }
                        }
                    }
                }
            }

            // ---- T22: Nitrosamine  N-N=O <-> N=N-OH ----
            if (graph.atomicNum[i] === 7 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j22 = nb_i[j];
                    if (graph.atomicNum[j22] === 7 && !graph.aromatic[j22]) {
                        nb_j = graph.neighbors[j22];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && graph.bondOrder(j22, nb_j[k]) === 2) {
                                cls = classId++;
                                tc([i, j22, nb_j[k]], cls, TW_NITROSAMINE);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T23: Vinyl amine  C=C-NH2 <-> C(-NH)-CH ----
            if (graph.atomicNum[i] === 7 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j23 = nb_i[j];
                    if (graph.atomicNum[j23] === 6 && graph.bondOrder(i, j23) === 1) {
                        nb_j = graph.neighbors[j23];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 6 && graph.bondOrder(j23, nb_j[k]) === 2) {
                                cls = classId++;
                                tc([i, j23, nb_j[k]], cls, TW_VINYL_AMINE);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T24: Oxime/nitroso  C=N-OH <-> C(=O)-NH (extended) ----
            // Similar to T3 but includes ring N
            if (graph.atomicNum[i] === 7 && graph.ring[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j24 = nb_i[j];
                    if (graph.atomicNum[j24] === 6 && graph.bondOrder(i, j24) === 2) {
                        nb_j = graph.neighbors[j24];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && graph.bondOrder(j24, nb_j[k]) === 1) {
                                cls = classId++;
                                tc([i, nb_j[k]], cls, TW_OXIME_NITROSO);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T25: Beta-keto acid  HOOC-CH2-C(=O) enolisation ----
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j25 = nb_i[j];
                    if (graph.atomicNum[j25] === 6 && graph.bondOrder(i, j25) === 2) {
                        nb_j = graph.neighbors[j25];
                        var hasOH25 = false;
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && graph.bondOrder(j25, nb_j[k]) === 1) { hasOH25 = true; break; }
                        }
                        if (hasOH25) {
                            for (k = 0; k < nb_j.length; k++) {
                                var bridge25 = nb_j[k];
                                if (bridge25 !== i && graph.atomicNum[bridge25] === 6 && !graph.aromatic[bridge25]) {
                                    nb_k = graph.neighbors[bridge25];
                                    for (m = 0; m < nb_k.length; m++) {
                                        var m25 = nb_k[m];
                                        if (m25 !== j25 && graph.atomicNum[m25] === 6) {
                                            nb_m = graph.neighbors[m25];
                                            for (o = 0; o < nb_m.length; o++) {
                                                if (nb_m[o] !== bridge25 && graph.atomicNum[nb_m[o]] === 8 && graph.bondOrder(m25, nb_m[o]) === 2) {
                                                    cls = classId++;
                                                    tc([i, j25, bridge25, m25, nb_m[o]], cls, TW_BETA_KETO_ACID);
                                                    break;
                                                }
                                            }
                                            if (graph.tautomerClass[i] !== -1) break;
                                        }
                                    }
                                    if (graph.tautomerClass[i] !== -1) break;
                                }
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T26: Enolate  C(=O)-C <-> C(-O^-)-C= (deprotonated enol) ----
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.formalCharge[i] === -1 && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j26 = nb_i[j];
                    if (graph.atomicNum[j26] === 6 && graph.bondOrder(i, j26) === 1) {
                        nb_j = graph.neighbors[j26];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 6) {
                                cls = classId++;
                                tc([i, j26, nb_j[k]], cls, TW_ENOLATE);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T27: Cyanuric acid / isocyanuric acid ----
            // Ring with alternating N-C(=O) pattern (6-membered, 3 C=O)
            if (graph.atomicNum[i] === 7 && graph.ring[i] && graph.minRingSize[i] === 6 && graph.tautomerClass[i] === -1) {
                var carbonylCount27 = 0;
                for (j = 0; j < nb_i.length; j++) {
                    var j27 = nb_i[j];
                    if (graph.atomicNum[j27] === 6 && graph.ring[j27]) {
                        nb_j = graph.neighbors[j27];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 8 && graph.bondOrder(j27, nb_j[k]) === 2) {
                                carbonylCount27++;
                                break;
                            }
                        }
                    }
                }
                if (carbonylCount27 >= 2) {
                    cls = classId++;
                    graph.tautomerClass[i] = cls; graph.tautomerWeight[i] = TW_CYANURIC;
                    for (j = 0; j < nb_i.length; j++) {
                        if (graph.atomicNum[nb_i[j]] === 6 && graph.ring[nb_i[j]] && graph.tautomerClass[nb_i[j]] === -1) {
                            graph.tautomerClass[nb_i[j]] = cls;
                            graph.tautomerWeight[nb_i[j]] = TW_CYANURIC;
                        }
                    }
                }
            }

            // ---- T28: 1,3-diimine  N=C-C=N <-> NH-C=C-NH ----
            if (graph.atomicNum[i] === 7 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j28 = nb_i[j];
                    if (graph.atomicNum[j28] === 6 && graph.bondOrder(i, j28) === 2) {
                        nb_j = graph.neighbors[j28];
                        for (k = 0; k < nb_j.length; k++) {
                            var k28 = nb_j[k];
                            if (k28 !== i && graph.atomicNum[k28] === 6) {
                                nb_k = graph.neighbors[k28];
                                for (m = 0; m < nb_k.length; m++) {
                                    if (nb_k[m] !== j28 && graph.atomicNum[nb_k[m]] === 7 && !graph.aromatic[nb_k[m]] && graph.bondOrder(k28, nb_k[m]) === 2) {
                                        cls = classId++;
                                        tc([i, j28, k28, nb_k[m]], cls, TW_DIIMINE);
                                        break;
                                    }
                                }
                                if (graph.tautomerClass[i] !== -1) break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T29: Thioketone enol  C=S-C <-> C(SH)=C ----
            if (graph.atomicNum[i] === 16 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j29 = nb_i[j];
                    if (graph.atomicNum[j29] === 6 && graph.bondOrder(i, j29) === 1) {
                        nb_j = graph.neighbors[j29];
                        for (k = 0; k < nb_j.length; k++) {
                            if (nb_j[k] !== i && graph.atomicNum[nb_j[k]] === 6 && graph.bondOrder(j29, nb_j[k]) === 2) {
                                cls = classId++;
                                tc([i, j29, nb_j[k]], cls, TW_THIOKETONE);
                                break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }

            // ---- T30: Extended enol  C=C-C=C-OH <-> C(-OH)-C=C-C=O ----
            // 1,5-proton shift through conjugated system
            if (graph.atomicNum[i] === 8 && !graph.aromatic[i] && graph.tautomerClass[i] === -1) {
                for (j = 0; j < nb_i.length; j++) {
                    var j30 = nb_i[j];
                    if (graph.atomicNum[j30] === 6 && graph.bondOrder(i, j30) === 1) {
                        nb_j = graph.neighbors[j30];
                        for (k = 0; k < nb_j.length; k++) {
                            var k30 = nb_j[k];
                            if (k30 !== i && graph.atomicNum[k30] === 6 && graph.bondOrder(j30, k30) === 2) {
                                nb_k = graph.neighbors[k30];
                                for (m = 0; m < nb_k.length; m++) {
                                    var m30 = nb_k[m];
                                    if (m30 !== j30 && graph.atomicNum[m30] === 6 && graph.bondOrder(k30, m30) === 1) {
                                        var nb_m30 = graph.neighbors[m30];
                                        for (o = 0; o < nb_m30.length; o++) {
                                            var o30 = nb_m30[o];
                                            if (o30 !== k30 && graph.atomicNum[o30] === 6 && graph.bondOrder(m30, o30) === 2) {
                                                var nb_o30 = graph.neighbors[o30];
                                                for (var p = 0; p < nb_o30.length; p++) {
                                                    if (nb_o30[p] !== m30 && graph.atomicNum[nb_o30[p]] === 8 && !graph.aromatic[nb_o30[p]] && graph.bondOrder(o30, nb_o30[p]) === 2) {
                                                        cls = classId++;
                                                        tc([i, j30, k30, m30, o30, nb_o30[p]], cls, TW_EXTENDED_ENOL);
                                                        break;
                                                    }
                                                }
                                                if (graph.tautomerClass[i] !== -1) break;
                                            }
                                        }
                                        if (graph.tautomerClass[i] !== -1) break;
                                    }
                                }
                                if (graph.tautomerClass[i] !== -1) break;
                            }
                        }
                        if (graph.tautomerClass[i] !== -1) break;
                    }
                }
            }
        }

        // Henderson-Hasselbalch pH adjustment (shift clamped to [-0.5, 0.5])
        if (Math.abs(graph._pH - 7.4) >= 0.01) {
            var rawShift = 0.1 * (graph._pH - 7.4);
            var shift = Math.max(-0.5, Math.min(0.5, rawShift));
            for (i = 0; i < n; i++) {
                if (graph.tautomerClass[i] === -1 || graph.tautomerWeight[i] === 1.0) continue;
                graph.tautomerWeight[i] = Math.max(0.1, Math.min(0.99, graph.tautomerWeight[i] + shift));
            }
        }
    }

    /** Ensure tautomer classes are computed. */
    SMSDGraph.prototype.ensureTautomers = function(pH) {
        this._pH = (pH !== undefined) ? pH : 7.4;
        if (!this.tautomerClass) {
            computeTautomerClasses(this);
        }
    };

    // ========================================================================
    // atomsCompatFast — quick compatibility check (atomic number + degree only)
    // Avoids full chemistry evaluation; suitable for candidate pre-filtering.
    // ========================================================================

    function atomsCompatFast(g1, i, g2, j) {
        if (g1.atomicNum[i] !== g2.atomicNum[j]) return false;
        if (g1.degree[i] > g2.degree[j]) return false;
        return true;
    }

    // ========================================================================
    // isAnilineN — check if atom i in graph g is an aniline nitrogen
    // (N bonded directly to an aromatic ring carbon, not in ring itself)
    // Aniline N should be excluded from positive ionisable classification.
    // ========================================================================

    function isAnilineN(g, i) {
        if (g.atomicNum[i] !== 7 || g.aromatic[i] || g.ring[i]) return false;
        var nb = g.neighbors[i];
        for (var k = 0; k < nb.length; k++) {
            if (g.atomicNum[nb[k]] === 6 && g.aromatic[nb[k]] && g.ring[nb[k]]) {
                return true;
            }
        }
        return false;
    }

    // ========================================================================
    // isEsterO — check if atom i in graph g is an ester oxygen
    // (O bonded to C which also has C=O, i.e. -O-C(=O)-)
    // Ester O should be excluded from negative ionisable classification.
    // ========================================================================

    function isEsterO(g, i) {
        if (g.atomicNum[i] !== 8) return false;
        var nb = g.neighbors[i];
        for (var k = 0; k < nb.length; k++) {
            var c = nb[k];
            if (g.atomicNum[c] !== 6) continue;
            if (g.bondOrder(i, c) !== 1) continue;
            var nbC = g.neighbors[c];
            for (var m = 0; m < nbC.length; m++) {
                if (nbC[m] !== i && g.atomicNum[nbC[m]] === 8 && g.bondOrder(c, nbC[m]) === 2) {
                    return true;
                }
            }
        }
        return false;
    }

    // ========================================================================
    // isPyridineN — check if atom i in graph g is a pyridine-type nitrogen
    // (aromatic N in ring with degree 2; H-bond acceptor, not donor)
    // ========================================================================

    function isPyridineN(g, i) {
        return g.atomicNum[i] === 7 && g.aromatic[i] && g.ring[i] && g.degree[i] === 2;
    }

    // ========================================================================
    // Export
    // ========================================================================

    window.SMSDGraph = {
        SMSDGraph: SMSDGraph,
        ChemOptions: ChemOptions,
        ATOMIC_NUMBERS: ATOMIC_NUMBERS,
        atomsCompat: atomsCompat,
        atomsCompatFast: atomsCompatFast,
        bondsCompat: bondsCompat,
        nlfOk: nlfOk,
        buildNLF1: buildNLF1,
        buildNLF2: buildNLF2,
        nlfLabel: nlfLabel,
        computeTautomerClasses: computeTautomerClasses,
        isAnilineN: isAnilineN,
        isEsterO: isEsterO,
        isPyridineN: isPyridineN
    };

})();
