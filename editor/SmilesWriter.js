/**
 * SmilesWriter.js — Molecule graph -> canonical SMILES string
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Generates canonical SMILES from a Molecule graph.
 * Features:
 *   - Morgan-algorithm canonical atom ordering
 *   - Aromatic perception using Hueckel 4n+2 rule on SSSR rings
 *   - Correct implicit-H count (valence - bond_order_sum - |charge|)
 *   - Bracket atoms only when required (non-organic, charge, isotope,
 *     explicit H different from default, map number, chirality)
 *   - Ring closure numbering (digits 1-9, then %10, %11, ...)
 *   - Stereochemistry output (@ / @@ tetrahedral, / \ E/Z double bonds)
 *   - Reaction SMILES (>>)
 */
(function(global) {
    'use strict';

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    var ORGANIC = { 'B': true, 'C': true, 'N': true, 'O': true,
                    'P': true, 'S': true, 'F': true, 'Cl': true,
                    'Br': true, 'I': true };

    // Standard valences for implicit-H
    var STANDARD_VALENCES = {
        'B': [3], 'C': [4], 'N': [3, 5], 'O': [2],
        'P': [3, 5], 'S': [2, 4, 6], 'F': [1], 'Cl': [1],
        'Br': [1], 'I': [1], 'H': [1], 'Si': [4],
        'Se': [2, 4, 6], 'As': [3, 5], 'Te': [2, 4, 6]
    };

    // Aromatic atoms that can be written in lowercase without brackets
    var AROMATIC_ORGANIC = { 'B': true, 'C': true, 'N': true, 'O': true,
                             'P': true, 'S': true };

    // Atomic numbers for canonical tie-breaking
    var ATOMIC_NUMBERS = {
        'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
        'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
        'K':19,'Ca':20,'Fe':26,'Cu':29,'Zn':30,'As':33,'Se':34,'Br':35,'I':53
    };

    // Aromatic valences for implicit-H calculation in the writer.
    // Must match the AROMATIC_VALENCE table in Molecule.js so that the
    // writer's expectedH matches what the reader's calcHydrogens returns.
    var AROMATIC_WRITER_VALENCE = {
        'C': 3, 'N': 2, 'O': 2, 'S': 2, 'B': 2, 'P': 2
    };

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Write SMILES from a Molecule object.
     * @param {Molecule} mol
     * @returns {string}
     */
    function write(mol) {
        if (!mol || mol.atoms.length === 0) return '';

        if (mol.reactionArrow) {
            return writeReaction(mol);
        }

        // Perceive aromaticity
        var aromaticSet = perceiveAromaticity(mol);

        // Multi-component
        var components = mol.getComponents();
        var parts = [];
        for (var i = 0; i < components.length; i++) {
            parts.push(writeComponent(mol, components[i], aromaticSet));
        }
        return parts.join('.');
    }

    /**
     * Write SMILES with the molecule name separated by a tab.
     * Always appends a tab and the name (empty string if no name).
     * @param {Molecule} mol
     * @returns {string}
     */
    function writeNamed(mol) {
        if (!mol || mol.atoms.length === 0) return '\t';

        if (mol.reactionArrow) {
            var rxn = writeReaction(mol);
            return rxn + '\t' + (mol.name || '');
        }

        var aromaticSet = perceiveAromaticity(mol);
        var components = mol.getComponents();
        var parts = [];
        for (var i = 0; i < components.length; i++) {
            parts.push(writeComponent(mol, components[i], aromaticSet));
        }
        var smiles = parts.join('.');
        return smiles + '\t' + (mol.name || '');
    }

    // -----------------------------------------------------------------------
    // Reaction SMILES
    // -----------------------------------------------------------------------

    function writeReaction(mol) {
        var arrow = mol.reactionArrow;
        var midX = (arrow.x1 + arrow.x2) / 2;
        var components = mol.getComponents();
        var reactants = [], products = [];

        // Perceive aromaticity
        var aromaticSet = perceiveAromaticity(mol);

        components.forEach(function(comp) {
            var avgX = 0;
            comp.forEach(function(id) {
                var a = mol.getAtom(id);
                if (a) avgX += a.x;
            });
            avgX /= comp.length;
            if (avgX < midX) {
                reactants.push(writeComponent(mol, comp, aromaticSet));
            } else {
                products.push(writeComponent(mol, comp, aromaticSet));
            }
        });

        return reactants.join('.') + '>>' + products.join('.');
    }

    // -----------------------------------------------------------------------
    // Morgan algorithm for canonical ordering
    // -----------------------------------------------------------------------

    /**
     * Compute Morgan-like canonical ranking for a set of atom IDs.
     * Returns an object mapping atomId -> rank (lower rank = earlier in SMILES).
     */
    function morganRanking(mol, atomIds) {
        var atomSet = {};
        atomIds.forEach(function(id) { atomSet[id] = true; });

        // Initial invariant: degree within this component
        var invariant = {};
        atomIds.forEach(function(id) {
            var neighbors = mol.getNeighbors(id);
            var deg = 0;
            for (var i = 0; i < neighbors.length; i++) {
                if (atomSet[neighbors[i]]) deg++;
            }
            invariant[id] = deg;
        });

        // Iterate: sum neighbor invariants until ranking stabilises
        var prevDistinct = 0;
        for (var iter = 0; iter < 20; iter++) {
            var newInvariant = {};
            atomIds.forEach(function(id) {
                var sum = invariant[id];
                var neighbors = mol.getNeighbors(id);
                for (var i = 0; i < neighbors.length; i++) {
                    if (atomSet[neighbors[i]]) sum += invariant[neighbors[i]];
                }
                newInvariant[id] = sum;
            });
            invariant = newInvariant;

            // Count distinct values
            var vals = {};
            atomIds.forEach(function(id) { vals[invariant[id]] = true; });
            var distinct = Object.keys(vals).length;
            if (distinct === prevDistinct) break;
            prevDistinct = distinct;
        }

        // Build tie-breaking key: [morganInvariant, atomicNumber, charge]
        var sortable = atomIds.map(function(id) {
            var atom = mol.getAtom(id);
            var atomicNum = ATOMIC_NUMBERS[atom.symbol] || 0;
            return {
                id: id,
                key: [invariant[id], atomicNum, atom.charge || 0]
            };
        });

        sortable.sort(function(a, b) {
            for (var i = 0; i < a.key.length; i++) {
                if (a.key[i] !== b.key[i]) return a.key[i] - b.key[i];
            }
            return a.id - b.id; // final tie-break on internal ID
        });

        var rank = {};
        for (var i = 0; i < sortable.length; i++) {
            rank[sortable[i].id] = i;
        }
        return rank;
    }

    // -----------------------------------------------------------------------
    // Write a single connected component
    // -----------------------------------------------------------------------

    function writeComponent(mol, atomIds, aromaticSet) {
        var atomSet = {};
        atomIds.forEach(function(id) { atomSet[id] = true; });

        var rank = morganRanking(mol, atomIds);

        // Start DFS from the atom with the highest canonical rank (best root)
        var startAtom = atomIds[0];
        var bestRank_ = rank[startAtom];
        for (var ri = 1; ri < atomIds.length; ri++) {
            if (rank[atomIds[ri]] > bestRank_) {
                bestRank_ = rank[atomIds[ri]];
                startAtom = atomIds[ri];
            }
        }

        // FIX: pass startAtom so ring closure detection uses the same root as the SMILES DFS traversal
        var rings = findRingClosures(mol, atomIds, rank, startAtom);

        // DFS traversal
        var visited = {};
        var result = '';

        function dfs(atomId, fromBondId) {
            visited[atomId] = true;
            var atom = mol.getAtom(atomId);

            // For stereocentres, we need to determine the correct @/@@ tag
            // based on the actual DFS traversal order of neighbors, not the
            // original parse order. This ensures correct round-trip behaviour.
            var chiralityOverride = '';
            if (atom.chirality) {
                chiralityOverride = _resolveChirality(mol, atom, atomSet, rank,
                    rings, visited, fromBondId);
            }

            result += atomToSmiles(mol, atom, aromaticSet, chiralityOverride);

            // Ring closure digits at this atom
            if (rings.openAt[atomId]) {
                rings.openAt[atomId].forEach(function(rc) {
                    var bondStr = bondTypeToSmiles(rc.bondType, rc.fromAromatic && rc.toAromatic);
                    result += bondStr + ringNumStr(rc.num);
                });
            }
            if (rings.closeAt[atomId]) {
                rings.closeAt[atomId].forEach(function(rc) {
                    var bondStr = bondTypeToSmiles(rc.bondType, rc.fromAromatic && rc.toAromatic);
                    result += bondStr + ringNumStr(rc.num);
                });
            }

            // Gather unvisited neighbors (sorted by canonical rank)
            var bonds = mol.getBondsOfAtom(atomId);
            var branches = [];
            for (var i = 0; i < bonds.length; i++) {
                var bond = bonds[i];
                if (bond.id === fromBondId) continue;
                var neighbor = bond.otherAtom(atomId);
                if (!atomSet[neighbor]) continue;
                if (visited[neighbor]) continue;
                if (rings.isRingBond[bond.id]) continue;
                branches.push({ bond: bond, neighbor: neighbor });
            }

            // Sort branches by canonical rank (highest rank last = main chain)
            branches.sort(function(a, b) {
                return rank[a.neighbor] - rank[b.neighbor];
            });

            if (branches.length === 0) return;

            for (var i = 0; i < branches.length; i++) {
                var br = branches[i];
                var fromAro = aromaticSet[atomId];
                var toAro = aromaticSet[br.neighbor];
                var bondStr = bondTypeToSmiles(br.bond.type, fromAro && toAro);
                var stereoStr = stereoToSmiles(br.bond, atomId);

                if (i < branches.length - 1) {
                    result += '(' + stereoStr + bondStr;
                    dfs(br.neighbor, br.bond.id);
                    result += ')';
                } else {
                    result += stereoStr + bondStr;
                    dfs(br.neighbor, br.bond.id);
                }
            }
        }

        dfs(startAtom, -1);
        return result;
    }

    // -----------------------------------------------------------------------
    // Chirality resolution for SMILES output
    // -----------------------------------------------------------------------

    /**
     * Determine the correct @/@@ for an atom in the DFS traversal order.
     *
     * SMILES chirality convention:
     *   '@' means anticlockwise when looking from the "from" atom (the implicit
     *   H or the atom we came from in the DFS).
     *
     * The original atom.chirality stores the parity relative to the original
     * SMILES parse order. We need to adjust if the DFS traversal visits
     * neighbors in a different order.
     *
     * The approach: build the original neighbor order, build the DFS neighbor
     * order, count the parity of swaps needed to go from one to the other,
     * and flip @/@@  if the parity is odd.
     */
    function _resolveChirality(mol, atom, atomSet, rank, rings, visited, fromBondId) {
        if (!atom.chirality) return '';

        // Get all neighbor IDs in their natural (adjacency list) order
        var bonds = mol.getBondsOfAtom(atom.id);
        var originalOrder = [];
        for (var i = 0; i < bonds.length; i++) {
            var nid = bonds[i].otherAtom(atom.id);
            if (atomSet[nid]) originalOrder.push(nid);
        }

        // Build the DFS output order: the order in which neighbors will appear
        // in the SMILES string for this atom.
        // 1. Ring closure atoms (already visited, in ring-closure order)
        // 2. Branch atoms (unvisited, sorted by canonical rank)
        var dfsOrder = [];

        // The "from" atom (parent in DFS) is implicitly first in SMILES convention.
        // If there's an implicit H, it takes the first position instead.
        var fromAtomId = -1;
        if (fromBondId >= 0) {
            var fromBond = mol.getBond(fromBondId);
            if (fromBond) fromAtomId = fromBond.otherAtom(atom.id);
        }

        // Collect ring-closure neighbors at this atom (they come before branches)
        var ringNeighbors = [];
        if (rings.openAt[atom.id]) {
            for (var r = 0; r < rings.openAt[atom.id].length; r++) {
                // Find the atom at the other end of this ring closure
                var rc = rings.openAt[atom.id][r];
                // Ring closures opened at this atom: the "close" end is another atom
                // We need to find which atom closes this ring number
                for (var aId in rings.closeAt) {
                    if (!rings.closeAt.hasOwnProperty(aId)) continue;
                    var closes = rings.closeAt[aId];
                    for (var ci = 0; ci < closes.length; ci++) {
                        if (closes[ci].num === rc.num) {
                            ringNeighbors.push(parseInt(aId));
                        }
                    }
                }
            }
        }
        if (rings.closeAt[atom.id]) {
            for (var r = 0; r < rings.closeAt[atom.id].length; r++) {
                var rc = rings.closeAt[atom.id][r];
                for (var aId in rings.openAt) {
                    if (!rings.openAt.hasOwnProperty(aId)) continue;
                    var opens = rings.openAt[aId];
                    for (var oi = 0; oi < opens.length; oi++) {
                        if (opens[oi].num === rc.num) {
                            ringNeighbors.push(parseInt(aId));
                        }
                    }
                }
            }
        }

        // Add ring closure neighbors
        for (var r = 0; r < ringNeighbors.length; r++) {
            dfsOrder.push(ringNeighbors[r]);
        }

        // Collect unvisited branch neighbors sorted by rank
        var branches = [];
        for (var i = 0; i < bonds.length; i++) {
            var bond = bonds[i];
            if (bond.id === fromBondId) continue;
            var neighbor = bond.otherAtom(atom.id);
            if (!atomSet[neighbor]) continue;
            if (visited[neighbor]) continue;
            if (rings.isRingBond[bond.id]) continue;
            branches.push(neighbor);
        }
        branches.sort(function(a, b) { return rank[a] - rank[b]; });
        for (var i = 0; i < branches.length; i++) {
            dfsOrder.push(branches[i]);
        }

        // Build the permutation from original order to DFS order.
        // If from-atom is present, it comes first in both.
        var origSeq = [];
        var dfsSeq = [];

        if (fromAtomId >= 0) {
            origSeq.push(fromAtomId);
            dfsSeq.push(fromAtomId);
        }
        for (var i = 0; i < originalOrder.length; i++) {
            if (originalOrder[i] !== fromAtomId) origSeq.push(originalOrder[i]);
        }
        for (var i = 0; i < dfsOrder.length; i++) {
            if (dfsOrder[i] !== fromAtomId) dfsSeq.push(dfsOrder[i]);
        }

        // Add implicit hydrogens to both sequences
        var hCount = mol.calcHydrogens(atom.id);
        for (var h = 0; h < hCount; h++) {
            origSeq.push(-(h + 1));
            dfsSeq.push(-(h + 1));
        }

        // Compute parity of the permutation from origSeq to dfsSeq
        if (origSeq.length !== dfsSeq.length || origSeq.length < 3) {
            return atom.chirality;
        }

        // Map origSeq indices to dfsSeq indices
        var perm = [];
        for (var i = 0; i < origSeq.length; i++) {
            var idx = -1;
            for (var j = 0; j < dfsSeq.length; j++) {
                if (dfsSeq[j] === origSeq[i]) { idx = j; break; }
            }
            if (idx < 0) return atom.chirality; // safety fallback
            perm.push(idx);
        }

        var parity = _permParity(perm);
        if (parity === 1) {
            // Odd permutation: flip chirality
            return atom.chirality === '@' ? '@@' : '@';
        }
        return atom.chirality;
    }

    /**
     * Compute permutation parity (0 = even, 1 = odd).
     */
    function _permParity(perm) {
        var n = perm.length;
        var arr = perm.slice();
        var swaps = 0;
        for (var i = 0; i < n; i++) {
            while (arr[i] !== i) {
                var target = arr[i];
                arr[i] = arr[target];
                arr[target] = target;
                swaps++;
            }
        }
        return swaps % 2;
    }

    // -----------------------------------------------------------------------
    // Ring closure detection (DFS back-edges with canonical ordering)
    // -----------------------------------------------------------------------

    function findRingClosures(mol, atomIds, rank, startAtomId) {
        var atomSet = {};
        atomIds.forEach(function(id) { atomSet[id] = true; });

        var visited = {};
        var openAt  = {};
        var closeAt = {};
        var isRingBond = {};
        var nextNum = 1;

        function dfs(atomId, parentBondId) {
            visited[atomId] = true;
            var bonds = mol.getBondsOfAtom(atomId);

            // Sort bonds by neighbor rank for deterministic ring numbering
            var sortedBonds = bonds.slice().sort(function(a, b) {
                var na = a.otherAtom(atomId);
                var nb = b.otherAtom(atomId);
                return (rank[na] || 0) - (rank[nb] || 0);
            });

            for (var i = 0; i < sortedBonds.length; i++) {
                var bond = sortedBonds[i];
                if (bond.id === parentBondId) continue;
                var neighbor = bond.otherAtom(atomId);
                if (!atomSet[neighbor]) continue;
                if (visited[neighbor]) {
                    if (!isRingBond[bond.id]) {
                        var num = nextNum++;
                        var atomA = mol.getAtom(neighbor);
                        var atomB = mol.getAtom(atomId);
                        if (!openAt[neighbor])  openAt[neighbor] = [];
                        openAt[neighbor].push({
                            num: num, bondType: bond.type,
                            fromAromatic: !!(atomA && atomA.aromatic),
                            toAromatic: !!(atomB && atomB.aromatic)
                        });
                        if (!closeAt[atomId]) closeAt[atomId] = [];
                        closeAt[atomId].push({
                            num: num, bondType: bond.type,
                            fromAromatic: !!(atomA && atomA.aromatic),
                            toAromatic: !!(atomB && atomB.aromatic)
                        });
                        isRingBond[bond.id] = true;
                    }
                } else {
                    dfs(neighbor, bond.id);
                }
            }
        }

        // FIX: use startAtomId (same root as SMILES DFS) instead of atomIds[0]
        dfs(startAtomId || atomIds[0], -1);
        return { openAt: openAt, closeAt: closeAt, isRingBond: isRingBond };
    }

    // -----------------------------------------------------------------------
    // Aromatic perception — Hueckel 4n+2 rule on SSSR rings
    // -----------------------------------------------------------------------

    /**
     * Returns a set (object) of atom IDs that are aromatic.
     * Uses the molecule's findRings to get SSSR, then checks each ring
     * for planarity candidates (C, N, O, S) and pi-electron count.
     */
    function perceiveAromaticity(mol) {
        var aromaticAtoms = {};

        // Atoms that can participate in aromatic rings
        var aromaticElements = { 'C': true, 'N': true, 'O': true, 'S': true };

        var rings = mol.findRings(8); // find rings up to size 8

        for (var r = 0; r < rings.length; r++) {
            var ring = rings[r].atoms;
            if (ring.length < 3) continue;

            // Check: all atoms must be potential aromatic elements
            var allAromatic = true;
            for (var i = 0; i < ring.length; i++) {
                var atom = mol.getAtom(ring[i]);
                if (!atom || !aromaticElements[atom.symbol]) {
                    allAromatic = false;
                    break;
                }
            }
            if (!allAromatic) continue;

            // Check: ring must have alternating single/double bonds, or
            // count pi electrons using Hueckel rule
            var piElectrons = 0;
            var valid = true;
            for (var i = 0; i < ring.length; i++) {
                var atomId = ring[i];
                var atom = mol.getAtom(atomId);
                var sym = atom.symbol;

                // Check ALL ring bonds for this atom for double bonds, not just
                // the bond to (i+1).  This avoids depending on ring traversal
                // order when Kekule bonds are present.
                var prevId = ring[(i - 1 + ring.length) % ring.length];
                var nextId = ring[(i + 1) % ring.length];
                var bondPrev = mol.getBondBetween(atomId, prevId);
                var bondNext = mol.getBondBetween(atomId, nextId);
                var hasDoubleBondInRing = (bondPrev && bondPrev.type === Molecule.BOND_DOUBLE) ||
                                          (bondNext && bondNext.type === Molecule.BOND_DOUBLE);

                if (sym === 'C') {
                    if (hasDoubleBondInRing) {
                        piElectrons += 1; // one electron from double bond
                    } else {
                        // Check if this C has an exocyclic double bond (no pi contribution)
                        var exoBonds = mol.getBondsOfAtom(atomId);
                        var hasExoDouble = false;
                        for (var j = 0; j < exoBonds.length; j++) {
                            var other = exoBonds[j].otherAtom(atomId);
                            if (ring.indexOf(other) < 0 && exoBonds[j].type === Molecule.BOND_DOUBLE) {
                                hasExoDouble = true;
                                break;
                            }
                        }
                        if (hasExoDouble) {
                            valid = false; break; // not aromatic
                        }
                        // Charged carbon with no in-ring double: lone pair / empty orbital.
                        // Cp- anion: [cH-] donates 2 electrons; tropylium [c+] donates 0.
                        var cCharge = atom.charge || 0;
                        if (cCharge === -1) {
                            piElectrons += 2; // anionic C: lone pair into pi system
                        } else if (cCharge === 1) {
                            piElectrons += 0; // cationic C: empty p-orbital
                        } else if (atom.aromatic) {
                            piElectrons += 1; // aromatic-flagged C (lowercase 'c')
                        } else {
                            // Neutral sp3 C with no in-ring double and no exocyclic
                            // double: not part of a pi system (e.g. cyclohexane,
                            // saturated steroid rings, cyclopentadiene CH2).
                            valid = false; break;
                        }
                    }
                } else if (sym === 'N') {
                    if (hasDoubleBondInRing) {
                        // Pyridine-type (Kekule): imine nitrogen contributes 1 pi electron
                        piElectrons += 1;
                    } else {
                        // Distinguish pyrrole-N (lone pair, 2 e-) from pyridine-N
                        // (1 e-) when bonds are stored as single (aromatic SMILES
                        // input or all-single-bond ring).
                        //
                        // Pyrrole-type N: has an explicit H (e.g. [nH]) or is
                        //   bonded to only 2 ring atoms with an explicit H count.
                        //   Lone pair donates 2 electrons to the pi system.
                        // Pyridine-type N: 2 ring bonds, no H, one lone pair not
                        //   in the pi system; contributes 1 pi electron from the
                        //   double bond it would have in Kekule form.
                        //
                        // Detection: if the atom was parsed as aromatic (from
                        // lowercase SMILES), an explicit H > 0 means pyrrole-type.
                        // If hydrogens == 0 (bracket [n] or organic 'n' with 2
                        // bonds), it is pyridine-type.  For non-aromatic input
                        // with all single bonds, count total degree + H vs
                        // standard N valence of 3.
                        var explH = atom.hydrogens; // >= 0 for bracket atoms, -1 for organic
                        var nRingBonds = 2; // always 2 for a ring atom
                        if (atom.aromatic) {
                            // Aromatic SMILES: 'n' means pyridine (no H),
                            // '[nH]' means pyrrole (has H).
                            // Pyridinium '[nH+]' has explH>0 but charge=+1 — its
                            // lone pair is consumed by H+; it is pyridine-type
                            // (1 pi e-), NOT pyrrole-type.
                            var nCharge = atom.charge || 0;
                            if (explH > 0 && nCharge <= 0) {
                                piElectrons += 2; // pyrrole-type
                            } else {
                                piElectrons += 1; // pyridine-type / pyridinium
                            }
                        } else {
                            // Non-aromatic input (all single bonds in ring):
                            // check total bond order + explicit/default H.
                            // If totalValence == 3 and no H -> pyridine-type
                            var bondSum = mol.bondOrderSum(atomId);
                            var totalNeighbors = mol.getNeighbors(atomId).length;
                            if (explH > 0 || (explH < 0 && totalNeighbors === 2 && bondSum === 2)) {
                                // 2 single bonds + implicit H -> pyrrole-type
                                piElectrons += 2;
                            } else {
                                piElectrons += 1;
                            }
                        }
                    }
                } else if (sym === 'O' || sym === 'S') {
                    // Furan/thiophene-type: lone pair contributes 2 electrons
                    piElectrons += 2;
                }
            }

            if (!valid) continue;

            // Hueckel rule: 4n+2 for n = 0, 1, 2, 3, ...
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

    // -----------------------------------------------------------------------
    // Atom to SMILES string
    // -----------------------------------------------------------------------

    function atomToSmiles(mol, atom, aromaticSet, chiralityOverride) {
        var sym = atom.symbol;
        var charge  = atom.charge || 0;
        var isotope = atom.isotope || 0;
        var mapNum  = atom.mapNumber || 0;
        var chirality = (chiralityOverride !== undefined && chiralityOverride !== '')
            ? chiralityOverride : (atom.chirality || '');
        var isAromatic = !!(aromaticSet && aromaticSet[atom.id]);

        // Compute the default implicit H count that a reader would calculate
        // for this atom (used to decide whether brackets are needed).
        // For aromatic atoms, use the aromatic valence to match what the
        // reader's calcHydrogens would produce.
        var bondSum = mol.bondOrderSum(atom.id);
        var defaultValence;
        if (isAromatic && AROMATIC_WRITER_VALENCE[sym] !== undefined) {
            // Check if bonds are stored as single (aromatic input)
            var bonds = mol.getBondsOfAtom(atom.id);
            var hasMultipleBond = false;
            for (var bi = 0; bi < bonds.length; bi++) {
                if (bonds[bi].type > Molecule.BOND_SINGLE) { hasMultipleBond = true; break; }
            }
            defaultValence = hasMultipleBond
                ? getDefaultValence(sym, bondSum, charge)
                : AROMATIC_WRITER_VALENCE[sym];
        } else {
            defaultValence = getDefaultValence(sym, bondSum, charge);
        }
        var expectedH = Math.max(0, defaultValence - bondSum - Math.abs(charge));
        var explicitH = (atom.hydrogens >= 0) ? atom.hydrogens : -1;

        // Determine if brackets are needed
        var needBrackets = false;
        if (!ORGANIC[sym] && !(isAromatic && AROMATIC_ORGANIC[sym])) needBrackets = true;
        if (charge !== 0)    needBrackets = true;
        if (isotope > 0)     needBrackets = true;
        if (mapNum > 0)      needBrackets = true;
        if (chirality)       needBrackets = true;
        if (explicitH >= 0 && explicitH !== expectedH) needBrackets = true;

        if (!needBrackets) {
            if (isAromatic) return sym.toLowerCase();
            return sym;
        }

        // Build bracketed representation
        var s = '[';
        if (isotope > 0) s += isotope;
        if (isAromatic) {
            s += sym.toLowerCase();
        } else {
            s += sym;
        }
        if (chirality) s += chirality;

        // H count inside brackets
        var hCount = (explicitH >= 0) ? explicitH : expectedH;
        if (hCount > 0) {
            s += 'H';
            if (hCount > 1) s += hCount;
        }

        // Charge
        if (charge > 0) {
            s += '+';
            if (charge > 1) s += charge;
        } else if (charge < 0) {
            s += '-';
            if (charge < -1) s += Math.abs(charge);
        }

        // Atom map
        if (mapNum > 0) s += ':' + mapNum;

        s += ']';
        return s;
    }

    // -----------------------------------------------------------------------
    // Bond type to SMILES
    // -----------------------------------------------------------------------

    function bondTypeToSmiles(type, bothAromatic) {
        if (type === Molecule.BOND_DOUBLE && !bothAromatic) return '=';
        if (type === Molecule.BOND_TRIPLE) return '#';
        // Single bond between aromatic atoms is implicit; between non-aromatic also implicit
        return '';
    }

    /**
     * Stereo bond markers (/ and \) for E/Z double bonds.
     */
    function stereoToSmiles(bond, fromAtomId) {
        if (!bond.stereo || bond.stereo === Molecule.STEREO_NONE) return '';
        // stereo 1 = '/', stereo 6 = '\'
        if (bond.stereo === 1) {
            return (bond.atom1 === fromAtomId) ? '/' : '\\';
        }
        if (bond.stereo === 6) {
            return (bond.atom1 === fromAtomId) ? '\\' : '/';
        }
        return '';
    }

    // -----------------------------------------------------------------------
    // Ring number formatting
    // -----------------------------------------------------------------------

    function ringNumStr(num) {
        if (num < 10) return '' + num;
        return '%' + num;
    }

    // -----------------------------------------------------------------------
    // Valence helper
    // -----------------------------------------------------------------------

    function getDefaultValence(symbol, bondSum, charge) {
        var vList = STANDARD_VALENCES[symbol];
        if (!vList) return 0;
        var target = bondSum + Math.abs(charge || 0);
        for (var i = 0; i < vList.length; i++) {
            if (vList[i] >= target) return vList[i];
        }
        return vList[vList.length - 1];
    }

    // -----------------------------------------------------------------------
    // Export
    // -----------------------------------------------------------------------

    global.SmilesWriter = { write: write, writeNamed: writeNamed };

})(typeof window !== 'undefined' ? window : this);
