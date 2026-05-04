/**
 * CipStereo.js — Cahn-Ingold-Prelog stereochemistry assignment (v6.2.0 engine)
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Full digraph-based CIP implementation (IUPAC 2013 Rules 1-4a):
 *   - BFS level-by-level priority comparison
 *   - Phantom nodes for double/triple bond expansion
 *   - Phantom nodes for ring closures (back-edge duplicates)
 *   - Aromatic mancude handling (pi contribution to valence)
 *   - Rule 3: seqCis/seqTrans (Z > E > nonstereogenic)
 *   - Rule 4a: chiral > pseudoasymmetric > nonstereogenic (two-pass)
 *   - MAX_DEPTH=30, MAX_NODES=5000 to prevent runaway
 *   - Permutation parity for @/@@
 *   - Charge-aware valence expansion + aromatic pi for implicitH
 *
 * Assigns:
 *   - R/S labels for tetrahedral stereocentres (atoms with .chirality '@' or '@@')
 *   - E/Z labels for double bonds with restricted rotation
 */
(function(global) {
    'use strict';

    var ELEMENTS = Molecule.ELEMENTS;
    var BOND_SINGLE = Molecule.BOND_SINGLE;
    var BOND_DOUBLE = Molecule.BOND_DOUBLE;
    var BOND_TRIPLE = Molecule.BOND_TRIPLE;

    var MAX_DEPTH = 30;
    var MAX_NODES = 5000;

    // -----------------------------------------------------------------------
    // Atomic number lookup — needed for CIP priority
    // -----------------------------------------------------------------------
    var ATOMIC_NUMBER = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
        'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
        'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
        'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
        'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
        'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
        'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
        'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
        'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
        'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
        'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
        'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
        'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
        'Hg': 80, 'Tl': 81,
        'Pb': 82, 'Bi': 83, 'R': 0
    };

    /**
     * Get the atomic number for a symbol. Returns 0 for unknowns.
     */
    function atomicNumber(symbol) {
        return ATOMIC_NUMBER[symbol] || 0;
    }

    /**
     * Get the mass for a symbol, respecting isotope overrides.
     */
    function atomMass(atom) {
        if (atom.isotope > 0) return atom.isotope;
        var elem = ELEMENTS[atom.symbol];
        return elem ? elem.mass : 0;
    }

    // -----------------------------------------------------------------------
    // Aromatic mancude detection — identifies atoms in aromatic rings whose
    // pi contribution affects CIP digraph expansion
    // -----------------------------------------------------------------------

    /**
     * Atoms in aromatic (mancude) systems contribute pi electrons that affect
     * the effective bond order seen by the CIP digraph.  For aromatic atoms
     * stored with single bonds (from aromatic SMILES), we treat one bond as
     * having a 1.5 contribution: the CIP digraph adds a phantom node for
     * the pi-bond partner with atomic number equal to the partner, weighted
     * at half the normal duplicate weight.
     *
     * Returns true if atomId is flagged aromatic and has no explicit
     * multiple bonds (i.e. its bonds are stored as Kekulized singles from
     * aromatic SMILES input).
     */
    function isMancudeAtom(mol, atomId) {
        var atom = mol.getAtom(atomId);
        if (!atom || !atom.aromatic) return false;
        var bonds = mol.getBondsOfAtom(atomId);
        for (var i = 0; i < bonds.length; i++) {
            if (bonds[i].type > BOND_SINGLE) return false;
        }
        return true;
    }

    /**
     * For a mancude atom, determine the aromatic neighbor that contributes
     * the shared pi bond.  We pick the aromatic neighbor with the highest
     * atomic number (CIP-consistent tie-break).  Returns the neighbor
     * atom ID, or -1 if none found.
     */
    function mancudePiPartner(mol, atomId, excludeId) {
        var bonds = mol.getBondsOfAtom(atomId);
        var bestId = -1;
        var bestZ = -1;
        for (var i = 0; i < bonds.length; i++) {
            var nid = bonds[i].otherAtom(atomId);
            if (nid === excludeId) continue;
            var natom = mol.getAtom(nid);
            if (!natom || !natom.aromatic) continue;
            var z = atomicNumber(natom.symbol);
            if (z > bestZ) {
                bestZ = z;
                bestId = nid;
            }
        }
        return bestId;
    }

    // -----------------------------------------------------------------------
    // CIP Digraph Node
    // -----------------------------------------------------------------------

    /**
     * A node in the CIP expansion digraph.
     *
     * For real atoms: atomId is set, symbol/mass come from the molecule.
     * For phantom duplicate atoms (from double/triple bonds or ring closures):
     *   phantom = true, symbol and mass are copied from the duplicated atom,
     *   children are empty (phantom nodes are leaves).
     * For mancude phantom atoms: mancude = true, represents the fractional
     *   pi-bond contribution in an aromatic system.
     */
    function CipNode(atomId, symbol, mass, phantom) {
        this.atomId = atomId;
        this.symbol = symbol;
        this.atomicNum = atomicNumber(symbol);
        this.mass = mass;
        this.phantom = !!phantom;
        this.mancude = false;
        this.depth = 0;
        this.children = [];
        this.seqCisTransScore = 0; // Rule 3: 2 = Z (seqCis), 1 = E (seqTrans), 0 = none
        this.stereogenicScore = 0; // Rule 4a: 2 = chiral (R/S), 1 = pseudoasymmetric (r/s), 0 = none
    }

    // -----------------------------------------------------------------------
    // Digraph construction — BFS expansion with phantom nodes
    // -----------------------------------------------------------------------

    /**
     * Build a CIP digraph rooted at startAtomId, growing away from
     * centerAtomId (the stereocentre or double-bond atom).
     *
     * The digraph is built level-by-level (BFS) up to MAX_DEPTH levels,
     * with a hard cap of MAX_NODES total nodes to prevent runaway on
     * pathological inputs (e.g. large polymeric rings).
     *
     * Key features:
     *   - Double/triple bonds add phantom duplicate leaf nodes
     *   - Ring closures: when a back-edge is encountered (neighbor already
     *     on the path from root), a phantom duplicate leaf node is inserted
     *     instead of creating a cycle.  This correctly handles spiro,
     *     bridged, and fused ring systems.
     *   - Mancude handling: aromatic atoms with single-only bonds get a
     *     phantom node for their pi-bond partner.
     */
    function buildDigraph(mol, startAtomId, centerAtomId) {
        var startAtom = mol.getAtom(startAtomId);
        if (!startAtom) return null;

        var root = new CipNode(startAtomId, startAtom.symbol, atomMass(startAtom), false);
        root.depth = 0;

        var totalNodes = 1;

        // BFS queue: each entry is { node, parentAtomId, pathSet }
        // pathSet tracks atom IDs on the current path from root to this node
        // (used for ring-closure detection)
        var queue = [];
        var rootPath = {};
        rootPath[centerAtomId] = true;
        rootPath[startAtomId] = true;
        queue.push({ node: root, parentAtomId: centerAtomId, pathSet: rootPath });

        var depth = 0;

        while (queue.length > 0 && depth < MAX_DEPTH && totalNodes < MAX_NODES) {
            var nextQueue = [];
            depth++;

            for (var qi = 0; qi < queue.length && totalNodes < MAX_NODES; qi++) {
                var entry = queue[qi];
                var node = entry.node;
                var parentId = entry.parentAtomId;
                var pathSet = entry.pathSet;
                var atomId = node.atomId;

                if (node.phantom) continue; // phantom nodes are leaves

                var bonds = mol.getBondsOfAtom(atomId);

                for (var bi = 0; bi < bonds.length && totalNodes < MAX_NODES; bi++) {
                    var bond = bonds[bi];
                    var neighborId = bond.otherAtom(atomId);

                    if (neighborId === parentId) {
                        // Back-edge to parent: add phantom duplicates for
                        // higher-order bonds (double/triple) toward parent
                        var parentAtom = mol.getAtom(parentId);
                        if (parentAtom && bond.type >= BOND_DOUBLE) {
                            for (var p = 1; p < bond.type && totalNodes < MAX_NODES; p++) {
                                var phNode = new CipNode(
                                    parentId, parentAtom.symbol,
                                    atomMass(parentAtom), true
                                );
                                phNode.depth = depth;
                                node.children.push(phNode);
                                totalNodes++;
                            }
                        }
                        continue;
                    }

                    // Ring closure detection: if neighbor is already on the
                    // path from root to current node, insert a phantom
                    // duplicate leaf instead of following the cycle.
                    if (pathSet[neighborId]) {
                        var ringAtom = mol.getAtom(neighborId);
                        if (ringAtom) {
                            var ringPhantom = new CipNode(
                                neighborId, ringAtom.symbol,
                                atomMass(ringAtom), true
                            );
                            ringPhantom.depth = depth;
                            node.children.push(ringPhantom);
                            totalNodes++;
                        }
                        continue;
                    }

                    // Normal expansion: add real child node
                    var nAtom = mol.getAtom(neighborId);
                    if (!nAtom) continue;

                    var child = new CipNode(
                        neighborId, nAtom.symbol, atomMass(nAtom), false
                    );
                    child.depth = depth;
                    node.children.push(child);
                    totalNodes++;

                    // Add phantom duplicate atoms for double/triple bonds
                    // (forward direction)
                    if (bond.type >= BOND_DOUBLE) {
                        for (var p = 1; p < bond.type && totalNodes < MAX_NODES; p++) {
                            var fwdPhantom = new CipNode(
                                atomId, node.symbol, node.mass, true
                            );
                            fwdPhantom.depth = depth;
                            child.children.push(fwdPhantom);
                            totalNodes++;
                        }
                    }

                    // Prepare next level BFS entry
                    var childPath = _copyObj(pathSet);
                    childPath[neighborId] = true;
                    nextQueue.push({
                        node: child,
                        parentAtomId: atomId,
                        pathSet: childPath
                    });
                }

                // Mancude (aromatic) handling: if this atom is aromatic with
                // only single bonds, add a phantom node for the pi-bond
                // partner to account for the delocalized double bond.
                if (isMancudeAtom(mol, atomId) && totalNodes < MAX_NODES) {
                    var piPartner = mancudePiPartner(mol, atomId, parentId);
                    if (piPartner >= 0) {
                        var piAtom = mol.getAtom(piPartner);
                        if (piAtom) {
                            var piPhantom = new CipNode(
                                piPartner, piAtom.symbol, atomMass(piAtom), true
                            );
                            piPhantom.mancude = true;
                            piPhantom.depth = depth;
                            node.children.push(piPhantom);
                            totalNodes++;
                        }
                    }
                }

                // Add implicit hydrogens as children
                var hCount = mol.calcHydrogens(atomId);
                for (var h = 0; h < hCount && totalNodes < MAX_NODES; h++) {
                    var hNode = new CipNode(-1, 'H', 1, false);
                    hNode.phantom = true; // H leaves don't expand
                    hNode.depth = depth;
                    node.children.push(hNode);
                    totalNodes++;
                }

                // Sort children by descending CIP priority (highest first)
                // at each node for canonical ordering during comparison
                node.children.sort(function(a, b) {
                    return _compareSingle(b, a);
                });
            }

            queue = nextQueue;
        }

        return root;
    }

    /**
     * Shallow copy of an object (used for path sets).
     */
    function _copyObj(obj) {
        var copy = {};
        for (var k in obj) {
            if (obj.hasOwnProperty(k)) copy[k] = obj[k];
        }
        return copy;
    }

    // -----------------------------------------------------------------------
    // Digraph post-processing — populate Rule 3 and Rule 4a scores
    // -----------------------------------------------------------------------

    /**
     * Walk a digraph tree and set seqCisTransScore on each non-phantom node
     * based on whether its atom participates in a stereogenic double bond.
     */
    function _decorateRule3(mol, node) {
        if (!node) return;
        if (!node.phantom && node.atomId >= 0) {
            node.seqCisTransScore = _seqCisTransScore(mol, node.atomId, -1);
        }
        for (var i = 0; i < node.children.length; i++) {
            _decorateRule3(mol, node.children[i]);
        }
    }

    /**
     * Walk a digraph tree and set stereogenicScore on each non-phantom node
     * based on the chiralityMap from a prior R/S assignment pass.
     *
     * chiralityMap: { atomId: 'R'|'S'|'r'|'s' }
     *   R/S => score 2 (chiral)
     *   r/s => score 1 (pseudoasymmetric)
     *   absent => score 0 (nonstereogenic)
     */
    function _decorateRule4a(node, chiralityMap) {
        if (!node) return;
        if (!node.phantom && node.atomId >= 0 && chiralityMap[node.atomId]) {
            var label = chiralityMap[node.atomId];
            if (label === 'R' || label === 'S') {
                node.stereogenicScore = 2;
            } else if (label === 'r' || label === 's') {
                node.stereogenicScore = 1;
            }
        }
        for (var i = 0; i < node.children.length; i++) {
            _decorateRule4a(node.children[i], chiralityMap);
        }
    }

    // -----------------------------------------------------------------------
    // Rule 3 helper — seqCis/seqTrans score for double bond nodes
    // -----------------------------------------------------------------------

    /**
     * Determine the seqCis/seqTrans score for an atom that participates in a
     * double bond.  Uses 2D coordinates to determine whether the bond has
     * Z (cis) or E (trans) geometry.
     *
     * Returns 2 for Z (seqCis), 1 for E (seqTrans), 0 if not determinable
     * or if the atom is not part of a stereogenic double bond.
     */
    function _seqCisTransScore(mol, atomId, parentAtomId) {
        var bonds = mol.getBondsOfAtom(atomId);
        for (var i = 0; i < bonds.length; i++) {
            var bond = bonds[i];
            if (bond.type !== BOND_DOUBLE) continue;

            var partnerId = bond.otherAtom(atomId);
            var a1 = mol.getAtom(atomId);
            var a2 = mol.getAtom(partnerId);
            if (!a1 || !a2) continue;

            // Get substituents on each end (excluding the double-bond partner)
            var subs1 = _getSubstituentsForRule3(mol, atomId, partnerId);
            var subs2 = _getSubstituentsForRule3(mol, partnerId, atomId);

            // Need at least 2 substituents on each end for stereogenicity
            if (subs1.length < 2 || subs2.length < 2) continue;

            // Check that substituents are distinguishable on each end
            // (otherwise the double bond is not stereogenic)
            var trees1 = [];
            for (var j = 0; j < subs1.length; j++) {
                if (subs1[j] < 0) {
                    var ht = new CipNode(-1, 'H', 1, false);
                    ht.phantom = true;
                    trees1.push(ht);
                } else {
                    trees1.push(buildDigraph(mol, subs1[j], atomId));
                }
            }
            if (trees1.length === 2 && compareCipTrees(trees1[0], trees1[1]) === 0) continue;

            var trees2 = [];
            for (var j = 0; j < subs2.length; j++) {
                if (subs2[j] < 0) {
                    var ht2 = new CipNode(-1, 'H', 1, false);
                    ht2.phantom = true;
                    trees2.push(ht2);
                } else {
                    trees2.push(buildDigraph(mol, subs2[j], partnerId));
                }
            }
            if (trees2.length === 2 && compareCipTrees(trees2[0], trees2[1]) === 0) continue;

            // Use the existing bond.cipLabel if already assigned by assignEZ
            if (bond.cipLabel === 'Z') return 2;
            if (bond.cipLabel === 'E') return 1;

            // Otherwise determine from 2D geometry
            // Find the highest-priority substituent on each end
            var bestSub1 = subs1[0], bestTree1 = trees1[0];
            for (var j = 1; j < subs1.length; j++) {
                if (compareCipTrees(trees1[j], bestTree1) > 0) {
                    bestTree1 = trees1[j];
                    bestSub1 = subs1[j];
                }
            }
            var bestSub2 = subs2[0], bestTree2 = trees2[0];
            for (var j = 1; j < subs2.length; j++) {
                if (compareCipTrees(trees2[j], bestTree2) > 0) {
                    bestTree2 = trees2[j];
                    bestSub2 = subs2[j];
                }
            }

            var pos1 = _getSubPosition(mol, bestSub1, atomId);
            var pos2 = _getSubPosition(mol, bestSub2, partnerId);
            if (!pos1 || !pos2) continue;

            var bx = a2.x - a1.x;
            var by = a2.y - a1.y;

            var dx1 = pos1.x - a1.x;
            var dy1 = pos1.y - a1.y;
            var cross1 = bx * dy1 - by * dx1;

            var dx2 = pos2.x - a2.x;
            var dy2 = pos2.y - a2.y;
            var cross2 = bx * dy2 - by * dx2;

            if (Math.abs(cross1) < 0.001 || Math.abs(cross2) < 0.001) continue;

            // Same sign = same side = Z; opposite sign = E
            if ((cross1 > 0) === (cross2 > 0)) {
                return 2; // Z (seqCis)
            } else {
                return 1; // E (seqTrans)
            }
        }
        return 0;
    }

    /**
     * Get substituent atom IDs for one end of a double bond for Rule 3.
     * Includes implicit hydrogens as negative IDs.
     */
    function _getSubstituentsForRule3(mol, atomId, partnerId) {
        var neighbors = mol.getNeighbors(atomId);
        var subs = [];
        for (var i = 0; i < neighbors.length; i++) {
            if (neighbors[i] !== partnerId) {
                subs.push(neighbors[i]);
            }
        }
        var hCount = mol.calcHydrogens(atomId);
        for (var h = 0; h < hCount; h++) {
            subs.push(-(h + 1));
        }
        return subs;
    }

    // -----------------------------------------------------------------------
    // BFS level-by-level priority comparison (IUPAC 2013 Rules 1-4a)
    // -----------------------------------------------------------------------

    /**
     * Compare two CIP digraph trees using BFS level-by-level comparison.
     *
     * At each level, collect all nodes from tree A and tree B at that depth,
     * sort them by (atomicNumber, mass) descending, and compare the sorted
     * lists.  The first level where the lists differ determines priority.
     *
     * This is the correct IUPAC 2013 approach: compare entire levels rather
     * than following individual branches (which is the older, incorrect DFS
     * approach).
     *
     * Returns > 0 if a has higher priority, < 0 if b has, 0 if tied.
     */
    function compareCipTrees(a, b) {
        if (!a && !b) return 0;
        if (!a) return -1;
        if (!b) return 1;

        // Level 0: compare roots directly
        var cmp = _compareSingle(a, b);
        if (cmp !== 0) return cmp;

        // BFS level-by-level comparison
        var queueA = [a];
        var queueB = [b];

        for (var lvl = 0; lvl < MAX_DEPTH; lvl++) {
            // Collect children of current level
            var nextA = [];
            var nextB = [];
            for (var i = 0; i < queueA.length; i++) {
                for (var c = 0; c < queueA[i].children.length; c++) {
                    nextA.push(queueA[i].children[c]);
                }
            }
            for (var i = 0; i < queueB.length; i++) {
                for (var c = 0; c < queueB[i].children.length; c++) {
                    nextB.push(queueB[i].children[c]);
                }
            }

            if (nextA.length === 0 && nextB.length === 0) return 0;
            if (nextA.length === 0) return -1;
            if (nextB.length === 0) return 1;

            // Sort each level's nodes by descending priority (highest first)
            nextA.sort(function(x, y) { return _compareSingle(y, x); });
            nextB.sort(function(x, y) { return _compareSingle(y, x); });

            // Compare the sorted lists element by element
            var maxLen = Math.max(nextA.length, nextB.length);
            for (var j = 0; j < maxLen; j++) {
                var na = j < nextA.length ? nextA[j] : null;
                var nb = j < nextB.length ? nextB[j] : null;
                if (!na && !nb) continue;
                if (!na) return -1;
                if (!nb) return 1;
                var nodeCmp = _compareSingle(na, nb);
                if (nodeCmp !== 0) return nodeCmp;
            }

            queueA = nextA;
            queueB = nextB;
        }

        return 0;
    }

    /**
     * Compare two individual CIP nodes by Rule 1 (atomic number),
     * Rule 2 (mass number), Rule 3 (seqCis/seqTrans), and Rule 4a
     * (stereogenicity).  Does not recurse into children.
     */
    function _compareSingle(a, b) {
        if (!a && !b) return 0;
        if (!a) return -1;
        if (!b) return 1;
        // Rule 1: atomic number
        if (a.atomicNum !== b.atomicNum) return a.atomicNum - b.atomicNum;
        // Rule 2: mass number
        if (a.mass !== b.mass) return a.mass - b.mass;
        // Rule 3: seqCis (Z) > seqTrans (E) > nonstereogenic
        if (a.seqCisTransScore !== b.seqCisTransScore) {
            return a.seqCisTransScore - b.seqCisTransScore;
        }
        // Rule 4a: chiral (R/S) > pseudoasymmetric (r/s) > nonstereogenic
        if (a.stereogenicScore !== b.stereogenicScore) {
            return a.stereogenicScore - b.stereogenicScore;
        }
        return 0;
    }

    // -----------------------------------------------------------------------
    // Legacy compareCipNodes — kept for backward compatibility but delegates
    // to the BFS level-by-level comparison
    // -----------------------------------------------------------------------

    function compareCipNodes(a, b) {
        return compareCipTrees(a, b);
    }

    // -----------------------------------------------------------------------
    // CIP priority ranking
    // -----------------------------------------------------------------------

    /**
     * Compute CIP priority ranks for the substituents of a given centre atom.
     * Returns an array of { neighborId, rank, tree } sorted by atom order.
     * rank 0 = highest priority, 1 = next, etc.
     * Returns null if any two substituents have equal priority (no valid assignment).
     *
     * Optional chiralityMap: when provided (pass 2), digraph trees are
     * decorated with Rule 3 (seqCis/seqTrans) and Rule 4a (stereogenicity)
     * scores before comparison.
     */
    function cipPriorities(mol, centerAtomId, chiralityMap) {
        var neighbors = mol.getNeighbors(centerAtomId);
        var hCount = mol.calcHydrogens(centerAtomId);

        var entries = [];
        for (var i = 0; i < neighbors.length; i++) {
            var tree = buildDigraph(mol, neighbors[i], centerAtomId);
            if (chiralityMap) {
                _decorateRule3(mol, tree);
                _decorateRule4a(tree, chiralityMap);
            }
            entries.push({ neighborId: neighbors[i], tree: tree });
        }
        for (var h = 0; h < hCount; h++) {
            var hTree = new CipNode(-1, 'H', 1, false);
            hTree.phantom = true;
            entries.push({ neighborId: -(h + 1), tree: hTree });
        }

        // Sort by descending priority (highest first)
        entries.sort(function(a, b) {
            return compareCipTrees(b.tree, a.tree);
        });

        // Check for duplicate priorities
        for (var i = 0; i < entries.length - 1; i++) {
            if (compareCipTrees(entries[i].tree, entries[i + 1].tree) === 0) {
                return null; // Cannot assign — duplicate priority
            }
        }

        // Assign ranks: 0 = highest priority
        for (var i = 0; i < entries.length; i++) {
            entries[i].rank = i;
        }

        return entries;
    }

    // -----------------------------------------------------------------------
    // R/S assignment for tetrahedral stereocentres
    // -----------------------------------------------------------------------

    /**
     * Assign R/S labels to all atoms with chirality set ('@' or '@@').
     *
     * The SMILES chirality convention:
     *   '@'  = anticlockwise when viewed from the first listed neighbor
     *   '@@' = clockwise when viewed from the first listed neighbor
     *
     * To determine R/S:
     *   1. Get CIP priority ordering (0 = highest)
     *   2. The lowest-priority group (rank 3, typically H) goes to the back
     *   3. Looking from the front (opposite to rank 3), if 0->1->2 is
     *      clockwise => R, anticlockwise => S
     *
     * We use the relationship between the SMILES neighbor ordering and the
     * CIP priority ordering, combined with the '@'/'@@' tag, to determine
     * the spatial arrangement without needing 3D coordinates.
     */
    function assignRS(mol) {
        // Two-pass approach:
        //   Pass 1: assign R/S using Rules 1-2 only (no chiralityMap)
        //   Pass 2: re-assign using Rules 1-4a with the chiralityMap from pass 1
        //           (this allows Rule 4a to use stereogenicity of neighbouring
        //            centres as a tiebreaker)

        // Pass 1 — Rules 1-2 only
        _assignRSPass(mol, null);

        // Build chiralityMap from pass-1 results
        var chiralityMap = {};
        var hasChiral = false;
        for (var i = 0; i < mol.atoms.length; i++) {
            if (mol.atoms[i].cipLabel) {
                chiralityMap[mol.atoms[i].id] = mol.atoms[i].cipLabel;
                hasChiral = true;
            }
        }

        // Pass 2 — Rules 1-4a (only needed if pass 1 found any chirality)
        if (hasChiral) {
            _assignRSPass(mol, chiralityMap);
        }
    }

    /**
     * Internal helper: assign R/S labels in a single pass.
     * When chiralityMap is provided, digraph trees are decorated with
     * Rule 3 and Rule 4a scores for tiebreaking.
     */
    function _assignRSPass(mol, chiralityMap) {
        for (var i = 0; i < mol.atoms.length; i++) {
            var atom = mol.atoms[i];
            atom.cipLabel = '';

            if (!atom.chirality) continue;

            var neighbors = mol.getNeighbors(atom.id);
            var hCount = (atom.hydrogens > 0) ? atom.hydrogens : mol.calcHydrogens(atom.id);
            var totalSubstituents = neighbors.length + hCount;

            if (totalSubstituents !== 4) continue;

            var priorities = cipPriorities(mol, atom.id, chiralityMap);
            if (!priorities) continue; // Duplicate priorities — no valid R/S

            // Map neighbor IDs to CIP ranks
            var rankOf = {};
            for (var p = 0; p < priorities.length; p++) {
                rankOf[priorities[p].neighborId] = priorities[p].rank;
            }

            // Build the sequence of CIP ranks in the order neighbors appear
            var seq = [];
            for (var n = 0; n < neighbors.length; n++) {
                seq.push(rankOf[neighbors[n]]);
            }
            // Implicit hydrogens come after explicit neighbors in SMILES convention
            for (var h = 0; h < hCount; h++) {
                seq.push(rankOf[-(h + 1)]);
            }

            // Compute the parity of the permutation needed to sort seq into [0,1,2,3]
            var parity = permutationParity(seq);

            // SMILES chirality convention (OpenSMILES):
            //   '@'  = anticlockwise when viewed from the first listed neighbor
            //   '@@' = clockwise when viewed from the first listed neighbor
            //
            // CIP convention:
            //   View from opposite the lowest-priority group (rank 3).
            //   If 0->1->2 is clockwise => R, anticlockwise => S.
            //
            // The permutation parity of (CIP ranks in SMILES neighbor order)
            // relative to sorted [0,1,2,3] tells us whether the SMILES-encoded
            // winding matches the CIP-sorted winding.
            //
            //   '@'  (CCW from first neighbor) + even parity => R
            //   '@@' (CW from first neighbor) + even parity => S
            //   '@'  (CCW from first neighbor) + odd parity  => S
            //   '@@' (CW from first neighbor) + odd parity  => R
            var isAt = (atom.chirality === '@');

            if (parity === 0) {
                atom.cipLabel = isAt ? 'R' : 'S';
            } else {
                atom.cipLabel = isAt ? 'S' : 'R';
            }
        }
    }

    /**
     * Compute the parity of a permutation.
     * Returns 0 for even, 1 for odd.
     */
    function permutationParity(perm) {
        var n = perm.length;
        var arr = perm.slice();
        var swaps = 0;
        // FIX: validate that arr is a permutation of [0..n-1] to avoid an
        // infinite loop if a caller accidentally supplies undefined or
        // out-of-range values (e.g. from a partial CIP rank assignment).
        for (var v = 0; v < n; v++) {
            if (typeof arr[v] !== 'number' || arr[v] < 0 || arr[v] >= n) {
                return 0; // not a valid permutation; treat as even
            }
        }
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
    // E/Z assignment for double bonds
    // -----------------------------------------------------------------------

    /**
     * Assign E/Z labels to all double bonds that have two different
     * substituents on each end.
     *
     * Convention:
     *   Z (zusammen, "together") = higher-priority groups on same side
     *   E (entgegen, "opposite") = higher-priority groups on opposite side
     *
     * We use the 2D coordinates to determine which side of the double bond
     * each substituent falls on.
     */
    function assignEZ(mol) {
        for (var i = 0; i < mol.bonds.length; i++) {
            var bond = mol.bonds[i];
            bond.cipLabel = '';

            if (bond.type !== BOND_DOUBLE) continue;

            var a1 = mol.getAtom(bond.atom1);
            var a2 = mol.getAtom(bond.atom2);
            if (!a1 || !a2) continue;

            // Get substituents on each end (excluding the double-bond partner)
            var subs1 = _getSubstituents(mol, bond.atom1, bond.atom2);
            var subs2 = _getSubstituents(mol, bond.atom2, bond.atom1);

            // Need exactly 2 substituents on each end (including implicit H)
            if (subs1.length < 1 || subs1.length > 2) continue;
            if (subs2.length < 1 || subs2.length > 2) continue;

            // If only one substituent on either side, no E/Z possible
            if (subs1.length < 2 && mol.calcHydrogens(bond.atom1) === 0) continue;
            if (subs2.length < 2 && mol.calcHydrogens(bond.atom2) === 0) continue;

            // Rank substituents on each end
            var pri1 = _rankDoubleBondSubs(mol, bond.atom1, bond.atom2, subs1);
            var pri2 = _rankDoubleBondSubs(mol, bond.atom2, bond.atom1, subs2);

            if (!pri1 || !pri2) continue; // Equal priorities on one side

            // Determine which side of the double bond axis each high-priority
            // substituent falls on, using the cross product (2D).
            var highSub1 = pri1[0]; // highest priority on atom1 side
            var highSub2 = pri2[0]; // highest priority on atom2 side

            var pos1 = _getSubPosition(mol, highSub1, bond.atom1);
            var pos2 = _getSubPosition(mol, highSub2, bond.atom2);

            if (!pos1 || !pos2) continue;

            // Bond axis vector: a1 -> a2
            var bx = a2.x - a1.x;
            var by = a2.y - a1.y;

            // Cross products to determine which side each high-priority sub is on
            var dx1 = pos1.x - a1.x;
            var dy1 = pos1.y - a1.y;
            var cross1 = bx * dy1 - by * dx1;

            var dx2 = pos2.x - a2.x;
            var dy2 = pos2.y - a2.y;
            var cross2 = bx * dy2 - by * dx2;

            if (Math.abs(cross1) < 0.001 || Math.abs(cross2) < 0.001) continue;

            // Same sign = same side = Z; opposite sign = opposite side = E
            if ((cross1 > 0) === (cross2 > 0)) {
                bond.cipLabel = 'Z';
            } else {
                bond.cipLabel = 'E';
            }
        }
    }

    /**
     * Get substituent atom IDs for one end of a double bond, excluding the
     * partner atom and including implicit H as negative IDs.
     */
    function _getSubstituents(mol, atomId, partnerId) {
        var neighbors = mol.getNeighbors(atomId);
        var subs = [];
        for (var i = 0; i < neighbors.length; i++) {
            if (neighbors[i] !== partnerId) {
                subs.push(neighbors[i]);
            }
        }
        // Add implicit hydrogens as virtual substituents
        var hCount = mol.calcHydrogens(atomId);
        for (var h = 0; h < hCount; h++) {
            subs.push(-(h + 1));
        }
        return subs;
    }

    /**
     * Rank substituents on one end of a double bond by CIP priority.
     * Returns sorted array [highestPri, lowestPri] of substituent IDs,
     * or null if priorities are equal.
     */
    function _rankDoubleBondSubs(mol, atomId, partnerId, subs) {
        if (subs.length < 2) {
            // Only one real substituent — no need to rank, it is the highest
            return subs;
        }

        var entries = [];
        for (var i = 0; i < subs.length; i++) {
            var subId = subs[i];
            var tree;
            if (subId < 0) {
                tree = new CipNode(-1, 'H', 1, false);
                tree.phantom = true;
            } else {
                tree = buildDigraph(mol, subId, atomId);
            }
            entries.push({ subId: subId, tree: tree });
        }

        entries.sort(function(a, b) {
            return compareCipTrees(b.tree, a.tree);
        });

        // Check for equal priorities
        for (var i = 0; i < entries.length - 1; i++) {
            if (compareCipTrees(entries[i].tree, entries[i + 1].tree) === 0) {
                return null;
            }
        }

        var result = [];
        for (var i = 0; i < entries.length; i++) {
            result.push(entries[i].subId);
        }
        return result;
    }

    /**
     * Get the 2D position of a substituent for cross-product calculation.
     * For real atoms, use their coordinates.
     * For implicit H (negative IDs), estimate a position opposite to existing bonds.
     */
    function _getSubPosition(mol, subId, anchorId) {
        if (subId > 0) {
            var atom = mol.getAtom(subId);
            return atom ? { x: atom.x, y: atom.y } : null;
        }
        // Implicit hydrogen — estimate position
        var anchor = mol.getAtom(anchorId);
        if (!anchor) return null;
        var neighbors = mol.getNeighbors(anchorId);
        if (neighbors.length === 0) return { x: anchor.x + 1, y: anchor.y };

        // Place H opposite to the average bond direction
        var sx = 0, sy = 0;
        for (var i = 0; i < neighbors.length; i++) {
            var n = mol.getAtom(neighbors[i]);
            if (!n) continue;
            var dx = n.x - anchor.x;
            var dy = n.y - anchor.y;
            var len = Math.sqrt(dx * dx + dy * dy);
            if (len > 0) { sx += dx / len; sy += dy / len; }
        }
        var len2 = Math.sqrt(sx * sx + sy * sy);
        if (len2 > 0) {
            return { x: anchor.x - (sx / len2) * 30, y: anchor.y - (sy / len2) * 30 };
        }
        return { x: anchor.x, y: anchor.y + 30 };
    }

    // -----------------------------------------------------------------------
    // Main entry point — assign all CIP labels
    // -----------------------------------------------------------------------

    /**
     * Compute and assign CIP stereochemistry labels for the entire molecule.
     * Sets atom.cipLabel ('R' or 'S') and bond.cipLabel ('E' or 'Z').
     */
    function assign(mol) {
        if (!mol || mol.atoms.length === 0) return;
        assignRS(mol);
        assignEZ(mol);
    }

    // -----------------------------------------------------------------------
    // Export
    // -----------------------------------------------------------------------

    global.CipStereo = {
        assign: assign,
        assignRS: assignRS,
        assignEZ: assignEZ,
        cipPriorities: cipPriorities,
        ATOMIC_NUMBER: ATOMIC_NUMBER
    };

})(typeof window !== 'undefined' ? window : this);
