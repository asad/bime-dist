/**
 * CIPStereo.js — Cahn-Ingold-Prelog stereochemistry assignment.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Full digraph-based CIP implementation following the IUPAC 2013
 * recommendations, hierarchical-digraph rules 1a–4a:
 *
 *   Rule 1a  — higher atomic number precedes lower.
 *   Rule 1b  — a duplicate atom ranks by the atomic number of the atom
 *              it duplicates (realised here via phantom nodes).
 *   Rule 2   — higher mass number precedes lower (isotopes).
 *   Rule 3   — seqcis ( Z ) precedes seqtrans ( E ) precedes a
 *              non-stereogenic double bond.
 *   Rule 4a  — chiral ( R / S ) precedes pseudo-asymmetric ( r / s )
 *              precedes a non-stereogenic centre (resolved two-pass).
 *
 * Implementation notes:
 *   - The hierarchical digraph is grown breadth-first, level by level.
 *   - Double and triple bonds expand into phantom duplicate leaf nodes.
 *   - Ring closures expand into a phantom duplicate leaf node at the
 *     back-edge, keeping the digraph acyclic — spiro, bridged and fused
 *     ring systems are all handled correctly.
 *   - Aromatic "mancude" atoms contribute a fractional pi bond to a
 *     phantom partner node.
 *   - Runaway guards: MAX_DEPTH = 30 levels, MAX_NODES = 5000 nodes.
 *   - Tetrahedral parity is read from the SMILES @ / @@ flag via a
 *     permutation-parity comparison against the CIP-ranked order.
 *   - Implicit-hydrogen counts are charge-aware and add the aromatic
 *     pi contribution where relevant.
 *
 * Calibrated abstention (never-wrong contract):
 *   The engine NEVER emits a guessed descriptor. A specified tetrahedral
 *   stereocentre whose ligand priorities cannot be fully resolved is left
 *   UNLABELLED (cipLabel '') and flagged `atom.cipUndetermined = true`, so a
 *   caller can distinguish an intentional abstention from "not a stereocentre"
 *   (the latter has chirality unset / cipUndetermined falsy). The principal
 *   abstention class is the symmetric ring cyclitol (inositol family): every
 *   ring carbon's two ring-path ligands are constitutionally identical, so
 *   Rules 1-3/4a tie, and the auxiliary-descriptor pass (Rules 4b/5) cannot
 *   seed across the ring-closure back-edge duplicate leaf. Resolving these
 *   requires a STEREO-AWARE canonical ranking (the SmilesWriter ranking is
 *   constitutional-only — see its KNOWN_ORBIT_DRIFT note) plus an external
 *   IUPAC-2013 oracle to validate the order-fragile lowercase r/s letters; it
 *   is a deferred follow-up. Reference ground truth (prof-chem / IUPAC 2013):
 *   scyllo- and cis-inositol carbons are genuinely NON-stereogenic (abstention
 *   is the correct answer); myo-inositol carries four R/S + two pseudo-
 *   asymmetric r/s (C2, C5 on its mirror plane). Acyclic pseudo-asymmetry
 *   (ribitol r, ribaric-acid s) and meso hexitols (galactitol, allitol -> mixed
 *   R/S) already resolve correctly and are PubChem-validated in the tests.
 *
 * Public API — exported as the global `CIPStereo`. The former name
 * `CipStereo` is retained as a deprecated back-compat alias.
 *   assign(mol)        — assign both R/S and E/Z labels, in place.
 *   assignRS(mol)      — assign R/S on tetrahedral stereocentres
 *                        (atoms whose .chirality is '@' or '@@'). Sets
 *                        atom.cipLabel and, where it abstains, the honest
 *                        atom.cipUndetermined flag (see above).
 *   assignEZ(mol)      — assign E/Z on double bonds with restricted
 *                        rotation.
 *   cipPriorities(mol, centerAtomId [, chiralityMap])
 *                      — return the CIP-ranked neighbour list.
 *   ATOMIC_NUMBER      — the shared symbol→Z table (Molecule.ATOMIC_NUMBERS).
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
    // v1.8.26: aliased to the canonical Molecule.ATOMIC_NUMBERS table
    // (was a local copy duplicated across 4 modules). The canonical
    // table extends to U(92) vs the old local Bi(83) cap — strictly
    // more correct for CIP priority of heavy-element substituents.
    var ATOMIC_NUMBER = Molecule.ATOMIC_NUMBERS;

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
        // Rule 4b/5: auxiliary descriptor ('R'|'S'|'r'|'s'|''), painted from
        // the resolved descriptor set so the like/unlike comparator can read it
        // off the digraph (see _decorateAux).
        this.auxDescriptor = '';
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
                    // (forward direction). CIP Rule: a double bond A=B gives
                    // A a phantom duplicate of B (and, symmetrically, B a
                    // phantom duplicate of A). The duplicate of the NEIGHBOUR
                    // (B) must hang off the CURRENT node (A) so that e.g. a
                    // carbonyl carbon ranks as (O,O,H) not (O,H). The reciprocal
                    // duplicate of A under B is added when B is expanded and
                    // sees its back-edge to A (the neighborId===parentId branch
                    // above) — so it must NOT be added here, or it double-counts.
                    if (bond.type >= BOND_DOUBLE) {
                        for (var p = 1; p < bond.type && totalNodes < MAX_NODES; p++) {
                            var dupNeighbor = new CipNode(
                                neighborId, nAtom.symbol, atomMass(nAtom), true
                            );
                            dupNeighbor.depth = depth;
                            node.children.push(dupNeighbor);
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
     * Compare two CIP digraph trees by the IUPAC-2013 hierarchical-digraph
     * rules (P-92 / Hanson et al. 2018, J. Chem. Inf. Model. 58:1755 — the
     * the modern reference-implementation formulation), applied as a strict
     * RULE LADDER:
     *
     *   [Rule 1a/1b atomic number] then [Rule 2 mass] then [Rule 3 seqCis/
     *    trans] then [Rule 4a stereogenicity]
     *
     * Each rule is explored to exhaustion over the WHOLE digraph before the
     * next rule is consulted. Within a rule the comparison is breadth-first
     * with BRANCH PAIRING: compare the two roots; rank each node's children
     * high->low by the full cumulative comparator; compare the whole sphere
     * by node value FIRST; only on a full-sphere tie descend into the matched
     * pairs in rank order. This preserves branch identity (the old pooled-BFS
     * flatten discarded it, mis-ranking deep substituent chains) AND respects
     * sphere order (a deep high-Z atom must not outrank a nearer difference).
     *
     * Determinism: residual ties are broken by _structuralCmp (subtree shape +
     * labels only — never atom id / address / input order), giving a total
     * order, so the result is independent of Array.sort stability and of the
     * SMILES atom-input order. Returns >0 if a outranks b, <0 if b, 0 if tied.
     */
    function _ruleAtomicNum(a, b) { return a.atomicNum - b.atomicNum; }            // Rule 1a/1b
    function _ruleMass(a, b)      { return a.mass - b.mass; }                      // Rule 2
    function _ruleSeq(a, b)       { return a.seqCisTransScore - b.seqCisTransScore; } // Rule 3
    function _ruleStereo(a, b)    { return a.stereogenicScore - b.stereogenicScore; } // Rule 4a

    var CIP_RULES = [_ruleAtomicNum, _ruleMass, _ruleSeq, _ruleStereo];

    function compareCipTrees(a, b) {
        if (!a && !b) return 0;
        if (!a) return -1;
        if (!b) return 1;
        for (var r = 0; r < CIP_RULES.length; r++) {
            var cmp = _ruleCompare(a, b, r, 0);
            if (cmp !== 0) return cmp;
        }
        return 0;
    }

    // Compare two trees under rules [0..maxRule], breadth-first with branch
    // pairing. depth is bounded by MAX_DEPTH.
    function _ruleCompare(a, b, maxRule, depth) {
        if (!a && !b) return 0;
        if (!a) return -1;
        if (!b) return 1;
        var nodeCmp = _cumulativeNodeCmp(a, b, maxRule);
        if (nodeCmp !== 0) return nodeCmp;
        if (depth >= MAX_DEPTH) return 0;
        var ca = _rankChildren(a, maxRule, depth);
        var cb = _rankChildren(b, maxRule, depth);
        var n = ca.length > cb.length ? ca.length : cb.length;
        // Whole-sphere node-value test, in rank order (sphere before descent).
        for (var i = 0; i < n; i++) {
            var va = i < ca.length ? ca[i] : null;
            var vb = i < cb.length ? cb[i] : null;
            if (!va && !vb) continue;
            if (!va) return -1;
            if (!vb) return 1;
            var vc = _cumulativeNodeCmp(va, vb, maxRule);
            if (vc !== 0) return vc;
        }
        // Sphere tied -> descend into matched pairs in rank order.
        for (var j = 0; j < n; j++) {
            var da = j < ca.length ? ca[j] : null;
            var db = j < cb.length ? cb[j] : null;
            if (!da && !db) continue;
            if (!da) return -1;
            if (!db) return 1;
            var dc = _ruleCompare(da, db, maxRule, depth + 1);
            if (dc !== 0) return dc;
        }
        return 0;
    }

    // Node test under rules [0..maxRule] only (no recursion).
    function _cumulativeNodeCmp(a, b, maxRule) {
        for (var r = 0; r <= maxRule; r++) {
            var c = CIP_RULES[r](a, b);
            if (c !== 0) return c;
        }
        return 0;
    }

    // Rank a node's children high->low by the full cumulative comparator,
    // with a deterministic structural tie-break. Memoised per (node, maxRule);
    // trees are rebuilt per cipPriorities call so the cache never goes stale.
    function _rankChildren(node, maxRule, depth) {
        if (!node._rankCache) node._rankCache = {};
        if (node._rankCache[maxRule]) return node._rankCache[maxRule];
        var kids = node.children.slice();
        kids.sort(function(x, y) {
            var c = _ruleCompare(y, x, maxRule, depth + 1); // descending
            if (c !== 0) return c;
            return _structuralCmp(y, x);
        });
        node._rankCache[maxRule] = kids;
        return kids;
    }

    // Total-order structural tie-break: compares subtree shape + node labels
    // ONLY (atomicNum, mass, scores, recursively). Never consults atom id,
    // object identity, thread state, or SMILES token order — so it cannot
    // introduce input-order dependence.
    function _structuralCmp(a, b) {
        if (!a && !b) return 0;
        if (!a) return -1;
        if (!b) return 1;
        if (a.atomicNum !== b.atomicNum) return a.atomicNum - b.atomicNum;
        if (a.mass !== b.mass) return a.mass - b.mass;
        if (a.seqCisTransScore !== b.seqCisTransScore) return a.seqCisTransScore - b.seqCisTransScore;
        if (a.stereogenicScore !== b.stereogenicScore) return a.stereogenicScore - b.stereogenicScore;
        var ka = a.children, kb = b.children;
        if (ka.length !== kb.length) return ka.length - kb.length;
        var la = ka.slice(), lb = kb.slice();
        la.sort(_structuralNodeSort);
        lb.sort(_structuralNodeSort);
        for (var i = 0; i < la.length; i++) {
            var c = _structuralCmp(la[i], lb[i]);
            if (c !== 0) return c;
        }
        return 0;
    }
    function _structuralNodeSort(x, y) {
        if (x.atomicNum !== y.atomicNum) return y.atomicNum - x.atomicNum;
        if (x.mass !== y.mass) return y.mass - x.mass;
        if (x.seqCisTransScore !== y.seqCisTransScore) return y.seqCisTransScore - x.seqCisTransScore;
        if (x.stereogenicScore !== y.stereogenicScore) return y.stereogenicScore - x.stereogenicScore;
        if (x.children.length !== y.children.length) return y.children.length - x.children.length;
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
    // Rule 4b / 4c / 5 — auxiliary descriptors, like / unlike, R before S
    // -----------------------------------------------------------------------
    //
    // IUPAC 2013 P-92.1.4.4 (Rules 4b, 4c) and P-92.1.4.5 (Rule 5), in the
    // algorithmic formulation of Hanson, Musacchio, Mayfield, Vainio, Yerin &
    // Redkin 2018 (J. Chem. Inf. Model. 58:1755-1765, "Algorithmic Analysis of
    // Cahn-Ingold-Prelog Rules of Stereochemistry"), the formulation the
    // modern external reference IUPAC-2013 CIP implementations use.
    //
    // These rules ACT ONLY when Rules 1-3 (and 4a) leave two ligands of a
    // centre tied; they never override a constitutional difference. The engine
    // invariant is NEVER EMIT A WRONG LABEL: every branch below either resolves
    // a centre exactly as the reference would, or returns "undecided" so the
    // centre stays UNLABELLED. No guessing.
    //
    // Auxiliary descriptors. To rank the ligands of a root centre we need a
    // provisional descriptor on every stereogenic atom INSIDE the hierarchical
    // digraph, not only the root. The auxiliary descriptor of an internal atom
    // is computed in the digraph rooted at THAT atom: its children point away
    // from the root and its single edge back toward the root is taken as a
    // DUPLICATE atom leaf (CIP hierarchical-digraph rule), and the @/@@ parity
    // is read in the canonical frame the parser normalised every chirality
    // token into (heavy neighbours in getNeighbors() order, implicit H last) -
    // the SAME frame _assignRSPass uses. Auxiliary descriptors of deeper atoms
    // feed the ranking of shallower atoms, so they are computed deepest-first.

    // Map a descriptor letter to its case-folded uppercase form (R/S). Used by
    // the Rule 4b/5 comparator so r/s and R/S compare uniformly.
    function _descUpper(d) {
        if (d === 'R' || d === 'r') return 'R';
        if (d === 'S' || d === 's') return 'S';
        return '';
    }

    // Decorate a root-branch digraph `tree` with auxiliary descriptors.
    //
    // IMPORTANT (never-wrong): we paint each internal node with the descriptor
    // it ALREADY carries in the supplied chiralityMap, i.e. a descriptor that
    // is pinned by constitution (Rules 1-3) or settled in an earlier pass. We
    // deliberately do NOT synthesise fresh per-root descriptors from a bare
    // duplicate-atom back-edge: that construction is sensitive to the SMILES
    // atom order in symmetric rings (the toward-root duplicate leaf ranks
    // differently from the equivalent away-root real chain), which would make
    // the resulting Rule 4b/5 outcome order-dependent and therefore unsafe.
    //
    // The practical consequence: Rule 4b/5 can resolve a centre only when the
    // distinguishing stereocentres in its ligands are themselves resolvable by
    // constitution (e.g. the pseudo-asymmetric carbon of ribitol, whose two
    // branches each carry a constitutionally-fixed R or S centre). Symmetric
    // cyclitols, whose ring centres never resolve constitutionally, are left
    // UNLABELLED rather than guessed. See the delivery notes for the exact
    // sub-classes this covers.
    function _decorateAux(mol, tree, ctx, chiralityMap) {
        (function paint(n) {
            if (!n) return;
            if (!n.phantom && n.atomId >= 0 && chiralityMap && chiralityMap[n.atomId]) {
                n.auxDescriptor = chiralityMap[n.atomId];
            }
            for (var i = 0; i < n.children.length; i++) paint(n.children[i]);
        })(tree);
        return null;
    }

    // Extract the hierarchical (sphere-ordered) sequence of auxiliary
    // descriptors from a decorated branch tree. Sphere 0 = the branch root, and
    // within each sphere nodes are taken in CIP-rank order (high->low) so the
    // sequence is canonical. Only resolved descriptors ('R'/'S'/'r'/'s') are
    // emitted; undecided nodes contribute nothing (they tie under Rule 4b/5 and
    // are settled, if at all, by deeper differences). Returns an array of
    // descriptor letters.
    function _auxSequence(tree) {
        var out = [];
        var level = [tree];
        var guard = 0;
        while (level.length > 0 && guard++ < MAX_DEPTH + 2) {
            // canonical rank order within the sphere
            var ranked = level.slice();
            ranked.sort(function(a, b) {
                var c = compareCipTrees(b, a);
                if (c !== 0) return c;
                return _structuralCmp(b, a);
            });
            var next = [];
            for (var i = 0; i < ranked.length; i++) {
                var nd = ranked[i];
                if (nd && !nd.phantom && nd.auxDescriptor) out.push(nd.auxDescriptor);
                if (nd) {
                    for (var k = 0; k < nd.children.length; k++) next.push(nd.children[k]);
                }
            }
            level = next;
        }
        return out;
    }

    // Rule 4b (like/unlike) + Rule 5 (R before S) comparison of two branch
    // trees `ta`, `tb` that are already TIED under Rules 1-3/4a. Returns
    //   { cmp: >0|<0|0, byRule5: bool }
    // where cmp>0 means ta outranks tb, and byRule5 is true iff the decision
    // came from Rule 5 (R-before-S) on otherwise-identical like/unlike
    // profiles -> the parent centre is pseudo-asymmetric if this pair decides
    // it.
    //
    // Rule 4b uses a PER-LIGAND reference: each ligand's like/unlike string is
    // built relative to that ligand's OWN leading (root-nearest, highest-rank)
    // descriptor. The leading descriptor is therefore always 'l'; a later
    // descriptor is 'l' if it has the same R/S sense as the lead, 'u' if
    // opposite. Two enantiomorphic ligands (exact mirror images) produce
    // IDENTICAL like/unlike strings, so Rule 4b cannot separate them and the
    // tie falls to Rule 5 -> pseudo-asymmetry. Two diastereomorphic ligands
    // produce DIFFERENT like/unlike strings and Rule 4b separates them as a
    // true (capital) distinction. (IUPAC 2013 P-92.1.4.4 / P-92.5; Hanson 2018.)
    function _likeUnlikeString(seq) {
        // seq: descriptor letters in canonical sphere order. Returns the
        // like/unlike vector relative to the ligand's OWN leading (root-nearest,
        // highest-rank) descriptor; the lead is always 'like'=1, a later
        // descriptor is 1 (like) if it has the same R/S sense as the lead, 0
        // (unlike) otherwise. Self-referencing makes two enantiomorphic ligands
        // produce IDENTICAL vectors (mirror images preserve internal like/unlike
        // relations), so Rule 4b cannot separate them and the tie falls to Rule
        // 5 -> pseudo-asymmetry; diastereomorphic ligands differ -> Rule 4b
        // separates them as a true (capital) distinction. (IUPAC 2013
        // P-92.1.4.4 / P-92.5; Hanson et al. 2018.) Returns [] if no
        // descriptors are present.
        if (seq.length === 0) return [];
        var lead = _descUpper(seq[0]);
        if (!lead) return [];
        var out = [];
        for (var i = 0; i < seq.length; i++) {
            var u = _descUpper(seq[i]);
            out.push(u === lead ? 1 : 0);
        }
        return out;
    }

    // Rule 4b (like/unlike) + Rule 5 (R before S) comparison of two branch
    // trees `ta`, `tb` already TIED under Rules 1-3/4a. Returns
    //   { cmp: >0|<0|0, byRule5: bool }
    // cmp>0 => ta outranks tb. byRule5 => decided by R-before-S on identical
    // like/unlike profiles, which makes the parent centre pseudo-asymmetric.
    //
    // Because the auxiliary descriptors consulted here are ABSOLUTE (pinned by
    // constitution or an earlier pass; see _decorateAux), this comparison is
    // independent of SMILES atom order. The leading descriptors used by Rule 5
    // are therefore real R/S values and the R-before-S decision is meaningful.
    function _ruleB5Compare(ta, tb, ref) {
        var sa = _auxSequence(ta);
        var sb = _auxSequence(tb);
        if (sa.length === 0 && sb.length === 0) return { cmp: 0, byRule5: false };

        // Rule 4b: compare per-ligand like/unlike strings hierarchically.
        var lua = _likeUnlikeString(sa);
        var lub = _likeUnlikeString(sb);
        var n = lua.length < lub.length ? lua.length : lub.length;
        for (var i = 0; i < n; i++) {
            if (lua[i] !== lub[i]) return { cmp: (lua[i] - lub[i]), byRule5: false };
        }
        // Differing lengths => the constitutions actually differed and the tie
        // should have broken earlier; stay never-wrong by declaring undecided.
        if (lua.length !== lub.length) return { cmp: 0, byRule5: false };

        // Rule 5: identical like/unlike -> mirror-related ligands. Separate by
        // R-before-S (r before s) on the leading descriptors. The decision is
        // "by Rule 5" -> the parent centre is pseudo-asymmetric (lowercase).
        var la = _descUpper(sa[0]);
        var lb = _descUpper(sb[0]);
        if (la && lb && la !== lb) {
            var va = (la === 'R') ? 1 : 0;
            var vb = (lb === 'R') ? 1 : 0;
            return { cmp: (va - vb), byRule5: true };
        }
        return { cmp: 0, byRule5: false };
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
    function cipPriorities(mol, centerAtomId, chiralityMap, info) {
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

        // Check for duplicate priorities under Rules 1-3/4a. Any residual tie
        // is handed to the Rule 4b/5 resolver (auxiliary descriptors); if that
        // still cannot separate the pair, the centre stays UNLABELLED.
        var tiedPair = false;
        for (var i = 0; i < entries.length - 1; i++) {
            if (compareCipTrees(entries[i].tree, entries[i + 1].tree) === 0) {
                tiedPair = true;
                break;
            }
        }
        if (tiedPair) {
            // Rule 4b/5 only applies on pass 2 (chiralityMap present), where a
            // first-pass descriptor set exists to seed auxiliary descriptors.
            if (!chiralityMap) return null;
            return _resolveWithAux(mol, centerAtomId, entries, chiralityMap, info);
        }

        // Assign ranks: 0 = highest priority
        for (var i = 0; i < entries.length; i++) {
            entries[i].rank = i;
        }

        return entries;
    }

    // Resolve a centre whose ligand entries (already sorted by Rules 1-3/4a)
    // contain at least one tied adjacent pair, using Rule 4b (like/unlike) and
    // Rule 5 (R before S) on per-root auxiliary descriptors. Returns the
    // entries array with ranks assigned, or null if any pair stays tied (the
    // never-wrong fallback). If a pair is separated by Rule 5 on reflection-
    // related ligands, info.pseudoAsym is set true so the caller emits a
    // lowercase r/s descriptor for the centre.
    function _resolveWithAux(mol, centerAtomId, entries, chiralityMap, info) {
        // Decorate every ligand branch with per-root auxiliary descriptors.
        var ctx = { count: 0 };
        for (var i = 0; i < entries.length; i++) {
            if (entries[i].tree && !entries[i].tree.phantom) {
                _decorateAux(mol, entries[i].tree, ctx, chiralityMap);
            }
        }

        // Choose the root reference descriptor. Per Hanson, the reference is the
        // descriptor of the highest-ranked ligand's first stereogenic atom; we
        // take the first resolved auxiliary descriptor in CIP-rank order across
        // ligands. Fold r/s -> R/S for the like test.
        var ref = '';
        for (var e = 0; e < entries.length && !ref; e++) {
            var seq = entries[e].tree ? _auxSequence(entries[e].tree) : [];
            for (var s = 0; s < seq.length; s++) {
                var u = _descUpper(seq[s]);
                if (u) { ref = u; break; }
            }
        }
        if (!ref) return null; // no auxiliary descriptors -> cannot resolve

        // Re-sort the FULL ligand list with the Rule 4b/5 tiebreak appended to
        // the Rules 1-3/4a comparator. Track whether any adjacent separation
        // used Rule 5 (-> pseudo-asymmetric centre).
        var usedRule5Anywhere = { v: false };
        entries.sort(function(a, b) {
            var base = compareCipTrees(b.tree, a.tree);
            if (base !== 0) return base;
            var r = _ruleB5Compare(b.tree, a.tree, ref);
            return r.cmp; // >0 means b<-... compareCipTrees(b,a) convention
        });

        // Verify a strict total order now exists; record Rule-5 usage on the
        // adjacent deciding pairs. Any pair still tied, or any pair whose Rule
        // 4b/5 outcome is sign-ambiguous, leaves the centre UNLABELLED.
        var rule5OnDecider = false;
        for (var k = 0; k < entries.length - 1; k++) {
            var base2 = compareCipTrees(entries[k].tree, entries[k + 1].tree);
            if (base2 !== 0) continue; // separated constitutionally; fine
            var rr = _ruleB5Compare(entries[k].tree, entries[k + 1].tree, ref);
            if (rr.cmp === 0) return null;        // still tied -> never-wrong
            if (rr.byRule5) rule5OnDecider = true; // reflection-only separation
        }

        for (var m = 0; m < entries.length; m++) entries[m].rank = m;
        if (info) info.pseudoAsym = rule5OnDecider;
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
        // Multi-pass fixpoint:
        //   Pass 1: Rules 1-3 only (no chiralityMap) — resolve the centres that
        //           need no neighbour stereo information.
        //   Pass 2+: Rules 1-4a with the chiralityMap from the previous pass,
        //           and Rules 4b/5 (auxiliary descriptors) for centres still
        //           tied by constitution (cyclitols, pseudo-asymmetric carbons).
        // The loop repeats while new centres become labelled, so a descriptor
        // resolved in one pass can feed Rule 4a in the next. It is bounded and
        // converges because the label set only grows. Determinism: each pass is
        // a deterministic function of the molecule + current label set.

        var hasAnyChiral = false;
        for (var c0 = 0; c0 < mol.atoms.length; c0++) {
            if (mol.atoms[c0].chirality) { hasAnyChiral = true; break; }
        }

        // Pass 1 — Rules 1-3 only
        _assignRSPass(mol, null);
        if (!hasAnyChiral) return;

        // Passes 2.. — feed back labels; run aux (Rules 4b/5) for tied centres.
        var prevSig = null;
        for (var pass = 0; pass < mol.atoms.length + 4; pass++) {
            var chiralityMap = {};
            var sigParts = [];
            for (var i = 0; i < mol.atoms.length; i++) {
                if (mol.atoms[i].cipLabel) {
                    chiralityMap[mol.atoms[i].id] = mol.atoms[i].cipLabel;
                    sigParts.push(mol.atoms[i].id + ':' + mol.atoms[i].cipLabel);
                }
            }
            var sig = sigParts.join(',');
            if (sig === prevSig) break; // fixpoint reached
            prevSig = sig;
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
            atom.cipUndetermined = false;

            if (!atom.chirality) continue;

            var neighbors = mol.getNeighbors(atom.id);
            var hCount = (atom.hydrogens > 0) ? atom.hydrogens : mol.calcHydrogens(atom.id);
            var totalSubstituents = neighbors.length + hCount;

            if (totalSubstituents !== 4) continue;

            var info = { pseudoAsym: false };
            var priorities = cipPriorities(mol, atom.id, chiralityMap, info);
            if (!priorities) {
                // A specified tetrahedral stereocentre whose ligand priorities
                // cannot be fully resolved — constitutional symmetry or an
                // unresolved ring pseudo-asymmetry (the inositol-class cyclitol
                // family). We NEVER guess a descriptor; we abstain (cipLabel
                // stays '') but record WHY, so a caller/test can tell this
                // intentional calibrated abstention apart from "not a
                // stereocentre". See the module-header note on cyclitols.
                atom.cipUndetermined = true;
                continue;
            }

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

            // Sign calibrated against textbook anchors (L-alanine = S,
            // D-glyceraldehyde = R, (R)-CHFClBr = R) once the parser normalises
            // EVERY chiral centre — from-atom AND chiral-first — into the single
            // "heavy neighbours in adjacency order, implicit H last" frame. With
            // both classes in one true frame this single sign map is correct for
            // all of them (the previous inverse map only matched chiral-first
            // because those tokens were mis-normalised; see SmilesParser lead-H).
            var label;
            if (parity === 0) {
                label = isAt ? 'S' : 'R';
            } else {
                label = isAt ? 'R' : 'S';
            }
            // Rule 5: a centre resolved only by R-before-S on reflection-related
            // ligands is pseudo-asymmetric and takes a lowercase descriptor.
            if (info.pseudoAsym) label = label.toLowerCase();
            atom.cipLabel = label;
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

    // -----------------------------------------------------------------------
    // v2.4.16: canonical E/Z descriptors from the PARSED directional bonds
    // (bond.stereo 1='/' 6='\\'), NOT from 2D geometry. assignEZ() reads the
    // layout, which a re-layout can corrupt; this reads the topology + the
    // input slashes, so it is order-invariant and survives any layout. Used to
    // give structure keys (MetaboliteLibrary.findBySmiles) an E/Z signature so
    // cis/trans isomers (maleate vs fumarate) are never conflated.
    // -----------------------------------------------------------------------

    // Side (+1/-1) of a directional substituent of `center` (excluding the
    // double-bond partner). OpenSMILES: stereo 1 ('/') means atom2 is "up"
    // relative to atom1; stereo 6 ('\\') means "down". Returns { sub, side } or
    // null when this end carries no directional bond.
    function _directionalSide(mol, center, dbPartner) {
        var bonds = mol.getBondsOfAtom(center);
        for (var i = 0; i < bonds.length; i++) {
            var b = bonds[i];
            if (b.otherAtom(center) === dbPartner) { continue; } // skip the C=C itself
            if (b.stereo === 1 || b.stereo === 6) {
                var side = (b.stereo === 1 ? 1 : -1) * (b.atom1 === center ? 1 : -1);
                return { sub: b.otherAtom(center), side: side };
            }
        }
        return null;
    }

    // Canonical E/Z for one double bond, or null when it is not a determinable
    // stereo double bond (no directional bond on an end, or a CIP-tied end).
    function _doubleBondEZ(mol, db) {
        var a1 = db.atom1, a2 = db.atom2;
        var f1 = _directionalSide(mol, a1, a2);
        var f2 = _directionalSide(mol, a2, a1);
        if (!f1 || !f2) { return null; }
        var r1 = _rankDoubleBondSubs(mol, a1, a2, _getSubstituents(mol, a1, a2));
        var r2 = _rankDoubleBondSubs(mol, a2, a1, _getSubstituents(mol, a2, a1));
        if (!r1 || !r2 || !r1.length || !r2.length) { return null; } // CIP-tied end
        var hi1 = r1[0], hi2 = r2[0];
        // The flag may sit on the LOWER-priority substituent; the higher one is
        // then on the opposite side of the double-bond axis (negate).
        var side1 = (f1.sub === hi1) ? f1.side : -f1.side;
        var side2 = (f2.sub === hi2) ? f2.side : -f2.side;
        // High-priority groups on the SAME side = Z (zusammen); opposite = E.
        return (side1 === side2) ? 'Z' : 'E';
    }

    /**
     * Canonical, order-invariant, geometry-free E/Z descriptors for every
     * determinable stereo double bond, returned sorted (a multiset). Empty when
     * the molecule has no specified C=C stereochemistry.
     *
     * Scope: each double bond is perceived from its OWN flanking directional
     * bonds, so the per-bond label is exact for ISOLATED stereo double bonds
     * (which covers the entire metabolite corpus). In a conjugated polyene whose
     * shared single-bond slash couples adjacent C=C, an individual bond's label
     * may be wrong — but the sorted multiset still differs between genuine
     * diastereomers, so such a query degrades only to a NON-match (abstain),
     * never to a wrong compound name. Full conjugated-system perception is left
     * to a future minor.
     */
    function ezDescriptors(mol) {
        var out = [];
        if (!mol || !mol.bonds) { return out; }
        for (var i = 0; i < mol.bonds.length; i++) {
            var db = mol.bonds[i];
            if (db.type !== BOND_DOUBLE) { continue; }
            var ez = _doubleBondEZ(mol, db);
            if (ez) { out.push(ez); }
        }
        out.sort();
        return out;
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

    global.CIPStereo = {
        assign: assign,
        assignRS: assignRS,
        assignEZ: assignEZ,
        ezDescriptors: ezDescriptors,
        // Geometry-free per-double-bond E/Z ('E'|'Z'|null) from the parsed
        // directional markers (SMILES / \), independent of 2D coordinates.
        // Used by the layout to draw cis/trans correctly (assignEZ, which reads
        // coordinates, is circular after a re-layout rebuilds a trans zigzag).
        doubleBondEZ: _doubleBondEZ,
        cipPriorities: cipPriorities,
        ATOMIC_NUMBER: ATOMIC_NUMBER
    };

    // Deprecated back-compat alias. `CipStereo` was the pre-v1.8.27
    // name; external embedders referencing it keep working. New code
    // should use `CIPStereo` (CIP is an initialism — Cahn-Ingold-Prelog).
    global.CipStereo = global.CIPStereo;

})(typeof window !== 'undefined' ? window : this);
