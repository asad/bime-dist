/**
 * SmilesWriter.js — Molecule graph -> canonical SMILES string
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Generates canonical SMILES from a Molecule graph.
 * Features:
 *   - Input-order-invariant canonical atom ranking (structural seed +
 *     iterative colour refinement + iterative tie-breaking)
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

    // v1.8.26: STANDARD_VALENCES + getDefaultValence now live in
    // Molecule.js (single canonical home — they were byte-identical
    // here and in SmilesParser.js). Local aliases keep call sites terse.
    var STANDARD_VALENCES = Molecule.STANDARD_VALENCES;
    var getDefaultValence = Molecule.getDefaultValence;

    // Aromatic atoms that can be written in lowercase without brackets
    var AROMATIC_ORGANIC = { 'B': true, 'C': true, 'N': true, 'O': true,
                             'P': true, 'S': true };

    // Atomic numbers for canonical tie-breaking.
    // v1.8.26: aliased to the canonical Molecule.ATOMIC_NUMBERS table.
    // The old local copy was a sparse subset (common elements only);
    // the canonical full table is a strict superset, so canonical
    // ordering is unchanged for those elements and now also defined
    // for the rest.
    var ATOMIC_NUMBERS = Molecule.ATOMIC_NUMBERS;

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
    // Canonical atom ranking
    // -----------------------------------------------------------------------

    /**
     * Compute an input-order-invariant canonical ranking for a set of atom
     * IDs. Returns an object mapping atomId -> rank (lower rank = earlier in
     * the SMILES DFS). The ranking is a total order — every atom ends up
     * with a distinct rank — but the rank values are not necessarily dense.
     *
     * Three stages, none of which consult the atom-creation order:
     *
     *   1. Seed — each atom gets a structural invariant tuple (degree,
     *      element, charge, explicit valence, attached-H, isotope, ring
     *      membership). Bond orders enter through the valence term: that is
     *      what separates a carbonyl O from a hydroxyl O. The old summed-
     *      connectivity invariant ignored bond order entirely, so the two oxygens
     *      of every -C(=O)O- collapsed to one rank and were then split by
     *      atom id — the root of the old input-order dependence.
     *   2. Colour refinement (iterative) — repeatedly re-rank each
     *      atom by [own rank, sorted multiset of neighbour ranks] until the
     *      partition stops splitting.
     *   3. Iterative tie-breaking — while atoms remain tied, take the lowest
     *      tied rank class, promote one representative, and re-refine.
     *      Refinement after a promotion propagates the new distinction
     *      through the whole graph, so every constitutional difference that
     *      1-WL can see is resolved. A class that still cannot be split is a
     *      genuine automorphism orbit, where the choice of representative
     *      does not change the emitted SMILES — so the deterministic atom-id
     *      pick used to break it is canonical-safe.
     *
     * The ranking is constitutional: it does not consult tetrahedral @ / @@
     * or E/Z descriptors. For stereocentres that sit in a constitutional
     * automorphism orbit this leaves the DFS root — and therefore the
     * emitted @ / @@ — dependent on the atom-id tie-break; folding
     * tetrahedral parity into the seed to resolve that (the stereo-aware
     * ranking) is a deferred follow-up. See the KNOWN_ORBIT_DRIFT note in
     * tests/test_round_trip.js.
     */
    function canonicalRank(mol, atomIds) {
        var atomSet = {};
        atomIds.forEach(function(id) { atomSet[id] = true; });

        // In-component neighbour lists — computed once, reused every round.
        var nbr = {};
        atomIds.forEach(function(id) {
            var all = mol.getNeighbors(id);
            var inComp = [];
            for (var i = 0; i < all.length; i++) {
                if (atomSet[all[i]]) inComp.push(all[i]);
            }
            nbr[id] = inComp;
        });

        // Ring membership — one findRings pass shared by every atom.
        var ringAtom = {};
        var rings = mol.findRings();
        for (var r = 0; r < rings.length; r++) {
            var ra = rings[r].atoms;
            for (var k = 0; k < ra.length; k++) ringAtom[ra[k]] = true;
        }

        // Lexicographic compare of two numeric-array keys. Shorter keys are
        // padded with -Infinity so the order is a well-defined total preorder.
        function keyCmp(ka, kb) {
            var n = Math.max(ka.length, kb.length);
            for (var i = 0; i < n; i++) {
                var va = i < ka.length ? ka[i] : -Infinity;
                var vb = i < kb.length ? kb[i] : -Infinity;
                if (va !== vb) return va < vb ? -1 : 1;
            }
            return 0;
        }

        // Turn a per-atom key function into ranks. Equal keys collapse to one
        // rank; the rank is the sorted-position of the first atom with that
        // key, so ranks are stable but not dense. There is deliberately NO
        // atom-id fallback — genuine ties keep an equal rank, and only the
        // tie-breaking stage may split them.
        //
        // rankByKey runs O(atomIds.length) times across refinement and
        // tie-breaking, so its sort buffer is allocated once and refilled
        // each call (not re-sliced) — canonicalRank is synchronous and
        // non-re-entrant, so the single shared buffer is safe.
        var _rankScratch = atomIds.slice();
        function rankByKey(keyOf) {
            var ids = _rankScratch;
            for (var c = 0; c < ids.length; c++) ids[c] = atomIds[c];
            var keys = {};
            for (var i = 0; i < ids.length; i++) keys[ids[i]] = keyOf(ids[i]);
            ids.sort(function(a, b) { return keyCmp(keys[a], keys[b]); });
            var rank = {};
            var cur = 0;
            for (var i = 0; i < ids.length; i++) {
                if (i > 0 && keyCmp(keys[ids[i - 1]], keys[ids[i]]) !== 0) cur = i;
                rank[ids[i]] = cur;
            }
            return rank;
        }

        function distinctCount(rank) {
            var seen = {};
            for (var i = 0; i < atomIds.length; i++) seen[rank[atomIds[i]]] = true;
            return Object.keys(seen).length;
        }

        // v1.8.34: stereo-aware tie-break descriptor. When stage-3 iterative
        // tie-breaking reaches a constitutional automorphism orbit that
        // contains stereocentres, atom.chirality lets the orbit be split
        // deterministically instead of by atom id (which is parse-order-
        // dependent — the root of the KNOWN_ORBIT_DRIFT class). The
        // descriptor re-expresses the stored @/@@ token (SmilesParser
        // normalises it into the heavy-atom adjacency frame) against the
        // CURRENT rank order of the atom's neighbours: two orbit
        // stereocentres of opposite handedness get different descriptors,
        // while a genuine stereo-automorphism orbit (same handedness, e.g.
        // a C2-symmetric pair) gets equal descriptors and correctly falls
        // through to the atom-id tie-break.
        //
        // Returns 1 or 2 for a usable descriptor, or 0 when none applies (no
        // chirality, not a tetrahedral centre, too few substituents, or a
        // neighbour-rank tie makes the canonical neighbour order ambiguous).
        // A 0 sends the caller to the atom-id fallback — i.e. never worse
        // than the pre-v1.8.34 behaviour.
        function stereoDescriptor(id, rk) {
            var atom = mol.getAtom(id);
            if (!atom || !atom.chirality) return 0;
            var bonds = mol.getBondsOfAtom(id);
            var adjSeq = [];
            for (var i = 0; i < bonds.length; i++) {
                var nid = bonds[i].otherAtom(id);
                if (atomSet[nid]) adjSeq.push(nid);
            }
            var hCount = mol.calcHydrogens(id);
            if (hCount > 1) return 0;
            if (hCount === 1) adjSeq.push(-1);          // implicit H, last
            if (adjSeq.length < 3) return 0;
            // Every heavy neighbour must hold a distinct rank — otherwise the
            // canonical neighbour order is ambiguous and the descriptor would
            // not be input-order-invariant.
            var seenR = {};
            for (var k = 0; k < adjSeq.length; k++) {
                if (adjSeq[k] === -1) continue;
                var rv = rk[adjSeq[k]];
                if (seenR[rv]) return 0;
                seenR[rv] = true;
            }
            var canonSeq = adjSeq.slice().sort(function(a, c) {
                if (a === -1) return 1;
                if (c === -1) return -1;
                return rk[a] - rk[c];
            });
            var perm = [];
            for (var p = 0; p < adjSeq.length; p++) {
                perm.push(canonSeq.indexOf(adjSeq[p]));
            }
            var base = (atom.chirality === '@') ? 1 : 2;
            return (_permParity(perm) === 1) ? (3 - base) : base;
        }

        // --- stage 1: structural seed -------------------------------------
        var seed = {};
        for (var s = 0; s < atomIds.length; s++) {
            var sid = atomIds[s];
            var atom = mol.getAtom(sid);
            var hCount = atom.hydrogens >= 0 ? atom.hydrogens
                                             : mol.calcHydrogens(sid);
            var key = [
                nbr[sid].length,                     // degree within component
                ATOMIC_NUMBERS[atom.symbol] || 0,    // element
                atom.charge || 0,                    // formal charge
                mol.bondOrderSum(sid),               // explicit valence (bond orders)
                hCount,                              // attached hydrogens
                atom.isotope || 0,                   // isotope label
                ringAtom[sid] ? 1 : 0                // ring membership
            ];
            seed[sid] = key;
        }
        var rank = rankByKey(function(id) { return seed[id]; });

        // --- stage 2: iterative colour refinement ----------------
        function refine(rk) {
            var prev = distinctCount(rk);
            for (var iter = 0; iter < atomIds.length + 1; iter++) {
                var cur = rk;
                var next = rankByKey(function(id) {
                    var ns = [];
                    var list = nbr[id];
                    for (var j = 0; j < list.length; j++) ns.push(cur[list[j]]);
                    ns.sort(function(a, b) { return a - b; });
                    return [cur[id]].concat(ns);
                });
                var d = distinctCount(next);
                rk = next;
                if (d === prev) break;          // partition is stable
                prev = d;
            }
            return rk;
        }
        rank = refine(rank);

        // --- stage 3: iterative tie-breaking ------------------------------
        // Each pass promotes one victim: every atom keeps rank*2 except the
        // victim, which gets rank*2 - 1 — an odd value no other atom can
        // hold (its former classmates sit at rank*2, lower classes at
        // <= rank*2 - 2, higher classes at >= rank*2 + 2). So after the
        // re-rank the victim is alone in its class and the distinct-rank
        // count strictly rises by at least one every pass. The loop
        // therefore runs at most (atomIds.length - 1) times; `guard` is
        // purely defensive.
        var guard = atomIds.length + 1;
        while (distinctCount(rank) < atomIds.length && guard-- > 0) {
            var bucket = {};
            for (var b = 0; b < atomIds.length; b++) {
                var rv = rank[atomIds[b]];
                (bucket[rv] || (bucket[rv] = [])).push(atomIds[b]);
            }
            var rankValues = Object.keys(bucket).map(Number)
                .sort(function(a, c) { return a - c; });
            // v1.8.34: pick the lowest tied class, but PREFER the lowest one
            // that contains a stereocentre. Splitting a stereo orbit by its
            // descriptor and letting refinement propagate the distinction
            // resolves the non-stereo orbits downstream of it (e.g. the ring
            // CH2 pairs of cyclohexane-1,2-diol, each next to a differently-
            // configured stereocentre) BEFORE they would be atom-id-tie-
            // broken — atom id is parse-order-dependent, the stereo
            // descriptor is not. A non-stereo orbit that still survives to be
            // atom-id-broken afterwards is then a genuine automorphism orbit,
            // where the atom-id pick is canonical-safe.
            var tiedRank = -1, tiedStereoRank = -1;
            for (var t = 0; t < rankValues.length; t++) {
                var cls = bucket[rankValues[t]];
                if (cls.length <= 1) continue;
                if (tiedRank < 0) tiedRank = rankValues[t];
                if (tiedStereoRank < 0) {
                    for (var cs = 0; cs < cls.length; cs++) {
                        if (mol.getAtom(cls[cs]).chirality) {
                            tiedStereoRank = rankValues[t];
                            break;
                        }
                    }
                }
                if (tiedStereoRank >= 0) break;
            }
            if (tiedRank < 0) break;
            if (tiedStereoRank >= 0) tiedRank = tiedStereoRank;

            // Promote one representative of the chosen tied class. After
            // refinement the class is a constitutional automorphism orbit.
            // v1.8.34: if it contains stereocentres, split it by a canonical
            // stereo descriptor — that is input-order-invariant, unlike the
            // atom-id fallback, so orbit stereocentres of opposite handedness
            // (cyclohexane-1,2-diol etc.) no longer drift across a round-trip.
            // A class with no stereocentre, or a genuine same-handedness
            // stereo orbit, gets equal descriptors and falls through to the
            // atom-id pick exactly as before — so every non-orbit-stereo
            // molecule is byte-identical to v1.8.33. Promotion = halve the
            // rank scale, then nudge the representative one step below its
            // former classmates.
            var tiedClass = bucket[tiedRank];
            var sdesc = null;
            for (var hs = 0; hs < tiedClass.length; hs++) {
                if (mol.getAtom(tiedClass[hs]).chirality) {
                    sdesc = {};
                    for (var sd = 0; sd < tiedClass.length; sd++) {
                        sdesc[tiedClass[sd]] = stereoDescriptor(tiedClass[sd], rank);
                    }
                    break;
                }
            }
            var victim = tiedClass.slice().sort(function(a, c) {
                if (sdesc && sdesc[a] !== sdesc[c]) return sdesc[a] - sdesc[c];
                return a - c;                       // final fallback: atom id
            })[0];
            var promoted = {};
            for (var p = 0; p < atomIds.length; p++) {
                var pid = atomIds[p];
                promoted[pid] = rank[pid] * 2 - (pid === victim ? 1 : 0);
            }
            rank = refine(rankByKey(function(id) { return [promoted[id]]; }));
        }

        return rank;
    }

    // -----------------------------------------------------------------------
    // Write a single connected component
    // -----------------------------------------------------------------------

    function writeComponent(mol, atomIds, aromaticSet) {
        var atomSet = {};
        atomIds.forEach(function(id) { atomSet[id] = true; });

        var rank = canonicalRank(mol, atomIds);

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

        var bonds = mol.getBondsOfAtom(atom.id);

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
                            ringNeighbors.push(parseInt(aId, 10));
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
                            ringNeighbors.push(parseInt(aId, 10));
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

        // Compute @/@@ as the parity of the permutation from the STORED frame
        // (origSeq — the neighbour order atom.chirality is normalised into) to
        // the EMITTED frame (dfsSeq — the order the written SMILES is read back
        // in). Each frame places its implicit H where that frame actually puts
        // it (stored: H LAST; emitted: H immediately after the from-atom, where
        // [C@H] writes it), which is what makes parse -> write a round-trip
        // fixed point for a stereocentre with an implicit hydrogen.
        //
        // v2.4.17: origSeq reproduces SmilesParser's getNeighbors()+H-last frame
        // with NO from-atom hoist. The previous hoist inverted parity for ~26%
        // of chiral molecules (naproxen / limonene / tartaric); v1.8.34's partial
        // revert of it merely traded limonene for morphine. Matching the parser
        // frame EXACTLY fixes them all (pinned by test_v2_4_17_writer_stereo).
        // Both frames are constructed below with per-block notes.
        var origSeq = [];
        var dfsSeq = [];
        var hCount = mol.calcHydrogens(atom.id);

        // dfsSeq — the EMITTED frame: the from-atom is written first, then the
        // implicit H exactly where [C@H] puts it (immediately after the from-
        // atom), then ring closures, then branches.
        if (fromAtomId >= 0) dfsSeq.push(fromAtomId);
        for (var h = 0; h < hCount; h++) dfsSeq.push(-(h + 1));
        for (var di = 0; di < dfsOrder.length; di++) {
            if (dfsOrder[di] !== fromAtomId) dfsSeq.push(dfsOrder[di]);
        }

        // origSeq — the STORED frame: SmilesParser normalises atom.chirality
        // into mol.getNeighbors() adjacency order with the implicit H LAST and
        // NO from-atom hoist (see SmilesParser.js "Normalise chirality tokens to
        // the adjacency frame", target = getNeighbors + ['H']). Reproduce that
        // frame EXACTLY — from the same getNeighbors() source, from-atom left in
        // its natural adjacency position — so the origSeq -> dfsSeq permutation
        // is the true stored->emitted map. (The previous code hoisted the from-
        // atom to the front, which flipped parity whenever the from-atom was not
        // already first in adjacency order — the naproxen/limonene/tartaric
        // inversion; v1.8.34's partial revert traded limonene for morphine.)
        var nbrs = mol.getNeighbors(atom.id);
        for (var ni = 0; ni < nbrs.length; ni++) {
            if (atomSet[nbrs[ni]]) origSeq.push(nbrs[ni]);
        }
        for (var h2 = 0; h2 < hCount; h2++) origSeq.push(-(h2 + 1));

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
        // FIX: validate that arr is a permutation of [0..n-1] to avoid an
        // infinite loop if a caller accidentally supplies out-of-range values.
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

    // v1.8.28: perceiveAromaticity now lives in Molecule.js — the
    // single canonical Hückel implementation. SmilesWriter and
    // SmartsMatch previously each carried a copy; the two had
    // diverged (audit BUG-H1). Local alias keeps call sites terse.
    //
    // v2.4.17: the Hückel re-perception is deliberately conservative and does
    // NOT re-derive every aromatic system BIME can PARSE (tetrazoles, azulene,
    // some fused azines, the nucleobase lactams). For an atom the parser flagged
    // aromatic whose bonds are still stored single (aromatic SMILES input that
    // was never Kekulised), trust that flag — otherwise the writer emitted it as
    // saturated uppercase-single and the round-trip silently changed the
    // molecular formula (uracil -> dihydrouracil, azulene C10H8 -> C10H18, the
    // sartan tetrazoles gaining 5 H). calcHydrogens already treats these atoms
    // via AROMATIC_VALENCE, so the re-parse recovers the exact formula.
    function perceiveAromaticity(mol) {
        var set = Molecule.perceiveAromaticity(mol) || {};
        var atoms = mol.atoms;
        for (var i = 0; i < atoms.length; i++) {
            var a = atoms[i];
            if (!a.aromatic || set[a.id]) { continue; }
            var bonds = mol.getBondsOfAtom(a.id);
            var kekulised = false;
            for (var b = 0; b < bonds.length; b++) {
                if (bonds[b].type > Molecule.BOND_SINGLE) { kekulised = true; break; }
            }
            if (!kekulised) { set[a.id] = true; }
        }
        return set;
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

        // v2.4.17: an atom emitted AROMATIC whose real implicit-H count exceeds
        // what a bare aromatic symbol implies must pin that H with brackets.
        // Otherwise a pyrrole-type N-H re-aromatised from Kekulé input (indole,
        // serotonin, the purine 5-ring) is written as bare `n` and the reader
        // drops the NH (bare aromatic n = pyridine-type, 0 H). calcHydrogens
        // gives the true count; brackets then emit it as [nH].
        if (isAromatic && explicitH < 0) {
            var arTrueH = mol.calcHydrogens(atom.id);
            if (arTrueH > expectedH) { explicitH = arTrueH; }
        }

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
        // OpenSMILES %nn supports two-digit ring numbers (10-99). Numbers
        // beyond 99 are extremely rare in practice (would need >99 unclosed
        // rings simultaneously); we emit %nn anyway since the modulus
        // approach truncates and there is no portable extension syntax.
        return '%' + num;
    }

    // (v1.8.26: getDefaultValence moved to Molecule.js — see the
    // `var getDefaultValence = Molecule.getDefaultValence;` alias above.)

    // -----------------------------------------------------------------------
    // Export
    // -----------------------------------------------------------------------

    global.SmilesWriter = { write: write, writeNamed: writeNamed };

})(typeof window !== 'undefined' ? window : this);
