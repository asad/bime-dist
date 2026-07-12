/**
 * SmartsMatch.js — subgraph-isomorphism-based substructure matching for SMARTS queries
 *
 * Copyright (c) 2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2026 Syed Asad Rahman
 * Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 *
 * Implements a subgraph-isomorphism algorithm for matching SMARTS query
 * patterns against target molecules. Evaluates atom and bond constraints
 * including ring membership, degree, valence, hydrogen count, charge,
 * aromaticity, recursive SMARTS, and logical operators.
 *
 * Usage:
 *   var matches = SmartsMatch.match(targetMol, queryMol);
 *   var hasHit  = SmartsMatch.hasMatch(targetMol, queryMol);
 *   SmartsMatch.highlightMatches(targetMol, queryMol);
 */
(function(global) {
    'use strict';

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Find all substructure matches of query in target.
     * @param {Molecule} target — the molecule to search in
     * @param {Molecule} query  — the query molecule (from SmartsParser)
     * @returns {Array} array of mappings [{queryAtomId: targetAtomId, ...}, ...]
     */
    function match(target, query) {
        if (!target || !query) return [];
        if (query.atoms.length === 0) return [];
        if (target.atoms.length === 0) return [];

        // Precompute target properties
        var ctx = buildContext(target, query);

        // Run matcher
        var allMappings = [];
        subIso(ctx, {}, {}, 0, allMappings);

        return allMappings;
    }

    /**
     * Check whether target contains at least one match for query.
     * @param {Molecule} target
     * @param {Molecule} query
     * @returns {boolean}
     */
    function hasMatch(target, query) {
        if (!target || !query) return false;
        if (query.atoms.length === 0) return true;
        if (target.atoms.length === 0) return false;

        var ctx = buildContext(target, query);
        return subIsoFirst(ctx, {}, {}, 0);
    }

    /**
     * Find all matches and highlight matching atoms in target.
     * Sets atom.highlighted = true and atom.bgColor to teal for matches.
     * @param {Molecule} target
     * @param {Molecule} query
     * @returns {number} number of matches found
     */
    function highlightMatches(target, query) {
        // Clear previous highlights
        for (var i = 0; i < target.atoms.length; i++) {
            target.atoms[i].highlighted = false;
            target.atoms[i].bgColor = null;
        }
        for (var i = 0; i < target.bonds.length; i++) {
            target.bonds[i].highlighted = false;
            target.bonds[i].bgColor = null;
        }

        var mappings = match(target, query);

        for (var mi = 0; mi < mappings.length; mi++) {
            var mapping = mappings[mi];
            for (var qId in mapping) {
                if (!mapping.hasOwnProperty(qId)) continue;
                var tId = mapping[qId];
                var atom = target.getAtom(tId);
                if (atom) {
                    atom.highlighted = true;
                    atom.bgColor = '#0d9488';
                }
            }

            // Highlight bonds between matched atom pairs
            for (var qId1 in mapping) {
                if (!mapping.hasOwnProperty(qId1)) continue;
                for (var qId2 in mapping) {
                    if (!mapping.hasOwnProperty(qId2)) continue;
                    if (parseInt(qId1) >= parseInt(qId2)) continue;
                    var tBond = target.getBondBetween(mapping[qId1], mapping[qId2]);
                    if (tBond) {
                        tBond.highlighted = true;
                        tBond.bgColor = '#0d9488';
                    }
                }
            }
        }

        return mappings.length;
    }

    // -----------------------------------------------------------------------
    // Context building — precompute ring info, degrees, etc.
    // -----------------------------------------------------------------------

    function buildContext(target, query) {
        var ctx = {
            target: target,
            query: query,
            queryAtomIds: [],
            targetAtomIds: [],
            targetRings: null,
            targetRingMembership: {},  // atomId -> [ringSize, ...]
            targetRingCount: {},       // atomId -> number of rings containing it
            targetAromaticAtoms: {},   // atomId -> boolean
            maxMappings: 1000          // safety limit
        };

        // Query atom IDs in order
        for (var i = 0; i < query.atoms.length; i++) {
            ctx.queryAtomIds.push(query.atoms[i].id);
        }

        // Target atom IDs
        for (var i = 0; i < target.atoms.length; i++) {
            ctx.targetAtomIds.push(target.atoms[i].id);
        }

        // Compute ring info for target
        ctx.targetRings = target.findRings(8);
        for (var ri = 0; ri < ctx.targetRings.length; ri++) {
            var ring = ctx.targetRings[ri];
            for (var ai = 0; ai < ring.atoms.length; ai++) {
                var aid = ring.atoms[ai];
                if (!ctx.targetRingMembership[aid]) ctx.targetRingMembership[aid] = [];
                ctx.targetRingMembership[aid].push(ring.size);
                ctx.targetRingCount[aid] = (ctx.targetRingCount[aid] || 0) + 1;
            }
        }

        // Perceive aromaticity on target
        ctx.targetAromaticAtoms = perceiveAromaticity(target);

        return ctx;
    }

    // -----------------------------------------------------------------------
    // matcher — find all matches
    // -----------------------------------------------------------------------

    /**
     * Recursive matching.
     * @param {Object} ctx — context
     * @param {Object} queryMap — queryAtomId -> targetAtomId (current partial mapping)
     * @param {Object} targetUsed — targetAtomId -> true (atoms already used)
     * @param {number} depth — current depth (index into queryAtomIds)
     * @param {Array} allMappings — results array
     */
    function subIso(ctx, queryMap, targetUsed, depth, allMappings) {
        if (allMappings.length >= ctx.maxMappings) return;

        if (depth === ctx.queryAtomIds.length) {
            // Complete mapping found
            var mapping = {};
            for (var k in queryMap) mapping[k] = queryMap[k];
            allMappings.push(mapping);
            return;
        }

        var queryAtomId = ctx.queryAtomIds[depth];
        var queryAtom = ctx.query.getAtom(queryAtomId);

        // Try each target atom as a candidate
        for (var ti = 0; ti < ctx.targetAtomIds.length; ti++) {
            var targetAtomId = ctx.targetAtomIds[ti];
            if (targetUsed[targetAtomId]) continue;

            var targetAtom = ctx.target.getAtom(targetAtomId);

            // Check atom compatibility
            if (!atomMatchesQuery(targetAtom, queryAtom, ctx)) continue;

            // Check bond consistency: all query bonds to already-mapped query atoms
            // must have corresponding bonds in target
            if (!checkBondConsistency(ctx, queryMap, queryAtomId, targetAtomId)) continue;

            // Extend mapping
            queryMap[queryAtomId] = targetAtomId;
            targetUsed[targetAtomId] = true;

            subIso(ctx, queryMap, targetUsed, depth + 1, allMappings);

            // Backtrack
            delete queryMap[queryAtomId];
            delete targetUsed[targetAtomId];

            if (allMappings.length >= ctx.maxMappings) return;
        }
    }

    /**
     * Variant that returns true on first match (early termination).
     */
    function subIsoFirst(ctx, queryMap, targetUsed, depth) {
        if (depth === ctx.queryAtomIds.length) return true;

        var queryAtomId = ctx.queryAtomIds[depth];
        var queryAtom = ctx.query.getAtom(queryAtomId);

        for (var ti = 0; ti < ctx.targetAtomIds.length; ti++) {
            var targetAtomId = ctx.targetAtomIds[ti];
            if (targetUsed[targetAtomId]) continue;

            var targetAtom = ctx.target.getAtom(targetAtomId);
            if (!atomMatchesQuery(targetAtom, queryAtom, ctx)) continue;
            if (!checkBondConsistency(ctx, queryMap, queryAtomId, targetAtomId)) continue;

            queryMap[queryAtomId] = targetAtomId;
            targetUsed[targetAtomId] = true;

            if (subIsoFirst(ctx, queryMap, targetUsed, depth + 1)) return true;

            delete queryMap[queryAtomId];
            delete targetUsed[targetAtomId];
        }

        return false;
    }

    // -----------------------------------------------------------------------
    // Bond consistency check
    // -----------------------------------------------------------------------

    /**
     * Check that for every query bond from queryAtomId to an already-mapped
     * query neighbor, there exists a corresponding bond in the target between
     * targetAtomId and the mapped target neighbor, and the bond types are compatible.
     */
    function checkBondConsistency(ctx, queryMap, queryAtomId, targetAtomId) {
        var queryBonds = ctx.query.getBondsOfAtom(queryAtomId);

        for (var i = 0; i < queryBonds.length; i++) {
            var qBond = queryBonds[i];
            var qNeighbor = qBond.otherAtom(queryAtomId);

            // Only check already-mapped neighbors
            if (queryMap[qNeighbor] === undefined) continue;

            var tNeighbor = queryMap[qNeighbor];
            var tBond = ctx.target.getBondBetween(targetAtomId, tNeighbor);

            // Must have a bond in target
            if (!tBond) return false;

            // Check bond type compatibility
            if (!bondMatchesQuery(tBond, qBond, ctx)) return false;
        }

        return true;
    }

    // -----------------------------------------------------------------------
    // Atom constraint evaluation
    // -----------------------------------------------------------------------

    /**
     * Check whether a target atom satisfies all query constraints on a query atom.
     * @param {Atom} targetAtom
     * @param {Atom} queryAtom
     * @param {Object} ctx — matching context
     * @returns {boolean}
     */
    function atomMatchesQuery(targetAtom, queryAtom, ctx) {
        var constraints = queryAtom.queryConstraints;
        if (!constraints || constraints.length === 0) {
            // No constraints = wildcard, matches anything
            return true;
        }

        for (var i = 0; i < constraints.length; i++) {
            if (!evaluateConstraint(constraints[i], targetAtom, ctx)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Evaluate a single atom constraint against a target atom.
     */
    function evaluateConstraint(constraint, targetAtom, ctx) {
        var result = false;

        switch (constraint.type) {
            case 'wildcard':
                result = true;
                break;

            case 'atomicNum':
                var targetNum = SmartsParser.ATOMIC_NUMBERS[targetAtom.symbol] || 0;
                result = constraint.values.indexOf(targetNum) >= 0;
                break;

            case 'aromatic':
                result = !!(ctx.targetAromaticAtoms[targetAtom.id] || targetAtom.aromatic);
                break;

            case 'aliphatic':
                result = !(ctx.targetAromaticAtoms[targetAtom.id] || targetAtom.aromatic);
                break;

            case 'ring':
                var inRing = !!(ctx.targetRingMembership[targetAtom.id] &&
                    ctx.targetRingMembership[targetAtom.id].length > 0);
                result = (constraint.value === true) ? inRing : !inRing;
                break;

            case 'ringCount':
                var rc = ctx.targetRingCount[targetAtom.id] || 0;
                result = (rc === constraint.value);
                break;

            case 'ringSize':
                var rings = ctx.targetRingMembership[targetAtom.id] || [];
                result = false;
                for (var ri = 0; ri < rings.length; ri++) {
                    if (rings[ri] === constraint.value) { result = true; break; }
                }
                break;

            case 'degree':
                result = (ctx.target.degree(targetAtom.id) === constraint.value);
                break;

            case 'totalConnections':
                // X = explicit bonds + implicit hydrogens
                var explDeg = ctx.target.degree(targetAtom.id);
                var hCount = ctx.target.calcHydrogens(targetAtom.id);
                result = ((explDeg + hCount) === constraint.value);
                break;

            case 'valence':
                // v = total bond order including implicit H
                var bondOrd = ctx.target.bondOrderSum(targetAtom.id);
                var hCount = ctx.target.calcHydrogens(targetAtom.id);
                result = ((bondOrd + hCount) === constraint.value);
                break;

            case 'hCount':
                var hc = ctx.target.calcHydrogens(targetAtom.id);
                result = (hc === constraint.value);
                break;

            case 'charge':
                result = (targetAtom.charge === constraint.value);
                break;

            case 'isotope':
                // FIX: SMARTS isotope constraint (from leading-digit `[13C]`).
                // Matches when the target atom's isotope equals the constraint
                // value.  Target atom isotope of 0 means "natural mass" — never
                // matches a specific isotope query.
                result = ((targetAtom.isotope || 0) === constraint.value);
                break;

            case 'recursive':
                // Recursive SMARTS: atom must match the first atom of the inner pattern
                result = evaluateRecursive(constraint.smarts, targetAtom, ctx);
                break;

            case 'or':
                // OR of alternatives: at least one must match
                result = false;
                if (constraint.alternatives) {
                    for (var ai = 0; ai < constraint.alternatives.length; ai++) {
                        if (evaluateConstraint(constraint.alternatives[ai], targetAtom, ctx)) {
                            result = true;
                            break;
                        }
                    }
                }
                break;

            default:
                // Unknown constraint type — assume match
                result = true;
        }

        // Apply negation
        if (constraint.negate) result = !result;

        return result;
    }

    /**
     * Evaluate a recursive SMARTS pattern.
     * The target atom must match the first atom of the inner SMARTS,
     * and the rest of the pattern must also match.
     */
    function evaluateRecursive(smartsStr, targetAtom, ctx) {
        if (!smartsStr || typeof SmartsParser === 'undefined') return false;

        var innerQuery = SmartsParser.parse(smartsStr);
        if (!innerQuery || innerQuery.atoms.length === 0) return false;

        // Build a new context for the recursive match
        var innerCtx = buildContext(ctx.target, innerQuery);

        // The first query atom must map to targetAtom
        var firstQueryAtomId = innerCtx.queryAtomIds[0];
        var firstQueryAtom = innerQuery.getAtom(firstQueryAtomId);

        // Check if targetAtom can match the first query atom
        if (!atomMatchesQuery(targetAtom, firstQueryAtom, innerCtx)) return false;

        // If single-atom query, the atom match is sufficient
        if (innerQuery.atoms.length === 1) return true;

        // For multi-atom queries, try matching with first atom pinned
        var queryMap = {};
        var targetUsed = {};
        queryMap[firstQueryAtomId] = targetAtom.id;
        targetUsed[targetAtom.id] = true;

        return subIsoFirst(innerCtx, queryMap, targetUsed, 1);
    }

    // -----------------------------------------------------------------------
    // Bond constraint evaluation
    // -----------------------------------------------------------------------

    /**
     * Check whether a target bond matches the query bond type.
     * @param {Bond} targetBond
     * @param {Bond} queryBond
     * @param {Object} ctx
     * @returns {boolean}
     */
    function bondMatchesQuery(targetBond, queryBond, ctx) {
        var queryType = queryBond.queryType || 'default';

        // Check if target bond is aromatic
        var isAromatic = false;
        if (ctx) {
            var atom1Aro = (ctx.targetAromaticAtoms && ctx.targetAromaticAtoms[targetBond.atom1]);
            var atom2Aro = (ctx.targetAromaticAtoms && ctx.targetAromaticAtoms[targetBond.atom2]);
            // Also check atom.aromatic flag (set by SMILES parser for lowercase atoms)
            if (!atom1Aro) {
                var a1 = ctx.target.getAtom(targetBond.atom1);
                if (a1 && a1.aromatic) atom1Aro = true;
            }
            if (!atom2Aro) {
                var a2 = ctx.target.getAtom(targetBond.atom2);
                if (a2 && a2.aromatic) atom2Aro = true;
            }
            isAromatic = !!(atom1Aro && atom2Aro);
        }

        switch (queryType) {
            case 'any':
                return true;

            case 'single':
                return targetBond.type === Molecule.BOND_SINGLE && !isAromatic;

            case 'double':
                return targetBond.type === Molecule.BOND_DOUBLE;

            case 'triple':
                return targetBond.type === Molecule.BOND_TRIPLE;

            case 'aromatic':
                return isAromatic;

            case 'not_single':
                return !(targetBond.type === Molecule.BOND_SINGLE && !isAromatic);

            case 'not_double':
                return targetBond.type !== Molecule.BOND_DOUBLE;

            case 'not_triple':
                return targetBond.type !== Molecule.BOND_TRIPLE;

            case 'not_aromatic':
                return !isAromatic;

            case 'none':
                return false;

            case 'default':
                // Default bond in SMARTS: matches single (1) or aromatic (4) only
                if (isAromatic) { return true; }
                return targetBond.type === Molecule.BOND_SINGLE;

            default:
                return true;
        }
    }

    // -----------------------------------------------------------------------
    // Aromatic perception (duplicated from SmilesWriter for self-containment)
    // -----------------------------------------------------------------------

    // v1.8.28: perceiveAromaticity now lives in Molecule.js — the
    // single canonical Hückel implementation. This module's old
    // local copy mis-perceived aromaticity on Kekulé-form rings
    // (it checked only the bond to atom i+1); the canonical copy
    // inspects both ring bonds. This closes audit BUG-H1.
    var perceiveAromaticity = Molecule.perceiveAromaticity;

    // -----------------------------------------------------------------------
    // Convenience: parse SMARTS string and match in one step
    // -----------------------------------------------------------------------

    /**
     * Parse a SMARTS string and match against a target molecule.
     * @param {Molecule} target
     * @param {string} smartsStr
     * @returns {Array} array of atom mappings
     */
    function matchSmarts(target, smartsStr) {
        if (typeof SmartsParser === 'undefined') return [];
        var query = SmartsParser.parse(smartsStr);
        if (!query || query.atoms.length === 0) return [];
        return match(target, query);
    }

    /**
     * Parse a SMARTS string and highlight matches in target.
     * @param {Molecule} target
     * @param {string} smartsStr
     * @returns {number} number of matches
     */
    function highlightSmarts(target, smartsStr) {
        if (typeof SmartsParser === 'undefined') return 0;
        var query = SmartsParser.parse(smartsStr);
        if (!query || query.atoms.length === 0) return 0;
        return highlightMatches(target, query);
    }

    // -----------------------------------------------------------------------
    // Exports
    // -----------------------------------------------------------------------

    global.SmartsMatch = {
        match: match,
        hasMatch: hasMatch,
        highlightMatches: highlightMatches,
        matchSmarts: matchSmarts,
        highlightSmarts: highlightSmarts,
        atomMatchesQuery: atomMatchesQuery,
        bondMatchesQuery: bondMatchesQuery
    };

})(typeof window !== 'undefined' ? window : this);
