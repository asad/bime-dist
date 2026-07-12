/**
 * SmartsWriter.js — Query Molecule -> SMARTS string
 *
 * Copyright (c) 2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2026 Syed Asad Rahman
 * Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 *
 * Converts a query molecule (produced by SmartsParser) back to a SMARTS string.
 * Uses DFS traversal with bracket notation for constraint atoms.
 */
(function(global) {
    'use strict';

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Write a SMARTS string from a query Molecule.
     * @param {Molecule} mol — query molecule with queryConstraints on atoms
     * @returns {string} SMARTS pattern
     */
    function write(mol) {
        if (!mol || mol.atoms.length === 0) return '';

        // Multi-component
        var components = mol.getComponents();
        var parts = [];
        for (var i = 0; i < components.length; i++) {
            parts.push(writeComponent(mol, components[i]));
        }
        return parts.join('.');
    }

    // -----------------------------------------------------------------------
    // Write a single connected component
    // -----------------------------------------------------------------------

    function writeComponent(mol, atomIds) {
        var atomSet = {};
        for (var i = 0; i < atomIds.length; i++) atomSet[atomIds[i]] = true;

        // Find ring closures via DFS
        var visited = {};
        var ringClosures = {};
        var ringBonds = {};
        var nextRing = 1;

        // First pass: detect ring bonds
        function detectRings(atomId, parentBondId) {
            visited[atomId] = true;
            var bonds = mol.getBondsOfAtom(atomId);
            for (var i = 0; i < bonds.length; i++) {
                var bond = bonds[i];
                if (bond.id === parentBondId) continue;
                var neighbor = bond.otherAtom(atomId);
                if (!atomSet[neighbor]) continue;
                if (visited[neighbor]) {
                    if (!ringBonds[bond.id]) {
                        var num = nextRing++;
                        if (!ringClosures[neighbor]) ringClosures[neighbor] = [];
                        ringClosures[neighbor].push({ num: num, bond: bond });
                        if (!ringClosures[atomId]) ringClosures[atomId] = [];
                        ringClosures[atomId].push({ num: num, bond: bond });
                        ringBonds[bond.id] = true;
                    }
                } else {
                    detectRings(neighbor, bond.id);
                }
            }
        }
        detectRings(atomIds[0], -1);

        // Second pass: DFS to build SMARTS string
        var visited2 = {};
        var result = '';

        function dfs(atomId, fromBondId) {
            visited2[atomId] = true;
            var atom = mol.getAtom(atomId);
            result += atomToSmarts(atom);

            // Ring closure digits
            if (ringClosures[atomId]) {
                for (var ri = 0; ri < ringClosures[atomId].length; ri++) {
                    var rc = ringClosures[atomId][ri];
                    if (visited2[rc.bond.otherAtom(atomId)] || rc.bond.id === fromBondId) {
                        // only emit from the opening side (first visit)
                    }
                    var bondStr = bondToSmarts(rc.bond);
                    var numStr = rc.num < 10 ? '' + rc.num : '%' + rc.num;
                    result += bondStr + numStr;
                }
            }

            // Gather unvisited neighbors
            var bonds = mol.getBondsOfAtom(atomId);
            var branches = [];
            for (var i = 0; i < bonds.length; i++) {
                var bond = bonds[i];
                if (bond.id === fromBondId) continue;
                var neighbor = bond.otherAtom(atomId);
                if (!atomSet[neighbor]) continue;
                if (visited2[neighbor]) continue;
                if (ringBonds[bond.id]) continue;
                branches.push({ bond: bond, neighbor: neighbor });
            }

            for (var i = 0; i < branches.length; i++) {
                var br = branches[i];
                var bondStr = bondToSmarts(br.bond);
                if (i < branches.length - 1) {
                    result += '(' + bondStr;
                    dfs(br.neighbor, br.bond.id);
                    result += ')';
                } else {
                    result += bondStr;
                    dfs(br.neighbor, br.bond.id);
                }
            }
        }

        dfs(atomIds[0], -1);
        return result;
    }

    // -----------------------------------------------------------------------
    // Atom to SMARTS string
    // -----------------------------------------------------------------------

    function atomToSmarts(atom) {
        var constraints = atom.queryConstraints;
        if (!constraints || constraints.length === 0) {
            // No constraints: wildcard
            return '[*]';
        }

        // Check if this is a simple element with no extra constraints
        if (canWriteUnbracketed(atom, constraints)) {
            var sym = atom.symbol;
            if (atom.aromatic) return sym.toLowerCase();
            return sym;
        }

        // Build bracket representation
        // FIX: isotope must be the first token inside the bracket WITHOUT a
        // ';' separator (SMARTS requires the bare leading mass), so split it
        // out from the rest of the constraint list.
        var isotopePrefix = '';
        var parts = [];
        for (var i = 0; i < constraints.length; i++) {
            var c = constraints[i];
            if (c.type === 'isotope' && !c.negate) {
                isotopePrefix = '' + c.value;
                continue;
            }
            var s = constraintToString(c);
            if (s) parts.push(s);
        }

        if (parts.length === 0 && !isotopePrefix) {
            if (atom.mapNumber > 0) return '[*:' + atom.mapNumber + ']';
            return '[*]';
        }
        var result = '[' + isotopePrefix + parts.join(';');
        if (atom.mapNumber > 0) result += ':' + atom.mapNumber;
        result += ']';
        return result;
    }

    /**
     * Check if an atom can be written without brackets (simple organic atom).
     */
    function canWriteUnbracketed(atom, constraints) {
        if (atom.mapNumber > 0) return false;
        // FIX: any isotope constraint forces brackets — `13c` is not legal SMARTS.
        if (atom.isotope > 0) return false;
        for (var ix = 0; ix < constraints.length; ix++) {
            if (constraints[ix].type === 'isotope') return false;
        }
        if (constraints.length > 2) return false;

        var organicSymbols = { 'B':1,'C':1,'N':1,'O':1,'P':1,'S':1,'F':1,'Cl':1,'Br':1,'I':1 };
        var aromaticOrganic = { 'B':1,'C':1,'N':1,'O':1,'P':1,'S':1 };

        // Must be just an atomicNum constraint (and optionally aromatic/aliphatic)
        var hasOnlyElement = true;
        var elementNum = -1;
        for (var i = 0; i < constraints.length; i++) {
            var c = constraints[i];
            if (c.type === 'atomicNum' && !c.negate && c.values.length === 1) {
                elementNum = c.values[0];
            } else if ((c.type === 'aromatic' || c.type === 'aliphatic') && !c.negate) {
                // OK — only safe to drop into unbracketed form when not negated.
                // A negated aromatic/aliphatic flag must be preserved as [!a]/[!A].
            } else {
                hasOnlyElement = false;
            }
        }

        if (!hasOnlyElement || elementNum < 0) return false;

        var sym = SmartsParser.SYMBOL_FOR_NUM[elementNum];
        if (!sym) return false;

        if (atom.aromatic) {
            return !!aromaticOrganic[sym];
        }
        return !!organicSymbols[sym];
    }

    /**
     * Convert a single constraint to SMARTS notation string.
     */
    function constraintToString(c) {
        var prefix = c.negate ? '!' : '';

        switch (c.type) {
            case 'wildcard':
                return prefix + '*';

            case 'atomicNum':
                if (c.values.length === 1) {
                    var sym = SmartsParser.SYMBOL_FOR_NUM[c.values[0]];
                    if (sym) return prefix + '#' + c.values[0];
                    return prefix + '#' + c.values[0];
                }
                // OR of multiple
                var parts = [];
                for (var i = 0; i < c.values.length; i++) {
                    parts.push(prefix + '#' + c.values[i]);
                }
                return parts.join(',');

            case 'aromatic':
                return prefix + 'a';

            case 'aliphatic':
                return prefix + 'A';

            case 'ring':
                return prefix + 'R';

            case 'ringCount':
                return prefix + 'R' + c.value;

            case 'ringSize':
                return prefix + 'r' + c.value;

            case 'degree':
                return prefix + 'D' + c.value;

            case 'totalConnections':
                return prefix + 'X' + c.value;

            case 'valence':
                return prefix + 'v' + c.value;

            case 'hCount':
                return prefix + 'H' + c.value;

            case 'charge':
                if (c.value > 0) return prefix + '+' + c.value;
                if (c.value < 0) return prefix + c.value;
                return prefix + '+0';

            case 'isotope':
                // FIX: SMARTS isotope round-trip.  SMARTS encodes isotope as a
                // leading mass number inside brackets (e.g. `[13C]`).  We emit
                // the digits inline; the bracket-wrap is added by atomToSmarts.
                // (Negation of an isotope is uncommon in SMARTS but supported
                // here for symmetry — `!13` would mean "not mass 13".)
                return prefix + c.value;

            case 'recursive':
                return prefix + '$(' + c.smarts + ')';

            case 'or':
                if (c.alternatives) {
                    var parts = [];
                    for (var i = 0; i < c.alternatives.length; i++) {
                        parts.push(constraintToString(c.alternatives[i]));
                    }
                    return parts.join(',');
                }
                return '';

            default:
                return '';
        }
    }

    // -----------------------------------------------------------------------
    // Bond to SMARTS string
    // -----------------------------------------------------------------------

    function bondToSmarts(bond) {
        var qt = bond.queryType || 'default';
        switch (qt) {
            case 'any':           return '~';
            case 'single':        return '-';
            case 'double':        return '=';
            case 'triple':        return '#';
            case 'aromatic':      return ':';
            case 'not_single':    return '!-';
            case 'not_double':    return '!=';
            case 'not_triple':    return '!#';
            case 'not_aromatic':  return '!:';
            case 'default':       return '';
            default:              return '';
        }
    }

    // -----------------------------------------------------------------------
    // Exports
    // -----------------------------------------------------------------------

    global.SmartsWriter = { write: write };

})(typeof window !== 'undefined' ? window : this);
