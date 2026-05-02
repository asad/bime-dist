/**
 * SmartsParser.js — SMARTS pattern string -> query Molecule with constraints
 *
 * Copyright (c) 2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2026 Syed Asad Rahman
 * Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 *
 * Parses SMARTS strings into query molecules where each atom and bond carries
 * constraint arrays for substructure matching. Supports:
 *   - Wildcard atoms: [*], [#6] (atomic number), [!#1] (NOT hydrogen)
 *   - Logical operators: [C,N] (OR), [C;R] (AND), [!C] (NOT)
 *   - Ring membership: [R], [R2], [r5]
 *   - Degree: [D2], [X3] (total connections)
 *   - Valence: [v3]
 *   - Hydrogen count: [H2]
 *   - Aromatic/aliphatic: [a], [A]
 *   - Charge: [+1], [-2]
 *   - Bond types: - = # : ~ !-
 *   - Recursive SMARTS: [$(*~[#7])]
 *   - Ring closures and branches
 */
(function(global) {
    'use strict';

    var BOND_LENGTH = Molecule.BOND_LENGTH;

    // Atomic numbers for element symbols
    var ATOMIC_NUMBERS = {
        'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
        'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,
        'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
        'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,
        'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
        'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,
        'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,
        'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,
        'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,
        'Pa':91,'U':92
    };

    // Reverse lookup: atomic number -> symbol
    var SYMBOL_FOR_NUM = {};
    for (var sym in ATOMIC_NUMBERS) {
        SYMBOL_FOR_NUM[ATOMIC_NUMBERS[sym]] = sym;
    }

    // Aromatic organic subset
    var AROMATIC_ORGANIC = { 'b':true, 'c':true, 'n':true, 'o':true, 'p':true, 's':true };

    // Organic subset (uppercase)
    var ORGANIC_SINGLE = { 'B':true, 'C':true, 'N':true, 'O':true, 'P':true, 'S':true, 'F':true, 'I':true };
    var ORGANIC_TWO = { 'Br':true, 'Cl':true };

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Parse a SMARTS pattern into a query Molecule.
     * Each atom will have a queryConstraints array.
     * Each bond will have a queryType string.
     * @param {string} smarts
     * @returns {Molecule} query molecule, or null on error
     */
    function parse(smarts) {
        if (!smarts || !smarts.trim()) return null;
        smarts = smarts.trim();

        var mol = new Molecule();
        mol.isQuery = true;
        mol.parseErrors = [];

        var result = parseFragment(smarts, mol);
        if (result.errors.length > 0) {
            mol.parseErrors = result.errors;
        }

        // Assign basic 2D coords via simple layout
        layoutQuery(mol);

        return mol;
    }

    // -----------------------------------------------------------------------
    // Fragment parser
    // -----------------------------------------------------------------------

    function parseFragment(smarts, mol) {
        var pos = 0;
        var len = smarts.length;
        var stack = [];
        var currentAtom = null;
        var ringOpenings = {};
        var pendingBondType = null;
        var errors = [];

        function peek() { return pos < len ? smarts[pos] : ''; }
        function next() { return pos < len ? smarts[pos++] : ''; }
        function isDigit(c) { return c >= '0' && c <= '9'; }

        while (pos < len) {
            var c = peek();

            // Branch open
            if (c === '(') {
                if (currentAtom === null) {
                    errors.push('Unexpected "(" at position ' + pos);
                }
                stack.push(currentAtom);
                next();
                continue;
            }

            // Branch close
            if (c === ')') {
                if (stack.length === 0) {
                    errors.push('Unbalanced ")" at position ' + pos);
                } else {
                    currentAtom = stack.pop();
                }
                next();
                continue;
            }

            // Bond symbols
            if (c === '~') {
                pendingBondType = 'any';
                next();
                continue;
            }
            if (c === '-') {
                next();
                pendingBondType = 'single';
                continue;
            }
            if (c === '=') {
                pendingBondType = 'double';
                next();
                continue;
            }
            if (c === '#') {
                // Could be triple bond or atomic number inside bracket
                // In SMARTS at top level, # is triple bond
                pendingBondType = 'triple';
                next();
                continue;
            }
            if (c === ':') {
                pendingBondType = 'aromatic';
                next();
                continue;
            }
            if (c === '/' || c === '\\') {
                // Geometric bond — treat as single for matching
                pendingBondType = 'single';
                next();
                continue;
            }
            if (c === '!') {
                // Negated bond: !- !:  etc. Peek before consuming so an invalid
                // follow-on character is reported without silently swallowing it.
                var bc = (pos + 1 < len) ? smarts[pos + 1] : '';
                if (bc === '-') { pendingBondType = 'not_single'; pos += 2; }
                else if (bc === '=') { pendingBondType = 'not_double'; pos += 2; }
                else if (bc === '#') { pendingBondType = 'not_triple'; pos += 2; }
                else if (bc === ':') { pendingBondType = 'not_aromatic'; pos += 2; }
                else if (bc === '~') { pendingBondType = 'none'; pos += 2; }
                else {
                    errors.push('Invalid negated bond at position ' + pos);
                    next(); // consume only the '!' so the next char is re-evaluated
                }
                continue;
            }

            // Ring closure
            if (c === '%' || isDigit(c)) {
                var ringStart = pos;
                var ringNum;
                if (c === '%') {
                    next();
                    var d1 = next();
                    var d2 = next();
                    if (!isDigit(d1) || !isDigit(d2)) {
                        errors.push('Invalid ring number at position ' + ringStart);
                        continue;
                    }
                    ringNum = parseInt(d1 + d2, 10);
                } else {
                    ringNum = parseInt(next(), 10);
                }

                if (currentAtom === null) {
                    errors.push('Ring closure with no current atom at position ' + ringStart);
                    continue;
                }

                if (ringOpenings[ringNum]) {
                    var opening = ringOpenings[ringNum];
                    var bt = pendingBondType || opening.bondType || 'default';
                    var bond = mol.addBond(opening.atomId, currentAtom, Molecule.BOND_SINGLE);
                    if (bond) bond.queryType = bt;
                    delete ringOpenings[ringNum];
                    pendingBondType = null;
                } else {
                    ringOpenings[ringNum] = {
                        atomId: currentAtom,
                        bondType: pendingBondType,
                        pos: ringStart
                    };
                    pendingBondType = null;
                }
                continue;
            }

            // Bracket atom (query atom specification)
            if (c === '[') {
                var bracketStart = pos;
                next(); // skip [
                var depth = 1;
                var bracketEnd = pos;
                while (bracketEnd < len && depth > 0) {
                    if (smarts[bracketEnd] === '[') depth++;
                    else if (smarts[bracketEnd] === ']') depth--;
                    if (depth > 0) bracketEnd++;
                }
                if (depth !== 0) {
                    errors.push('Unclosed bracket at position ' + bracketStart);
                    break;
                }
                var bracketStr = smarts.substring(pos, bracketEnd);
                pos = bracketEnd + 1;

                var constraints = parseBracketAtom(bracketStr, errors);
                var sym = constraints._symbol || '*';
                var atom = mol.addAtom(sym, 0, 0);
                atom.queryConstraints = constraints.constraints;
                atom.aromatic = constraints._aromatic || false;
                atom.isQuery = true;
                if (constraints._hCount >= 0) {
                    atom.hydrogens = constraints._hCount;
                }
                if (constraints._mapNumber > 0) {
                    atom.mapNumber = constraints._mapNumber;
                }
                // FIX: persist isotope on query atom so downstream consumers
                // (SmartsWriter, debug renderers) can see what mass was matched.
                if (constraints._isotope > 0) {
                    atom.isotope = constraints._isotope;
                }

                if (currentAtom !== null) {
                    var bt = pendingBondType || 'default';
                    var bond = mol.addBond(currentAtom, atom.id, Molecule.BOND_SINGLE);
                    if (bond) bond.queryType = bt;
                }
                currentAtom = atom.id;
                pendingBondType = null;
                continue;
            }

            // Organic/aromatic atom (unbracketed)
            if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
                var sym = null;
                var aromatic = false;
                var atomNum = 0;

                if (c >= 'a' && c <= 'z') {
                    if (AROMATIC_ORGANIC[c]) {
                        sym = c.toUpperCase();
                        aromatic = true;
                        atomNum = ATOMIC_NUMBERS[sym] || 0;
                        next();
                    } else {
                        errors.push('Unknown aromatic atom "' + c + '" at position ' + pos);
                        next();
                        continue;
                    }
                } else {
                    var twoChar = c + (pos + 1 < len ? smarts[pos + 1] : '');
                    if (ORGANIC_TWO[twoChar]) {
                        sym = twoChar;
                        atomNum = ATOMIC_NUMBERS[sym] || 0;
                        pos += 2;
                    } else if (ORGANIC_SINGLE[c]) {
                        sym = c;
                        atomNum = ATOMIC_NUMBERS[sym] || 0;
                        next();
                    } else {
                        errors.push('Unknown atom "' + c + '" at position ' + pos);
                        next();
                        continue;
                    }
                }

                var atom = mol.addAtom(sym, 0, 0);
                atom.aromatic = aromatic;
                atom.isQuery = true;
                // For unbracketed atoms, constraint is just matching the element
                atom.queryConstraints = [{
                    type: 'atomicNum',
                    values: [atomNum],
                    negate: false
                }];
                if (aromatic) {
                    atom.queryConstraints.push({ type: 'aromatic', value: true });
                } else {
                    atom.queryConstraints.push({ type: 'aliphatic', value: true });
                }

                if (currentAtom !== null) {
                    var bt = pendingBondType || 'default';
                    var bond = mol.addBond(currentAtom, atom.id, Molecule.BOND_SINGLE);
                    if (bond) bond.queryType = bt;
                }
                currentAtom = atom.id;
                pendingBondType = null;
                continue;
            }

            // Dot separator
            if (c === '.') {
                currentAtom = null;
                pendingBondType = null;
                next();
                continue;
            }

            errors.push('Unexpected character "' + c + '" at position ' + pos);
            next();
        }

        // Check unclosed rings
        var unclosed = Object.keys(ringOpenings);
        for (var i = 0; i < unclosed.length; i++) {
            errors.push('Unclosed ring ' + unclosed[i]);
        }

        if (stack.length > 0) {
            errors.push('Unbalanced parentheses: ' + stack.length + ' unclosed');
        }

        return { errors: errors };
    }

    // -----------------------------------------------------------------------
    // Bracket atom parser — parses constraint expressions inside [...]
    // -----------------------------------------------------------------------

    function parseBracketAtom(str, errors) {
        var constraints = [];
        var symbol = '*';
        var aromatic = false;
        var hCount = -1;
        var mapNumber = 0;
        var isotope = 0;
        var pos = 0;
        var len = str.length;

        function isDigit(c) { return c >= '0' && c <= '9'; }
        function peek() { return pos < len ? str[pos] : ''; }
        function next() { return pos < len ? str[pos++] : ''; }

        // FIX: SMARTS leading-digit isotope.  `[13C]` / `[2H]` etc. SMARTS
        // bracket atoms can begin with an isotope mass (per OpenSMARTS spec).
        // Previously the digits fell through `parseAtomPrimitive`'s unknown-
        // character path and were silently discarded, so `[13C]` matched any C.
        // We strip the leading digits here (before the per-primitive loop) and
        // emit an `isotope` constraint for the matcher.
        var isoStr = '';
        while (pos < len && isDigit(str[pos])) isoStr += str[pos++];
        if (isoStr) {
            isotope = parseInt(isoStr, 10);
            constraints.push({ type: 'isotope', value: isotope, negate: false });
        }

        // Handle recursive SMARTS: $(...) — parse the inner part
        function parseRecursive() {
            // We're positioned right after '$('
            var depth = 1;
            var start = pos;
            while (pos < len && depth > 0) {
                if (str[pos] === '(') depth++;
                else if (str[pos] === ')') depth--;
                if (depth > 0) pos++;
            }
            var inner = str.substring(start, pos);
            if (pos < len) pos++; // skip closing )
            return inner;
        }

        // Parse comma-separated OR groups and semicolon-separated AND groups
        // SMARTS bracket atoms: [expr1,expr2;expr3]
        // Comma = OR (lower precedence), Semicolon = AND (higher precedence)
        // We parse the entire bracket content as a series of primitive expressions

        while (pos < len) {
            var c = peek();

            // Skip semicolons and commas (we handle them as part of constraint logic)
            if (c === ';') {
                // AND — constraints are implicitly ANDed in our array
                next();
                continue;
            }
            if (c === ',') {
                // OR — we need to group with the previous constraint
                next();
                var orConstraints = parseAtomPrimitive(str, pos, errors);
                pos = orConstraints._pos;
                // Merge with last constraint if types match
                if (constraints.length > 0 && orConstraints.constraints.length > 0) {
                    var last = constraints[constraints.length - 1];
                    var orC = orConstraints.constraints[0];
                    if (last.type === 'atomicNum' && orC.type === 'atomicNum' && !last.negate && !orC.negate) {
                        // Merge OR: combine values
                        for (var vi = 0; vi < orC.values.length; vi++) {
                            if (last.values.indexOf(orC.values[vi]) < 0) {
                                last.values.push(orC.values[vi]);
                            }
                        }
                    } else {
                        // Can't merge, add as OR constraint
                        constraints.push({ type: 'or', alternatives: [last, orC] });
                        constraints.splice(constraints.length - 2, 1); // remove the old last
                    }
                }
                if (orConstraints._symbol && orConstraints._symbol !== '*') {
                    // Don't override primary symbol with OR alternatives
                }
                if (orConstraints._aromatic) aromatic = true;
                continue;
            }

            // Atom map number (:N) — must come at end of bracket, before closing ]
            if (c === ':') {
                next(); // skip ':'
                var mapStr = '';
                while (pos < len && str[pos] >= '0' && str[pos] <= '9') mapStr += str[pos++];
                if (mapStr) {
                    mapNumber = parseInt(mapStr, 10);
                }
                continue;
            }

            var prim = parseAtomPrimitive(str, pos, errors);
            pos = prim._pos;
            for (var ci = 0; ci < prim.constraints.length; ci++) {
                constraints.push(prim.constraints[ci]);
            }
            if (prim._symbol && prim._symbol !== '*') symbol = prim._symbol;
            if (prim._aromatic) aromatic = true;
            if (prim._hCount >= 0) hCount = prim._hCount;
        }

        return {
            constraints: constraints,
            _symbol: symbol,
            _aromatic: aromatic,
            _hCount: hCount,
            _mapNumber: mapNumber,
            _isotope: isotope
        };
    }

    /**
     * Parse a single atom primitive expression from position pos in str.
     * Returns { constraints: [...], _pos: newPos, _symbol, _aromatic, _hCount }
     */
    function parseAtomPrimitive(str, pos, errors) {
        var constraints = [];
        var symbol = '*';
        var aromatic = false;
        var hCount = -1;
        var len = str.length;

        function isDigit(c) { return c >= '0' && c <= '9'; }

        // Consume one primitive
        if (pos >= len) {
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        var c = str[pos];
        var negate = false;

        // NOT operator
        if (c === '!') {
            negate = true;
            pos++;
            c = pos < len ? str[pos] : '';
        }

        // Wildcard
        if (c === '*') {
            pos++;
            constraints.push({ type: 'wildcard', negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: '*', _aromatic: aromatic, _hCount: hCount };
        }

        // Atomic number: #N
        if (c === '#') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            var num = parseInt(numStr, 10);
            constraints.push({ type: 'atomicNum', values: [num], negate: negate });
            symbol = SYMBOL_FOR_NUM[num] || '*';
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Ring membership: R or Rn
        if (c === 'R') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            if (numStr) {
                constraints.push({ type: 'ringCount', value: parseInt(numStr, 10), negate: negate });
            } else {
                constraints.push({ type: 'ring', value: true, negate: negate });
            }
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Ring size: rN
        if (c === 'r') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            if (numStr) {
                constraints.push({ type: 'ringSize', value: parseInt(numStr, 10), negate: negate });
            } else {
                constraints.push({ type: 'ring', value: true, negate: negate });
            }
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Degree: Dn
        if (c === 'D') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            constraints.push({ type: 'degree', value: numStr ? parseInt(numStr, 10) : 1, negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Total connections: Xn
        if (c === 'X') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            constraints.push({ type: 'totalConnections', value: numStr ? parseInt(numStr, 10) : 1, negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Valence: vN
        if (c === 'v') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            constraints.push({ type: 'valence', value: numStr ? parseInt(numStr, 10) : 1, negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Hydrogen count: Hn
        if (c === 'H') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            var hVal = numStr ? parseInt(numStr, 10) : 1;
            constraints.push({ type: 'hCount', value: hVal, negate: negate });
            hCount = hVal;
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Aromatic: a
        if (c === 'a') {
            pos++;
            aromatic = !negate;
            constraints.push({ type: 'aromatic', value: true, negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Aliphatic: A
        if (c === 'A') {
            pos++;
            constraints.push({ type: 'aliphatic', value: true, negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Charge: +n or -n
        if (c === '+') {
            pos++;
            var numStr = '';
            while (pos < len && isDigit(str[pos])) numStr += str[pos++];
            // Support '++' / '+++' style multi-charge
            var multiCharge = 1;
            while (pos < len && str[pos] === '+') { multiCharge++; pos++; }
            var charge = numStr ? parseInt(numStr, 10) : multiCharge;
            constraints.push({ type: 'charge', value: charge, negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }
        // Negative charge: '-' followed by digits, end-of-bracket, or another
        // bracket-level separator (';', ',', ':', '+'/'-' for multi-charges).
        // Previously only digits/']' were accepted, so '[C-;R]' silently lost
        // the charge constraint.
        if (c === '-') {
            var nextCh = (pos + 1 < len) ? str[pos + 1] : '';
            if (pos + 1 >= len || isDigit(nextCh) || nextCh === ']' ||
                nextCh === ';' || nextCh === ',' || nextCh === ':' ||
                nextCh === '+' || nextCh === '-') {
                pos++;
                var numStr = '';
                while (pos < len && isDigit(str[pos])) numStr += str[pos++];
                // Support '--' / '---' style multi-charge (rare but legal SMARTS)
                var multiCharge = -1;
                while (pos < len && str[pos] === '-') { multiCharge--; pos++; }
                var charge = numStr ? -parseInt(numStr, 10) : multiCharge;
                constraints.push({ type: 'charge', value: charge, negate: negate });
                return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
            }
        }

        // Recursive SMARTS: $(...)
        if (c === '$') {
            pos++;
            if (pos < len && str[pos] === '(') {
                pos++; // skip (
                var depth = 1;
                var start = pos;
                while (pos < len && depth > 0) {
                    if (str[pos] === '(') depth++;
                    else if (str[pos] === ')') depth--;
                    if (depth > 0) pos++;
                }
                var inner = str.substring(start, pos);
                if (pos < len) pos++; // skip )
                constraints.push({ type: 'recursive', smarts: inner, negate: negate });
            }
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // Element symbol (aromatic lowercase or uppercase)
        if (c >= 'a' && c <= 'z') {
            // Aromatic element
            var twoLower = str.substring(pos, pos + 2);
            if (twoLower === 'se' || twoLower === 'as' || twoLower === 'te') {
                symbol = twoLower.charAt(0).toUpperCase() + twoLower.charAt(1);
                aromatic = !negate;
                pos += 2;
            } else {
                symbol = c.toUpperCase();
                aromatic = !negate;
                pos++;
            }
            var atomNum = ATOMIC_NUMBERS[symbol] || 0;
            constraints.push({ type: 'atomicNum', values: [atomNum], negate: negate });
            if (!negate) constraints.push({ type: 'aromatic', value: true });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        if (c >= 'A' && c <= 'Z') {
            // Regular element
            var twoChar = str.substring(pos, pos + 2);
            // Check for two-char elements
            if (pos + 1 < len && str[pos + 1] >= 'a' && str[pos + 1] <= 'z') {
                if (ATOMIC_NUMBERS[twoChar]) {
                    symbol = twoChar;
                    pos += 2;
                } else {
                    symbol = c;
                    pos++;
                }
            } else {
                symbol = c;
                pos++;
            }
            var atomNum = ATOMIC_NUMBERS[symbol] || 0;
            constraints.push({ type: 'atomicNum', values: [atomNum], negate: negate });
            return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
        }

        // If we get here, unknown character — skip it
        pos++;
        return { constraints: constraints, _pos: pos, _symbol: symbol, _aromatic: aromatic, _hCount: hCount };
    }

    // -----------------------------------------------------------------------
    // Simple query layout (BFS tree)
    // -----------------------------------------------------------------------

    function layoutQuery(mol) {
        if (mol.atoms.length === 0) return;

        var placed = {};
        var queue = [];

        mol.atoms[0].x = 0;
        mol.atoms[0].y = 0;
        placed[mol.atoms[0].id] = true;
        queue.push(mol.atoms[0].id);

        while (queue.length > 0) {
            var atomId = queue.shift();
            var atom = mol.getAtom(atomId);
            var neighbors = mol.getNeighbors(atomId);
            var placedNeighbors = [];
            var unplaced = [];
            for (var i = 0; i < neighbors.length; i++) {
                if (placed[neighbors[i]]) placedNeighbors.push(neighbors[i]);
                else unplaced.push(neighbors[i]);
            }

            if (unplaced.length === 0) continue;

            var startAngle;
            if (placedNeighbors.length > 0) {
                var avgAngle = 0;
                for (var j = 0; j < placedNeighbors.length; j++) {
                    var n = mol.getAtom(placedNeighbors[j]);
                    avgAngle += Math.atan2(n.y - atom.y, n.x - atom.x);
                }
                avgAngle /= placedNeighbors.length;
                startAngle = avgAngle + Math.PI;
            } else {
                startAngle = 0;
            }

            var spread = Math.PI / 3;
            var totalSpread = (unplaced.length - 1) * spread;
            var firstAngle = startAngle - totalSpread / 2;

            for (var i = 0; i < unplaced.length; i++) {
                var nAtom = mol.getAtom(unplaced[i]);
                var angle = firstAngle + i * spread;
                nAtom.x = atom.x + BOND_LENGTH * Math.cos(angle);
                nAtom.y = atom.y + BOND_LENGTH * Math.sin(angle);
                placed[unplaced[i]] = true;
                queue.push(unplaced[i]);
            }
        }

        // Center
        if (mol.atoms.length > 0) {
            var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
            for (var i = 0; i < mol.atoms.length; i++) {
                var a = mol.atoms[i];
                if (a.x < minX) minX = a.x;
                if (a.y < minY) minY = a.y;
                if (a.x > maxX) maxX = a.x;
                if (a.y > maxY) maxY = a.y;
            }
            var cx = (minX + maxX) / 2;
            var cy = (minY + maxY) / 2;
            for (var i = 0; i < mol.atoms.length; i++) {
                mol.atoms[i].x -= cx;
                mol.atoms[i].y -= cy;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Exports
    // -----------------------------------------------------------------------

    var SmartsParser = {
        parse: parse,
        ATOMIC_NUMBERS: ATOMIC_NUMBERS,
        SYMBOL_FOR_NUM: SYMBOL_FOR_NUM
    };

    global.SmartsParser = SmartsParser;

})(typeof window !== 'undefined' ? window : this);
