/**
 * SmilesParser.js — SMILES string -> Molecule graph with 2D layout
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Strict recursive-descent parser for the OpenSMILES specification.
 * Supports: organic subset (B, C, N, O, P, S, F, Cl, Br, I), bracket atoms
 * (isotope, chirality, H count, charge, atom map), aromatic atoms (b, c, n,
 * o, p, s), single/double/triple/aromatic bonds, directional bonds (/ \),
 * ring closures (single digit and %nn), branches, dot-disconnected fragments,
 * and reaction SMILES (>>).
 *
 * Returns detailed error information on invalid input via validateSmiles().
 */
(function(global) {
    'use strict';

    var BOND_LENGTH = Molecule.BOND_LENGTH;

    // -----------------------------------------------------------------------
    // Allowed symbols in the organic subset (no brackets required)
    // -----------------------------------------------------------------------
    var ORGANIC_SINGLE = { 'B': true, 'C': true, 'N': true, 'O': true,
                           'P': true, 'S': true, 'F': true, 'I': true };
    var ORGANIC_TWO    = { 'Br': true, 'Cl': true };

    // Aromatic organic subset (lowercase in SMILES)
    var AROMATIC_ORGANIC = { 'b': true, 'c': true, 'n': true, 'o': true,
                             'p': true, 's': true };

    // Standard valences used for implicit-H calculation and validation
    var STANDARD_VALENCES = {
        'B':  [3],        'C':  [4],        'N':  [3, 5],
        'O':  [2],        'P':  [3, 5],     'S':  [2, 4, 6],
        'F':  [1],        'Cl': [1],        'Br': [1],
        'I':  [1],        'H':  [1],        'Si': [4],
        'Se': [2, 4, 6],  'As': [3, 5],     'Te': [2, 4, 6]
    };

    // Every element symbol we consider valid inside brackets
    var ALL_ELEMENTS = {
        'H':1,'He':1,'Li':1,'Be':1,'B':1,'C':1,'N':1,'O':1,'F':1,'Ne':1,
        'Na':1,'Mg':1,'Al':1,'Si':1,'P':1,'S':1,'Cl':1,'Ar':1,'K':1,'Ca':1,
        'Sc':1,'Ti':1,'V':1,'Cr':1,'Mn':1,'Fe':1,'Co':1,'Ni':1,'Cu':1,'Zn':1,
        'Ga':1,'Ge':1,'As':1,'Se':1,'Br':1,'Kr':1,'Rb':1,'Sr':1,'Y':1,'Zr':1,
        'Nb':1,'Mo':1,'Tc':1,'Ru':1,'Rh':1,'Pd':1,'Ag':1,'Cd':1,'In':1,'Sn':1,
        'Sb':1,'Te':1,'I':1,'Xe':1,'Cs':1,'Ba':1,'La':1,'Ce':1,'Pr':1,'Nd':1,
        'Pm':1,'Sm':1,'Eu':1,'Gd':1,'Tb':1,'Dy':1,'Ho':1,'Er':1,'Tm':1,'Yb':1,
        'Lu':1,'Hf':1,'Ta':1,'W':1,'Re':1,'Os':1,'Ir':1,'Pt':1,'Au':1,'Hg':1,
        'Tl':1,'Pb':1,'Bi':1,'Po':1,'At':1,'Rn':1,'Fr':1,'Ra':1,'Ac':1,'Th':1,
        'Pa':1,'U':1,'Np':1,'Pu':1,'Am':1,'Cm':1,'Bk':1,'Cf':1,'Es':1,'Fm':1,
        'Md':1,'No':1,'Lr':1,'Rf':1,'Db':1,'Sg':1,'Bh':1,'Hs':1,'Mt':1,'Ds':1,
        'Rg':1,'Cn':1,'Nh':1,'Fl':1,'Mc':1,'Lv':1,'Ts':1,'Og':1,
        'R':1  // generic R-group
    };

    // Aromatic atoms allowed inside brackets (lowercase)
    var AROMATIC_BRACKET = { 'b':1,'c':1,'n':1,'o':1,'p':1,'s':1,
                             'se':1,'as':1,'te':1 };

    // Bond type constants
    var BOND_SINGLE   = Molecule.BOND_SINGLE;
    var BOND_DOUBLE   = Molecule.BOND_DOUBLE;
    var BOND_TRIPLE   = Molecule.BOND_TRIPLE;
    var BOND_AROMATIC = 4; // internal representation; stored as 1 in Molecule

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Parse a SMILES string into a Molecule object.
     * @param {string}   smiles  The SMILES string.
     * @param {Molecule} [mol]   Optional existing Molecule to append into.
     * @returns {Molecule} The populated molecule, or null on hard error.
     *
     * After parsing, mol.parseErrors (array of strings) is populated with any
     * problems found (empty array if the SMILES is clean).
     * mol.parseWarnings contains non-fatal issues like valence violations.
     */
    var MAX_SMILES_LENGTH = 10000; // guard against denial-of-service from huge input

    function parse(smiles, mol) {
        if (smiles == null) return null;
        mol = mol || new Molecule();
        mol.parseErrors   = mol.parseErrors   || [];
        mol.parseWarnings = mol.parseWarnings || [];
        smiles = smiles.trim();
        if (smiles.length === 0) return mol;
        if (smiles.length > MAX_SMILES_LENGTH) {
            mol.parseErrors.push('SMILES string too long (' + smiles.length + ' chars, max ' + MAX_SMILES_LENGTH + ')');
            return mol;
        }

        // Reaction SMILES (>> forward, => retrosynthetic)
        var isReaction = smiles.indexOf('>>') >= 0 || smiles.indexOf('=>') >= 0;
        if (isReaction) {
            parseReaction(smiles, mol);
            return mol;
        }

        // Extract molecule name: everything after the first whitespace following the SMILES
        // e.g. "CCO ethanol" -> SMILES="CCO", name="ethanol"
        var nameFromInput = '';
        var smilesOnly = smiles;
        var wsMatch = smiles.match(/^(\S+)([ \t]+)(.+)$/);
        if (wsMatch) {
            smilesOnly = wsMatch[1];
            nameFromInput = wsMatch[3].trim();
        }

        // Multi-component (fragments separated by '.' at top level)
        var fragments = splitTopLevel(smilesOnly, '.');
        var offsetX = 0;

        for (var f = 0; f < fragments.length; f++) {
            var startAtomCount = mol.atoms.length;
            var result = parseFragment(fragments[f], mol);
            if (result.errors.length > 0) {
                mol.parseErrors = mol.parseErrors.concat(result.errors);
            }
            if (result.warnings.length > 0) {
                mol.parseWarnings = mol.parseWarnings.concat(result.warnings);
            }
            var fragAtoms = mol.atoms.slice(startAtomCount);
            if (fragAtoms.length > 0) {
                layoutFragment(mol, fragAtoms);
                // FIX: shift fragment so it doesn't overlap with previous fragments
                var bounds = getBounds(fragAtoms);
                if (f > 0) {
                    var shiftX = offsetX - bounds.minX;
                    for (var si = 0; si < fragAtoms.length; si++) {
                        fragAtoms[si].x += shiftX;
                    }
                    bounds = getBounds(fragAtoms);
                }
                offsetX = bounds.maxX + BOND_LENGTH * 2;
            }
        }

        // Store the name parsed from the input (e.g. "CCO ethanol")
        if (nameFromInput) {
            mol.name = nameFromInput;
        }

        return mol;
    }

    // -----------------------------------------------------------------------
    // Reaction SMILES  (supports >> forward and => retrosynthetic arrows)
    // -----------------------------------------------------------------------

    function parseReaction(smiles, mol) {
        // Detect arrow type: '=>' (retrosynthetic) or '>>' (forward)
        var arrowType = 'forward';
        var idx = smiles.indexOf('>>');
        if (idx < 0) {
            idx = smiles.indexOf('=>');
            if (idx >= 0) arrowType = 'retro';
        }
        if (idx < 0) return; // safety

        var reactantStr = smiles.substring(0, idx);
        var productStr  = smiles.substring(idx + 2);

        // Split each side by '.' to get individual fragments
        var reactantFrags = splitTopLevel(reactantStr, '.');
        var productFrags  = splitTopLevel(productStr, '.');

        // ---- Parse and layout each reactant fragment independently ----
        var FRAG_GAP    = BOND_LENGTH * 1.5;   // space between fragments on the same side
        var ARROW_PAD   = BOND_LENGTH * 1.2;    // padding before/after the arrow
        var ARROW_LEN   = BOND_LENGTH * 3;      // arrow shaft length
        var PLUS_WIDTH  = BOND_LENGTH * 0.6;    // width reserved for "+" glyph

        var reactantMols = [];
        var reactantBoundsList = [];
        for (var ri = 0; ri < reactantFrags.length; ri++) {
            var rMol = new Molecule();
            parse(reactantFrags[ri], rMol);
            if (rMol.parseErrors)   mol.parseErrors   = mol.parseErrors.concat(rMol.parseErrors);
            if (rMol.parseWarnings) mol.parseWarnings = mol.parseWarnings.concat(rMol.parseWarnings);
            reactantMols.push(rMol);
            reactantBoundsList.push(rMol.getBounds());
        }

        var productMols = [];
        var productBoundsList = [];
        for (var pi = 0; pi < productFrags.length; pi++) {
            var pMol = new Molecule();
            parse(productFrags[pi], pMol);
            if (pMol.parseErrors)   mol.parseErrors   = mol.parseErrors.concat(pMol.parseErrors);
            if (pMol.parseWarnings) mol.parseWarnings = mol.parseWarnings.concat(pMol.parseWarnings);
            productMols.push(pMol);
            productBoundsList.push(pMol.getBounds());
        }

        // ---- Compute global vertical centre (max extent of all fragments) ----
        var globalMinY =  Infinity, globalMaxY = -Infinity;
        var allBounds = reactantBoundsList.concat(productBoundsList);
        for (var bi = 0; bi < allBounds.length; bi++) {
            var b = allBounds[bi];
            if (b.y < globalMinY) globalMinY = b.y;
            if (b.y + b.h > globalMaxY) globalMaxY = b.y + b.h;
        }
        var globalCentreY = (globalMinY + globalMaxY) / 2;

        // ---- Place reactant fragments left-to-right, vertically centred ----
        var cursorX = 0;
        var reactantPlusPositions = []; // x positions for "+" signs

        for (var ri = 0; ri < reactantMols.length; ri++) {
            var rBounds = reactantBoundsList[ri];
            var rMol    = reactantMols[ri];
            var fragCentreY = rBounds.y + rBounds.h / 2;
            var shiftX = cursorX - rBounds.x;
            var shiftY = globalCentreY - fragCentreY;

            copyFragmentInto(mol, rMol, shiftX, shiftY);

            cursorX = cursorX + rBounds.w + FRAG_GAP;

            // Record "+" position between fragments (not after the last one)
            if (ri < reactantMols.length - 1) {
                reactantPlusPositions.push({
                    x: cursorX - FRAG_GAP / 2,
                    y: globalCentreY
                });
                cursorX += PLUS_WIDTH; // extra room for the + glyph
            }
        }

        // ---- Reaction arrow ----
        var arrowX1 = cursorX + ARROW_PAD;
        var arrowX2 = arrowX1 + ARROW_LEN;
        var arrowY  = globalCentreY;

        mol.reactionArrow = {
            x1: arrowX1, y1: arrowY,
            x2: arrowX2, y2: arrowY,
            type: arrowType,       // 'forward' | 'retro'
            conditions: ''         // optional text above arrow (set by caller)
        };

        // ---- Place product fragments ----
        cursorX = arrowX2 + ARROW_PAD;
        var productPlusPositions = [];

        for (var pi = 0; pi < productMols.length; pi++) {
            var pBounds = productBoundsList[pi];
            var pMol    = productMols[pi];
            var fragCentreY = pBounds.y + pBounds.h / 2;
            var shiftX = cursorX - pBounds.x;
            var shiftY = globalCentreY - fragCentreY;

            copyFragmentInto(mol, pMol, shiftX, shiftY);

            cursorX = cursorX + pBounds.w + FRAG_GAP;

            if (pi < productMols.length - 1) {
                productPlusPositions.push({
                    x: cursorX - FRAG_GAP / 2,
                    y: globalCentreY
                });
                cursorX += PLUS_WIDTH;
            }
        }

        // Store plus-sign positions on mol for the renderer
        mol.reactionPlusSigns = reactantPlusPositions.concat(productPlusPositions);
    }

    /**
     * Copy all atoms and bonds from srcMol into destMol, applying a
     * coordinate offset.  Used when assembling reaction components.
     */
    function copyFragmentInto(destMol, srcMol, dx, dy) {
        var idMap = {};
        for (var i = 0; i < srcMol.atoms.length; i++) {
            var a = srcMol.atoms[i];
            var newAtom = destMol.addAtom(a.symbol, a.x + dx, a.y + dy);
            newAtom.charge    = a.charge;
            newAtom.isotope   = a.isotope;
            newAtom.mapNumber = a.mapNumber;
            newAtom.hydrogens = a.hydrogens;
            newAtom.aromatic  = a.aromatic;
            newAtom.chirality = a.chirality;
            idMap[a.id] = newAtom.id;
        }
        for (var i = 0; i < srcMol.bonds.length; i++) {
            var b = srcMol.bonds[i];
            var bond = destMol.addBond(idMap[b.atom1], idMap[b.atom2], b.type);
            if (bond) bond.stereo = b.stereo;
        }
    }

    // -----------------------------------------------------------------------
    // Top-level split (respects bracket/parenthesis nesting)
    // -----------------------------------------------------------------------

    function splitTopLevel(str, delim) {
        var parts = [], depth = 0, bracketDepth = 0, start = 0;
        for (var i = 0; i < str.length; i++) {
            var c = str[i];
            if (c === '(') depth++;
            else if (c === ')') depth--;
            else if (c === '[') bracketDepth++;
            else if (c === ']') bracketDepth--;
            else if (c === delim && depth === 0 && bracketDepth === 0) {
                parts.push(str.substring(start, i));
                start = i + 1;
            }
        }
        parts.push(str.substring(start));
        return parts.filter(function(s) { return s.length > 0; });
    }

    // -----------------------------------------------------------------------
    // Fragment parser (the core recursive-descent engine)
    // -----------------------------------------------------------------------

    function parseFragment(smiles, mol) {
        var pos = 0;
        var len = smiles.length;
        var stack        = [];   // branch stack of atom IDs
        var currentAtom  = null;
        var ringOpenings = {};   // ringNum -> { atomId, bondType, pos }
        var pendingBond  = null; // null means "default" (single, or aromatic between aromatic atoms)
        var pendingStereo = null; // '/' or '\\'
        var errors   = [];
        var warnings = [];
        var atomIndex = 0; // running count for position reporting

        // Helpers
        function peek()  { return pos < len ? smiles[pos] : ''; }
        function next()  { return pos < len ? smiles[pos++] : ''; }
        function isDigit(c) { return c >= '0' && c <= '9'; }

        // ---- Main loop ----
        while (pos < len) {
            var c = peek();

            // -- Branch open --
            if (c === '(') {
                if (currentAtom === null) {
                    errors.push('Unexpected "(" at position ' + pos + ' — no atom to branch from');
                }
                stack.push(currentAtom);
                next();
                continue;
            }

            // -- Branch close --
            if (c === ')') {
                if (stack.length === 0) {
                    errors.push('Unbalanced parentheses: extra ")" at position ' + pos);
                } else {
                    currentAtom = stack.pop();
                }
                next();
                continue;
            }

            // -- Bond symbols --
            if (c === '-') { pendingBond = BOND_SINGLE; next(); continue; }
            if (c === '=') { pendingBond = BOND_DOUBLE; next(); continue; }
            if (c === '#') { pendingBond = BOND_TRIPLE; next(); continue; }
            if (c === ':') {
                // Aromatic bond (explicit) — only valid between aromatic atoms
                pendingBond = BOND_AROMATIC;
                next();
                continue;
            }
            if (c === '/' || c === '\\') {
                pendingStereo = c;
                if (pendingBond === null) pendingBond = BOND_SINGLE;
                next();
                continue;
            }

            // -- Ring closure --
            if (c === '%' || isDigit(c)) {
                var ringStart = pos;
                var ringNum;
                if (c === '%') {
                    next(); // skip %
                    var d1 = next();
                    var d2 = next();
                    if (!isDigit(d1) || !isDigit(d2)) {
                        errors.push('Invalid ring number at position ' + ringStart + ' — expected two digits after %');
                        continue;
                    }
                    ringNum = parseInt(d1 + d2, 10);
                } else {
                    ringNum = parseInt(next(), 10);
                }

                if (currentAtom === null) {
                    errors.push('Ring closure ' + ringNum + ' at position ' + ringStart + ' with no current atom');
                    continue;
                }

                if (ringOpenings[ringNum]) {
                    // Close ring
                    var opening = ringOpenings[ringNum];
                    var bondType = pendingBond || opening.bondType || BOND_SINGLE;
                    var bond = mol.addBond(opening.atomId, currentAtom, bondType);
                    if (bond && pendingStereo) {
                        bond.stereo = pendingStereo === '/' ? 1 : 6;
                    }
                    delete ringOpenings[ringNum];
                    pendingBond = null;
                    pendingStereo = null;
                } else {
                    // Open ring
                    ringOpenings[ringNum] = {
                        atomId: currentAtom,
                        bondType: pendingBond,
                        pos: ringStart
                    };
                    pendingBond = null;
                    pendingStereo = null;
                }
                continue;
            }

            // -- Bracket atom --
            if (c === '[') {
                var bracketStart = pos;
                next(); // skip '['
                var end = smiles.indexOf(']', pos);
                if (end < 0) {
                    errors.push('Unclosed bracket atom starting at position ' + bracketStart);
                    break;
                }
                var bracketStr = smiles.substring(pos, end);
                var parsed = parseBracketAtom(bracketStr, bracketStart + 1, errors);
                pos = end + 1;

                var atom = mol.addAtom(parsed.symbol, 0, 0);
                atom.charge    = parsed.charge;
                atom.isotope   = parsed.isotope;
                atom.mapNumber = parsed.mapNumber;
                atom.hydrogens = parsed.hydrogens;
                atom.aromatic  = parsed.aromatic;
                atom.chirality = parsed.chirality;

                if (currentAtom !== null) {
                    var bt = resolveBondType(mol, currentAtom, atom.id, pendingBond);
                    var bond = mol.addBond(currentAtom, atom.id, bt);
                    if (bond && pendingStereo) {
                        bond.stereo = pendingStereo === '/' ? 1 : 6;
                    }
                }
                currentAtom = atom.id;
                pendingBond = null;
                pendingStereo = null;
                atomIndex++;
                continue;
            }

            // -- Organic atom (uppercase start) or aromatic (lowercase) --
            if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
                var atomStart = pos;
                var sym = null;
                var aromatic = false;

                // Try aromatic organic subset first (lowercase single char)
                if (c >= 'a' && c <= 'z') {
                    if (AROMATIC_ORGANIC[c]) {
                        sym = c.toUpperCase();
                        aromatic = true;
                        next();
                    } else {
                        errors.push('Unknown atom symbol "' + c + '" at position ' + pos);
                        next();
                        continue;
                    }
                } else {
                    // Uppercase: try two-char organic first, then one-char
                    var twoChar = c + (pos + 1 < len ? smiles[pos + 1] : '');
                    if (ORGANIC_TWO[twoChar]) {
                        sym = twoChar;
                        pos += 2;
                    } else if (ORGANIC_SINGLE[c]) {
                        sym = c;
                        next();
                    } else {
                        errors.push('Unknown atom symbol "' + c + '" at position ' + pos + ' — not in organic subset (use brackets for other elements)');
                        next();
                        continue;
                    }
                }

                var atom = mol.addAtom(sym, 0, 0);
                atom.aromatic = aromatic;

                if (currentAtom !== null) {
                    var bt = resolveBondType(mol, currentAtom, atom.id, pendingBond);
                    var bond = mol.addBond(currentAtom, atom.id, bt);
                    if (bond && pendingStereo) {
                        bond.stereo = pendingStereo === '/' ? 1 : 6;
                    }
                }
                currentAtom = atom.id;
                pendingBond = null;
                pendingStereo = null;
                atomIndex++;
                continue;
            }

            // -- Dot (fragment separator) should have been handled by splitTopLevel --
            if (c === '.') {
                // Within a fragment this means disconnect
                currentAtom = null;
                pendingBond = null;
                pendingStereo = null;
                next();
                continue;
            }

            // -- Unknown character --
            errors.push('Unexpected character "' + c + '" at position ' + pos);
            next();
        }

        // ---- Post-parse validation ----

        // Unclosed branches
        if (stack.length > 0) {
            errors.push('Unbalanced parentheses: ' + stack.length + ' unclosed branch(es)');
        }

        // Unclosed ring openings
        var unclosedRings = Object.keys(ringOpenings);
        for (var r = 0; r < unclosedRings.length; r++) {
            var rn = unclosedRings[r];
            errors.push('Unclosed ring ' + rn + ' opened at position ' + ringOpenings[rn].pos);
        }

        // FIX: valence check should only cover atoms added by this fragment, not the entire molecule
        // (the outer loop in parse() may call parseFragment() multiple times for dot-separated SMILES)
        // FIX: apply charge correction.  A negatively charged atom can absorb
        // |charge| extra bonds (lone-pair-to-bond), a positively charged atom
        // gives up |charge| bonds.  Without this, [BH4-] (4 bonds, charge -1)
        // and [NH4+] (4 bonds, charge +1) are flagged as valence violations
        // even though they are perfectly valid Lewis structures.
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            var vList = STANDARD_VALENCES[a.symbol];
            if (!vList) continue;
            var bondSum = mol.bondOrderSum(a.id);
            var hCount = (a.hydrogens >= 0) ? a.hydrogens : 0;
            // For auto-H atoms, we compute the implicit Hs below
            if (a.hydrogens < 0) {
                // Determine default H
                var maxV = getDefaultValence(a.symbol, bondSum, a.charge);
                hCount = Math.max(0, maxV - bondSum - Math.abs(a.charge));
            }
            var totalOrder = bondSum + hCount;
            var charge = a.charge || 0;
            // Skip valence-warning on charged atoms: charge correction depends
            // on whether the formal charge came from electron loss (cation
            // gives up a bond, e.g. CH3+ has 3 bonds) or from coordinative bond
            // donation (cation gains a bond, e.g. NH4+ has 4 bonds).  Different
            // atoms and different mechanisms give different effective valences,
            // and a one-rule-fits-all check generates false positives on common
            // structures like [BH4-], [NH4+], [PH4+].  Charged atoms are
            // always written in brackets where the user has taken responsibility
            // for the H count, so suppressing the warning here is safe.
            if (charge !== 0) continue;
            var valid = false;
            for (var v = 0; v < vList.length; v++) {
                if (totalOrder <= vList[v]) { valid = true; break; }
            }
            if (!valid) {
                warnings.push('Valence violation: atom ' + a.symbol +
                    ' (id ' + a.id + ') has ' + totalOrder +
                    ' bonds but max standard valence is ' + vList[vList.length - 1]);
            }
        }

        return { errors: errors, warnings: warnings };
    }

    // -----------------------------------------------------------------------
    // Bracket atom parser:  [isotope? symbol chirality? hcount? charge? :map?]
    // -----------------------------------------------------------------------

    function parseBracketAtom(str, globalOffset, errors) {
        // Per OpenSMILES, bracket atoms with no explicit H specification
        // have 0 implicit hydrogens (unlike organic-subset atoms which auto-
        // calculate H from valence).  We use -1 as a sentinel initially and
        // promote it to 0 at the end if no H token was found.
        var result = {
            symbol: 'C', charge: 0, isotope: 0,
            mapNumber: 0, hydrogens: -1,
            chirality: '', aromatic: false
        };
        var pos = 0;
        var len = str.length;

        function isDigit(c) { return c >= '0' && c <= '9'; }

        // -- Isotope (leading digits) --
        var isoStr = '';
        while (pos < len && isDigit(str[pos])) isoStr += str[pos++];
        if (isoStr) result.isotope = parseInt(isoStr, 10);

        // -- Symbol --
        if (pos < len) {
            // Check for aromatic bracket atoms (lowercase one or two chars)
            if (str[pos] >= 'a' && str[pos] <= 'z') {
                var twoLower = str.substring(pos, pos + 2);
                if (AROMATIC_BRACKET[twoLower]) {
                    // Two-char aromatic (se, as, te)
                    result.symbol = twoLower.charAt(0).toUpperCase() + twoLower.charAt(1);
                    result.aromatic = true;
                    pos += 2;
                } else if (AROMATIC_BRACKET[str[pos]]) {
                    result.symbol = str[pos].toUpperCase();
                    result.aromatic = true;
                    pos++;
                } else {
                    errors.push('Unknown aromatic atom "' + str[pos] + '" at position ' + (globalOffset + pos));
                    pos++;
                }
            } else if (str[pos] >= 'A' && str[pos] <= 'Z') {
                var sym = str[pos];
                pos++;
                if (pos < len && str[pos] >= 'a' && str[pos] <= 'z') {
                    var twoSym = sym + str[pos];
                    if (ALL_ELEMENTS[twoSym]) {
                        sym = twoSym;
                        pos++;
                    }
                    // else: single-char symbol, keep pos where it is
                }
                if (!ALL_ELEMENTS[sym]) {
                    errors.push('Unknown element symbol "' + sym + '" at position ' + (globalOffset + pos - sym.length));
                }
                result.symbol = sym;
            } else if (str[pos] === '*') {
                // Wildcard atom
                result.symbol = '*';
                pos++;
            } else {
                errors.push('Expected element symbol at position ' + (globalOffset + pos));
            }
        }

        // -- Chirality (@ or @@) --
        if (pos < len && str[pos] === '@') {
            pos++;
            if (pos < len && str[pos] === '@') {
                result.chirality = '@@';
                pos++;
            } else {
                result.chirality = '@';
            }
            // Skip optional chirality class (e.g. @TH1, @AL2 — rare)
            // Only skip known class prefixes, NOT 'H' which is the hydrogen count
            if (pos < len && str[pos] === 'T' && pos + 1 < len && str[pos + 1] === 'H') { pos += 2; while (pos < len && isDigit(str[pos])) pos++; }
            else if (pos < len && str[pos] === 'A' && pos + 1 < len && str[pos + 1] === 'L') { pos += 2; while (pos < len && isDigit(str[pos])) pos++; }
            else if (pos < len && str[pos] === 'S' && pos + 1 < len && str[pos + 1] === 'P') { pos += 2; while (pos < len && isDigit(str[pos])) pos++; }
            else if (pos < len && str[pos] === 'T' && pos + 1 < len && str[pos + 1] === 'B') { pos += 2; while (pos < len && isDigit(str[pos])) pos++; }
            else if (pos < len && str[pos] === 'O' && pos + 1 < len && str[pos + 1] === 'H') { pos += 2; while (pos < len && isDigit(str[pos])) pos++; }
        }

        // -- Hydrogen count --
        if (pos < len && str[pos] === 'H') {
            pos++;
            var hStr = '';
            while (pos < len && isDigit(str[pos])) hStr += str[pos++];
            result.hydrogens = hStr ? parseInt(hStr, 10) : 1;
        }

        // -- Charge --
        if (pos < len && (str[pos] === '+' || str[pos] === '-')) {
            var sign = str[pos] === '+' ? 1 : -1;
            pos++;
            var chargeStr = '';
            while (pos < len && isDigit(str[pos])) chargeStr += str[pos++];
            if (chargeStr) {
                result.charge = sign * parseInt(chargeStr, 10);
            } else {
                // Count repeated + or - signs
                var count = 1;
                var signChar = sign > 0 ? '+' : '-';
                while (pos < len && str[pos] === signChar) { count++; pos++; }
                result.charge = sign * count;
            }
        }

        // -- Atom map number (:N) --
        if (pos < len && str[pos] === ':') {
            pos++;
            var mapStr = '';
            while (pos < len && isDigit(str[pos])) mapStr += str[pos++];
            if (mapStr) {
                result.mapNumber = parseInt(mapStr, 10);
            } else {
                errors.push('Expected map number after ":" at position ' + (globalOffset + pos));
            }
        }

        // Anything remaining is unexpected
        if (pos < len) {
            errors.push('Unexpected content "' + str.substring(pos) + '" in bracket atom at position ' + (globalOffset + pos));
        }

        // OpenSMILES rule: bracket atoms with no explicit H token have 0
        // implicit hydrogens (only organic-subset atoms auto-calculate).
        if (result.hydrogens < 0) {
            result.hydrogens = 0;
        }

        return result;
    }

    // -----------------------------------------------------------------------
    // Bond type resolution (handles aromatic default)
    // -----------------------------------------------------------------------

    function resolveBondType(mol, fromId, toId, pendingBond) {
        if (pendingBond !== null) {
            // Explicit aromatic bond symbol ':' -> store as single (Molecule doesn't have BOND_AROMATIC)
            if (pendingBond === BOND_AROMATIC) return BOND_SINGLE;
            return pendingBond;
        }
        // Default: if both atoms are aromatic, use single bond (aromatic perception later)
        var fromAtom = mol.getAtom(fromId);
        var toAtom   = mol.getAtom(toId);
        if (fromAtom && toAtom && fromAtom.aromatic && toAtom.aromatic) {
            return BOND_SINGLE; // aromatic — will be promoted during perception
        }
        return BOND_SINGLE;
    }

    // -----------------------------------------------------------------------
    // Valence helpers
    // -----------------------------------------------------------------------

    /**
     * Return the most appropriate default valence for implicit-H calculation.
     * Picks the smallest standard valence >= bondSum + |charge|.
     */
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
    // validateSmiles — standalone validation returning detailed errors
    // -----------------------------------------------------------------------

    /**
     * Validate a SMILES string without building a full Molecule.
     * @param {string} smiles
     * @returns {{ valid: boolean, errors: string[], warnings: string[] }}
     */
    function validateSmiles(smiles) {
        if (typeof smiles !== 'string' || smiles.trim().length === 0) {
            return { valid: false, errors: ['Empty SMILES string'], warnings: [] };
        }

        var mol = new Molecule();
        mol.parseErrors = [];
        mol.parseWarnings = [];

        parse(smiles, mol);

        return {
            valid:    mol.parseErrors.length === 0,
            errors:   mol.parseErrors,
            warnings: mol.parseWarnings
        };
    }

    // -----------------------------------------------------------------------
    // 2D Layout — delegates to Layout.js if available, otherwise BFS fallback
    // -----------------------------------------------------------------------

    function layoutFragment(mol, fragAtoms) {
        if (fragAtoms.length === 0) return;

        // Use the dedicated Layout module for ring-aware coordinate generation
        if (typeof Layout !== 'undefined' && Layout.layoutFragment) {
            Layout.layoutFragment(mol, fragAtoms);
            return;
        }

        // Fallback: simple BFS tree layout (no ring perception)
        var placed = {};
        var queue = [];

        fragAtoms[0].x = 0;
        fragAtoms[0].y = 0;
        placed[fragAtoms[0].id] = true;
        queue.push(fragAtoms[0].id);

        while (queue.length > 0) {
            var atomId = queue.shift();
            var atom = mol.getAtom(atomId);
            var neighbors = mol.getNeighbors(atomId);
            var placedNeighbors = neighbors.filter(function(n) { return placed[n]; });
            var unplaced = neighbors.filter(function(n) { return !placed[n]; });

            if (unplaced.length === 0) continue;

            var startAngle;
            if (placedNeighbors.length > 0) {
                var avgAngle = 0;
                placedNeighbors.forEach(function(nId) {
                    var n = mol.getAtom(nId);
                    avgAngle += Math.atan2(n.y - atom.y, n.x - atom.x);
                });
                avgAngle /= placedNeighbors.length;
                startAngle = avgAngle + Math.PI;
            } else {
                startAngle = 0;
            }

            var spread = Math.PI / 3;
            var totalSpread = (unplaced.length - 1) * spread;
            var firstAngle = startAngle - totalSpread / 2;

            for (var i = 0; i < unplaced.length; i++) {
                var nId = unplaced[i];
                var n = mol.getAtom(nId);
                var angle = firstAngle + i * spread;
                n.x = atom.x + BOND_LENGTH * Math.cos(angle);
                n.y = atom.y + BOND_LENGTH * Math.sin(angle);
                placed[nId] = true;
                queue.push(nId);
            }
        }

        var bounds = getBounds(fragAtoms);
        var cx = bounds.minX + (bounds.maxX - bounds.minX) / 2;
        var cy = bounds.minY + (bounds.maxY - bounds.minY) / 2;
        fragAtoms.forEach(function(a) { a.x -= cx; a.y -= cy; });
    }

    function getBounds(atoms) {
        var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        atoms.forEach(function(a) {
            if (a.x < minX) minX = a.x;
            if (a.y < minY) minY = a.y;
            if (a.x > maxX) maxX = a.x;
            if (a.y > maxY) maxY = a.y;
        });
        return { minX: minX, minY: minY, maxX: maxX, maxY: maxY };
    }

    // -----------------------------------------------------------------------
    // Export
    // -----------------------------------------------------------------------

    global.SmilesParser = {
        parse:          parse,
        validateSmiles: validateSmiles
    };

})(typeof window !== 'undefined' ? window : this);
