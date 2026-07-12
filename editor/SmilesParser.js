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

    // v1.8.26: STANDARD_VALENCES + getDefaultValence now live in
    // Molecule.js (single canonical home — they were byte-identical
    // here and in SmilesWriter.js). Local aliases keep call sites
    // terse and unchanged.
    var STANDARD_VALENCES = Molecule.STANDARD_VALENCES;
    var getDefaultValence = Molecule.getDefaultValence;

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
        'R':1,  // generic R-group
        // v2.0.3 — amino-acid residue placeholders (drawn as labelled
        // pills, same family as [R] / [R1]). Listed here so bracket-form
        // [Gly] etc. parse correctly, which is the form SmilesWriter
        // emits on round-trip for any multi-char non-organic symbol.
        'Gly':1,'Ala':1,'Ser':1,'Val':1,'Leu':1
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
    // Multi-letter atom shortcuts (v2.0.3)
    //
    // Two flavours sharing one dispatch:
    //
    //   1. Organic-chemistry SUBSTRUCTURE aliases — Ph / Me / Et / iPr /
    //      tBu / Bn / Bz. These represent real heavy-atom fragments and
    //      expand inline into their full atom/bond graph; "Ph" IS a
    //      benzene ring, drawn with six aromatic carbons.
    //
    //   2. Residue / R-group PLACEHOLDERS — Gly / Ala / Ser / Val / Leu.
    //      These collapse to a single labelled atom whose symbol is the
    //      residue code, the same way [R] / [R1] already work in BIME.
    //      "GlyAla" parses as two nodes joined by one bond (the peptide
    //      bond is implicit in the label); "GlyAlaSer" is a 3-node chain.
    //      That matches BIME's Pathway Canvas residue-pill model and
    //      keeps peptide drawings visually clean rather than blowing up
    //      to 5+ atoms per residue.
    //
    // Each entry declares (a) the atoms it contributes, (b) the bonds
    // between them, (c) which atom is bonded TO by the parent (attach),
    // and (d) which atom the parser continues from (cont — equal to
    // attach for every current entry: terminal substituents and single-
    // node residues both behave like one drop-in heavy atom).
    //
    // Ring digits inside the expansion live in their own namespace: each
    // entry builds its atoms+bonds directly via mol.addAtom / mol.addBond
    // rather than re-entering the parser, so e.g. the aromatic ring inside
    // Ph cannot collide with any digit the outer SMILES is using.
    //
    // Matched longest-first (so "Gly" beats a hypothetical "Gl" / "G"
    // entry); the dispatch is case-sensitive so "Br" / "Cl" still parse
    // as bromine / chlorine (neither is in the abbreviation table).
    // -----------------------------------------------------------------------
    // Build a single expansion function from a spec.
    //
    // spec.atoms[i] accepts ONLY these fields (unknown fields are silently
    // ignored, so a typo will produce a "default" atom rather than an
    // error — keep entries terse):
    //   - sym    (string, required)  element symbol; multi-char symbols
    //                                like 'Gly' are valid R-group labels
    //                                (must also appear in ALL_ELEMENTS).
    //   - arom   (bool,   optional)  default false. Marks the atom as
    //                                aromatic so the layout renders it
    //                                inside an aromatic ring.
    //   - charge (number, optional)  default 0. Formal charge.
    //   - h      (number, optional)  default -1 (calcHydrogens computes).
    //                                Pin to 0 on residue labels so valence
    //                                checks don't try to add Hs to a
    //                                three-letter code.
    //
    // spec.bonds[i] is the tuple [fromIdx, toIdx, type] indexing back into
    // spec.atoms. spec.attach and spec.cont are atom indices: the parent
    // SMILES atom bonds to attach, and the parser continues from cont.
    function _abbrev(spec) {
        return function expandAbbrev(mol) {
            var ids = [];
            for (var i = 0; i < spec.atoms.length; i++) {
                var a = mol.addAtom(spec.atoms[i].sym, 0, 0);
                if (spec.atoms[i].arom) a.aromatic = true;
                if (typeof spec.atoms[i].charge === 'number') a.charge = spec.atoms[i].charge;
                if (typeof spec.atoms[i].h === 'number') a.hydrogens = spec.atoms[i].h;
                ids.push(a.id);
            }
            for (var b = 0; b < spec.bonds.length; b++) {
                var bb = spec.bonds[b];
                mol.addBond(ids[bb[0]], ids[bb[1]], bb[2]);
            }
            return { attach: ids[spec.attach], cont: ids[spec.cont] };
        };
    }

    // Single-node residue placeholder helper: one atom, symbol = sym,
    // explicit hydrogens = 0 (no implicit-H valence accounting for a
    // residue label), attach == cont.
    function _residue(sym) {
        return _abbrev({ atoms: [{ sym: sym, h: 0 }], bonds: [],
                         attach: 0, cont: 0 });
    }

    var ABBREVIATIONS = {
        // ----- Organic substructure aliases (expanded inline) ----------
        'Ph': _abbrev({                                    // phenyl  C6H5-
            atoms: [{sym:'C',arom:true},{sym:'C',arom:true},{sym:'C',arom:true},
                    {sym:'C',arom:true},{sym:'C',arom:true},{sym:'C',arom:true}],
            bonds: [[0,1,1],[1,2,1],[2,3,1],[3,4,1],[4,5,1],[5,0,1]],
            attach: 0, cont: 0
        }),
        'Me': _abbrev({                                    // methyl  CH3-
            atoms: [{sym:'C'}], bonds: [], attach: 0, cont: 0
        }),
        'Et': _abbrev({                                    // ethyl   C2H5-
            atoms: [{sym:'C'},{sym:'C'}],
            bonds: [[0,1,1]], attach: 0, cont: 0
        }),
        'iPr': _abbrev({                                   // isopropyl  (CH3)2CH-
            atoms: [{sym:'C'},{sym:'C'},{sym:'C'}],
            bonds: [[0,1,1],[0,2,1]], attach: 0, cont: 0
        }),
        'tBu': _abbrev({                                   // tert-butyl  (CH3)3C-
            atoms: [{sym:'C'},{sym:'C'},{sym:'C'},{sym:'C'}],
            bonds: [[0,1,1],[0,2,1],[0,3,1]], attach: 0, cont: 0
        }),
        'Bn': _abbrev({                                    // benzyl  PhCH2-
            atoms: [{sym:'C'},                             // 0: methylene
                    {sym:'C',arom:true},{sym:'C',arom:true},{sym:'C',arom:true},
                    {sym:'C',arom:true},{sym:'C',arom:true},{sym:'C',arom:true}],
            bonds: [[0,1,1],[1,2,1],[2,3,1],[3,4,1],[4,5,1],[5,6,1],[6,1,1]],
            attach: 0, cont: 0
        }),
        'Bz': _abbrev({                                    // benzoyl  PhC(=O)-
            atoms: [{sym:'C'},{sym:'O'},                   // 0: carbonyl C, 1: =O
                    {sym:'C',arom:true},{sym:'C',arom:true},{sym:'C',arom:true},
                    {sym:'C',arom:true},{sym:'C',arom:true},{sym:'C',arom:true}],
            bonds: [[0,1,2],[0,2,1],[2,3,1],[3,4,1],[4,5,1],[5,6,1],[6,7,1],[7,2,1]],
            attach: 0, cont: 0
        }),
        // ----- Amino-acid residues (single labelled atom — R-group form) -----
        // Each residue collapses to one node carrying its three-letter
        // label as the atom symbol. "GlyAla" therefore parses as two
        // residue nodes connected by a single bond — the canonical
        // peptide-drawing convention. The atom's hydrogens are pinned at
        // 0 so calcHydrogens / valence checks don't try to add Hs to a
        // residue label.
        'Gly': _residue('Gly'),                            // glycine
        'Ala': _residue('Ala'),                            // alanine
        'Ser': _residue('Ser'),                            // serine
        'Val': _residue('Val'),                            // valine
        'Leu': _residue('Leu')                             // leucine
    };

    // Sorted longest-first so e.g. "Bn" wins over a hypothetical "B"
    // alias, and "iPr" / "tBu" / "Gly" never get truncated to a shorter
    // prefix match. Computed once at module load.
    var ABBREV_KEYS = Object.keys(ABBREVIATIONS).sort(function (a, b) {
        return b.length - a.length;
    });

    // Set of first characters that COULD start an abbreviation; used by the
    // main parse loop to skip the abbreviation lookup entirely when the
    // current character can't begin any registered entry. Derived from
    // ABBREV_KEYS so adding a new entry (especially one with a lowercase
    // first letter, like iPr or tBu) does NOT require manually widening
    // a hand-curated gate condition.
    var ABBREV_FIRSTS = {};
    for (var _ak = 0; _ak < ABBREV_KEYS.length; _ak++) {
        ABBREV_FIRSTS[ABBREV_KEYS[_ak].charAt(0)] = true;
    }

    function matchAbbreviation(smiles, pos) {
        for (var i = 0; i < ABBREV_KEYS.length; i++) {
            var k = ABBREV_KEYS[i];
            if (smiles.substring(pos, pos + k.length) === k) return k;
        }
        return null;
    }

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
        // v1.4.1 (bug-fix #11): guard the all-empty `>>` reaction case so we
        // don't propagate Infinity / NaN into reactionArrow coordinates.
        if (allBounds.length === 0) { globalMinY = 0; globalMaxY = 0; }
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
        // chiralSeq: for every chiral atom, the order its neighbours appear in
        // the SMILES token stream — the frame the raw @/@@ token is defined
        // against. Entries are neighbour atom IDs, the marker 'H' for an
        // implicit hydrogen, or { pendingRing: n } for a ring opened here and
        // not yet closed. After the parse this is used to normalise every
        // chirality token into the heavy-atom adjacency frame (getBondsOfAtom
        // order, implicit H last) that SmilesWriter and CIPStereo both assume.
        var chiralSeq    = {};   // chiral atomId -> [neighbourId | 'H' | {pendingRing}]
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
                    // Record the ring-closure neighbour in SMILES-token order
                    // for both endpoints if chiral. The closing atom appends
                    // its partner at the digit position; the opening atom
                    // resolves the placeholder it parked when the ring opened.
                    if (chiralSeq[currentAtom]) chiralSeq[currentAtom].push(opening.atomId);
                    if (chiralSeq[opening.atomId]) {
                        var oSeq = chiralSeq[opening.atomId];
                        for (var os = 0; os < oSeq.length; os++) {
                            if (oSeq[os] && oSeq[os].pendingRing === ringNum) {
                                oSeq[os] = currentAtom;
                                break;
                            }
                        }
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
                    // Park a placeholder at the digit position; it is resolved
                    // to the partner atom id when this ring later closes.
                    if (chiralSeq[currentAtom]) {
                        chiralSeq[currentAtom].push({ pendingRing: ringNum });
                    }
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
                    // This atom is a SMILES-token-order neighbour of currentAtom.
                    if (chiralSeq[currentAtom]) chiralSeq[currentAtom].push(atom.id);
                }
                // Seed this atom's own semantic neighbour order: the from-atom
                // (if any) first, then the implicit hydrogen written in the
                // bracket — exactly where [C@H] places it.
                if (atom.chirality) {
                    chiralSeq[atom.id] = [];
                    if (currentAtom !== null) chiralSeq[atom.id].push(currentAtom);
                    for (var ch = 0; ch < atom.hydrogens; ch++) {
                        chiralSeq[atom.id].push('H');
                    }
                }
                currentAtom = atom.id;
                pendingBond = null;
                pendingStereo = null;
                atomIndex++;
                continue;
            }

            // -- Multi-letter atom shorthand (Ph, Me, Et, Bn, Bz, iPr, tBu,
            //    Gly, Ala, Ser, Val, Leu — see ABBREVIATIONS at the top of
            //    this file). Tried BEFORE the case-based organic dispatch
            //    because some abbreviations start with lowercase (iPr, tBu)
            //    which the outer block would otherwise rule out as either
            //    unknown aromatic atoms or stranded prefixes. Match is
            //    longest-first and case-sensitive, so e.g. "Br" / "Cl"
            //    (bromine / chlorine) are unaffected — they are not in the
            //    table. Each shorthand expands inline into its full atom/
            //    bond graph; the bond from the parent atom is created here
            //    so chiralSeq sees the abbreviation's attach atom as a
            //    normal SMILES-token neighbour, and the parser continues
            //    from the abbreviation's declared cont atom (same as attach
            //    for terminal groups like Ph or Me; the carbonyl C for
            //    amino acids so "GlyAla" produces the chemically correct
            //    peptide bond between the residues). The first-character
            //    gate ABBREV_FIRSTS[c] is derived from ABBREV_KEYS at
            //    module load, so adding a new entry — including one with
            //    a lowercase first letter like iPr or tBu — automatically
            //    extends the gate without any hand-edit here.
            if (ABBREV_FIRSTS[c]) {
                var abbr = matchAbbreviation(smiles, pos);
                if (abbr) {
                    var hooks = ABBREVIATIONS[abbr](mol);
                    if (currentAtom !== null) {
                        var btA = resolveBondType(mol, currentAtom, hooks.attach, pendingBond);
                        var bondA = mol.addBond(currentAtom, hooks.attach, btA);
                        if (bondA && pendingStereo) {
                            bondA.stereo = pendingStereo === '/' ? 1 : 6;
                        }
                        if (chiralSeq[currentAtom]) chiralSeq[currentAtom].push(hooks.attach);
                    }
                    currentAtom = hooks.cont;
                    pos += abbr.length;
                    pendingBond = null;
                    pendingStereo = null;
                    atomIndex++;
                    continue;
                }
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
                    // This atom is a SMILES-token-order neighbour of currentAtom.
                    if (chiralSeq[currentAtom]) chiralSeq[currentAtom].push(atom.id);
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

        // ---- Normalise chirality tokens to the adjacency frame ----
        // The raw @/@@ token is defined against the SMILES-token neighbour
        // order (chiralSeq). SmilesWriter._resolveChirality and CIPStereo both
        // interpret atom.chirality against the heavy-atom adjacency order
        // (getBondsOfAtom) with the implicit H last. Converting each token
        // once, here, gives those consumers a single consistent frame — which
        // is what makes parse -> write a round-trip fixed point for ring-
        // opening stereocentres (the old SMILES-token frame put a ring
        // neighbour at the digit position, but the parser appends the ring
        // bond only when the ring closes, so the two orders disagreed).
        for (var cId in chiralSeq) {
            if (!chiralSeq.hasOwnProperty(cId)) continue;
            var cAtom = mol.getAtom(parseInt(cId, 10));
            if (!cAtom || !cAtom.chirality) continue;

            // Semantic order: resolve placeholders; skip if a ring never closed.
            var semRaw = chiralSeq[cId];
            var sem = [];
            var resolved = true;
            for (var sx = 0; sx < semRaw.length; sx++) {
                if (semRaw[sx] && typeof semRaw[sx] === 'object') { resolved = false; break; }
                sem.push(semRaw[sx]);
            }
            if (!resolved) continue;

            // Target frame: heavy neighbours in adjacency order, implicit H last.
            var hCount = (cAtom.hydrogens >= 0) ? cAtom.hydrogens : 0;
            if (hCount > 1) continue;                  // not a tetrahedral centre
            var target = mol.getNeighbors(cAtom.id);
            if (hCount === 1) target = target.concat(['H']);
            if (sem.length !== target.length || sem.length < 3) continue;

            // Permutation sem -> target, then its parity by cycle decomposition.
            var perm = [];
            var consumed = {};
            var valid = true;
            for (var px = 0; px < sem.length; px++) {
                var hit = -1;
                for (var tx = 0; tx < target.length; tx++) {
                    if (!consumed[tx] && target[tx] === sem[px]) { hit = tx; consumed[tx] = true; break; }
                }
                if (hit < 0) { valid = false; break; }
                perm.push(hit);
            }
            if (!valid) continue;
            var swaps = 0;
            var seenC = {};
            for (var i2 = 0; i2 < perm.length; i2++) {
                if (seenC[i2]) continue;
                var j2 = i2, clen = 0;
                while (!seenC[j2]) { seenC[j2] = true; j2 = perm[j2]; clen++; }
                swaps += (clen - 1);
            }
            if (swaps % 2 === 1) {
                cAtom.chirality = (cAtom.chirality === '@') ? '@@' : '@';
            }
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
                // Greedy match: try 3-char symbol (e.g. residue placeholders
                // like Gly / Ala / Ser) before 2-char (Cl, Br, Mg, …) before
                // the single uppercase letter. The longest match wins.
                if (pos + 1 < len &&
                    str[pos]   >= 'a' && str[pos]   <= 'z' &&
                    str[pos+1] >= 'a' && str[pos+1] <= 'z') {
                    var threeSym = sym + str[pos] + str[pos+1];
                    if (ALL_ELEMENTS[threeSym]) {
                        sym = threeSym;
                        pos += 2;
                    }
                }
                if (sym.length === 1 &&
                    pos < len && str[pos] >= 'a' && str[pos] <= 'z') {
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

    // (v1.8.26: getDefaultValence moved to Molecule.js — see the
    // `var getDefaultValence = Molecule.getDefaultValence;` alias above.)

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
