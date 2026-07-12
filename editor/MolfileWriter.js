/**
 * MolfileWriter.js — Molecule graph → MOL V2000/V3000 file format
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 */
(function(global) {
    'use strict';

    function write(mol, v3000) {
        if (v3000) return writeV3000(mol);
        return writeV2000(mol);
    }

    function writeV2000(mol) {
        var lines = [];

        // FIX: V2000 format supports max 999 atoms/bonds; fall back to V3000 for larger molecules
        if (mol.atoms.length > 999 || mol.bonds.length > 999 || hasLargeAtomMap(mol)) {
            return writeV3000(mol);
        }

        // Header.  FIX: round-trip mol.program (line 2) and mol.comment (line 3)
        //          when set.  Otherwise emit a BIME stamp and a blank line.
        //          Spec caps each header line at 80 chars.
        lines.push(truncate(mol.name || '', 80));
        lines.push(truncate(mol.program || ('  BIME    ' + timestamp() + '2D'), 80));
        lines.push(truncate(mol.comment || '', 80));

        // Counts line
        var nAtoms = mol.atoms.length;
        var nBonds = mol.bonds.length;
        lines.push(pad(nAtoms, 3) + pad(nBonds, 3) + '  0  0  0  0  0  0  0  0999 V2000');

        // Atom block
        var atomIdxMap = {};
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            atomIdxMap[a.id] = i + 1;
            var aam = a.mapNumber || 0;
            lines.push(
                padF(a.x / 30, 10, 4) + padF(-a.y / 30, 10, 4) + padF(0, 10, 4) + ' ' +
                padR(a.symbol, 3) + ' 0' +
                pad(chargeToMol(a.charge), 3) + '  0  0  0  0  0  0  0' + pad(aam, 3) + '  0'
            );
        }

        // Bond block
        for (var i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            // FIX: MOL bond-stereo field 4 has different semantics for single vs
            // double bonds.  Single bonds: 1=wedge-up, 6=wedge-down (3D config).
            // Double bonds: 0=none, 3=cis-or-trans-either (E/Z encoded by 2D
            // geometry, not by this field).  BIME stores SMILES `/` and `\` as
            // bond.stereo=1/6 on directional single bonds, which are correctly
            // emitted as wedges.  But for double bonds we must NOT emit 1/6
            // (those would be read as ill-formed in any V2000 parser); emit 0.
            var stereo;
            if (b.type === Molecule.BOND_DOUBLE) {
                stereo = 0;
            } else {
                stereo = b.stereo === Molecule.STEREO_WEDGE ? 1 : (b.stereo === Molecule.STEREO_DASH ? 6 : 0);
            }
            lines.push(
                pad(atomIdxMap[b.atom1], 3) + pad(atomIdxMap[b.atom2], 3) +
                pad(b.type, 3) + pad(stereo, 3) + '  0  0  0'
            );
        }

        // Properties
        // FIX: V2000 M CHG/ISO lines allow max 8 entries per line; split into chunks
        var charged = mol.atoms.filter(function(a) { return a.charge !== 0; });
        for (var ci = 0; ci < charged.length; ci += 8) {
            var chunk = charged.slice(ci, ci + 8);
            var line = 'M  CHG' + pad(chunk.length, 3);
            chunk.forEach(function(a) {
                line += pad(atomIdxMap[a.id], 4) + pad(a.charge, 4);
            });
            lines.push(line);
        }

        // Isotopes
        var isotoped = mol.atoms.filter(function(a) { return a.isotope > 0; });
        for (var ii = 0; ii < isotoped.length; ii += 8) {
            var chunk = isotoped.slice(ii, ii + 8);
            var line = 'M  ISO' + pad(chunk.length, 3);
            chunk.forEach(function(a) {
                line += pad(atomIdxMap[a.id], 4) + pad(a.isotope, 4);
            });
            lines.push(line);
        }

        // FIX: write radical multiplicity (M  RAD lines).  Without this, an
        // atom with radical=1/2/3 round-trips through MOL silently as a
        // closed-shell species, breaking radical fidelity.
        var radicaled = mol.atoms.filter(function(a) { return a.radical > 0; });
        for (var ri = 0; ri < radicaled.length; ri += 8) {
            var chunk = radicaled.slice(ri, ri + 8);
            var line = 'M  RAD' + pad(chunk.length, 3);
            chunk.forEach(function(a) {
                line += pad(atomIdxMap[a.id], 4) + pad(a.radical, 4);
            });
            lines.push(line);
        }

        lines.push('M  END');
        return lines.join('\n');
    }

    function writeV3000(mol) {
        var lines = [];
        // FIX: round-trip mol.program/comment header lines like V2000 does.
        lines.push(truncate(mol.name || '', 80));
        lines.push(truncate(mol.program || ('  BIME    ' + timestamp() + '2D'), 80));
        lines.push(truncate(mol.comment || '', 80));
        lines.push('  0  0  0  0  0  0  0  0  0  0999 V3000');
        lines.push('M  V30 BEGIN CTAB');
        lines.push('M  V30 COUNTS ' + mol.atoms.length + ' ' + mol.bonds.length + ' 0 0 0');

        lines.push('M  V30 BEGIN ATOM');
        var atomIdxMap = {};
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            atomIdxMap[a.id] = i + 1;
            var line = 'M  V30 ' + (i + 1) + ' ' + a.symbol + ' ' +
                (a.x / 30).toFixed(4) + ' ' + (-a.y / 30).toFixed(4) + ' 0';
            if (a.charge !== 0) line += ' CHG=' + a.charge;
            if (a.isotope > 0) line += ' MASS=' + a.isotope;
            if (a.mapNumber > 0) line += ' AAM=' + a.mapNumber;
            // FIX: emit V3000 RAD= per-atom radical multiplicity so radicals
            // round-trip through V3000.  Spec: 1=singlet, 2=doublet, 3=triplet.
            if (a.radical > 0) line += ' RAD=' + a.radical;
            lines.push(line);
        }
        lines.push('M  V30 END ATOM');

        lines.push('M  V30 BEGIN BOND');
        for (var i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            // FIX: V3000 CFG= is only meaningful on single bonds (1=wedge-up,
            // 3=wedge-down).  Suppress it on double/triple bonds where the
            // wedge-up/down semantic does not apply (E/Z lives in coordinates).
            var cfg = '';
            if (b.type === Molecule.BOND_SINGLE) {
                cfg = b.stereo === Molecule.STEREO_WEDGE ? ' CFG=1'
                    : (b.stereo === Molecule.STEREO_DASH ? ' CFG=3' : '');
            }
            lines.push('M  V30 ' + (i + 1) + ' ' + b.type + ' ' + atomIdxMap[b.atom1] + ' ' + atomIdxMap[b.atom2] + cfg);
        }
        lines.push('M  V30 END BOND');
        lines.push('M  V30 END CTAB');
        lines.push('M  END');
        return lines.join('\n');
    }

    function writeSDF(mol) {
        var molblock = writeV2000(mol);
        var sdf = molblock;
        // FIX: emit > <NAME> and > <COMMENT> data fields when set, each
        //      followed by the spec-required blank line.  Multiple fields
        //      can stack before the $$$$ terminator.
        if (mol.name) {
            sdf += '\n> <NAME>\n' + sanitizeDataFieldValue(mol.name, 1000) + '\n';
        }
        if (mol.comment) {
            sdf += '\n> <COMMENT>\n' + sanitizeDataFieldValue(mol.comment, 1000) + '\n';
        }
        sdf += '\n$$$$\n';
        return sdf;
    }

    function writeRXN(mol, v3000) {
        if (!mol || !mol.reactionArrow) return '';
        var sides = splitReactionComponents(mol);
        var useV3000 = !!v3000 || hasLargeAtomMap(mol);
        var lines = [
            useV3000 ? '$RXN V3000' : '$RXN',
            truncate(mol.name || '', 80),
            truncate(mol.program || ('  BIME    ' + timestamp() + '2D'), 80),
            truncate(mol.comment || '', 80),
            useV3000
                ? 'M  V30 COUNTS ' + sides.reactants.length + ' ' + sides.products.length + ' 0'
                : pad(sides.reactants.length, 3) + pad(sides.products.length, 3)
        ];
        sides.reactants.concat(sides.products).forEach(function(component) {
            lines.push('$MOL');
            lines.push(write(component.mol, useV3000));
        });
        return lines.join('\n');
    }

    function splitReactionComponents(mol) {
        var arrow = mol.reactionArrow;
        var midX = ((arrow.x1 || 0) + (arrow.x2 || 0)) / 2;
        var midY = ((arrow.y1 || 0) + (arrow.y2 || 0)) / 2;
        var dx = (arrow.x2 || 0) - (arrow.x1 || 0);
        var dy = (arrow.y2 || 0) - (arrow.y1 || 0);
        var len = Math.sqrt(dx * dx + dy * dy) || 1;
        dx /= len;
        dy /= len;
        var components = mol.getComponents ? mol.getComponents() : [];
        var reactants = [];
        var products = [];
        components.forEach(function(ids) {
            var c = componentToMol(mol, ids);
            if (!c || c.ids.length === 0) return;
            var side = ((c.avgX - midX) * dx) + ((c.avgY - midY) * dy);
            if (side < 0) reactants.push(c);
            else products.push(c);
        });
        reactants.sort(function(a, b) { return a.avgX - b.avgX; });
        products.sort(function(a, b) { return a.avgX - b.avgX; });
        return { reactants: reactants, products: products };
    }

    function componentToMol(src, ids) {
        var idSet = {};
        var out = new global.Molecule();
        out.name = src.name || '';
        out.program = src.program || '';
        out.comment = src.comment || '';
        var idMap = {};
        var avgX = 0;
        var avgY = 0;
        for (var i = 0; i < ids.length; i++) idSet[ids[i]] = true;
        ids.forEach(function(id) {
            var a = src.getAtom(id);
            if (!a) return;
            var copy = out.addAtom(a.symbol, a.x, a.y);
            copy.charge = a.charge || 0;
            copy.isotope = a.isotope || 0;
            copy.mapNumber = a.mapNumber || 0;
            copy.hydrogens = a.hydrogens;
            copy.aromatic = !!a.aromatic;
            copy.chirality = a.chirality || '';
            copy.radical = a.radical || 0;
            idMap[a.id] = copy.id;
            avgX += a.x;
            avgY += a.y;
        });
        src.bonds.forEach(function(b) {
            if (!idSet[b.atom1] || !idSet[b.atom2]) return;
            var bond = out.addBond(idMap[b.atom1], idMap[b.atom2], b.type);
            if (bond) {
                bond.stereo = b.stereo || 0;
                bond.cipLabel = b.cipLabel || '';
            }
        });
        return {
            mol: out,
            ids: ids,
            avgX: ids.length ? avgX / ids.length : 0,
            avgY: ids.length ? avgY / ids.length : 0
        };
    }

    function hasLargeAtomMap(mol) {
        if (!mol || !mol.atoms) return false;
        for (var i = 0; i < mol.atoms.length; i++) {
            if ((mol.atoms[i].mapNumber || 0) > 999) return true;
        }
        return false;
    }

    // Helpers
    function pad(num, width) {
        var s = '' + num;
        while (s.length < width) s = ' ' + s;
        return s;
    }

    // FIX: clip header strings to spec maximum (80 chars) and strip any
    //      embedded newlines / SDF delimiters that would otherwise corrupt
    //      the line layout when the MOL block is embedded in an SDF.
    function truncate(str, maxLen) {
        if (str === null || str === undefined) return '';
        var s = ('' + str).replace(/\$\$\$\$/g, '').replace(/[\r\n]+/g, ' ');
        if (s.length > maxLen) s = s.substring(0, maxLen);
        return s;
    }

    function sanitizeDataFieldValue(str, maxLen) {
        if (str === null || str === undefined) return '';
        var s = ('' + str).replace(/\$\$\$\$/g, '').replace(/[\r\n]+/g, ' ');
        if (s.length > maxLen) s = s.substring(0, maxLen);
        return s;
    }

    function padR(str, width) {
        while (str.length < width) str = str + ' ';
        return str;
    }

    function padF(num, width, decimals) {
        var s = num.toFixed(decimals);
        while (s.length < width) s = ' ' + s;
        return s;
    }

    function chargeToMol(charge) {
        if (charge === 0) return 0;
        if (charge === 1) return 3;
        if (charge === 2) return 2;
        if (charge === 3) return 1;
        if (charge === -1) return 5;
        if (charge === -2) return 6;
        if (charge === -3) return 7;
        return 0;
    }

    function timestamp() {
        var d = new Date();
        return pad(d.getMonth() + 1, 2).replace(/ /g, '0') +
               pad(d.getDate(), 2).replace(/ /g, '0') +
               pad(d.getFullYear() % 100, 2).replace(/ /g, '0') +
               pad(d.getHours(), 2).replace(/ /g, '0') +
               pad(d.getMinutes(), 2).replace(/ /g, '0');
    }

    global.MolfileWriter = { write: write, writeSDF: writeSDF, writeRXN: writeRXN };

})(typeof window !== 'undefined' ? window : this);
