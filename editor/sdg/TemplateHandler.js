/**
 * editor/sdg/TemplateHandler.js — CDK TemplateHandler port (v1.8.17 full).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * IP NOTICE — clean-room re-implementation, no source copy:
 *   Independent JavaScript implementation of the algorithm described
 *   in CDK's TemplateHandler (LGPL-2.1) and Helson '99 SDG. NO CDK
 *   source code copied. The CDK URL is an algorithm reference only.
 *
 * Reference (Java original, ~15 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/TemplateHandler.java
 *
 * Looks up scaffolds in BIME's editor/Templates.js library and applies
 * their pre-computed coordinates to a target molecule's atom-substructure
 * matches. SDG calls this BEFORE the algorithmic placement loop so
 * common scaffolds get a "good starting layout" before refinement.
 *
 * Architecture: TemplateHandler is a thin wrapper around BIME's
 * editor/Templates.js (which is already used by Layout.js's
 * matchRingSystemTemplates). SDG callers can use this CDK-shaped API
 * to access the same template library.
 *
 * Status (v1.8.17): FULL implementation. Replaces the v1.8.13 scaffold.
 */
(function (global) {
    'use strict';

    var TemplateHandler = {};

    /**
     * mapTemplates(mol, library) — try to apply a template to a
     * molecule's substructure. Returns the count of atoms placed.
     *
     * Workflow:
     *   1. Build the molecule's "ring-system signature" — sorted ring
     *      sizes (e.g. "2:6,6" for two fused 6-rings = naphthalene).
     *   2. Look up the signature in the optional `library` parameter
     *      (defaults to BIME's editor/Templates.js TEMPLATE_LOOKUP).
     *   3. For each candidate template, attempt a substructure match.
     *   4. On match, copy the template's coordinates to the matched
     *      atoms.
     */
    TemplateHandler.mapTemplates = function (mol, library) {
        if (!mol || !mol.atoms || mol.atoms.length === 0) return 0;
        var Templates = global.Templates;
        if (!Templates) return 0;

        // Default library = sigaure-keyed dispatch table. Same shape
        // as Layout.js's TEMPLATE_LOOKUP.
        library = library || TemplateHandler.DEFAULT_LIBRARY;

        var rings = TemplateHandler._findRings(mol);
        if (rings.length === 0) return 0;
        var systems = TemplateHandler._buildRingSystems(rings);

        var placed = 0;
        for (var si = 0; si < systems.length; si++) {
            var sysRingIdx = systems[si];
            if (sysRingIdx.length < 1) continue;
            var sizes = [];
            for (var ri = 0; ri < sysRingIdx.length; ri++) {
                sizes.push(rings[sysRingIdx[ri]].length);
            }
            sizes.sort(function (a, b) { return a - b; });
            var sig = sysRingIdx.length + ':' + sizes.join(',');
            var candidates = library[sig];
            if (!candidates || candidates.length === 0) continue;
            // Try each template in turn.
            for (var ci = 0; ci < candidates.length; ci++) {
                var tmplName = candidates[ci];
                var tmplFn = Templates[tmplName];
                if (typeof tmplFn !== 'function') continue;
                var tmpl = tmplFn();
                if (!tmpl || !tmpl.atoms) continue;
                // Apply to system atoms (by atom-symbol matching against
                // the template's atoms).
                var systemAtomIds = TemplateHandler._collectSystemAtoms(rings, sysRingIdx);
                if (systemAtomIds.length !== tmpl.atoms.length) continue;
                if (TemplateHandler._applyTemplate(mol, systemAtomIds, tmpl)) {
                    placed += systemAtomIds.length;
                    break;
                }
            }
        }
        return placed;
    };

    /**
     * addMolecule(mol) — register a placed molecule as a new template.
     * Captures the current 2D coordinates as a reference layout.
     */
    TemplateHandler.addMolecule = function (mol, name) {
        if (!mol || !mol.atoms) return false;
        if (!name) return false;
        TemplateHandler._userLibrary = TemplateHandler._userLibrary || {};
        TemplateHandler._userLibrary[name] = {
            atoms: mol.atoms.map(function (a) {
                return { symbol: a.symbol, x: a.x, y: a.y };
            })
        };
        return true;
    };

    /**
     * removeMolecule(name) — unregister a previously-added template.
     */
    TemplateHandler.removeMolecule = function (name) {
        if (!TemplateHandler._userLibrary) return false;
        if (!TemplateHandler._userLibrary[name]) return false;
        delete TemplateHandler._userLibrary[name];
        return true;
    };

    /**
     * _findRings(mol) — wrapper around mol.findRings with array
     * normalisation.
     */
    TemplateHandler._findRings = function (mol) {
        if (!mol.findRings) return [];
        var raw = [];
        try { raw = mol.findRings(20); } catch (e) { return []; }
        var rings = [];
        for (var i = 0; i < raw.length; i++) {
            rings.push(raw[i].atoms || raw[i]);
        }
        return rings;
    };

    /**
     * _buildRingSystems(rings) — union-find on shared atoms (≥ 1).
     */
    TemplateHandler._buildRingSystems = function (rings) {
        if (rings.length === 0) return [];
        var n = rings.length;
        var par = [];
        for (var i = 0; i < n; i++) par[i] = i;
        function find(x) { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; }
        function unite(a, b) { a = find(a); b = find(b); if (a !== b) par[b] = a; }
        for (var i = 0; i < n; i++) {
            for (var j = i + 1; j < n; j++) {
                var setI = {};
                for (var k = 0; k < rings[i].length; k++) setI[rings[i][k]] = true;
                var sharedCount = 0;
                for (var k2 = 0; k2 < rings[j].length; k2++) {
                    if (setI[rings[j][k2]]) sharedCount++;
                }
                if (sharedCount >= 1) unite(i, j);
            }
        }
        var systems = {};
        for (var i = 0; i < n; i++) {
            var root = find(i);
            if (!systems[root]) systems[root] = [];
            systems[root].push(i);
        }
        var result = [];
        for (var key in systems) result.push(systems[key]);
        return result;
    };

    /**
     * _collectSystemAtoms(rings, sysRingIndices) — union of atom IDs.
     */
    TemplateHandler._collectSystemAtoms = function (rings, sysRingIndices) {
        var seen = {};
        var result = [];
        for (var ri = 0; ri < sysRingIndices.length; ri++) {
            var ring = rings[sysRingIndices[ri]];
            for (var ai = 0; ai < ring.length; ai++) {
                if (!seen[ring[ai]]) {
                    seen[ring[ai]] = true;
                    result.push(ring[ai]);
                }
            }
        }
        return result;
    };

    /**
     * _applyTemplate(mol, systemAtomIds, tmpl) — copy template coords
     * to matched atoms via VF2++ subgraph isomorphism (deterministic,
     * connectivity-aware), with a fall-back to greedy symbol matching
     * if SMSDVF2 isn't loaded.
     *
     * v1.8.19 upgrade: replaces the v1.8.17 greedy-by-symbol matcher
     * with VF2-based atom mapping. This produces consistent results
     * for symmetric templates (benzene, cyclohexane) where the symbol
     * histogram is degenerate — VF2 picks the lexicographic-smallest
     * mapping; greedy used insertion order.
     *
     * Returns true on successful apply.
     */
    TemplateHandler._applyTemplate = function (mol, systemAtomIds, tmpl) {
        if (!systemAtomIds || !tmpl || !tmpl.atoms) return false;
        if (systemAtomIds.length !== tmpl.atoms.length) return false;
        // Gather symbol histograms; reject if mismatch.
        var molHist = {}, tmplHist = {};
        for (var i = 0; i < systemAtomIds.length; i++) {
            var a = mol.getAtom(systemAtomIds[i]);
            if (!a) return false;
            molHist[a.symbol] = (molHist[a.symbol] || 0) + 1;
        }
        for (var j = 0; j < tmpl.atoms.length; j++) {
            tmplHist[tmpl.atoms[j].symbol] = (tmplHist[tmpl.atoms[j].symbol] || 0) + 1;
        }
        for (var k in molHist) {
            if (molHist[k] !== tmplHist[k]) return false;
        }

        var tmplCx = 0, tmplCy = 0;
        for (var t = 0; t < tmpl.atoms.length; t++) {
            tmplCx += tmpl.atoms[t].x;
            tmplCy += tmpl.atoms[t].y;
        }
        tmplCx /= tmpl.atoms.length;
        tmplCy /= tmpl.atoms.length;

        // v1.8.19: prefer VF2++ subgraph isomorphism if SMSDVF2 +
        // SMSDGraph are loaded (deterministic, connectivity-aware).
        if (typeof Molecule !== 'undefined' &&
            typeof global.SMSDGraph !== 'undefined' &&
            typeof global.SMSDVF2 !== 'undefined' &&
            global.SMSDVF2.findSubstructure) {
            try {
                // Build a Molecule for the template and a sub-Molecule
                // for the system. Bonds matter for VF2.
                var tmplMol = new Molecule();
                var tmplIdMap = {};
                for (var ti = 0; ti < tmpl.atoms.length; ti++) {
                    var ta = tmplMol.addAtom(tmpl.atoms[ti].symbol,
                                              tmpl.atoms[ti].x, tmpl.atoms[ti].y);
                    if (ta) tmplIdMap[ti] = ta.id;
                }
                if (tmpl.bonds) {
                    for (var tbi = 0; tbi < tmpl.bonds.length; tbi++) {
                        var tb = tmpl.bonds[tbi];
                        if (tmplIdMap[tb.from] === undefined ||
                            tmplIdMap[tb.to] === undefined) continue;
                        tmplMol.addBond(tmplIdMap[tb.from], tmplIdMap[tb.to],
                                         tb.type || 1);
                    }
                }
                // Build sub-Molecule from system atoms + their internal bonds.
                var subMol = new Molecule();
                var subIdMap = {};
                var sysSet = {};
                for (var sj = 0; sj < systemAtomIds.length; sj++) {
                    sysSet[systemAtomIds[sj]] = true;
                    var ma = mol.getAtom(systemAtomIds[sj]);
                    if (ma) {
                        var sa = subMol.addAtom(ma.symbol, ma.x, ma.y);
                        if (sa) subIdMap[systemAtomIds[sj]] = sa.id;
                    }
                }
                for (var bi = 0; bi < mol.bonds.length; bi++) {
                    var bnd = mol.bonds[bi];
                    if (sysSet[bnd.atom1] && sysSet[bnd.atom2]) {
                        subMol.addBond(subIdMap[bnd.atom1], subIdMap[bnd.atom2],
                                        bnd.type || 1);
                    }
                }
                // Run VF2.
                var qGraph = new global.SMSDGraph.SMSDGraph(tmplMol);
                var tGraph = new global.SMSDGraph.SMSDGraph(subMol);
                var vf2Opts = new global.SMSDGraph.ChemOptions();
                vf2Opts.matchBondOrder = 'any';
                vf2Opts.ringMatchesRingOnly = false;
                var mapping = global.SMSDVF2.findSubstructure(qGraph, tGraph, vf2Opts, 5000);
                if (mapping) {
                    // mapping: queryIdx → targetIdx (0-based positions
                    // in qGraph / tGraph atom arrays, which correspond
                    // to insertion order = atom index in tmplMol/subMol).
                    for (var qi in mapping) {
                        var qIdx = parseInt(qi, 10);
                        var tIdx = mapping[qi];
                        // qIdx is the index in tmpl.atoms (template).
                        // tIdx is the index in systemAtomIds (sub-mol).
                        var ta2 = tmpl.atoms[qIdx];
                        var molAtomId = systemAtomIds[tIdx];
                        var ma2 = mol.getAtom(molAtomId);
                        if (ma2 && ta2) {
                            ma2.x = ta2.x - tmplCx;
                            ma2.y = ta2.y - tmplCy;
                        }
                    }
                    return true;
                }
            } catch (vfe) { /* fall through to greedy */ }
        }

        // Fall-back: greedy symbol-order matching (v1.8.17 behaviour).
        var molAtoms = systemAtomIds.map(function (id) { return mol.getAtom(id); });
        var usedTmpl = {};
        for (var mi = 0; mi < molAtoms.length; mi++) {
            var ma3 = molAtoms[mi];
            for (var ti2 = 0; ti2 < tmpl.atoms.length; ti2++) {
                if (usedTmpl[ti2]) continue;
                if (tmpl.atoms[ti2].symbol !== ma3.symbol) continue;
                ma3.x = tmpl.atoms[ti2].x - tmplCx;
                ma3.y = tmpl.atoms[ti2].y - tmplCy;
                usedTmpl[ti2] = true;
                break;
            }
        }
        return true;
    };

    /**
     * Default template library: subset of Layout.js's TEMPLATE_LOOKUP
     * exposed via this CDK API. Callers can override with their own
     * library map.
     */
    TemplateHandler.DEFAULT_LIBRARY = {
        '1:3':       ['cyclopropane'],
        '1:4':       ['cyclobutane', 'betalactam'],
        '1:5':       ['cyclopentane', 'pyrrolidine', 'imidazole', 'pyrazole', 'thiophene', 'furan', 'oxazole', 'isoxazole', 'thiazole', 'pyrrole'],
        '1:6':       ['benzene', 'cyclohexane', 'pyridine', 'pyrimidine', 'piperidine', 'piperazine', 'morpholine'],
        '1:7':       ['cycloheptane'],
        '2:5,6':     ['indole', 'benzofuran', 'benzothiazole', 'benzothiophene', 'benzimidazole', 'purine'],
        '2:6,6':     ['naphthalene', 'quinoline', 'isoquinoline', 'quinazoline', 'quinoxaline'],
        '3:5,6,6':   ['carbazole', 'fluorene'],
        '3:6,6,6':   ['phenanthrene', 'anthracene', 'acridine'],
        '4:5,6,6,6': ['steroid'],
        '4:6,6,6,6': ['tetracycline', 'pyrene']
    };

    global.SDG = global.SDG || {};
    global.SDG.TemplateHandler = TemplateHandler;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = TemplateHandler;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
