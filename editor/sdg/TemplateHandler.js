/**
 * editor/sdg/TemplateHandler.js — native template handler (v1.8.17 full).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Looks up scaffolds in BIME's editor/Templates.js library and applies
 * their pre-computed coordinates to a target molecule's atom-substructure
 * matches. SDG calls this BEFORE the algorithmic placement loop so
 * common scaffolds get a "good starting layout" before refinement.
 *
 * Architecture: TemplateHandler is a thin wrapper around BIME's
 * editor/Templates.js (which is already used by Layout.js's
 * matchRingSystemTemplates). SDG callers can use this API
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
     * fitTemplateToMapping(mol, mapping, tmpl, options)
     *
     * Computes the best rigid 2D transform from template coordinates to
     * existing molecule coordinates for mapped atom pairs. This is the 2D
     * Procrustes/Kabsch case: translate both point sets to their centroids,
     * then solve the least-squares rotation analytically.
     *
     * If the molecule coordinates are degenerate (common during initial SDG
     * placement, where all atoms start at 0,0), the transform falls back to
     * centered template coordinates at the target centroid.
     */
    TemplateHandler.fitTemplateToMapping = function (mol, mapping, tmpl, options) {
        options = options || {};
        if (!mol || !mapping || !tmpl || !tmpl.atoms) return null;

        var pairs = [];
        for (var ti in mapping) {
            if (!mapping.hasOwnProperty(ti)) continue;
            var tmplIdx = parseInt(ti, 10);
            var targetAtom = mol.getAtom(mapping[ti]);
            var sourceAtom = tmpl.atoms[tmplIdx];
            if (!targetAtom || !sourceAtom ||
                !isFinite(targetAtom.x) || !isFinite(targetAtom.y) ||
                !isFinite(sourceAtom.x) || !isFinite(sourceAtom.y)) {
                continue;
            }
            pairs.push({ source: sourceAtom, target: targetAtom });
        }
        if (pairs.length === 0) return null;

        var scx = 0, scy = 0, tcx = 0, tcy = 0;
        for (var i = 0; i < pairs.length; i++) {
            scx += pairs[i].source.x;
            scy += pairs[i].source.y;
            tcx += pairs[i].target.x;
            tcy += pairs[i].target.y;
        }
        scx /= pairs.length;
        scy /= pairs.length;
        tcx /= pairs.length;
        tcy /= pairs.length;

        var sourceVar = 0;
        var targetVar = 0;
        var dot = 0;
        var cross = 0;
        for (var j = 0; j < pairs.length; j++) {
            var sx = pairs[j].source.x - scx;
            var sy = pairs[j].source.y - scy;
            var tx = pairs[j].target.x - tcx;
            var ty = pairs[j].target.y - tcy;
            sourceVar += sx * sx + sy * sy;
            targetVar += tx * tx + ty * ty;
            dot += sx * tx + sy * ty;
            cross += sx * ty - sy * tx;
        }

        var minTargetSpread = options.minTargetSpread || 1e-6;
        var cos = 1;
        var sin = 0;
        var scale = 1;
        var degenerate = sourceVar <= 1e-9 || targetVar <= minTargetSpread;
        if (!degenerate) {
            var norm = Math.sqrt(dot * dot + cross * cross);
            if (norm > 1e-9) {
                cos = dot / norm;
                sin = cross / norm;
                if (options.allowScale === true) {
                    scale = norm / sourceVar;
                }
            }
        }

        return {
            sourceCx: scx,
            sourceCy: scy,
            targetCx: tcx,
            targetCy: tcy,
            cos: cos,
            sin: sin,
            scale: scale,
            degenerate: degenerate,
            pairCount: pairs.length
        };
    };

    TemplateHandler.transformTemplatePoint = function (point, transform) {
        var sx = (point.x - transform.sourceCx) * transform.scale;
        var sy = (point.y - transform.sourceCy) * transform.scale;
        return {
            x: transform.targetCx + sx * transform.cos - sy * transform.sin,
            y: transform.targetCy + sx * transform.sin + sy * transform.cos
        };
    };

    /**
     * applyTemplateBestFit(mol, mapping, tmpl, placed, options)
     *
     * Applies mapped template coordinates using fitTemplateToMapping. The
     * optional `placed` object is updated with every mapped molecule atom.
     */
    TemplateHandler.applyTemplateBestFit = function (mol, mapping, tmpl, placed, options) {
        var transform = TemplateHandler.fitTemplateToMapping(mol, mapping, tmpl, options);
        if (!transform) return null;
        for (var ti in mapping) {
            if (!mapping.hasOwnProperty(ti)) continue;
            var atom = mol.getAtom(mapping[ti]);
            var tmplAtom = tmpl.atoms[parseInt(ti, 10)];
            if (!atom || !tmplAtom) continue;
            var p = TemplateHandler.transformTemplatePoint(tmplAtom, transform);
            atom.x = p.x;
            atom.y = p.y;
            if (placed) placed[mapping[ti]] = true;
        }
        return transform;
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
     * to matched atoms via deterministic, connectivity-aware substructure
     * mapping, with a fall-back to greedy symbol matching if SMSDVF2 isn't
     * loaded.
     *
     * v1.8.19 upgrade: replaces the v1.8.17 greedy-by-symbol matcher
     * with graph-aware atom mapping. This produces consistent results
     * for symmetric templates (benzene, cyclohexane) where the symbol
     * histogram is degenerate — the mapper picks the lexicographic-smallest
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

        // v1.8.19: prefer deterministic substructure mapping if SMSDVF2 +
        // SMSDGraph are loaded.
        if (typeof Molecule !== 'undefined' &&
            typeof global.SMSDGraph !== 'undefined' &&
            typeof global.SMSDVF2 !== 'undefined' &&
            global.SMSDVF2.findSubstructure) {
            try {
                // Build a Molecule for the template and a sub-Molecule
                // for the system. Bonds matter for topology-aware mapping.
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
                        var from = tb.a1 !== undefined ? tb.a1 : (tb.from !== undefined ? tb.from : tb[0]);
                        var to = tb.a2 !== undefined ? tb.a2 : (tb.to !== undefined ? tb.to : tb[1]);
                        if (tmplIdMap[from] === undefined ||
                            tmplIdMap[to] === undefined) continue;
                        tmplMol.addBond(tmplIdMap[from], tmplIdMap[to],
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
                // Run matcher.
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
                    var tmplToMol = {};
                    for (var qi in mapping) {
                        var qIdx = parseInt(qi, 10);
                        var tIdx = mapping[qi];
                        tmplToMol[qIdx] = systemAtomIds[tIdx];
                    }
                    TemplateHandler.applyTemplateBestFit(mol, tmplToMol, tmpl, null);
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
                var p = TemplateHandler.transformTemplatePoint(tmpl.atoms[ti2], {
                    sourceCx: tmplCx,
                    sourceCy: tmplCy,
                    targetCx: 0,
                    targetCy: 0,
                    cos: 1,
                    sin: 0,
                    scale: 1
                });
                ma3.x = p.x;
                ma3.y = p.y;
                usedTmpl[ti2] = true;
                break;
            }
        }
        return true;
    };

    /**
     * Default template library: subset of Layout.js's TEMPLATE_LOOKUP
     * exposed via this API. Callers can override with their own
     * library map.
     */
    TemplateHandler.DEFAULT_LIBRARY = {
        '1:3':       ['cyclopropane'],
        '1:4':       ['cyclobutane', 'beta_lactam'],
        '1:5':       ['cyclopentane', 'pyrrolidine', 'imidazole', 'pyrazole', 'thiophene', 'furan', 'oxazole', 'isoxazole', 'thiazole', 'pyrrole'],
        '1:6':       ['benzene', 'cyclohexane', 'pyridine', 'pyrimidine', 'piperidine', 'piperazine', 'morpholine'],
        '1:7':       ['cycloheptane'],
        '2:5,6':     ['indole', 'benzofuran', 'benzothiazole', 'benzothiophene', 'benzimidazole', 'purine'],
        '2:6,6':     ['naphthalene', 'quinoline', 'isoquinoline', 'quinazoline', 'quinoxaline'],
        '3:5,6,6':   ['carbazole', 'fluorene'],
        '3:6,6,6':   ['phenanthrene', 'anthracene', 'acridine'],
        '4:5,6,6,6': ['steroid'],
        // v2.0.35: pyrene restored with correct peri-condensed template.
        '4:6,6,6,6': ['tetracycline', 'pyrene']
    };

    global.SDG = global.SDG || {};
    global.SDG.TemplateHandler = TemplateHandler;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = TemplateHandler;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
