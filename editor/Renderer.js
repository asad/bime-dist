/**
 * Renderer.js — SVG renderer for molecular structures
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Renders a Molecule graph as SVG elements: atom labels, bonds, reaction arrows,
 * atom-atom mapping numbers, selection highlights, and hover effects.
 */
(function(global) {
    'use strict';

    var ELEMENTS = Molecule.ELEMENTS;
    var BOND_LENGTH = Molecule.BOND_LENGTH;

    // Rendering constants
    var FONT_SIZE = 13;
    var FONT_FAMILY = '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    var FONT_MONO = 'SFMono-Regular, Consolas, "Liberation Mono", Menlo, monospace';
    var LABEL_PAD = 4;
    var DOUBLE_BOND_GAP = 3;
    var TRIPLE_BOND_GAP = 4;
    var WEDGE_WIDTH = 5;
    var SELECT_COLOR = '#0d9488';
    var HIGHLIGHT_COLOR = '#99f6e4';
    var MAP_FONT_SIZE = 12;
    var MAP_COLOR = '#0d9488';
    var HYDROGEN_FONT_SIZE = 10;
    var ARROW_HEAD_SIZE = 12;
    var ARROW_STROKE_WIDTH = 2.2;
    var PLUS_FONT_SIZE = 18;
    var PLUS_COLOR = '#555';
    var CIP_FONT_SIZE = 9;
    var CIP_COLOR = '#888';

    // v1.4.0: default palette for mol-mol component pair highlighting.
    // Pair 1 reuses the BIME teal brand colour; remaining seven cycle through
    // colour-blind-friendly hues. Rotates if a reaction has more pairs.
    var PAIR_PALETTE = [
        '#0d9488', // teal (BIME brand)
        '#a855f7', // purple
        '#f97316', // orange
        '#06b6d4', // cyan
        '#84cc16', // lime
        '#ec4899', // pink
        '#eab308', // yellow
        '#3b82f6'  // blue
    ];

    // v1.4.0: RDT-style per-atom soft halo palette. Ten distinct pale
    // tints suitable as background halos behind atom symbols (CPK colours
    // remain intact, halo is purely additive). Cycled by paletteIndex when
    // a reaction has more component pairs than the palette length.
    var COMPONENT_PAIR_PALETTE = [
        '#a5d8d2',  // pale teal
        '#c7e8b9',  // pale green
        '#fbd5a5',  // pale orange
        '#d8c5e7',  // pale purple
        '#a8d8e8',  // pale cyan
        '#fbc4c1',  // pale pink
        '#e8d8a5',  // pale yellow
        '#c5e2c7',  // pale mint
        '#e5b8d8',  // pale rose
        '#b8c5e5'   // pale lavender
    ];
    // Neutral grey halo for unpaired components (paletteIndex === -1).
    var COMPONENT_PAIR_NEUTRAL = '#cccccc';

    function Renderer(container, width, height) {
        this.container = typeof container === 'string' ? document.getElementById(container) : container;
        this.width = parseInt(width) || 500;
        this.height = parseInt(height) || 400;
        this.molecule = null;
        this.showHydrogens = true;
        this.showMapNumbers = true;
        this.showMolName = false;    // only show name when explicitly set
        this.depictMode = false;
        this.showCipLabels = false;  // legacy: show both R/S and E/Z
        this.showCipRS = false;     // show R/S at stereocentres
        this.showCipEZ = false;     // show E/Z at double bonds
        this.aromaticStyle = 'circle'; // 'circle' or 'kekule'
        // v1.4.0: optional component-pair highlight overlay. When non-empty,
        // each entry { reactantAtomIds, productAtomIds, color, pairIndex,
        // paletteIndex } is rendered as RDT-style soft per-atom halos
        // behind atoms/bonds. Populated by MolEditor after RDT.mapReaction.
        // Cleared on user "Clear" or new molecule load.
        this.componentPairs = [];
        this.showComponentPairs = true;
        // v1.4.0: master switch for per-atom halo rendering ("Color atoms"
        // checkbox). Default ON: halos paint when componentPairs is set.
        // Setting to false suppresses halos but leaves componentPairs and
        // map numbers intact.
        this.colorAtoms = true;

        // Viewport transform (pan + zoom)
        this.offsetX = 0;
        this.offsetY = 0;
        this.scale = 1;

        // Create SVG element
        this.svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        this.svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
        this.svg.setAttribute('width', '100%');
        this.svg.setAttribute('height', '100%');
        this.svg.setAttribute('viewBox', '0 0 ' + this.width + ' ' + this.height);
        this.svg.style.display = 'block';
        this.svg.style.background = 'var(--color-bg, white)';
        this.svg.style.cursor = 'crosshair';
        this.svg.style.userSelect = 'none';

        // Accessibility: expose this SVG as an image to assistive technology
        // and provide a <title>/<desc> updated on every render(). WCAG 1.1.1.
        this.svg.setAttribute('role', 'img');
        this.svg.setAttribute('aria-labelledby', 'bime-svg-title bime-svg-desc');
        this._a11yTitle = document.createElementNS('http://www.w3.org/2000/svg', 'title');
        this._a11yTitle.setAttribute('id', 'bime-svg-title');
        this._a11yTitle.textContent = 'Empty molecule canvas';
        this._a11yDesc = document.createElementNS('http://www.w3.org/2000/svg', 'desc');
        this._a11yDesc.setAttribute('id', 'bime-svg-desc');
        this._a11yDesc.textContent = 'Use the toolbar above to draw a molecule.';
        this.svg.appendChild(this._a11yTitle);
        this.svg.appendChild(this._a11yDesc);

        // Layers (back to front)
        this.bgLayer = this._createGroup('bg-layer');
        this.bondLayer = this._createGroup('bond-layer');
        this.atomLayer = this._createGroup('atom-layer');
        this.overlayLayer = this._createGroup('overlay-layer');

        this.svg.appendChild(this.bgLayer);
        this.svg.appendChild(this.bondLayer);
        this.svg.appendChild(this.atomLayer);
        this.svg.appendChild(this.overlayLayer);

        if (this.container) {
            this.container.appendChild(this.svg);
        }
    }

    Renderer.prototype._createGroup = function(id) {
        var g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        g.setAttribute('id', id);
        return g;
    };

    Renderer.prototype.setMolecule = function(mol) {
        this.molecule = mol;
    };

    Renderer.prototype.setSize = function(w, h) {
        this.width = parseInt(w) || this.width;
        this.height = parseInt(h) || this.height;
        this.svg.setAttribute('viewBox', '0 0 ' + this.width + ' ' + this.height);
    };

    // =========================================================================
    // Main render
    // =========================================================================
    Renderer.prototype.render = function() {
        this._clear(this.bgLayer);
        this._clear(this.bondLayer);
        this._clear(this.atomLayer);
        this._clear(this.overlayLayer);

        if (!this.molecule) return;

        var mol = this.molecule;

        // Cache ring info for bond rendering decisions
        // FIX: also check that getRingInfo exists to avoid TypeError if Layout is partially loaded
        this._ringInfo = (typeof Layout !== 'undefined' && Layout.getRingInfo) ? Layout.getRingInfo(mol) : [];

        // v1.4.0: draw mol-mol soft per-atom halos BEHIND bonds (in bgLayer)
        // when set by MolEditor after Auto-map. The "Color atoms" UI toggle
        // sets `colorAtoms` to false to suppress without losing pair data.
        if (this.colorAtoms !== false && this.showComponentPairs &&
            this.componentPairs && this.componentPairs.length > 0) {
            this._drawComponentPairs();
        }

        // Draw aromatic ring circles (behind bonds) — only in 'circle' style
        if (this.aromaticStyle !== 'kekule') {
            for (var ri = 0; ri < this._ringInfo.length; ri++) {
                var ring = this._ringInfo[ri];
                if (ring.aromatic) {
                    this._drawAromaticCircle(ring);
                }
            }
        }

        // Draw bonds (behind atoms)
        for (var i = 0; i < mol.bonds.length; i++) {
            this._drawBond(mol.bonds[i]);
        }

        // Draw atoms
        for (var i = 0; i < mol.atoms.length; i++) {
            this._drawAtom(mol.atoms[i]);
        }

        // Draw CIP stereochemistry labels
        var showRS = this.showCipRS || this.showCipLabels;
        var showEZ = this.showCipEZ || this.showCipLabels;
        if ((showRS || showEZ) && typeof CipStereo !== 'undefined') {
            if (showRS) CipStereo.assignRS(mol);
            if (showEZ) CipStereo.assignEZ(mol);
            this._drawCipLabels(mol, showRS, showEZ);
        }

        // Draw reaction arrow
        if (mol.reactionArrow) {
            this._drawReactionArrow(mol.reactionArrow);
        }

        // Draw "+" signs between reaction components
        if (mol.reactionPlusSigns && mol.reactionPlusSigns.length > 0) {
            for (var pi = 0; pi < mol.reactionPlusSigns.length; pi++) {
                this._drawPlusSign(mol.reactionPlusSigns[pi]);
            }
        }

        // Draw atom-atom mapping connection lines between highlighted map partners
        this._drawMapConnections(mol);

        // Draw molecule name below structure
        if (this.showMolName && mol.name) {
            this._drawMolName(mol);
        }

        // Update accessible name/description for screen readers (WCAG 1.1.1).
        this._updateA11yLabels(mol);
    };

    /**
     * Update the SVG <title> and <desc> so screen-reader users get a
     * meaningful summary of the current molecule. Called from render().
     */
    Renderer.prototype._updateA11yLabels = function(mol) {
        if (!this._a11yTitle || !this._a11yDesc) return;
        var atomCount = mol.atoms ? mol.atoms.length : 0;
        var bondCount = mol.bonds ? mol.bonds.length : 0;
        var name = (mol.name || '').trim();
        var isReaction = !!(mol.reactionArrow);
        var titleStr;
        if (atomCount === 0) {
            titleStr = 'Empty molecule canvas';
        } else if (name) {
            titleStr = (isReaction ? 'Reaction: ' : 'Molecule: ') + name;
        } else {
            titleStr = isReaction
                ? 'Drawn reaction with ' + atomCount + ' atoms'
                : 'Drawn molecule with ' + atomCount + ' atoms';
        }
        var descStr = atomCount === 0
            ? 'Use the toolbar above to draw a molecule, or load one from the library.'
            : ('Structure has ' + atomCount + ' atom' + (atomCount === 1 ? '' : 's') +
               ' and ' + bondCount + ' bond' + (bondCount === 1 ? '' : 's') + '.');
        this._a11yTitle.textContent = titleStr;
        this._a11yDesc.textContent = descStr;
    };

    /**
     * Draw a dashed inner circle to represent aromaticity inside a ring.
     */
    Renderer.prototype._drawAromaticCircle = function(ringInfo) {
        var cx = this._tx(ringInfo.center.x);
        var cy = this._ty(ringInfo.center.y);

        // Compute average distance from center to ring atoms (= polygon radius)
        var mol = this.molecule;
        // FIX: ringInfo may carry stale atom IDs after a delete; guard against
        // a missing atom rather than NPE-ing on a.x / a.y.
        var totalDist = 0;
        var counted = 0;
        for (var i = 0; i < ringInfo.atoms.length; i++) {
            var a = mol.getAtom(ringInfo.atoms[i]);
            if (!a) continue;
            var dx = this._tx(a.x) - cx;
            var dy = this._ty(a.y) - cy;
            totalDist += Math.sqrt(dx * dx + dy * dy);
            counted++;
        }
        if (counted === 0) return; // nothing to render
        var avgRadius = totalDist / counted;
        // Inner circle at ~60% of radius
        var innerRadius = avgRadius * 0.6;

        var circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        circle.setAttribute('cx', cx);
        circle.setAttribute('cy', cy);
        circle.setAttribute('r', innerRadius);
        circle.setAttribute('fill', 'none');
        circle.setAttribute('stroke', '#999');
        circle.setAttribute('stroke-width', '1');
        circle.setAttribute('stroke-dasharray', '3,2');
        this.bondLayer.appendChild(circle);
    };

    // =========================================================================
    // v1.4.0: Mol-mol (component-component) per-atom halo highlight
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    //
    // RDT-inspired per-atom soft halos: every atom in a paired reactant
    // component gets a light coloured circular halo behind its symbol, and
    // every atom in the matching product component gets a halo in the same
    // colour. Atoms in unpaired components (paletteIndex === -1) get a
    // neutral grey halo. The halos sit on bgLayer so they render BENEATH
    // bonds/atom labels and survive ImageExport.toSVG/toPNG/toPDF
    // (export paths serialise the live DOM).
    // =========================================================================
    /**
     * Draw soft circular halos behind every atom in every component pair.
     * Each entry of `this.componentPairs` is:
     *   { reactantAtomIds: [..], productAtomIds: [..], color: '#...',
     *     pairIndex: N, paletteIndex: N }
     * paletteIndex === -1 marks unpaired components (neutral grey halo).
     */
    Renderer.prototype._drawComponentPairs = function() {
        var mol = this.molecule;
        if (!mol) return;
        var HALO_RADIUS = 14;
        var HALO_OPACITY = 0.5;
        for (var pi = 0; pi < this.componentPairs.length; pi++) {
            var pair = this.componentPairs[pi];
            if (!pair) continue;
            var color = this._haloColorForPair(pair, pi);
            this._drawAtomHalos(pair.reactantAtomIds, color, HALO_RADIUS, HALO_OPACITY, pair.paletteIndex);
            this._drawAtomHalos(pair.productAtomIds,  color, HALO_RADIUS, HALO_OPACITY, pair.paletteIndex);
        }
    };

    /**
     * Resolve the halo fill colour for a component-pair entry. Honours an
     * explicit `pair.color` override, then falls back to the palette indexed
     * by paletteIndex. paletteIndex === -1 (unpaired components) yields the
     * neutral grey halo. paletteIndex >= palette.length wraps modulo.
     */
    Renderer.prototype._haloColorForPair = function(pair, fallbackIdx) {
        if (pair && pair.color) return pair.color;
        var idx = (pair && typeof pair.paletteIndex === 'number') ? pair.paletteIndex : fallbackIdx;
        if (idx === undefined || idx === null) idx = fallbackIdx;
        if (idx < 0) return COMPONENT_PAIR_NEUTRAL;
        return COMPONENT_PAIR_PALETTE[idx % COMPONENT_PAIR_PALETTE.length];
    };

    /**
     * Append one halo <circle> per atom id behind the atom symbol.
     * Halos go in bgLayer so bonds and labels render on top.
     */
    Renderer.prototype._drawAtomHalos = function(atomIds, color, radius, opacity, paletteIndex) {
        if (!atomIds || atomIds.length === 0) return;
        var mol = this.molecule;
        for (var i = 0; i < atomIds.length; i++) {
            var atom = mol.getAtom(atomIds[i]);
            if (!atom) continue;
            var cx = this._tx(atom.x);
            var cy = this._ty(atom.y);
            var halo = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            halo.setAttribute('cx', cx);
            halo.setAttribute('cy', cy);
            halo.setAttribute('r', radius);
            halo.setAttribute('fill', color);
            halo.setAttribute('opacity', String(opacity));
            halo.setAttribute('class', 'bime-pair-halo');
            halo.setAttribute('data-atom-id', String(atom.id));
            if (paletteIndex !== undefined && paletteIndex !== null) {
                halo.setAttribute('data-palette-index', String(paletteIndex));
            }
            this.bgLayer.appendChild(halo);
        }
    };

    // =========================================================================
    // v1.5.2: Search-match halo (library Search panel hit click)
    // -------------------------------------------------------------------------
    // The Search panel calls MolEditor.highlightSearchMatch(mode, querySmiles)
    // which flags every matched atom/bond with `searchMatchKind` plus a
    // timestamp `searchMatchTs`. The renderer turns those flags into a
    // polished two-tier halo:
    //
    //   - outer soft pale-teal disc/underlay (r=11 / line w=7, ~0.32 opacity)
    //   - inner brand-teal ring outline / inner stroke (~0.85 opacity)
    //
    // On first paint we also fire a one-shot SMIL reveal animation that
    // brightens both layers from the resting state to a "find result"-style
    // flash (#5eead4) and settles back over ~1.0 s. The animation is gated
    // on `prefers-reduced-motion: reduce`: when a user has that preference
    // set, we skip the <animate> children and paint the resting state.
    //
    // BRAND_COLOR `#0d9488` is intentionally re-used as the inner stroke
    // for WCAG 2.1 AA against `#fff` (~4.96:1) and against `#0f172a` dark
    // canvas (~5.36:1). The pale fill `#99f6e4` is purely decorative —
    // atom symbols are still rendered on top in their CPK colour.
    // =========================================================================
    var SEARCH_MATCH_FILL = '#99f6e4';     // pale teal halo fill (decorative)
    var SEARCH_MATCH_RING = '#0d9488';     // BIME brand teal — WCAG-AA stroke
    var SEARCH_MATCH_FLASH = '#5eead4';    // brighter teal for reveal flash
    var SEARCH_MATCH_REVEAL_MS = 1000;     // total reveal duration
    var SEARCH_MATCH_REVEAL_WINDOW_MS = 250; // only animate hits painted within
                                              //  this window of the flag ts so
                                              //  pan/zoom redraws don't replay

    /**
     * Detect `prefers-reduced-motion: reduce`. Cached per-renderer; safe in
     * non-browser test environments (returns false).
     */
    Renderer.prototype._prefersReducedMotion = function() {
        if (typeof window === 'undefined' || !window.matchMedia) return false;
        try {
            return window.matchMedia('(prefers-reduced-motion: reduce)').matches === true;
        } catch (e) { return false; }
    };

    /**
     * True when a freshly-flagged hit (matchTs within REVEAL_WINDOW_MS of
     * Date.now()) should play the reveal animation. Returns false in tests
     * (no Date.now branch) so unit tests don't depend on timing.
     */
    Renderer.prototype._shouldAnimateSearchMatch = function(matchTs) {
        if (this._prefersReducedMotion()) return false;
        if (!matchTs) return false;
        if (typeof Date === 'undefined' || !Date.now) return false;
        var dt = Date.now() - matchTs;
        return dt >= 0 && dt <= SEARCH_MATCH_REVEAL_WINDOW_MS;
    };

    /**
     * Append a SMIL <animate> child that brightens `attrName` from
     * `flashValue` to `restValue` over SEARCH_MATCH_REVEAL_MS, using
     * fill="freeze" so the resting state is what stays after the reveal.
     */
    Renderer.prototype._appendRevealAnim = function(parent, attrName, flashValue, restValue) {
        var anim = document.createElementNS('http://www.w3.org/2000/svg', 'animate');
        anim.setAttribute('attributeName', attrName);
        anim.setAttribute('from', flashValue);
        anim.setAttribute('to', restValue);
        anim.setAttribute('dur', (SEARCH_MATCH_REVEAL_MS / 1000) + 's');
        anim.setAttribute('begin', '0s');
        anim.setAttribute('fill', 'freeze');
        anim.setAttribute('calcMode', 'spline');
        anim.setAttribute('keyTimes', '0;1');
        anim.setAttribute('keySplines', '0.16 1 0.3 1'); // ease-out (cubic)
        parent.appendChild(anim);
    };

    /**
     * Append a SVG <title> child so screen readers announce matched
     * atoms/bonds as "Search match (substructure)" etc.
     */
    Renderer.prototype._appendMatchA11yTitle = function(parent, kind) {
        var t = document.createElementNS('http://www.w3.org/2000/svg', 'title');
        t.textContent = 'Search match' + (kind ? ' (' + kind + ')' : '');
        parent.appendChild(t);
    };

    /**
     * Paint the polished v1.5.2 atom search-match halo at (cx, cy).
     * Layers (back→front, all in bgLayer):
     *   1. r=11 disc filled SEARCH_MATCH_FILL @ ~0.32 opacity
     *   2. r=12 ring stroked SEARCH_MATCH_RING @ ~0.85 opacity, width 1.6
     * If first-paint, both layers also carry a SMIL <animate> that flashes
     * to SEARCH_MATCH_FLASH and settles back. Honours reduced-motion.
     */
    Renderer.prototype._drawSearchMatchAtomHalo = function(cx, cy, kind, matchTs) {
        var animate = this._shouldAnimateSearchMatch(matchTs);

        // Outer soft pale-teal disc.
        var disc = this._circle(cx, cy, 11, SEARCH_MATCH_FILL, 'none', 0.32);
        disc.setAttribute('class', 'bime-search-match-fill');
        disc.setAttribute('data-match-kind', String(kind || ''));
        if (animate) {
            this._appendRevealAnim(disc, 'r', '4', '11');
            this._appendRevealAnim(disc, 'opacity', '0.7', '0.32');
        }
        this._appendMatchA11yTitle(disc, kind);
        this.bgLayer.appendChild(disc);

        // Inner brand-teal ring — sharper, gives "circled" cue.
        var ring = this._circle(cx, cy, 12, 'none', SEARCH_MATCH_RING, 0.85);
        ring.setAttribute('stroke-width', '1.6');
        ring.setAttribute('class', 'bime-search-match-ring');
        ring.setAttribute('data-match-kind', String(kind || ''));
        if (animate) {
            this._appendRevealAnim(ring, 'stroke', SEARCH_MATCH_FLASH, SEARCH_MATCH_RING);
            this._appendRevealAnim(ring, 'stroke-width', '3.2', '1.6');
        }
        this.bgLayer.appendChild(ring);
    };

    /**
     * Paint the v1.5.2 bond search-match glow between (x1,y1)→(x2,y2):
     *   1. soft pale-teal underlay line  (width 7, opacity 0.32, round cap)
     *   2. brighter brand-teal inner line (width 2.4, opacity 0.7, round cap)
     * On first paint, fires a one-shot SMIL reveal that flashes to
     * SEARCH_MATCH_FLASH and settles. Honours reduced-motion.
     */
    Renderer.prototype._drawSearchMatchBondGlow = function(x1, y1, x2, y2, kind, matchTs) {
        var animate = this._shouldAnimateSearchMatch(matchTs);

        var under = this._line(x1, y1, x2, y2, SEARCH_MATCH_FILL, 7, 'round');
        under.style.opacity = '0.32';
        under.setAttribute('class', 'bime-search-match-bond-under');
        under.setAttribute('data-match-kind', String(kind || ''));
        if (animate) {
            this._appendRevealAnim(under, 'opacity', '0.7', '0.32');
            this._appendRevealAnim(under, 'stroke-width', '11', '7');
        }
        this._appendMatchA11yTitle(under, kind);
        this.bgLayer.appendChild(under);

        var inner = this._line(x1, y1, x2, y2, SEARCH_MATCH_RING, 2.4, 'round');
        inner.style.opacity = '0.7';
        inner.setAttribute('class', 'bime-search-match-bond-inner');
        inner.setAttribute('data-match-kind', String(kind || ''));
        if (animate) {
            this._appendRevealAnim(inner, 'stroke', SEARCH_MATCH_FLASH, SEARCH_MATCH_RING);
            this._appendRevealAnim(inner, 'stroke-width', '4', '2.4');
        }
        this.bgLayer.appendChild(inner);
    };

    /**
     * v1.4.0: helper — return a flat array of colour strings, one per
     * componentPairs entry plus a trailing neutral-grey for paletteIndex
     * === -1. Useful for tests and external palette inspection.
     *
     * v1.4.1 (bug-fix #6 documentation): the returned array is always
     * `componentPairs.length + 1` long: the leading N entries map the input
     * list, the trailing entry is the neutral grey halo used for
     * paletteIndex = -1 (unpaired components / spectator molecules). This
     * asymmetry is deliberate — callers iterating the list pair-wise should
     * stop at `out.length - 1` if they only want pair-colours, or use
     * `out[out.length - 1]` to look up the unpaired-component colour.
     */
    Renderer.componentPairColors = function(componentPairs) {
        var out = [];
        if (!componentPairs) return out;
        for (var i = 0; i < componentPairs.length; i++) {
            var p = componentPairs[i];
            var idx = (p && typeof p.paletteIndex === 'number') ? p.paletteIndex : i;
            if (idx < 0) {
                out.push(COMPONENT_PAIR_NEUTRAL);
            } else {
                out.push(COMPONENT_PAIR_PALETTE[idx % COMPONENT_PAIR_PALETTE.length]);
            }
        }
        out.push(COMPONENT_PAIR_NEUTRAL);
        return out;
    };

    // Expose the palettes so tests / external callers can inspect defaults.
    Renderer.PAIR_PALETTE = PAIR_PALETTE.slice();
    Renderer.COMPONENT_PAIR_PALETTE = COMPONENT_PAIR_PALETTE.slice();
    Renderer.COMPONENT_PAIR_NEUTRAL = COMPONENT_PAIR_NEUTRAL;
    // v1.5.2: expose the search-match halo palette for tests / theming.
    Renderer.SEARCH_MATCH_FILL = SEARCH_MATCH_FILL;
    Renderer.SEARCH_MATCH_RING = SEARCH_MATCH_RING;
    Renderer.SEARCH_MATCH_FLASH = SEARCH_MATCH_FLASH;
    Renderer.SEARCH_MATCH_REVEAL_MS = SEARCH_MATCH_REVEAL_MS;

    // =========================================================================
    // Bond rendering
    // =========================================================================
    Renderer.prototype._drawBond = function(bond) {
        var mol = this.molecule;
        var a1 = mol.getAtom(bond.atom1);
        var a2 = mol.getAtom(bond.atom2);
        if (!a1 || !a2) return;

        var x1 = this._tx(a1.x), y1 = this._ty(a1.y);
        var x2 = this._tx(a2.x), y2 = this._ty(a2.y);

        // Shorten bonds to avoid overlapping atom labels
        var dx = x2 - x1, dy = y2 - y1;
        var len = Math.sqrt(dx * dx + dy * dy);
        if (len < 1) return;
        var nx = dx / len, ny = dy / len;

        var trim1 = this._labelRadius(a1);
        var trim2 = this._labelRadius(a2);
        x1 += nx * trim1; y1 += ny * trim1;
        x2 -= nx * trim2; y2 -= ny * trim2;

        var color = bond.highlighted ? HIGHLIGHT_COLOR : (bond.selected ? SELECT_COLOR : '#333');
        var strokeWidth = bond.selected ? 2.5 : 1.8;

        // Background highlight.
        // v1.5.2: when the bond was flagged by `highlightSearchMatch` (i.e.
        // bond.searchMatchKind is set), paint a richer two-tier glow:
        //   1. a soft pale-teal underlay (width 7, opacity ~0.32)
        //   2. a brighter brand-teal inner stroke (width 2.5, opacity ~0.7)
        // and (if the user does not prefer-reduced-motion) fire a one-shot
        // 1.0 s SMIL reveal that brightens both layers from the centre out.
        // Legacy callers that set bond.bgColor without searchMatchKind keep
        // the prior single-layer flood.
        if (bond.bgColor) {
            var smk = bond.searchMatchKind;
            if (smk) {
                this._drawSearchMatchBondGlow(x1, y1, x2, y2, smk, bond.searchMatchTs);
            } else {
                var bg = this._line(x1, y1, x2, y2, bond.bgColor, 6, 'round');
                bg.style.opacity = '0.4';
                this.bgLayer.appendChild(bg);
            }
        }

        // Query bond: draw as dashed line with teal colour
        if (bond.isQuery || bond.queryType) {
            var qColor = '#0d9488';
            var qLine = this._line(x1, y1, x2, y2, qColor, 1.8);
            qLine.setAttribute('stroke-dasharray', '5,3');
            qLine.style.opacity = '0.8';
            this.bondLayer.appendChild(qLine);
            // Label for non-default query bond types
            if (bond.queryType && bond.queryType !== 'default' && bond.queryType !== 'single') {
                var midX = (x1 + x2) / 2;
                var midY = (y1 + y2) / 2 - 6;
                var labelMap = {
                    'any': '~', 'double': '=', 'triple': '#', 'aromatic': ':',
                    'not_single': '!-', 'not_double': '!=', 'not_triple': '!#'
                };
                var lbl = labelMap[bond.queryType] || '';
                if (lbl) {
                    var txt = this._text(midX, midY, lbl, qColor, 9);
                    txt.setAttribute('font-family', 'monospace');
                    this.bondLayer.appendChild(txt);
                }
            }
            // Hit target
            var hit = this._line(x1, y1, x2, y2, 'transparent', 10);
            hit.setAttribute('data-bond-id', bond.id);
            hit.style.cursor = 'pointer';
            this.overlayLayer.appendChild(hit);
            return; // skip normal bond rendering
        }

        if (bond.stereo === Molecule.STEREO_WEDGE) {
            this._drawWedgeBond(x1, y1, x2, y2, color, this.bondLayer);
        } else if (bond.stereo === Molecule.STEREO_DASH) {
            this._drawDashBond(x1, y1, x2, y2, color, this.bondLayer);
        } else if (bond.type === Molecule.BOND_SINGLE) {
            this.bondLayer.appendChild(this._line(x1, y1, x2, y2, color, strokeWidth));
        } else if (bond.type === Molecule.BOND_DOUBLE) {
            this._drawDoubleBond(bond, x1, y1, x2, y2, nx, ny, color, strokeWidth);
        } else if (bond.type === Molecule.BOND_TRIPLE) {
            var px = -ny * TRIPLE_BOND_GAP, py = nx * TRIPLE_BOND_GAP;
            this.bondLayer.appendChild(this._line(x1, y1, x2, y2, color, strokeWidth));
            this.bondLayer.appendChild(this._line(x1 + px, y1 + py, x2 + px, y2 + py, color, strokeWidth * 0.7));
            this.bondLayer.appendChild(this._line(x1 - px, y1 - py, x2 - px, y2 - py, color, strokeWidth * 0.7));
        }

        // Hit target (invisible wider line for easier clicking)
        var hit = this._line(x1, y1, x2, y2, 'transparent', 10);
        hit.setAttribute('data-bond-id', bond.id);
        hit.style.cursor = 'pointer';
        this.overlayLayer.appendChild(hit);
    };

    /**
     * Draw a double bond with ring-aware offset.
     * In-ring double bonds: one line centered, second offset toward ring center.
     * Chain double bonds: symmetric (centered) double bond.
     */
    Renderer.prototype._drawDoubleBond = function(bond, x1, y1, x2, y2, nx, ny, color, strokeWidth) {
        var ringInfo = this._findRingForBond(bond);

        if (ringInfo) {
            // Offset double bond: main bond on center line, second bond toward ring center
            var rcx = this._tx(ringInfo.center.x);
            var rcy = this._ty(ringInfo.center.y);
            var mx = (x1 + x2) / 2;
            var my = (y1 + y2) / 2;
            // Direction from bond midpoint toward ring center
            var toCx = rcx - mx;
            var toCy = rcy - my;
            // Project onto perpendicular direction
            var perpX = -ny;
            var perpY = nx;
            var dot = toCx * perpX + toCy * perpY;
            var sign = dot >= 0 ? 1 : -1;

            var offsetX = perpX * DOUBLE_BOND_GAP * 1.5 * sign;
            var offsetY = perpY * DOUBLE_BOND_GAP * 1.5 * sign;

            // Main bond
            this.bondLayer.appendChild(this._line(x1, y1, x2, y2, color, strokeWidth));
            // Shorter offset bond (trimmed by 15% on each end)
            var trim = 0.15;
            var bdx = x2 - x1, bdy = y2 - y1;
            var blen = Math.sqrt(bdx * bdx + bdy * bdy);
            var sx = x1 + offsetX + nx * blen * trim;
            var sy = y1 + offsetY + ny * blen * trim;
            var ex = x2 + offsetX - nx * blen * trim;
            var ey = y2 + offsetY - ny * blen * trim;
            this.bondLayer.appendChild(this._line(sx, sy, ex, ey, color, strokeWidth * 0.8));
        } else {
            // Symmetric double bond (chain)
            var px = -ny * DOUBLE_BOND_GAP, py = nx * DOUBLE_BOND_GAP;
            this.bondLayer.appendChild(this._line(x1 + px, y1 + py, x2 + px, y2 + py, color, strokeWidth));
            this.bondLayer.appendChild(this._line(x1 - px, y1 - py, x2 - px, y2 - py, color, strokeWidth));
        }
    };

    /**
     * Find the ring (if any) that contains this bond.
     */
    Renderer.prototype._findRingForBond = function(bond) {
        if (!this._ringInfo) return null;
        for (var i = 0; i < this._ringInfo.length; i++) {
            var ring = this._ringInfo[i];
            var atoms = ring.atoms;
            var idx1 = atoms.indexOf(bond.atom1);
            var idx2 = atoms.indexOf(bond.atom2);
            if (idx1 >= 0 && idx2 >= 0) {
                var diff = Math.abs(idx1 - idx2);
                if (diff === 1 || diff === atoms.length - 1) {
                    return ring;
                }
            }
        }
        return null;
    };

    Renderer.prototype._drawWedgeBond = function(x1, y1, x2, y2, color, layer) {
        var dx = x2 - x1, dy = y2 - y1;
        var len = Math.sqrt(dx * dx + dy * dy);
        if (len < 0.5) return;
        var nx = -dy / len, ny = dx / len;
        var w = WEDGE_WIDTH;

        // Smooth filled wedge with gradient-like appearance
        var points = x1 + ',' + y1 + ' ' +
            (x2 + nx * w) + ',' + (y2 + ny * w) + ' ' +
            (x2 - nx * w) + ',' + (y2 - ny * w);
        var poly = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
        poly.setAttribute('points', points);
        poly.setAttribute('fill', color);
        poly.setAttribute('stroke', color);
        poly.setAttribute('stroke-width', '0.5');
        poly.setAttribute('stroke-linejoin', 'round');
        layer.appendChild(poly);
    };

    Renderer.prototype._drawDashBond = function(x1, y1, x2, y2, color, layer) {
        var dx = x2 - x1, dy = y2 - y1;
        var len = Math.sqrt(dx * dx + dy * dy);
        if (len < 0.5) return;
        var nx = -dy / len, ny = dx / len;
        var segments = 7;
        // Draw tapered dashes from narrow (start) to wide (end)
        for (var i = 1; i <= segments; i++) {
            var t = i / segments;
            var cx = x1 + dx * t, cy = y1 + dy * t;
            var w = WEDGE_WIDTH * t;
            var dash = this._line(cx - nx * w, cy - ny * w, cx + nx * w, cy + ny * w, color, 1.5);
            dash.setAttribute('stroke-linecap', 'round');
            layer.appendChild(dash);
        }
    };

    // =========================================================================
    // Atom rendering
    // =========================================================================
    Renderer.prototype._drawAtom = function(atom) {
        var x = this._tx(atom.x);
        var y = this._ty(atom.y);
        var elem = atom.getElement();
        var mol = this.molecule;

        // Background highlight.
        // v1.5.2: when this atom was flagged by `highlightSearchMatch` (so
        // atom.searchMatchKind is one of 'substructure' | 'mcs' | 'exact'),
        // draw a polished two-tier "search match" halo:
        //   - outer soft pale-teal disc  (r=11, opacity ~0.32)
        //   - inner brand-teal ring stroke (r=12, opacity ~0.85)
        // plus a one-shot SMIL reveal that brightens then settles. Honours
        // `prefers-reduced-motion: reduce` (no animation, painted at the
        // resting state). Legacy callers that set bgColor directly without a
        // searchMatchKind retain the original soft flood.
        if (atom.bgColor) {
            var atomSmk = atom.searchMatchKind;
            if (atomSmk) {
                this._drawSearchMatchAtomHalo(x, y, atomSmk, atom.searchMatchTs);
            } else {
                var circle = this._circle(x, y, 12, atom.bgColor, 'none');
                circle.style.opacity = '0.4';
                this.bgLayer.appendChild(circle);
            }
        }

        // Selection highlight
        if (atom.selected) {
            this.bgLayer.appendChild(this._circle(x, y, 12, SELECT_COLOR, 'none', 0.2));
        }

        // Hover highlight — stronger teal for map-highlighted atoms
        if (atom.highlighted && atom.mapHighlighted && atom.mapNumber > 0) {
            this.bgLayer.appendChild(this._circle(x, y, 15, SELECT_COLOR, 'none', 0.3));
            this.bgLayer.appendChild(this._circle(x, y, 15, 'none', SELECT_COLOR, 0.6));
            this.bgLayer.lastChild.setAttribute('stroke-width', '2');
        } else if (atom.highlighted && !atom.searchMatchKind) {
            // v1.5.2: search-match atoms already painted a polished halo
            // above (`_drawSearchMatchAtomHalo`); skip the flat HIGHLIGHT_COLOR
            // disc to avoid stacking two halos on the same atom.
            this.bgLayer.appendChild(this._circle(x, y, 14, HIGHLIGHT_COLOR, 'none', 0.5));
        }

        // Query atom: dashed border circle
        if (atom.isQuery || atom.queryConstraints) {
            var qCircle = this._circle(x, y, 13, 'none', '#0d9488');
            qCircle.setAttribute('stroke-dasharray', '4,2');
            qCircle.setAttribute('stroke-width', '1.5');
            qCircle.style.opacity = '0.7';
            this.bgLayer.appendChild(qCircle);
        }

        // Only draw label for non-carbon atoms, or carbon with charge/isotope/no bonds
        var showLabel = atom.symbol !== 'C' || atom.charge !== 0 || atom.isotope > 0 || mol.degree(atom.id) === 0;

        if (showLabel) {
            // White background behind label
            var textWidth = this._estimateTextWidth(atom.symbol);
            this.atomLayer.appendChild(this._rect(x - textWidth / 2 - LABEL_PAD, y - FONT_SIZE / 2 - 2,
                textWidth + LABEL_PAD * 2, FONT_SIZE + 4, 'var(--color-bg, white)', 'none'));

            // Atom label
            var label = this._text(x, y + FONT_SIZE * 0.35, atom.symbol, elem.color, FONT_SIZE, 'bold');
            this.atomLayer.appendChild(label);

            // Charge
            if (atom.charge !== 0) {
                var chargeStr = atom.charge > 0 ? (atom.charge === 1 ? '+' : atom.charge + '+') :
                    (atom.charge === -1 ? '−' : Math.abs(atom.charge) + '−');
                var chargeLabel = this._text(x + textWidth / 2 + 2, y - FONT_SIZE * 0.2, chargeStr, elem.color, 9);
                this.atomLayer.appendChild(chargeLabel);
            }

            // Hydrogens — position based on bond directions
            if (this.showHydrogens) {
                var hCount = mol.calcHydrogens(atom.id);
                if (hCount > 0) {
                    var hStr = 'H';
                    var hCountStr = hCount > 1 ? '' + hCount : '';
                    var hDir = this._hydrogenDirection(atom, mol);
                    var hTextW = this._estimateTextWidth(hStr) + (hCountStr ? this._estimateTextWidth(hCountStr) * 0.6 : 0);

                    var hx, hy, hAnchor;
                    if (hDir === 'left') {
                        // Place H to the left of the atom symbol
                        hx = x - textWidth / 2 - (atom.charge !== 0 ? 10 : 2) - hTextW / 2;
                        hy = y + FONT_SIZE * 0.35;
                        hAnchor = 'end';
                    } else if (hDir === 'up') {
                        hx = x;
                        hy = y - FONT_SIZE * 0.6;
                        hAnchor = 'middle';
                    } else if (hDir === 'down') {
                        hx = x;
                        hy = y + FONT_SIZE * 1.2;
                        hAnchor = 'middle';
                    } else {
                        // Default: right
                        hx = x + textWidth / 2 + (atom.charge !== 0 ? 10 : 2) + hTextW / 2;
                        hy = y + FONT_SIZE * 0.35;
                        hAnchor = 'start';
                    }

                    // White background behind H label
                    if (hDir === 'left' || hDir === 'right') {
                        var bgX = hDir === 'right' ?
                            x + textWidth / 2 + (atom.charge !== 0 ? 10 : 2) :
                            x - textWidth / 2 - (atom.charge !== 0 ? 10 : 2) - hTextW;
                        this.atomLayer.appendChild(this._rect(bgX - 1, y - FONT_SIZE / 2 - 2,
                            hTextW + 2, FONT_SIZE + 4, 'var(--color-bg, white)', 'none'));
                    }

                    var hLabel = this._text(hx, hy, hStr + hCountStr, '#666', HYDROGEN_FONT_SIZE);
                    if (hAnchor !== 'middle') hLabel.setAttribute('text-anchor', hAnchor);
                    this.atomLayer.appendChild(hLabel);
                }
            }

            // Isotope
            if (atom.isotope > 0) {
                var isoLabel = this._text(x - textWidth / 2 - 4, y - FONT_SIZE * 0.2, '' + atom.isotope, elem.color, 9);
                isoLabel.setAttribute('text-anchor', 'end');
                this.atomLayer.appendChild(isoLabel);
            }
        }

        // Atom-atom mapping number — larger teal label below the atom
        if (this.showMapNumbers && atom.mapNumber > 0) {
            var isMapHL = atom.mapHighlighted && atom.highlighted;
            var mapFS = isMapHL ? 14 : MAP_FONT_SIZE;
            var mapY = y + FONT_SIZE + 6;
            // Background pill for readability
            var mapText = '' + atom.mapNumber;
            var mapW = mapText.length * mapFS * 0.65 + 6;
            var pillColor = isMapHL ? SELECT_COLOR : 'var(--color-bg, white)';
            var mapBg = this._rect(x - mapW / 2, mapY - mapFS * 0.7,
                mapW, mapFS + 4, pillColor, 'none');
            mapBg.style.opacity = isMapHL ? '0.2' : '0.85';
            mapBg.setAttribute('rx', '3');
            mapBg.setAttribute('ry', '3');
            this.atomLayer.appendChild(mapBg);
            var mapLabel = this._text(x, mapY + mapFS * 0.3, mapText, MAP_COLOR, mapFS, 'bold');
            mapLabel.setAttribute('font-family', FONT_MONO);
            this.atomLayer.appendChild(mapLabel);
        }

        // MCS mapping number — orange label above the atom
        if (this.showMapNumbers && atom.mcsMapNumber > 0) {
            var mcsFS = MAP_FONT_SIZE;
            var mcsY = y - FONT_SIZE - 2;
            var mcsText = '' + atom.mcsMapNumber;
            var mcsW = mcsText.length * mcsFS * 0.65 + 6;
            var mcsBg = this._rect(x - mcsW / 2, mcsY - mcsFS * 0.7,
                mcsW, mcsFS + 4, '#99f6e4', 'none');
            mcsBg.style.opacity = '0.7';
            mcsBg.setAttribute('rx', '3');
            mcsBg.setAttribute('ry', '3');
            this.atomLayer.appendChild(mcsBg);
            var mcsLabel = this._text(x, mcsY + mcsFS * 0.3, mcsText, '#0f766e', mcsFS, 'bold');
            mcsLabel.setAttribute('font-family', FONT_MONO);
            this.atomLayer.appendChild(mcsLabel);
        }

        // Hit target
        var hit = this._circle(x, y, 14, 'transparent', 'none');
        hit.setAttribute('data-atom-id', atom.id);
        hit.style.cursor = 'pointer';
        this.overlayLayer.appendChild(hit);
    };

    /**
     * Determine the best direction for placing implicit H labels.
     * Looks at bond directions to find the most open side.
     * Returns 'right', 'left', 'up', or 'down'.
     */
    Renderer.prototype._hydrogenDirection = function(atom, mol) {
        var neighbors = mol.getNeighbors(atom.id);
        if (neighbors.length === 0) return 'right';

        // Compute average direction of bonds
        var sumX = 0, sumY = 0;
        for (var i = 0; i < neighbors.length; i++) {
            var n = mol.getAtom(neighbors[i]);
            if (!n) continue;
            var dx = n.x - atom.x;
            var dy = n.y - atom.y;
            var len = Math.sqrt(dx * dx + dy * dy);
            if (len > 0) { sumX += dx / len; sumY += dy / len; }
        }

        // Place H opposite to the average bond direction
        var angle = Math.atan2(-sumY, -sumX); // opposite direction
        // Snap to cardinal direction
        if (angle > -Math.PI / 4 && angle <= Math.PI / 4) return 'right';
        if (angle > Math.PI / 4 && angle <= 3 * Math.PI / 4) return 'down';
        if (angle > -3 * Math.PI / 4 && angle <= -Math.PI / 4) return 'up';
        return 'left';
    };

    Renderer.prototype._labelRadius = function(atom) {
        if (atom.symbol === 'C' && atom.charge === 0 && atom.isotope === 0 && this.molecule.degree(atom.id) > 0) return 0;
        return this._estimateTextWidth(atom.symbol) / 2 + LABEL_PAD;
    };

    // =========================================================================
    // Reaction arrow (forward >> and retrosynthetic =>)
    // =========================================================================
    Renderer.prototype._drawReactionArrow = function(arrow) {
        var x1 = this._tx(arrow.x1), y1 = this._ty(arrow.y1);
        var x2 = this._tx(arrow.x2), y2 = this._ty(arrow.y2);
        var isRetro = (arrow.type === 'retro');
        var arrowColor = '#333';

        var dx = x2 - x1, dy = y2 - y1;
        var len = Math.sqrt(dx * dx + dy * dy);
        if (len < 1) return;
        var nx = dx / len, ny = dy / len;
        var px = -ny, py = nx;
        var hs = ARROW_HEAD_SIZE;

        if (isRetro) {
            // ---- Retrosynthetic arrow: open (hollow) arrowhead, no tail ----
            // Shaft line (stops short of head so it doesn't poke through)
            var shaftEnd = len - hs;
            var sx2 = x1 + nx * shaftEnd, sy2 = y1 + ny * shaftEnd;
            var shaft = this._line(x1, y1, sx2, sy2, arrowColor, ARROW_STROKE_WIDTH);
            this.overlayLayer.appendChild(shaft);

            // Open arrowhead (two lines forming a V)
            var hw = hs * 0.5; // half-width
            var tipX = x2, tipY = y2;
            var wing1X = tipX - nx * hs + px * hw;
            var wing1Y = tipY - ny * hs + py * hw;
            var wing2X = tipX - nx * hs - px * hw;
            var wing2Y = tipY - ny * hs - py * hw;

            var headLine1 = this._line(wing1X, wing1Y, tipX, tipY, arrowColor, ARROW_STROKE_WIDTH);
            var headLine2 = this._line(wing2X, wing2Y, tipX, tipY, arrowColor, ARROW_STROKE_WIDTH);
            this.overlayLayer.appendChild(headLine1);
            this.overlayLayer.appendChild(headLine2);
        } else {
            // ---- Forward arrow: solid filled triangular arrowhead ----
            // Shaft line (stops at base of arrowhead)
            var shaftEnd = len - hs * 0.85;
            var sx2 = x1 + nx * shaftEnd, sy2 = y1 + ny * shaftEnd;
            var shaft = this._line(x1, y1, sx2, sy2, arrowColor, ARROW_STROKE_WIDTH);
            this.overlayLayer.appendChild(shaft);

            // Filled triangular arrowhead
            var hw = hs * 0.45; // half-width
            var tipX = x2, tipY = y2;
            var baseX = tipX - nx * hs, baseY = tipY - ny * hs;
            var points = tipX + ',' + tipY + ' ' +
                (baseX + px * hw) + ',' + (baseY + py * hw) + ' ' +
                (baseX - px * hw) + ',' + (baseY - py * hw);
            var head = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
            head.setAttribute('points', points);
            head.setAttribute('fill', arrowColor);
            head.setAttribute('stroke', arrowColor);
            head.setAttribute('stroke-width', '0.5');
            head.setAttribute('stroke-linejoin', 'round');
            this.overlayLayer.appendChild(head);
        }

        // ---- Conditions text above arrow (if provided) ----
        if (arrow.conditions) {
            var midX = (x1 + x2) / 2;
            var midY = (y1 + y2) / 2 - 12;
            var condLabel = this._text(midX, midY, arrow.conditions, '#666', 10);
            condLabel.setAttribute('font-style', 'italic');
            this.overlayLayer.appendChild(condLabel);
        }
    };

    // =========================================================================
    // "+" sign between multi-component reactants / products
    // =========================================================================
    Renderer.prototype._drawPlusSign = function(pos) {
        var x = this._tx(pos.x);
        var y = this._ty(pos.y);
        var plusLabel = this._text(x, y + PLUS_FONT_SIZE * 0.35, '+', PLUS_COLOR, PLUS_FONT_SIZE, 'bold');
        this.overlayLayer.appendChild(plusLabel);
    };

    // =========================================================================
    // CIP stereochemistry label rendering
    // =========================================================================

    /**
     * Draw CIP labels: (R)/(S) near stereocentres, (E)/(Z) near double bonds.
     * Labels are rendered in small italic grey text on the overlay layer.
     */
    Renderer.prototype._drawCipLabels = function(mol, showRS, showEZ) {
        var i, atom, bond, x, y, label;
        if (showRS === undefined) showRS = true;
        if (showEZ === undefined) showEZ = true;

        // R/S labels at stereocentres
        for (i = 0; i < mol.atoms.length; i++) {
            atom = mol.atoms[i];
            if (!atom.cipLabel || !showRS) continue;

            x = this._tx(atom.x);
            y = this._ty(atom.y);

            // Position label slightly below and to the right of the atom
            var offsetX = 10;
            var offsetY = -10;

            label = this._text(x + offsetX, y + offsetY, '(' + atom.cipLabel + ')',
                CIP_COLOR, CIP_FONT_SIZE);
            label.setAttribute('font-style', 'italic');
            label.setAttribute('text-anchor', 'start');
            label.setAttribute('font-family', FONT_FAMILY);
            this.overlayLayer.appendChild(label);
        }

        // E/Z labels at double bonds
        for (i = 0; i < mol.bonds.length; i++) {
            bond = mol.bonds[i];
            if (!bond.cipLabel || !showEZ) continue;

            var a1 = mol.getAtom(bond.atom1);
            var a2 = mol.getAtom(bond.atom2);
            if (!a1 || !a2) continue;

            // Position label at the midpoint of the bond, offset perpendicular
            var mx = this._tx((a1.x + a2.x) / 2);
            var my = this._ty((a1.y + a2.y) / 2);

            // Offset perpendicular to the bond
            var dx = a2.x - a1.x;
            var dy = a2.y - a1.y;
            var len = Math.sqrt(dx * dx + dy * dy);
            if (len > 0) {
                var perpX = -dy / len * 10 * this.scale;
                var perpY = dx / len * 10 * this.scale;
                mx += perpX;
                my += perpY;
            }

            label = this._text(mx, my, '(' + bond.cipLabel + ')',
                CIP_COLOR, CIP_FONT_SIZE);
            label.setAttribute('font-style', 'italic');
            label.setAttribute('font-family', FONT_FAMILY);
            this.overlayLayer.appendChild(label);
        }
    };

    // =========================================================================
    // Atom-atom mapping connection lines
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /**
     * Draw dashed teal bezier curves connecting atoms that share a highlighted
     * mapNumber across the reaction (e.g., reactant atom :1 <-> product atom :1).
     * Only renders when at least two atoms with the same mapNumber are highlighted.
     */
    Renderer.prototype._drawMapConnections = function(mol) {
        if (!mol || !mol.reactionArrow) return;

        // Group highlighted map-atoms by their mapNumber
        var mapGroups = {};
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            if (a.mapHighlighted && a.highlighted && a.mapNumber > 0) {
                if (!mapGroups[a.mapNumber]) mapGroups[a.mapNumber] = [];
                mapGroups[a.mapNumber].push(a);
            }
        }

        // Draw curved dashed connections for each map group with 2+ atoms
        for (var mn in mapGroups) {
            if (!mapGroups.hasOwnProperty(mn)) continue;
            var atoms = mapGroups[mn];
            if (atoms.length < 2) continue;

            // Connect each pair (usually just 2 — one reactant, one product)
            for (var j = 0; j < atoms.length - 1; j++) {
                for (var k = j + 1; k < atoms.length; k++) {
                    var a1 = atoms[j], a2 = atoms[k];
                    var x1 = this._tx(a1.x), y1 = this._ty(a1.y);
                    var x2 = this._tx(a2.x), y2 = this._ty(a2.y);

                    // Bezier control point: arch upward over the reaction arrow
                    var midX = (x1 + x2) / 2;
                    var midY = (y1 + y2) / 2;
                    var dist = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                    var archHeight = Math.max(30, dist * 0.3);
                    // Arch above the midpoint (negative Y = upward in SVG)
                    var cpY = midY - archHeight;

                    var path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                    path.setAttribute('d', 'M' + x1 + ',' + y1 +
                        ' Q' + midX + ',' + cpY + ' ' + x2 + ',' + y2);
                    path.setAttribute('fill', 'none');
                    path.setAttribute('stroke', SELECT_COLOR);
                    path.setAttribute('stroke-width', '1.5');
                    path.setAttribute('stroke-dasharray', '5,3');
                    path.setAttribute('stroke-linecap', 'round');
                    path.style.opacity = '0.6';
                    this.overlayLayer.appendChild(path);
                }
            }
        }
    };

    // =========================================================================
    // Molecule name rendering
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    var MOL_NAME_FONT_SIZE = 14;
    var MOL_NAME_PAD = 20;

    /**
     * Draw the molecule name centred below the structure bounding box.
     * Renders a semi-transparent white pill background behind the text
     * for readability, then bold text in the configured colour.
     */
    Renderer.prototype._drawMolName = function(mol) {
        var bounds = mol.getBounds();
        // Centre x of the molecule in screen coords
        var cx = this._tx(bounds.x + bounds.w / 2);
        // Bottom of molecule + padding in screen coords
        var bottomY = this._ty(bounds.y + bounds.h) + MOL_NAME_PAD * this.scale;
        var nameText = mol.name;

        // Estimate text width for background pill
        var textW = nameText.length * MOL_NAME_FONT_SIZE * 0.55;
        var pillW = textW + 16;
        var pillH = MOL_NAME_FONT_SIZE + 8;

        // Semi-transparent white pill background
        var bg = this._rect(cx - pillW / 2, bottomY - pillH / 2,
            pillW, pillH, 'rgba(255,255,255,0.85)', 'none');
        bg.setAttribute('rx', '4');
        bg.setAttribute('ry', '4');
        this.overlayLayer.appendChild(bg);

        // Name text
        var label = this._text(cx, bottomY + MOL_NAME_FONT_SIZE * 0.35,
            nameText, 'var(--color-text, #333)', MOL_NAME_FONT_SIZE, 'bold');
        label.setAttribute('class', 'bime-mol-name');
        this.overlayLayer.appendChild(label);
    };

    /**
     * Return the screen-space bounding box of the rendered molecule name.
     * Used by MolEditor for hit-testing double-click-to-edit on the name.
     * Returns { x, y, w, h } or null if name is not shown.
     */
    Renderer.prototype.getMolNameBounds = function() {
        if (!this.showMolName || !this.molecule || !this.molecule.name) return null;
        var bounds = this.molecule.getBounds();
        var cx = this._tx(bounds.x + bounds.w / 2);
        var bottomY = this._ty(bounds.y + bounds.h) + MOL_NAME_PAD * this.scale;
        var textW = this.molecule.name.length * MOL_NAME_FONT_SIZE * 0.55;
        var pillW = textW + 16;
        var pillH = MOL_NAME_FONT_SIZE + 8;
        return { x: cx - pillW / 2, y: bottomY - pillH / 2, w: pillW, h: pillH };
    };

    // =========================================================================
    // Preview elements (for drawing in progress)
    // =========================================================================
    Renderer.prototype.drawPreviewBond = function(x1, y1, x2, y2) {
        this._clear(this.overlayLayer);
        this.overlayLayer.appendChild(this._line(
            this._tx(x1), this._ty(y1), this._tx(x2), this._ty(y2), SELECT_COLOR, 1.5));
    };

    Renderer.prototype.drawPreviewAtom = function(x, y, symbol) {
        var px = this._tx(x), py = this._ty(y);
        this.overlayLayer.appendChild(this._circle(px, py, 12, SELECT_COLOR, 'none', 0.2));
        if (symbol !== 'C') {
            this.overlayLayer.appendChild(this._text(px, py + FONT_SIZE * 0.35, symbol, SELECT_COLOR, FONT_SIZE, 'bold'));
        }
    };

    Renderer.prototype.clearPreview = function() {
        this._clear(this.overlayLayer);
    };

    // =========================================================================
    // Coordinate transforms
    // =========================================================================
    Renderer.prototype._tx = function(x) { return (x + this.offsetX) * this.scale; };
    Renderer.prototype._ty = function(y) { return (y + this.offsetY) * this.scale; };
    Renderer.prototype.screenToMol = function(sx, sy) {
        var rect = this.svg.getBoundingClientRect();
        var vx = (sx - rect.left) / rect.width * this.width;
        var vy = (sy - rect.top) / rect.height * this.height;
        return { x: vx / this.scale - this.offsetX, y: vy / this.scale - this.offsetY };
    };

    Renderer.prototype.centerMolecule = function() {
        if (!this.molecule || this.molecule.isEmpty()) return;
        var bounds = this.molecule.getBounds();

        // Auto-scale reactions to fit viewport width with padding
        if (this.molecule.reactionArrow) {
            var pad = BOND_LENGTH * 2;
            var scaleX = (this.width  - pad) / (bounds.w > 0 ? bounds.w : 1);
            var scaleY = (this.height - pad) / (bounds.h > 0 ? bounds.h : 1);
            var fitScale = Math.min(scaleX, scaleY, 1.6); // cap at 1.6x to avoid giant rendering
            fitScale = Math.max(fitScale, 0.3);            // floor so it stays readable
            this.scale = fitScale;
        }

        this.offsetX = (this.width / this.scale - bounds.w) / 2 - bounds.x;
        this.offsetY = (this.height / this.scale - bounds.h) / 2 - bounds.y;
    };

    // =========================================================================
    // SVG helpers
    // =========================================================================
    Renderer.prototype._clear = function(g) {
        while (g.firstChild) g.removeChild(g.firstChild);
    };

    Renderer.prototype._line = function(x1, y1, x2, y2, stroke, width, cap) {
        var l = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        l.setAttribute('x1', x1); l.setAttribute('y1', y1);
        l.setAttribute('x2', x2); l.setAttribute('y2', y2);
        l.setAttribute('stroke', stroke || '#333');
        l.setAttribute('stroke-width', width || 1.5);
        l.setAttribute('stroke-linecap', cap || 'round');
        return l;
    };

    Renderer.prototype._circle = function(cx, cy, r, fill, stroke, opacity) {
        var c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        c.setAttribute('cx', cx); c.setAttribute('cy', cy); c.setAttribute('r', r);
        c.setAttribute('fill', fill || 'none');
        c.setAttribute('stroke', stroke || 'none');
        if (opacity !== undefined) c.style.opacity = opacity;
        return c;
    };

    Renderer.prototype._rect = function(x, y, w, h, fill, stroke) {
        var r = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        r.setAttribute('x', x); r.setAttribute('y', y);
        r.setAttribute('width', w); r.setAttribute('height', h);
        r.setAttribute('fill', fill || 'none');
        r.setAttribute('stroke', stroke || 'none');
        return r;
    };

    Renderer.prototype._text = function(x, y, text, fill, fontSize, weight) {
        var t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        t.setAttribute('x', x); t.setAttribute('y', y);
        t.setAttribute('fill', fill || '#333');
        t.setAttribute('font-size', (fontSize || FONT_SIZE) + 'px');
        t.setAttribute('font-family', FONT_FAMILY);
        t.setAttribute('text-anchor', 'middle');
        if (weight) t.setAttribute('font-weight', weight);
        t.textContent = text;
        return t;
    };

    Renderer.prototype._estimateTextWidth = function(text) {
        return text.length * FONT_SIZE * 0.6;
    };

    Renderer.prototype.getSVGString = function() {
        return this.svg.outerHTML;
    };

    global.Renderer = Renderer;

})(typeof window !== 'undefined' ? window : this);
