/**
 * tests/test_v2_3_7_accessibility.js — accessibility hardening (v2.3.7).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.3.7 closes two WCAG failures a UX review flagged. No chemistry or algorithm
 * change. This suite pins both fixes:
 *
 *   P1-5  Use of Colour (WCAG 1.4.1): three encodings that were carried by
 *         colour ALONE gain a redundant, greyscale-readable cue.
 *           (a) PATHWAY NODE TYPE — renderPathwayNode now stamps a small
 *               uppercase letter badge (M metabolite / R residue / Co cofactor /
 *               Rxn reaction) on top of the (kept) coloured fill/stroke, so a
 *               greyscale or colour-blind viewer can tell a residue from a
 *               metabolite from a cofactor.
 *           (b) REACTION-ARROW TYPE — the four arrow types are shape-distinct:
 *               forward (one head at end), reverse (one head at start),
 *               reversible (two offset half-harpoons, no full heads) and
 *               resonance (full heads at BOTH ends, on a DASHED shaft so it is
 *               not confused with a solid forward/reverse arrow). The arrow
 *               colour never encoded type (it only flips on selection), so the
 *               fix is purely shape/pattern.
 *           (c) ATOM-TRACE STATE — the four moiety verdicts now carry DISTINCT
 *               border PATTERNS (intact solid / fragmented dashed / partial
 *               dotted / absent double) plus a status glyph+word tag in the cell
 *               title, instead of intact+fragmented both-solid and
 *               partial+absent both-dashed differing only by hue.
 *
 *   P1-6  the editor draw surface is REACHABLE + NAMED + DESCRIBED (WCAG 2.1.1
 *         keyboard, 4.1.2 name/role/value). MolEditor._drawArea was an unlabeled,
 *         unfocusable <div>; it now carries role=img, tabindex=0 and an aria-label
 *         that is a LIVE description of what is drawn ("empty canvas", else the
 *         molecule name + atom/bond tally + "reaction scheme" when an arrow is
 *         present), refreshed from the existing changed() signal. A focus-visible
 *         ring is injected so keyboard focus is visible. Reachable + described
 *         only — full keyboard atom PLACEMENT is a later (Phase B) item.
 *
 * Source-shape assertions (the workbench glue + MolEditor DOM construction need
 * the whole page / a real renderer to execute, exactly as test_v2_3_0 /
 * test_v2_3_4 / test_v2_3_6 note) PLUS executable sections: the workbench cue
 * mappings are modelled the way test_v2_3_6 section E models control flow, and
 * the MolEditor live-description methods are exercised FOR REAL against a real
 * Molecule via the prototype (MolEditor.js loads cleanly under the shim).
 *
 * Plain Node, no real DOM beyond the shim.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

// MolEditor + its Renderer dependency load under the shim (the constructor only
// touches the DOM when instantiated; the prototype methods we exercise below
// take a hand-built `this`). Pull the real classes off the global.
require('../editor/Renderer.js');
require('../editor/MolEditor.js');
var Molecule = globalThis.Molecule;
var MolEditor = globalThis.MolEditor;

var runner = shim.makeRunner('Accessibility: non-colour cues + reachable canvas (v2.3.7)');
var test = runner.test;
console.log('Accessibility: non-colour cues + reachable canvas (v2.3.7)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var CSS = readRepoFile('css/workbench.css');
var ME = readRepoFile('editor/MolEditor.js');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing the source-shape-pinned suites use.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// ===========================================================================
// A. P1-5a — pathway node TYPE letter badge (renderPathwayNode).
// ===========================================================================
test('A1 the badge helpers exist as top-level functions', function () {
    assert.ok(WB.indexOf('function pathwayNodeTypeBadgeText(') !== -1,
        'pathwayNodeTypeBadgeText defined');
    assert.ok(WB.indexOf('function appendPathwayNodeTypeBadge(') !== -1,
        'appendPathwayNodeTypeBadge defined');
});

test('A2 renderPathwayNode is NOT relocated and delegates to the badge helper', function () {
    // renderPathwayNode must still exist as a top-level function (source-shape
    // contract) and must call the new badge helper before appending to the
    // viewport, so EVERY node gets the non-colour cue.
    var fn = wbFn('renderPathwayNode');
    assert.ok(fn, 'renderPathwayNode located as a top-level function');
    assert.ok(fn.indexOf('appendPathwayNodeTypeBadge(g, node)') !== -1,
        'renderPathwayNode stamps the type badge');
    var badgeIdx = fn.indexOf('appendPathwayNodeTypeBadge(g, node)');
    var appendIdx = fn.lastIndexOf('viewport.appendChild(g)');
    assert.ok(badgeIdx !== -1 && appendIdx !== -1 && badgeIdx < appendIdx,
        'the badge is added before the group is appended to the viewport');
});

test('A3 each node type maps to a DISTINCT non-colour letter badge', function () {
    // Model the badge-text decision exactly as pathwayNodeTypeBadgeText does, so
    // a future edit that collapses two types to the same tag fails here.
    var TAGS = { metabolite: 'M', residue: 'R', cofactor: 'Co', reaction: 'Rxn' };
    function badgeText(node) {
        if (!node) { return ''; }
        var kind = node.kind;
        if (!kind && node.structure && node.structure.type === 'reaction') { kind = 'reaction'; }
        if (kind === 'reaction' || (node.structure && node.structure.type === 'reaction')) { return 'Rxn'; }
        if (Object.prototype.hasOwnProperty.call(TAGS, kind)) { return TAGS[kind]; }
        return 'M';
    }
    var metabolite = badgeText({ kind: 'metabolite' });
    var dflt = badgeText({});                       // untyped/white node
    var residue = badgeText({ kind: 'residue' });
    var cofactor = badgeText({ kind: 'cofactor' });
    var reaction = badgeText({ kind: 'reaction' });
    var rxnByStructure = badgeText({ structure: { type: 'reaction' } });

    assert.strictEqual(metabolite, 'M', 'metabolite => M');
    assert.strictEqual(dflt, 'M', 'untyped white node => M (still cued)');
    assert.strictEqual(residue, 'R', 'residue => R');
    assert.strictEqual(cofactor, 'Co', 'cofactor => Co');
    assert.strictEqual(reaction, 'Rxn', 'reaction => Rxn');
    assert.strictEqual(rxnByStructure, 'Rxn', 'reaction-by-structure => Rxn');

    // The four meaningful types must be pairwise distinct (a greyscale viewer
    // tells residue from metabolite from cofactor from reaction).
    var set = {};
    [residue, cofactor, reaction, metabolite].forEach(function (t) { set[t] = true; });
    assert.strictEqual(Object.keys(set).length, 4,
        'residue / cofactor / reaction / metabolite badges are all distinct');
    // The helper source agrees with the model.
    assert.ok(WB.indexOf("residue: 'R'") !== -1 && WB.indexOf("cofactor: 'Co'") !== -1 &&
        WB.indexOf("reaction: 'Rxn'") !== -1 && WB.indexOf("metabolite: 'M'") !== -1,
        'PATHWAY_NODE_TYPE_BADGE map carries the four distinct tags');
});

test('A4 CSS styles the badge as a greyscale-readable chip (not colour-only)', function () {
    assert.ok(CSS.indexOf('.wb-pathway-node-badge') !== -1, 'badge class is styled');
    var idx = CSS.indexOf('.wb-pathway-node .wb-pathway-node-badge rect');
    assert.ok(idx !== -1, 'badge rect rule present');
    var rule = CSS.slice(idx, idx + 160);
    // Slate stroke + pale fill => visible chip outline in greyscale; the drop
    // shadow inherited from .wb-pathway-node rect is cancelled.
    assert.ok(rule.indexOf('stroke:') !== -1, 'badge chip has a border (greyscale-visible)');
    assert.ok(rule.indexOf('filter: none') !== -1, 'badge rect cancels the node drop-shadow');
});

// ===========================================================================
// B. P1-5b — reaction-arrow types are shape-distinguishable, not colour-coded.
// ===========================================================================
test('B1 drawPathwayEdgeArrow is NOT relocated and the resonance shaft is dashed', function () {
    var fn = wbFn('drawPathwayEdgeArrow');
    assert.ok(fn, 'drawPathwayEdgeArrow located as a top-level function');
    // The shaft helper takes a `dash` flag and applies a dasharray.
    assert.ok(/function line\([^)]*dash\)/.test(fn), 'the shaft line helper accepts a dash flag');
    assert.ok(fn.indexOf("attrs['stroke-dasharray']") !== -1,
        'a dashed shaft variant exists');
    assert.ok(fn.indexOf("line(pts.x1, pts.y1, pts.x2, pts.y2, type === 'resonance')") !== -1,
        'resonance draws a DASHED shaft; forward/reverse/reversible stay solid');
});

test('B2 reversible draws harpoons (no full heads) and resonance draws two full heads', function () {
    var fn = wbFn('drawPathwayEdgeArrow');
    // reversible: early-return branch that appends harpoons.
    assert.ok(/if \(type === 'reversible'\)/.test(fn), 'reversible has its own branch');
    assert.ok(fn.indexOf('harpoon(') !== -1, 'reversible uses single-barb harpoons');
    // forward/resonance head at end; reverse/resonance head at start.
    assert.ok(fn.indexOf("type === 'forward' || type === 'resonance'") !== -1,
        'forward + resonance get a head at the end');
    assert.ok(fn.indexOf("type === 'reverse' || type === 'resonance'") !== -1,
        'reverse + resonance get a head at the start');
});

test('B3 the four arrow types have pairwise-distinct shape fingerprints', function () {
    // Model the shape each type draws (heads at start/end, harpoons, dashed
    // shaft) the way drawPathwayEdgeArrow does, and assert all four differ.
    function fingerprint(type) {
        if (type === 'reversible') {
            return { headStart: false, headEnd: false, harpoons: true, dashed: false };
        }
        return {
            headStart: (type === 'reverse' || type === 'resonance'),
            headEnd: (type === 'forward' || type === 'resonance'),
            harpoons: false,
            dashed: (type === 'resonance')
        };
    }
    var types = ['forward', 'reverse', 'reversible', 'resonance'];
    var seen = {};
    for (var i = 0; i < types.length; i++) {
        var key = JSON.stringify(fingerprint(types[i]));
        assert.ok(!seen[key], types[i] + ' has a UNIQUE shape fingerprint (not shared)');
        seen[key] = true;
    }
    // Spot the previously-ambiguous pair: forward vs resonance now differ by BOTH
    // a second head AND the dashed shaft.
    var fwd = fingerprint('forward');
    var res = fingerprint('resonance');
    assert.notStrictEqual(fwd.headStart, res.headStart, 'forward vs resonance differ at the start head');
    assert.notStrictEqual(fwd.dashed, res.dashed, 'forward vs resonance differ by shaft pattern');
});

// ===========================================================================
// C. P1-5c — atom-trace states get distinct border patterns + a status tag.
// ===========================================================================
test('C1 CSS gives the four moiety verdicts DISTINCT border patterns', function () {
    function borderStyle(stateClass) {
        var idx = CSS.indexOf('.wb-atom-trace-cell.' + stateClass + ' {');
        assert.ok(idx !== -1, stateClass + ' rule present');
        var block = CSS.slice(idx, CSS.indexOf('}', idx));
        var m = block.match(/border-style:\s*([a-z]+)/);
        assert.ok(m, stateClass + ' declares a border-style');
        return m[1];
    }
    var intact = borderStyle('is-moiety-intact');
    var fragmented = borderStyle('is-moiety-fragmented');
    var partial = borderStyle('is-moiety-partial-loss');
    var absent = borderStyle('is-moiety-absent');
    assert.strictEqual(intact, 'solid', 'intact => solid');
    assert.strictEqual(fragmented, 'dashed', 'fragmented => dashed');
    assert.strictEqual(partial, 'dotted', 'partial-loss => dotted');
    assert.strictEqual(absent, 'double', 'absent => double');
    var styles = {};
    [intact, fragmented, partial, absent].forEach(function (s) { styles[s] = true; });
    assert.strictEqual(Object.keys(styles).length, 4,
        'all four trace states have a distinct border pattern (readable without colour)');
});

test('C2 the moiety status-tag helper exists and maps each verdict to a glyph+word', function () {
    assert.ok(WB.indexOf('function atomTraceMoietyTag(') !== -1, 'atomTraceMoietyTag defined');
    // Model the mapping the source declares.
    var TAGS = {
        'intact': '✓ intact',
        'fragmented': '⚠ fragmented',
        'partial-loss': '◐ partial',
        'absent': '✕ absent'
    };
    function tag(status) {
        return Object.prototype.hasOwnProperty.call(TAGS, status) ? TAGS[status] : '';
    }
    assert.strictEqual(tag('intact'), '✓ intact');
    assert.strictEqual(tag('fragmented'), '⚠ fragmented');
    assert.strictEqual(tag('partial-loss'), '◐ partial');
    assert.strictEqual(tag('absent'), '✕ absent');
    assert.strictEqual(tag('unreachable'), '', 'unreachable carries no tag (no border either)');
    // Each tag includes a word (so it is screen-reader meaningful, not glyph-only)
    // and all four are distinct.
    var vals = ['intact', 'fragmented', 'partial-loss', 'absent'].map(tag);
    var distinct = {};
    vals.forEach(function (v) {
        distinct[v] = true;
        assert.ok(/[a-z]/.test(v), 'tag "' + v + '" includes a readable word');
    });
    assert.strictEqual(Object.keys(distinct).length, 4, 'all four tags are distinct');
    // Source carries the same vocabulary.
    assert.ok(WB.indexOf("'intact': '✓ intact'") !== -1 &&
        WB.indexOf("'absent': '✕ absent'") !== -1,
        'ATOM_TRACE_MOIETY_TAGS map present in source');
});

test('C3 showAtomTraceStrip stamps the status tag into the cell title (not colour-only)', function () {
    var fn = WB.substring(WB.indexOf('function showAtomTraceStrip'));
    fn = fn.slice(0, 16000);
    assert.ok(fn.indexOf('atomTraceMoietyTag(') !== -1,
        'showAtomTraceStrip computes a status tag');
    assert.ok(fn.indexOf('wb-atom-trace-moiety-tag') !== -1,
        'the tag is rendered with its own class in the cell title');
});

// ===========================================================================
// D. P1-6 — the draw surface is reachable + named + described.
// ===========================================================================
test('D1 _drawArea gets role=img, tabindex=0 and an aria-label', function () {
    // Slice the draw-area construction block.
    var idx = ME.indexOf('this._drawArea = document.createElement');
    assert.ok(idx !== -1, '_drawArea is constructed');
    // Widen past the rationale comment to the setAttribute calls.
    var block = ME.slice(idx, idx + 2400);
    assert.ok(block.indexOf("setAttribute('role', 'img')") !== -1,
        'draw area gets role=img (documented: not application until Phase B keyboard support)');
    assert.ok(block.indexOf("setAttribute('tabindex', '0')") !== -1,
        'draw area is in the focus order (tabindex=0)');
    assert.ok(block.indexOf("setAttribute('aria-label', 'Molecule drawing canvas") !== -1,
        'draw area has a Molecule-drawing-canvas aria-label');
    assert.ok(block.indexOf("'bime-draw-area'") !== -1, 'draw area carries the focus-style hook class');
});

test('D2 a visible focus-visible ring is injected for the draw area', function () {
    assert.ok(ME.indexOf('.bime-draw-area:focus-visible') !== -1,
        'a focus-visible rule for the draw area is injected');
    var idx = ME.indexOf('.bime-draw-area:focus-visible');
    var rule = ME.slice(idx, idx + 80);
    assert.ok(rule.indexOf('box-shadow') !== -1,
        'the focus-visible rule paints a visible ring (mirrors .bime-btn:focus-visible)');
});

test('D3 changed() refreshes the live description from the existing signal', function () {
    var idx = ME.indexOf('MolEditor.prototype.changed = function');
    var fn = ME.slice(idx, ME.indexOf('};', idx));
    assert.ok(fn.indexOf('this._updateCanvasDescription()') !== -1,
        'changed() updates the canvas description (rides the existing change signal)');
    // It must also still fire the original callback (no regression).
    assert.ok(fn.indexOf("_fireCallback('AfterStructureModified')") !== -1,
        'changed() still fires AfterStructureModified');
});

test('D4 canvasDescription() reports "empty canvas" for a fresh molecule (REAL method)', function () {
    var mol = new Molecule();
    var self = Object.create(MolEditor.prototype);
    self.molecule = mol;
    assert.strictEqual(self.canvasDescription(), 'Molecule drawing canvas, empty canvas',
        'no atoms => empty-canvas description');
});

test('D5 canvasDescription() reports the atom/bond tally + name + reaction (REAL method)', function () {
    var mol = new Molecule();
    var a1 = mol.addAtom('C', 0, 0);
    var a2 = mol.addAtom('C', 1, 0);
    var a3 = mol.addAtom('O', 2, 0);
    mol.addBond(a1.id, a2.id, 1);
    mol.addBond(a2.id, a3.id, 1);
    var self = Object.create(MolEditor.prototype);
    self.molecule = mol;
    assert.strictEqual(self.canvasDescription(),
        'Molecule drawing canvas: 3 atoms, 2 bonds',
        'tally is pluralised correctly');

    mol.name = 'Ethanol';
    assert.strictEqual(self.canvasDescription(),
        'Molecule drawing canvas: Ethanol, 3 atoms, 2 bonds',
        'molecule name is included when present');

    mol.reactionArrow = { x1: 0, y1: 0, x2: 1, y2: 0 };
    assert.strictEqual(self.canvasDescription(),
        'Molecule drawing canvas: Ethanol, 3 atoms, 2 bonds, reaction scheme',
        'a reaction arrow is announced as a reaction scheme');
});

test('D6 canvasDescription() singularises a single atom / single bond (REAL method)', function () {
    var mol = new Molecule();
    var a1 = mol.addAtom('C', 0, 0);
    var a2 = mol.addAtom('N', 1, 0);
    mol.addBond(a1.id, a2.id, 3);
    var self = Object.create(MolEditor.prototype);
    self.molecule = mol;
    assert.strictEqual(self.canvasDescription(),
        'Molecule drawing canvas: 2 atoms, 1 bond',
        '1 bond is singular');

    var mol1 = new Molecule();
    mol1.addAtom('C', 0, 0);
    var self1 = Object.create(MolEditor.prototype);
    self1.molecule = mol1;
    assert.strictEqual(self1.canvasDescription(),
        'Molecule drawing canvas: 1 atom, 0 bonds',
        '1 atom is singular, 0 bonds is plural');
});

test('D7 _updateCanvasDescription() writes the description onto the draw-area aria-label (REAL method)', function () {
    var mol = new Molecule();
    var a1 = mol.addAtom('C', 0, 0);
    var a2 = mol.addAtom('C', 1, 0);
    mol.addBond(a1.id, a2.id, 1);
    var rec = { _attrs: {}, setAttribute: function (k, v) { this._attrs[k] = v; } };
    var self = Object.create(MolEditor.prototype);
    self.molecule = mol;
    self._drawArea = rec;
    self._updateCanvasDescription();
    assert.strictEqual(rec._attrs['aria-label'],
        'Molecule drawing canvas: 2 atoms, 1 bond',
        'aria-label is updated to the live description');

    // Updates again as the structure changes (drives the changed() contract).
    mol.addAtom('O', 2, 0);
    self._updateCanvasDescription();
    assert.strictEqual(rec._attrs['aria-label'],
        'Molecule drawing canvas: 3 atoms, 1 bond',
        'aria-label tracks subsequent edits');
});

test('D8 _updateCanvasDescription() is a no-op when the draw area is absent (no throw)', function () {
    var self = Object.create(MolEditor.prototype);
    self.molecule = new Molecule();
    self._drawArea = null;
    assert.doesNotThrow(function () { self._updateCanvasDescription(); },
        'guards against a missing draw area (teardown / headless)');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
