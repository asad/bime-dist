/**
 * tests/test_v2_4_8_editor_affordances.js — the dedicated Editor view exposes the
 * chemistry display toggles + the pan tool (v2.4.8).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Follow-on to the v2.4.6 reaction-group fix: the same continuum redesign that
 * hid the Reaction group also force-hid five toolbar-row-2 affordances in the
 * editor — Pan Mol (panmol), R/S (ciprs), E/Z (cipez), H (toggleh), Ar
 * (togglearo) — to declutter the in-place FOCUS LENS. But in the dedicated v2.4.0
 * EDITOR view that hide stranded real features:
 *   - H and Ar had NO other route at all (no command, no menu);
 *   - R/S and E/Z were reachable only through the command palette, never a button;
 *   - Pan Mol is the editor's ONLY pan (empty-drag box-selects; the wheel zooms).
 * v2.4.8 re-shows all five in the dedicated Editor view ONLY, via an override
 * scoped to `.wb-main.wb-view-editor` — strictly more specific than the hide and
 * unable to match the lens (which carries `.wb-lens-open`, never `.wb-view-editor`),
 * so the focus-lens/pathway decluttering is untouched. The `actions` + `export`
 * rail groups deliberately stay hidden (covered by Editor Output + Export Presets).
 *
 * Source-shape contract on css/workbench.css + editor/MolEditor.js (the same level
 * as test_v2_4_6 / test_v2_4_0 §C), with CSS comments stripped so prose can't
 * spoof a rule-presence regex.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Editor view exposes display toggles + pan (v2.4.8)');
var test = runner.test;
console.log('Editor view exposes display toggles + pan (v2.4.8)');

function readRepoFile(rel) { return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8'); }

// Strip CSS block comments first (this file's own fix comment names selectors).
var CSS = readRepoFile('css/workbench.css').replace(/\/\*[\s\S]*?\*\//g, '');
var MOL = readRepoFile('editor/MolEditor.js');

var TOOLS = ['panmol', 'ciprs', 'cipez', 'toggleh', 'togglearo'];

function selectorWeight(sel) {
    return (sel.match(/\.[A-Za-z0-9_-]+/g) || []).length
        + (sel.match(/\[[^\]]+\]/g) || []).length
        + (sel.match(/#[A-Za-z0-9_-]+/g) || []).length;
}

// ---------------------------------------------------------------------------
// A. editor/MolEditor.js — the five tools are real (the fix surfaces real chrome).
// ---------------------------------------------------------------------------
test('A1 MolEditor defines the R/S, E/Z, H, Ar display toggles and the Pan Mol tool', function () {
    assert.ok(/id:\s*'ciprs',\s*label:\s*'R\/S'/.test(MOL), "R/S toggle (id:'ciprs') exists");
    assert.ok(/id:\s*'cipez',\s*label:\s*'E\/Z'/.test(MOL), "E/Z toggle (id:'cipez') exists");
    assert.ok(/id:\s*'toggleh',\s*label:\s*'H'/.test(MOL), "H toggle (id:'toggleh') exists");
    assert.ok(/id:\s*'togglearo',\s*label:\s*'Ar'/.test(MOL), "Ar toggle (id:'togglearo') exists");
    assert.ok(/id:\s*'panmol',\s*label:\s*'Pan Mol'/.test(MOL), "Pan Mol tool (id:'panmol') exists");
});

// ---------------------------------------------------------------------------
// B. css/workbench.css — the continuum hide baseline is kept (lens + pathway).
// ---------------------------------------------------------------------------
test('B1 the baseline still HIDES all five row-2 affordances in the editor chrome', function () {
    // The 5-selector hide block ending in display:none !important.
    var re = /\.wb-editor-wrap\s+\.bime-toolbar-row-2\s+\.bime-btn\[data-tool-id="panmol"\]\s*,[\s\S]*?\[data-tool-id="togglearo"\]\s*\{[^}]*display:\s*none\s*!important/;
    assert.ok(re.test(CSS), 'the panmol..togglearo hide rule (display:none !important) is kept');
});

// ---------------------------------------------------------------------------
// C. css/workbench.css — the v2.4.8 scoped un-hide (the fix).
// ---------------------------------------------------------------------------
test('C1 the Editor view RE-SHOWS all five affordances (display:flex !important)', function () {
    TOOLS.forEach(function (id) {
        var re = new RegExp('\\.wb-main\\.wb-view-editor\\s+\\.wb-editor-wrap\\s+\\.bime-toolbar-row-2\\s+\\.bime-btn\\[data-tool-id="' + id + '"\\]');
        assert.ok(re.test(CSS), 'the Editor-view override includes ' + id);
    });
    // The override block resolves to display:flex !important (matches visible siblings).
    var blockRe = /\.wb-main\.wb-view-editor[^{]*\[data-tool-id="togglearo"\]\s*\{[^}]*display:\s*flex\s*!important/;
    assert.ok(blockRe.test(CSS), 'the override sets display:flex !important');
});

test('C2 the override is keyed on .wb-view-editor and is strictly more specific than the hide', function () {
    var hideSel = '.wb-editor-wrap .bime-toolbar-row-2 .bime-btn[data-tool-id="togglearo"]';
    var overSel = '.wb-main.wb-view-editor .wb-editor-wrap .bime-toolbar-row-2 .bime-btn[data-tool-id="togglearo"]';
    assert.ok(CSS.indexOf(hideSel) !== -1, 'hide selector present');
    assert.ok(CSS.indexOf(overSel) !== -1, 'override selector present');
    assert.ok(selectorWeight(overSel) > selectorWeight(hideSel),
        'override is more specific (' + selectorWeight(overSel) + ' > ' + selectorWeight(hideSel) + ')');
    assert.ok(overSel.indexOf(hideSel) === overSel.length - hideSel.length,
        'the override is the hide selector with .wb-main.wb-view-editor prepended');
});

// ---------------------------------------------------------------------------
// D. SCOPE GUARDS — only these five; actions/export stay hidden; lens unaffected.
// ---------------------------------------------------------------------------
test('D1 the Editor view does NOT re-show the actions or export rail groups', function () {
    var reAct = /\.wb-view-editor[^{}]*\[data-group="actions"\][^{}]*\{[^}]*display:\s*flex/;
    var reExp = /\.wb-view-editor[^{}]*\[data-group="export"\][^{}]*\{[^}]*display:\s*flex/;
    assert.strictEqual(reAct.test(CSS), false, 'actions group stays hidden (Validate is in Editor Output)');
    assert.strictEqual(reExp.test(CSS), false, 'export group stays hidden (Editor Output + Export Presets cover it)');
});

test('D2 the lens contract holds: no .wb-lens-open rule re-shows these affordances', function () {
    TOOLS.forEach(function (id) {
        var re = new RegExp('\\.wb-lens-open[^{}]*\\[data-tool-id="' + id + '"\\][^{}]*\\{[^}]*display:\\s*flex');
        assert.strictEqual(re.test(CSS), false, 'no .wb-lens-open rule re-shows ' + id + ' (lens stays decluttered)');
    });
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
