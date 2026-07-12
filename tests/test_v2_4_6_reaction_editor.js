/**
 * tests/test_v2_4_6_reaction_editor.js — the Reaction group is reachable in the
 * dedicated Editor view (v2.4.6).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * THE BUG: the continuum redesign (v2.3.x) force-hid the molecule editor's left-
 * rail Reaction group — `.wb-editor-wrap .bime-left-rail .bime-group[data-group=
 * "reaction"] { display: none !important }` — because back then the editor only
 * ever appeared as a transient in-place FOCUS LENS over a single structure node
 * inside a pathway, where a reaction arrow is a pathway-level concept, not a
 * node-level one. v2.4.0 then added a DEDICATED standalone Editor view (the
 * `.wb-view-editor` hero surface) whose entire purpose is drawing molecules AND
 * reactions (educts -> products, split by the -> arrow) — but the old hide still
 * fired there, so the Reaction group (Arrow / Atom-map / Auto-map / role Colors /
 * mapping Pairs) was invisible. A user could not see any reaction-creation arrow.
 *
 * THE FIX (CSS-only, surgical): re-show the Reaction group ONLY in the dedicated
 * Editor view, via an override keyed on `.wb-main.wb-view-editor` that is strictly
 * MORE specific than the hide AND can never match the lens (which carries
 * `.wb-lens-open`, never `.wb-view-editor`). `actions` + `export` deliberately
 * stay hidden — the workbench's own Editor Output row already exposes Validate +
 * SVG/PNG/Add-to-pathway, so they have alternative access; the Reaction group was
 * the ONLY rail group with no other route, which is exactly why it was the bug.
 *
 * This suite pins all of that as a source-shape contract (no browser DOM), the
 * same level test_v2_4_0 section C uses:
 *   A. editor/MolEditor.js — the Reaction group + Arrow TOOL still exist in the
 *      data model, so the fix surfaces a REAL tool (not dead chrome).
 *   B. css/workbench.css — the continuum hide baseline is KEPT (lens + pathway).
 *   C. css/workbench.css — the v2.4.6 scoped un-hide exists, is keyed on
 *      `.wb-view-editor`, uses !important, and is strictly more specific than the
 *      hide (so it wins in the Editor view, loses everywhere else).
 *   D. SCOPE GUARDS — the override re-shows ONLY reaction (NOT actions/export),
 *      and the lens/pathway states keep reaction hidden by construction.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Reaction group reachable in the Editor view (v2.4.6)');
var test = runner.test;
console.log('Reaction group reachable in the Editor view (v2.4.6)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// Strip CSS block comments before any rule-presence regex: this file's own fix
// comment QUOTES selectors (it names `.wb-lens-open` to explain why the override
// can't match the lens), and a loose `[^{}]*` regex would otherwise match that
// prose instead of a real rule. We assert on RULES, never on comment text.
var CSS = readRepoFile('css/workbench.css').replace(/\/\*[\s\S]*?\*\//g, '');
var MOL = readRepoFile('editor/MolEditor.js');

// Count the compound (class/attribute/id) selectors in a CSS selector string —
// a faithful-enough specificity proxy: both the hide and the override end in the
// SAME `.bime-group[data-group="reaction"]` tail and BOTH use !important, so the
// one with MORE leading compound selectors wins. (No pseudo-classes are used in
// either selector, so counting classes + attrs + ids is exact here.)
function selectorWeight(sel) {
    var classes = (sel.match(/\.[A-Za-z0-9_-]+/g) || []).length;
    var attrs = (sel.match(/\[[^\]]+\]/g) || []).length;
    var ids = (sel.match(/#[A-Za-z0-9_-]+/g) || []).length;
    return classes + attrs + ids;
}

// ---------------------------------------------------------------------------
// A. editor/MolEditor.js — the Reaction group + Arrow tool are REAL.
// ---------------------------------------------------------------------------
test('A1 MolEditor still defines a left-placed Reaction group', function () {
    assert.ok(/id:\s*'reaction',\s*label:\s*'Reaction',\s*placement:\s*'left'/.test(MOL),
        "TOOLBAR_GROUPS carries { id:'reaction', label:'Reaction', placement:'left' }");
});

test('A2 the Reaction group contains the Arrow tool that draws the -> arrow', function () {
    // The fix surfaces a tool that actually creates the reaction arrow; pin it.
    assert.ok(/id:\s*'reaction',\s*label:\s*'Arrow',\s*icon:\s*'react',\s*type:\s*'tool'/.test(MOL),
        "the Arrow tool { id:'reaction', label:'Arrow', icon:'react', type:'tool' } exists");
    assert.ok(MOL.indexOf('Reaction arrow') !== -1 && MOL.indexOf('reactants and products') !== -1,
        'its tooltip describes drawing a reaction arrow between reactants and products');
});

// ---------------------------------------------------------------------------
// B. css/workbench.css — the continuum hide baseline is kept.
// ---------------------------------------------------------------------------
test('B1 the baseline still HIDES reaction/actions/export in the editor rail', function () {
    // This is the lens + pathway contract: in those states the editor chrome must
    // stay decluttered. The override in section C only lifts it for .wb-view-editor.
    var re = /\.wb-editor-wrap\s+\.bime-left-rail\s+\.bime-group\[data-group="reaction"\]\s*,[\s\S]*?\[data-group="actions"\]\s*,[\s\S]*?\[data-group="export"\]\s*\{[^}]*display:\s*none\s*!important/;
    assert.ok(re.test(CSS),
        'the 3-group hide rule (reaction, actions, export) -> display:none !important is kept');
});

// ---------------------------------------------------------------------------
// C. css/workbench.css — the v2.4.6 scoped un-hide (the fix).
// ---------------------------------------------------------------------------
test('C1 the Editor view RE-SHOWS the Reaction group (display:flex !important)', function () {
    var re = /\.wb-main\.wb-view-editor\s+\.wb-editor-wrap\s+\.bime-left-rail\s+\.bime-group\[data-group="reaction"\]\s*\{[^}]*display:\s*flex\s*!important/;
    assert.ok(re.test(CSS),
        '.wb-main.wb-view-editor ... .bime-group[data-group="reaction"] { display: flex !important } exists');
});

test('C2 the override is keyed on .wb-view-editor (so it cannot match the lens)', function () {
    // The focus lens carries .wb-lens-open, NEVER .wb-view-editor; keying the
    // un-hide on .wb-view-editor is what keeps the lens hide intact (section D).
    var i = CSS.indexOf('.wb-main.wb-view-editor .wb-editor-wrap .bime-left-rail .bime-group[data-group="reaction"]');
    assert.ok(i !== -1, 'the override selector is anchored on .wb-main.wb-view-editor');
    var block = CSS.slice(i, CSS.indexOf('}', i) + 1);
    assert.ok(/!important/.test(block),
        'the override uses !important (required to beat the display:none !important hide)');
});

test('C3 the override is strictly MORE specific than the hide (so it wins in the Editor view)', function () {
    var hideSel = '.wb-editor-wrap .bime-left-rail .bime-group[data-group="reaction"]';
    var overSel = '.wb-main.wb-view-editor .wb-editor-wrap .bime-left-rail .bime-group[data-group="reaction"]';
    assert.ok(CSS.indexOf(hideSel) !== -1, 'hide selector present');
    assert.ok(CSS.indexOf(overSel) !== -1, 'override selector present');
    assert.ok(selectorWeight(overSel) > selectorWeight(hideSel),
        'override has higher specificity (' + selectorWeight(overSel) + ' > ' + selectorWeight(hideSel) + ')');
    // And the override's tail IS the hide selector (it only PREPENDS .wb-main.wb-view-editor).
    assert.ok(overSel.indexOf(hideSel) === overSel.length - hideSel.length,
        'the override is the hide selector with .wb-main.wb-view-editor prepended (a strict refinement)');
});

// ---------------------------------------------------------------------------
// D. SCOPE GUARDS — only reaction is lifted; actions/export stay hidden, and the
//    lens/pathway states keep reaction hidden by construction.
// ---------------------------------------------------------------------------
test('D1 the Editor view does NOT re-show actions or export (they have other access)', function () {
    // Validate lives in the Editor Output row; SVG/PNG/etc. too. Re-showing the
    // rail copies would duplicate. Assert no .wb-view-editor un-hide for them.
    var reAct = /\.wb-view-editor[^{}]*\[data-group="actions"\][^{}]*\{[^}]*display:\s*flex/;
    var reExp = /\.wb-view-editor[^{}]*\[data-group="export"\][^{}]*\{[^}]*display:\s*flex/;
    assert.strictEqual(reAct.test(CSS), false, 'the actions group is NOT re-shown in the Editor view');
    assert.strictEqual(reExp.test(CSS), false, 'the export group is NOT re-shown in the Editor view');
});

test('D2 the un-hide is reaction-ONLY (exactly one .wb-view-editor reaction-group override)', function () {
    var m = CSS.match(/\.wb-main\.wb-view-editor\s+\.wb-editor-wrap\s+\.bime-left-rail\s+\.bime-group\[data-group="reaction"\]/g) || [];
    assert.strictEqual(m.length, 1, 'exactly one Editor-view reaction-group override (no scope creep)');
});

test('D3 the lens contract is intact: no .wb-lens-open un-hide leaks reaction in', function () {
    // The lens must keep the reaction group hidden. There must be NO rule that
    // re-shows the reaction group while .wb-lens-open is set.
    var reLens = /\.wb-lens-open[^{}]*\[data-group="reaction"\][^{}]*\{[^}]*display:\s*flex/;
    assert.strictEqual(reLens.test(CSS), false,
        'no .wb-lens-open rule re-shows the reaction group (the lens stays decluttered)');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
