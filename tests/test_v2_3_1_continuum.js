/**
 * tests/test_v2_3_1_continuum.js — focus-lens continuum (v2.3.1, merge Stage 3).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Stage 3 retires the Molecule/Pathway MODE TOGGLE. The pathway canvas is now
 * ONE continuous surface — the permanent workspace — and structure editing
 * happens through the Stage-2 focus lens (the molecule editor floats over a
 * node) rather than a full-view swap.
 *
 * setWorkbenchMode + syncWorkbenchLayoutMode survive as COMPATIBILITY SHIMS so
 * their callers (and the source-shape tests that pin `setWorkbenchMode('molecule')`)
 * keep working. This suite pins the continuum shape of js/workbench.js +
 * css/workbench.css without needing a browser DOM:
 *   (a) both shims are still defined and callable;
 *   (b) the pinned single-molecule entry points still carry the literal
 *       setWorkbenchMode('molecule') (the legacy fallback the older suites count);
 *   (c) editPathwayNodeStructure routes to the focus lens (openStructureLens);
 *   (d) the single-molecule workflow drops a fresh node + opens the lens; and
 *   (e) the continuum CSS marker (the .wb-lens-open editor-overlay rule) exists,
 *       while the dead full-screen-editor swap is gone.
 *
 * Plain Node, no DOM (read files as text), mirroring test_v2_3_0 sections D/E.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Focus-lens continuum (v2.3.1)');
var test = runner.test;
console.log('Focus-lens continuum (v2.3.1)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var CSS = readRepoFile('css/workbench.css');
var HTML = readRepoFile('workbench.html');

// Slice a top-level function body from `function name(` up to the next
// top-level `\nfunction ` — the same slicing test_v2_0_67 / test_v2_3_0 use, so
// a refactor that splits a pinned body fails here too.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// ---------------------------------------------------------------------------
// A. The mode functions survive as compatibility shims.
// ---------------------------------------------------------------------------
test('A1 setWorkbenchMode + syncWorkbenchLayoutMode are still defined', function () {
    assert.ok(WB.indexOf('function setWorkbenchMode(') !== -1,
        'setWorkbenchMode symbol is preserved (compat shim)');
    assert.ok(WB.indexOf('function syncWorkbenchLayoutMode(') !== -1,
        'syncWorkbenchLayoutMode symbol is preserved (compat shim)');
});

test('A2 the dispatcher still routes set-workbench-mode through the shim', function () {
    assert.ok(WB.indexOf("action === 'set-workbench-mode'") !== -1
        && WB.indexOf("setWorkbenchMode(target.getAttribute('data-mode'))") !== -1,
        'set-workbench-mode action stays resolvable via the dispatcher');
    // The action + data-mode wiring must remain in the markup even though the
    // toggle is hidden (release_integrity also pins this).
    assert.ok(HTML.indexOf('data-wb-action="set-workbench-mode"') !== -1
        && HTML.indexOf('data-mode="pathway"') !== -1,
        'set-workbench-mode + data-mode attributes remain in workbench.html');
});

test('A3 setWorkbenchMode no longer performs the legacy full-view swap', function () {
    var fn = wbFn('setWorkbenchMode');
    assert.ok(fn, 'setWorkbenchMode located');
    // The retired swap hid the pathway section when returning to "molecule".
    assert.strictEqual(fn.indexOf("setWorkbenchSectionOpen('wb-pathway-section', false)"), -1,
        'the shim no longer hides the pathway workspace (no mode swap)');
    // The canvas is the workspace in every mode now.
    assert.ok(fn.indexOf("openWorkbenchSection('wb-pathway-section')") !== -1,
        'the shim keeps the canvas workspace mounted');
});

// ---------------------------------------------------------------------------
// B. Pinned single-molecule entry points keep the legacy substring.
// ---------------------------------------------------------------------------
test('B1 addCurrentMoleculeToPathway still contains setWorkbenchMode(\'molecule\')', function () {
    var fn = wbFn('addCurrentMoleculeToPathway');
    assert.ok(fn, 'addCurrentMoleculeToPathway located');
    assert.ok(fn.indexOf("setWorkbenchMode('molecule')") !== -1,
        'legacy fallback substring preserved (test_v2_0_52 counts it)');
});

test('B2 editPathwayNodeStructure still contains setWorkbenchMode(\'molecule\') + load(source)', function () {
    var fn = wbFn('editPathwayNodeStructure');
    assert.ok(fn, 'editPathwayNodeStructure located');
    assert.ok(fn.indexOf("setWorkbenchMode('molecule')") !== -1,
        'legacy fallback substring preserved (test_v2_0_60 counts it)');
    assert.ok(fn.indexOf('pathwayNodeStructureInput(node)') !== -1,
        'still chooses the saved RXN/MOL/SMILES payload (pinned)');
    assert.ok(fn.indexOf('load(source)') !== -1,
        'still loads through the normal molecular-input path (pinned)');
});

// ---------------------------------------------------------------------------
// C. editPathwayNodeStructure routes to the focus lens (continuum).
// ---------------------------------------------------------------------------
test('C1 editPathwayNodeStructure reaches openStructureLens(', function () {
    var fn = wbFn('editPathwayNodeStructure');
    assert.ok(fn.indexOf('openStructureLens(') !== -1,
        'edit-in-place opens the focus lens on the continuum');
    // It opens the lens on the current selection, returning before the legacy
    // full-editor load path runs.
    assert.ok(fn.indexOf('openStructureLens(_pathway.selectedId)') !== -1,
        'lens opens on the selected node id');
    var lensIdx = fn.indexOf('openStructureLens(');
    var loadIdx = fn.indexOf('load(source)');
    assert.ok(lensIdx !== -1 && loadIdx !== -1 && lensIdx < loadIdx,
        'lens route precedes the legacy load() fallback');
});

// ---------------------------------------------------------------------------
// D. The single-molecule "draw one structure" flow opens the lens on a fresh
//    node (so "draw one structure -> read its SMILES" survives the continuum).
// ---------------------------------------------------------------------------
test('D1 the empty-editor "Draw Structure" branch drops a fresh node into the lens', function () {
    var helper = wbFn('newPathwayNodeInLens');
    assert.ok(helper, 'newPathwayNodeInLens helper is defined');
    assert.ok(helper.indexOf('addPathwayNode(') !== -1,
        'reuses node creation for the fresh standalone structure');
    assert.ok(helper.indexOf('openStructureLens(') !== -1,
        'opens the focus lens on the fresh node');
    // The Draw-Structure entry delegates to it when the editor is empty.
    var add = wbFn('addCurrentMoleculeToPathway');
    assert.ok(add.indexOf('newPathwayNodeInLens()') !== -1,
        'empty-editor branch delegates to the lens entry point');
});

// ---------------------------------------------------------------------------
// E. The continuum CSS marker exists; the dead editor-hero swap is gone.
// ---------------------------------------------------------------------------
test('E1 the .wb-lens-open editor-overlay rule is the only editor-visible state', function () {
    assert.ok(CSS.indexOf('.wb-main.wb-lens-open .wb-editor-wrap') !== -1,
        'the focus-lens overlay rule (continuum marker) exists');
    // The editor chrome is out of flow except as the lens overlay.
    assert.ok(/\.wb-editor-wrap\s*\{[^}]*display:\s*none/.test(CSS),
        'the editor chrome is display:none in the base continuum layout');
});

test('E2 the legacy full-screen-editor swap layout is retired', function () {
    // The old swap hid the editor/output/inspector behind .is-pathway-open and
    // built an "editor" grid area. Neither should survive.
    assert.strictEqual(CSS.indexOf('grid-area: editor'), -1,
        'no dead "editor" grid area remains');
    assert.strictEqual(CSS.indexOf('"editor"'), -1,
        'no grid-template-areas still place an editor hero');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
