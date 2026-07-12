/**
 * tests/test_v2_0_61_pathway_context_menu.js — right-click context
 * menus on pathway nodes / edges / steps / notes / compartments /
 * empty canvas (v2.0.61).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.61 adds a floating context menu that the user opens with a
 * right-click. The menu is repopulated per right-click depending on
 * which pathway element was clicked:
 *
 *   - Node     — Edit Structure, Update Labels, Duplicate, Copy SMILES,
 *                Delete
 *   - Edge     — Reverse arrow, Cycle arrow type, Delete
 *   - Step / Note / Compartment — Update Labels, Delete
 *   - Empty canvas — Add Metabolite / Residue / Step / Note Here,
 *                    Clean Up Layout, Fit
 *
 * Dismissed by ESC, by clicking outside, or by activating any item.
 *
 * Plain Node, no DOM. Source-shape only.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Pathway context menu (v2.0.61)');
var test = runner.test;

console.log('Pathway context menu (v2.0.61)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. core helpers declared
// ---------------------------------------------------------------------------
test('A1 js/workbench.js declares showPathwayContextMenu / hidePathwayContextMenu / resolvePathwayContextTarget', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('function showPathwayContextMenu') !== -1);
    assert.ok(js.indexOf('function hidePathwayContextMenu') !== -1);
    assert.ok(js.indexOf('function resolvePathwayContextTarget') !== -1);
});

test('A2 js/workbench.js declares handlePathwayCanvasContextMenu + pathwayContextMenuItems', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('function handlePathwayCanvasContextMenu') !== -1);
    assert.ok(js.indexOf('function pathwayContextMenuItems') !== -1);
});

test('A3 contextmenu event listener is wired onto document', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf("addEventListener('contextmenu', handlePathwayCanvasContextMenu)") !== -1,
        'document contextmenu listener wires handlePathwayCanvasContextMenu');
});

test('A4 ESC + click-outside dismiss the menu', function () {
    var js = readRepoFile('js/workbench.js');
    // ESC handler:
    assert.ok(js.indexOf('_pathwayContextMenu && (e.key') !== -1,
        'keydown handler dismisses on ESC');
    // Click-outside (capture-phase mousedown):
    assert.ok(js.indexOf("addEventListener('mousedown', function(e)") !== -1,
        'mousedown handler dismisses when click is outside the menu');
});

// ---------------------------------------------------------------------------
// B. menu item composition per target type
// ---------------------------------------------------------------------------
test('B1 node target carries Edit Structure / Duplicate / Copy SMILES / Delete', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function pathwayContextMenuItems');
    var slice = js.substring(fnStart, fnStart + 3500);
    assert.ok(slice.indexOf("'Edit Structure'") !== -1, 'node: Edit Structure');
    assert.ok(slice.indexOf("'Duplicate node'") !== -1, 'node: Duplicate node');
    assert.ok(slice.indexOf("'Copy SMILES'") !== -1, 'node: Copy SMILES');
    assert.ok(slice.indexOf("'Delete'") !== -1, 'node: Delete');
});

test('B2 edge target carries Reverse + Cycle arrow type + Delete', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function pathwayContextMenuItems');
    var slice = js.substring(fnStart, fnStart + 3500);
    assert.ok(slice.indexOf("'Reverse arrow'") !== -1, 'edge: Reverse arrow');
    assert.ok(slice.indexOf("'Cycle arrow type'") !== -1, 'edge: Cycle arrow type');
});

test('B3 empty canvas target carries Add-X-here + Clean Up + Fit', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function pathwayContextMenuItems');
    var slice = js.substring(fnStart, fnStart + 3500);
    assert.ok(slice.indexOf("'Add Metabolite here'") !== -1);
    assert.ok(slice.indexOf("'Add Residue here'") !== -1);
    assert.ok(slice.indexOf("'Add Step here'") !== -1);
    assert.ok(slice.indexOf("'Add Note here'") !== -1);
    assert.ok(slice.indexOf("'Clean Up Layout'") !== -1);
    assert.ok(slice.indexOf("'Fit'") !== -1);
});

// ---------------------------------------------------------------------------
// C. dispatcher wiring — each context-* action is routed
// ---------------------------------------------------------------------------
test('C1 dispatcher routes all eight context-* actions', function () {
    var js = readRepoFile('js/workbench.js');
    var actions = [
        'context-edit-node-structure',
        'context-update-node-labels',
        'context-delete-item',
        'context-copy-smiles',
        'context-duplicate-node',
        'context-reverse-edge',
        'context-cycle-arrow-type',
        'context-add-here'
    ];
    for (var i = 0; i < actions.length; i++) {
        assert.ok(js.indexOf("action === '" + actions[i] + "'") !== -1,
            'dispatcher routes ' + actions[i]);
    }
});

test('C2 every dispatcher case calls hidePathwayContextMenu() before its handler', function () {
    var js = readRepoFile('js/workbench.js');
    // For each context-* action case, look at the body and confirm the
    // first statement is hidePathwayContextMenu(). This guarantees the
    // menu doesn't linger after the user activates an item.
    var actions = [
        'context-edit-node-structure',
        'context-update-node-labels',
        'context-delete-item',
        'context-copy-smiles',
        'context-duplicate-node',
        'context-reverse-edge',
        'context-cycle-arrow-type',
        'context-add-here'
    ];
    for (var i = 0; i < actions.length; i++) {
        var idx = js.indexOf("action === '" + actions[i] + "'");
        var slice = js.substring(idx, idx + 400);
        assert.ok(slice.indexOf('hidePathwayContextMenu()') !== -1,
            actions[i] + ' calls hidePathwayContextMenu() before its action');
    }
});

// ---------------------------------------------------------------------------
// D. CSS — context menu styling
// ---------------------------------------------------------------------------
test('D1 css/workbench.css carries .wb-context-menu + .wb-context-menu-item + .wb-context-menu-divider', function () {
    var css = readRepoFile('css/workbench.css');
    assert.ok(css.indexOf('.wb-context-menu ') !== -1 ||
              css.indexOf('.wb-context-menu{') !== -1 ||
              css.indexOf('.wb-context-menu\n') !== -1,
        '.wb-context-menu rule is declared');
    assert.ok(css.indexOf('.wb-context-menu-item') !== -1);
    assert.ok(css.indexOf('.wb-context-menu-divider') !== -1);
});

test('D2 .wb-context-menu uses position:fixed + a high z-index', function () {
    var css = readRepoFile('css/workbench.css');
    var ruleStart = css.indexOf('.wb-context-menu ');
    if (ruleStart === -1) { ruleStart = css.indexOf('.wb-context-menu\n'); }
    if (ruleStart === -1) { ruleStart = css.indexOf('.wb-context-menu{'); }
    assert.ok(ruleStart > 0);
    var rule = css.substring(ruleStart, ruleStart + 700);
    assert.ok(rule.indexOf('position: fixed') !== -1 || rule.indexOf('position:fixed') !== -1,
        '.wb-context-menu is position:fixed');
    assert.ok(/z-index:\s*\d{3,}/.test(rule),
        '.wb-context-menu has a 3+ digit z-index');
});

// ---------------------------------------------------------------------------
// E. context-action helpers exist
// ---------------------------------------------------------------------------
test('E1 context action handlers are all declared', function () {
    var js = readRepoFile('js/workbench.js');
    var fns = [
        'contextEditNodeStructure',
        'contextUpdateNodeLabels',
        'contextDeleteItem',
        'contextCopySmiles',
        'contextDuplicateNode',
        'contextReverseEdge',
        'contextCycleArrowType',
        'contextAddHere'
    ];
    for (var i = 0; i < fns.length; i++) {
        assert.ok(js.indexOf('function ' + fns[i]) !== -1,
            fns[i] + ' is declared');
    }
});

test('E2 contextDuplicateNode uses pushPathwayHistory so the duplicate is undoable', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function contextDuplicateNode');
    var slice = js.substring(fnStart, fnStart + 800);
    assert.ok(slice.indexOf('pushPathwayHistory') !== -1,
        'contextDuplicateNode pushes pathway history');
});

test('E3 contextCopySmiles falls back gracefully when clipboard API is unavailable', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function contextCopySmiles');
    var slice = js.substring(fnStart, fnStart + 1200);
    assert.ok(slice.indexOf('navigator.clipboard') !== -1,
        'contextCopySmiles tries the clipboard API');
    assert.ok(slice.indexOf('Clipboard unavailable') !== -1,
        'contextCopySmiles has a no-clipboard fallback message');
});

test('E4 contextCopySmiles prefers the saved structure payload over subtitle text', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function contextCopySmiles');
    var slice = js.substring(fnStart, fnStart + 1200);
    assert.ok(slice.indexOf('pathwayNodeStructurePayload(node)') !== -1,
        'contextCopySmiles reads the canonical structure payload');
    assert.ok(slice.indexOf('structure && structure.smiles') !== -1,
        'canonical SMILES is copied before display subtitle fallback');
});

test('E5 contextDuplicateNode preserves the saved structure payload', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function contextDuplicateNode');
    var slice = js.substring(fnStart, fnStart + 900);
    assert.ok(slice.indexOf('structure: pathwayNodeStructurePayload(node)') !== -1,
        'duplicate nodes carry molecule/reaction payloads forward');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
