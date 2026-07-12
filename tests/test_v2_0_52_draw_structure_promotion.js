/**
 * tests/test_v2_0_52_draw_structure_promotion.js — "Draw Structure"
 * affordance surfaced in the pathway Tools row + smart-empty handling
 * (v2.0.52).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.17 turned Pathway mode into a TRUE single-view swap — the molecule
 * editor (`.wb-editor-wrap`) becomes `display: none` while the pathway
 * canvas takes over the layout. That left the user with no visible
 * drawing surface in pathway-only mode, and the existing "Add Structure"
 * action button was parked in a command grid below the canvas where it
 * was easy to miss.
 *
 * v2.0.52 surfaces the same action ("add-pathway-current") in the pathway
 * Tools toolbar — adjacent to Select / Metabolite / Residue / Arrow /
 * Curly / Step / Note / Compartment — and renames the surfaced affordance
 * to "Draw Structure" so the intent reads correctly when the editor is
 * hidden. The companion behaviour: when the editor's molecule has zero
 * atoms, the click auto-switches back to Molecule mode so the user has
 * a drawing surface, with a status hint pointing back at the same
 * action.
 *
 * Plain Node, no DOM. The tests assert the markup + workbench.js source
 * shape rather than executing the click handler (which would need a
 * full DOM + bundle harness).
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Draw Structure promotion (v2.0.52)');
var test = runner.test;

console.log('Draw Structure promotion (v2.0.52)');

// Helpers ----------------------------------------------------------------

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. markup: the "Draw Structure" affordance lives in the Tools toolbar
// ---------------------------------------------------------------------------
test('A1 workbench.html — Tools toolbar contains a "Draw Structure" button bound to add-pathway-current', function() {
    var html = readRepoFile('workbench.html');
    // Locate the Tools toolbar block.
    var toolsStart = html.indexOf('wb-pathway-tools" role="toolbar"');
    assert.ok(toolsStart > 0, 'pathway tools toolbar block is present');
    var toolsEnd = html.indexOf('</div>', toolsStart);
    var toolsBlock = html.substring(toolsStart, toolsEnd);
    // The promoted button must be inside the Tools block.
    assert.ok(toolsBlock.indexOf('Draw Structure') !== -1,
        'Tools toolbar carries the "Draw Structure" label');
    assert.ok(toolsBlock.indexOf('data-wb-action="add-pathway-current"') !== -1,
        'Tools toolbar wires the new button to add-pathway-current');
});

test('A2 workbench.html — the promoted button uses the solid primary-action styling, not the tool-toggle styling', function() {
    var html = readRepoFile('workbench.html');
    // The "Draw Structure" button is unique in the file by its label; locate
    // it directly and read the enclosing <button ...> open tag.
    var drawIdx = html.indexOf('>Draw Structure<');
    assert.ok(drawIdx > 0, '"Draw Structure" button label is present');
    var openIdx = html.lastIndexOf('<button', drawIdx);
    assert.ok(openIdx > 0 && openIdx < drawIdx, 'matching <button> open tag is present');
    var openTag = html.substring(openIdx, drawIdx);
    // Must be a primary-action class.
    assert.ok(openTag.indexOf('wb-btn-solid') !== -1,
        'promoted button uses wb-btn-solid (primary action); openTag was: ' + openTag);
    assert.ok(openTag.indexOf('wb-pathway-stamp-action') !== -1,
        'promoted button carries the v2.0.52 stamp-action marker class');
    // The button must NOT be a wb-pathway-tool toggle (no aria-pressed
    // semantic) — those are click-on-canvas tools, this is a one-click
    // action. Test the exact class token to avoid the substring being
    // matched as part of "wb-pathway-tool" inside a different attribute.
    var classMatch = /class="([^"]+)"/.exec(openTag);
    assert.ok(classMatch, 'open tag has a class attribute');
    var classList = classMatch[1].split(/\s+/);
    assert.strictEqual(classList.indexOf('wb-pathway-tool'), -1,
        'promoted button is NOT a wb-pathway-tool toggle');
});

test('A3 workbench.html — the command-grid add-pathway-current button is preserved (redundant placement is intentional)', function() {
    var html = readRepoFile('workbench.html');
    // v2.0.52 promotes the affordance to the Tools row but leaves the
    // command-grid button in place — users who look in either spot find it.
    var occurrences = (html.match(/data-wb-action="add-pathway-current"/g) || []).length;
    assert.strictEqual(occurrences, 2,
        'add-pathway-current is wired in BOTH the Tools row and the command grid (v2.0.52)');
    // v2.3.5 (P1-3): the same action now reads "Draw Structure" EVERYWHERE
    // (the command-grid copy was unified from the old "Add Structure" so one
    // action no longer wears two labels). Both placements carry the label.
    var labelOccurrences = (html.match(/>Draw Structure</g) || []).length;
    assert.strictEqual(labelOccurrences, 2,
        'the unified "Draw Structure" label appears on BOTH add-pathway-current buttons');
});

// ---------------------------------------------------------------------------
// B. CSS: the new affordance is styled distinctly from the toggle tools
// ---------------------------------------------------------------------------
test('B1 css/workbench.css — adds the .wb-pathway-tools-sep visual separator', function() {
    var css = readRepoFile('css/workbench.css');
    assert.ok(css.indexOf('.wb-pathway-tools-sep') !== -1,
        'wb-pathway-tools-sep class is declared');
});

test('B2 css/workbench.css — adds the .wb-pathway-stamp-action emphasis style', function() {
    var css = readRepoFile('css/workbench.css');
    assert.ok(css.indexOf('.wb-pathway-stamp-action') !== -1,
        'wb-pathway-stamp-action class is declared (visual emphasis for the primary action)');
});

// ---------------------------------------------------------------------------
// C. js/workbench.js: smart-empty auto-switch to Molecule mode
// ---------------------------------------------------------------------------
test('C1 addCurrentMoleculeToPathway — empty editor branch calls setWorkbenchMode("molecule")', function() {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function addCurrentMoleculeToPathway()');
    assert.ok(fnStart > 0, 'addCurrentMoleculeToPathway is defined');
    // Use a generous slice so we capture the full empty-branch body.
    var fnSlice = js.substring(fnStart, fnStart + 1500);
    assert.ok(fnSlice.indexOf("setWorkbenchMode('molecule')") !== -1,
        'empty editor branch auto-switches back to Molecule mode (v2.0.52)');
});

test('C2 addCurrentMoleculeToPathway — empty branch shows the v2.0.52 hint pointing back at "Draw Structure"', function() {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function addCurrentMoleculeToPathway()');
    var fnSlice = js.substring(fnStart, fnStart + 1500);
    // Hint mentions "Draw Structure" so the user knows the same action will
    // pick the new molecule up once they've drawn it.
    assert.ok(fnSlice.indexOf('Draw Structure') !== -1,
        'empty-branch status hint references the "Draw Structure" button by name');
});

test('C3 addCurrentMoleculeToPathway — pre-v2.0.52 "Draw or load a structure first" hint is replaced', function() {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function addCurrentMoleculeToPathway()');
    var fnSlice = js.substring(fnStart, fnStart + 1500);
    // The old early-return hint shouldn't ship anymore — the new flow is
    // auto-switch + clearer message. Guards against a regression where the
    // old message accidentally lingers after the refactor.
    assert.strictEqual(fnSlice.indexOf('Draw or load a structure first'), -1,
        'pre-v2.0.52 status hint is gone (replaced by auto-switch + hint)');
});

// ---------------------------------------------------------------------------
// D. wiring sanity — the dispatcher already routes add-pathway-current
// ---------------------------------------------------------------------------
test('D1 js/workbench.js — dispatcher routes add-pathway-current to addCurrentMoleculeToPathway', function() {
    var js = readRepoFile('js/workbench.js');
    // Single canonical wire: action 'add-pathway-current' -> function call.
    var idx = js.indexOf("action === 'add-pathway-current'");
    assert.ok(idx > 0, 'dispatcher case is present');
    var slice = js.substring(idx, idx + 200);
    assert.ok(slice.indexOf('addCurrentMoleculeToPathway()') !== -1,
        'dispatcher invokes addCurrentMoleculeToPathway()');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
