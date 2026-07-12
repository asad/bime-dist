/**
 * tests/test_v2_0_60_edit_pathway_structure.js — edit-in-place
 * structure for pathway nodes (v2.0.60).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background. v2.0.52 promoted "Draw Structure" into the pathway Tools
 * row, but the action only ADDED a fresh node — there was no flow to
 * REVISE an existing node's structure (the user-reported scenario: a
 * pathway node landed with a wrong / missing structure and needed to
 * be corrected without losing its position or edges).
 *
 * v2.0.60 turns "Draw Structure" into an edit-in-place action when a
 * pathway node is selected, AND adds a new "Edit Structure" action
 * that loads the selected node's SMILES back into the molecule
 * editor for revision.
 *
 * Plain Node, no DOM. Tests assert source shape rather than executing
 * the click handler.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Edit pathway structure (v2.0.60)');
var test = runner.test;

console.log('Edit pathway structure (v2.0.60)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. addCurrentMoleculeToPathway updates the selected node in place
// ---------------------------------------------------------------------------
test('A1 addCurrentMoleculeToPathway branches on _pathway.selectedType==="node"', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function addCurrentMoleculeToPathway');
    assert.ok(fnStart > 0, 'addCurrentMoleculeToPathway is defined');
    var slice = js.substring(fnStart, fnStart + 4000);
    assert.ok(slice.indexOf("_pathway.selectedType === 'node'") !== -1,
        'function checks for a selected node');
    assert.ok(slice.indexOf("findPathwayItem('node', _pathway.selectedId)") !== -1,
        'function looks up the selected node by id');
});

test('A2 in-place branch applies the full structure payload to the selected node', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function addCurrentMoleculeToPathway');
    var slice = js.substring(fnStart, fnStart + 4000);
    assert.ok(slice.indexOf('captureEditorStructureForPathway()') !== -1,
        'captures a canonical editor structure payload');
    assert.ok(slice.indexOf('applyPathwayStructurePayload(node, payload)') !== -1,
        'in-place branch updates the selected node via one payload adapter');
    assert.ok(slice.indexOf('renderPathwayCanvas()') !== -1,
        'in-place branch re-renders the canvas');
    assert.ok(slice.indexOf("pathwayStatus('Updated '") !== -1,
        'in-place branch emits a status message naming the update');
});

test('A2b payload stores full molecule/reaction text separate from display subtitle', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('function normalizePathwayStructurePayload') !== -1,
        'pathway nodes keep a structured payload');
    assert.ok(js.indexOf('smiles: smiles') !== -1,
        'payload stores full SMILES separately from subtitle');
    assert.ok(js.indexOf('rxn: rxnfile') !== -1,
        'payload stores RXN text for reactions');
    assert.ok(js.indexOf('PATHWAY_STRUCTURE_TEXT_LIMIT') !== -1,
        'payload has a local file-size guard');
});

test('A3 in-place branch returns early — does NOT also stamp a new node', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function addCurrentMoleculeToPathway');
    var slice = js.substring(fnStart, fnStart + 4000);
    // The branch ends with `return;` before the addPathwayNode fall-through.
    var inPlaceIdx = slice.indexOf("_pathway.selectedType === 'node'");
    var addPathwayIdx = slice.indexOf('addPathwayNode(center.x, center.y');
    assert.ok(inPlaceIdx > 0 && addPathwayIdx > inPlaceIdx,
        'in-place branch appears before the fall-through addPathwayNode call');
    // Between the two, there must be at least one `return;` (the early
    // return inside the in-place branch).
    var between = slice.substring(inPlaceIdx, addPathwayIdx);
    assert.ok(between.indexOf('return;') !== -1,
        'in-place branch returns before falling through to addPathwayNode');
});

// ---------------------------------------------------------------------------
// B. editPathwayNodeStructure: load selected node's saved chemistry into the editor
// ---------------------------------------------------------------------------
test('B1 editPathwayNodeStructure is defined', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('function editPathwayNodeStructure') !== -1,
        'editPathwayNodeStructure is declared');
});

test('B2 editPathwayNodeStructure guards on selection + loads saved structure into the editor', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function editPathwayNodeStructure');
    var slice = js.substring(fnStart, fnStart + 3000);
    // Guard against no-selection.
    assert.ok(slice.indexOf("_pathway.selectedType !== 'node'") !== -1,
        'guards on a node being selected');
    assert.ok(slice.indexOf('pathwayNodeStructureInput(node)') !== -1,
        'chooses the saved RXN/MOL/SMILES payload');
    assert.ok(slice.indexOf('load(source)') !== -1,
        'loads through the normal workbench molecular-input path');
    assert.strictEqual(slice.indexOf('editor.loadSmiles'), -1,
        'does not depend on a non-existent editor.loadSmiles helper');
    // Switches back to Molecule mode so the editor canvas is visible.
    assert.ok(slice.indexOf("setWorkbenchMode('molecule')") !== -1,
        'switches to Molecule mode so the editor is visible');
});

test('B3 node structure input prefers reaction RXN, then MOL, then SMILES', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function pathwayNodeStructureInput');
    assert.ok(fnStart > 0, 'pathwayNodeStructureInput is defined');
    var slice = js.substring(fnStart, fnStart + 900);
    assert.ok(slice.indexOf("structure.type === 'reaction' && structure.rxn") !== -1,
        'reaction RXN is preferred when available');
    assert.ok(slice.indexOf('if (structure.mol) return structure.mol') !== -1,
        'MOL is the next editor-load format');
    assert.ok(slice.indexOf('if (structure.smiles) return structure.smiles') !== -1,
        'SMILES remains the portable fallback');
});

test('B4 node structure input never parses display subtitles as chemistry', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function pathwayNodeStructureInput');
    assert.ok(fnStart > 0, 'pathwayNodeStructureInput is defined');
    var slice = js.substring(fnStart, fnStart + 900);
    assert.strictEqual(slice.indexOf('node.subtitle'), -1,
        'display subtitles such as names or accessions are not treated as molecular input');
    assert.ok(slice.indexOf("return '';") !== -1,
        'nodes without saved chemistry fail closed instead of loading a wrong fragment');
});

test('B5 built-in pathway examples store editable structure payloads on nodes', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function loadExamplePathway');
    assert.ok(fnStart > 0, 'loadExamplePathway is defined');
    var slice = js.substring(fnStart, fnStart + 2600);
    assert.ok(slice.indexOf('normalizePathwayStructurePayload') !== -1,
        'example nodes build a real structure payload from metabolite metadata');
    assert.ok(slice.indexOf('smiles: meta.smiles') !== -1,
        'example nodes persist the source SMILES used for the thumbnail');
    assert.ok(slice.indexOf('structure: structure') !== -1,
        'example nodes pass the payload into addPathwayNode');
});

test('B6 double-clicking a pathway node opens structure editing', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf("addEventListener('dblclick', handlePathwayCanvasDoubleClick)") !== -1,
        'double-click handler is bound');
    var fnStart = js.indexOf('function handlePathwayCanvasDoubleClick');
    assert.ok(fnStart > 0, 'double-click handler is declared');
    var slice = js.substring(fnStart, fnStart + 900);
    assert.ok(slice.indexOf("target.type !== 'node'") !== -1,
        'double-click editing is limited to pathway nodes');
    assert.ok(slice.indexOf('editPathwayNodeStructure()') !== -1,
        'double-click delegates to the same edit flow as the context menu');
});

// ---------------------------------------------------------------------------
// C. dispatcher wiring
// ---------------------------------------------------------------------------
test('C1 dispatcher routes edit-pathway-structure to editPathwayNodeStructure', function () {
    var js = readRepoFile('js/workbench.js');
    var idx = js.indexOf("action === 'edit-pathway-structure'");
    assert.ok(idx > 0, 'dispatcher case is present');
    var slice = js.substring(idx, idx + 200);
    assert.ok(slice.indexOf('editPathwayNodeStructure()') !== -1,
        'dispatcher invokes editPathwayNodeStructure()');
});

// ---------------------------------------------------------------------------
// D. workbench.html — buttons + tooltips
// ---------------------------------------------------------------------------
test('D1 workbench.html has Edit-Structure buttons in BOTH the Tools row and the Structure command group', function () {
    var html = readRepoFile('workbench.html');
    var occurrences = (html.match(/data-wb-action="edit-pathway-structure"/g) || []).length;
    assert.strictEqual(occurrences, 2,
        'edit-pathway-structure wired in BOTH the Tools row and the command grid');
    assert.ok(html.indexOf('>Edit Structure<') !== -1,
        'visible "Edit Structure" button label is present');
});

test('D2 Draw Structure button tooltip mentions in-place update on selected node', function () {
    var html = readRepoFile('workbench.html');
    var drawIdx = html.indexOf('>Draw Structure<');
    assert.ok(drawIdx > 0);
    var openIdx = html.lastIndexOf('<button', drawIdx);
    var openTag = html.substring(openIdx, drawIdx);
    assert.ok(openTag.toLowerCase().indexOf('updates that node') !== -1 ||
              openTag.toLowerCase().indexOf('update') !== -1,
        'Draw Structure tooltip references the in-place update flow');
});

test('D3 the command-grid add-pathway-current button is preserved and documents the in-place behaviour', function () {
    var html = readRepoFile('workbench.html');
    // v2.3.5 (P1-3): the command-grid label was unified from "Add Structure"
    // to "Draw Structure" so one action does not read two different ways. The
    // SECOND "Draw Structure" instance is the command-grid button (the first is
    // the Tools-row stamp action). It keeps a tooltip describing the in-place
    // update flow.
    var firstIdx = html.indexOf('>Draw Structure<');
    var addIdx = html.indexOf('>Draw Structure<', firstIdx + 1);
    assert.ok(addIdx > firstIdx,
        'a second (command-grid) "Draw Structure" button is present');
    var openIdx = html.lastIndexOf('<button', addIdx);
    var openTag = html.substring(openIdx, addIdx);
    assert.ok(openTag.indexOf('updates that node') !== -1 ||
              openTag.indexOf('node selected') !== -1,
        'command-grid Draw Structure tooltip documents the in-place behaviour');
});

test('D4 command-grid "Update" button is relabelled to "Update Labels"', function () {
    var html = readRepoFile('workbench.html');
    // The v2.0.52-era button was just labelled "Update", which was
    // ambiguous (update what — structure? label?). v2.0.60 disambiguates.
    assert.ok(html.indexOf('>Update Labels<') !== -1,
        'Update button relabelled to "Update Labels"');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
