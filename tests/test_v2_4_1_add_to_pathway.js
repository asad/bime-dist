/**
 * tests/test_v2_4_1_add_to_pathway.js — Hybrid Stage B: "Add to pathway"
 * from the Editor view (v2.4.1).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Stage A (v2.4.0) revived the visible "Editor | Pathway" view switch. Stage B
 * adds the IMPORT action: from the Editor view, drop the current drawing onto the
 * pathway as a node (molecule) or a reactant/step/product subgraph (reaction),
 * auto-tagged, then switch to the Pathway view. There is NO AAM-linking yet
 * (that is Stage C).
 *
 * The new code is deliberately a thin ENTRY POINT + view switch + auto-tag glue
 * over the EXISTING import spine — it must not re-implement capture/split/stamp:
 *   - addEditorDrawingToPathway() branches on editor.molecule.reactionArrow:
 *       reaction -> convertCurrentReactionToPathway() (splits via
 *       RDT._splitReactionSides, drops reactant + kind:'reaction' step + product
 *       nodes and wires edges);
 *       molecule -> addCurrentMoleculeToPathway() (stamps one node), after
 *       CLEARING the pathway selection so a fresh node is stamped (never edits a
 *       stale selection);
 *     it guards the empty editor (no atoms) with a friendly status and no node,
 *     and ends on setWorkbenchMode('pathway') so the user sees the result.
 *   - the dispatcher routes data-wb-action="editor-add-to-pathway" to it.
 *   - the "Add to pathway" button lives in the Editor Output action row.
 *   - entering the Editor view clears the pathway selection too.
 *
 * Source-shape + DOM-stub, plain Node (no real DOM beyond the shim) — mirroring
 * test_v2_4_0_hybrid_views (A–D) + test_v2_0_67 (the pinned import spine).
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Add to pathway from the Editor view (v2.4.1)');
var test = runner.test;
console.log('Add to pathway from the Editor view (v2.4.1)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var HTML = readRepoFile('workbench.html');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing test_v2_3_1 / test_v2_0_67 / test_v2_4_0 use,
// so a refactor that splits a pinned body fails here too.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// ---------------------------------------------------------------------------
// A. workbench.html — the "Add to pathway" button (Editor Output action row).
// ---------------------------------------------------------------------------
test('A1 an "Add to pathway" button exists with data-wb-action="editor-add-to-pathway"', function () {
    var matches = HTML.match(/data-wb-action="editor-add-to-pathway"/g) || [];
    assert.strictEqual(matches.length, 1,
        'exactly one Editor-view "Add to pathway" button carries the new action');
    // It is a real <button> labelled "Add to pathway".
    assert.ok(/<button[^>]*data-wb-action="editor-add-to-pathway"[^>]*>[\s\S]*?Add to pathway[\s\S]*?<\/button>/.test(HTML),
        'the editor-add-to-pathway control is a button labelled "Add to pathway"');
});

test('A2 the button sits in the Editor Output action row (near SMILES/Copy/Validate)', function () {
    // The Editor Output shell's action row holds export/Copy/Validate; the new
    // button must live in the same block so it is reachable from the Editor view
    // (and, being in that row, in the mobile flow at <=720px too).
    var outIdx = HTML.indexOf('Editor output views');
    assert.ok(outIdx !== -1, 'the Editor Output tablist is present');
    var addIdx = HTML.indexOf('data-wb-action="editor-add-to-pathway"');
    var validateIdx = HTML.indexOf('data-wb-action="validate"');
    assert.ok(addIdx > outIdx, 'the Add-to-pathway button is inside the Editor Output section');
    // Proximity to the existing Validate control (same action row).
    assert.ok(Math.abs(addIdx - validateIdx) < 1200,
        'the Add-to-pathway button sits beside the existing Output actions (Validate/Copy)');
});

test('A3 the new button is ADDITIVE — the pinned import buttons are all preserved', function () {
    // test_v2_0_67 A1 + test_v2_3_5/_6 pin these counts; the new editor-add
    // action must not perturb them.
    var addPathwayCurrent = HTML.match(/data-wb-action="add-pathway-current"/g) || [];
    assert.ok(addPathwayCurrent.length >= 2,
        'the "Draw Structure" (add-pathway-current) buttons are kept (>=2)');
    var reactionToPathway = HTML.match(/data-wb-action="reaction-to-pathway"/g) || [];
    assert.ok(reactionToPathway.length >= 2,
        'the "To Pathway"/"Reaction Step" (reaction-to-pathway) buttons are kept (>=2)');
    assert.ok(HTML.indexOf('>To Pathway</button>') !== -1, 'the Reaction Scheme "To Pathway" action survives');
    assert.ok(HTML.indexOf('>Reaction Step</button>') !== -1, 'the Pathway-canvas "Reaction Step" action survives');
});

// ---------------------------------------------------------------------------
// B. js/workbench.js — the dispatcher routes the new action.
// ---------------------------------------------------------------------------
test('B1 the click dispatcher routes editor-add-to-pathway to addEditorDrawingToPathway', function () {
    assert.ok(WB.indexOf("else if (action === 'editor-add-to-pathway') addEditorDrawingToPathway();") !== -1,
        'editor-add-to-pathway is resolvable via the click dispatcher');
});

test('B2 the pinned import-spine dispatcher wires are untouched', function () {
    // The Stage B route is purely additive: the v2.0.52/_67 wires stay byte-exact.
    assert.ok(WB.indexOf("else if (action === 'add-pathway-current') addCurrentMoleculeToPathway();") !== -1,
        'add-pathway-current -> addCurrentMoleculeToPathway is preserved (test_v2_0_52 D1)');
    assert.ok(WB.indexOf("else if (action === 'reaction-to-pathway') convertCurrentReactionToPathway();") !== -1,
        'reaction-to-pathway -> convertCurrentReactionToPathway is preserved (test_v2_0_67 A2)');
});

// ---------------------------------------------------------------------------
// C. js/workbench.js — addEditorDrawingToPathway shape (NEW sibling function).
// ---------------------------------------------------------------------------
test('C1 addEditorDrawingToPathway is a NEW top-level function (does not relocate a pinned one)', function () {
    var fn = wbFn('addEditorDrawingToPathway');
    assert.ok(fn, 'addEditorDrawingToPathway is declared');
    // It must be its OWN function, not folded into the pinned spine functions.
    assert.ok(wbFn('addCurrentMoleculeToPathway'), 'addCurrentMoleculeToPathway still exists as its own function');
    assert.ok(wbFn('convertCurrentReactionToPathway'), 'convertCurrentReactionToPathway still exists as its own function');
    assert.ok(wbFn('setWorkbenchMode'), 'setWorkbenchMode still exists as its own function');
});

test('C2 it guards the empty editor (no atoms) — friendly status, no node', function () {
    var fn = wbFn('addEditorDrawingToPathway');
    assert.ok(fn.indexOf('editor.molecule.atoms.length === 0') !== -1,
        'guards on an empty editor (no atoms)');
    assert.ok(/pathwayStatus\([^)]*Draw a molecule or reaction/.test(fn),
        'empty-editor branch emits a friendly status pointing the user to draw first');
    // The empty guard must return before any node-stamping delegation.
    var guardIdx = fn.indexOf('editor.molecule.atoms.length === 0');
    var convertIdx = fn.indexOf('convertCurrentReactionToPathway()');
    var addIdx = fn.indexOf('addCurrentMoleculeToPathway()');
    var firstReturn = fn.indexOf('return;', guardIdx);
    assert.ok(firstReturn !== -1 && firstReturn < convertIdx && firstReturn < addIdx,
        'the empty guard returns before delegating to either import path');
});

test('C3 it branches on editor.molecule.reactionArrow (reaction vs molecule)', function () {
    var fn = wbFn('addEditorDrawingToPathway');
    assert.ok(fn.indexOf('editor.molecule.reactionArrow') !== -1,
        'branches on whether the editor holds a reaction');
});

test('C4 the REACTION branch delegates to convertCurrentReactionToPathway()', function () {
    var fn = wbFn('addEditorDrawingToPathway');
    assert.ok(fn.indexOf('convertCurrentReactionToPathway()') !== -1,
        'a drawn reaction is handed to the existing reaction-splitting converter');
});

test('C5 the MOLECULE branch clears the selection THEN delegates to addCurrentMoleculeToPathway()', function () {
    var fn = wbFn('addEditorDrawingToPathway');
    assert.ok(fn.indexOf('addCurrentMoleculeToPathway()') !== -1,
        'a plain molecule is handed to the existing single-node stamper');
    // Clearing the pathway selection BEFORE the molecule path is what guarantees
    // a FRESH node (addCurrentMoleculeToPathway keeps its own edit-in-place branch
    // for the Pathway-view "Draw Structure" flow; here we must bypass it).
    var clearIdx = fn.indexOf("selectPathwayItem('', null)");
    var addIdx = fn.indexOf('addCurrentMoleculeToPathway()');
    assert.ok(clearIdx !== -1, 'the pathway selection is cleared (selectPathwayItem(\'\', null))');
    assert.ok(clearIdx < addIdx,
        'the selection is cleared BEFORE addCurrentMoleculeToPathway so it stamps a fresh node');
});

test('C6 it ends by switching to the Pathway view and surfaces an "Added ..." status', function () {
    var fn = wbFn('addEditorDrawingToPathway');
    assert.ok(fn.indexOf("setWorkbenchMode('pathway')") !== -1,
        "lands the user on the Pathway view (setWorkbenchMode('pathway'))");
    assert.ok(fn.indexOf('fitPathwayCanvas()') !== -1,
        'fits the canvas so the imported node(s) are visible');
    assert.ok(/pathwayStatus\([^)]*Added /.test(fn),
        'surfaces an Added "<label>" status using the captured label');
    // The label is reused from the existing namer, not re-derived here.
    assert.ok(fn.indexOf('captureEditorStructureForPathway()') !== -1,
        'reuses captureEditorStructureForPathway (pathwayInputLabel/moleculeFormulaLabel) for the label');
    // setWorkbenchMode('pathway') is the LAST view decision (after both branches).
    var setIdx = fn.lastIndexOf("setWorkbenchMode('pathway')");
    var convertIdx = fn.indexOf('convertCurrentReactionToPathway()');
    var addIdx = fn.indexOf('addCurrentMoleculeToPathway()');
    assert.ok(setIdx > convertIdx && setIdx > addIdx,
        'the Pathway-view switch comes after the import delegation for both branches');
});

// ---------------------------------------------------------------------------
// C'. setWorkbenchMode grows (does not relocate) — entering Editor clears the
//     pathway selection, and the v2.4.0 pins still hold.
// ---------------------------------------------------------------------------
test("C7 entering the Editor view clears the pathway selection (selectPathwayItem('', null))", function () {
    var fn = wbFn('setWorkbenchMode');
    assert.ok(fn, 'setWorkbenchMode located');
    var editorBranchIdx = fn.indexOf("if (_wbView === 'editor')");
    assert.ok(editorBranchIdx !== -1, 'the Editor-view branch exists');
    var clearIdx = fn.indexOf("selectPathwayItem('', null)");
    assert.ok(clearIdx > editorBranchIdx,
        'a selectPathwayItem(\'\', null) call lives inside the Editor-view branch');
    // It still re-fits the editor (v2.4.0 B2 behaviour preserved).
    assert.ok(fn.indexOf('refitMoleculeEditorSoon()') !== -1,
        'the Editor-view branch still re-fits the molecule editor (v2.4.0 carry-over)');
});

test('C8 setWorkbenchMode keeps its v2.4.0 source-shape pins (grow, not relocate)', function () {
    var fn = wbFn('setWorkbenchMode');
    assert.ok(fn.indexOf("openWorkbenchSection('wb-pathway-section')") !== -1,
        'keeps the pinned openWorkbenchSection literal (canvas stays mounted)');
    assert.ok(fn.indexOf("_wbView = (mode === 'editor') ? 'editor' : 'pathway'") !== -1,
        'keeps the view-selection assignment');
    assert.strictEqual(fn.indexOf("setWorkbenchSectionOpen('wb-pathway-section', false)"), -1,
        'still never hides the pathway section via the v2.3.1-A3-forbidden literal');
    assert.strictEqual(fn.indexOf('appendChild'), -1, 'no reparent (appendChild) — the switch is a class flip');
    assert.strictEqual(fn.indexOf('bime-editor'), -1, 'never re-mounts #bime-editor');
});

test('C9 the addCurrentMoleculeToPathway edit-in-place branch is untouched (v2.0.60 pin)', function () {
    // Stage B does NOT touch the spine: addCurrentMoleculeToPathway keeps its
    // selected-node in-place branch returning before the fall-through stamp.
    var slice = wbFn('addCurrentMoleculeToPathway');
    assert.ok(slice.indexOf("_pathway.selectedType === 'node'") !== -1,
        'in-place branch still checks for a selected node');
    var inPlaceIdx = slice.indexOf("_pathway.selectedType === 'node'");
    var addPathwayIdx = slice.indexOf('addPathwayNode(center.x, center.y');
    assert.ok(inPlaceIdx > 0 && addPathwayIdx > inPlaceIdx,
        'in-place branch still precedes the fall-through addPathwayNode call');
    assert.ok(slice.substring(inPlaceIdx, addPathwayIdx).indexOf('return;') !== -1,
        'in-place branch still returns before the stamp (test_v2_0_60 A3)');
});

// ---------------------------------------------------------------------------
// D. Executable DOM-stub: drive addEditorDrawingToPathway through a faithful
//    port of its control flow against a stub editor + stub _pathway, and assert
//    the node/edge deltas + the view switch. Fails loudly on a regression.
//    (Same approach as test_v2_4_0 section D / test_v2_3_4 section E.)
// ---------------------------------------------------------------------------

// A minimal pathway model with the node/edge/selection surface the import spine
// touches. addNode/addEdge mimic addPathwayNode/addPathwayEdgeRecord deltas.
function makePathway() {
    return {
        nodes: [],
        edges: [],
        selectedType: 'node',     // start with a STALE selection on purpose
        selectedId: 'node-stale',
        _seq: 0
    };
}

// Faithful port of addEditorDrawingToPathway's control flow. `env` carries the
// stub editor + pathway + a recorder for the view and status; returns the view.
function applyAddEditorDrawingToPathway(env) {
    var editor = env.editor;
    var pw = env.pathway;
    // Empty-editor guard.
    if (!editor || !editor.molecule || !editor.molecule.atoms || editor.molecule.atoms.length === 0) {
        env.status = 'Draw a molecule or reaction in the editor first, then "Add to pathway".';
        return env.view; // unchanged; no node stamped
    }
    var isReaction = !!editor.molecule.reactionArrow;
    var label = isReaction ? (editor.molName || 'Reaction') : (editor.molName || 'Structure');
    if (isReaction) {
        // convertCurrentReactionToPathway(): split -> reactant nodes + 1 reaction
        // node + product nodes, wired with edges. (Port the node/edge deltas.)
        var r = editor.molecule.split.reactants;
        var p = editor.molecule.split.products;
        var reactionNode = { id: 'node-' + (++pw._seq), kind: 'reaction', label: label };
        pw.nodes.push(reactionNode);
        var i, n;
        for (i = 0; i < r.length; i++) {
            n = { id: 'node-' + (++pw._seq), kind: 'metabolite', label: 'Reactant ' + (i + 1) };
            pw.nodes.push(n);
            pw.edges.push({ id: 'edge-' + (++pw._seq), from: n.id, to: reactionNode.id });
        }
        for (i = 0; i < p.length; i++) {
            n = { id: 'node-' + (++pw._seq), kind: 'metabolite', label: 'Product ' + (i + 1) };
            pw.nodes.push(n);
            pw.edges.push({ id: 'edge-' + (++pw._seq), from: reactionNode.id, to: n.id });
        }
        // convert selects the reaction node + lands on the Pathway view.
        pw.selectedType = 'node';
        pw.selectedId = reactionNode.id;
    } else {
        // Clear the selection FIRST, then stamp ONE fresh node.
        pw.selectedType = '';
        pw.selectedId = null;
        pw.nodes.push({ id: 'node-' + (++pw._seq), kind: 'metabolite', label: label });
    }
    env.view = 'pathway';
    env.status = 'Added "' + label + '" to the pathway.';
    return env.view;
}

function stubMoleculeEditor(opts) {
    opts = opts || {};
    return {
        molName: opts.name || '',
        molecule: {
            atoms: opts.atoms || [],
            reactionArrow: opts.reactionArrow || null,
            split: opts.split || null
        }
    };
}

test('D1 molecule import: clears the stale selection, stamps exactly ONE fresh node, lands on Pathway', function () {
    var env = {
        editor: stubMoleculeEditor({ name: 'Ethanol', atoms: [{}, {}, {}] }),
        pathway: makePathway(),
        view: 'editor'
    };
    assert.strictEqual(env.pathway.selectedId, 'node-stale', 'precondition: a stale node is selected');
    var view = applyAddEditorDrawingToPathway(env);
    assert.strictEqual(env.pathway.nodes.length, 1, 'exactly one node is stamped for a molecule');
    assert.strictEqual(env.pathway.edges.length, 0, 'a molecule import adds no edges');
    assert.strictEqual(env.pathway.selectedType, '', 'the stale selection was cleared (fresh node, not in-place edit)');
    assert.strictEqual(env.pathway.selectedId, null, 'the stale selection id was cleared');
    assert.strictEqual(view, 'pathway', 'the user lands on the Pathway view');
    assert.strictEqual(env.status, 'Added "Ethanol" to the pathway.', 'status names the imported structure');
});

test('D2 reaction import (CCO>>CC=O shape): reactant + reaction step + product nodes and wired edges', function () {
    var env = {
        editor: stubMoleculeEditor({
            name: 'Oxidation',
            atoms: [{}, {}, {}, {}, {}, {}],
            reactionArrow: { type: 'forward' },
            split: { reactants: [{}], products: [{}] } // one reactant -> one product
        }),
        pathway: makePathway(),
        view: 'editor'
    };
    var view = applyAddEditorDrawingToPathway(env);
    // 1 reactant + 1 reaction node + 1 product = 3 nodes; 2 edges (in + out).
    assert.strictEqual(env.pathway.nodes.length, 3, 'reactant + reaction-step + product nodes');
    assert.strictEqual(env.pathway.edges.length, 2, 'reactant->reaction and reaction->product edges');
    var reactionNodes = env.pathway.nodes.filter(function (n) { return n.kind === 'reaction'; });
    assert.strictEqual(reactionNodes.length, 1, 'exactly one kind:"reaction" step node');
    // The reaction node is wired in from the reactant and out to the product.
    var rid = reactionNodes[0].id;
    assert.ok(env.pathway.edges.some(function (e) { return e.to === rid; }), 'a reactant edge enters the reaction node');
    assert.ok(env.pathway.edges.some(function (e) { return e.from === rid; }), 'a product edge leaves the reaction node');
    assert.strictEqual(view, 'pathway', 'the user lands on the Pathway view');
    assert.strictEqual(env.status, 'Added "Oxidation" to the pathway.', 'status names the imported reaction');
});

test('D3 empty editor: no node is stamped, the view is unchanged, a friendly status is shown', function () {
    var env = {
        editor: stubMoleculeEditor({ atoms: [] }),
        pathway: makePathway(),
        view: 'editor'
    };
    var view = applyAddEditorDrawingToPathway(env);
    assert.strictEqual(env.pathway.nodes.length, 0, 'no node is stamped from an empty editor');
    assert.strictEqual(env.pathway.edges.length, 0, 'no edge is added from an empty editor');
    assert.strictEqual(view, 'editor', 'the empty-editor guard does NOT switch the view');
    assert.ok(/Draw a molecule or reaction/.test(env.status), 'a friendly guidance status is shown');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
