/**
 * tests/test_v2_3_3_continuum_polish.js — continuum POLISH + CLEANUP (v2.3.3, merge Stage 5).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Stage 5 is the FINAL polish of the editor/pathway focus-lens continuum (the
 * Molecule/Pathway mode toggle was retired in Stage 3; the focus lens shipped in
 * Stage 2; focus-scoped undo in Stage 4). This suite pins the Stage-5 source
 * shape WITHOUT a browser DOM, mirroring test_v2_3_1 / test_v2_3_2:
 *   (a) the hidden mode-bar buttons no longer carry stale `aria-pressed`, while
 *       the `set-workbench-mode` + `data-mode` wiring stays so the action and the
 *       compatibility shims remain resolvable;
 *   (b) onPathwayKeyDown has an Enter -> openStructureLens branch, guarded by a
 *       single NODE selection and a no-lens-open condition (the keyboard twin of
 *       double-click), placed AFTER the Ctrl+Z/Y lens-routing block;
 *   (c) positionStructureLens flips the overlay to the opposite side of the node
 *       on edge overflow (LEFT of the node on right overflow, ABOVE on bottom
 *       overflow) instead of dragging the panel back over itself; and
 *   (d) the stale "on the next mode toggle" comment is gone (no mode toggle
 *       exists post-Stage-3).
 *
 * Plain Node, no DOM (read files as text), mirroring test_v2_3_1 / test_v2_3_2.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Continuum polish + cleanup (v2.3.3)');
var test = runner.test;
console.log('Continuum polish + cleanup (v2.3.3)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var HTML = readRepoFile('workbench.html');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing test_v2_0_67 / test_v2_3_0 / test_v2_3_1 /
// test_v2_3_2 use, so a refactor that splits a pinned body fails here too.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// Slice the hidden mode-bar markup block (the `wb-modebar` <div> through its
// closing tag) so assertions are about THAT block, not unrelated markup.
function modebarBlock() {
    var start = HTML.indexOf('wb-modebar');
    if (start === -1) { return ''; }
    // Back up to the opening '<div' of the modebar container.
    var open = HTML.lastIndexOf('<div', start);
    var slice = HTML.slice(open === -1 ? start : open);
    var end = slice.indexOf('</div>');
    return end === -1 ? slice : slice.slice(0, end + '</div>'.length);
}

// ---------------------------------------------------------------------------
// A. The hidden mode-bar buttons dropped the stale aria-pressed, but keep the
//    set-workbench-mode + data-mode wiring (so the action + shims still resolve).
// ---------------------------------------------------------------------------
// v2.4.0 (Hybrid Stage A): the owner chose the HYBRID architecture, which
// SUPERSEDES the pure-continuum "hidden mode bar" invariant this assertion
// originally encoded. The modebar is now the VISIBLE "Editor | Pathway" view
// switch (a deliberate addition on top of the kept focus lens), so two parts of
// the v2.3.3 nit no longer hold and are intentionally updated here (the rest of
// the assertion — that the set-workbench-mode action stays wired — is unchanged):
//   1. the "Molecule" mode value is repurposed to the standalone "Editor" view
//      (data-mode="molecule" -> data-mode="editor"); data-mode="pathway" stays;
//   2. the buttons are role=tab and now legitimately carry the active-state
//      attribute (aria-selected, plus the aria-pressed synonym the dark-mode CSS
//      keys on) — that is the point of a visible tablist, not the stale leftover
//      the v2.3.3 nit removed from the then-hidden buttons.
// All OTHER v2.3.3 assertions (Enter-to-lens, overlay flip, no "mode toggle"
// prose) and every other suite stay green. See tests/test_v2_4_0_hybrid_views.js.
test('A1 the wb-modebar block keeps set-workbench-mode + the Editor/Pathway data-mode values (Hybrid view switch)', function () {
    var block = modebarBlock();
    assert.ok(block, 'wb-modebar block located');
    // The action wiring MUST stay (release_integrity + test_v2_3_1 + test_v2_4_0
    // pin it, and the pinned occurrence counts depend on it).
    assert.ok(block.indexOf('data-wb-action="set-workbench-mode"') !== -1,
        'set-workbench-mode action stays in the (now visible) view-switch bar');
    // v2.4.0: the two views are Editor + Pathway (the old "molecule" value is the
    // standalone Editor view now); data-mode="pathway" remains (pinned elsewhere).
    assert.ok(block.indexOf('data-mode="editor"') !== -1
        && block.indexOf('data-mode="pathway"') !== -1,
        'both view data-mode values (editor + pathway) are present in the view-switch bar');
});

test('A2 the pinned occurrence counts are intact (no markup churn from the nit)', function () {
    function count(hay, needle) {
        var n = 0, i = 0;
        while ((i = hay.indexOf(needle, i)) !== -1) { n++; i += needle.length; }
        return n;
    }
    assert.strictEqual(count(HTML, 'data-wb-action="set-workbench-mode"'), 2,
        'set-workbench-mode wired exactly twice (the two Editor | Pathway view tabs)');
    assert.strictEqual(count(HTML, 'data-mode="pathway"'), 1, 'data-mode="pathway" x1');
    assert.strictEqual(count(HTML, 'data-wb-action="add-pathway-current"'), 2, 'add-pathway-current x2');
    assert.strictEqual(count(HTML, 'data-wb-action="edit-pathway-structure"'), 2, 'edit-pathway-structure x2');
    assert.strictEqual(count(HTML, 'data-wb-action="reaction-to-pathway"'), 2, 'reaction-to-pathway x2');
});

// ---------------------------------------------------------------------------
// B. Enter-to-focus parity: a selected NODE + Enter opens the focus lens, but
//    ONLY when no lens is already open, and AFTER the Ctrl+Z/Y lens routing.
// ---------------------------------------------------------------------------
test('B1 onPathwayKeyDown has an Enter -> openStructureLens branch guarded by a node selection + no-lens condition', function () {
    var fn = wbFn('onPathwayKeyDown');
    assert.ok(fn, 'onPathwayKeyDown located');
    var enterIdx = fn.indexOf("e.key === 'Enter'");
    assert.ok(enterIdx !== -1, 'an Enter chord is handled');
    // Isolate the Enter branch so the guards/route are asserted about IT.
    var block = fn.slice(enterIdx, enterIdx + 320);
    assert.ok(block.indexOf('openStructureLens(') !== -1,
        'Enter opens the focus lens (the keyboard twin of double-click)');
    assert.ok(block.indexOf("_pathway.selectedType === 'node'") !== -1,
        'guarded on a single NODE selection');
    assert.ok(block.indexOf('_pathway.selectedId') !== -1,
        'guarded on an actual selected id');
    assert.ok(block.indexOf('_lens.isOpen()') !== -1,
        'guarded on the NO-lens-open condition (never re-opens while focused)');
    assert.ok(block.indexOf('preventDefault') !== -1,
        'the Enter branch preventDefaults so the key is consumed');
});

test('B2 the Enter branch sits AFTER the Ctrl+Z/Y lens-routing block (undo branches undisturbed)', function () {
    var fn = wbFn('onPathwayKeyDown');
    var enterIdx = fn.indexOf("e.key === 'Enter'");
    // The lens-routing Ctrl+Z/Y block (Stage 4) routes to editor._doAction and
    // must precede the new Enter branch; pathway undo/redo must still be present.
    var doActionIdx = fn.indexOf('editor._doAction(');
    var undoIdx = fn.indexOf('undoPathway()');
    var redoIdx = fn.indexOf('redoPathway()');
    assert.ok(doActionIdx !== -1 && undoIdx !== -1 && redoIdx !== -1,
        'the Ctrl+Z/Y lens route + pathway undo/redo branches are all still present');
    assert.ok(doActionIdx < enterIdx,
        'the lens-routing Ctrl+Z/Y block precedes the Enter branch');
    assert.ok(undoIdx < enterIdx && redoIdx < enterIdx,
        'the pathway undo/redo branches precede the Enter branch (left undisturbed)');
});

// ---------------------------------------------------------------------------
// C. positionStructureLens flips the overlay to the opposite side of the node on
//    edge overflow instead of dragging the panel back over the node.
// ---------------------------------------------------------------------------
test('C1 positionStructureLens flips LEFT-of-node on right overflow and ABOVE-node on bottom overflow', function () {
    var fn = wbFn('positionStructureLens');
    assert.ok(fn, 'positionStructureLens located');
    // The flip anchors relative to the node's own edge (anchor.left/anchor.top),
    // NOT the viewport edge — that is what keeps the panel off the node.
    assert.ok(fn.indexOf('anchor.left -') !== -1,
        'right-overflow flip anchors to the LEFT of the node (anchor.left - panel width)');
    assert.ok(fn.indexOf('anchor.top -') !== -1,
        'bottom-overflow flip anchors ABOVE the node (anchor.top - panel height)');
    // The flip is gated on the same overflow test that used to drag the panel
    // back over the node.
    assert.ok(fn.indexOf('vw - MARGIN') !== -1 && fn.indexOf('vh - MARGIN') !== -1,
        'the flip is gated on the right/bottom viewport-overflow tests');
    // The min-size + viewport clamp survives as the FINAL safety after the flip.
    var flipIdx = fn.indexOf('anchor.left -');
    var clampIdx = fn.indexOf('left = vw - MARGIN - w');
    assert.ok(clampIdx !== -1, 'the viewport clamp is still present as the final safety');
    assert.ok(flipIdx < clampIdx,
        'the flip is attempted BEFORE the final viewport clamp');
});

// ---------------------------------------------------------------------------
// D. The stale "mode toggle" comment is gone (no mode toggle exists post-Stage-3).
// ---------------------------------------------------------------------------
test('D1 no residual "mode toggle" comment remains in js/workbench.js', function () {
    assert.strictEqual(WB.indexOf('on the next mode toggle'), -1,
        'the stale "on the next mode toggle" comment was reworded (no toggle exists now)');
    assert.strictEqual(WB.indexOf('mode toggle'), -1,
        'no "mode toggle" language survives anywhere in workbench.js');
});

test('D2 collapseStructureLens still clears the inline overlay positioning (the comment\'s subject is intact)', function () {
    var fn = wbFn('collapseStructureLens');
    assert.ok(fn, 'collapseStructureLens located');
    // The reworded comment still documents the wrap-style clear it sits above.
    assert.ok(fn.indexOf("wrap.style.position = ''") !== -1,
        'the inline-position clear (the comment refers to it) is intact');
    assert.ok(fn.indexOf('once the lens collapses') !== -1,
        'the comment now references the lens collapse, not a mode toggle');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
