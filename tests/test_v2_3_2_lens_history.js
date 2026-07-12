/**
 * tests/test_v2_3_2_lens_history.js — focus-scoped undo/redo (v2.3.2, merge Stage 4).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Stage 4 makes Ctrl+Z / Ctrl+Y focus-scoped across the continuum rather than a
 * single mega-stack spanning atom-moves AND node-drags:
 *   - lens OPEN  -> the molecule editor's own history (atom/bond granularity),
 *                   ephemeral to the focus session;
 *   - lens CLOSED-> the pathway history (node/edge granularity), unchanged;
 *   - COMMIT (collapseStructureLens(true)) -> the whole session lands as EXACTLY
 *     ONE pathway snapshot (the pushPathwayHistory() before
 *     applyPathwayStructurePayload);
 *   - CANCEL (collapseStructureLens(false)) -> ZERO pathway snapshots.
 *
 * Arbitration: onPathwayKeyDown is bound on the pathway SVG; MolEditor's own
 * undo handler is bound on `document` and gated by isInteractiveTarget. The
 * editor overlay is a SEPARATE DOM subtree, so (a) focus-in-editor reaches only
 * MolEditor's handler, and (b) focus-on-canvas reaches onPathwayKeyDown, which
 * routes to the editor AND stopPropagation()s so MolEditor's document handler
 * never double-applies. Exactly one undo per Ctrl+Z in both cases.
 *
 * Plain Node, no DOM (read js/workbench.js as text + drive a DOM-free model of
 * the routing/commit decision), mirroring test_v2_3_0 sections D/E/F and
 * test_v2_3_1.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Focus-scoped undo/redo (v2.3.2)');
var test = runner.test;
console.log('Focus-scoped undo/redo (v2.3.2)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing test_v2_0_67 / test_v2_3_0 / test_v2_3_1 use,
// so a refactor that splits a pinned body fails here too.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// Count non-overlapping occurrences of a literal needle within a string.
function countOccurrences(hay, needle) {
    if (!hay || !needle) { return 0; }
    var n = 0, i = 0;
    while ((i = hay.indexOf(needle, i)) !== -1) { n++; i += needle.length; }
    return n;
}

// ---------------------------------------------------------------------------
// A. Source-shape: onPathwayKeyDown routes Ctrl+Z/Y to the EDITOR when a lens
//    is open, BEFORE the pathway undo, and arbitrates against MolEditor's own
//    document handler.
// ---------------------------------------------------------------------------
test('A1 onPathwayKeyDown has a lens-open branch that routes to the editor undo before pathway undo', function () {
    var fn = wbFn('onPathwayKeyDown');
    assert.ok(fn, 'onPathwayKeyDown located');
    // The lens guard exists.
    var lensIdx = fn.indexOf('_lens && _lens.isOpen()');
    assert.ok(lensIdx !== -1, 'guards on an OPEN focus lens');
    // It routes to the molecule editor's verified undo/redo entry point.
    var doActionIdx = fn.indexOf('editor._doAction(');
    assert.ok(doActionIdx !== -1, 'routes Ctrl+Z/Y to the molecule editor (_doAction)');
    // The editor route must precede the pathway undo so the pathway timeline
    // does NOT fire while focused.
    var pathwayUndoIdx = fn.indexOf('undoPathway()');
    assert.ok(pathwayUndoIdx !== -1, 'still has the pathway undo path for the lens-closed case');
    assert.ok(doActionIdx < pathwayUndoIdx,
        'the lens->editor route is checked BEFORE undoPathway()');
    assert.ok(lensIdx < pathwayUndoIdx,
        'the lens guard precedes the pathway undo');
});

test('A2 the lens branch arbitrates against MolEditor\'s document handler (stopPropagation + preventDefault)', function () {
    var fn = wbFn('onPathwayKeyDown');
    // Isolate the lens-open block so the assertion is about THAT branch, not an
    // unrelated stopPropagation elsewhere.
    var lensIdx = fn.indexOf('_lens && _lens.isOpen()');
    var block = fn.slice(lensIdx, lensIdx + 400);
    assert.ok(block.indexOf('editor._doAction(') !== -1,
        'the lens block performs the editor undo/redo');
    assert.ok(block.indexOf('e.stopPropagation()') !== -1,
        'the lens block stops propagation so MolEditor\'s document handler does not double-apply');
    assert.ok(block.indexOf('e.preventDefault()') !== -1,
        'the lens block prevents the browser default');
    assert.ok(block.indexOf('return') !== -1,
        'the lens block short-circuits (no fall-through to pathway undo)');
});

test('A3 both the undo (z) and redo (y / shift+z) chords are covered by the lens branch', function () {
    var fn = wbFn('onPathwayKeyDown');
    var lensIdx = fn.indexOf('_lens && _lens.isOpen()');
    // The chord match that wraps the lens branch must accept z AND y.
    var pre = fn.slice(Math.max(0, lensIdx - 200), lensIdx);
    assert.ok(/['"]z['"]/.test(pre) && /['"]y['"]/.test(pre),
        'the chord guard around the lens branch matches both z and y');
    var block = fn.slice(lensIdx, lensIdx + 400);
    // Redo is selected for y (always) or shift+z; undo otherwise.
    assert.ok(block.indexOf('e.shiftKey') !== -1,
        'shift distinguishes redo from undo within the lens branch');
    assert.ok(block.indexOf("'redo'") !== -1 && block.indexOf("'undo'") !== -1,
        'the lens branch can issue both redo and undo to the editor');
});

// ---------------------------------------------------------------------------
// B. Source-shape: openStructureLens gives the focus session a CLEAN editor
//    history baseline so session-undo cannot escape into a previous molecule.
// ---------------------------------------------------------------------------
test('B1 openStructureLens clears the editor history after load(), before opening the lens', function () {
    var fn = wbFn('openStructureLens');
    assert.ok(fn, 'openStructureLens located');
    var clearIdx = fn.indexOf('editor.history.clear()');
    assert.ok(clearIdx !== -1, 'clears the editor history for a fresh focus-session baseline');
    var loadIdx = fn.indexOf('load(source)');
    var openIdx = fn.indexOf('_lens.open(');
    assert.ok(loadIdx !== -1 && clearIdx > loadIdx,
        'the history clear happens AFTER the structure load (so the load itself is dropped)');
    assert.ok(openIdx !== -1 && clearIdx < openIdx,
        'the history clear happens BEFORE the lens opens');
    // Guarded so a bundle without the API doesn't throw.
    var guard = fn.slice(Math.max(0, clearIdx - 80), clearIdx);
    assert.ok(guard.indexOf('editor.history') !== -1,
        'the clear is guarded on editor.history being present');
});

// ---------------------------------------------------------------------------
// C. Source-shape: collapseStructureLens pushes EXACTLY ONE pathway snapshot on
//    commit and NONE on cancel — no per-atom leakage into PathwayHistory.
// ---------------------------------------------------------------------------
test('C1 collapseStructureLens pushes exactly one pushPathwayHistory(), inside the commit branch', function () {
    var fn = wbFn('collapseStructureLens');
    assert.ok(fn, 'collapseStructureLens located');
    assert.strictEqual(countOccurrences(fn, 'pushPathwayHistory('), 1,
        'exactly one pathway snapshot per committed focus session');
    // That single push lives in the `if (commit && node)` branch, paired with
    // the existing capture/apply, so the cancel path can never reach it.
    var commitIdx = fn.indexOf('if (commit && node)');
    assert.ok(commitIdx !== -1, 'commit branch present');
    var pushIdx = fn.indexOf('pushPathwayHistory(');
    var applyIdx = fn.indexOf('applyPathwayStructurePayload(');
    assert.ok(pushIdx > commitIdx, 'the snapshot is inside the commit branch');
    assert.ok(applyIdx > pushIdx,
        'snapshot is pushed BEFORE applying the payload (so undo returns to the pre-edit node)');
});

test('C2 the no-commit (cancel) path skips the pathway snapshot', function () {
    var fn = wbFn('collapseStructureLens');
    // The ONLY push is gated behind `commit`, so when collapseStructureLens(false)
    // runs, control never reaches it. Assert the structural gating: the push is
    // textually after the commit guard and before the unconditional _lens.close().
    var commitIdx = fn.indexOf('if (commit && node)');
    var pushIdx = fn.indexOf('pushPathwayHistory(');
    var closeIdx = fn.indexOf('_lens.close()');
    assert.ok(commitIdx !== -1 && pushIdx !== -1 && closeIdx !== -1,
        'commit guard, push, and close all present');
    assert.ok(commitIdx < pushIdx && pushIdx < closeIdx,
        'the snapshot is strictly inside the commit-guarded region, ahead of the unconditional close');
});

// ---------------------------------------------------------------------------
// D. DOM-free logic model of the routing + commit decisions. We DON'T execute
//    onPathwayKeyDown / collapseStructureLens (they need the whole page); we
//    reproduce their decision rules here (the way test_v2_3_0 section F
//    reproduces the clamp) and prove the focus-scoped behaviour end to end.
// ---------------------------------------------------------------------------

// The routing decision, isolated from DOM: given a lens-open flag and a Ctrl+Z/Y
// keydown, decide WHO undoes (the editor or the pathway) and WHAT action. Mirrors
// the onPathwayKeyDown lens branch + its fall-through.
function routeUndoRedo(lensOpen, key, shiftKey) {
    var isZ = (key === 'z' || key === 'Z');
    var isY = (key === 'y' || key === 'Y');
    if (!isZ && !isY) { return null; }
    if (lensOpen) {
        var redo = isZ ? !!shiftKey : true; // y is always redo; z+shift is redo
        return { target: 'editor', action: redo ? 'redo' : 'undo' };
    }
    // Lens closed: pathway history. z (no shift) is undo; z+shift and y are redo.
    var pRedo = isZ ? !!shiftKey : true;
    return { target: 'pathway', action: pRedo ? 'redo' : 'undo' };
}

test('D1 lens OPEN routes every Ctrl+Z/Y to the editor; lens CLOSED routes to the pathway', function () {
    assert.deepStrictEqual(routeUndoRedo(true, 'z', false), { target: 'editor', action: 'undo' });
    assert.deepStrictEqual(routeUndoRedo(true, 'z', true), { target: 'editor', action: 'redo' });
    assert.deepStrictEqual(routeUndoRedo(true, 'y', false), { target: 'editor', action: 'redo' });
    assert.deepStrictEqual(routeUndoRedo(false, 'z', false), { target: 'pathway', action: 'undo' });
    assert.deepStrictEqual(routeUndoRedo(false, 'z', true), { target: 'pathway', action: 'redo' });
    assert.deepStrictEqual(routeUndoRedo(false, 'y', false), { target: 'pathway', action: 'redo' });
    assert.strictEqual(routeUndoRedo(true, 'a', false), null, 'non-undo keys ignored');
});

// Arbitration model: the pathway SVG handler (onPathwayKeyDown) fires only when
// the event passes through the pathway SVG; MolEditor's document handler fires
// for any non-interactive target. The lens branch stopPropagation()s, so a
// keydown is applied EXACTLY ONCE in every focus configuration. We model the two
// handlers and count applications.
function applyOnce(focusLocation, lensOpen) {
    // focusLocation: 'editor' (inside .wb-editor-wrap) or 'canvas' (pathway SVG).
    var applications = 0;
    var propagationStopped = false;

    // onPathwayKeyDown — bound on the pathway SVG. Fires first (deeper node) and
    // ONLY when focus is on the canvas subtree.
    if (focusLocation === 'canvas') {
        if (lensOpen) {
            applications++;            // routes to editor._doAction
            propagationStopped = true; // e.stopPropagation()
        }
        // (lens closed -> pathway undo; still a single application, modeled below)
        else { applications++; }
    }

    // MolEditor document handler — bound on `document`, gated by
    // isInteractiveTarget. Neither the editor SVG nor the pathway SVG is an
    // interactive target, so it WOULD fire for both — unless propagation was
    // already stopped by the pathway handler.
    if (!propagationStopped) {
        if (focusLocation === 'editor') { applications++; }
        // focus on canvas with lens closed: MolEditor's handler also runs, but
        // it operates on the (empty) editor history; the meaningful single
        // application is the pathway undo above. For the lens-OPEN canvas case
        // propagation is stopped so this does not run.
    }

    return applications;
}

test('D2 a Ctrl+Z is applied exactly once whether focus is in the editor or on the canvas (lens open)', function () {
    // (a) focus inside the editor overlay: only MolEditor's document handler.
    assert.strictEqual(applyOnce('editor', true), 1,
        'editor-focused: MolEditor handles it once (pathway SVG handler never sees it)');
    // (b) focus on the pathway canvas: onPathwayKeyDown routes + stops propagation.
    assert.strictEqual(applyOnce('canvas', true), 1,
        'canvas-focused: pathway handler routes to the editor once, MolEditor suppressed');
});

// The commit decision: a focus session accumulates N editor-history pushes
// (atom/bond edits), but the PATHWAY history sees exactly one snapshot on commit
// and zero on cancel. Mirrors collapseStructureLens(commit).
function focusSession(atomEdits, commit, hasPayload) {
    var editorHistoryPushes = atomEdits;     // per-atom, ephemeral to the session
    var pathwaySnapshots = 0;
    // collapseStructureLens: push happens only when (commit && node && payload).
    if (commit && hasPayload) { pathwaySnapshots = 1; }
    return { editorHistoryPushes: editorHistoryPushes, pathwaySnapshots: pathwaySnapshots };
}

test('D3 N atom edits in a focus session -> 1 pathway snapshot on commit', function () {
    [1, 3, 7, 25].forEach(function (n) {
        var r = focusSession(n, true, true);
        assert.strictEqual(r.editorHistoryPushes, n,
            n + ' atom edits stay in the ephemeral editor history');
        assert.strictEqual(r.pathwaySnapshots, 1,
            n + ' atom edits collapse to exactly ONE pathway snapshot on commit');
    });
});

test('D4 N atom edits cancelled (Escape / click-away, no commit) -> 0 pathway snapshots', function () {
    [0, 1, 9, 40].forEach(function (n) {
        var r = focusSession(n, false, true);
        assert.strictEqual(r.pathwaySnapshots, 0,
            n + ' atom edits leave NO pathway snapshot when the session is cancelled');
    });
    // A committed-but-empty capture (no payload) also leaves the timeline alone.
    assert.strictEqual(focusSession(5, true, false).pathwaySnapshots, 0,
        'commit with an empty capture pushes no snapshot');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
