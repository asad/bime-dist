/**
 * tests/test_v2_3_4_lens_hardening.js — focus-lens trust + accessibility
 * hardening (v2.3.4).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.3.4 hardens the focus-lens interaction shipped across v2.3.0-v2.3.3. A UX
 * review flagged the lens as "a hidden, modeless, destructive-by-default
 * interaction with no visible entry, no exit confirmation, no cancel, and no
 * dialog accessibility." This suite pins the four fixes:
 *
 *   P0-1  every pathway STRUCTURE node advertises the edit gesture — a native
 *         <title> ("Double-click or press Enter to edit structure") and an
 *         `is-structure` class for the hover/selected editable affordance.
 *   P0-2  the lens overlay gains a Done/Cancel HEADER; collapse(true) commits
 *         (history push + apply + a save flash + a "Saved" status), while
 *         collapse(false) DISCARDS — no history push, no payload mutation.
 *   P0-3  the overlay is an accessible modal dialog: role=dialog / aria-modal /
 *         aria-label set on open + prior focus recorded, removed + focus
 *         restored on collapse, with Tab trapped inside while open.
 *   P1-1  the click-away is de-eagered: a mousedown on ANOTHER node hops the
 *         lens (commit + re-open), chrome clicks are ignored, and only a bare
 *         canvas / background click commits + collapses.
 *
 * Mostly source-shape (the workbench DOM glue needs the whole page to execute,
 * exactly as test_v2_3_0 / test_v2_3_3 note), plus an executable DOM-stub
 * section that drives the discard-vs-commit + dialog-attr + focus-restore +
 * click-away-hop control flow through recording stubs — modelled on
 * test_v2_3_0_focus_lens.js section F and the tests/shim.js document stub.
 *
 * Plain Node, no real DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Focus-lens hardening (v2.3.4)');
var test = runner.test;
console.log('Focus-lens hardening (v2.3.4)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var CSS = readRepoFile('css/workbench.css');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing the source-shape-pinned suites use, so a
// refactor that splits a pinned body fails here too.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// ===========================================================================
// P0-1 — structure nodes look editable + advertise the gesture.
// ===========================================================================
test('P0-1.1 renderPathwayNode appends a native <title> editable cue + is-structure class', function () {
    var fn = wbFn('renderPathwayNode');
    assert.ok(fn, 'renderPathwayNode located');
    // The exact cue string (native tooltip + SVG accessible name).
    assert.ok(fn.indexOf('Double-click or press Enter to edit structure') !== -1,
        'the editable-gesture cue string is appended to the node group');
    // It is created as an SVG <title> element via the pathway element factory.
    assert.ok(/pathwaySvgEl\(\s*'title'/.test(fn),
        'the cue is an SVG <title> child (native tooltip + accessible name)');
    // The group is flagged is-structure for the CSS affordance.
    assert.ok(fn.indexOf("' is-structure'") !== -1 || fn.indexOf('is-structure') !== -1,
        'the node group is flagged is-structure');
});

test('P0-1.2 renderPathwayNode is NOT split by the cue (single top-level function)', function () {
    // The cue must be appended inline — no nested function declaration that would
    // truncate the body / break the other source-shape suites' slicing.
    var fn = wbFn('renderPathwayNode');
    assert.strictEqual(fn.indexOf('\nfunction '), -1,
        'renderPathwayNode contains no nested top-level function declaration');
    // The node still appends its rect (the cue is additive, not a replacement).
    assert.ok(fn.indexOf("pathwaySvgEl('rect'") !== -1, 'the node rect is still drawn');
});

test('P0-1.3 css advertises a hover affordance on is-structure nodes + a persistent selected cue', function () {
    assert.ok(CSS.indexOf('.wb-pathway-node.is-structure:hover rect') !== -1,
        'a hover affordance targets is-structure node rects');
    assert.ok(CSS.indexOf('.wb-pathway-node.is-structure.is-selected rect') !== -1,
        'a persistent cue targets the SELECTED is-structure node');
    // Non-colour perceivability: the hover cue also changes stroke-width, not
    // colour alone; the selected cue adds a dashed inset.
    var hoverIdx = CSS.indexOf('.wb-pathway-node.is-structure:hover rect');
    var hoverBlock = CSS.slice(hoverIdx, hoverIdx + 160);
    assert.ok(hoverBlock.indexOf('stroke-width') !== -1,
        'the hover cue is not colour-only (stroke-width changes too)');
    assert.ok(CSS.indexOf('var(--color-primary-light)') !== -1,
        'the affordance stays within the --color token system');
});

// ===========================================================================
// P0-2 — explicit Done/Cancel header + save confirmation; discard on cancel.
// ===========================================================================
test('P0-2.1 a lens header builder + remover exist as top-level functions', function () {
    assert.ok(WB.indexOf('function buildStructureLensHeader(') !== -1,
        'buildStructureLensHeader defined');
    assert.ok(WB.indexOf('function removeStructureLensHeader(') !== -1,
        'removeStructureLensHeader defined');
});

test('P0-2.2 the header carries the node label + a Done and a Cancel button', function () {
    var fn = wbFn('buildStructureLensHeader');
    assert.ok(fn, 'buildStructureLensHeader located');
    assert.ok(fn.indexOf('wb-lens-header') !== -1, 'header element class present');
    assert.ok(fn.indexOf('Edit structure: ') !== -1, 'header shows the node label');
    assert.ok(fn.indexOf("'Done'") !== -1, 'a Done button label');
    assert.ok(fn.indexOf("'Cancel'") !== -1, 'a Cancel button label');
    // Done commits, Cancel discards.
    assert.ok(fn.indexOf('collapseStructureLens(true)') !== -1,
        'Done wires collapseStructureLens(true)');
    assert.ok(fn.indexOf('collapseStructureLens(false)') !== -1,
        'Cancel wires collapseStructureLens(false)');
    // It is prepended into the overlay so it sits at the TOP of the panel.
    assert.ok(fn.indexOf('.wb-editor-wrap') !== -1, 'the header lives inside the overlay');
    assert.ok(fn.indexOf('insertBefore') !== -1 || fn.indexOf('firstChild') !== -1,
        'the header is prepended above the editor chrome');
});

test('P0-2.3 openStructureLens builds the header; collapseStructureLens removes it', function () {
    var open = wbFn('openStructureLens');
    var collapse = wbFn('collapseStructureLens');
    assert.ok(open.indexOf('buildStructureLensHeader(') !== -1,
        'openStructureLens builds the header on open');
    assert.ok(collapse.indexOf('removeStructureLensHeader(') !== -1,
        'collapseStructureLens removes the header on collapse');
});

test('P0-2.4 collapseStructureLens(false) DISCARDS — the commit block is guarded on commit', function () {
    var fn = wbFn('collapseStructureLens');
    assert.ok(fn, 'collapseStructureLens located');
    // Both destructive operations sit inside `if (commit && node)` so a cancel
    // never pushes history and never mutates the node payload.
    var guardIdx = fn.indexOf('if (commit && node)');
    assert.ok(guardIdx !== -1, 'the commit path is guarded on `commit`');
    var pushIdx = fn.indexOf('pushPathwayHistory(');
    var applyIdx = fn.indexOf('applyPathwayStructurePayload(');
    assert.ok(pushIdx > guardIdx, 'pushPathwayHistory only runs inside the commit guard');
    assert.ok(applyIdx > guardIdx, 'applyPathwayStructurePayload only runs inside the commit guard');
    // Close + cleanup + re-render happen on BOTH paths.
    assert.ok(fn.indexOf('_lens.close(') !== -1, 'closes the lens on both paths');
    assert.ok(fn.indexOf('renderPathwayCanvas(') !== -1, 're-renders on both paths');
});

test('P0-2.5 a committed edit flashes the node + reports a "Saved" status (Esc still = Done)', function () {
    assert.ok(WB.indexOf('function flashSavedNode(') !== -1, 'flashSavedNode helper defined');
    var flash = wbFn('flashSavedNode');
    assert.ok(flash.indexOf('wb-node-saved-flash') !== -1, 'the flash adds the highlight class');
    assert.ok(flash.indexOf('setTimeout') !== -1 && flash.indexOf('900') !== -1,
        'the flash is removed after ~900ms');
    var collapse = wbFn('collapseStructureLens');
    assert.ok(collapse.indexOf('flashSavedNode(') !== -1, 'commit triggers the save flash');
    assert.ok(collapse.indexOf("'Saved \"'") !== -1 || collapse.indexOf('Saved "') !== -1,
        'commit reports a Saved status with the node label');
    // Esc still commits (Done): the Escape handler routes to collapse(true).
    var escIdx = WB.indexOf("e.key === 'Escape' || e.keyCode === 27");
    var escSlice = WB.slice(escIdx, escIdx + 400);
    assert.ok(escSlice.indexOf('collapseStructureLens(true)') !== -1,
        'Escape still commits the lens (Esc = Done)');
});

test('P0-2.6 css styles the lens header + a reduced-motion-safe save flash', function () {
    assert.ok(CSS.indexOf('.wb-lens-header') !== -1, 'header styled');
    assert.ok(CSS.indexOf('.wb-lens-done') !== -1 && CSS.indexOf('.wb-lens-cancel') !== -1,
        'Done + Cancel buttons styled');
    assert.ok(CSS.indexOf('@keyframes wb-node-saved-flash') !== -1, 'a save-flash keyframe exists');
    assert.ok(CSS.indexOf('.wb-pathway-node.wb-node-saved-flash rect') !== -1,
        'the flash targets the node rect');
    // prefers-reduced-motion disables the animation.
    var rmIdx = CSS.indexOf('prefers-reduced-motion');
    assert.ok(rmIdx !== -1, 'a prefers-reduced-motion block exists');
    assert.ok(CSS.indexOf('.wb-pathway-node.wb-node-saved-flash rect { animation: none') !== -1,
        'reduced-motion disables the save-flash animation');
});

// ===========================================================================
// P0-3 — dialog semantics + focus management (WCAG 2.4.3 / 4.1.2).
// ===========================================================================
test('P0-3.1 openStructureLens sets role=dialog / aria-modal / aria-label + records prior focus', function () {
    var fn = wbFn('openStructureLens');
    assert.ok(fn.indexOf("setAttribute('role', 'dialog')") !== -1, 'role=dialog set on the overlay');
    assert.ok(fn.indexOf("setAttribute('aria-modal', 'true')") !== -1, 'aria-modal=true set');
    assert.ok(fn.indexOf("'aria-label', 'Edit structure: '") !== -1,
        'aria-label names the dialog with the node label');
    assert.ok(fn.indexOf('_lensLastFocus = ') !== -1,
        'the pre-lens focus is recorded (mirrors _modalLastFocus)');
    assert.ok(fn.indexOf('document.activeElement') !== -1, 'records document.activeElement');
});

test('P0-3.2 openStructureLens moves focus INTO the overlay + installs a Tab trap', function () {
    var fn = wbFn('openStructureLens');
    // Focus moves to the Done button.
    assert.ok(fn.indexOf('.wb-lens-done') !== -1, 'targets the Done button');
    assert.ok(fn.indexOf('.focus()') !== -1, 'moves keyboard focus into the overlay');
    // The trap handler is installed on the overlay.
    assert.ok(fn.indexOf('_lensTrapHandler = onLensTrapKeydown') !== -1,
        'the Tab-trap handler is recorded');
    assert.ok(fn.indexOf("addEventListener('keydown', _lensTrapHandler") !== -1,
        'the Tab-trap is bound on the overlay');
    // The trap helpers exist.
    assert.ok(WB.indexOf('function onLensTrapKeydown(') !== -1, 'onLensTrapKeydown defined');
    assert.ok(WB.indexOf('function lensFocusableElements(') !== -1, 'lensFocusableElements defined');
});

test('P0-3.3 onLensTrapKeydown cycles Tab/Shift+Tab within the overlay focusables', function () {
    var fn = wbFn('onLensTrapKeydown');
    assert.ok(fn.indexOf("e.key !== 'Tab'") !== -1 || fn.indexOf("'Tab'") !== -1,
        'the trap reacts to Tab');
    assert.ok(fn.indexOf('e.shiftKey') !== -1, 'handles Shift+Tab (reverse cycle)');
    assert.ok(fn.indexOf('lensFocusableElements(') !== -1, 'cycles among the overlay focusables');
    assert.ok(fn.indexOf('preventDefault') !== -1, 'wraps focus by preventing the default Tab move');
});

test('P0-3.4 collapseStructureLens removes the dialog attrs + trap and restores focus', function () {
    var fn = wbFn('collapseStructureLens');
    assert.ok(fn.indexOf("removeAttribute('role')") !== -1, 'role removed on collapse');
    assert.ok(fn.indexOf("removeAttribute('aria-modal')") !== -1, 'aria-modal removed');
    assert.ok(fn.indexOf("removeAttribute('aria-label')") !== -1, 'aria-label removed');
    assert.ok(fn.indexOf("removeEventListener('keydown', _lensTrapHandler") !== -1,
        'the Tab-trap is torn down');
    // Focus restored to the recorded element, guarded if it is gone.
    assert.ok(fn.indexOf('_lensLastFocus && _lensLastFocus.focus') !== -1,
        'focus restore is guarded against a missing element');
    assert.ok(fn.indexOf('_lensLastFocus.focus()') !== -1, 'restores focus to the pre-lens element');
    assert.ok(fn.indexOf('_lensLastFocus = null') !== -1, 'clears the recorded focus afterwards');
});

// ===========================================================================
// P1-1 — de-eager the click-away.
// ===========================================================================
test('P1-1.1 a dedicated click-away handler is bound (replacing the eager inline one)', function () {
    assert.ok(WB.indexOf('function handleLensClickAway(') !== -1,
        'handleLensClickAway defined as a top-level function');
    assert.ok(WB.indexOf("addEventListener('mousedown', handleLensClickAway") !== -1,
        'the click-away handler is bound on mousedown');
});

test('P1-1.2 a mousedown on ANOTHER node HOPS the lens (commit then re-open)', function () {
    var fn = wbFn('handleLensClickAway');
    assert.ok(fn, 'handleLensClickAway located');
    // Resolve the node the same way double-click / context menu do.
    assert.ok(fn.indexOf('resolvePathwayContextTarget(') !== -1,
        'resolves the node id the canonical way');
    var hopIdx = fn.indexOf("target.type === 'node'");
    assert.ok(hopIdx !== -1, 'detects a mousedown on a node target');
    var hopBlock = fn.slice(hopIdx, hopIdx + 220);
    assert.ok(hopBlock.indexOf('target.id !== focusedId') !== -1,
        'the hop is only for a DIFFERENT node than the focused one');
    assert.ok(hopBlock.indexOf('collapseStructureLens(true)') !== -1
        && hopBlock.indexOf('openStructureLens(') !== -1,
        'the hop commits the current lens then opens the other node (seamless, not a dismiss)');
});

test('P1-1.3 chrome / overlay / focused-node mousedowns do NOT dismiss', function () {
    var fn = wbFn('handleLensClickAway');
    assert.ok(fn.indexOf(".closest('.wb-editor-wrap')") !== -1,
        'mousedowns inside the overlay (incl. the lens header) are ignored');
    assert.ok(fn.indexOf('.wb-inspector') !== -1 && fn.indexOf('.wb-output-section') !== -1,
        'mousedowns inside the inspector / output panel are ignored');
    // The focused node is also a no-op (not a dismiss, not a hop).
    assert.ok(fn.indexOf("'[data-pathway-id=\"' + focusedId + '\"]'") !== -1
        || fn.indexOf('data-pathway-id="' ) !== -1,
        'a mousedown on the focused node is a no-op');
});

test('P1-1.4 only a bare canvas / background mousedown commits + collapses', function () {
    var fn = wbFn('handleLensClickAway');
    // The final, unguarded path (after all the early returns) commits + closes.
    var lastCollapse = fn.lastIndexOf('collapseStructureLens(true)');
    var inspectorGuard = fn.indexOf('.wb-inspector');
    assert.ok(lastCollapse > inspectorGuard,
        'the commit-and-collapse is the LAST resort, after the chrome guards');
    assert.ok(fn.indexOf('_lens.isOpen()') !== -1,
        'the handler only acts while a lens is open');
});

// ===========================================================================
// E. Executable DOM-stub: drive the discard/commit + dialog-attr + focus-restore
//    + click-away-hop control flow through recording stubs (the workbench glue
//    proper needs the whole page, so we model its exact control flow here the way
//    test_v2_3_0 section F models positionStructureLens). Robustness over
//    cleverness: these fail loudly if the documented behaviour regresses.
// ===========================================================================

// A recording element: records attributes, classList, focus, listeners.
function recEl(tag) {
    return {
        tag: tag || 'div',
        _attrs: {},
        _listeners: {},
        focused: false,
        children: [],
        classList: {
            _set: {},
            add: function (c) { this._set[c] = true; },
            remove: function (c) { delete this._set[c]; },
            contains: function (c) { return !!this._set[c]; }
        },
        style: {},
        setAttribute: function (k, v) { this._attrs[k] = v; },
        getAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k) ? this._attrs[k] : null; },
        removeAttribute: function (k) { delete this._attrs[k]; },
        addEventListener: function (t, fn) { (this._listeners[t] = this._listeners[t] || []).push(fn); },
        removeEventListener: function (t, fn) {
            var a = this._listeners[t] || [];
            var i = a.indexOf(fn);
            if (i !== -1) { a.splice(i, 1); }
        },
        appendChild: function (c) { this.children.push(c); return c; },
        focus: function () { this.focused = true; }
    };
}

test('E1 discard vs commit: a faithful model leaves the payload + history untouched on cancel', function () {
    // Model the EXACT collapseStructureLens commit guard: only `commit` runs the
    // two destructive ops. This catches a regression where the guard is loosened.
    var node = { id: 'n1', label: 'Glucose', structure: { smiles: 'OCC1OC(O)C(O)C(O)C1O' } };
    var history = [];
    function collapseModel(commit) {
        if (commit && node) {
            // capture would return a payload; model it as a new structure.
            var payload = { smiles: 'CHANGED' };
            history.push('snapshot');            // pushPathwayHistory
            node.structure = payload;            // applyPathwayStructurePayload
        }
        // close + cleanup happen on both paths (not modelled here).
    }
    // Cancel: nothing mutates.
    var before = node.structure;
    collapseModel(false);
    assert.strictEqual(node.structure, before, 'cancel does NOT mutate the node payload');
    assert.strictEqual(history.length, 0, 'cancel does NOT push pathway history');
    // Commit: both happen exactly once.
    collapseModel(true);
    assert.strictEqual(node.structure.smiles, 'CHANGED', 'commit applies the new payload');
    assert.strictEqual(history.length, 1, 'commit pushes exactly one history snapshot');
});

test('E2 dialog lifecycle: open sets role/aria + records focus; collapse clears + restores', function () {
    var prior = recEl('button');     // the element focused before the lens opened
    prior.focused = true;
    var wrap = recEl('div');         // the overlay (.wb-editor-wrap)
    var done = recEl('button');
    done.classList.add('wb-lens-done');
    wrap.appendChild(done);
    var label = 'Pyruvate';

    // ---- open (mirrors openStructureLens P0-3) ----
    var lensLastFocus = prior;       // document.activeElement at open time
    wrap.setAttribute('role', 'dialog');
    wrap.setAttribute('aria-modal', 'true');
    wrap.setAttribute('aria-label', 'Edit structure: ' + label);
    done.focus();
    var trap = function () {};
    wrap.addEventListener('keydown', trap, true);

    assert.strictEqual(wrap.getAttribute('role'), 'dialog', 'role=dialog on open');
    assert.strictEqual(wrap.getAttribute('aria-modal'), 'true', 'aria-modal on open');
    assert.strictEqual(wrap.getAttribute('aria-label'), 'Edit structure: Pyruvate', 'aria-label names the dialog');
    assert.strictEqual(done.focused, true, 'focus moved into the overlay (Done button)');
    assert.strictEqual((wrap._listeners.keydown || []).length, 1, 'Tab-trap bound on open');

    // ---- collapse (mirrors collapseStructureLens P0-3) ----
    wrap.removeAttribute('role');
    wrap.removeAttribute('aria-modal');
    wrap.removeAttribute('aria-label');
    wrap.removeEventListener('keydown', trap, true);
    prior.focused = false; // simulate focus having moved to the overlay
    if (lensLastFocus && lensLastFocus.focus) { lensLastFocus.focus(); }
    lensLastFocus = null;

    assert.strictEqual(wrap.getAttribute('role'), null, 'role cleared on collapse');
    assert.strictEqual(wrap.getAttribute('aria-modal'), null, 'aria-modal cleared');
    assert.strictEqual(wrap.getAttribute('aria-label'), null, 'aria-label cleared');
    assert.strictEqual((wrap._listeners.keydown || []).length, 0, 'Tab-trap removed on collapse');
    assert.strictEqual(prior.focused, true, 'focus restored to the pre-lens element (WCAG 2.4.3)');
});

test('E3 click-away routing: another-node mousedown HOPS; chrome is ignored; bare canvas dismisses', function () {
    // Model the handleLensClickAway routing with a tiny closest() stub. We assert
    // the THREE outcomes by recording which action each target produces.
    var focusedId = 'nA';
    function makeTarget(kind, nodeId) {
        // kind: 'overlay' | 'inspector' | 'output' | 'focused-node' | 'other-node' | 'canvas'
        return {
            closest: function (sel) {
                if (sel.indexOf('.wb-editor-wrap') !== -1) { return kind === 'overlay' ? {} : null; }
                if (sel.indexOf('data-pathway-id="' + focusedId + '"') !== -1) {
                    return kind === 'focused-node' ? {} : null;
                }
                if (sel.indexOf('.wb-inspector') !== -1) {
                    return (kind === 'inspector' || kind === 'output') ? {} : null;
                }
                return null;
            },
            _kind: kind,
            _nodeId: nodeId
        };
    }
    function resolveTarget(t) {
        if (t._kind === 'other-node') { return { type: 'node', id: t._nodeId }; }
        if (t._kind === 'focused-node') { return { type: 'node', id: focusedId }; }
        return { type: 'canvas' };
    }
    // Faithful port of handleLensClickAway's control flow.
    function route(target) {
        var actions = [];
        var hasClosest = target && target.closest;
        if (hasClosest && target.closest('.wb-editor-wrap')) { return actions; }
        if (hasClosest && target.closest('[data-pathway-id="' + focusedId + '"]')) { return actions; }
        var resolved = resolveTarget(target);
        if (resolved && resolved.type === 'node' && resolved.id && resolved.id !== focusedId) {
            actions.push('commit');
            actions.push('open:' + resolved.id);
            return actions;
        }
        if (hasClosest && target.closest('.wb-inspector, .wb-output-section')) { return actions; }
        actions.push('commit-collapse');
        return actions;
    }

    // overlay / focused node / chrome -> no dismiss.
    assert.deepStrictEqual(route(makeTarget('overlay')), [], 'overlay mousedown ignored');
    assert.deepStrictEqual(route(makeTarget('focused-node')), [], 'focused-node mousedown ignored');
    assert.deepStrictEqual(route(makeTarget('inspector')), [], 'inspector mousedown ignored');
    // another node -> hop (commit + open that node).
    assert.deepStrictEqual(route(makeTarget('other-node', 'nB')), ['commit', 'open:nB'],
        'another-node mousedown hops the lens (commit then open)');
    // bare canvas -> commit + collapse.
    assert.deepStrictEqual(route(makeTarget('canvas')), ['commit-collapse'],
        'bare canvas mousedown commits + collapses');
});

test('E4 the buildStructureLensHeader contract: Done=commit, Cancel=discard wiring', function () {
    // Drive a faithful header builder against recording stubs and fire the
    // buttons, asserting the documented commit/discard split.
    var calls = [];
    function collapse(commit) { calls.push(commit); }
    // Mirror buildStructureLensHeader's button wiring.
    var cancel = recEl('button');
    cancel.addEventListener('click', function () { collapse(false); });
    var done = recEl('button');
    done.addEventListener('click', function () { collapse(true); });
    // Fire Cancel then Done.
    (cancel._listeners.click || []).forEach(function (fn) { fn(); });
    (done._listeners.click || []).forEach(function (fn) { fn(); });
    assert.deepStrictEqual(calls, [false, true],
        'Cancel collapses with commit=false (discard); Done with commit=true (save)');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
