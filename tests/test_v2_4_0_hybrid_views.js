/**
 * tests/test_v2_4_0_hybrid_views.js — Hybrid Stage A: the Editor | Pathway view
 * switch (v2.4.0).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.3.x is the focus-lens CONTINUUM: no mode toggle; the pathway canvas is the
 * single permanent surface; double-clicking a structure node floats the molecule
 * editor as an in-place focus-lens overlay. v2.4.0 ADDS a dedicated standalone
 * Editor view + an explicit "Editor | Pathway" switch, while KEEPING the in-place
 * focus lens. mol+reaction stay one editor; pathway is a separate view. This is
 * the view-switch SKELETON only (no import-to-pathway, no AAM — later stages).
 *
 * The switch is SHOW/HIDE, never a reparent: Renderer.js positions atoms per
 * element (there is no single viewport <g>), so #bime-editor must never be
 * remounted. setWorkbenchMode stamps `.wb-view-editor` / `.wb-view-pathway` on
 * `.wb-main`; CSS owns the visibility. CRITICAL: the focus lens stays orthogonal
 * — its `.wb-main.wb-lens-open .wb-editor-wrap` overlay rule is `!important` and
 * keyed on `.wb-lens-open`, so the lens still opens over a node in the Pathway
 * view. This suite pins all of that without a browser DOM:
 *   A. workbench.html — the modebar is a VISIBLE tablist with Editor/Pathway
 *      data-mode buttons that keep data-wb-action="set-workbench-mode".
 *   B. js/workbench.js — `_wbView` defaults to 'pathway'; setWorkbenchMode still
 *      carries the pinned `openWorkbenchSection('wb-pathway-section')` literal,
 *      now also stamps the `.wb-view-*` classes, never hides the section via the
 *      v2.3.1-A3-forbidden `setWorkbenchSectionOpen('wb-pathway-section', false)`,
 *      and refits the editor when entering the Editor view.
 *   C. css/workbench.css — `.wb-view-editor .wb-editor-wrap{display:flex}` +
 *      `.wb-view-editor #wb-pathway-section{display:none}`; the `.wb-lens-open`
 *      `!important` overlay survives (lens protected); the base
 *      `.wb-editor-wrap{display:none}` is kept; the modebar is visible with a
 *      non-colour active-tab cue.
 *   D. an executable DOM-stub flips set-workbench-mode through a faithful port of
 *      setWorkbenchMode's control flow and asserts the `.wb-main` view class +
 *      the pathway-section visibility toggle, AND that `.wb-lens-open` coexists
 *      with `.wb-view-pathway` (the lens survives in the Pathway view).
 *
 * Source-shape + DOM-stub, plain Node (no real DOM beyond the shim) — mirroring
 * test_v2_3_1_continuum (A–E) + test_v2_3_4 section E.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Hybrid views: Editor | Pathway switch (v2.4.0)');
var test = runner.test;
console.log('Hybrid views: Editor | Pathway switch (v2.4.0)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var CSS = readRepoFile('css/workbench.css');
var HTML = readRepoFile('workbench.html');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing test_v2_3_1 / test_v2_0_67 use, so a refactor
// that splits a pinned body fails here too.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

// Slice the .wb-modebar markup block (opening tag through its closing </div>).
function modebarBlock() {
    var i = HTML.indexOf('class="wb-modebar"');
    if (i === -1) { return ''; }
    var open = HTML.lastIndexOf('<div', i);
    var close = HTML.indexOf('</div>', i);
    return open === -1 || close === -1 ? '' : HTML.slice(open, close + 6);
}

// ---------------------------------------------------------------------------
// A. workbench.html — the modebar is a VISIBLE Editor | Pathway tablist.
// ---------------------------------------------------------------------------
test('A1 the modebar is a real tablist (role=tablist) and no longer hidden', function () {
    var bar = modebarBlock();
    assert.ok(bar, '.wb-modebar block located in workbench.html');
    assert.ok(bar.indexOf('role="tablist"') !== -1,
        'the modebar is exposed as a tablist (it is the view switch now)');
    // v2.3.1 hid it with `hidden aria-hidden="true"`; the Hybrid reveals it. The
    // opening <div ...> tag must carry neither the `hidden` attr nor be
    // aria-hidden (the label span keeps its own aria-hidden — scope to the tag).
    var openTag = bar.slice(0, bar.indexOf('>') + 1);
    assert.strictEqual(/\shidden(\s|>|=)/.test(openTag), false,
        'the modebar is no longer hidden (the Editor | Pathway switch is visible)');
    assert.strictEqual(openTag.indexOf('aria-hidden="true"'), -1,
        'the modebar container is not aria-hidden');
    assert.strictEqual(bar.indexOf('wb-continuum-chrome'), -1,
        'the v2.3.1 continuum-chrome hide marker is gone');
});

test('A2 both tabs keep data-wb-action="set-workbench-mode" with role=tab + tabindex 0', function () {
    var bar = modebarBlock();
    var tabs = bar.match(/<button[^>]*class="wb-mode-btn"[^>]*>/g) || [];
    assert.strictEqual(tabs.length, 2, 'exactly two view tabs (Editor + Pathway)');
    for (var i = 0; i < tabs.length; i++) {
        assert.ok(tabs[i].indexOf('data-wb-action="set-workbench-mode"') !== -1,
            'tab ' + i + ' keeps the pinned set-workbench-mode action (dispatcher + compat shims)');
        assert.ok(tabs[i].indexOf('role="tab"') !== -1, 'tab ' + i + ' has role=tab');
        assert.ok(tabs[i].indexOf('tabindex="0"') !== -1, 'tab ' + i + ' is keyboard-reachable (tabindex 0)');
        assert.ok(tabs[i].indexOf('aria-selected') !== -1, 'tab ' + i + ' carries aria-selected');
    }
});

test('A3 the tabs are labelled Editor (data-mode=editor) and Pathway (data-mode=pathway)', function () {
    var bar = modebarBlock();
    assert.ok(/data-mode="editor"[^>]*>\s*Editor\s*<\/button>/.test(bar)
        || (bar.indexOf('data-mode="editor"') !== -1 && bar.indexOf('>Editor<') !== -1),
        'an "Editor" tab carries data-mode="editor"');
    assert.ok(/data-mode="pathway"[^>]*>\s*Pathway\s*<\/button>/.test(bar)
        || (bar.indexOf('data-mode="pathway"') !== -1 && bar.indexOf('>Pathway<') !== -1),
        'a "Pathway" tab carries data-mode="pathway"');
    // data-mode="pathway" must remain in the page (release_integrity + v2.3.1 A2
    // both pin it); the old data-mode="molecule" button is repurposed to editor.
    assert.ok(HTML.indexOf('data-mode="pathway"') !== -1,
        'data-mode="pathway" remains in workbench.html (pinned by other suites)');
});

test('A4 the dispatcher still routes set-workbench-mode through setWorkbenchMode', function () {
    assert.ok(WB.indexOf("action === 'set-workbench-mode'") !== -1
        && WB.indexOf("setWorkbenchMode(target.getAttribute('data-mode'))") !== -1,
        'set-workbench-mode stays resolvable via the click dispatcher');
});

// ---------------------------------------------------------------------------
// B. js/workbench.js — _wbView + grown setWorkbenchMode + tab sync.
// ---------------------------------------------------------------------------
test('B1 the default boot view _wbView is gated by PATHWAY_ENABLED (Pathway when enabled, Editor when hidden)', function () {
    // v2.4.18: the default boot view follows the pathway feature flag — 'pathway'
    // when the pathway editor ships, 'editor' for the core-only public build.
    assert.ok(/var\s+PATHWAY_ENABLED\s*=\s*(?:true|false)\s*;/.test(WB),
        'PATHWAY_ENABLED feature flag declared as a module global');
    assert.ok(/var\s+_wbView\s*=\s*PATHWAY_ENABLED\s*\?\s*'pathway'\s*:\s*'editor'\s*;/.test(WB),
        "_wbView = PATHWAY_ENABLED ? 'pathway' : 'editor' (flag-gated default boot view)");
});

test('B2 setWorkbenchMode keeps the pinned openWorkbenchSection literal AND grows the view switch', function () {
    var fn = wbFn('setWorkbenchMode');
    assert.ok(fn, 'setWorkbenchMode located');
    // Pinned by v2.3.1 A3 + a source-shape grep: the canvas stays mounted.
    assert.ok(fn.indexOf("openWorkbenchSection('wb-pathway-section')") !== -1,
        'keeps the pinned openWorkbenchSection(\'wb-pathway-section\') literal');
    // New: it now picks the standalone view by stamping the .wb-view-* classes.
    assert.ok(fn.indexOf("_wbView = (mode === 'editor') ? 'editor' : 'pathway'") !== -1,
        "sets _wbView from the mode ('editor' vs everything-else -> pathway)");
    assert.ok(fn.indexOf("'wb-view-editor'") !== -1 && fn.indexOf("'wb-view-pathway'") !== -1,
        'stamps both .wb-view-editor and .wb-view-pathway on .wb-main');
    assert.ok(fn.indexOf('classList.toggle') !== -1,
        'uses classList.toggle to add one view class and remove the other');
    // Entering the Editor view re-fits the just-shown editor.
    assert.ok(fn.indexOf('refitMoleculeEditorSoon()') !== -1,
        'refits the molecule editor when it becomes the hero (Editor view)');
});

test('B3 setWorkbenchMode does NOT hide the pathway section via the v2.3.1-A3-forbidden literal', function () {
    var fn = wbFn('setWorkbenchMode');
    assert.strictEqual(fn.indexOf("setWorkbenchSectionOpen('wb-pathway-section', false)"), -1,
        'the pathway section is hidden by the .wb-view-* CSS class, NOT by setWorkbenchSectionOpen(...,false)');
});

test('B4 syncWorkbenchModeButtons reflects the active view via aria-selected (role=tab)', function () {
    var fn = wbFn('syncWorkbenchModeButtons');
    assert.ok(fn, 'syncWorkbenchModeButtons located');
    assert.ok(fn.indexOf("setAttribute('aria-selected'") !== -1,
        'sets aria-selected on the tabs (the primary state for role=tab)');
    assert.ok(fn.indexOf("getAttribute('data-mode')") !== -1,
        'keys the selected tab off its data-mode (the view name)');
});

test('B5 syncWorkbenchLayoutMode defers the tab state to _wbView (not the section-open heuristic)', function () {
    var fn = wbFn('syncWorkbenchLayoutMode');
    assert.ok(fn, 'syncWorkbenchLayoutMode located');
    assert.ok(fn.indexOf('syncWorkbenchModeButtons(_wbView)') !== -1,
        'reflects the current _wbView so opening other drawers cannot steal the Editor tab');
});

test('B6 the boot path stamps the default Pathway view class on .wb-main', function () {
    // The default boot view is Pathway: .wb-main gets .wb-view-pathway from first
    // paint so the editor stays out of flow (base display:none) on load.
    assert.ok(WB.indexOf("classList.add('wb-view-pathway')") !== -1,
        'boot stamps .wb-view-pathway on .wb-main (default view)');
    assert.ok(WB.indexOf("classList.remove('wb-view-editor')") !== -1,
        'boot clears any .wb-view-editor');
});

// ---------------------------------------------------------------------------
// C. css/workbench.css — the Editor view rules + the protected lens overlay.
// ---------------------------------------------------------------------------
test('C1 the Editor view shows the editor as the in-flow hero (display:flex)', function () {
    assert.ok(/\.wb-main\.wb-view-editor\s+\.wb-editor-wrap\s*\{[^}]*display:\s*flex/.test(CSS),
        '.wb-main.wb-view-editor .wb-editor-wrap { display: flex } makes the editor the hero');
});

test('C2 the Editor view hides the pathway canvas section', function () {
    assert.ok(/\.wb-main\.wb-view-editor\s+#wb-pathway-section\s*\{[^}]*display:\s*none/.test(CSS),
        '.wb-main.wb-view-editor #wb-pathway-section { display: none }');
});

test('C3 CRITICAL: the focus-lens overlay rule still WINS (!important, keyed on .wb-lens-open)', function () {
    // The lens overlay is the ONLY editor-visible state in the Pathway view, and
    // it must beat both the base display:none AND any view rule. It is keyed on
    // .wb-lens-open (orthogonal to .wb-view-*) and uses !important.
    // Anchor on the RULE (selector + opening brace), not the prose that mentions
    // the same selector in the comments above it.
    var i = CSS.indexOf('.wb-main.wb-lens-open .wb-editor-wrap {');
    assert.ok(i !== -1, 'the .wb-main.wb-lens-open .wb-editor-wrap overlay rule exists');
    var block = CSS.slice(i, CSS.indexOf('}', i) + 1);
    assert.ok(/display:\s*flex\s*!important/.test(block),
        'the lens overlay forces display:flex !important so it beats the view rule -> the lens survives in Pathway view');
});

test('C4 the base editor chrome is display:none (continuum + lens contract, v2.3.1 E1)', function () {
    assert.ok(/\.wb-editor-wrap\s*\{[^}]*display:\s*none/.test(CSS),
        'base .wb-editor-wrap { display: none } kept (editor out of flow unless Editor view or lens)');
});

test('C5 the modebar is a VISIBLE tab bar (display:flex), no longer display:none', function () {
    var i = CSS.indexOf('.wb-modebar {');
    assert.ok(i !== -1, 'the .wb-modebar base rule exists');
    var block = CSS.slice(i, CSS.indexOf('}', i) + 1);
    assert.ok(/display:\s*flex/.test(block), '.wb-modebar is display:flex (the view switch is visible)');
    assert.strictEqual(/display:\s*none/.test(block), false,
        '.wb-modebar is no longer display:none (v2.3.1 hid it; the Hybrid reveals it)');
});

test('C6 the active tab carries a NON-COLOUR cue (weight + underline), not colour alone (v2.3.7)', function () {
    var i = CSS.indexOf('.wb-mode-btn[aria-selected="true"]');
    assert.ok(i !== -1, 'an active-tab rule keyed on aria-selected exists');
    var block = CSS.slice(i, CSS.indexOf('}', i) + 1);
    assert.ok(/font-weight:\s*[6-9]\d\d/.test(block),
        'active tab is bolder (a weight cue readable in greyscale)');
    assert.ok(/border-bottom-color/.test(block) || /border-bottom:/.test(block),
        'active tab gets an underline (a shape cue readable without colour)');
});

// ---------------------------------------------------------------------------
// D. Executable DOM-stub: drive the set-workbench-mode flip through a faithful
//    port of setWorkbenchMode's control flow (the workbench glue proper needs
//    the whole page, so we model its exact control flow here the way
//    test_v2_3_4 section E / test_v2_3_0 section F do). Fails loudly on regress.
// ---------------------------------------------------------------------------

// A recording element with a working classList (records add/remove/contains).
function recEl() {
    return {
        _attrs: {},
        classList: {
            _set: {},
            add: function (c) { this._set[c] = true; },
            remove: function (c) { delete this._set[c]; },
            toggle: function (c, on) { if (on) { this._set[c] = true; } else { delete this._set[c]; } },
            contains: function (c) { return !!this._set[c]; }
        },
        setAttribute: function (k, v) { this._attrs[k] = v; },
        getAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k) ? this._attrs[k] : null; }
    };
}

// Faithful port of setWorkbenchMode's view-switching control flow. `wbView` is
// the module-global stand-in; we return the new value so the test can thread it.
function applySetWorkbenchMode(mode, main, tabs, calls) {
    calls = calls || {};
    calls.openedSection = (calls.openedSection || 0) + 1; // openWorkbenchSection('wb-pathway-section')
    var wbView = (mode === 'editor') ? 'editor' : 'pathway';
    if (main) {
        main.classList.toggle('wb-view-editor', wbView === 'editor');
        main.classList.toggle('wb-view-pathway', wbView !== 'editor');
    }
    // syncWorkbenchModeButtons(wbView)
    for (var i = 0; i < tabs.length; i++) {
        var on = tabs[i].getAttribute('data-mode') === wbView ? 'true' : 'false';
        tabs[i].setAttribute('aria-selected', on);
        tabs[i].setAttribute('aria-pressed', on);
    }
    if (wbView === 'editor') { calls.refit = (calls.refit || 0) + 1; }
    return wbView;
}

function makeTabs() {
    var ed = recEl(); ed.setAttribute('data-mode', 'editor');
    var pw = recEl(); pw.setAttribute('data-mode', 'pathway');
    return { editor: ed, pathway: pw, all: [ed, pw] };
}

test('D1 flipping to Editor stamps .wb-view-editor, hides pathway, selects the Editor tab, refits', function () {
    var main = recEl();
    var tabs = makeTabs();
    var calls = {};
    var view = applySetWorkbenchMode('editor', main, tabs.all, calls);
    assert.strictEqual(view, 'editor', 'mode "editor" -> Editor view');
    assert.ok(main.classList.contains('wb-view-editor'), '.wb-main gains .wb-view-editor');
    assert.ok(!main.classList.contains('wb-view-pathway'), '.wb-main loses .wb-view-pathway');
    assert.strictEqual(tabs.editor.getAttribute('aria-selected'), 'true', 'Editor tab aria-selected=true');
    assert.strictEqual(tabs.pathway.getAttribute('aria-selected'), 'false', 'Pathway tab aria-selected=false');
    assert.strictEqual(calls.openedSection, 1, 'the pathway section stays MOUNTED (openWorkbenchSection ran)');
    assert.strictEqual(calls.refit, 1, 'the hero editor is re-fitted on entry');
});

test('D2 flipping back to Pathway restores .wb-view-pathway, re-shows the canvas, selects the Pathway tab', function () {
    var main = recEl();
    var tabs = makeTabs();
    applySetWorkbenchMode('editor', main, tabs.all, {});   // start in Editor
    var calls = {};
    var view = applySetWorkbenchMode('pathway', main, tabs.all, calls);
    assert.strictEqual(view, 'pathway', 'mode "pathway" -> Pathway view');
    assert.ok(main.classList.contains('wb-view-pathway'), '.wb-main gains .wb-view-pathway');
    assert.ok(!main.classList.contains('wb-view-editor'), '.wb-main loses .wb-view-editor (canvas re-shown)');
    assert.strictEqual(tabs.pathway.getAttribute('aria-selected'), 'true', 'Pathway tab aria-selected=true');
    assert.strictEqual(tabs.editor.getAttribute('aria-selected'), 'false', 'Editor tab aria-selected=false');
    assert.ok(!calls.refit, 'no hero re-fit when no lens is open in the Pathway view');
});

test('D3 legacy callers (e.g. setWorkbenchMode(\'molecule\')) land on the default Pathway view', function () {
    var main = recEl();
    var tabs = makeTabs();
    var view = applySetWorkbenchMode('molecule', main, tabs.all, {});
    assert.strictEqual(view, 'pathway', "any non-'editor' mode (incl. the legacy 'molecule') -> Pathway view");
    assert.ok(main.classList.contains('wb-view-pathway'), 'legacy caller leaves .wb-main on the Pathway view');
});

test('D4 CRITICAL: .wb-lens-open coexists with .wb-view-pathway (the focus lens survives in Pathway view)', function () {
    // In the Pathway view a double-click opens the focus lens: .wb-main carries
    // BOTH .wb-view-pathway AND .wb-lens-open. The classes are independent, and
    // the lens overlay rule (!important, keyed on .wb-lens-open) beats the base
    // display:none — so the editor paints as the floating lens over the node.
    var main = recEl();
    var tabs = makeTabs();
    applySetWorkbenchMode('pathway', main, tabs.all, {});
    // opening the lens only adds .wb-lens-open (it never touches the view class).
    main.classList.add('wb-lens-open');
    assert.ok(main.classList.contains('wb-view-pathway'), 'still in the Pathway view while the lens is open');
    assert.ok(main.classList.contains('wb-lens-open'), '.wb-lens-open is set (lens open over a node)');
    // The CSS contract that makes this paint correctly: the !important overlay
    // exists and the view rule has NO !important, so the overlay wins.
    var lensIdx = CSS.indexOf('.wb-main.wb-lens-open .wb-editor-wrap {');
    var lensBlock = CSS.slice(lensIdx, CSS.indexOf('}', lensIdx) + 1);
    assert.ok(/display:\s*flex\s*!important/.test(lensBlock),
        'the lens overlay forces display:flex !important -> wins over the view rule + base display:none');
    var viewIdx = CSS.indexOf('.wb-main.wb-view-editor .wb-editor-wrap {');
    var viewBlock = CSS.slice(viewIdx, CSS.indexOf('}', viewIdx) + 1);
    assert.strictEqual(/!important/.test(viewBlock), false,
        'the Editor-view rule uses NO !important, so it never out-ranks the lens overlay');
});

test('D5 the lens does not reparent the editor: the switch only ever flips classes on .wb-main', function () {
    // Regression guard for the architecture invariant: Renderer.js positions
    // atoms per element, so #bime-editor must not be remounted. The view switch
    // is a pure class flip — assert setWorkbenchMode touches no appendChild /
    // removeChild / insertBefore (no reparent of the editor host).
    var fn = wbFn('setWorkbenchMode');
    assert.strictEqual(fn.indexOf('appendChild'), -1, 'setWorkbenchMode does not appendChild (no reparent)');
    assert.strictEqual(fn.indexOf('removeChild'), -1, 'setWorkbenchMode does not removeChild (no reparent)');
    assert.strictEqual(fn.indexOf('insertBefore'), -1, 'setWorkbenchMode does not insertBefore (no reparent)');
    assert.strictEqual(fn.indexOf('bime-editor'), -1,
        'setWorkbenchMode never re-mounts #bime-editor — it only stamps the .wb-view-* class');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
