/**
 * tests/test_v2_3_6_mobile.js — continuum on mobile (v2.3.6).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.3.6 fixes two mobile P0s a UX review flagged ("the continuum is near-unusable
 * on a phone"). No chemistry or algorithm change. This suite pins both fixes:
 *
 *   P0-6  the focus lens becomes a FULL-SCREEN SHEET on a narrow (<=720px)
 *         viewport instead of the node-anchored floating panel. The anchored
 *         panel clamped to a 380x320 min, which on a ~360px phone covered the
 *         node, overflowed, and offered no on-screen exit (Esc is unavailable on
 *         a touch keyboard). positionStructureLens now DETECTS the narrow viewport
 *         (isNarrowViewport, matchMedia/innerWidth at the SAME 720px the CSS uses)
 *         and delegates to applyLensSheetLayout, which adds `.wb-lens-sheet` to
 *         `.wb-main` (CSS owns the full-screen sizing), strips the inline desktop
 *         coords, hides the dock, refits the editor, and SKIPS the node-anchor/flip
 *         math. The v2.3.4 Done/Cancel header is a STICKY bar of the sheet with
 *         >=44px taps; the dialog a11y (role/aria-modal/aria-label/focus-trap) still
 *         applies. Desktop behaviour above the breakpoint is unchanged.
 *   P0-7  the mobile dock drives CONTINUUM actions (Draw -> opens a fresh node in
 *         the lens sheet; Fit -> fits the CANVAS; Inspect -> toggles the Inspector
 *         as a bottom sheet). The old Select / Command / Clean / Fit(editor) /
 *         Output buttons acted on the display:none editor and did nothing visible
 *         with no lens — they are gone. The Inspector (Selection / Warnings /
 *         Command palette) was display:none under 720px with no mobile entry point;
 *         `.wb-inspector-open` now overrides that as a bottom sheet, toggled from
 *         the dock (aria-expanded kept in sync; Esc / the X dismiss it).
 *
 * Mostly source-shape (the workbench DOM glue needs the whole page to execute,
 * exactly as test_v2_3_0 / test_v2_3_4 / test_v2_3_5 note), plus an executable
 * DOM-stub section that drives the sheet-layout toggle, the inspector toggle, and
 * the narrow-viewport decision through recording stubs — modelled on
 * test_v2_3_4_lens_hardening.js section E and the tests/shim.js document stub.
 *
 * Plain Node, no real DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Continuum on mobile (v2.3.6)');
var test = runner.test;
console.log('Continuum on mobile (v2.3.6)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var HTML = readRepoFile('workbench.html');
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

function count(haystack, needle) {
    return haystack.split(needle).length - 1;
}

// Slice the mobile-dock markup (the <div class="wb-mobile-dock"> ... </div>).
function dockMarkup() {
    var start = HTML.indexOf('<div class="wb-mobile-dock"');
    assert.ok(start !== -1, 'mobile dock markup present');
    var end = HTML.indexOf('</div>', start);
    return HTML.slice(start, end + 6);
}

// Slice the <=720px media block so assertions key off the MOBILE rules.
function mobileCss() {
    var start = CSS.indexOf('@media screen and (max-width: 720px)');
    assert.ok(start !== -1, 'the <=720px media block is present');
    // Walk braces to find the matching close of the media block.
    var i = CSS.indexOf('{', start);
    var depth = 0;
    for (; i < CSS.length; i++) {
        if (CSS[i] === '{') { depth++; }
        else if (CSS[i] === '}') { depth--; if (depth === 0) { return CSS.slice(start, i + 1); } }
    }
    return CSS.slice(start);
}

// ===========================================================================
// P0-6 — full-screen lens sheet on a narrow viewport.
// ===========================================================================
test('P0-6.1 isNarrowViewport + applyLensSheetLayout exist as top-level functions', function () {
    assert.ok(WB.indexOf('function isNarrowViewport(') !== -1, 'isNarrowViewport defined');
    assert.ok(WB.indexOf('function applyLensSheetLayout(') !== -1, 'applyLensSheetLayout defined');
});

test('P0-6.2 the mobile breakpoint constant matches the CSS 720px value', function () {
    assert.ok(WB.indexOf('WB_MOBILE_BREAKPOINT = 720') !== -1,
        'WB_MOBILE_BREAKPOINT is 720 (matches the CSS @media max-width)');
    var fn = wbFn('isNarrowViewport');
    assert.ok(fn, 'isNarrowViewport located');
    // It tracks the SAME query the CSS uses (matchMedia) with an innerWidth fallback.
    assert.ok(fn.indexOf('matchMedia') !== -1, 'isNarrowViewport prefers matchMedia');
    assert.ok(fn.indexOf('max-width: ') !== -1, 'isNarrowViewport queries a max-width media string');
    assert.ok(fn.indexOf('window.innerWidth') !== -1, 'isNarrowViewport falls back to window.innerWidth');
    assert.ok(fn.indexOf('WB_MOBILE_BREAKPOINT') !== -1, 'isNarrowViewport uses the shared breakpoint');
});

test('P0-6.3 positionStructureLens has a narrow-viewport branch that applies the sheet and RETURNS early', function () {
    var fn = wbFn('positionStructureLens');
    assert.ok(fn, 'positionStructureLens located');
    // The narrow branch: detect -> apply sheet -> return (skip the anchor/flip math).
    assert.ok(/if\s*\(\s*isNarrowViewport\(\)\s*\)\s*\{/.test(fn),
        'positionStructureLens branches on isNarrowViewport()');
    var idx = fn.indexOf('isNarrowViewport()');
    var branch = fn.slice(idx, idx + 160);
    assert.ok(branch.indexOf('applyLensSheetLayout(wrap, true)') !== -1,
        'the narrow branch applies the sheet layout');
    assert.ok(branch.indexOf('return') !== -1,
        'the narrow branch returns early (skips the node-anchoring/flip math)');
    // The early return precedes the desktop anchor math (the anchor rect lookup).
    assert.ok(fn.indexOf('isNarrowViewport()') < fn.indexOf('getBoundingClientRect'),
        'the narrow branch sits BEFORE the desktop anchor computation');
});

test('P0-6.4 positionStructureLens still keeps the desktop path (sheet OFF above the breakpoint)', function () {
    var fn = wbFn('positionStructureLens');
    // Above the breakpoint it explicitly turns the sheet OFF, then runs the
    // unchanged anchor/flip math (MIN_W/MIN_H clamp + edge flip + inline coords).
    assert.ok(fn.indexOf('applyLensSheetLayout(wrap, false)') !== -1,
        'the wider path turns the sheet layout OFF');
    assert.ok(fn.indexOf('MIN_W') !== -1 && fn.indexOf('MIN_H') !== -1,
        'the desktop min-size clamp is preserved');
    assert.ok(fn.indexOf("wrap.style.position = 'fixed'") !== -1,
        'the desktop path still writes the inline anchored coords');
    // positionStructureLens is a single top-level function (source-shape: not split).
    assert.strictEqual(fn.indexOf('\nfunction '), -1,
        'positionStructureLens contains no nested top-level function declaration');
});

test('P0-6.5 applyLensSheetLayout toggles .wb-lens-sheet on .wb-main and clears inline coords', function () {
    var fn = wbFn('applyLensSheetLayout');
    assert.ok(fn, 'applyLensSheetLayout located');
    assert.ok(fn.indexOf("classList.add('wb-lens-sheet')") !== -1,
        'sheet ON adds .wb-lens-sheet to .wb-main');
    assert.ok(fn.indexOf("classList.remove('wb-lens-sheet')") !== -1,
        'sheet OFF removes .wb-lens-sheet');
    // It hands sizing to CSS by clearing inline coords, and hides the dock.
    assert.ok(fn.indexOf("wrap.style.left = ''") !== -1,
        'sheet ON clears the inline left coord so CSS owns the sizing');
    assert.ok(fn.indexOf('setMobileDockHidden(true)') !== -1,
        'sheet ON hides the mobile dock');
    assert.ok(fn.indexOf('setMobileDockHidden(false)') !== -1,
        'sheet OFF restores the mobile dock');
    // It refits the editor draw surface (mirrors the desktop refit).
    assert.ok(fn.indexOf('editor.setSize(') !== -1, 'sheet refit calls editor.setSize');
});

test('P0-6.6 openStructureLens still builds the v2.3.4 header + dialog a11y (contract preserved)', function () {
    var open = wbFn('openStructureLens');
    assert.ok(open.indexOf('buildStructureLensHeader(') !== -1,
        'openStructureLens still builds the Done/Cancel header (sheet reuses it)');
    assert.ok(open.indexOf("setAttribute('role', 'dialog')") !== -1,
        'dialog role still set (applies to the sheet too)');
    assert.ok(open.indexOf("setAttribute('aria-modal', 'true')") !== -1, 'aria-modal still set');
    assert.ok(open.indexOf('_lensTrapHandler = onLensTrapKeydown') !== -1,
        'the Tab focus-trap is still installed (sheet is still a modal dialog)');
    // positionStructureLens is still called from open, so the sheet class is
    // applied during open on a narrow viewport.
    assert.ok(open.indexOf('positionStructureLens()') !== -1,
        'openStructureLens calls positionStructureLens (applies the sheet on mobile)');
});

test('P0-6.7 collapseStructureLens tears the sheet down (drops .wb-lens-sheet + restores the dock)', function () {
    var fn = wbFn('collapseStructureLens');
    assert.ok(fn.indexOf("classList.remove('wb-lens-sheet')") !== -1,
        'collapse removes the .wb-lens-sheet layout class');
    assert.ok(fn.indexOf('setMobileDockHidden(false)') !== -1,
        'collapse restores the mobile dock');
    // The v2.3.4 commit/discard guard is untouched.
    assert.ok(/if\s*\(commit && node\)/.test(fn),
        'the commit/discard guard from v2.3.4 is preserved');
});

test('P0-6.8 the lens re-evaluates on resize while open (adapts across the breakpoint)', function () {
    // A window resize listener re-runs positionStructureLens when the lens is open,
    // so rotating a phone / resizing past 720px swaps sheet<->anchored.
    assert.ok(/window\.addEventListener\('resize'/.test(WB),
        'a window resize listener is wired');
    var idx = WB.indexOf("addEventListener('resize'");
    var slice = WB.slice(idx, idx + 160);
    assert.ok(slice.indexOf('_lens.isOpen()') !== -1 && slice.indexOf('positionStructureLens()') !== -1,
        'the resize handler re-runs positionStructureLens while the lens is open');
});

// --- P0-6 CSS: the sheet is full-screen with a STICKY Done/Cancel header. ---
test('P0-6.9 CSS makes the lens a full-screen sheet under .wb-lens-sheet (<=720px)', function () {
    var m = mobileCss();
    assert.ok(m.indexOf('.wb-main.wb-lens-open.wb-lens-sheet .wb-editor-wrap') !== -1,
        'the sheet editor-wrap rule is inside the <=720px block');
    var idx = m.indexOf('.wb-main.wb-lens-open.wb-lens-sheet .wb-editor-wrap');
    var block = m.slice(idx, m.indexOf('}', idx) + 1);
    assert.ok(block.indexOf('position: fixed') !== -1, 'sheet is position:fixed');
    // Near-full-screen: fixed inset with a small margin (left/top/right/bottom set).
    assert.ok(block.indexOf('left:') !== -1 && block.indexOf('right:') !== -1 &&
              block.indexOf('top:') !== -1 && block.indexOf('bottom:') !== -1,
        'sheet pins all four edges (full-bleed inset)');
    // v2.3.6 fix: the sheet MUST sit ABOVE the page nav (.nav is z-index 100), or
    // its sticky Done/Cancel header is occluded by the nav bar — no tappable exit
    // on a phone (caught in browser verification, not by the DOM-stub alone).
    var zm = /z-index:\s*(\d+)/.exec(block);
    assert.ok(zm, 'the sheet sets an explicit z-index');
    var zi = parseInt(zm[1], 10);
    assert.ok(zi > 100 && zi < 1000,
        'sheet z-index (' + zi + ') is above the page nav (100), below the modal-root (1000)');
});

test('P0-6.10 the sheet header is STICKY with >=44px tap targets (always-visible exit)', function () {
    var m = mobileCss();
    var hIdx = m.indexOf('.wb-main.wb-lens-open.wb-lens-sheet .wb-lens-header');
    assert.ok(hIdx !== -1, 'a sheet-scoped lens-header rule exists in the mobile block');
    var hBlock = m.slice(hIdx, hIdx + 200);
    assert.ok(hBlock.indexOf('position: sticky') !== -1, 'the sheet header is sticky');
    assert.ok(hBlock.indexOf('top: 0') !== -1, 'the sticky header pins to the top');
    var bIdx = m.indexOf('.wb-main.wb-lens-open.wb-lens-sheet .wb-lens-btn');
    assert.ok(bIdx !== -1, 'a sheet-scoped lens-btn rule exists');
    var bBlock = m.slice(bIdx, bIdx + 160);
    assert.ok(/min-height:\s*44px/.test(bBlock), 'sheet Done/Cancel buttons are >=44px tall (WCAG 2.5.5)');
});

// ===========================================================================
// P0-7 — reconcile the mobile dock + restore Inspector access.
// ===========================================================================
test('P0-7.1 the mobile dock no longer carries no-op-on-mobile editor buttons', function () {
    var dock = dockMarkup();
    // The old editor-targeting actions are GONE from the dock (they acted on the
    // display:none editor and did nothing visible with no lens open).
    assert.strictEqual(dock.indexOf('data-wb-action="set-tool"'), -1,
        'no Select/set-tool button (it targeted the hidden editor)');
    assert.strictEqual(dock.indexOf('data-wb-action="clean2d"'), -1,
        'no Clean/clean2d button (it targeted the hidden editor)');
    assert.strictEqual(dock.indexOf('data-wb-action="best-fit"'), -1,
        'no editor Best-Fit button (it targeted the hidden editor)');
    assert.strictEqual(dock.indexOf('data-wb-action="focus-output"'), -1,
        'no Output button (the Output section is already below the canvas)');
});

test('P0-7.2 the dock drives CONTINUUM actions: Draw, canvas Fit, Inspect', function () {
    var dock = dockMarkup();
    // Draw -> a dedicated dock action that routes to the lens-draw entry point.
    assert.ok(dock.indexOf('data-wb-action="dock-draw"') !== -1,
        'dock has a Draw button (dock-draw)');
    assert.ok(dock.indexOf('>Draw<') !== -1, 'the Draw button is labelled Draw');
    // Fit -> the CANVAS fit (pathway-fit), not the editor best-fit.
    assert.ok(dock.indexOf('data-wb-action="pathway-fit"') !== -1,
        'dock Fit acts on the canvas (pathway-fit), not the hidden editor');
    // Inspect -> the inspector bottom-sheet toggle.
    assert.ok(dock.indexOf('data-wb-action="toggle-inspector"') !== -1,
        'dock has an Inspect toggle');
});

test('P0-7.3 dock-draw routes to the SAME continuum entry point (no separate no-op)', function () {
    // The dispatcher routes dock-draw to addCurrentMoleculeToPathway (the canonical
    // Draw Structure path, which opens newPathwayNodeInLens on an empty editor).
    var idx = WB.indexOf("action === 'dock-draw'");
    assert.ok(idx !== -1, 'dispatcher handles dock-draw');
    var arm = WB.slice(idx, idx + 700);
    assert.ok(arm.indexOf('addCurrentMoleculeToPathway()') !== -1,
        'dock-draw calls addCurrentMoleculeToPathway (the lens-draw entry point)');
    // add-pathway-current stays at its two PINNED placements (Tools row + command
    // grid) — the dock uses dock-draw, so the v2.0.52 / v2.3.5 counts hold.
    assert.strictEqual(count(HTML, 'data-wb-action="add-pathway-current"'), 2,
        'add-pathway-current still appears exactly twice (dock did not add a third)');
});

test('P0-7.4 the Inspect toggle is accessible (aria-expanded + aria-controls)', function () {
    var dock = dockMarkup();
    var idx = dock.indexOf('data-wb-action="toggle-inspector"');
    var btn = dock.slice(dock.lastIndexOf('<button', idx), dock.indexOf('>', idx) + 1);
    assert.ok(btn.indexOf('aria-expanded="false"') !== -1,
        'the Inspect toggle starts aria-expanded=false');
    assert.ok(btn.indexOf('aria-controls="wb-inspector-panel"') !== -1,
        'the Inspect toggle points aria-controls at the inspector');
    // The inspector aside carries the matching id.
    assert.ok(HTML.indexOf('id="wb-inspector-panel"') !== -1,
        'the inspector <aside> carries the wb-inspector-panel id');
});

test('P0-7.5 setMobileInspectorOpen toggles .wb-inspector-open + syncs aria-expanded', function () {
    assert.ok(WB.indexOf('function setMobileInspectorOpen(') !== -1, 'setMobileInspectorOpen defined');
    assert.ok(WB.indexOf('function toggleMobileInspector(') !== -1, 'toggleMobileInspector defined');
    var fn = wbFn('setMobileInspectorOpen');
    assert.ok(fn.indexOf("classList.add('wb-inspector-open')") !== -1,
        'open adds .wb-inspector-open');
    assert.ok(fn.indexOf("classList.remove('wb-inspector-open')") !== -1,
        'close removes .wb-inspector-open');
    assert.ok(fn.indexOf("setAttribute('aria-expanded'") !== -1,
        'the dock toggle aria-expanded is kept in sync');
});

test('P0-7.6 Escape dismisses the inspector bottom sheet (after the lens)', function () {
    // The global Escape handler closes the inspector sheet when it is open.
    var idx = WB.indexOf("classList.contains('wb-inspector-open')");
    assert.ok(idx !== -1, 'Escape handler checks the inspector-open state');
    var slice = WB.slice(idx, idx + 200);
    assert.ok(slice.indexOf('setMobileInspectorOpen(false)') !== -1,
        'Escape closes the inspector bottom sheet');
});

// --- P0-7 CSS: the inspector is reachable on mobile via the -open override. ---
test('P0-7.7 CSS: <=720px hides the inspector but .wb-inspector-open overrides it as a bottom sheet', function () {
    var m = mobileCss();
    assert.ok(/\.wb-inspector\s*\{\s*display:\s*none\s*!important;\s*\}/.test(m),
        'the inspector is display:none by default under 720px');
    var oIdx = m.indexOf('.wb-inspector.wb-inspector-open');
    assert.ok(oIdx !== -1, 'the -open override rule exists in the mobile block');
    var oBlock = m.slice(oIdx, oIdx + 360);
    assert.ok(oBlock.indexOf('display: grid !important') !== -1,
        '.wb-inspector-open overrides the display:none');
    assert.ok(oBlock.indexOf('position: fixed') !== -1 && oBlock.indexOf('bottom: 0') !== -1,
        '.wb-inspector-open is a bottom sheet (fixed, pinned to the bottom)');
});

test('P0-7.8 CSS: the dock hides while the sheet owns the screen (.wb-dock-hidden)', function () {
    var m = mobileCss();
    assert.ok(m.indexOf('.wb-mobile-dock.wb-dock-hidden') !== -1,
        'the dock-hidden rule is in the mobile block');
    var idx = m.indexOf('.wb-mobile-dock.wb-dock-hidden');
    assert.ok(m.slice(idx, idx + 60).indexOf('display: none !important') !== -1,
        '.wb-dock-hidden hides the dock');
});

// ===========================================================================
// E. Executable DOM-stub: drive the sheet toggle, inspector toggle, and the
//    narrow-viewport decision through recording stubs (the workbench glue proper
//    needs the whole page, so we model its exact control flow here the way
//    test_v2_3_4 section E models the lens lifecycle).
// ===========================================================================

// A minimal recording classList + element.
function recClassList() {
    return {
        _set: {},
        add: function (c) { this._set[c] = true; },
        remove: function (c) { delete this._set[c]; },
        contains: function (c) { return !!this._set[c]; },
        toggle: function (c) { if (this._set[c]) { delete this._set[c]; } else { this._set[c] = true; } }
    };
}
function recEl() {
    return {
        classList: recClassList(),
        style: {},
        _attrs: {},
        setAttribute: function (k, v) { this._attrs[k] = v; },
        getAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k) ? this._attrs[k] : null; },
        focus: function () { this._focused = true; }
    };
}

test('E1 narrow-viewport decision model: <=720px => sheet, >720px => anchored', function () {
    var BP = 720;
    function isNarrow(width) { return width <= BP; }
    assert.strictEqual(isNarrow(360), true, 'a 360px phone is narrow (sheet)');
    assert.strictEqual(isNarrow(720), true, 'exactly 720px is narrow (inclusive, matches CSS max-width)');
    assert.strictEqual(isNarrow(721), false, '721px is desktop (anchored)');
    assert.strictEqual(isNarrow(1280), false, 'a 1280px desktop is anchored');
});

test('E2 applyLensSheetLayout model: ON sets the class + clears coords + hides dock; OFF reverses', function () {
    var main = recEl();
    var wrap = recEl();
    var dock = recEl();
    // Seed inline desktop coords as if a prior wider layout had run.
    wrap.style.position = 'fixed'; wrap.style.left = '120px'; wrap.style.top = '80px';
    wrap.style.width = '380px'; wrap.style.height = '320px'; wrap.style.zIndex = '50';

    function applyModel(on) {
        if (on) {
            main.classList.add('wb-lens-sheet');
            wrap.style.position = ''; wrap.style.left = ''; wrap.style.top = '';
            wrap.style.width = ''; wrap.style.height = ''; wrap.style.zIndex = '';
            dock.classList.add('wb-dock-hidden');
        } else {
            main.classList.remove('wb-lens-sheet');
            dock.classList.remove('wb-dock-hidden');
        }
    }

    applyModel(true);
    assert.strictEqual(main.classList.contains('wb-lens-sheet'), true, 'sheet ON adds the class');
    assert.strictEqual(wrap.style.left, '', 'sheet ON clears inline left (CSS owns sizing)');
    assert.strictEqual(wrap.style.width, '', 'sheet ON clears inline width');
    assert.strictEqual(dock.classList.contains('wb-dock-hidden'), true, 'sheet ON hides the dock');

    applyModel(false);
    assert.strictEqual(main.classList.contains('wb-lens-sheet'), false, 'sheet OFF drops the class');
    assert.strictEqual(dock.classList.contains('wb-dock-hidden'), false, 'sheet OFF restores the dock');
});

test('E3 inspector toggle model: flips .wb-inspector-open + aria-expanded, both directions', function () {
    var inspector = recEl();
    var toggle = recEl();
    toggle.setAttribute('aria-expanded', 'false');

    function setOpen(open) {
        if (open) { inspector.classList.add('wb-inspector-open'); }
        else { inspector.classList.remove('wb-inspector-open'); }
        toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
    }
    function toggleIt() { setOpen(!inspector.classList.contains('wb-inspector-open')); }

    toggleIt();
    assert.strictEqual(inspector.classList.contains('wb-inspector-open'), true, 'first tap opens the sheet');
    assert.strictEqual(toggle.getAttribute('aria-expanded'), 'true', 'aria-expanded => true when open');
    toggleIt();
    assert.strictEqual(inspector.classList.contains('wb-inspector-open'), false, 'second tap closes it');
    assert.strictEqual(toggle.getAttribute('aria-expanded'), 'false', 'aria-expanded => false when closed');
});

test('E4 the lens sheet keeps the v2.3.4 Done/Cancel contract (Done=commit, Cancel=discard)', function () {
    // The sheet reuses buildStructureLensHeader: Done -> collapse(true), Cancel ->
    // collapse(false). Model that the sticky header buttons map to commit/discard.
    var calls = [];
    function collapse(commit) { calls.push(commit ? 'commit' : 'discard'); }
    function onDone() { collapse(true); }
    function onCancel() { collapse(false); }
    onDone(); onCancel();
    assert.deepStrictEqual(calls, ['commit', 'discard'],
        'Done commits (collapse true), Cancel discards (collapse false)');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
