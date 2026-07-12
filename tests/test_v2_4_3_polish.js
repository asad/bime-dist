/**
 * tests/test_v2_4_3_polish.js — Hybrid Stage D: polish + cleanup (v2.4.3).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * The Hybrid is functionally complete (v2.4.0 Editor|Pathway view switch;
 * v2.4.1 "Add to pathway" import; v2.4.2 AAM-link so atom-tracing chains across
 * user-built pathways). Stage D is the finishing polish. No chemistry or
 * algorithm change. This suite pins four small, additive items:
 *
 *   1. Unnamed-molecule auto-tag -> formula fallback. captureEditorStructureForPathway
 *      routes a PLAIN MOLECULE with no editor/mol name through the existing
 *      moleculeFormulaLabel(...) helper BEFORE the literal 'Structure', so an
 *      unnamed ethanol lands as 'C2H6O' rather than the generic "Structure".
 *      The reaction path (named-or-'Reaction') and the named-molecule path are
 *      untouched.
 *   2. First-run discoverability. The empty-canvas CTA (drawEmptyCanvasCta) keeps
 *      its v2.3.5 double-click + Example copy AND adds a second terse line naming
 *      the Editor tab as a third way in (draw there, then Add to pathway). Still
 *      pointer-events:none + aria-hidden + --color tokens.
 *   3. docs.html no longer hardcodes the brittle "1,415" test count (it drifts
 *      EVERY release); the prose is now count-free.
 *   4. Mobile coherence. The Editor|Pathway tablist is reachable at <=720px with
 *      >=44px targets; the standalone Editor view's "Add to pathway" button is NOT
 *      hidden by the mobile output-action trimming; and the dock (a Pathway-canvas
 *      affordance) is hidden in the Editor view via a dedicated, lens-orthogonal
 *      flag so Draw/Fit are not redundant/confusing there — while Draw still opens
 *      the lens in the Pathway view (v2.3.6 dock contract preserved).
 *
 * Source-shape + DOM-stub, plain Node (no real DOM beyond the shim) — mirroring
 * test_v2_4_0 / test_v2_4_1 / test_v2_3_5 / test_v2_3_6.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Hybrid polish + cleanup (v2.4.3)');
var test = runner.test;
console.log('Hybrid polish + cleanup (v2.4.3)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var WB = readRepoFile('js/workbench.js');
var CSS = readRepoFile('css/workbench.css');
var HTML = readRepoFile('workbench.html');
var DOCS = readRepoFile('docs.html');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing test_v2_4_0 / test_v2_4_1 / test_v2_0_67 use,
// so a refactor that splits a pinned body fails here too.
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

// Slice the <=720px media block so mobile assertions key off the MOBILE rules.
function mobileCss() {
    var start = CSS.indexOf('@media screen and (max-width: 720px)');
    assert.ok(start !== -1, 'the <=720px media block is present');
    var i = CSS.indexOf('{', start);
    var depth = 0;
    for (; i < CSS.length; i++) {
        if (CSS[i] === '{') { depth++; }
        else if (CSS[i] === '}') { depth--; if (depth === 0) { return CSS.slice(start, i + 1); } }
    }
    return CSS.slice(start);
}

// ===========================================================================
// 1. Unnamed plain molecule -> molecular-formula label fallback.
// ===========================================================================
test('1.1 captureEditorStructureForPathway routes the unnamed-molecule fallback through moleculeFormulaLabel before \'Structure\'', function () {
    var fn = wbFn('captureEditorStructureForPathway');
    assert.ok(fn, 'captureEditorStructureForPathway located');
    // The plain-molecule branch derives a formula via the EXISTING helper, with
    // 'Structure' kept only as the final fallback. The empty-string second arg
    // ensures moleculeFormulaLabel does NOT itself short-circuit to a fallback
    // before the literal 'Structure'.
    assert.ok(/moleculeFormulaLabel\(mol,\s*''\)\s*\|\|\s*'Structure'/.test(fn),
        'unnamed plain molecule falls back to moleculeFormulaLabel(mol, \'\') before \'Structure\'');
    // The reaction path is unchanged: named-or-'Reaction' (NOT routed through formula).
    assert.ok(fn.indexOf("title = name || 'Reaction'") !== -1,
        'the reaction path keeps the named-or-\'Reaction\' fallback (no formula)');
    // The named path is honoured first for both molecule + reaction.
    assert.ok(/var name = \(editor\.getMolName && editor\.getMolName\(\)\) \|\| mol\.name \|\| ''/.test(fn),
        'the editor/mol name is honoured first (named molecules + reactions unaffected)');
});

test('1.2 captureEditorStructureForPathway stays a single top-level function (source-shape pin)', function () {
    var fn = wbFn('captureEditorStructureForPathway');
    assert.strictEqual(fn.indexOf('\nfunction '), -1,
        'captureEditorStructureForPathway is grown in place, never split/relocated');
    // The pinned payload literals (v2.0.60 A2b) still live in the function.
    assert.ok(fn.indexOf('smiles: smiles') !== -1 && fn.indexOf('rxn: rxnfile') !== -1,
        'the structured payload literals are preserved');
});

// Executable DOM-stub: a faithful port of the label-fallback control flow,
// proving the fallback ORDER (name -> formula -> 'Structure' for a molecule;
// name -> 'Reaction' for a reaction).
function labelFallbackModel(opts) {
    // opts: { isReaction, name, formula } — `formula` is what moleculeFormulaLabel
    // would return for an empty-string fallback ('' means "not derivable").
    var moleculeFormulaLabel = function (mol, fallback) {
        return opts.formula || fallback || '';
    };
    var name = opts.name || '';
    var title;
    if (opts.isReaction) {
        title = name || 'Reaction';
    } else {
        title = name || moleculeFormulaLabel(null, '') || 'Structure';
    }
    return title;
}

test('1.3 model: an unnamed CCO tags as its formula (C2H6O), not "Structure"', function () {
    assert.strictEqual(
        labelFallbackModel({ isReaction: false, name: '', formula: 'C2H6O' }),
        'C2H6O',
        'an unnamed plain molecule lands as its molecular formula');
});

test('1.4 model: a named molecule keeps its name; a formula-less molecule still falls back to "Structure"', function () {
    assert.strictEqual(
        labelFallbackModel({ isReaction: false, name: 'Ethanol', formula: 'C2H6O' }),
        'Ethanol',
        'a named molecule keeps its name (formula never overrides a real name)');
    assert.strictEqual(
        labelFallbackModel({ isReaction: false, name: '', formula: '' }),
        'Structure',
        'with no name and no derivable formula, the final fallback is still "Structure"');
});

test('1.5 model: the reaction path is unchanged (named-or-"Reaction", never formula)', function () {
    assert.strictEqual(
        labelFallbackModel({ isReaction: true, name: '', formula: 'C2H6O' }),
        'Reaction',
        'an unnamed reaction stays "Reaction" — the formula fallback is molecule-only');
    assert.strictEqual(
        labelFallbackModel({ isReaction: true, name: 'Esterification', formula: 'C2H6O' }),
        'Esterification',
        'a named reaction keeps its name');
});

// ===========================================================================
// 2. First-run CTA also surfaces the Editor tab.
// ===========================================================================
test('2.1 drawEmptyCanvasCta keeps its v2.3.5 copy AND names the Editor tab', function () {
    var fn = wbFn('drawEmptyCanvasCta');
    assert.ok(fn, 'drawEmptyCanvasCta located');
    // The v2.3.5 line is preserved verbatim (pinned by test_v2_3_5 P0-5.4/.5/E1).
    assert.ok(fn.indexOf('Double-click anywhere to draw a structure — or load an Example below') !== -1,
        'the original double-click + Example line is preserved verbatim');
    // The new line surfaces the Editor tab + the Add-to-pathway import flow.
    assert.ok(/Editor tab to draw/i.test(fn),
        'the CTA now mentions opening the Editor tab to draw');
    assert.ok(/Add to pathway/i.test(fn),
        'the CTA names the "Add to pathway" import flow as the follow-up');
});

test('2.2 the added CTA line stays on-brand: muted token, no emoji, in the same group', function () {
    var fn = wbFn('drawEmptyCanvasCta');
    // The new text node reuses the muted --color token + the CTA line class.
    var idx = fn.indexOf('Or open the Editor tab');
    assert.ok(idx !== -1, 'the Editor-tab line is present');
    // The fill token and class appear in the helper (defense-in-depth: the whole
    // helper uses --color-text-muted + wb-pathway-cta-line for its text nodes).
    assert.ok(fn.indexOf('var(--color-text-muted)') !== -1,
        'CTA text uses the --color-text-muted token');
    assert.ok(count(fn, "'class': 'wb-pathway-cta-line'") >= 2,
        'the new line reuses the wb-pathway-cta-line class (>=2 lines now)');
    // No emoji surrogate pairs anywhere in the helper (matches test_v2_3_5 P0-5.5).
    assert.ok(!/[\uD800-\uDBFF]/.test(fn), 'CTA copy contains no emoji');
    // Still pointer-events:none + aria-hidden (the group-level a11y contract).
    assert.ok(fn.indexOf("'pointer-events': 'none'") !== -1 && fn.indexOf("'aria-hidden': 'true'") !== -1,
        'the CTA group stays pointer-events:none + aria-hidden');
});

test('2.3 the default boot view stays gated by PATHWAY_ENABLED (CTA + ?pathway= protected when enabled)', function () {
    // v2.4.18: the default view now follows the pathway feature flag. Re-pin the
    // gated declaration so a stray edit here is caught beside the CTA change.
    assert.ok(/var\s+_wbView\s*=\s*PATHWAY_ENABLED\s*\?\s*'pathway'\s*:\s*'editor'\s*;/.test(WB),
        "_wbView = PATHWAY_ENABLED ? 'pathway' : 'editor' (flag-gated default boot view)");
});

// ===========================================================================
// 3. docs.html no longer hardcodes the brittle "1,415" test count.
// ===========================================================================
test('3.1 docs.html does not hardcode the stale "1,415" test count', function () {
    assert.strictEqual(DOCS.indexOf('1,415'), -1,
        'docs.html must not hardcode 1,415 (the real count drifts every release)');
    // The previously-derived "4,280 test executions" number was equally brittle.
    assert.strictEqual(DOCS.indexOf('4,280'), -1,
        'docs.html must not hardcode the derived 4,280 total either');
});

test('3.2 docs.html still describes the multi-tier suite (source / bundle / min-bundle in lockstep)', function () {
    // Non-brittle phrasing: the Testing section keeps describing all four tiers
    // without committing to an exact, drift-prone integer.
    var lower = DOCS.toLowerCase();
    assert.ok(lower.indexOf('source') !== -1 && lower.indexOf('bundle') !== -1,
        'the Testing section still names the source + bundle tiers');
    assert.ok(lower.indexOf('minified-bundle') !== -1 || lower.indexOf('min-bundle') !== -1,
        'the Testing section still names the minified-bundle tier');
    assert.ok(lower.indexOf('regression') !== -1 || lower.indexOf('lockstep') !== -1,
        'the prose frames it as a comprehensive regression battery (count-free)');
});

// ===========================================================================
// 4. Mobile: the tablist + dock are coherent at <=720px.
// ===========================================================================
test('4.1 the Editor|Pathway tablist has >=44px touch targets in the mobile block', function () {
    var m = mobileCss();
    var i = m.indexOf('.wb-mode-btn');
    assert.ok(i !== -1, 'the mobile block tunes the mode-btn tabs');
    var block = m.slice(i, i + 120);
    assert.ok(/min-height:\s*44px/.test(block),
        'the Editor|Pathway tabs are >=44px tall on mobile (WCAG 2.5.5)');
    // The modebar is display:flex at every width (v2.4.0 C5 pins the base rule);
    // re-assert it is NOT display:none in the mobile block.
    var mb = m.indexOf('.wb-modebar');
    if (mb !== -1) {
        var mbBlock = m.slice(mb, m.indexOf('}', mb) + 1);
        assert.strictEqual(/display:\s*none/.test(mbBlock), false,
            'the modebar is not hidden on mobile (the view switch stays reachable)');
    }
});

test('4.2 the Editor-view "Add to pathway" button is NOT hidden by the mobile output trimming', function () {
    // The mobile block trims export/validate/clear-editor from the Output action
    // row, but the Add-to-pathway control must stay reachable on a phone.
    var m = mobileCss();
    assert.strictEqual(m.indexOf('[data-wb-action="editor-add-to-pathway"]'), -1,
        'no mobile rule hides editor-add-to-pathway');
    // Sanity: the trimming that DOES exist targets only export/validate/clear.
    assert.ok(m.indexOf('[data-wb-action="validate"]') !== -1,
        'the mobile output trimming still targets validate (unchanged)');
    // The button still exists exactly once in the page (v2.4.1 A1 count holds).
    assert.strictEqual(count(HTML, 'data-wb-action="editor-add-to-pathway"'), 1,
        'the Add-to-pathway button is present exactly once (unchanged)');
});

test('4.3 setWorkbenchMode hides the dock in the Editor view via a lens-orthogonal flag', function () {
    var fn = wbFn('setWorkbenchMode');
    assert.ok(fn, 'setWorkbenchMode located');
    assert.ok(fn.indexOf("setMobileDockEditorHidden(_wbView === 'editor')") !== -1,
        'setWorkbenchMode hides the dock iff the Editor view is active');
    // The helper toggles its OWN flag (NOT the lens-sheet .wb-dock-hidden), so
    // the two state machines never fight.
    assert.ok(WB.indexOf('function setMobileDockEditorHidden(') !== -1,
        'setMobileDockEditorHidden is a dedicated helper');
    var helper = wbFn('setMobileDockEditorHidden');
    assert.ok(helper.indexOf("classList.add('wb-dock-editor-hidden')") !== -1
        && helper.indexOf("classList.remove('wb-dock-editor-hidden')") !== -1,
        'the helper toggles the dedicated .wb-dock-editor-hidden class');
    assert.strictEqual(helper.indexOf('wb-dock-hidden\''), -1,
        'the Editor-view flag does NOT touch the lens-sheet .wb-dock-hidden flag');
    // setWorkbenchMode stays a pure class-flip (v2.4.0 D5): no reparenting.
    assert.strictEqual(fn.indexOf('appendChild'), -1, 'setWorkbenchMode does not appendChild');
    assert.strictEqual(fn.indexOf('removeChild'), -1, 'setWorkbenchMode does not removeChild');
});

test('4.4 CSS: .wb-dock-editor-hidden hides the dock in the mobile block (distinct from .wb-dock-hidden)', function () {
    var m = mobileCss();
    var i = m.indexOf('.wb-mobile-dock.wb-dock-editor-hidden');
    assert.ok(i !== -1, 'the Editor-view dock-hidden rule is inside the <=720px block');
    var block = m.slice(i, i + 80);
    assert.ok(block.indexOf('display: none !important') !== -1,
        '.wb-dock-editor-hidden hides the dock');
    // The v2.3.6 lens-sheet rule still exists and is independent.
    assert.ok(m.indexOf('.wb-mobile-dock.wb-dock-hidden') !== -1,
        'the v2.3.6 lens-sheet .wb-dock-hidden rule is preserved (orthogonal)');
});

test('4.5 the v2.3.6 dock contract is intact: Draw still opens the lens in the Pathway view', function () {
    // Stage D must not break the v2.3.6 dock: Draw -> dock-draw -> the lens-draw
    // entry point; Fit -> pathway-fit; Inspect -> the inspector toggle.
    var dockStart = HTML.indexOf('<div class="wb-mobile-dock"');
    assert.ok(dockStart !== -1, 'mobile dock markup present');
    var dock = HTML.slice(dockStart, HTML.indexOf('</div>', dockStart) + 6);
    assert.ok(dock.indexOf('data-wb-action="dock-draw"') !== -1, 'dock keeps the Draw button');
    assert.ok(dock.indexOf('data-wb-action="pathway-fit"') !== -1, 'dock keeps the canvas Fit');
    assert.ok(dock.indexOf('data-wb-action="toggle-inspector"') !== -1, 'dock keeps the Inspect toggle');
    var idx = WB.indexOf("action === 'dock-draw'");
    assert.ok(idx !== -1, 'dispatcher still handles dock-draw');
    var arm = WB.slice(idx, idx + 700);
    assert.ok(arm.indexOf('addCurrentMoleculeToPathway()') !== -1,
        'dock-draw still routes to the lens-draw entry point (unchanged)');
});

// E1: an executable model of the dock-flag interaction proving the Editor-view
// flag and the lens-sheet flag are independent (collapsing the lens cannot
// reveal the dock while the Editor view owns the screen, and vice versa).
function recDock() {
    return {
        _set: {},
        classList: {
            _set: {},
            add: function (c) { this._set[c] = true; },
            remove: function (c) { delete this._set[c]; },
            contains: function (c) { return !!this._set[c]; }
        }
    };
}

test('E1 model: the Editor-view dock flag and the lens-sheet dock flag are orthogonal', function () {
    var dock = recDock();
    function setSheetHidden(h) { if (h) dock.classList.add('wb-dock-hidden'); else dock.classList.remove('wb-dock-hidden'); }
    function setEditorHidden(h) { if (h) dock.classList.add('wb-dock-editor-hidden'); else dock.classList.remove('wb-dock-editor-hidden'); }
    function dockVisible() { return !dock.classList.contains('wb-dock-hidden') && !dock.classList.contains('wb-dock-editor-hidden'); }

    // Enter the Editor view: dock hidden by the Editor flag.
    setEditorHidden(true);
    assert.strictEqual(dockVisible(), false, 'Editor view hides the dock');
    // Open + collapse a lens while in the Editor view: the lens flag toggles, but
    // collapsing it must NOT reveal the dock (the Editor flag still holds it).
    setSheetHidden(true);
    setSheetHidden(false);
    assert.strictEqual(dockVisible(), false,
        'collapsing a lens cannot reveal the dock while the Editor view holds it');
    // Return to the Pathway view: the dock comes back.
    setEditorHidden(false);
    assert.strictEqual(dockVisible(), true, 'Pathway view restores the dock');
    // A lens sheet in the Pathway view still hides + restores the dock on its own.
    setSheetHidden(true);
    assert.strictEqual(dockVisible(), false, 'a lens sheet hides the dock in the Pathway view');
    setSheetHidden(false);
    assert.strictEqual(dockVisible(), true, 'collapsing the lens restores the dock in the Pathway view');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
