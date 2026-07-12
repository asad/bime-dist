/**
 * tests/test_v2_3_8_toolbar_clarity.js — toolbar clarity + live lens output
 * (v2.3.8).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Fifth UX-hardening release after v2.3.4–2.3.7. No chemistry or algorithm
 * change. This suite pins both fixes a UX review flagged:
 *
 *   P1-2  TOOLBAR CLARITY. The pathway controls were a flat wall of
 *         near-identical text pills with NO iconography, so a chemist could not
 *         pre-attentively find a tool (contrast the all-icons molecule-editor
 *         toolbar). Each TOOLS-row button (Select / Metabolite / Residue /
 *         Arrow / Curly / Step / Note / Compartment / Draw Structure / Edit
 *         Structure) now carries a LEADING decorative inline-SVG icon
 *         (aria-hidden="true"), and the key VIEW actions (Zoom -/+, Fit, Clean
 *         Up) do too. The icons are additive: no button's data-wb-action /
 *         data-tool / id changed, no visible label TEXT changed, and the NUMBER
 *         of buttons is unchanged. A PRIMARY/secondary hierarchy is introduced:
 *         the constantly-used Select (home tool) + Draw Structure read as
 *         primary (.wb-pathway-tool-primary), the once-in-a-while actions stay
 *         quiet.
 *
 *   P1-4  LIVE OUTPUT BESIDE THE LENS. While a focus lens is open the Inspector
 *         / Editor Output dock BELOW the canvas (off-screen — especially on the
 *         mobile full-screen sheet), so a user drawing in the lens could not see
 *         the SMILES or warnings their edit produced without closing the lens
 *         (which commits). openStructureLens now builds a COMPACT live readout
 *         strip inside `.wb-editor-wrap` (aria-live="polite") showing the current
 *         SMILES (truncated) + a warning-COUNT badge, fed by the SAME
 *         change->output signal Editor Output uses (updateLensReadout at the tail
 *         of updateOutputNow — no extra recompute). collapseStructureLens tears
 *         it down (mirrors the v2.3.4 header).
 *
 * Source-shape assertions (the workbench glue needs the whole page to execute,
 * exactly as test_v2_3_0 / test_v2_3_4 / test_v2_3_6 / test_v2_3_7 note) PLUS an
 * executable DOM-stub section that drives the readout build / update / teardown
 * control flow through recording stubs (modelled on test_v2_3_4 section E).
 *
 * Plain Node, no real DOM beyond the shim.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Toolbar clarity + live lens output (v2.3.8)');
var test = runner.test;
console.log('Toolbar clarity + live lens output (v2.3.8)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var HTML = readRepoFile('workbench.html');
var CSS = readRepoFile('css/workbench.css');
var WB = readRepoFile('js/workbench.js');

// Slice a top-level function body from `function name(` up to the next top-level
// `\nfunction ` — the same slicing the source-shape-pinned suites use.
function wbFn(name) {
    var start = WB.indexOf('function ' + name + '(');
    if (start === -1) { return ''; }
    var rest = WB.slice(start);
    var next = rest.indexOf('\nfunction ');
    return next === -1 ? rest : rest.slice(0, next);
}

function count(haystack, needle) {
    return (haystack.match(new RegExp(needle.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'g')) || []).length;
}

// The ten Tools-row buttons, identified by their stable label TEXT + the action
// that must NOT have changed.
var TOOLS_ROW = [
    { label: 'Select',          attr: 'data-tool="select"' },
    { label: 'Metabolite',      attr: 'data-tool="metabolite"' },
    { label: 'Residue',         attr: 'data-tool="residue"' },
    { label: 'Arrow',           attr: 'data-tool="arrow"' },
    { label: 'Curly',           attr: 'data-tool="curly"' },
    { label: 'Step',            attr: 'data-tool="step"' },
    { label: 'Note',            attr: 'data-tool="note"' },
    { label: 'Compartment',     attr: 'data-tool="compartment"' }
];

// Return the open-tag-through-label slice of the button whose visible label is
// `label` AND that sits inside the Tools toolbar block (the first one).
function toolsRowBlock() {
    var start = HTML.indexOf('wb-pathway-tools" role="toolbar"');
    assert.ok(start > 0, 'pathway Tools toolbar block present');
    var end = HTML.indexOf('</div>', start);
    return HTML.substring(start, end);
}

// ===========================================================================
// A. P1-2 — every Tools-row button gains a leading inline-SVG icon, additively.
// ===========================================================================
test('A1 the Tools toolbar still holds exactly the same buttons (count unchanged)', function () {
    var block = toolsRowBlock();
    // 8 click-on-canvas tools + the 2 stamp actions (Draw / Edit Structure) = 10.
    assert.strictEqual(count(block, '<button'), 10,
        'Tools row has exactly 10 buttons (8 tools + Draw/Edit Structure) — number unchanged');
    assert.strictEqual(count(block, 'data-wb-action="set-pathway-tool"'), 8,
        'the 8 set-pathway-tool toggles are intact');
});

test('A2 each set-pathway-tool button carries a LEADING aria-hidden inline-SVG icon', function () {
    for (var i = 0; i < TOOLS_ROW.length; i++) {
        var t = TOOLS_ROW[i];
        var idx = HTML.indexOf(t.attr);
        assert.ok(idx > 0, t.label + ' tool button present (' + t.attr + ')');
        // From the data-tool attr to the closing label `>Label<` — capture the
        // button's inner markup (everything between the open-tag `>` and the
        // label text).
        var gt = HTML.indexOf('>', idx);
        var labelMark = HTML.indexOf('>' + t.label + '<', idx);
        assert.ok(labelMark > gt, t.label + ' label `>' + t.label + '<` follows the open tag');
        var inner = HTML.substring(gt + 1, labelMark + 1); // includes the trailing `>` of </svg>
        assert.ok(inner.indexOf('<svg') !== -1, t.label + ' button has an inline <svg> icon');
        assert.ok(/<svg[^>]*aria-hidden="true"/.test(inner),
            t.label + ' icon is aria-hidden="true" (decorative)');
        assert.ok(inner.indexOf('</svg>') !== -1, t.label + ' icon is a complete inline SVG');
    }
});

test('A3 the two Tools-row stamp actions (Draw / Edit Structure) also get a leading icon', function () {
    // Draw Structure (Tools row): the solid stamp-action; find its open tag.
    var drawIdx = HTML.indexOf('wb-pathway-stamp-action wb-pathway-tool-primary');
    assert.ok(drawIdx > 0, 'the Tools-row Draw Structure stamp button is present');
    var drawGt = HTML.indexOf('>', drawIdx);
    var drawLabel = HTML.indexOf('>Draw Structure<', drawIdx);
    var drawInner = HTML.substring(drawGt + 1, drawLabel + 1);
    assert.ok(/<svg[^>]*aria-hidden="true"/.test(drawInner),
        'Draw Structure (Tools row) has an aria-hidden inline-SVG icon');
    // Edit Structure (Tools row): the non-solid stamp action right after it.
    var editIdx = HTML.indexOf('data-wb-action="edit-pathway-structure"');
    var editGt = HTML.indexOf('>', editIdx);
    var editLabel = HTML.indexOf('>Edit Structure<', editIdx);
    var editInner = HTML.substring(editGt + 1, editLabel + 1);
    assert.ok(/<svg[^>]*aria-hidden="true"/.test(editInner),
        'Edit Structure (Tools row) has an aria-hidden inline-SVG icon');
});

test('A4 the key VIEW actions (Zoom -/+, Fit, Clean Up) gain a leading icon too', function () {
    var views = [
        { attr: 'data-wb-action="pathway-zoom" data-dir="out"', label: 'Zoom -' },
        { attr: 'data-wb-action="pathway-zoom" data-dir="in"',  label: 'Zoom +' },
        { attr: 'data-wb-action="pathway-fit"',                 label: 'Fit' },
        { attr: 'data-wb-action="pathway-layout"',              label: 'Clean Up' }
    ];
    for (var i = 0; i < views.length; i++) {
        var v = views[i];
        var idx = HTML.indexOf(v.attr);
        assert.ok(idx > 0, v.label + ' view action present');
        var gt = HTML.indexOf('>', idx);
        var labelMark = HTML.indexOf('>' + v.label + '<', idx);
        var inner = HTML.substring(gt + 1, labelMark + 1);
        assert.ok(/<svg[^>]*aria-hidden="true"/.test(inner),
            v.label + ' view action has an aria-hidden inline-SVG icon');
    }
});

test('A5 icons are ADDITIVE — no data-action, data-tool, or label TEXT changed', function () {
    // The literal labels + action/count invariants the pinned suites assert.
    assert.strictEqual(count(HTML, '>Draw Structure<'), 2,
        'both Draw Structure buttons keep the literal label (icon did not displace it)');
    assert.strictEqual(count(HTML, '>Edit Structure<'), 2, 'both Edit Structure labels intact');
    assert.strictEqual(count(HTML, 'data-wb-action="add-pathway-current"'), 2,
        'add-pathway-current still appears exactly twice');
    assert.strictEqual(count(HTML, 'data-wb-action="edit-pathway-structure"'), 2,
        'edit-pathway-structure still appears exactly twice');
    assert.strictEqual(count(HTML, 'data-wb-action="set-pathway-tool"'), 8,
        'the 8 set-pathway-tool actions are unchanged in count');
    // The tool data-tool tokens are all still present.
    var tools = ['select', 'metabolite', 'residue', 'arrow', 'curly', 'step', 'note', 'compartment'];
    for (var i = 0; i < tools.length; i++) {
        assert.ok(HTML.indexOf('data-tool="' + tools[i] + '"') !== -1,
            'data-tool="' + tools[i] + '" preserved');
    }
    // No emoji / external <img> assets crept in for icons — inline SVG only.
    var block = toolsRowBlock();
    assert.strictEqual(block.indexOf('<img'), -1, 'no external <img> icons in the Tools row');
});

// ===========================================================================
// B. P1-2 — a PRIMARY / secondary hierarchy distinction exists.
// ===========================================================================
test('B1 a primary-tool class exists and marks Select + Draw Structure', function () {
    assert.ok(HTML.indexOf('wb-pathway-tool-primary') !== -1,
        'the .wb-pathway-tool-primary distinction is used in the markup');
    // Select (the home tool) is primary.
    var selIdx = HTML.indexOf('data-tool="select"');
    var selOpen = HTML.substring(HTML.lastIndexOf('<button', selIdx), HTML.indexOf('>', selIdx));
    assert.ok(selOpen.indexOf('wb-pathway-tool-primary') !== -1,
        'Select (home tool) reads as primary');
    // Draw Structure (the constantly-used stamp) is primary AND keeps the
    // existing solid/stamp-action treatment (reuse, not replace).
    var drawIdx = HTML.indexOf('data-wb-action="add-pathway-current"');
    var drawOpen = HTML.substring(HTML.lastIndexOf('<button', drawIdx), HTML.indexOf('>', drawIdx));
    assert.ok(drawOpen.indexOf('wb-pathway-tool-primary') !== -1,
        'Draw Structure reads as primary');
    assert.ok(drawOpen.indexOf('wb-btn-solid') !== -1 &&
              drawOpen.indexOf('wb-pathway-stamp-action') !== -1,
        'Draw Structure keeps the existing solid + stamp-action patterns');
});

test('B2 the primary class is styled in CSS', function () {
    assert.ok(CSS.indexOf('.wb-pathway-tool-primary') !== -1,
        '.wb-pathway-tool-primary has a CSS rule (the primary emphasis)');
    // The Tools-row icons are sized by the shared `.wb-btn svg` rule (no per-icon
    // width/height inline), so the icon flips colour with the button.
    assert.ok(CSS.indexOf('.wb-btn svg') !== -1,
        'the shared .wb-btn svg sizing rule (icons inherit currentColor) is present');
});

// ===========================================================================
// C. P1-4 — openStructureLens builds the live readout; collapse tears it down.
// ===========================================================================
test('C1 the readout helpers exist as top-level functions', function () {
    assert.ok(WB.indexOf('function buildStructureLensReadout(') !== -1,
        'buildStructureLensReadout defined');
    assert.ok(WB.indexOf('function removeStructureLensReadout(') !== -1,
        'removeStructureLensReadout defined');
    assert.ok(WB.indexOf('function updateLensReadout(') !== -1,
        'updateLensReadout defined');
});

test('C2 buildStructureLensReadout builds an aria-live strip with SMILES + a warning count', function () {
    var fn = wbFn('buildStructureLensReadout');
    assert.ok(fn, 'buildStructureLensReadout located');
    assert.ok(fn.indexOf('wb-lens-readout') !== -1, 'the readout element class is present');
    assert.ok(fn.indexOf("setAttribute('aria-live', 'polite')") !== -1,
        'the readout is an aria-live="polite" region (screen-reader announces updates)');
    assert.ok(fn.indexOf('wb-lens-readout-smiles') !== -1, 'a SMILES element is built');
    assert.ok(fn.indexOf('wb-lens-readout-warn') !== -1, 'a warning-count element is built');
    // It appends to the editor wrap (so it lives inside the lens overlay/sheet).
    assert.ok(fn.indexOf(".wb-editor-wrap") !== -1 || fn.indexOf("wrap.appendChild") !== -1,
        'the strip is mounted inside .wb-editor-wrap');
});

test('C3 updateLensReadout reuses the EXISTING cached output signal (no new recompute)', function () {
    var fn = wbFn('updateLensReadout');
    assert.ok(fn, 'updateLensReadout located');
    // SMILES from the cached outputTextFor('smiles'); warnings from the cached
    // getMoleculeInsights() — the same two signals Editor Output already reads.
    assert.ok(fn.indexOf("outputTextFor('smiles')") !== -1,
        'SMILES comes from the cached outputTextFor (no fresh SMILES write)');
    assert.ok(fn.indexOf('getMoleculeInsights()') !== -1,
        'warnings come from the cached getMoleculeInsights (no fresh layout eval)');
    assert.ok(fn.indexOf('warnings.length') !== -1, 'it shows a warning COUNT');
    // It is a no-op unless a lens is open.
    assert.ok(fn.indexOf('_lens') !== -1 && fn.indexOf('isOpen()') !== -1,
        'updateLensReadout is a no-op when no lens is open');
});

test('C4 updateOutputNow drives the readout on the same per-change RAF path', function () {
    var fn = wbFn('updateOutputNow');
    assert.ok(fn, 'updateOutputNow located');
    assert.ok(fn.indexOf('updateLensReadout()') !== -1,
        'updateOutputNow calls updateLensReadout (hooked into the existing change->output path)');
});

test('C5 openStructureLens builds the readout; collapseStructureLens removes it', function () {
    var open = wbFn('openStructureLens');
    assert.ok(open.indexOf('buildStructureLensReadout(') !== -1,
        'openStructureLens builds the live readout strip');
    // Built alongside / after the v2.3.4 header (both inside the same overlay).
    assert.ok(open.indexOf('buildStructureLensHeader(') !== -1,
        'the v2.3.4 Done/Cancel header is still built (not regressed)');
    var collapse = wbFn('collapseStructureLens');
    assert.ok(collapse.indexOf('removeStructureLensReadout()') !== -1,
        'collapseStructureLens tears the readout strip down');
    assert.ok(collapse.indexOf('removeStructureLensHeader()') !== -1,
        'the v2.3.4 header teardown is still present (not regressed)');
});

test('C6 the readout is styled (desktop overlay + mobile sheet) with --color tokens, no emoji', function () {
    assert.ok(CSS.indexOf('.wb-lens-readout') !== -1, 'the readout strip is styled');
    assert.ok(CSS.indexOf('.wb-lens-readout-smiles') !== -1, 'the SMILES line is styled');
    assert.ok(CSS.indexOf('.wb-lens-readout-warn') !== -1, 'the warning badge is styled');
    // The has-warnings state uses the warning token (not a colour-only literal).
    var warnIdx = CSS.indexOf('.wb-lens-readout-warn.has-warnings');
    assert.ok(warnIdx !== -1, 'a has-warnings state exists');
    var warnBlock = CSS.slice(warnIdx, CSS.indexOf('}', warnIdx) + 1);
    assert.ok(warnBlock.indexOf('var(--color-warning)') !== -1,
        'the warning state uses the --color-warning token');
    // It must also style under the mobile sheet (Output is most out of reach there).
    assert.ok(CSS.indexOf('.wb-main.wb-lens-open.wb-lens-sheet .wb-lens-readout') !== -1,
        'the readout has a mobile-sheet rule (sticky bottom bar)');
});

test('C7 the warning glyph is a text/SVG mark, not an emoji', function () {
    var fn = wbFn('updateLensReadout');
    // We use the U+26A0 WARNING SIGN (a plain symbol char) + the word "warning(s)";
    // no emoji presentation selector / pictographic emoji is used.
    assert.ok(fn.indexOf('warning') !== -1, 'the badge spells the word "warning(s)"');
    // No common emoji codepoints in the readout strings.
    assert.strictEqual(/[\u{1F300}-\u{1FAFF}\u{2700}-\u{27BF}]/u.test(fn), false,
        'no pictographic emoji in the readout copy');
});

// ===========================================================================
// D. Executable DOM-stub — drive the readout build / update / teardown control
//    flow through recording stubs (the workbench glue proper needs the whole
//    page, so we model its exact control flow here the way test_v2_3_4 section E
//    models the lens lifecycle).
// ===========================================================================

// A tiny recording element + document, just rich enough for the readout glue.
function makeEl(tag) {
    return {
        tagName: tag,
        className: '',
        type: '',
        title: '',
        textContent: '',
        attrs: {},
        _disabled: false,
        children: [],
        parentNode: null,
        _click: null,
        setAttribute: function (k, v) { this.attrs[k] = v; if (k === 'disabled') { this._disabled = true; } },
        getAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this.attrs, k) ? this.attrs[k] : null; },
        removeAttribute: function (k) { delete this.attrs[k]; if (k === 'disabled') { this._disabled = false; } },
        appendChild: function (c) { c.parentNode = this; this.children.push(c); return c; },
        removeChild: function (c) { var i = this.children.indexOf(c); if (i !== -1) { this.children.splice(i, 1); c.parentNode = null; } },
        addEventListener: function (ev, fn) { if (ev === 'click') { this._click = fn; } },
        classList: {
            _set: {},
            add: function (c) { this._set[c] = true; },
            remove: function (c) { delete this._set[c]; },
            toggle: function (c, on) { if (on) { this._set[c] = true; } else { delete this._set[c]; } },
            contains: function (c) { return !!this._set[c]; }
        },
        querySelector: function (sel) { return findBySelector(this, sel); }
    };
}

function findBySelector(root, sel) {
    // Only the class selectors the readout glue uses.
    var cls = sel.replace(/^\./, '');
    var stack = root.children.slice();
    while (stack.length) {
        var n = stack.shift();
        if ((' ' + n.className + ' ').indexOf(' ' + cls + ' ') !== -1) { return n; }
        if (n.children) { stack = stack.concat(n.children); }
    }
    return null;
}

// A faithful model of buildStructureLensReadout + updateLensReadout +
// removeStructureLensReadout, parameterised by a (smiles, warnings) provider so
// we can drive "edits" without a real editor/molecule.
function modelReadout(wrap, signalProvider) {
    function truncateMiddle(s, max) {
        s = String(s == null ? '' : s);
        if (s.length <= max) { return s; }
        var keep = max - 1;
        var head = Math.ceil(keep / 2);
        var tail = Math.floor(keep / 2);
        return s.slice(0, head) + '…' + s.slice(s.length - tail);
    }
    function remove() {
        var ex = wrap.querySelector('.wb-lens-readout');
        if (ex && ex.parentNode) { ex.parentNode.removeChild(ex); }
    }
    function update() {
        var bar = wrap.querySelector('.wb-lens-readout');
        if (!bar) { return; }
        var smEl = bar.querySelector('.wb-lens-readout-smiles');
        var wnEl = bar.querySelector('.wb-lens-readout-warn');
        var sig = signalProvider();
        var empty = !sig.atoms;
        if (smEl) {
            if (empty) { smEl.textContent = 'No structure yet'; smEl.setAttribute('disabled', 'disabled'); }
            else { smEl.textContent = sig.smiles ? truncateMiddle(sig.smiles, 48) : '(no SMILES)'; smEl.removeAttribute('disabled'); }
        }
        if (wnEl) {
            var n = sig.warnings.length;
            wnEl.textContent = n === 0 ? '0 warnings' : ('⚠ ' + n + ' warning' + (n === 1 ? '' : 's'));
            wnEl.classList.toggle('has-warnings', n > 0);
        }
    }
    function build() {
        remove();
        var bar = makeEl('div');
        bar.className = 'wb-lens-readout';
        bar.setAttribute('aria-live', 'polite');
        var sm = makeEl('button'); sm.type = 'button'; sm.className = 'wb-lens-readout-smiles'; sm.textContent = '—';
        var wn = makeEl('span'); wn.className = 'wb-lens-readout-warn'; wn.textContent = '0 warnings';
        bar.appendChild(sm); bar.appendChild(wn);
        wrap.appendChild(bar);
        update();
        return bar;
    }
    return { build: build, update: update, remove: remove };
}

test('D1 build creates an aria-live strip with a SMILES button + a warning badge', function () {
    var wrap = makeEl('div');
    var state = { atoms: 6, smiles: 'c1ccccc1', warnings: [] };
    var rd = modelReadout(wrap, function () { return state; });
    var bar = rd.build();
    assert.strictEqual(bar.getAttribute('aria-live'), 'polite', 'strip is aria-live=polite');
    var sm = wrap.querySelector('.wb-lens-readout-smiles');
    var wn = wrap.querySelector('.wb-lens-readout-warn');
    assert.ok(sm, 'SMILES element present');
    assert.ok(wn, 'warning element present');
    assert.strictEqual(sm.textContent, 'c1ccccc1', 'SMILES filled from the signal on build');
    assert.strictEqual(wn.textContent, '0 warnings', 'zero warnings reads "0 warnings"');
});

test('D2 update tracks edits: SMILES + warning count refresh, has-warnings toggles', function () {
    var wrap = makeEl('div');
    var state = { atoms: 6, smiles: 'c1ccccc1', warnings: [] };
    var rd = modelReadout(wrap, function () { return state; });
    rd.build();
    // Simulate the user drawing a structure that raises 2 warnings.
    state.smiles = 'CC(=O)O';
    state.warnings = ['a', 'b'];
    rd.update();
    var sm = wrap.querySelector('.wb-lens-readout-smiles');
    var wn = wrap.querySelector('.wb-lens-readout-warn');
    assert.strictEqual(sm.textContent, 'CC(=O)O', 'SMILES updated on edit');
    assert.strictEqual(wn.textContent, '⚠ 2 warnings', 'warning count updated with the alert glyph');
    assert.ok(wn.classList.contains('has-warnings'), 'has-warnings set when count > 0');
    // One-warning singular grammar.
    state.warnings = ['only-one'];
    rd.update();
    assert.strictEqual(wn.textContent, '⚠ 1 warning', 'singular "warning" for a count of 1');
});

test('D3 a long SMILES is truncated in the MIDDLE with an ellipsis', function () {
    var wrap = makeEl('div');
    var longSmiles = 'C'.repeat(80);
    var state = { atoms: 80, smiles: longSmiles, warnings: [] };
    var rd = modelReadout(wrap, function () { return state; });
    rd.build();
    var sm = wrap.querySelector('.wb-lens-readout-smiles');
    assert.ok(sm.textContent.length < longSmiles.length, 'long SMILES is shortened');
    assert.ok(sm.textContent.indexOf('…') !== -1, 'an ellipsis marks the elision');
});

test('D4 an empty editor reads "No structure yet" and disables the copy button', function () {
    var wrap = makeEl('div');
    var state = { atoms: 0, smiles: '', warnings: [] };
    var rd = modelReadout(wrap, function () { return state; });
    rd.build();
    var sm = wrap.querySelector('.wb-lens-readout-smiles');
    assert.strictEqual(sm.textContent, 'No structure yet', 'empty editor placeholder');
    assert.strictEqual(sm.getAttribute('disabled'), 'disabled', 'copy disabled with no structure');
});

test('D5 remove() tears the strip down (idempotent)', function () {
    var wrap = makeEl('div');
    var rd = modelReadout(wrap, function () { return { atoms: 6, smiles: 'C', warnings: [] }; });
    rd.build();
    assert.ok(wrap.querySelector('.wb-lens-readout'), 'strip present after build');
    rd.remove();
    assert.strictEqual(wrap.querySelector('.wb-lens-readout'), null, 'strip gone after remove');
    rd.remove(); // idempotent — no throw
    assert.strictEqual(wrap.querySelector('.wb-lens-readout'), null, 'still gone (idempotent)');
});

test('D6 rebuild does not double the strip (stale one removed first)', function () {
    var wrap = makeEl('div');
    var rd = modelReadout(wrap, function () { return { atoms: 6, smiles: 'C', warnings: [] }; });
    rd.build();
    rd.build();
    var n = 0;
    for (var i = 0; i < wrap.children.length; i++) {
        if ((' ' + wrap.children[i].className + ' ').indexOf(' wb-lens-readout ') !== -1) { n++; }
    }
    assert.strictEqual(n, 1, 'exactly one readout strip after a rebuild');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
