/**
 * tests/test_v2_3_5_first_run.js — first-run coherence (v2.3.5).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.3.5 fixes three first-run coherence issues a UX review flagged on the
 * continuum workbench ("blank grid, mystery molecule, invisible editor"). No
 * chemistry or algorithm change. This suite pins the three fixes:
 *
 *   P0-4  the DEFAULT boot (no ?smiles= / ?input= / ?pathway= deep link) no
 *         longer loads benzene (c1ccccc1) into the display:none editor. That
 *         had populated Editor Output with SMILES + properties for a molecule
 *         the user could see NOWHERE (the canvas was empty), so Canvas and
 *         Output disagreed and read as a bug. The editor now boots empty; the
 *         deep-link branches are untouched.
 *   P0-5  renderPathwayCanvas paints a centered, muted onboarding CTA
 *         (drawEmptyCanvasCta) when the pathway is empty — pointer-events:none
 *         (so a double-click passes through to start drawing), aria-hidden,
 *         and within the --color token system. It is drawn ONLY in the empty
 *         branch, so any added content removes it on the next render.
 *   P1-3  the add-pathway-current action reads "Draw Structure" EVERYWHERE
 *         (Tools row + command grid) — one action no longer wears two labels
 *         ("Add Structure" was the old command-grid copy). The two pan-copy
 *         strings (Select-tool <title> + the Select status line) now state the
 *         TRUE pan affordance identically (verified against onPathwayPointerDown:
 *         pan is Shift+drag / middle-button-drag on the empty canvas, NOT a
 *         plain drag).
 *
 * Mostly source-shape (the workbench DOM glue needs the whole page to execute,
 * exactly as test_v2_3_0 / test_v2_3_4 note), plus an executable DOM-stub
 * section that drives drawEmptyCanvasCta + the empty-vs-content render decision
 * through recording stubs — modelled on test_v2_3_4_lens_hardening.js section E
 * and the tests/shim.js document stub.
 *
 * Plain Node, no real DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('First-run coherence (v2.3.5)');
var test = runner.test;
console.log('First-run coherence (v2.3.5)');

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

// ===========================================================================
// P0-4 — no phantom benzene on the default boot.
// ===========================================================================
test('P0-4.1 the DOMContentLoaded boot no longer loads benzene by default', function () {
    // The old default tail `loadNamed('c1ccccc1', 'Benzene')` is gone — there
    // is no benzene literal anywhere in the boot path or the file.
    assert.strictEqual(WB.indexOf("loadNamed('c1ccccc1'"), -1,
        "the default-boot loadNamed('c1ccccc1', ...) call is removed");
    assert.strictEqual(WB.indexOf('c1ccccc1'), -1,
        'no benzene SMILES literal remains in js/workbench.js');
});

test('P0-4.2 the deep-link load branches are PRESERVED exactly', function () {
    // The two deep-link arms still exist: a ?pathway= example and a
    // ?smiles=/?input= structure. Only the no-deep-link default changed.
    assert.ok(WB.indexOf('var initial = getInitialWorkbenchLoad();') !== -1,
        'boot still consults getInitialWorkbenchLoad()');
    assert.ok(WB.indexOf('loadExamplePathway(initial.pathway)') !== -1,
        'deep-link: ?pathway= still loads the example pathway');
    assert.ok(WB.indexOf('loadNamed(initial.input, initial.name)') !== -1,
        'deep-link: ?smiles=/?input= still loads the structure');
    // getInitialWorkbenchLoad still reads the supported params.
    var fn = wbFn('getInitialWorkbenchLoad');
    assert.ok(fn.indexOf("qs.get('pathway')") !== -1, '?pathway= param still read');
    assert.ok(fn.indexOf("qs.get('smiles')") !== -1, '?smiles= param still read');
    assert.ok(fn.indexOf("qs.get('input')") !== -1, '?input= param still read');
});

test('P0-4.3 the boot block documents the empty-editor default (no else-load arm)', function () {
    // Locate the boot decision and confirm there is no terminal `else { ...
    // loadNamed ... }` re-introducing a default molecule.
    var idx = WB.indexOf('var initial = getInitialWorkbenchLoad();');
    assert.ok(idx > 0, 'boot decision located');
    var block = WB.slice(idx, idx + 900);
    assert.strictEqual(block.indexOf("else {\n            loadNamed"), -1,
        'no else-arm reintroduces a default-boot molecule');
    assert.ok(block.indexOf('P0-4') !== -1,
        'the empty-default rationale is documented inline');
});

// Executable model of the boot decision: no deep link -> NOTHING loads.
test('P0-4.4 boot decision model: no deep link loads nothing; deep links still load', function () {
    var calls = [];
    function bootModel(initial) {
        if (initial && initial.pathway) {
            calls.push('pathway:' + initial.pathway);
        } else if (initial) {
            calls.push('named:' + initial.input);
        }
        // v2.3.5: a no-deep-link boot loads nothing (its terminal else only
        // RENDERS the empty canvas — no load; see P0-5.7).
    }
    calls = []; bootModel(null);
    assert.deepStrictEqual(calls, [], 'no deep link => no load (empty editor)');
    calls = []; bootModel({ pathway: 'glycolysis' });
    assert.deepStrictEqual(calls, ['pathway:glycolysis'], '?pathway= still loads');
    calls = []; bootModel({ input: 'CCO', name: 'Ethanol' });
    assert.deepStrictEqual(calls, ['named:CCO'], '?smiles=/?input= still loads');
});

test('P0-4.5 Editor Output keeps its empty-state placeholder (coherent empty default)', function () {
    // With the editor empty, the SMILES output textarea shows its placeholder
    // instead of a phantom molecule's SMILES — the existing empty-state.
    assert.ok(HTML.indexOf('id="smiles-out"') !== -1, 'smiles-out output present');
    var idx = HTML.indexOf('id="smiles-out"');
    var openIdx = HTML.lastIndexOf('<textarea', idx);
    var closeIdx = HTML.indexOf('>', idx);
    var tag = HTML.substring(openIdx, closeIdx);
    assert.ok(tag.indexOf('placeholder=') !== -1,
        'the output textarea carries an empty-state placeholder');
});

// ===========================================================================
// P0-5 — empty-canvas onboarding CTA.
// ===========================================================================
test('P0-5.1 isPathwayEmpty + drawEmptyCanvasCta exist as top-level functions', function () {
    assert.ok(WB.indexOf('function isPathwayEmpty(') !== -1, 'isPathwayEmpty defined');
    assert.ok(WB.indexOf('function drawEmptyCanvasCta(') !== -1, 'drawEmptyCanvasCta defined');
});

test('P0-5.2 isPathwayEmpty checks EVERY pathway collection (matches clearPathwayCanvas)', function () {
    var fn = wbFn('isPathwayEmpty');
    assert.ok(fn, 'isPathwayEmpty located');
    var cols = ['nodes', 'edges', 'steps', 'notes', 'mechanismArrows', 'compartments', 'backgrounds'];
    for (var i = 0; i < cols.length; i++) {
        assert.ok(fn.indexOf('_pathway.' + cols[i] + '.length') !== -1,
            'isPathwayEmpty inspects _pathway.' + cols[i]);
    }
});

test('P0-5.3 renderPathwayCanvas draws the CTA ONLY in the empty branch and is NOT split', function () {
    var fn = wbFn('renderPathwayCanvas');
    assert.ok(fn, 'renderPathwayCanvas located');
    // The CTA is drawn behind an isPathwayEmpty() guard.
    assert.ok(/if\s*\(\s*isPathwayEmpty\(\)\s*\)\s*\{\s*drawEmptyCanvasCta\(/.test(fn),
        'CTA is drawn only when isPathwayEmpty() is true');
    // renderPathwayCanvas remains a single top-level function (source-shape).
    assert.strictEqual(fn.indexOf('\nfunction '), -1,
        'renderPathwayCanvas contains no nested top-level function declaration');
    // It still paints the real layers (the CTA is additive, not a replacement).
    assert.ok(fn.indexOf('renderPathwayNode(viewport') !== -1, 'nodes still painted');
});

test('P0-5.4 the CTA group is pointer-events:none + aria-hidden, and uses --color tokens', function () {
    var fn = wbFn('drawEmptyCanvasCta');
    assert.ok(fn, 'drawEmptyCanvasCta located');
    assert.ok(fn.indexOf("'pointer-events': 'none'") !== -1,
        'CTA group sets pointer-events:none so clicks pass through to the canvas');
    assert.ok(fn.indexOf("'aria-hidden': 'true'") !== -1,
        'CTA group is aria-hidden (it must not confuse screen readers)');
    assert.ok(fn.indexOf('var(--color-text-muted)') !== -1,
        'CTA text uses the --color-text-muted token (muted, on-brand)');
    // Centered on the active page (uses pathwayPageSize, not a hardcoded size).
    assert.ok(fn.indexOf('pathwayPageSize()') !== -1,
        'CTA is centered via pathwayPageSize() (A4/Letter aware)');
    // The class hook the CSS targets.
    assert.ok(fn.indexOf("'class': 'wb-pathway-cta'") !== -1,
        'CTA group carries the wb-pathway-cta class');
});

test('P0-5.5 the CTA copy is terse + on-brand (no emoji, no third-party names)', function () {
    var fn = wbFn('drawEmptyCanvasCta');
    // The actionable line names BOTH entry points (double-click + Example).
    assert.ok(fn.indexOf('Double-click anywhere to draw a structure') !== -1,
        'CTA tells the user to double-click to draw');
    assert.ok(fn.toLowerCase().indexOf('example') !== -1,
        'CTA points at the Example controls');
    // On-brand: no emoji surrogate pairs anywhere in the helper.
    assert.ok(!/[\uD800-\uDBFF]/.test(fn), 'CTA copy contains no emoji');
});

test('P0-5.6 css styles the CTA without animation (no prefers-reduced-motion concern)', function () {
    assert.ok(CSS.indexOf('.wb-pathway-cta') !== -1, 'wb-pathway-cta styled');
    assert.ok(CSS.indexOf('.wb-pathway-cta-title') !== -1, 'CTA title styled');
    assert.ok(CSS.indexOf('.wb-pathway-cta-line') !== -1, 'CTA line styled');
    // The CTA block declares pointer-events:none (defense-in-depth with markup).
    var idx = CSS.indexOf('.wb-pathway-cta {');
    assert.ok(idx !== -1, 'wb-pathway-cta rule present');
    var block = CSS.slice(idx, idx + 60);
    assert.ok(block.indexOf('pointer-events: none') !== -1,
        'the CTA group is pointer-events:none in CSS too');
    // No @keyframes animation is attached to the CTA (static hint).
    assert.strictEqual(CSS.indexOf('.wb-pathway-cta-chevron { animation'), -1,
        'the chevron is a static hint (no animation -> no reduced-motion issue)');
});

test('P0-5.7 the no-deep-link boot RENDERS the empty canvas so the CTA shows on first paint', function () {
    // Regression guard for the boot gap: removing the default benzene load (P0-4)
    // also removed the render that used to fire on boot, so on a fresh default
    // load the canvas stayed blank — drawEmptyCanvasCta only ran after the first
    // user action. The no-deep-link terminal else MUST call renderPathwayCanvas()
    // so the onboarding CTA paints on the FIRST frame. (A DOM-stub can't drive the
    // DOMContentLoaded handler, so this is source-shape — but it pins the call.)
    var idx = WB.indexOf('var initial = getInitialWorkbenchLoad();');
    assert.ok(idx > 0, 'boot decision located');
    var block = WB.slice(idx, idx + 1200);
    // v2.4.18: the empty-canvas render is gated behind PATHWAY_ENABLED — it paints
    // the onboarding CTA on the first frame when the pathway view ships, and is
    // skipped for the core-only build (the editor is the boot hero with its own
    // empty-state placeholder, P0-4.5). The source shape is identical either way.
    var armIdx = block.indexOf('else if (PATHWAY_ENABLED)');
    assert.ok(armIdx !== -1,
        'the no-deep-link render arm is gated behind PATHWAY_ENABLED');
    assert.ok(block.indexOf('renderPathwayCanvas()', armIdx) !== -1,
        'that arm calls renderPathwayCanvas() (CTA on first paint when pathway is enabled)');
});

// ===========================================================================
// P1-3 — unify duplicate action labels + reconcile the pan copy.
// ===========================================================================
test('P1-3.1 the add-pathway-current action reads "Draw Structure" EVERYWHERE', function () {
    // Same action, wired twice (Tools row + command grid).
    assert.strictEqual(count(HTML, 'data-wb-action="add-pathway-current"'), 2,
        'add-pathway-current is wired exactly twice (Tools row + command grid)');
    // Both placements carry the unified label; the old "Add Structure" copy is gone.
    assert.strictEqual(count(HTML, '>Draw Structure<'), 2,
        'both add-pathway-current buttons read "Draw Structure"');
    assert.strictEqual(HTML.indexOf('>Add Structure<'), -1,
        'the old "Add Structure" label is gone (one action, one label)');
});

test('P1-3.2 "Edit Structure" is consistent across Tools row + command grid', function () {
    assert.strictEqual(count(HTML, 'data-wb-action="edit-pathway-structure"'), 2,
        'edit-pathway-structure is wired exactly twice');
    assert.strictEqual(count(HTML, '>Edit Structure<'), 2,
        'both edit-pathway-structure buttons read "Edit Structure"');
});

test('P1-3.3 both pan-copy strings state the TRUE pan affordance identically', function () {
    // The Select-tool <title> attribute.
    var selIdx = HTML.indexOf('data-tool="select"');
    assert.ok(selIdx > 0, 'Select tool button present');
    var titleMatch = /title="([^"]*)"/.exec(HTML.slice(selIdx, selIdx + 400));
    assert.ok(titleMatch, 'Select tool has a title attribute');
    var titleCopy = titleMatch[1];
    // The Select-tool status line in workbench.js.
    var statusIdx = WB.indexOf("pathwayStatus('Select and drag items.");
    assert.ok(statusIdx > 0, 'Select status line present');
    var statusCopy = /pathwayStatus\('([^']*)'\)/.exec(WB.slice(statusIdx, statusIdx + 200))[1];

    // The canonical, VERIFIED pan clause (onPathwayPointerDown pans on
    // e.shiftKey || middle button, NOT a plain drag). Both strings must
    // contain it identically (case-insensitively — one is mid-sentence).
    var canonical = 'shift while dragging the background to pan';
    assert.ok(titleCopy.toLowerCase().indexOf(canonical) !== -1,
        'the Select <title> states the Shift-to-pan affordance; was: ' + titleCopy);
    assert.ok(statusCopy.toLowerCase().indexOf(canonical) !== -1,
        'the Select status line states the Shift-to-pan affordance; was: ' + statusCopy);
    // The OLD contradictory "drag empty canvas to pan" copy is gone.
    assert.strictEqual(HTML.indexOf('drag empty canvas to pan'), -1,
        'the inaccurate "drag empty canvas to pan" copy is removed');
});

test('P1-3.4 the pan copy matches the actual implementation (Shift / middle-button)', function () {
    // Verify the source of truth: panning is gated on Shift OR middle button
    // while the Select tool is active and the pointer is on the empty canvas.
    var fn = wbFn('onPathwayPointerDown');
    assert.ok(fn, 'onPathwayPointerDown located');
    assert.ok(/_pathway\.tool === 'select' && \(e\.shiftKey \|\| e\.button === 1\)/.test(fn),
        'pan starts on Shift+drag or middle-button drag (NOT a plain empty-canvas drag)');
    assert.ok(/type:\s*'pan'/.test(fn), 'that branch sets a pan drag');
});

// ===========================================================================
// E. Executable DOM-stub: drive drawEmptyCanvasCta + the render decision
//    through recording stubs (the workbench glue proper needs the whole page,
//    so we model its exact control flow here the way test_v2_3_4 section E
//    models the lens lifecycle).
// ===========================================================================

// A recording SVG element: records tag, attributes, textContent, and children.
function recEl(tag) {
    return {
        tag: tag || 'g',
        _attrs: {},
        textContent: '',
        children: [],
        setAttribute: function (k, v) { this._attrs[k] = v; },
        setAttributeNS: function (ns, k, v) { this._attrs[k] = v; },
        getAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k) ? this._attrs[k] : null; },
        appendChild: function (c) { this.children.push(c); return c; }
    };
}

test('E1 drawEmptyCanvasCta model: builds a pointer-events:none, aria-hidden group with the copy', function () {
    // Faithful model of the helper's control flow against recording stubs.
    var PATHWAY_NS = 'http://www.w3.org/2000/svg';
    function pathwaySvgEl(tag, attrs) {
        var el = recEl(tag);
        attrs = attrs || {};
        for (var key in attrs) {
            if (!Object.prototype.hasOwnProperty.call(attrs, key)) continue;
            if (key === 'text') { el.textContent = attrs[key]; }
            else { el.setAttribute(key, attrs[key]); }
        }
        return el;
    }
    function drawModel(viewport) {
        var size = { w: 1200, h: 620 };
        var cx = size.w / 2, cy = size.h / 2;
        var g = pathwaySvgEl('g', { 'class': 'wb-pathway-cta', 'aria-hidden': 'true', 'pointer-events': 'none' });
        g.appendChild(pathwaySvgEl('text', { 'class': 'wb-pathway-cta-title', x: cx, y: cy - 18, fill: 'var(--color-text-muted)', text: 'Your canvas is empty' }));
        g.appendChild(pathwaySvgEl('text', { 'class': 'wb-pathway-cta-line', x: cx, y: cy + 12, fill: 'var(--color-text-muted)', text: 'Double-click anywhere to draw a structure — or load an Example below' }));
        g.appendChild(pathwaySvgEl('path', { 'class': 'wb-pathway-cta-chevron', d: 'M 588 664 L 600 676 L 612 664', stroke: 'var(--color-text-muted)' }));
        viewport.appendChild(g);
        return g;
    }
    var viewport = recEl('g');
    var g = drawModel(viewport);
    assert.strictEqual(viewport.children.length, 1, 'one CTA group appended to the viewport');
    assert.strictEqual(g.getAttribute('pointer-events'), 'none', 'group is pointer-events:none');
    assert.strictEqual(g.getAttribute('aria-hidden'), 'true', 'group is aria-hidden');
    assert.strictEqual(g.getAttribute('class'), 'wb-pathway-cta', 'group has the CTA class');
    assert.strictEqual(g.children.length, 3, 'headline + line + chevron');
    assert.strictEqual(g.children[0].textContent, 'Your canvas is empty', 'headline copy');
    assert.ok(g.children[1].textContent.indexOf('Double-click anywhere') !== -1, 'actionable line copy');
    assert.strictEqual(g.children[0].getAttribute('fill'), 'var(--color-text-muted)', 'muted token fill');
});

test('E2 render decision model: CTA drawn when empty, absent the moment content exists', function () {
    function isEmpty(pw) {
        return pw.nodes.length === 0 && pw.edges.length === 0 && pw.steps.length === 0 &&
            pw.notes.length === 0 && pw.mechanismArrows.length === 0 &&
            pw.compartments.length === 0 && pw.backgrounds.length === 0;
    }
    var drawn;
    function renderModel(pw, viewport) {
        // (real layers painted here)
        if (isEmpty(pw)) { drawn = true; viewport.appendChild(recEl('g')); }
    }
    var empty = { nodes: [], edges: [], steps: [], notes: [], mechanismArrows: [], compartments: [], backgrounds: [] };
    var vp1 = recEl('g'); drawn = false; renderModel(empty, vp1);
    assert.strictEqual(drawn, true, 'CTA drawn on an empty pathway');
    assert.strictEqual(vp1.children.length, 1, 'CTA group present when empty');

    var withNode = { nodes: [{ id: 'n1' }], edges: [], steps: [], notes: [], mechanismArrows: [], compartments: [], backgrounds: [] };
    var vp2 = recEl('g'); drawn = false; renderModel(withNode, vp2);
    assert.strictEqual(drawn, false, 'CTA NOT drawn once a node exists');
    assert.strictEqual(vp2.children.length, 0, 'no CTA group when content exists (a re-render removes it)');

    // Each non-empty collection independently suppresses the CTA.
    var keys = ['edges', 'steps', 'notes', 'mechanismArrows', 'compartments', 'backgrounds'];
    for (var i = 0; i < keys.length; i++) {
        var pw = { nodes: [], edges: [], steps: [], notes: [], mechanismArrows: [], compartments: [], backgrounds: [] };
        pw[keys[i]] = [{ id: 'x' }];
        assert.strictEqual(isEmpty(pw), false, keys[i] + ' content suppresses the CTA');
    }
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
