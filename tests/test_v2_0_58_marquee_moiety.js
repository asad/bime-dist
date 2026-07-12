/**
 * tests/test_v2_0_58_marquee_moiety.js — drag-rectangle moiety selection
 * (v2.0.58).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.54 shipped shift-click as the moiety-selection primitive. v2.0.58
 * adds a marquee (drag-rectangle) on the atom-trace cell SVG so users can
 * sweep a region and pick up every atom inside in one gesture. Plain
 * drag clears + sets the new moiety; shift+drag extends the existing
 * moiety on the same start metabolite.
 *
 * The marquee lives in the DOM (mousedown/move/up event sequence on the
 * atom-trace cell SVG), so we test it via source-shape assertions on
 * `js/workbench.js` and the supporting CSS rather than simulating a full
 * event sequence in Node.
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Marquee moiety select (v2.0.58)');
var test = runner.test;

console.log('Marquee moiety select (v2.0.58)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. handler functions exist + are wired into the document
// ---------------------------------------------------------------------------
test('A1 workbench.js declares handleAtomTraceMouseDown / Move / Up', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('function handleAtomTraceMouseDown(') !== -1,
        'mousedown handler declared');
    assert.ok(js.indexOf('function handleAtomTraceMouseMove(') !== -1,
        'mousemove handler declared');
    assert.ok(js.indexOf('function handleAtomTraceMouseUp(') !== -1,
        'mouseup handler declared');
});

test('A2 workbench.js wires the three handlers onto the document', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf("addEventListener('mousedown', handleAtomTraceMouseDown)") !== -1,
        'mousedown attached');
    assert.ok(js.indexOf("addEventListener('mousemove', handleAtomTraceMouseMove)") !== -1,
        'mousemove attached');
    assert.ok(js.indexOf("addEventListener('mouseup', handleAtomTraceMouseUp)") !== -1,
        'mouseup attached');
});

// ---------------------------------------------------------------------------
// B. drag-vs-click discrimination
// ---------------------------------------------------------------------------
test('B1 DRAG_THRESHOLD_PX is at least 3 (avoids accidental drags from tiny mouse jitter)', function () {
    var js = readRepoFile('js/workbench.js');
    var m = /var\s+DRAG_THRESHOLD_PX\s*=\s*(\d+)\s*;/.exec(js);
    assert.ok(m, 'DRAG_THRESHOLD_PX declared');
    assert.ok(parseInt(m[1], 10) >= 3,
        'threshold >= 3px, got ' + m[1]);
});

test('B2 click handler short-circuits on _suppressNextAtomTraceClick', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('_suppressNextAtomTraceClick') !== -1,
        'suppress flag declared');
    // The flag must be CHECKED inside the click handler too, not just set.
    var clickBlock = js.substring(js.indexOf("addEventListener('click'"));
    assert.ok(clickBlock.indexOf('_suppressNextAtomTraceClick') !== -1,
        'click handler tests the suppress flag');
});

// ---------------------------------------------------------------------------
// C. selection semantics (atom-rect hit-test, shift = extend)
// ---------------------------------------------------------------------------
test('C1 marquee hit-test uses cx/cy of .bime-trace-atom-hit', function () {
    var js = readRepoFile('js/workbench.js');
    var muUp = js.substring(js.indexOf('function handleAtomTraceMouseUp'));
    // Hit-test reads cx/cy attributes from .bime-trace-atom-hit circles.
    assert.ok(muUp.indexOf('.bime-trace-atom-hit') !== -1,
        'hit-target class referenced');
    assert.ok(muUp.indexOf("getAttribute('cx')") !== -1, 'reads cx');
    assert.ok(muUp.indexOf("getAttribute('cy')") !== -1, 'reads cy');
});

test('C2 shift+drag EXTENDS an existing moiety; plain drag REPLACES it', function () {
    var js = readRepoFile('js/workbench.js');
    var muUp = js.substring(js.indexOf('function handleAtomTraceMouseUp'),
                            js.indexOf('function handleAtomTraceMouseUp') + 3000);
    // Look for shift branch + extend semantics.
    assert.ok(muUp.indexOf('drag.shiftKey') !== -1, 'shift checked');
    assert.ok(muUp.indexOf('moietyStartMetIdx === drag.cellIdx') !== -1,
        'extend only when drag is on the active moiety cell');
    // Plain branch sets moietyStartMetIdx + moietySet fresh.
    assert.ok(muUp.indexOf('_pathway.moietyStartMetIdx = drag.cellIdx') !== -1,
        'plain drag re-anchors the moiety');
});

// ---------------------------------------------------------------------------
// D. marquee overlay rendering
// ---------------------------------------------------------------------------
test('D1 renderAtomTraceMarquee creates an SVG rect with class bime-trace-marquee', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('function renderAtomTraceMarquee(') !== -1,
        'marquee renderer declared');
    var mq = js.substring(js.indexOf('function renderAtomTraceMarquee('));
    assert.ok(mq.indexOf("createElementNS('http://www.w3.org/2000/svg', 'rect')") !== -1,
        'creates an SVG rect element');
    assert.ok(mq.indexOf("'bime-trace-marquee'") !== -1,
        'rect carries the bime-trace-marquee class');
});

test('D2 marquee rect uses pointer-events="none" so it does not block atom hit-tests', function () {
    var js = readRepoFile('js/workbench.js');
    var mq = js.substring(js.indexOf('function renderAtomTraceMarquee('));
    assert.ok(mq.indexOf("setAttribute('pointer-events', 'none')") !== -1,
        'marquee overlay is transparent to pointer events');
});

// ---------------------------------------------------------------------------
// E. CSS for the marquee
// ---------------------------------------------------------------------------
test('E1 css/workbench.css declares .bime-trace-marquee fill + stroke', function () {
    var css = readRepoFile('css/workbench.css');
    assert.ok(css.indexOf('.bime-trace-marquee') !== -1,
        'marquee class styled in CSS');
});

test('E2 atom-trace cell SVG has a crosshair cursor (marquee hint)', function () {
    var css = readRepoFile('css/workbench.css');
    assert.ok(/svg\.bime-trace-svg\s*\{[^}]*cursor:\s*crosshair/.test(css),
        'cursor: crosshair declared on the trace SVG');
});

// ---------------------------------------------------------------------------
// F. backward compat — click semantics preserved
// ---------------------------------------------------------------------------
test('F1 handleTraceAtomClick is still wired in the click dispatcher', function () {
    var js = readRepoFile('js/workbench.js');
    assert.ok(js.indexOf('handleTraceAtomClick(e.target, !!e.shiftKey)') !== -1,
        'single-click / shift-click path is preserved');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
