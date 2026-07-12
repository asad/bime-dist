/**
 * tests/test_v2_0_59_cofactor_lineage_ui.js — workbench wiring of v2.0.56
 * cofactor lineage into the atom-trace strip (v2.0.59).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.56 shipped the AtomTrace.parseCofactorString +
 * mapStepWithCofactors + tracePathway-cofactors plumbing but never
 * wired it into showAtomTraceStrip — so loading the glycolysis demo
 * still produced atom-trace strips with no cofactor-origin halos.
 *
 * v2.0.59 closes that wiring gap:
 *
 *  - showAtomTraceStrip parses each pathway-step `cofactors` string via
 *    AtomTrace.parseCofactorString into a per-edge spec array, and
 *    passes it to AtomTrace.tracePathway as options.cofactors.
 *  - The render loop reads result.cofactorLabels, lifts it into the
 *    full-strip frame, and applies an `is-cofactor-origin` class plus
 *    an SVG <title> tooltip naming the cofactor source for every atom
 *    whose ancestry includes a cofactor IN.
 *  - css/workbench.css carries the .bime-trace-atom.is-cofactor-origin
 *    rule (purple halo, layering below the moiety teal and primary
 *    single-trace amber).
 *
 * Plain Node, no DOM. Tests assert source shape rather than executing
 * the render path.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Cofactor lineage UI (v2.0.59)');
var test = runner.test;

console.log('Cofactor lineage UI (v2.0.59)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. js/workbench.js — parseCofactorString + traceCofactors wiring
// ---------------------------------------------------------------------------
test('A1 showAtomTraceStrip builds a traceCofactors array via parseCofactorString', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function showAtomTraceStrip');
    assert.ok(fnStart > 0, 'showAtomTraceStrip is defined');
    var slice = js.substring(fnStart, fnStart + 15000);
    assert.ok(slice.indexOf('AtomTrace.parseCofactorString') !== -1,
        'function calls AtomTrace.parseCofactorString');
    assert.ok(slice.indexOf('traceCofactors') !== -1,
        'function builds a traceCofactors array');
});

test('A2 traceOpts is extended with the cofactors spec', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function showAtomTraceStrip');
    var slice = js.substring(fnStart, fnStart + 15000);
    assert.ok(slice.indexOf('traceOpts.cofactors') !== -1,
        'traceOpts.cofactors is assigned from traceCofactors');
});

// ---------------------------------------------------------------------------
// B. cofactorLabelsByCell render
// ---------------------------------------------------------------------------
test('B1 showAtomTraceStrip lifts result.cofactorLabels into cofactorLabelsByCell', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function showAtomTraceStrip');
    var slice = js.substring(fnStart, fnStart + 15000);
    assert.ok(slice.indexOf('cofactorLabelsByCell') !== -1,
        'cofactorLabelsByCell variable is declared');
    assert.ok(slice.indexOf('result.cofactorLabels') !== -1,
        'reads result.cofactorLabels from tracePathway output');
});

test('B2 render loop applies is-cofactor-origin class to atoms with cofactor records', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function showAtomTraceStrip');
    var slice = js.substring(fnStart, fnStart + 15000);
    assert.ok(slice.indexOf("'is-cofactor-origin'") !== -1,
        'render loop adds the is-cofactor-origin class');
});

test('B3 render loop creates an SVG <title> tooltip naming the cofactor source', function () {
    var js = readRepoFile('js/workbench.js');
    var fnStart = js.indexOf('function showAtomTraceStrip');
    var slice = js.substring(fnStart, fnStart + 15000);
    // The title element is created via createElementNS('http://www.w3.org/2000/svg', 'title')
    assert.ok(slice.indexOf("createElementNS('http://www.w3.org/2000/svg', 'title')") !== -1,
        'SVG title element is constructed');
    assert.ok(slice.indexOf('Came from cofactor:') !== -1,
        'tooltip text includes "Came from cofactor:"');
});

// ---------------------------------------------------------------------------
// C. css/workbench.css — .is-cofactor-origin halo
// ---------------------------------------------------------------------------
test('C1 css/workbench.css declares the .is-cofactor-origin halo style', function () {
    var css = readRepoFile('css/workbench.css');
    assert.ok(css.indexOf('.bime-trace-atom.is-cofactor-origin') !== -1,
        '.bime-trace-atom.is-cofactor-origin rule is declared');
});

test('C2 cofactor halo defers to the primary single-trace amber and moiety teal', function () {
    var css = readRepoFile('css/workbench.css');
    // Layering: when an atom is BOTH cofactor-origin AND is-traced, the
    // amber single-trace halo wins; same for is-moiety-source -> teal.
    assert.ok(css.indexOf('.is-cofactor-origin.is-traced') !== -1,
        'layering rule for cofactor + single-trace present');
    assert.ok(css.indexOf('.is-cofactor-origin.is-moiety-source') !== -1,
        'layering rule for cofactor + moiety-source present');
});

// ---------------------------------------------------------------------------
// D. v2.0.56 primitives still exported (regression guard)
// ---------------------------------------------------------------------------
test('D1 AtomTrace exports parseCofactorString + COFACTOR_SMILES (v2.0.56 surface)', function () {
    require('../editor/AtomTrace.js');
    assert.strictEqual(typeof AtomTrace.parseCofactorString, 'function');
    assert.strictEqual(typeof AtomTrace.COFACTOR_SMILES, 'object');
    assert.ok(AtomTrace.COFACTOR_SMILES['atp'], 'ATP SMILES present');
});

test('D2 AtomTrace.version is at least 2.0.59 (forward-compat check)', function () {
    require('../editor/AtomTrace.js');
    assert.ok(typeof AtomTrace.version === 'string');
    var parts = AtomTrace.version.split('.').map(Number);
    var atLeast = (parts[0] > 2) || (parts[0] === 2 && (parts[1] > 0 || parts[2] >= 59));
    assert.ok(atLeast, 'AtomTrace.version >= 2.0.59, got ' + AtomTrace.version);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
