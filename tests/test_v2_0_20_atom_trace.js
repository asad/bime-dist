/**
 * tests/test_v2_0_20_atom_trace.js — cross-pathway atom tracing (v2.0.20).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Verifies AtomTrace chains per-reaction RDT atom–atom maps into per-atom paths
 * across a pathway. The scientific check: trace the carbons of D-glucose
 * through the (open, repertoire-derived) glycolysis backbone — all six reach
 * fructose-1,6-bisphosphate, then the aldolase cleavage forks the C6 skeleton
 * into two trioses so only the three carbons that flow into glyceraldehyde-3-P
 * continue down to pyruvate.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');
require('../editor/MetaboliteLibrary.js');
require('../editor/AtomTrace.js');

var runner = shim.makeRunner('Atom Trace (v2.0.20)');
var test = runner.test;

console.log('Atom Trace (v2.0.20)');

// linear glycolysis backbone (the aldolase split's other branch, DHAP, is a
// separate step and is intentionally not in this linear chain).
var GLYCOLYSIS = ['D-Glucose', 'Glucose 6-phosphate', 'Fructose 6-phosphate',
    'Fructose 1,6-bisphosphate', 'Glyceraldehyde 3-phosphate', '1,3-Bisphosphoglycerate',
    '3-Phosphoglycerate', '2-Phosphoglycerate', 'Phosphoenolpyruvate', 'Pyruvate'
].map(function(n) { return MetaboliteLibrary.find(n).smiles; });
var ALDOLASE_PRODUCT_INDEX = 4; // GAP — first metabolite after the F1,6BP cleavage
var PYRUVATE_INDEX = 9;

// ---------------------------------------------------------------------------
// A. mapStep
// ---------------------------------------------------------------------------
test('A1 mapStep maps a simple oxidation backbone (CCO -> CC=O)', function() {
    var step = AtomTrace.mapStep('CCO', 'CC=O');
    assert.strictEqual(step.mappedCount, 3, 'all 3 heavy atoms map');
    // index map covers reactant indices 0,1,2
    assert.ok(step.fromToIndex[0] !== undefined && step.fromToIndex[1] !== undefined &&
              step.fromToIndex[2] !== undefined, 'every reactant atom has an image');
});

// ---------------------------------------------------------------------------
// B. tiny chain
// ---------------------------------------------------------------------------
test('B1 two-step chain produces a path of length 3 for a conserved atom', function() {
    var r = AtomTrace.tracePathway(['CCO', 'CC=O', 'CC(=O)O']);
    assert.strictEqual(r.metaboliteCount, 3);
    assert.strictEqual(r.startAtoms, 3);
    // at least one carbon should be conserved through both steps
    var full = r.traces.filter(function(t) { return t.element === 'C' && t.path.length === 3; });
    assert.ok(full.length >= 1, 'a carbon survives both steps');
});

// ---------------------------------------------------------------------------
// C. glycolysis carbon fate (the headline trace)
// ---------------------------------------------------------------------------
var GLY = AtomTrace.tracePathway(GLYCOLYSIS);
var glucoseCarbons = GLY.traces.filter(function(t) { return t.element === 'C'; });

test('C1 D-glucose contributes exactly 6 carbon traces', function() {
    assert.strictEqual(glucoseCarbons.length, 6);
});

test('C2 all 6 glucose carbons survive to fructose-1,6-bisphosphate (pre-aldolase)', function() {
    var reachF16BP = glucoseCarbons.filter(function(t) {
        return t.path.some(function(p) { return p.metabolite === 3; });
    });
    assert.strictEqual(reachF16BP.length, 6, 'hexose phase conserves all carbons');
});

test('C3 the aldolase cleavage forks the skeleton — only 3 carbons enter GAP', function() {
    var reachGAP = glucoseCarbons.filter(function(t) {
        return t.path.some(function(p) { return p.metabolite === ALDOLASE_PRODUCT_INDEX; });
    });
    assert.strictEqual(reachGAP.length, 3, 'half the carbons flow into glyceraldehyde-3-P');
    // the other 3 are lost exactly at the aldolase step
    var lostAtAldolase = glucoseCarbons.filter(function(t) { return t.lostAtStep === ALDOLASE_PRODUCT_INDEX; });
    assert.strictEqual(lostAtAldolase.length, 3, 'the other three carbons fork off at aldolase');
});

test('C4 three glucose carbons trace all the way to pyruvate', function() {
    var reachPyruvate = glucoseCarbons.filter(function(t) {
        return t.path.some(function(p) { return p.metabolite === PYRUVATE_INDEX && p.element === 'C'; });
    });
    assert.strictEqual(reachPyruvate.length, 3, 'three carbons of glucose become a pyruvate');
    reachPyruvate.forEach(function(t) {
        assert.strictEqual(t.reachesEnd, true, 'a through-carbon reaches the final metabolite');
    });
});

test('C5 every step of the backbone atom-maps (status mapped)', function() {
    GLY.steps.forEach(function(s) {
        assert.strictEqual(s.status, 'mapped', 'step ' + s.from + '->' + s.to + ' maps');
        assert.ok(s.mapped > 0, 'step ' + s.from + '->' + s.to + ' has mapped atoms');
    });
});

// ---------------------------------------------------------------------------
// D. renderInteractiveSvg — interactive strip-renderer with addressable atoms
// ---------------------------------------------------------------------------
test('D1 renderInteractiveSvg emits one addressable group per heavy atom', function() {
    var m = new Molecule();
    SmilesParser.parse('CCO', m);
    var svg = AtomTrace.renderInteractiveSvg(m, { width: 160, height: 90 });
    assert.ok(svg.indexOf('<svg') === 0, 'starts with <svg');
    var groups = (svg.match(/class="bime-trace-atom"/g) || []).length;
    assert.strictEqual(groups, 3, 'one trace-atom group per atom');
    assert.ok(svg.indexOf('data-atom-index="0"') > 0, 'first atom indexed');
    assert.ok(svg.indexOf('data-atom-index="2"') > 0, 'last atom indexed');
});

test('D2 renderInteractiveSvg handles an empty molecule without throwing', function() {
    var svg = AtomTrace.renderInteractiveSvg(new Molecule());
    assert.ok(svg.indexOf('<svg') === 0, 'still emits a valid svg root');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
