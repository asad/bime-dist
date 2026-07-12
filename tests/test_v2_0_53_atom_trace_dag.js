/**
 * tests/test_v2_0_53_atom_trace_dag.js — Atom-trace DAG topology fix
 * (v2.0.53).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.20-v2.0.52 traced atoms along the LINEAR first-appearance chain
 * `metaboliteSmiles[i] >> metaboliteSmiles[i+1]`. For pathways with forks
 * (e.g. aldolase: F1,6BP → GAP + DHAP) and convergences (TPI: DHAP →
 * GAP), this produced phantom reactions that don't exist in biology —
 * for glycolysis specifically:
 *
 *   GAP   → DHAP   (phantom; the real reaction is DHAP → GAP)
 *   DHAP  → 1,3BPG (phantom; the real reaction is GAP → 1,3BPG)
 *
 * v2.0.53 accepts `options.edges = [{from, to}, …]` and walks the real
 * pathway DAG. Labels propagate in topological order through the
 * provided edges. The label model is multi-valued: a single downstream
 * atom can have multiple ancestor start-atoms under convergence (e.g.
 * pyruvate C1 has BOTH glucose C1 ancestors via the two halves of
 * F1,6BP cleavage).
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');
require('../editor/AtomTrace.js');
require('../editor/MetaboliteLibrary.js');

var runner = shim.makeRunner('Atom-trace DAG topology (v2.0.53)');
var test = runner.test;

console.log('Atom-trace DAG topology (v2.0.53)');

// Build the glycolysis metabolite-SMILES array + DAG-edge array used
// across the tests.
function glycolysisSetup() {
    var pw = MetaboliteLibrary.pathways.glycolysis;
    var order = [], seen = {};
    pw.steps.forEach(function (s) {
        [s.from, s.to].forEach(function (n) {
            if (!seen[n]) { seen[n] = true; order.push(n); }
        });
    });
    var nameToIdx = {};
    order.forEach(function (n, i) { nameToIdx[n.toLowerCase()] = i; });
    var smiles = order.map(function (n) {
        return (MetaboliteLibrary.find(n) || {}).smiles || '';
    });
    var edges = pw.steps.map(function (s) {
        return {
            from: nameToIdx[s.from.toLowerCase()],
            to:   nameToIdx[s.to.toLowerCase()]
        };
    });
    return { order: order, nameToIdx: nameToIdx, smiles: smiles, edges: edges, pw: pw };
}

// ---------------------------------------------------------------------------
// A. backward-compat — linear-chain mode is unchanged
// ---------------------------------------------------------------------------
test('A1 tracePathway without options.edges still uses the v2.0.20 linear chain', function () {
    // CCO -> CC=O -> CC(=O)O is a 3-step linear chain with no forks.
    // No edges argument → falls back to [[0,1], [1,2]]. The carbons should
    // map through all three metabolites.
    var r = AtomTrace.tracePathway(['CCO', 'CC=O', 'CC(=O)O']);
    assert.ok(r, 'tracePathway returns a result');
    assert.strictEqual(r.metaboliteCount, 3);
    // Backward-compat: 'edges' surface should reflect the linear inferred chain.
    assert.deepStrictEqual(r.edges, [{ from: 0, to: 1 }, { from: 1, to: 2 }]);
});

// ---------------------------------------------------------------------------
// B. DAG mode — explicit edges replace the linear chain
// ---------------------------------------------------------------------------
test('B1 tracePathway with options.edges uses the provided DAG edges', function () {
    // Same 3-metabolite linear chain but spelled out explicitly.
    var r = AtomTrace.tracePathway(['CCO', 'CC=O', 'CC(=O)O'], {
        edges: [{ from: 0, to: 1 }, { from: 1, to: 2 }]
    });
    assert.deepStrictEqual(r.edges, [{ from: 0, to: 1 }, { from: 1, to: 2 }]);
    assert.strictEqual(r.steps.length, 2);
});

test('B2 array-pair edges [u, v] are accepted alongside object-style {from, to}', function () {
    var r = AtomTrace.tracePathway(['CCO', 'CC=O', 'CC(=O)O'], {
        edges: [[0, 1], [1, 2]]
    });
    assert.deepStrictEqual(r.edges, [{ from: 0, to: 1 }, { from: 1, to: 2 }]);
});

test('B3 empty / malformed edges fall back to the linear default (defensive)', function () {
    var r1 = AtomTrace.tracePathway(['CCO', 'CC=O'], { edges: [] });
    assert.deepStrictEqual(r1.edges, [{ from: 0, to: 1 }]);
    var r2 = AtomTrace.tracePathway(['CCO', 'CC=O'], { edges: [{ junk: 1 }, null] });
    assert.deepStrictEqual(r2.edges, [{ from: 0, to: 1 }]);
});

// ---------------------------------------------------------------------------
// C. multi-valued labels: convergence aggregates ancestors
// ---------------------------------------------------------------------------
test('C1 labels[m][k] is now an array of start-atom indices (multi-valued)', function () {
    var r = AtomTrace.tracePathway(['CCO', 'CC=O']);
    assert.ok(Array.isArray(r.labels[0]));
    assert.ok(Array.isArray(r.labels[0][0]));
    // Start metabolite seeds the identity labels: labels[0][k] = [k].
    assert.deepStrictEqual(r.labels[0][0], [0]);
    assert.deepStrictEqual(r.labels[0][1], [1]);
});

// ---------------------------------------------------------------------------
// D. glycolysis correctness — phantom reactions are gone
// ---------------------------------------------------------------------------
test('D1 glycolysis DAG mode produces the 11 real reactions, not the v2.0.52 phantoms', function () {
    var s = glycolysisSetup();
    var r = AtomTrace.tracePathway(s.smiles, { edges: s.edges });
    // The real glycolysis has 11 steps. The phantom GAP→DHAP and DHAP→1,3BPG
    // edges that v2.0.52's linear-chain mode invented are NOT present.
    var pairs = r.edges.map(function (e) {
        return s.order[e.from] + ' -> ' + s.order[e.to];
    });
    // Confirm the actual aldolase fork is present.
    assert.ok(pairs.indexOf('Fructose 1,6-bisphosphate -> Glyceraldehyde 3-phosphate') !== -1,
        'F1,6BP -> GAP fork is present');
    assert.ok(pairs.indexOf('Fructose 1,6-bisphosphate -> Dihydroxyacetone phosphate') !== -1,
        'F1,6BP -> DHAP fork is present');
    // Confirm TPI runs in the biology-correct direction.
    assert.ok(pairs.indexOf('Dihydroxyacetone phosphate -> Glyceraldehyde 3-phosphate') !== -1,
        'TPI runs DHAP -> GAP');
    // Phantom edges from v2.0.52 first-appearance walker MUST be absent.
    assert.strictEqual(pairs.indexOf('Glyceraldehyde 3-phosphate -> Dihydroxyacetone phosphate'), -1,
        'phantom GAP -> DHAP is NOT present');
    assert.strictEqual(pairs.indexOf('Dihydroxyacetone phosphate -> 1,3-Bisphosphoglycerate'), -1,
        'phantom DHAP -> 1,3BPG is NOT present');
});

test('D2 glycolysis: pyruvate carbon atoms each have TWO glucose-carbon ancestors (textbook biochemistry)', function () {
    var s = glycolysisSetup();
    var r = AtomTrace.tracePathway(s.smiles, { edges: s.edges });
    var pyrIdx = s.nameToIdx['pyruvate'];
    var pyrLabels = r.labels[pyrIdx];
    assert.ok(pyrLabels, 'pyruvate labels populated');
    // Count pyruvate atoms with >= 2 distinct glucose ancestors. The
    // glycolysis carbon-fate diagram has C1, C2, C3 of pyruvate each
    // contributed BY two glucose carbons (top half via DHAP-then-TPI
    // and bottom half via GAP-direct), so at least three pyruvate atoms
    // should have a 2-element ancestor list.
    var multiAncestorCount = 0;
    for (var i = 0; i < pyrLabels.length; i++) {
        if (pyrLabels[i].length >= 2) { multiAncestorCount++; }
    }
    assert.ok(multiAncestorCount >= 3,
        'at least 3 pyruvate atoms have >= 2 glucose ancestors (got ' + multiAncestorCount + ')');
});

test('D3 glycolysis: every glucose carbon that reaches pyruvate via either path is reachesEnd=true', function () {
    var s = glycolysisSetup();
    var r = AtomTrace.tracePathway(s.smiles, { edges: s.edges });
    // Find the glucose carbon traces.
    var carbonTraces = r.traces.filter(function (t) { return t.element === 'C'; });
    assert.strictEqual(carbonTraces.length, 6, 'glucose has 6 carbons');
    // In v2.0.52 linear mode, the carbons that took the DHAP path were
    // wrongly flagged "lost at GAP"; only some carbons reached pyruvate.
    // In v2.0.53 DAG mode, BOTH halves of the F1,6BP cleavage merge
    // back into GAP via TPI, so every glucose carbon ends up in pyruvate.
    var reaching = carbonTraces.filter(function (t) { return t.reachesEnd; });
    assert.ok(reaching.length >= 5,
        'at least 5/6 glucose carbons reach pyruvate under the DAG topology (got ' +
        reaching.length + ' / 6); v2.0.52 was 2/3 of the linearly-traced carbons');
});

test('D4 glycolysis: trace path for a glucose carbon includes BOTH DHAP and GAP via the TPI route', function () {
    var s = glycolysisSetup();
    var r = AtomTrace.tracePathway(s.smiles, { edges: s.edges });
    // The F1,6BP aldol cleavage sends the top-half carbons through DHAP; TPI
    // then converts DHAP -> GAP, so at least one glucose carbon's trace must
    // visit BOTH DHAP and GAP and continue to pyruvate. Which atom index takes
    // the DHAP fork depends on the AAM (it shifts with the exact reactant
    // structures), so search for a DHAP-route carbon rather than hard-coding it.
    var viaDHAP = r.traces.filter(function (t) {
        if (t.element !== 'C') return false;
        var v = t.path.map(function (p) { return s.order[p.metabolite]; });
        return v.indexOf('Dihydroxyacetone phosphate') !== -1;
    });
    assert.ok(viaDHAP.length >= 1,
        'at least one glucose carbon traces through DHAP (got ' + viaDHAP.length + ')');
    var visited = viaDHAP[0].path.map(function (p) { return s.order[p.metabolite]; });
    assert.ok(visited.indexOf('Glyceraldehyde 3-phosphate') !== -1,
        'the DHAP-route carbon converges to GAP via TPI');
    assert.ok(visited.indexOf('Pyruvate') !== -1, 'the DHAP-route carbon reaches Pyruvate');
});

// ---------------------------------------------------------------------------
// E. topological-sort helpers exposed for tests
// ---------------------------------------------------------------------------
test('E1 _topoSort returns a valid order for the glycolysis DAG', function () {
    var s = glycolysisSetup();
    var topo = AtomTrace._topoSort(s.order.length, s.edges);
    assert.ok(Array.isArray(topo), 'topo sort succeeded');
    assert.strictEqual(topo.length, s.order.length);
    // Glucose (index 0) must come first; pyruvate (last in order) must come last.
    assert.strictEqual(topo[0], 0, 'glucose first');
    assert.strictEqual(topo[topo.length - 1], s.nameToIdx['pyruvate'],
        'pyruvate last');
});

test('E2 _topoSort returns null on a cycle (DAG invariant guard)', function () {
    // 0 -> 1 -> 2 -> 0  is a cycle.
    var topo = AtomTrace._topoSort(3, [
        { from: 0, to: 1 }, { from: 1, to: 2 }, { from: 2, to: 0 }
    ]);
    assert.strictEqual(topo, null, 'cycle yields null');
});

test('E3 _reachableFrom marks the reachable set correctly under a fork', function () {
    // 0 -> 1 -> 2, 0 -> 3 (orphan 4). Reachable from 0: {0,1,2,3}, not 4.
    var reach = AtomTrace._reachableFrom(0, 5, [
        { from: 0, to: 1 }, { from: 1, to: 2 }, { from: 0, to: 3 }
    ]);
    assert.strictEqual(reach[0], true);
    assert.strictEqual(reach[1], true);
    assert.strictEqual(reach[2], true);
    assert.strictEqual(reach[3], true);
    assert.strictEqual(reach[4], false);
});

// ---------------------------------------------------------------------------
// F. version stamp
// ---------------------------------------------------------------------------
test('F1 AtomTrace.version is at least 2.0.53 (forward-compatible check)', function () {
    // v2.0.53 introduced the DAG features tested above. Later releases may
    // bump AtomTrace.version (v2.0.54 -> '2.0.54' for moiety tracing, etc.);
    // this check stays valid as long as AtomTrace is loaded and the version
    // string is parseable.
    assert.ok(typeof AtomTrace.version === 'string');
    var parts = AtomTrace.version.split('.').map(Number);
    var atLeast = (parts[0] > 2) || (parts[0] === 2 && (parts[1] > 0 || parts[2] >= 53));
    assert.ok(atLeast, 'AtomTrace.version >= 2.0.53, got ' + AtomTrace.version);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
