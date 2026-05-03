/**
 * tests/test_aam.js — Atom-Atom Mapping (AAM) regressions for BIME v1.1.1.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Extends test_rdt.js with strategy/determinism/timeout/edge-case coverage of
 * the four pure-JS RDT pairing strategies (MIN, MAX, MIXTURE, RING).
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('AAM (RDT extended)');
var test = runner.test;

console.log('AAM (RDT extended)');

function parse(rxnSmiles) {
    var m = SmilesParser.parse(rxnSmiles);
    assert.strictEqual(m.parseErrors.length, 0,
        'parse error in ' + rxnSmiles + ': ' + m.parseErrors.join('; '));
    return m;
}

function mapCount(res) { return Object.keys(res.mapping || {}).length; }

// =========================================================================
// (1) Per-strategy mapping size on three reference reactions
// =========================================================================
// Reference 1: ethanol oxidation (3 heavy atoms each side)
test('MIN on CCO>>CC=O maps all 3 heavy atoms', function() {
    assert.strictEqual(mapCount(RDT.runMinPairwise(parse('CCO>>CC=O'))), 3);
});
test('MAX on CCO>>CC=O maps all 3 heavy atoms', function() {
    assert.strictEqual(mapCount(RDT.runMaxPairwise(parse('CCO>>CC=O'))), 3);
});
test('MIXTURE on CCO>>CC=O maps all 3 heavy atoms', function() {
    assert.strictEqual(mapCount(RDT.runMixturePairwise(parse('CCO>>CC=O'))), 3);
});
test('RING on CCO>>CC=O maps all 3 heavy atoms (falls back to MAX-style on ring-free)', function() {
    assert.strictEqual(mapCount(RDT.runRingPairwise(parse('CCO>>CC=O'))), 3);
});

// Reference 2: SN2 on ethyl bromide. RDT's pairwise pass only maps the two
// carbons because each molecule's halogen/hydroxide is on a *different*
// component on the other side; pinning the count protects against silent
// drift in component-pair selection.
test('MIN on CCBr.[OH-]>>CCO.[Br-] maps 2 atoms (the two carbons)', function() {
    assert.strictEqual(mapCount(RDT.runMinPairwise(parse('CCBr.[OH-]>>CCO.[Br-]'))), 2);
});
test('MAX on CCBr.[OH-]>>CCO.[Br-] maps 2 atoms', function() {
    assert.strictEqual(mapCount(RDT.runMaxPairwise(parse('CCBr.[OH-]>>CCO.[Br-]'))), 2);
});
test('MIXTURE on CCBr.[OH-]>>CCO.[Br-] maps 2 atoms', function() {
    assert.strictEqual(mapCount(RDT.runMixturePairwise(parse('CCBr.[OH-]>>CCO.[Br-]'))), 2);
});
test('RING on CCBr.[OH-]>>CCO.[Br-] maps 2 atoms (no rings, falls back)', function() {
    assert.strictEqual(mapCount(RDT.runRingPairwise(parse('CCBr.[OH-]>>CCO.[Br-]'))), 2);
});

// Reference 3: aromatic chlorination — RING/MAX/MIXTURE map all 7 heavy atoms;
// MIN intentionally trims to the smallest matching component.
test('MIN on benzene chlorination maps 1 atom (smallest-match policy: the chlorine)', function() {
    assert.strictEqual(mapCount(RDT.runMinPairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]'))), 1);
});
test('MAX on benzene chlorination maps 7 atoms (6 ring + 1 Cl)', function() {
    assert.strictEqual(mapCount(RDT.runMaxPairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]'))), 7);
});
test('MIXTURE on benzene chlorination maps 7 atoms', function() {
    assert.strictEqual(mapCount(RDT.runMixturePairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]'))), 7);
});
test('RING on benzene chlorination maps 7 atoms (6-ring + Cl filler)', function() {
    assert.strictEqual(mapCount(RDT.runRingPairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]'))), 7);
});

// =========================================================================
// (2) Strategy comparison — all 4 produce *some* mapping for canonical reactions
// =========================================================================
test('All 4 strategies map heavy atoms on Diels-Alder (MIN=2, MAX=4, MIXTURE=4, RING=3)', function() {
    var sm = 'C=CC=C.C=C>>C1CC=CCC1';
    // Per-strategy counts pinned to current BIME policy. RING falls between
    // MIN and MAX because it routes the cyclohexene ring to its own MCS pass.
    assert.strictEqual(mapCount(RDT.runMinPairwise(parse(sm))), 2);
    assert.strictEqual(mapCount(RDT.runMaxPairwise(parse(sm))), 4);
    assert.strictEqual(mapCount(RDT.runMixturePairwise(parse(sm))), 4);
    assert.strictEqual(mapCount(RDT.runRingPairwise(parse(sm))), 3);
});

test('All 4 strategies produce non-empty mappings on amide formation', function() {
    var sm = 'CC(=O)Cl.NCC>>CC(=O)NCC.Cl';
    assert.ok(mapCount(RDT.runMinPairwise(parse(sm))) > 0);
    assert.ok(mapCount(RDT.runMaxPairwise(parse(sm))) > 0);
    assert.ok(mapCount(RDT.runMixturePairwise(parse(sm))) > 0);
    assert.ok(mapCount(RDT.runRingPairwise(parse(sm))) > 0);
});

// =========================================================================
// (3) Determinism — 3 runs of mapReaction give identical structural fingerprint
// =========================================================================
test('mapReaction on esterification: 3 runs give identical signature', function() {
    function sig() {
        var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
        var res = RDT.mapReaction(r);
        var t = res.bondChanges.map(function(e) { return e.type; }).sort();
        return res.strategy + '|' + Object.keys(res.mapping).length + '|' + t.join(',');
    }
    var a = sig(), b = sig(), c = sig();
    assert.strictEqual(a, b);
    assert.strictEqual(b, c);
});

// =========================================================================
// (4) timeoutMs honoured — short timeout returns a well-formed result
// =========================================================================
// Note: with warmed-up V8 small reactions can finish in << 1 ms, so we cannot
// require timedOut === true. We assert the result shape is valid and that the
// timedOut flag is a boolean (i.e. the timer was consulted).
test('Short timeoutMs (1) returns a well-formed result without crashing', function() {
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    var res = RDT.mapReaction(r, { timeoutMs: 1 });
    assert.strictEqual(typeof res.timedOut, 'boolean');
    assert.ok(res.mapping !== null && typeof res.mapping === 'object');
    assert.ok(Array.isArray(res.bondChanges));
});

// =========================================================================
// (5) Pre-mapped atoms preserved
// =========================================================================
test('Pre-mapped atoms keep their map numbers as 1<->1, 2<->2', function() {
    var r = parse('[CH3:1][OH:2]>>[CH2:1]=[O:2]');
    RDT.mapReaction(r);
    var counts = {};
    for (var i = 0; i < r.atoms.length; i++) {
        var n = r.atoms[i].mapNumber || 0;
        if (n > 0) { counts[n] = (counts[n] || 0) + 1; }
    }
    // Each map number must appear exactly twice (one reactant, one product).
    assert.strictEqual(counts[1], 2);
    assert.strictEqual(counts[2], 2);
});

// =========================================================================
// (6) Empty reactant / product / both — graceful (warning, not error)
// =========================================================================
test('Empty reactant side: returns empty mapping with a warning', function() {
    var m = new Molecule();
    m.reactionArrow = { x1: -100, y1: 0, x2: 0, y2: 0 };
    var a = m.addAtom('C', 50, 0);
    var b = m.addAtom('O', 80, 0);
    m.addBond(a.id, b.id, 1);
    var res = RDT.mapReaction(m);
    assert.deepStrictEqual(res.mapping, {});
    assert.ok(res.warnings && res.warnings.length > 0);
});

test('Empty product side: returns empty mapping with a warning', function() {
    var m = new Molecule();
    m.addAtom('C', 0, 0);
    m.addAtom('C', 30, 0);
    m.addBond(m.atoms[0].id, m.atoms[1].id, 1);
    m.reactionArrow = { x1: 100, y1: 0, x2: 200, y2: 0 };
    var res = RDT.mapReaction(m);
    assert.deepStrictEqual(res.mapping, {});
    assert.ok(res.warnings && res.warnings.length > 0);
});

// =========================================================================
// (7) Disconnected reactant / product
// =========================================================================
test('Disconnected reactants (CC.NN>>CCNN) map exactly 2 atoms (one CC pair)', function() {
    var r = parse('CC.NN>>CCNN');
    var res = RDT.mapReaction(r);
    assert.strictEqual(mapCount(res), 2);
});

test('Disconnected products (CCOCC>>CCO.CC) map all 3 reactant heavy atoms', function() {
    var r = parse('CCOCC>>CCO.CC');
    var res = RDT.mapReaction(r);
    assert.strictEqual(mapCount(res), 3);
});

// =========================================================================
// (8) Reaction with explicit hydrogens
// =========================================================================
test('[H][H]>>[H][H] maps both hydrogen atoms (hydrogens treated as atoms)', function() {
    var r = parse('[H][H]>>[H][H]');
    var res = RDT.mapReaction(r);
    assert.strictEqual(mapCount(res), 2);
});

test('Reduction CC=O.[H][H]>>CCO maps the full C-C-O backbone (3 atoms each side)', function() {
    var r = parse('CC=O.[H][H]>>CCO');
    var res = RDT.mapReaction(r);
    // The 3 heavy atoms must be mapped on both sides (= 6 mapped atoms total).
    var heavyMapped = 0;
    for (var i = 0; i < r.atoms.length; i++) {
        if (r.atoms[i].symbol !== 'H' && r.atoms[i].mapNumber > 0) { heavyMapped++; }
    }
    assert.strictEqual(heavyMapped, 6);
});

// =========================================================================
// (9) options.strategies restricts the search set
// =========================================================================
test('options.strategies = [MIN] pins MIN as the winner', function() {
    var res = RDT.mapReaction(parse('CCO>>CC=O'), { strategies: ['MIN'] });
    assert.strictEqual(res.strategy, 'MIN');
});

test('options.strategies = [RING] pins RING as the winner', function() {
    var res = RDT.mapReaction(parse('CCO>>CC=O'), { strategies: ['RING'] });
    assert.strictEqual(res.strategy, 'RING');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
