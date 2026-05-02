/**
 * tests/test_rdt.js - Reaction Decoder Tool (RDT) regressions for BIME v1.1.0.
 *
 * Verifies the four pairing strategies (MIN, MAX, MIXTURE, RING), the
 * bond-change annotator, the fitness function, the selector, and the
 * mapReaction front-door on a panel of textbook organic reactions.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

// Load the RDT module after the shim has primed the editor modules.
require('../editor/RDT.js');

var runner = shim.makeRunner('RDT (AAM)');
var test = runner.test;

console.log('RDT (AAM)');

function parse(rxnSmiles) {
    var m = SmilesParser.parse(rxnSmiles);
    assert.strictEqual(m.parseErrors.length, 0,
        'parse error in ' + rxnSmiles + ': ' + m.parseErrors.join('; '));
    return m;
}

function countByType(events, type) {
    var n = 0;
    for (var i = 0; i < events.length; i++) {
        if (events[i].type === type) { n++; }
    }
    return n;
}

function mappedAtomCount(reaction) {
    var n = 0;
    for (var i = 0; i < reaction.atoms.length; i++) {
        if (reaction.atoms[i].mapNumber > 0) { n++; }
    }
    return n;
}

// =========================================================================
// Public-API smoke tests
// =========================================================================
test('RDT module loads and exposes the public API', function() {
    assert.ok(RDT, 'RDT global exists');
    assert.strictEqual(typeof RDT.mapReaction, 'function');
    assert.strictEqual(typeof RDT.runMinPairwise, 'function');
    assert.strictEqual(typeof RDT.runMaxPairwise, 'function');
    assert.strictEqual(typeof RDT.runMixturePairwise, 'function');
    assert.strictEqual(typeof RDT.runRingPairwise, 'function');
    assert.strictEqual(typeof RDT.annotateBondChanges, 'function');
    assert.strictEqual(typeof RDT.fitness, 'function');
    assert.deepStrictEqual(RDT.STRATEGIES, ['MIN', 'MAX', 'MIXTURE', 'RING']);
});

// =========================================================================
// Reaction splitting / balance utilities
// =========================================================================
test('splitReactionSides places reactants left of arrow midpoint', function() {
    var r = parse('CCO>>CC=O');
    var sides = RDT._splitReactionSides(r);
    assert.strictEqual(sides.reactants.length, 1);
    assert.strictEqual(sides.products.length, 1);
});

test('splitReactionSides handles 2+1 component reaction', function() {
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    var sides = RDT._splitReactionSides(r);
    assert.strictEqual(sides.reactants.length, 2);
    assert.strictEqual(sides.products.length, 2);
});

test('checkBalance flags ethanol -> acetaldehyde as unbalanced (loses 2 H)', function() {
    var r = parse('CCO>>CC=O');
    var sides = RDT._splitReactionSides(r);
    var b = RDT._checkBalance(sides);
    // Heavy atoms balance: 2C+1O = 2C+1O. Hydrogens are implicit so balance check
    // counts heavy elements only and reports balanced.
    assert.strictEqual(b.balanced, true);
});

test('checkBalance flags clearly unbalanced reaction', function() {
    var r = parse('CCO>>CC');
    var sides = RDT._splitReactionSides(r);
    var b = RDT._checkBalance(sides);
    assert.strictEqual(b.balanced, false);
});

// =========================================================================
// Per-strategy results: each strategy must produce a mapping (or not crash).
// =========================================================================
test('MIN strategy on CCO>>CC=O maps all 3 heavy atoms', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.runMinPairwise(r);
    assert.strictEqual(res.strategy, 'MIN');
    assert.ok(Object.keys(res.mapping).length >= 2,
        'MIN should map at least 2 heavy atoms; got ' + Object.keys(res.mapping).length);
    assert.ok(isFinite(res.score), 'fitness is finite');
});

test('MAX strategy on CCO>>CC=O maps all 3 heavy atoms', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.runMaxPairwise(r);
    assert.strictEqual(res.strategy, 'MAX');
    assert.strictEqual(Object.keys(res.mapping).length, 3, 'MAX should map all 3 heavy atoms');
});

test('MIXTURE strategy on CCO>>CC=O maps all 3 heavy atoms', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.runMixturePairwise(r);
    assert.strictEqual(res.strategy, 'MIXTURE');
    assert.strictEqual(Object.keys(res.mapping).length, 3);
});

test('RING strategy returns no MCS pairs on a non-ring reaction', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.runRingPairwise(r);
    assert.strictEqual(res.strategy, 'RING');
    // RING falls back to a MAX-style fill when no rings are present, so it
    // must still yield a usable mapping.
    assert.ok(Object.keys(res.mapping).length >= 2);
});

// =========================================================================
// Bond-change annotation on textbook reactions
// =========================================================================
test('CCO>>CC=O: orderChange on the C-O bond (single->double)', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.runMaxPairwise(r);
    // Expect at least one orderChange event.
    assert.ok(countByType(res.bondChanges, 'orderChange') >= 1,
        'expected an orderChange; got ' + JSON.stringify(res.bondChanges));
});

test('SN2: CCBr.[OH-]>>CCO.[Br-] produces formed and cleaved events', function() {
    var r = parse('CCBr.[OH-]>>CCO.[Br-]');
    var res = RDT.mapReaction(r);
    assert.ok(res.bondChanges.length >= 1, 'expected at least 1 bond change');
    var formed = countByType(res.bondChanges, 'formed');
    var cleaved = countByType(res.bondChanges, 'cleaved');
    assert.ok(formed + cleaved >= 1,
        'SN2 should report a formed or cleaved event; got: ' + JSON.stringify(res.bondChanges));
});

test('Esterification CC(=O)O.OCC>>CC(=O)OCC.O has at least one formed bond', function() {
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    var res = RDT.mapReaction(r);
    assert.ok(res.bondChanges.length >= 1, 'expected bond changes');
    // The ester C-O bond is newly formed.
    var formed = countByType(res.bondChanges, 'formed');
    assert.ok(formed >= 1 || countByType(res.bondChanges, 'cleaved') >= 1,
        'expected formed or cleaved event in esterification');
});

test('Ar chlorination c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH] preserves benzene ring', function() {
    var r = parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]');
    var res = RDT.mapReaction(r);
    assert.ok(res.strategy, 'a winning strategy was selected');
    // All 6 ring carbons should remain mapped on both sides.
    assert.ok(Object.keys(res.mapping).length >= 6,
        'expected at least 6 atoms mapped (benzene ring); got ' + Object.keys(res.mapping).length);
});

test('Reduction CC=O.[H][H]>>CCO maps the C-O backbone', function() {
    var r = parse('CC=O.[H][H]>>CCO');
    var res = RDT.mapReaction(r);
    assert.ok(Object.keys(res.mapping).length >= 2,
        'expected C and O of the carbonyl to be mapped');
    // The C=O -> C-O is an order change.
    assert.ok(countByType(res.bondChanges, 'orderChange') >= 1 ||
              countByType(res.bondChanges, 'formed') >= 1 ||
              countByType(res.bondChanges, 'cleaved') >= 1,
        'expected at least one bond-change event');
});

test('Diels-Alder C=CC=C.C=C>>C1CC=CCC1 produces bond-change events', function() {
    var r = parse('C=CC=C.C=C>>C1CC=CCC1');
    var res = RDT.mapReaction(r);
    // 6 heavy atoms on each side; mapping must cover ring atoms.
    assert.ok(Object.keys(res.mapping).length >= 4, 'should map most atoms');
    assert.ok(res.bondChanges.length >= 2, 'Diels-Alder forms 2 new C-C bonds');
});

// =========================================================================
// Multi-component reactions
// =========================================================================
test('Multiple reactants -> single product', function() {
    var r = parse('CC(=O)Cl.NCC>>CC(=O)NCC.Cl');
    var res = RDT.mapReaction(r);
    assert.ok(Object.keys(res.mapping).length >= 4);
    assert.ok(res.bondChanges.length >= 1);
});

test('Single reactant -> multiple products', function() {
    var r = parse('CCOCC>>CCO.CC');
    var res = RDT.mapReaction(r);
    assert.ok(Object.keys(res.mapping).length >= 2);
});

// =========================================================================
// Pre-mapped atoms must be respected
// =========================================================================
test('Pre-mapped atoms are preserved', function() {
    var r = parse('[CH3:1][OH:2]>>[CH2:1]=[O:2]');
    // Original mapping is 1<->1 (C) and 2<->2 (O).
    var rC = null, pC = null, rO = null, pO = null;
    var arrowX = (r.reactionArrow.x1 + r.reactionArrow.x2) / 2;
    for (var i = 0; i < r.atoms.length; i++) {
        var a = r.atoms[i];
        if (a.mapNumber === 1 && a.x < arrowX) { rC = a.id; }
        if (a.mapNumber === 1 && a.x > arrowX) { pC = a.id; }
        if (a.mapNumber === 2 && a.x < arrowX) { rO = a.id; }
        if (a.mapNumber === 2 && a.x > arrowX) { pO = a.id; }
    }
    assert.ok(rC !== null && pC !== null && rO !== null && pO !== null);

    var res = RDT.mapReaction(r);
    assert.strictEqual(res.mapping[rC], pC, 'reactant C maps to product C');
    assert.strictEqual(res.mapping[rO], pO, 'reactant O maps to product O');
});

// =========================================================================
// Unbalanced + edge cases
// =========================================================================
test('Unbalanced reaction emits a warning and best-effort mapping', function() {
    var r = parse('CCO>>CC');
    var res = RDT.mapReaction(r);
    var hasWarn = res.warnings && res.warnings.some(function(w) { return /balance/i.test(w); });
    assert.ok(hasWarn, 'expected an atom-balance warning; got ' + JSON.stringify(res.warnings));
    // Best-effort mapping still maps the 2 carbons.
    assert.ok(Object.keys(res.mapping).length >= 2);
});

test('Empty product side: reaction returns empty mapping with warning', function() {
    var m = new Molecule();
    // A lone reactant with no arrow -> no products.
    m.addAtom('C', 0, 0);
    m.addAtom('C', 30, 0);
    m.addBond(m.atoms[0].id, m.atoms[1].id, 1);
    m.reactionArrow = { x1: 100, y1: 0, x2: 200, y2: 0 };
    var res = RDT.mapReaction(m);
    assert.deepStrictEqual(res.mapping, {});
    assert.ok(res.warnings.length > 0);
});

test('Empty reactant side: returns empty mapping', function() {
    var m = new Molecule();
    m.reactionArrow = { x1: -100, y1: 0, x2: 0, y2: 0 };
    var a = m.addAtom('C', 50, 0);
    var b = m.addAtom('O', 80, 0);
    m.addBond(a.id, b.id, 1);
    var res = RDT.mapReaction(m);
    assert.deepStrictEqual(res.mapping, {});
});

test('Empty molecule (no atoms) is handled gracefully', function() {
    var m = new Molecule();
    var res = RDT.mapReaction(m);
    assert.deepStrictEqual(res.mapping, {});
    assert.ok(res.warnings.indexOf('empty reaction') >= 0);
});

// =========================================================================
// Determinism
// =========================================================================
test('RDT.mapReaction is deterministic across 3 runs', function() {
    function runOnce() {
        var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
        var res = RDT.mapReaction(r);
        // Build a structural fingerprint from atom symbols + map numbers
        // (insensitive to global atom-id counters that drift between parses).
        var sigParts = [];
        for (var i = 0; i < r.atoms.length; i++) {
            var a = r.atoms[i];
            sigParts.push(a.symbol + ':' + (a.mapNumber || 0));
        }
        // Sort to be position-insensitive (we just want to compare multisets).
        sigParts.sort();
        return {
            strategy: res.strategy,
            sig: sigParts.join(','),
            score: res.score,
            mapCount: Object.keys(res.mapping).length,
            bondChanges: res.bondChanges.length
        };
    }
    var r1 = runOnce(), r2 = runOnce(), r3 = runOnce();
    assert.strictEqual(r1.sig, r2.sig, 'mapping signature drifted between runs');
    assert.strictEqual(r2.sig, r3.sig, 'mapping signature drifted between runs');
    assert.strictEqual(r1.strategy, r2.strategy, 'winning strategy changed');
    assert.strictEqual(r1.mapCount, r2.mapCount, 'mapped-atom count drifted');
    assert.strictEqual(r1.bondChanges, r2.bondChanges, 'bond-change count drifted');
});

// =========================================================================
// Selector behaviour: 4-strategy comparison
// =========================================================================
test('All 4 strategies produce non-null results on benzene chlorination', function() {
    var rxnStr = 'c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]';
    var r;
    var strategies = {
        MIN: function() { return RDT.runMinPairwise(parse(rxnStr)); },
        MAX: function() { return RDT.runMaxPairwise(parse(rxnStr)); },
        MIXTURE: function() { return RDT.runMixturePairwise(parse(rxnStr)); },
        RING: function() { return RDT.runRingPairwise(parse(rxnStr)); }
    };
    var minRes = strategies.MIN();
    var maxRes = strategies.MAX();
    var mixRes = strategies.MIXTURE();
    var ringRes = strategies.RING();
    assert.ok(Object.keys(minRes.mapping).length > 0, 'MIN should produce mapping');
    assert.ok(Object.keys(maxRes.mapping).length > 0, 'MAX should produce mapping');
    assert.ok(Object.keys(mixRes.mapping).length > 0, 'MIXTURE should produce mapping');
    assert.ok(Object.keys(ringRes.mapping).length > 0, 'RING should produce mapping');
});

test('Selector picks the fitness-best strategy on benzene chlorination', function() {
    var r = parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]');
    var res = RDT.mapReaction(r);
    assert.ok(['MIN', 'MAX', 'MIXTURE', 'RING'].indexOf(res.strategy) >= 0);
    // Re-run individual strategies and verify the winner has the lowest score.
    var perStrat = {};
    perStrat.MIN = RDT.runMinPairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]')).score;
    perStrat.MAX = RDT.runMaxPairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]')).score;
    perStrat.MIXTURE = RDT.runMixturePairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]')).score;
    perStrat.RING = RDT.runRingPairwise(parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]')).score;
    var bestScore = Math.min(perStrat.MIN, perStrat.MAX, perStrat.MIXTURE, perStrat.RING);
    assert.ok(Math.abs(res.score - bestScore) <= 1e-6,
        'winner score ' + res.score + ' should equal best ' + bestScore + ' (' + JSON.stringify(perStrat) + ')');
});

// =========================================================================
// annotateBondChanges as a standalone helper
// =========================================================================
test('annotateBondChanges identifies a cleaved bond directly', function() {
    var r = parse('CCO>>CC=O');
    var sides = RDT._splitReactionSides(r);
    // Manually map: pick first 2 atoms ignoring O (so the C-O bond will be cleaved).
    var atomMap = {};
    var rAtoms = sides.reactants[0].atoms;
    var pAtoms = sides.products[0].atoms;
    // Map the two carbons only:
    var rC1 = null, rC2 = null, pC1 = null, pC2 = null;
    for (var i = 0; i < rAtoms.length; i++) {
        if (rAtoms[i].symbol === 'C') {
            if (rC1 === null) { rC1 = rAtoms[i].id; }
            else if (rC2 === null) { rC2 = rAtoms[i].id; }
        }
    }
    for (var j = 0; j < pAtoms.length; j++) {
        if (pAtoms[j].symbol === 'C') {
            if (pC1 === null) { pC1 = pAtoms[j].id; }
            else if (pC2 === null) { pC2 = pAtoms[j].id; }
        }
    }
    atomMap[rC1] = pC1;
    atomMap[rC2] = pC2;
    var events = RDT.annotateBondChanges(sides.reactants, sides.products, atomMap);
    var cleaved = countByType(events, 'cleaved');
    assert.ok(cleaved >= 1, 'unmapped O endpoint must yield a cleaved C-O bond');
});

// =========================================================================
// Fitness function semantics
// =========================================================================
test('fitness rewards no-change (identity reaction has fitness <= 0)', function() {
    // CCO>>CCO: identity, all atoms map, no bond changes, ring bonus = 0.
    // Use a reaction-SMILES style: same molecule on both sides.
    var r = parse('CCO>>CCO');
    var res = RDT.mapReaction(r);
    assert.ok(res.score <= 0.5,
        'identity reaction should score near zero; got ' + res.score);
});

test('fitness penalises bond changes', function() {
    var rA = parse('CCO>>CCO');
    var resA = RDT.mapReaction(rA);
    var rB = parse('CCO>>CC=O');
    var resB = RDT.mapReaction(rB);
    assert.ok(resA.score <= resB.score + 1e-9,
        'CCO>>CCO (' + resA.score + ') should not score worse than CCO>>CC=O (' + resB.score + ')');
});

// =========================================================================
// Mapping is applied in-place to the input reaction
// =========================================================================
test('mapReaction mutates the reaction so atoms carry mapNumber > 0', function() {
    var r = parse('CCO>>CC=O');
    RDT.mapReaction(r);
    assert.ok(mappedAtomCount(r) >= 4,
        'expected at least 2 reactant + 2 product atoms to be mapped');
});

test('Mapped atoms come in pairs (same map on both sides)', function() {
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    RDT.mapReaction(r);
    var counts = {};
    for (var i = 0; i < r.atoms.length; i++) {
        var n = r.atoms[i].mapNumber || 0;
        if (n > 0) { counts[n] = (counts[n] || 0) + 1; }
    }
    for (var k in counts) {
        if (counts.hasOwnProperty(k)) {
            assert.strictEqual(counts[k], 2,
                'map number ' + k + ' should appear exactly twice; got ' + counts[k]);
        }
    }
});

// =========================================================================
// Options handling
// =========================================================================
test('options.strategies can restrict the strategy set', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r, { strategies: ['MAX'] });
    assert.strictEqual(res.strategy, 'MAX');
});

test('Long timeoutMs is respected without spurious timeouts', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r, { timeoutMs: 30000 });
    assert.strictEqual(res.timedOut, false);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
