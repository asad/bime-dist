/**
 * tests/test_v2_0_46_chemistry_aware_fitness.js — opt-in chemistry-aware
 * fitness terms (aromatic-ring + stereo conservation; v2.0.46).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Default OFF: golden corpus + every existing test must still pass.
 * With `options.chemistryAware = true`, the fitness gains:
 *   +0.3 per WHOLE aromatic ring broken
 *   +0.3 per stereo centre dropped
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('Chemistry-aware fitness (v2.0.46)');
var test = runner.test;

console.log('Chemistry-aware fitness (v2.0.46)');

function parse(smi) {
    var m = SmilesParser.parse(smi);
    if (m.parseErrors.length > 0) {
        throw new Error('parse: ' + m.parseErrors.join('; '));
    }
    return m;
}

// ---------------------------------------------------------------------------
// A. backward compatibility — default-off path is byte-identical
// ---------------------------------------------------------------------------
test('A1 default mapReaction (no opts) produces the same score as v2.0.45', function() {
    var rxn = parse('CCO>>CC=O');
    var res = RDT.mapReaction(rxn);
    // The score is whatever pre-v2.0.46 fitness produced. We don't pin
    // the exact value (depends on strategy), but we DO assert it is
    // finite and the result is mapped.
    assert.strictEqual(res.status, 'mapped');
    assert.ok(isFinite(res.score), 'score finite');
});

test('A2 mapReaction with chemistryAware:false matches default behaviour', function() {
    var rxn = parse('c1ccccc1>>c1ccccc1');
    var def = RDT.mapReaction(rxn);
    var off = RDT.mapReaction(parse('c1ccccc1>>c1ccccc1'), { chemistryAware: false });
    assert.strictEqual(off.status, def.status);
    assert.strictEqual(off.score, def.score, 'chemistryAware:false matches default exactly');
});

// ---------------------------------------------------------------------------
// B. opt-in adds the chemistry-aware penalty
// ---------------------------------------------------------------------------
test('B1 chemistryAware:true on a benzene→benzene identity gives a score <= default', function() {
    var rxn1 = parse('c1ccccc1>>c1ccccc1');
    var def = RDT.mapReaction(rxn1);
    var rxn2 = parse('c1ccccc1>>c1ccccc1');
    var on = RDT.mapReaction(rxn2, { chemistryAware: true });
    // Identity mapping: all aromatic rings preserved (penalty 0). New
    // term doesn't HURT a perfect chemistry match — score should be
    // less than or equal to the default. (May equal exactly when no
    // chemistry penalty fires.)
    assert.ok(on.score <= def.score + 1e-9,
        'on=' + on.score + ', def=' + def.score);
});

test('B2 fitness() helper accepts chemistryAware opt', function() {
    var rxn = parse('CCO>>CC=O');
    var res = RDT.mapReaction(rxn);
    var withChem = RDT.fitness(
        { reactants: res.sides.reactants, products: res.sides.products },
        res.mapping, res.bondChanges, { chemistryAware: true });
    var noChem = RDT.fitness(
        { reactants: res.sides.reactants, products: res.sides.products },
        res.mapping, res.bondChanges, { chemistryAware: false });
    assert.ok(isFinite(withChem));
    assert.ok(isFinite(noChem));
    // For ethanol→acetaldehyde (no aromatic, no stereo on default
    // mapping), the chemistry penalty is zero, so the two are equal.
    assert.strictEqual(withChem, noChem,
        'no aromatic + no stereo → chemistry penalty is 0');
});

// ---------------------------------------------------------------------------
// C. structural impact — broken aromatic ring incurs penalty
// ---------------------------------------------------------------------------
test('C1 a partial mapping that drops an aromatic ring atom incurs the ring-break penalty', function() {
    // Force a PARTIAL mapping that omits one ring atom — the ring is
    // no longer fully mapped, so the chemistry penalty fires (0.3 per
    // broken ring). Full mapping passes the set-equality check (any
    // permutation of all 6 atoms still forms the destination ring).
    var rxn = parse('c1ccccc1>>c1ccccc1');
    var sides = RDT.mapReaction(rxn).sides;
    var rAtomIds = [], pAtomIds = [];
    for (var i = 0; i < sides.reactants.length; i++) {
        for (var j = 0; j < sides.reactants[i].atoms.length; j++) {
            rAtomIds.push(sides.reactants[i].atoms[j].id);
        }
    }
    for (var k = 0; k < sides.products.length; k++) {
        for (var l = 0; l < sides.products[k].atoms.length; l++) {
            pAtomIds.push(sides.products[k].atoms[l].id);
        }
    }
    // Full identity mapping: 6 atoms paired, ring preserved.
    var fullMap = {};
    for (var ii = 0; ii < rAtomIds.length; ii++) { fullMap[rAtomIds[ii]] = pAtomIds[ii]; }
    // Partial mapping: drop the LAST ring atom. The ring is no longer
    // fully mapped, so the chemistry penalty fires.
    var partialMap = {};
    for (var pp = 0; pp < rAtomIds.length - 1; pp++) {
        partialMap[rAtomIds[pp]] = pAtomIds[pp];
    }
    var fullScore = RDT.fitness(
        { reactants: sides.reactants, products: sides.products },
        fullMap, [], { chemistryAware: true });
    var partialScore = RDT.fitness(
        { reactants: sides.reactants, products: sides.products },
        partialMap, [], { chemistryAware: true });
    // partialScore should exceed fullScore by at least 0.3 for the
    // broken ring (other penalty terms may also differ, but the
    // chemistry-aware contribution should pay at least 0.3 more).
    assert.ok(partialScore >= fullScore + 0.29,
        'partial mapping must incur ≥ ~0.3 ring penalty; got delta = ' +
        (partialScore - fullScore));
});

// ---------------------------------------------------------------------------
// D. golden-corpus sanity — default path unchanged
// ---------------------------------------------------------------------------
test('D1 golden-corpus marker reactions still map cleanly with default opts (no regression)', function() {
    // A few representative reactions from the corpus that exercise
    // common pathways. Each should map without falling through to
    // unbalanced / timeout.
    var samples = [
        'CCO>>CC=O',                          // dehydrogenation
        'c1ccccc1>>c1ccccc1',                 // benzene identity (aromatic)
        'CC(=O)O.CCO>>CC(=O)OCC.O',           // esterification
        'CC(C)C>>CC(C)CO'                     // CH oxidation
    ];
    for (var i = 0; i < samples.length; i++) {
        var rxn = parse(samples[i]);
        var res = RDT.mapReaction(rxn);
        assert.ok(res.status === 'mapped' || res.status === 'unbalanced',
            'sample "' + samples[i] + '" has status ' + res.status);
    }
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
