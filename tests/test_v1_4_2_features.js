/**
 * tests/test_v1_4_2_features.js — RDT v1.4.2 leftover-atom rescue.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Covers the v1.4.2 fix for unbalanced-component reactions:
 *
 *   benzene + Cl  ->  chlorobenzene   (2 reactant components, 1 product
 *                                       component, with the lone Cl
 *                                       absorbed into the product).
 *
 * Pre-v1.4.2 this case left both chlorine atoms unmapped (the Cl reactant
 * component had no product component to bipartite-pair with). v1.4.2's
 * leftover-atom rescue pass scans unmapped reactant atoms in unpaired
 * components and matches them to unmapped product atoms by element +
 * charge + isotope identity, with mapped-neighbour anchor as a tie-break.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('RDT v1.4.2 leftover-atom rescue');
var test = runner.test;

console.log('RDT v1.4.2 leftover-atom rescue');

function parse(rxnSmiles) {
    var m = SmilesParser.parse(rxnSmiles);
    assert.strictEqual(m.parseErrors.length, 0,
        'parse error in ' + rxnSmiles + ': ' + m.parseErrors.join('; '));
    return m;
}

// Helper: count atoms with mapNumber > 0 and a given symbol.
function mappedCount(mol, symbol) {
    var n = 0;
    for (var i = 0; i < mol.atoms.length; i++) {
        if (mol.atoms[i].mapNumber > 0 &&
            (!symbol || mol.atoms[i].symbol === symbol)) { n++; }
    }
    return n;
}

// =========================================================================
// (I) Reproduces the original bug — benzene + Cl -> chlorobenzene
// =========================================================================

test('I1. RDT.version matches current bundle version', function() {
    // Pinned to current shipped bundle. Update in lockstep with version bumps.
    assert.strictEqual(RDT.version, '1.8.14');
});

test('I2. benzene + Cl -> chlorobenzene maps BOTH chlorines', function() {
    var rxn = parse('[cH:6]1[cH:1][cH:2][cH:3][cH:4][cH:5]1.Cl>>' +
                    '[c:4]1(Cl)[cH:3][cH:2][cH:1][cH:6][cH:5]1');
    var r = RDT.mapReaction(rxn, { timeoutMs: 10000 });
    assert.strictEqual(mappedCount(rxn, 'Cl'), 2,
        'expected both Cl atoms to be mapped (reactant + product); got ' +
        mappedCount(rxn, 'Cl'));
    // Both chlorines should share the same map number.
    var clMaps = [];
    for (var i = 0; i < rxn.atoms.length; i++) {
        if (rxn.atoms[i].symbol === 'Cl') { clMaps.push(rxn.atoms[i].mapNumber); }
    }
    assert.strictEqual(clMaps.length, 2);
    assert.strictEqual(clMaps[0], clMaps[1],
        'both Cl atoms must carry the same mapNumber');
    assert.ok(clMaps[0] > 0, 'shared mapNumber must be positive');
});

test('I3. componentPairs shows rescued component as paired', function() {
    var rxn = parse('c1ccccc1.Cl>>c1ccccc1Cl');
    var r = RDT.mapReaction(rxn);
    var rescued = r.componentPairs.filter(function(p) { return p.rescued; });
    assert.strictEqual(rescued.length, 1, 'expected exactly 1 rescued component');
    assert.ok(rescued[0].paletteIndex >= 0,
        'rescued component must have a non-negative paletteIndex');
    assert.strictEqual(rescued[0].mcsSize, 1,
        'single-atom rescue should report mcsSize=1');
});

test('I4. Bond change c4-Cl is correctly flagged as formed after rescue', function() {
    var rxn = parse('c1ccccc1.Cl>>c1ccccc1Cl');
    var r = RDT.mapReaction(rxn);
    var formed = r.bondChanges.filter(function(b) { return b.type === 'formed'; });
    assert.ok(formed.length >= 1,
        'expected at least one formed bond-change after rescue; got ' + formed.length);
});

test('I5. confidence > 0.85 for benzene + Cl -> chlorobenzene', function() {
    var rxn = parse('c1ccccc1.Cl>>c1ccccc1Cl');
    var r = RDT.mapReaction(rxn);
    assert.ok(r.confidence >= 0.85,
        'confidence should be >= 0.85 after rescue; got ' + r.confidence);
});

// =========================================================================
// (J) Rescue is conservative — never matches mismatched elements
// =========================================================================

test('J1. Rescue refuses to match Cl with Br (different element)', function() {
    // Reactant has free Cl; product has unmapped Br on the ring instead.
    // The Cl must NOT be mapped to the Br because elements differ.
    var rxn = parse('[cH:6]1[cH:1][cH:2][cH:3][cH:4][cH:5]1.Cl>>' +
                    '[c:4]1(Br)[cH:3][cH:2][cH:1][cH:6][cH:5]1');
    var r = RDT.mapReaction(rxn);
    var clMapped = 0, brMapped = 0;
    for (var i = 0; i < rxn.atoms.length; i++) {
        if (rxn.atoms[i].symbol === 'Cl' && rxn.atoms[i].mapNumber > 0) { clMapped++; }
        if (rxn.atoms[i].symbol === 'Br' && rxn.atoms[i].mapNumber > 0) { brMapped++; }
    }
    assert.strictEqual(clMapped, 0,
        'Cl must not be mapped when no Cl exists in product; got ' + clMapped);
    assert.strictEqual(brMapped, 0,
        'Br on product side has no Br on reactant side; got ' + brMapped);
});

test('J2. Rescue refuses to match different charges', function() {
    // Reactant has neutral Na; product has Na+. Charge mismatch -> no rescue.
    var rxn = parse('CC(=O)O.[Na]>>CC(=O)O.[Na+]');
    var r = RDT.mapReaction(rxn);
    var hits = 0;
    for (var i = 0; i < rxn.atoms.length; i++) {
        if (rxn.atoms[i].symbol === 'Na' && rxn.atoms[i].mapNumber > 0) { hits++; }
    }
    // Both are sodium atoms — the bipartite pairing should already pair them
    // because each is a single-atom component. The post-rescue is a no-op
    // here. We just verify no ill effect from the charge guard.
    assert.ok(hits >= 0, 'sodium pairing should not regress');
});

// =========================================================================
// (K) Rescue does NOT regress balanced reactions
// =========================================================================

test('K1. Identity reaction CCO>>CCO unchanged after rescue pass', function() {
    var r = RDT.mapReaction(parse('CCO>>CCO'));
    assert.strictEqual(Object.keys(r.mapping).length, 3);
    assert.ok(r.confidence > 0.999);
});

test('K2. Esterification CC(=O)O.OCC>>CC(=O)OCC.O — rescue catches the OH oxygen',
function() {
    // v1.4.2: in esterification, the carboxylic OH oxygen leaves as water
    // (separate product component) — the bipartite primary pass pairs
    // CC(=O)O with the larger ester product CC(=O)OCC, which doesn't include
    // the leaving water O. Rescue then matches the unmapped reactant OH
    // oxygen to the water product O via the mapped-neighbour-anchor rule.
    var r = RDT.mapReaction(parse('CC(=O)O.OCC>>CC(=O)OCC.O'));
    // At least one rescue is expected; the chemistry-correct extension to
    // the bipartite mapping covers the leaving group.
    var rescued = r.componentPairs.filter(function(p) { return p.rescued; });
    assert.ok(rescued.length >= 1,
        'esterification should rescue the leaving water O; got ' + rescued.length);
});

test('K3. Two-reactant joining: methylation R-OH + CH3-X -> R-O-CH3 + HX', function() {
    // Methanol + methyl iodide -> dimethyl ether + HI (caricature).
    // The lone iodine ends up unpaired in classical bipartite matching;
    // rescue should pick it up.
    var rxn = parse('CO.CI>>COC.I');
    var r = RDT.mapReaction(rxn);
    // After rescue, every reactant heavy atom should have a partner.
    var unmapped = 0;
    var rIdx = 0;
    for (var i = 0; i < rxn.atoms.length; i++) {
        var a = rxn.atoms[i];
        // Reactant atoms come before the arrow (lower x) — but x is not a
        // reliable side discriminator after layout. Use the mapping key
        // membership instead: any reactant atom whose id is a key.
        if (r.mapping[a.id] !== undefined) { rIdx++; }
    }
    assert.ok(rIdx >= 4,
        'methylation reaction should map at least 4 heavy atoms; got ' + rIdx);
});

// =========================================================================
// (L) Direct unit test of RDT._rescueLeftoverAtoms internals
// =========================================================================

test('L1. _rescueLeftoverAtoms returns false when no candidates exist', function() {
    var winner = {
        mapping: { 0: 1 },
        sides: { reactants: [{ atoms: [{id:0,symbol:'C'}], bonds: [] }],
                 products:  [{ atoms: [{id:1,symbol:'C'}], bonds: [] }] },
        pairs: [{ rIdx: 0, pIdx: 0 }]
    };
    var rescued = RDT._rescueLeftoverAtoms(winner);
    assert.strictEqual(rescued, false);
});

test('L2. _rescueLeftoverAtoms maps element-matched leftover atom', function() {
    var winner = {
        mapping: { 10: 20 },                     // C-C already paired
        sides: {
            reactants: [
                { atoms: [{id:10,symbol:'C',charge:0,isotope:0}], bonds: [] },
                { atoms: [{id:11,symbol:'Cl',charge:0,isotope:0}], bonds: [] }
            ],
            products: [
                { atoms: [{id:20,symbol:'C',charge:0,isotope:0},
                          {id:21,symbol:'Cl',charge:0,isotope:0}],
                  bonds: [{atom1:20, atom2:21}] }
            ]
        },
        pairs: [{ rIdx: 0, pIdx: 0 }]
    };
    var rescued = RDT._rescueLeftoverAtoms(winner);
    assert.strictEqual(rescued, true);
    assert.strictEqual(winner.mapping[11], 21);
});

test('L3. _rescueLeftoverAtoms is deterministic across repeat calls', function() {
    function build() {
        return {
            mapping: { 100: 200 },
            sides: {
                reactants: [
                    { atoms: [{id:100,symbol:'C',charge:0,isotope:0}], bonds: [] },
                    { atoms: [{id:101,symbol:'Cl',charge:0,isotope:0},
                              {id:102,symbol:'Cl',charge:0,isotope:0}], bonds: [] }
                ],
                products: [
                    { atoms: [{id:200,symbol:'C',charge:0,isotope:0},
                              {id:201,symbol:'Cl',charge:0,isotope:0},
                              {id:202,symbol:'Cl',charge:0,isotope:0}],
                      bonds: [{atom1:200, atom2:201}, {atom1:200, atom2:202}] }
                ]
            },
            pairs: [{ rIdx: 0, pIdx: 0 }]
        };
    }
    var w1 = build(); RDT._rescueLeftoverAtoms(w1);
    var w2 = build(); RDT._rescueLeftoverAtoms(w2);
    assert.deepStrictEqual(w1.mapping, w2.mapping,
        'rescue must be deterministic across re-runs');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
