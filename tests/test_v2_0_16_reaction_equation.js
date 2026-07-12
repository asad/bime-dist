/**
 * tests/test_v2_0_16_reaction_equation.js — Reaction-equation core (v2.0.16).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Exercises the pure, bundled reaction-equation helpers added to Molecule.js
 * that back the workbench's "+ between fragments" glyphs and the
 * coefficient-collapsed stoichiometric equation in the Editor Output:
 *
 *   Molecule.prototype.elementCounts()            — Hill element counts (+ implicit H)
 *   Molecule.formulaString(counts)                — Hill-notation formula string
 *   Molecule.prototype.computeReactionPlusSigns() — "+" positions between same-side fragments
 *   Molecule.reactionEquation(rMols, pMols, key)  — coefficient equation + balance
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Reaction Equation (v2.0.16)');
var test = runner.test;

console.log('Reaction Equation (v2.0.16)');

function molFromSmiles(smi) {
    var m = new Molecule();
    SmilesParser.parse(smi, m);
    return m;
}

// Canonical-SMILES grouping key, mirroring what the workbench passes so that
// identical species collapse and isomers stay distinct.
function canonKey(m) {
    try { return SmilesWriter.write(m); } catch (e) { return Molecule.formulaString(m.elementCounts()); }
}

// ---------------------------------------------------------------------------
// A. elementCounts + formulaString (Hill notation)
// ---------------------------------------------------------------------------
test('A1 elementCounts: water O has 2 implicit H', function() {
    var c = molFromSmiles('O').elementCounts();
    assert.strictEqual(c.O, 1);
    assert.strictEqual(c.H, 2);
});

test('A2 formulaString: water is H2O (C-less => alphabetical, H included)', function() {
    assert.strictEqual(Molecule.formulaString({ O: 1, H: 2 }), 'H2O');
});

test('A3 formulaString: ethanol is C2H6O (Hill: C, H, then alpha)', function() {
    var c = molFromSmiles('CCO').elementCounts();
    assert.strictEqual(Molecule.formulaString(c), 'C2H6O');
});

test('A4 formulaString: count of 1 omits the digit', function() {
    assert.strictEqual(Molecule.formulaString({ C: 1, H: 4 }), 'CH4');
});

test('A5 formulaString: empty counts -> empty string', function() {
    assert.strictEqual(Molecule.formulaString({}), '');
});

// ---------------------------------------------------------------------------
// B. computeReactionPlusSigns — count, side classification, ordering
// ---------------------------------------------------------------------------

// Build a molecule with N reactant fragments left of x=0 and M product
// fragments right of it, each fragment a single well-separated carbon, plus a
// reaction arrow crossing the midline.
function buildReactionFragments(reactantXs, productXs) {
    var mol = new Molecule();
    var i;
    for (i = 0; i < reactantXs.length; i++) { mol.addAtom('C', reactantXs[i], 0); }
    for (i = 0; i < productXs.length; i++) { mol.addAtom('C', productXs[i], 0); }
    return mol;
}

test('B1 no arrow => no plus signs', function() {
    var mol = buildReactionFragments([-100, -50], [100]);
    assert.deepStrictEqual(mol.computeReactionPlusSigns(null), []);
});

test('B2 single reactant -> single product => no plus signs', function() {
    var mol = buildReactionFragments([-100], [100]);
    mol.reactionArrow = { x1: -10, y1: 0, x2: 10, y2: 0 };
    assert.strictEqual(mol.computeReactionPlusSigns().length, 0);
});

test('B3 two reactants -> one product => exactly one plus (reactant side)', function() {
    var mol = buildReactionFragments([-100, -50], [100]);
    mol.reactionArrow = { x1: -10, y1: 0, x2: 10, y2: 0 };
    var plus = mol.computeReactionPlusSigns();
    assert.strictEqual(plus.length, 1);
    // Placed between the two reactant carbons (x in (-100, -50)).
    assert.ok(plus[0].x > -100 && plus[0].x < -50, 'plus x between reactant fragments: ' + plus[0].x);
});

test('B4 two reactants -> two products => two plus signs (one per side)', function() {
    var mol = buildReactionFragments([-100, -50], [50, 100]);
    mol.reactionArrow = { x1: -10, y1: 0, x2: 10, y2: 0 };
    var plus = mol.computeReactionPlusSigns();
    assert.strictEqual(plus.length, 2);
    var xs = plus.map(function(p) { return p.x; }).sort(function(a, b) { return a - b; });
    assert.ok(xs[0] > -100 && xs[0] < -50, 'reactant-side plus position');
    assert.ok(xs[1] > 50 && xs[1] < 100, 'product-side plus position');
});

test('B5 three reactants -> one product => two plus on reactant side', function() {
    var mol = buildReactionFragments([-150, -100, -50], [100]);
    mol.reactionArrow = { x1: -10, y1: 0, x2: 10, y2: 0 };
    assert.strictEqual(mol.computeReactionPlusSigns().length, 2);
});

test('B6 classification follows arrow midpoint, not sign of x', function() {
    // Shift everything positive; arrow midpoint at x=125 splits 2 vs 1.
    var mol = buildReactionFragments([100, 150], [300]);
    mol.reactionArrow = { x1: 200, y1: 0, x2: 250, y2: 0 }; // midpoint 225
    var plus = mol.computeReactionPlusSigns();
    assert.strictEqual(plus.length, 1);
    assert.ok(plus[0].x > 100 && plus[0].x < 150);
});

// ---------------------------------------------------------------------------
// C. reactionEquation — coefficient grouping, balance, arrow rendering
// ---------------------------------------------------------------------------
test('C1 combustion-style: 2 H2 + O2 -> 2 H2O collapses with coefficients', function() {
    var r = [molFromSmiles('[H][H]'), molFromSmiles('[H][H]'), molFromSmiles('O=O')];
    var p = [molFromSmiles('O'), molFromSmiles('O')];
    var eq = Molecule.reactionEquation(r, p, canonKey);
    assert.strictEqual(eq.reactantText, '2 H2 + O2');
    assert.strictEqual(eq.productText, '2 H2O');
    assert.strictEqual(eq.equation, '2 H2 + O2 → 2 H2O');
    assert.strictEqual(eq.balanced, true);
});

test('C2 single species each side, coefficient 1 omitted', function() {
    var eq = Molecule.reactionEquation([molFromSmiles('CCO')], [molFromSmiles('CC=O')], canonKey);
    assert.strictEqual(eq.reactantText, 'C2H6O');
    assert.strictEqual(eq.productText, 'C2H4O');
    // Lost H2 on oxidation => not balanced.
    assert.strictEqual(eq.balanced, false);
});

test('C3 balanced esterification model is reported balanced', function() {
    // CH3COOH + CH3OH -> CH3COOCH3 + H2O
    var r = [molFromSmiles('CC(=O)O'), molFromSmiles('CO')];
    var p = [molFromSmiles('CC(=O)OC'), molFromSmiles('O')];
    var eq = Molecule.reactionEquation(r, p, canonKey);
    assert.strictEqual(eq.balanced, true);
    assert.ok(eq.equation.indexOf('→') > 0, 'equation uses a Unicode arrow');
});

test('C4 distinct species sharing a formula stay separate under canonical key', function() {
    // Two C2H6O isomers: ethanol (CCO) and dimethyl ether (COC).
    var r = [molFromSmiles('CCO'), molFromSmiles('COC')];
    var eq = Molecule.reactionEquation(r, [], canonKey);
    // Same formula, different canonical SMILES => two coeff-1 terms, not "2 C2H6O".
    assert.strictEqual(eq.reactants.length, 2);
    assert.strictEqual(eq.reactantText, 'C2H6O + C2H6O');
});

test('C5 default keyFn groups by formula (isomers merge)', function() {
    var r = [molFromSmiles('CCO'), molFromSmiles('COC')];
    var eq = Molecule.reactionEquation(r, []); // no keyFn => formula grouping
    assert.strictEqual(eq.reactants.length, 1);
    assert.strictEqual(eq.reactantText, '2 C2H6O');
});

test('C6 empty side renders a placeholder, not a crash', function() {
    var eq = Molecule.reactionEquation([molFromSmiles('O')], []);
    assert.strictEqual(eq.productText, '');
    assert.strictEqual(eq.equation, 'H2O → ?');
});

test('C7 redox half-reaction [Fe+2] -> [Fe+3]: element-balanced but NOT balanced (charge)', function() {
    var eq = Molecule.reactionEquation([molFromSmiles('[Fe+2]')], [molFromSmiles('[Fe+3]')], canonKey);
    assert.strictEqual(eq.elementsBalanced, true, 'same element Fe on both sides');
    assert.strictEqual(eq.chargeBalanced, false, '+2 vs +3 net charge');
    assert.strictEqual(eq.balanced, false, 'overall balance requires charge conservation too');
    assert.strictEqual(eq.netCharge.reactants, 2);
    assert.strictEqual(eq.netCharge.products, 3);
});

test('C8 charge-neutral salt metathesis stays balanced', function() {
    // NaCl + (implicit) -> Na+ . Cl- ; both sides net charge 0, same elements.
    var eq = Molecule.reactionEquation(
        [molFromSmiles('[Na+]'), molFromSmiles('[Cl-]')],
        [molFromSmiles('[Na+]'), molFromSmiles('[Cl-]')], canonKey);
    assert.strictEqual(eq.chargeBalanced, true);
    assert.strictEqual(eq.balanced, true);
    assert.strictEqual(eq.netCharge.reactants, 0);
    assert.strictEqual(eq.netCharge.products, 0);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
