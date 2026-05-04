/**
 * tests/test_aam_round_trip.js — Reaction SMILES with atom-atom mapping
 * (AAM) round-trip regression for BIME v1.0.2.
 *
 * Each case is a textbook organic transformation expressed in standard
 * Daylight reaction-SMILES syntax (`reactants>>products`, optionally
 * `reactants>agents>products`). Atom-map labels (`[X:n]`) annotate which
 * heavy atoms correspond across the arrow.
 *
 * The round-trip contract under test:
 *   1. `SmilesParser.parse()` accepts the reaction string with no errors,
 *      records the reactionArrow, and assigns the integer mapNumber to
 *      each atom that carries a map.
 *   2. `SmilesWriter.write(mol)` emits a string that itself contains `>>`
 *      and re-parses to the same atom count and the same multiset of map
 *      numbers.
 *
 * Atom *order* within each reactant/product fragment is allowed to drift
 * (the writer chooses its own canonical traversal); the *set* of maps
 * must survive intact.
 *
 * No file I/O, no network, runs in <50 ms.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('AAM Round-Trip');
var test = runner.test;

console.log('AAM Round-Trip');

function P(smi) { return SmilesParser.parse(smi); }
function W(mol) { return SmilesWriter.write(mol); }

function mapMultiset(mol) {
    return mol.atoms
        .map(function(a) { return a.mapNumber || 0; })
        .filter(function(n) { return n > 0; })
        .sort(function(a, b) { return a - b; });
}

// Reaction SMILES with full atom-atom mapping.
// Each test asserts:
//   parse → no errors
//   parse → reactionArrow recorded
//   parse → expected atom count
//   parse → multiset of map numbers matches expected
//   write → output contains >>
//   re-parse → atom count and map multiset preserved
var CASES = [
    {
        name: 'esterification (acetic acid + methanol)',
        smi: '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][OH:6]>>[CH3:1][C:2](=[O:3])[O:6][CH3:5].[OH2:4]',
        atoms: 12,
        // Each atom number 1..6 appears once on each side: multiset has each map twice.
        expectedMaps: [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
    },
    {
        name: 'acetylene hydration (Markovnikov)',
        smi: '[CH:1]#[CH:2].[OH2:3]>>[CH3:1][CH:2]=[O:3]',
        atoms: 6,
        expectedMaps: [1, 1, 2, 2, 3, 3]
    },
    {
        name: 'amide formation (methylamine + acetic acid)',
        smi: '[NH2:1][CH3:2].[CH3:3][C:4](=[O:5])[OH:6]>>[CH3:3][C:4](=[O:5])[NH:1][CH3:2].[OH2:6]',
        atoms: 12,
        expectedMaps: [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
    },
    {
        name: 'Claisen condensation (two methyl acetate)',
        smi: '[CH3:1][C:2](=[O:3])[O:4][CH3:5].[CH3:6][C:7](=[O:8])[O:9][CH3:10]>>[CH3:1][C:2](=[O:3])[CH2:6][C:7](=[O:8])[O:9][CH3:10].[CH3:5][OH:4]',
        atoms: 20,
        expectedMaps: [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10]
    },
    {
        name: 'SN2 (methyl bromide + hydroxide)',
        smi: '[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]',
        atoms: 6,
        expectedMaps: [1, 1, 2, 2, 3, 3]
    },
    {
        name: 'Diels–Alder (butadiene + ethylene)',
        smi: '[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]>>[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:5][CH2:6]1',
        atoms: 12,
        expectedMaps: [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
    },
    {
        name: 'partial AAM (only some atoms mapped)',
        smi: '[CH3:1][OH:2].C(=O)O>>[CH3:1][O:2]C(=O)C',
        atoms: 10,
        expectedMaps: [1, 1, 2, 2]
    },
    {
        name: 'unmapped reaction (oxidation: ethanol -> acetaldehyde)',
        smi: 'CCO>>CC=O',
        atoms: 6,
        expectedMaps: []
    },
    {
        name: 'aromatic AAM (benzene chlorination)',
        smi: '[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1.[Cl:7][Cl:8]>>[c:1]1([Cl:7])[cH:2][cH:3][cH:4][cH:5][cH:6]1.[ClH:8]',
        atoms: 16,
        expectedMaps: [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]
    },
    {
        name: 'salt fragment with AAM (sodium hydroxide + methyl chloride)',
        smi: '[Na+:1].[OH-:2].[CH3:3][Cl:4]>>[Na+:1].[Cl-:4].[CH3:3][OH:2]',
        atoms: 8,
        expectedMaps: [1, 1, 2, 2, 3, 3, 4, 4]
    }
];

CASES.forEach(function(c) {
    test(c.name + ': parses with arrow and full map multiset', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0,
                           'parse errors: ' + m.parseErrors.join('; '));
        assert.ok(m.reactionArrow, 'expected reactionArrow to be set');
        assert.strictEqual(m.atoms.length, c.atoms,
                           'atom count: expected ' + c.atoms + ', got ' + m.atoms.length);
        assert.deepStrictEqual(mapMultiset(m), c.expectedMaps,
                               'atom-map multiset mismatch');
    });

    test(c.name + ': round-trips through writer (arrow + maps preserved)', function() {
        var m1 = P(c.smi);
        var written = W(m1);
        assert.ok(written.indexOf('>>') >= 0,
                  'expected >> in writer output: ' + written);
        var m2 = P(written);
        assert.strictEqual(m2.parseErrors.length, 0,
                           're-parse errors: ' + m2.parseErrors.join('; '));
        assert.strictEqual(m2.atoms.length, c.atoms,
                           're-parsed atom count drift');
        assert.deepStrictEqual(mapMultiset(m2), c.expectedMaps,
                               're-parsed atom-map multiset mismatch');
    });
});

// ---------------------------------------------------------------------------
// Edge cases on AAM
// ---------------------------------------------------------------------------
test('atom map :0 is treated as no map (bracket atom without :n)', function() {
    var m = P('[CH4:0]');
    // BIME stores 0 as "no map number"; mapNumber should be 0/undefined.
    assert.ok(!m.atoms[0].mapNumber, 'expected no map number for :0');
});

test('large atom map number :999 parses and round-trips', function() {
    var m = P('[CH4:999]');
    assert.strictEqual(m.atoms[0].mapNumber, 999);
    var s = W(m);
    var m2 = P(s);
    assert.strictEqual(m2.atoms[0].mapNumber, 999);
});

test('reaction with two reactants and one product preserves atom count', function() {
    var m = P('[CH3:1][OH:2].[CH3:3][OH:4]>>[CH3:1][O:2][CH3:3].[OH2:4]');
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 8);
    assert.ok(m.reactionArrow);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
