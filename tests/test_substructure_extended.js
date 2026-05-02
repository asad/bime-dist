/**
 * tests/test_substructure_extended.js — extended SUB regressions for v1.1.1.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Adds drug + functional-group + recursive-SMARTS + disconnected-target +
 * ring-constraint + negative-case + multi-mapping + performance test cases for
 * SmartsMatch.match against real molecules.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SUB (SMARTS extended)');
var test = runner.test;

console.log('SUB (SMARTS extended)');

function P(smi) { return SmilesParser.parse(smi); }
function Q(sma) { return SmartsParser.parse(sma); }
function nMatches(targetSmi, smarts) {
    return SmartsMatch.match(P(targetSmi), Q(smarts)).length;
}

// =========================================================================
// (1) Phenyl substructure in toluene / ethylbenzene / diphenylmethane
// Note: SmartsMatch returns symmetry-equivalent mappings. Each benzene admits
// 12 (6 rotations × 2 reflections), so we assert >= 12 for one ring and 24 for
// two rings.
// =========================================================================
test('phenyl c1ccccc1 in toluene → >=12 matches (one benzene ring with symmetry)', function() {
    assert.ok(nMatches('Cc1ccccc1', 'c1ccccc1') >= 12);
});
test('phenyl c1ccccc1 in ethylbenzene → >=12 matches', function() {
    assert.ok(nMatches('CCc1ccccc1', 'c1ccccc1') >= 12);
});
test('phenyl c1ccccc1 in diphenylmethane → 24 matches (two rings)', function() {
    assert.strictEqual(nMatches('c1ccc(Cc2ccccc2)cc1', 'c1ccccc1'), 24);
});

// =========================================================================
// (2) Carboxylic acid C(=O)O in NSAIDs
// =========================================================================
test('carboxylic acid C(=O)O in aspirin → 2 matches (acid + ester)', function() {
    // C(=O)O matches both the carboxylic acid AND the acetate ester
    // because the SMARTS is loose. This documents BIME's permissive policy.
    assert.strictEqual(nMatches('CC(=O)Oc1ccccc1C(=O)O', 'C(=O)O'), 2);
});
test('carboxylic acid C(=O)O in ibuprofen → 1 match', function() {
    assert.strictEqual(nMatches('CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'C(=O)O'), 1);
});
test('carboxylic acid C(=O)O in naproxen → 1 match', function() {
    assert.strictEqual(nMatches('COc1ccc2cc(C(C)C(=O)O)ccc2c1', 'C(=O)O'), 1);
});

// =========================================================================
// (3) Amide bond C(=O)N in paracetamol / phenacetin
// =========================================================================
test('amide C(=O)N in paracetamol → 1 match', function() {
    assert.strictEqual(nMatches('CC(=O)Nc1ccc(O)cc1', 'C(=O)N'), 1);
});
test('amide C(=O)N in phenacetin → 1 match', function() {
    assert.strictEqual(nMatches('CCOc1ccc(NC(C)=O)cc1', 'C(=O)N'), 1);
});

// =========================================================================
// (4) Aromatic-N pattern [n] in heterocycles
// =========================================================================
test('aromatic N [n] in pyridine → 1 match', function() {
    assert.strictEqual(nMatches('c1ccncc1', '[n]'), 1);
});
test('aromatic N [n] in imidazole → 2 matches', function() {
    assert.strictEqual(nMatches('c1cnc[nH]1', '[n]'), 2);
});
test('aromatic N [n] in pyrimidine → 2 matches', function() {
    assert.strictEqual(nMatches('c1ccncn1', '[n]'), 2);
});

// =========================================================================
// (5) Sulfonamide pattern S(=O)(=O)N
// =========================================================================
test('sulfonamide S(=O)(=O)N in CS(=O)(=O)N → >=1 match', function() {
    // Match count is 2 due to the symmetric S=O dual encoding.
    assert.ok(nMatches('CS(=O)(=O)N', 'S(=O)(=O)N') >= 1);
});

// =========================================================================
// (6) Recursive SMARTS [$([CX3]=[OX1])] — carbonyl
// =========================================================================
test('recursive SMARTS [$([CX3]=[OX1])] in acetone → 1 match (the C=O carbon)', function() {
    assert.strictEqual(nMatches('CC(C)=O', '[$([CX3]=[OX1])]'), 1);
});

// =========================================================================
// (7) SMARTS in a disconnected target — only the matching component returns hits
// =========================================================================
test('phenyl c1ccccc1 in (benzene . ethane) target → matches only the benzene', function() {
    // SmartsMatch returns mappings only for atoms that actually match.
    // Disconnected ethane component yields zero, so total = 12 (one ring symmetry set).
    assert.strictEqual(nMatches('c1ccccc1.CC', 'c1ccccc1'), 12);
});
test('phenyl c1ccccc1 in pure-aliphatic disconnected target → 0 matches', function() {
    assert.strictEqual(nMatches('CC.CC', 'c1ccccc1'), 0);
});

// =========================================================================
// (8) Ring constraint [R]
// =========================================================================
test('[R] in benzene → 6 matches (every C is a ring atom)', function() {
    assert.strictEqual(nMatches('c1ccccc1', '[R]'), 6);
});
test('[R] in ethane → 0 matches (no ring atoms)', function() {
    assert.strictEqual(nMatches('CC', '[R]'), 0);
});

// =========================================================================
// (9) Negative case: nitrogen query against pure-carbon target
// =========================================================================
test('N in benzene → 0 matches', function() {
    assert.strictEqual(nMatches('c1ccccc1', 'N'), 0);
});

// =========================================================================
// (10) Multi-mapping: [#6] in cyclohexane → 6 matches (one per atom)
// =========================================================================
test('[#6] in cyclohexane → 6 matches', function() {
    assert.strictEqual(nMatches('C1CCCCC1', '[#6]'), 6);
});

// =========================================================================
// (11) Performance: substructure search on a 50-atom molecule completes <100ms
// =========================================================================
test('CCC in 50-atom alkane completes in < 100ms', function() {
    var t0 = Date.now();
    var bigSmi = '';
    for (var i = 0; i < 50; i++) bigSmi += 'C';
    var n = nMatches(bigSmi, 'CCC');
    var elapsed = Date.now() - t0;
    assert.ok(n > 0, 'expected >0 matches');
    assert.ok(elapsed < 100, 'expected <100ms; got ' + elapsed + 'ms');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
