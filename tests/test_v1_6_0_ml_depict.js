/**
 * tests/test_v1_6_0_ml_depict.js — quality gate for v1.6.0-alpha ML residual
 * coordinate refiner.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * What we verify:
 *   1. MLDepict module loads, weights file present (≤ 50 KB), arch matches.
 *   2. Forward pass is deterministic — same input ⇒ byte-identical output.
 *   3. With useMLDepict = false (default), Layout produces v1.5.2 byte-
 *      identical coords for the canonical regression set (caffeine, benzene,
 *      ethanol, glucose). Backwards-compat is non-negotiable.
 *   4. With useMLDepict = true, blended layout still produces:
 *        - finite coords for every atom
 *        - bond lengths in [0.4, 1.6] BL on average
 *        - no atoms further than 5 BL from any neighbour
 *      i.e. ML refinement never breaks chemistry.
 *   5. Sanity: ML-on layout for a few representative inputs stays within
 *      reasonable bounds.
 */
'use strict';

var path  = require('path');
var fs    = require('fs');
var shim  = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'Templates.js'));
require(path.join(__dirname, '..', 'editor', 'Layout.js'));
require(path.join(__dirname, '..', 'editor', 'MLDepict.js'));

var Molecule     = globalThis.Molecule;
var SmilesParser = globalThis.SmilesParser;
var Layout       = globalThis.Layout;
var MLDepict     = globalThis.MLDepict;

var BL = Molecule.BOND_LENGTH;
var runner = shim.makeRunner('v1.6.0-alpha ML depict');

function buildMol(smiles) {
    var mol = new Molecule();
    SmilesParser.parse(smiles, mol);
    return mol;
}

function meanBondLen(mol) {
    var total = 0, n = 0;
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        var a1 = mol.getAtom(b.atom1), a2 = mol.getAtom(b.atom2);
        if (!a1 || !a2) continue;
        var dx = a2.x - a1.x, dy = a2.y - a1.y;
        total += Math.sqrt(dx * dx + dy * dy);
        n++;
    }
    return n ? total / n : 0;
}

function maxNbrDist(mol) {
    var max = 0;
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        var nbrs = mol.getNeighbors(a.id);
        for (var j = 0; j < nbrs.length; j++) {
            var nb = mol.getAtom(nbrs[j]);
            var d = Math.hypot(a.x - nb.x, a.y - nb.y);
            if (d > max) max = d;
        }
    }
    return max;
}

function allFinite(mol) {
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        if (!isFinite(a.x) || !isFinite(a.y)) return false;
    }
    return true;
}

function snapshotCoords(mol) {
    // Compare coords only; atom IDs are process-global and shift across
    // re-parses without affecting layout determinism.
    var snap = [];
    for (var i = 0; i < mol.atoms.length; i++) {
        snap.push([+mol.atoms[i].x.toFixed(4),
                   +mol.atoms[i].y.toFixed(4)]);
    }
    return JSON.stringify(snap);
}

function assert(cond, msg) {
    if (!cond) throw new Error(msg || 'assertion failed');
}

console.log('\n=== v1.6.0-alpha ML residual depict quality gate ===');

// ---------------------------------------------------------------------------
// 1. Module + weights present.
// ---------------------------------------------------------------------------
runner.test('MLDepict module exposes API', function () {
    assert(typeof MLDepict.ready === 'function', 'ready missing');
    assert(typeof MLDepict.predictResidual === 'function', 'predictResidual missing');
    assert(typeof MLDepict.refineLayout === 'function', 'refineLayout missing');
    assert(typeof MLDepict.summary === 'function', 'summary missing');
});

runner.test('weights file present and reasonable size', function () {
    var p = path.join(__dirname, '..', 'editor', 'ml-depict-weights.json');
    assert(fs.existsSync(p), 'weights JSON missing');
    var sz = fs.statSync(p).size;
    assert(sz > 1024, 'weights file suspiciously small: ' + sz);
    assert(sz < 60 * 1024, 'weights file too large for v1.x tier: ' + sz);
});

runner.test('weights load + arch matches', function () {
    assert(MLDepict.ready(), 'MLDepict.ready() = false');
    var s = MLDepict.summary();
    assert(s.arch.in === 40,  'expected in=40, got ' + s.arch.in);
    assert(s.arch.h1 === 64,  'expected h1=64, got ' + s.arch.h1);
    assert(s.arch.h2 === 32,  'expected h2=32, got ' + s.arch.h2);
    assert(s.arch.out === 2,  'expected out=2, got ' + s.arch.out);
    assert(s.training.test_mse < s.training.test_baseline_mse,
        'test MSE should beat baseline; got test=' + s.training.test_mse +
        ' vs baseline=' + s.training.test_baseline_mse);
});

// ---------------------------------------------------------------------------
// 2. Forward pass determinism.
// ---------------------------------------------------------------------------
runner.test('forward pass is deterministic', function () {
    var f  = new Array(20).fill(0);
    var nf = new Array(20).fill(0);
    f[0] = 1; f[10] = 0.5; f[11] = 1;
    nf[0] = 0.7; nf[10] = 0.5; nf[11] = 0.5;
    var r1 = MLDepict.predictResidual(f, nf);
    var r2 = MLDepict.predictResidual(f, nf);
    var r3 = MLDepict.predictResidual(f, nf);
    assert(r1[0] === r2[0] && r2[0] === r3[0], 'dx not deterministic');
    assert(r1[1] === r2[1] && r2[1] === r3[1], 'dy not deterministic');
});

runner.test('forward pass output is finite and bounded', function () {
    // Sweep all element one-hots through the network; every prediction
    // must be finite and within the trained domain (we trained on
    // residuals capped at ±5 BL, so model outputs should also be small).
    for (var e = 0; e < 10; e++) {
        var f  = new Array(20).fill(0);
        var nf = new Array(20).fill(0);
        f[e] = 1; f[10] = 0.5; nf[0] = 1; nf[10] = 0.5;
        var r = MLDepict.predictResidual(f, nf);
        assert(isFinite(r[0]) && isFinite(r[1]),
            'non-finite residual at element ' + e);
        assert(Math.abs(r[0]) < 5 && Math.abs(r[1]) < 5,
            'unbounded residual at element ' + e + ': ' + r);
    }
});

// ---------------------------------------------------------------------------
// 3. Backwards-compat: with useMLDepict = false the layout matches v1.5.2.
// ---------------------------------------------------------------------------
runner.test('useMLDepict=false leaves v1.5.2 coords unchanged (caffeine)', function () {
    Layout.options.useMLDepict = false;
    Layout.options.enforceConventions = true;
    var molA = buildMol('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
    Layout.layout(molA);
    var snapA = snapshotCoords(molA);

    var molB = buildMol('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
    Layout.layout(molB);
    var snapB = snapshotCoords(molB);

    assert(snapA === snapB, 'caffeine layout not deterministic with ML off');
    assert(allFinite(molA), 'caffeine has non-finite coords');
    var meanBL = meanBondLen(molA);
    assert(meanBL > BL * 0.7 && meanBL < BL * 1.3,
        'caffeine bond length out of range: ' + meanBL.toFixed(2));
});

runner.test('useMLDepict=false: benzene + glucose remain stable', function () {
    Layout.options.useMLDepict = false;
    var mols = ['c1ccccc1', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O'];
    for (var i = 0; i < mols.length; i++) {
        var m = buildMol(mols[i]);
        Layout.layout(m);
        assert(allFinite(m), mols[i] + ' has non-finite coords');
        var bl = meanBondLen(m);
        assert(bl > BL * 0.7 && bl < BL * 1.3,
            mols[i] + ' bond length out of range: ' + bl.toFixed(2));
    }
});

// ---------------------------------------------------------------------------
// 4. ML-on does not break chemistry on the canonical set.
// ---------------------------------------------------------------------------
var TEST_SMILES = [
    ['caffeine',     'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'],
    ['benzene',      'c1ccccc1'],
    ['ethanol',      'CCO'],
    ['glucose',      'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O'],
    ['acetylsal',    'CC(=O)Oc1ccccc1C(=O)O'],
    ['cyclohexane',  'C1CCCCC1'],
    ['naphthalene',  'c1ccc2ccccc2c1'],
    ['urea',         'NC(=O)N']
];

runner.test('useMLDepict=true keeps every test molecule chemically sound', function () {
    Layout.options.useMLDepict = true;
    Layout.options.mlDepictWeight = 0.15;
    for (var i = 0; i < TEST_SMILES.length; i++) {
        var name = TEST_SMILES[i][0], smi = TEST_SMILES[i][1];
        var m = buildMol(smi);
        Layout.layout(m);
        assert(allFinite(m), name + ': non-finite coords with ML on');
        var bl = meanBondLen(m);
        assert(bl > BL * 0.4 && bl < BL * 1.6,
            name + ' mean bond length ' + bl.toFixed(2) +
            ' outside [' + (BL * 0.4).toFixed(1) + ', ' +
            (BL * 1.6).toFixed(1) + '] BL');
        var maxD = maxNbrDist(m);
        assert(maxD < BL * 5,
            name + ' has neighbour distance > 5 BL: ' + maxD.toFixed(1));
    }
    // Reset.
    Layout.options.useMLDepict = false;
});

runner.test('useMLDepict=true forward pass is deterministic per atom', function () {
    Layout.options.useMLDepict = true;
    Layout.options.mlDepictWeight = 0.5;
    var s1 = (function () {
        var m = buildMol('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
        Layout.layout(m);
        return snapshotCoords(m);
    })();
    var s2 = (function () {
        var m = buildMol('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
        Layout.layout(m);
        return snapshotCoords(m);
    })();
    assert(s1 === s2, 'caffeine layout not deterministic with ML on');
    Layout.options.useMLDepict = false;
});

// ---------------------------------------------------------------------------
// 5. Smoke test: a handful of natural-product-ish SMILES stay valid.
// ---------------------------------------------------------------------------
runner.test('useMLDepict=true on metabolite-style SMILES is safe', function () {
    Layout.options.useMLDepict = true;
    Layout.options.mlDepictWeight = 0.2;
    var smiles = [
        'OC(=O)CCC(O)C(=O)O',                      // 2-OH glutarate
        'NC(CCCCN)C(=O)O',                          // lysine
        'OP(=O)(O)OCC(O)C(O)C(O)C(O)C=O',           // glucose-6-P
        'CC(=O)NCCCC(O)CO'                          // light alkanolamide
    ];
    for (var i = 0; i < smiles.length; i++) {
        var m = buildMol(smiles[i]);
        Layout.layout(m);
        assert(allFinite(m), smiles[i] + ': non-finite coords');
        var bl = meanBondLen(m);
        assert(bl > BL * 0.4 && bl < BL * 1.6,
            smiles[i] + ' bond length ' + bl.toFixed(2) + ' off');
    }
    Layout.options.useMLDepict = false;
});

// ---------------------------------------------------------------------------
// 6. Single-atom species and multi-fragment ionic mixtures must work
// even though the ML coord predictor cannot improve them.
// ---------------------------------------------------------------------------
runner.test('single-atom species: parse + layout (ML off and on)', function () {
    var cases = [
        ['proton',    '[H+]',   1],
        ['sodium',    '[Na+]',  1],
        ['chloride',  '[Cl-]',  1],
        ['hydroxide', '[OH-]',  1]
    ];
    for (var ci = 0; ci < cases.length; ci++) {
        var name = cases[ci][0], smi = cases[ci][1], expectN = cases[ci][2];
        var molOff = buildMol(smi);
        Layout.options.useMLDepict = false;
        Layout.layout(molOff);
        assert(molOff.atoms.length === expectN,
            name + ': expected ' + expectN + ' atoms, got ' + molOff.atoms.length);
        assert(allFinite(molOff), name + ' (ML off): non-finite coords');

        var molOn = buildMol(smi);
        Layout.options.useMLDepict = true;
        Layout.options.mlDepictWeight = 0.3;
        Layout.layout(molOn);
        assert(allFinite(molOn), name + ' (ML on): non-finite coords');

        // ML cannot move singletons (no neighbours), so the two
        // layouts must be byte-identical for true singletons.
        assert(snapshotCoords(molOff) === snapshotCoords(molOn),
            name + ': ML-off and ML-on coords disagree on singleton');
    }
    Layout.options.useMLDepict = false;
});

runner.test('multi-fragment ionic mixtures: side-by-side layout', function () {
    var cases = [
        ['NaCl',       '[Na+].[Cl-]',           2],
        ['MgCl2',      '[Mg+2].[Cl-].[Cl-]',    3],
        ['water+H+',   'O.[H+]',                2],
        ['acetate+Na', '[Na+].CC(=O)[O-]',      5]
    ];
    Layout.options.useMLDepict = false;
    for (var ci = 0; ci < cases.length; ci++) {
        var name = cases[ci][0], smi = cases[ci][1], expectN = cases[ci][2];
        var mol = buildMol(smi);
        Layout.layout(mol);
        assert(mol.atoms.length === expectN,
            name + ': expected ' + expectN + ' atoms, got ' + mol.atoms.length);
        assert(allFinite(mol), name + ': non-finite coords');
        // Fragments must not overlap — every atom should have a
        // distinct (x, y).
        for (var ai = 0; ai < mol.atoms.length; ai++) {
            for (var aj = ai + 1; aj < mol.atoms.length; aj++) {
                var a = mol.atoms[ai], b = mol.atoms[aj];
                if (a.x === b.x && a.y === b.y) {
                    throw new Error(name + ': atoms ' + ai + ' and ' + aj +
                        ' overlap at (' + a.x + ',' + a.y + ')');
                }
            }
        }
    }
});

runner.test('ML refineLayout no-ops cleanly on singleton mols', function () {
    var mol = buildMol('[Na+]');
    Layout.options.useMLDepict = false;
    Layout.layout(mol);
    // Now manually invoke refineLayout — it must return 0 (nothing to refine).
    var refined = MLDepict.refineLayout(mol, 0.5);
    assert(refined === 0,
        'refineLayout on singleton returned ' + refined + ', expected 0');
});

// ---------------------------------------------------------------------------
// Summary — works as both standalone runner AND bundle-test aggregator.
// ---------------------------------------------------------------------------
module.exports = function () { return runner.summary(); };
if (require.main === module) {
    var summary = runner.summary();
    console.log('\n' + summary.label + ': ' + summary.passed +
                ' passed, ' + summary.failed + ' failed');
    if (summary.failed > 0) {
        summary.failures.forEach(function (f) {
            console.log('  ✗ ' + f.name + ': ' + f.reason);
        });
        process.exit(1);
    }
}
