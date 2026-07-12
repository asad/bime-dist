/**
 * tests/test_v2_4_12_cyclitol_cip.js — CIP cyclitol calibrated abstention +
 * reference-validated stereo coverage (v2.4.12).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * What this guards
 * ----------------
 * The CIP R/S engine (editor/CIPStereo.js) leaves the symmetric ring cyclitols
 * (the inositol family) UNLABELLED: every ring carbon's two ring-path ligands
 * are constitutionally identical, so Rules 1-3/4a tie, and the auxiliary-
 * descriptor pass cannot seed across the ring-closure back-edge. Resolving
 * those would need a STEREO-AWARE canonical ranking (the SmilesWriter ranking
 * is constitutional-only by design — its KNOWN_ORBIT_DRIFT note) plus an
 * external IUPAC-2013 oracle to validate the order-fragile lowercase r/s
 * letters; that is a deferred follow-up. Per the never-wrong contract the engine
 * abstains rather than guess.
 *
 * v2.4.12 makes that abstention HONEST and INTROSPECTABLE: an abstained
 * stereocentre is flagged `atom.cipUndetermined = true` (cipLabel stays ''), so
 * a caller/test can tell an intentional calibrated abstention apart from "not a
 * stereocentre". This suite pins:
 *   A. The acyclic resolver is CORRECT against PubChem — including both meso
 *      hexitols (galactitol, allitol -> mixed R/S, never all-same). This is the
 *      regression armour that proves the shared comparator is sound before any
 *      future ring extension reuses it.
 *   B. Acyclic pseudo-asymmetric centres still resolve (ribitol r, ribaric s)
 *      with zero undetermined flags.
 *   C. Ring cyclitols (myo/scyllo/cis-inositol) abstain on EVERY ring centre:
 *      no fabricated descriptor, and each specified centre flagged
 *      cipUndetermined — order-invariant across equivalent SMILES.
 *   D. The marker is precise: resolved centres and plain non-stereocentres are
 *      never flagged.
 *   E. Source-shape pins on editor/CIPStereo.js.
 *
 * Reference ground truth (prof-chem / IUPAC 2013 P-91/P-92; recorded here for
 * the future stereo-aware follow-up, NOT asserted as engine output):
 *   - scyllo-inositol & cis-inositol: every carbon is genuinely NON-stereogenic
 *     (achiral molecule) -> "unlabelled" is the chemically correct answer.
 *   - myo-inositol: four true R/S (reflection-paired, opposite) + two pseudo-
 *     asymmetric lowercase r/s at C2 and C5 (the mirror-plane carbons).
 * PubChem CIDs for the hexitol gold values: galactitol 11850, allitol 643385,
 * D-mannitol 6251, D-glucitol/sorbitol 5780.
 *
 * Pure Node (real engine via shim.loadAll), no DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var SmilesParser = global.SmilesParser;
var CIPStereo = global.CIPStereo;

var runner = shim.makeRunner('CIP cyclitol calibrated abstention (v2.4.12)');
var test = runner.test;
console.log('CIP cyclitol calibrated abstention (v2.4.12)');

var SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'CIPStereo.js'), 'utf8');

// Parse + assign, then summarise the specified tetrahedral stereocentres.
function summary(smi) {
    var m = SmilesParser.parse(smi);
    assert.ok(m && m.atoms.length, 'parsed ' + smi);
    CIPStereo.assign(m);
    var labels = [], undetermined = 0, specified = 0, strayFlag = 0;
    for (var i = 0; i < m.atoms.length; i++) {
        var a = m.atoms[i];
        if (a.chirality) {
            specified++;
            if (a.cipLabel) labels.push(a.cipLabel);
            if (a.cipUndetermined) undetermined++;
        } else if (a.cipUndetermined) {
            // A non-specified atom must never carry the abstention flag.
            strayFlag++;
        }
    }
    labels.sort();
    return { labels: labels, multiset: labels.join(','), undetermined: undetermined,
             specified: specified, strayFlag: strayFlag };
}

// ---------------------------------------------------------------------------
// A. The acyclic resolver matches PubChem — including the meso hexitols.
//    (Refutes the "meso hexitol is mislabelled" false alarm and guards it.)
// ---------------------------------------------------------------------------
var HEXITOLS = [
    // name, PubChem isomeric SMILES, expected multiset, PubChem IUPAC config
    ['galactitol (meso)', 'C([C@H]([C@@H]([C@@H]([C@H](CO)O)O)O)O)O', 'R,R,S,S', '(2R,3S,4R,5S)'],
    ['allitol (meso)',    'C([C@H]([C@H]([C@H]([C@H](CO)O)O)O)O)O',   'R,R,S,S', '(2R,3R,4S,5S)'],
    ['D-mannitol',        'C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O', 'R,R,R,R', '(2R,3R,4R,5R)'],
    ['D-glucitol',        'C([C@H]([C@H]([C@@H]([C@H](CO)O)O)O)O)O',  'R,R,R,S', '(2R,3R,4R,5S)']
];
HEXITOLS.forEach(function (row) {
    test('A:' + row[0] + ' multiset = ' + row[2] + ' ' + row[3] + ' (PubChem)', function () {
        var s = summary(row[1]);
        assert.strictEqual(s.multiset, row[2],
            row[0] + ' must be ' + row[2] + ' to match PubChem ' + row[3] + '; got ' + s.multiset);
        assert.strictEqual(s.undetermined, 0, row[0] + ' fully resolves (no abstention)');
    });
});
test('A5 the meso hexitols are MIXED R/S, never all-same (a meso compound cannot be all-R/all-S)', function () {
    assert.strictEqual(summary(HEXITOLS[0][1]).multiset, 'R,R,S,S', 'galactitol is meso -> 2R+2S');
    assert.strictEqual(summary(HEXITOLS[1][1]).multiset, 'R,R,S,S', 'allitol is meso -> 2R+2S');
});

// ---------------------------------------------------------------------------
// B. Acyclic pseudo-asymmetric centres still resolve (no regression), 0 flags.
// ---------------------------------------------------------------------------
test('B1 ribitol stays R,S,r and ribaric acid stays R,S,s with zero abstention', function () {
    var ribitol = summary('OC[C@@H](O)[C@H](O)[C@@H](O)CO');
    var ribaric = summary('O=C(O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O');
    assert.strictEqual(ribitol.multiset, 'R,S,r', 'ribitol multiset R,S,r');
    assert.strictEqual(ribaric.multiset, 'R,S,s', 'ribaric acid multiset R,S,s');
    assert.strictEqual(ribitol.undetermined, 0, 'ribitol fully determined');
    assert.strictEqual(ribaric.undetermined, 0, 'ribaric acid fully determined');
});

// ---------------------------------------------------------------------------
// C. Ring cyclitols abstain on every ring centre — never fabricated, flagged,
//    order-invariant. (myo/scyllo/cis-inositol.)
// ---------------------------------------------------------------------------
var INOSITOLS = {
    'myo-inositol': [
        'O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O',
        'O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    ],
    'scyllo-inositol': [
        'O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O'
    ],
    'cis-inositol': [
        'O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O'
    ]
};
Object.keys(INOSITOLS).forEach(function (name) {
    test('C:' + name + ' abstains on all 6 ring centres (no fabricated R/S, all flagged)', function () {
        INOSITOLS[name].forEach(function (smi) {
            var s = summary(smi);
            assert.strictEqual(s.specified, 6, name + ' has 6 specified ring stereocentres (' + smi + ')');
            assert.strictEqual(s.labels.length, 0,
                name + ' must emit NO descriptor (never-wrong); got [' + s.multiset + '] for ' + smi);
            assert.strictEqual(s.undetermined, 6,
                name + ' flags every abstained centre cipUndetermined; got ' + s.undetermined + ' for ' + smi);
        });
    });
});
test('C-order myo-inositol abstention is identical across two equivalent SMILES', function () {
    var a = summary(INOSITOLS['myo-inositol'][0]);
    var b = summary(INOSITOLS['myo-inositol'][1]);
    assert.strictEqual(a.labels.length, 0); assert.strictEqual(b.labels.length, 0);
    assert.strictEqual(a.undetermined, 6); assert.strictEqual(b.undetermined, 6);
});

// ---------------------------------------------------------------------------
// D. The marker is precise — never set on a resolved centre or a non-centre.
// ---------------------------------------------------------------------------
test('D1 resolved stereocentres and plain non-stereocentres are never flagged', function () {
    // D-glucose: every centre resolves -> no flags anywhere.
    var glu = summary('O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO');
    assert.strictEqual(glu.undetermined, 0, 'glucose: resolved centres carry no abstention flag');
    assert.strictEqual(glu.strayFlag, 0, 'glucose: no non-specified atom is flagged');
    // Ethanol: no stereocentre at all -> nothing flagged.
    var eth = summary('CCO');
    assert.strictEqual(eth.specified, 0, 'ethanol has no specified stereocentre');
    assert.strictEqual(eth.strayFlag, 0, 'ethanol: no atom carries cipUndetermined');
});
test('D2 inositol carries the flag ONLY on its specified ring carbons, not its O/H', function () {
    var s = summary('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O');
    assert.strictEqual(s.strayFlag, 0, 'no oxygen / non-centre is flagged undetermined');
    assert.strictEqual(s.undetermined, 6, 'exactly the six ring carbons are flagged');
});

// ---------------------------------------------------------------------------
// E. Source-shape pins on editor/CIPStereo.js.
// ---------------------------------------------------------------------------
test('E1 the abstention path sets the honest cipUndetermined flag', function () {
    assert.ok(/atom\.cipUndetermined = true;/.test(SRC),
        'the !priorities abstention branch records atom.cipUndetermined = true');
    assert.ok(/atom\.cipUndetermined = false;/.test(SRC),
        'cipUndetermined is reset per pass (no stale flag across re-assignment)');
});
test('E2 the calibrated-abstention contract + cyclitol ground truth are documented', function () {
    assert.ok(/Calibrated abstention/.test(SRC), 'header documents the never-wrong calibrated abstention');
    assert.ok(/cyclitol|inositol/i.test(SRC), 'header names the cyclitol abstention class');
    assert.ok(/non-stereogenic/.test(SRC) && /pseudo-\s*\n?\s*\*?\s*asymmetric|pseudo-asymmetric/i.test(SRC),
        'header records the scyllo/cis non-stereogenic + myo pseudo-asymmetric ground truth');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
