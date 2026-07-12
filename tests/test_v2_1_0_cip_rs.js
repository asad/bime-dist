/**
 * tests/test_v2_1_0_cip_rs.js — metabolite stereochemistry data corrections
 * (v2.1.0) + a documented CIP-engine limitation.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists
 * ---------------
 * A prof-chem audit for v2.1.0 found TWO separate things:
 *
 *  (1) Six metabolites in editor/MetaboliteLibrary.js carried the WRONG
 *      enantiomer for their name. These are objective configuration facts
 *      (independent of any software): an L-α-amino acid is (S) at the
 *      α-carbon; D-glyceraldehyde-3-phosphate is (R); the β-oxidation
 *      L-3-hydroxyacyl-CoA intermediate is (S). The corrected SMILES tokens
 *      are pinned below so the data cannot silently regress.
 *
 *  (2) The CIP R/S DISPLAY engine (editor/CIPStereo.js) was itself
 *      input-order-dependent in v2.1.0 — the same molecule written two ways
 *      could get opposite descriptors. That engine bug was FIXED in v2.1.1
 *      (reference-validated; see tests/test_v2_1_1_cip_engine.js). This suite is
 *      retained as the DATA regression guard for the corrected tokens — it
 *      pins the structural SMILES, which is orthogonal to the engine fix.
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/MetaboliteLibrary.js');
var ML = global.MetaboliteLibrary;

var runner = shim.makeRunner('Metabolite stereochemistry (v2.1.0)');
var test = runner.test;

console.log('Metabolite stereochemistry (v2.1.0)');

// The corrected SMILES token for each metabolite whose stored configuration
// was the wrong enantiomer before v2.1.0, with the objective configuration
// it encodes (prof-chem-verified; textbook for the amino acids and D-GAP).
var CORRECTED = [
    { name: 'L-Ornithine',            smiles: 'NCCC[C@H](N)C(O)=O',              config: '(S) α-carbon (L-amino acid)' },
    { name: 'L-Citrulline',           smiles: 'NC(=O)NCCC[C@H](N)C(O)=O',        config: '(S) α-carbon (L-amino acid)' },
    { name: 'L-Arginine',             smiles: 'N=C(N)NCCC[C@H](N)C(O)=O',        config: '(S) α-carbon (L-amino acid)' },
    { name: 'L-Argininosuccinate',    smiles: 'N=C(NCCC[C@H](N)C(O)=O)N[C@@H](CC(O)=O)C(O)=O', config: '(S,S) both α-carbons (L)' },
    { name: 'Glyceraldehyde 3-phosphate', smiles: 'O=C[C@H](O)COP(O)(O)=O',      config: '(R) — D-GAP' },
    { name: 'L-3-Hydroxybutyryl-CoA', smiles: 'C[C@H](O)CC(=O)SC',               config: '(S) — L / β-oxidation intermediate' }
];

// ---------------------------------------------------------------------------
// A. The corrected stereo tokens are the ones shipped (data regression guard).
// ---------------------------------------------------------------------------
CORRECTED.forEach(function (row) {
    test('A:' + row.name + ' has the corrected ' + row.config + ' token',
        function () {
            var m = ML.find(row.name);
            assert.ok(m, row.name + ' must be in the library');
            assert.strictEqual(m.smiles, row.smiles,
                row.name + ' SMILES must be the corrected token ' + row.smiles
                    + ' (encoding ' + row.config + ')');
        });
});

// ---------------------------------------------------------------------------
// B. Each corrected token still parses cleanly and keeps exactly the
// stereocentre(s) it should (structural sanity — the fix flipped a token,
// it did not drop or add a centre).
// ---------------------------------------------------------------------------
test('B1 each corrected metabolite parses with its expected stereocentres',
    function () {
        var expectCentres = {
            'L-Ornithine': 1, 'L-Citrulline': 1, 'L-Arginine': 1,
            'L-Argininosuccinate': 2, 'Glyceraldehyde 3-phosphate': 1,
            'L-3-Hydroxybutyryl-CoA': 1
        };
        CORRECTED.forEach(function (row) {
            var mol = SmilesParser.parse(ML.find(row.name).smiles);
            assert.strictEqual(mol.parseErrors.length, 0,
                row.name + ' must parse cleanly');
            var n = mol.atoms.filter(function (a) { return a.chirality; }).length;
            assert.strictEqual(n, expectCentres[row.name],
                row.name + ' should have ' + expectCentres[row.name]
                    + ' stereocentre(s)');
        });
    });

// ---------------------------------------------------------------------------
// C. Regression guard for the v2.0.76 ketone body (already correct).
// ---------------------------------------------------------------------------
test('C1 D-β-Hydroxybutyrate keeps its (R) token', function () {
    assert.strictEqual(ML.find('D-β-Hydroxybutyrate').smiles,
        'C[C@@H](O)CC(O)=O', 'the (R) ketone body token is unchanged');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
