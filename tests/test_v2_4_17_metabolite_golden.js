/**
 * tests/test_v2_4_17_metabolite_golden.js — golden-structure assertions for the
 * teaching-critical phosphosugar / glycolysis / PPP metabolites (v2.4.17).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists: before v2.4.17 the MetaboliteLibrary tests only checked
 * SELF-consistency (findBySmiles(entry) === entry). That is satisfied even when
 * an entry's structure is wrong — it just resolves the wrong structure to its
 * own (wrong) name. Six entries shipped with wrong structures under correct
 * names (a teaching hazard): Glucose 6-phosphate was drawn as glucose-4-P;
 * ribulose-5-P / xylulose-5-P had the ketone adjacent to the phosphate (a
 * 1-phosphate connectivity); ribose-5-P was the C3-epimer (xylose-5-P); and
 * 3-phosphoglycerate / 1,3-BPG were the L (unnatural) enantiomer.
 *
 * This test pins each corrected entry to an AUTHORITATIVE fingerprint — the
 * molecular formula (Hill) and the sorted CIP R/S multiset read by the (sound)
 * CIP reader — cross-checked against the authoritative compound record. It deliberately
 * verifies through the PARSER + CIP READER, never the SMILES writer.
 *
 * Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'CIPStereo.js'));
require(path.join(__dirname, '..', 'editor', 'MetaboliteLibrary.js'));

var SmilesParser = globalThis.SmilesParser;
var CIPStereo = globalThis.CIPStereo;
var Molecule = globalThis.Molecule;
var MetaboliteLibrary = globalThis.MetaboliteLibrary;

var runner = shim.makeRunner('Metabolite golden structures (v2.4.17)');
var test = runner.test;
console.log('Metabolite golden structures (v2.4.17)');

function formulaOf(m) {
    var counts = (typeof m.elementCounts === 'function') ? m.elementCounts() : null;
    return counts ? Molecule.formulaString(counts) : '?';
}
function cipMultiset(m) {
    CIPStereo.assignRS(m);
    return m.atoms.map(function (a) { return a.cipLabel; })
        .filter(Boolean).sort().join(',');
}

// name in library -> { kegg, cid, formula, cip (sorted R/S multiset) }
// cip/formula independently verified against the authoritative compound record by prof-chem AND
// re-derived here from the library's own stored SMILES via the CIP reader.
var GOLD = {
    'Glucose 6-phosphate':     { kegg: 'C00668', cid: 5958,    formula: 'C6H13O9P',  cip: 'R,R,S,S' },
    '1,3-Bisphosphoglycerate': { kegg: 'C00236', cid: 683,     formula: 'C3H8O10P2', cip: 'R' },
    '3-Phosphoglycerate':      { kegg: 'C00197', cid: 724,     formula: 'C3H7O7P',   cip: 'R' },
    'D-Ribulose 5-phosphate':  { kegg: 'C00199', cid: 439184,  formula: 'C5H11O8P',  cip: 'R,R' },
    'D-Ribose 5-phosphate':    { kegg: 'C00117', cid: 77982,   formula: 'C5H11O8P',  cip: 'R,R,R' },
    'D-Xylulose 5-phosphate':  { kegg: 'C00231', cid: 5459919, formula: 'C5H11O8P',  cip: 'R,S' }
};

Object.keys(GOLD).forEach(function (name) {
    var gold = GOLD[name];
    test(name + ' has the authoritative formula + CIP (' + gold.kegg + ' / CID ' + gold.cid + ')', function () {
        var entry = MetaboliteLibrary.find(name);
        assert.ok(entry && entry.smiles, name + ' is in the library with a SMILES');
        assert.strictEqual(entry.kegg, gold.kegg, name + ' KEGG id');
        var m = SmilesParser.parse(entry.smiles);
        assert.strictEqual((m.parseErrors || []).length, 0, name + ' parses cleanly');
        assert.strictEqual(formulaOf(m), gold.formula, name + ' molecular formula');
        assert.strictEqual(cipMultiset(m), gold.cip, name + ' CIP R/S multiset (absolute config)');
    });
});

// Explicit anti-regression: the specific wrong structures that shipped must
// NEVER come back. These are the exact pre-v2.4.17 (incorrect) SMILES.
var FORBIDDEN = {
    'Glucose 6-phosphate':     'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1OP(O)(O)=O', // phosphate on ring C4
    'D-Ribose 5-phosphate':    'O=C[C@H](O)[C@@H](O)[C@H](O)COP(O)(O)=O',             // C3-epimer (xylose-5-P)
    '3-Phosphoglycerate':      'O[C@@H](COP(O)(O)=O)C(O)=O'                            // L enantiomer
};
test('the previously-shipped WRONG structures are not present', function () {
    Object.keys(FORBIDDEN).forEach(function (name) {
        var entry = MetaboliteLibrary.find(name);
        assert.notStrictEqual(entry.smiles, FORBIDDEN[name],
            name + ' must not carry its old (incorrect) SMILES again');
    });
});

// Sanity anchor: the trusted D-Glucose entry the G6P skeleton is built from.
test('D-Glucose anchor still reads (2R,3S,4R,5R) ring multiset', function () {
    var g = MetaboliteLibrary.find('D-Glucose');
    var m = SmilesParser.parse(g.smiles);
    assert.strictEqual(formulaOf(m), 'C6H12O6', 'D-Glucose formula');
    assert.strictEqual(cipMultiset(m), 'R,R,S,S', 'D-Glucose ring CIP multiset');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
