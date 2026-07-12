/**
 * tests/test_v2_4_17_writer_stereo.js — the SMILES writer must preserve absolute
 * configuration through parse -> write -> parse (v2.4.17).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background: SmilesWriter._resolveChirality built its "stored" neighbour frame
 * by HOISTING the from-atom to the front, but SmilesParser normalises
 * atom.chirality into mol.getNeighbors() adjacency order with NO hoist. The
 * mismatch inverted @/@@ for ~26% of chiral molecules whose stereocentre
 * neighbours are reordered by the DFS (ring / adjacent centres) — e.g. it
 * exported (S)-naproxen as (R), R-limonene as (S), and (R,R)-tartaric as meso.
 * The previous round-trip test (test_round_trip.js KNOWN_ORBIT_DRIFT) only
 * asserted atom/bond COUNT survival, so the epimerisation went uncaught.
 *
 * A second facet: RDT.extractComponents rebuilds a component's adjacency in
 * sorted-id order, which permutes stereocentre neighbours; it now re-normalises
 * chirality so reaction components (used for figure compound labels and RXN
 * export) keep their configuration too.
 *
 * These tests pin CIP-multiset preservation, verified through the (sound) CIP
 * reader, for both single-molecule and reaction-component write paths.
 *
 * Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'RDT.js'));
require(path.join(__dirname, '..', 'editor', 'CIPStereo.js'));
require(path.join(__dirname, '..', 'editor', 'MetaboliteLibrary.js'));

var SmilesParser = globalThis.SmilesParser;
var SmilesWriter = globalThis.SmilesWriter;
var CIPStereo = globalThis.CIPStereo;
var RDT = globalThis.RDT;
var MetaboliteLibrary = globalThis.MetaboliteLibrary;

var runner = shim.makeRunner('Writer stereo preservation (v2.4.17)');
var test = runner.test;
console.log('Writer stereo preservation (v2.4.17)');

function cipMS(m) { CIPStereo.assignRS(m); return m.atoms.map(function (a) { return a.cipLabel; }).filter(Boolean).sort().join(','); }
function roundtripMS(smi) { return cipMS(SmilesParser.parse(SmilesWriter.write(SmilesParser.parse(smi)))); }

// Molecules that USED to invert on write (single stereocentre => a flip is a
// genuine enantiomer inversion, not benign symmetric-orbit drift).
var WAS_FLIPPING = {
    'R-limonene':       'CC(=C)[C@@H]1CCC(=CC1)C',
    'naproxen':         'C[C@H](C(=O)O)c1ccc2cc(OC)ccc2c1',
    'L-tartaric':       'O[C@H]([C@@H](O)C(=O)O)C(=O)O',
    'quinine':          'C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@@H](O)c1ccnc2ccc(OC)cc12',
    'D-glyceraldehyde-3-P': 'O=C[C@H](O)COP(O)(O)=O'
};
// Controls that were always correct — must STAY correct.
var CONTROLS = {
    'L-alanine':        'N[C@@H](C)C(O)=O',
    '(R)-2-butanol':    'C[C@@H](O)CC',
    'D-glucose':        'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
    'L-malate':         'OC(=O)[C@@H](O)CC(O)=O'
};

test('A previously-inverting molecules now preserve their CIP configuration', function () {
    Object.keys(WAS_FLIPPING).forEach(function (name) {
        var before = cipMS(SmilesParser.parse(WAS_FLIPPING[name]));
        var after = roundtripMS(WAS_FLIPPING[name]);
        assert.ok(before.length > 0, name + ' has stereocentres to check');
        assert.strictEqual(after, before, name + ' CIP multiset must survive parse->write->parse');
    });
});

test('B control chiral molecules stay correct (no new inversions)', function () {
    Object.keys(CONTROLS).forEach(function (name) {
        assert.strictEqual(roundtripMS(CONTROLS[name]), cipMS(SmilesParser.parse(CONTROLS[name])),
            name + ' must remain configuration-stable');
    });
});

test('C the explicit naproxen/limonene enantiomer identity is exact', function () {
    // (S)-naproxen must not become (R); R-limonene must not become (S).
    assert.strictEqual(roundtripMS('C[C@H](C(=O)O)c1ccc2cc(OC)ccc2c1'), 'S', '(S)-naproxen stays S');
    assert.strictEqual(roundtripMS('CC(=C)[C@@H]1CCC(=CC1)C'), 'R', 'R-limonene stays R');
});

test('D reaction components keep their configuration (RDT.extractComponents re-normalises)', function () {
    // The glucose reactant of glucose -> G6P, written from the RDT component,
    // must still resolve to D-Glucose (its CIP multiset must be intact).
    var GLU = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O';
    var G6P = 'O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O';
    var res = RDT.mapReaction(SmilesParser.parse(GLU + '>>' + G6P), {});
    var comp = res.sides.reactants[0];
    var written = SmilesWriter.write(comp);
    assert.strictEqual(cipMS(SmilesParser.parse(written)), 'R,R,S,S',
        'reaction-component glucose keeps (R,R,S,S)');
    var entry = MetaboliteLibrary.findBySmiles(written);
    assert.ok(entry && entry.name === 'D-Glucose',
        'reaction-component glucose resolves to D-Glucose (got ' + (entry ? entry.name : 'null') + ')');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
