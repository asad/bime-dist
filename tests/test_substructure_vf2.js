/**
 * tests/test_substructure_vf2.js — SMSDVF2 substructure regressions.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('SMSD VF2');
var test = runner.test;

console.log('SMSD VF2 substructure');

function G(smi) {
    var mol = SmilesParser.parse(smi);
    return new SMSDGraph.SMSDGraph(mol);
}
function opts() { return new SMSDGraph.ChemOptions(); }

// --- Self-match ---
test('benzene contains benzene (self-match)', function() {
    var bz = G('c1ccccc1');
    var bz2 = G('c1ccccc1');
    var maps = SMSDVF2.findAllSubstructures(bz, bz2, opts(), 100);
    assert.ok(maps.length >= 1, 'expected ≥1 mapping, got ' + maps.length);
});

// --- Subgraph ---
test('CC is a substructure of CCC', function() {
    var cc = G('CC');
    var ccc = G('CCC');
    assert.strictEqual(SMSDVF2.isSubstructure(cc, ccc, opts()), true);
});

test('findAllSubstructures CC in CCC returns ≥1 mapping', function() {
    var cc = G('CC');
    var ccc = G('CCC');
    var maps = SMSDVF2.findAllSubstructures(cc, ccc, opts(), 100);
    assert.ok(maps.length >= 1, 'expected ≥1, got ' + maps.length);
});

// --- Empty query (B-3 documented behavior) ---
test('empty-query findSubstructure returns the empty mapping {}', function() {
    var emptyG = new SMSDGraph.SMSDGraph(new Molecule());
    var bz = G('c1ccccc1');
    var r = SMSDVF2.findSubstructure(emptyG, bz, opts());
    // BIME contract: empty query is trivially a substructure of anything
    // (returns {} not null). See SMSDVF2.findSubstructure line 733.
    assert.deepStrictEqual(r, {});
});

test('empty-query isSubstructure returns true', function() {
    var emptyG = new SMSDGraph.SMSDGraph(new Molecule());
    var bz = G('c1ccccc1');
    assert.strictEqual(SMSDVF2.isSubstructure(emptyG, bz, opts()), true);
});

// --- Disconnected molecule ---
test('O.O in O.O — at least one mapping found', function() {
    var oo = G('O.O');
    var oo2 = G('O.O');
    var maps = SMSDVF2.findAllSubstructures(oo, oo2, opts(), 100);
    assert.ok(maps.length >= 1, 'expected ≥1 mapping, got ' + maps.length);
});

// --- Isolated atom query ---
test('single C in benzene — 6 placements', function() {
    var single = G('C');
    var bz = G('c1ccccc1');
    var maps = SMSDVF2.findAllSubstructures(single, bz, opts(), 100);
    // 6 carbons × any orientation, but the symmetric C atom maps to each → ≥6
    assert.ok(maps.length >= 6, 'expected ≥6 mappings, got ' + maps.length);
});

// --- No-match cases ---
test('N (alone) not a substructure of benzene (only C atoms)', function() {
    var n = G('N');
    var bz = G('c1ccccc1');
    assert.strictEqual(SMSDVF2.isSubstructure(n, bz, opts()), false);
});

test('benzene not a substructure of cyclohexane (aromaticity differs)', function() {
    var bz = G('c1ccccc1');
    var cy = G('C1CCCCC1');
    // BIME's atom-type matching treats aromatic C ≠ aliphatic C by default
    // (matchAtomType=true). This locks in that policy.
    assert.strictEqual(SMSDVF2.isSubstructure(bz, cy, opts()), false);
});

test('larger query (CCC) cannot match smaller target (CC) → null', function() {
    var ccc = G('CCC');
    var cc = G('CC');
    assert.strictEqual(SMSDVF2.findSubstructure(ccc, cc, opts()), null);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
