/**
 * tests/test_chem_torture.js — Diverse polycyclic / charged / heterocyclic
 * regression coverage for BIME v1.0.2.
 *
 * Each case is a publicly-known textbook structure (PubChem / ChEMBL /
 * IUPAC-style references). The expected counts are derived from BIME's
 * v1.0.2 SmilesParser policy and serve as regression anchors:
 *   - polycyclic aromatic perception (anthracene, phenanthrene, pyrene,
 *     azulene)
 *   - charged species (zwitterions, multi-charge ions, aromatic cations)
 *   - sugars / nucleotides (ATP-like fragment, glucopyranose, glycoside
 *     dimer)
 *   - bridged / spiro / cage systems (norbornane, adamantane, cubane)
 *   - macrocycles (12-membered lactone, cyclam)
 *   - bio-relevant motifs (peptide bond, sulfonamide, phosphate,
 *     guanidinium, boronic acid)
 *
 * No file I/O, no network, runs in <50 ms.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Chem Torture');
var test = runner.test;

console.log('Chem Torture');

function P(smi) { return SmilesParser.parse(smi); }
function aromatic(m) { return m.atoms.filter(function(a) { return a.aromatic; }).length; }
function totalCharge(m) { return m.atoms.reduce(function(s, a) { return s + a.charge; }, 0); }

// ---------------------------------------------------------------------------
// Polycyclic aromatics (Hückel 4n+2 across fused rings)
// ---------------------------------------------------------------------------
[
    { name: 'anthracene',   smi: 'c1ccc2cc3ccccc3cc2c1',     atoms: 14, bonds: 16, arom: 14 },
    { name: 'phenanthrene', smi: 'c1ccc2ccc3ccccc3c2c1',     atoms: 14, bonds: 16, arom: 14 },
    { name: 'pyrene',       smi: 'c1cc2ccc3cccc4ccc(c1)c2c34', atoms: 16, bonds: 19, arom: 16 },
    { name: 'azulene',      smi: 'c1cc2cccccc2c1',           atoms: 10, bonds: 11, arom: 10 },
    { name: 'phenanthroline', smi: 'c1ccc2ncc3ncccc3c2c1',   atoms: 14, bonds: 16, arom: 14 }
].forEach(function(c) {
    test(c.name + ': polycyclic aromatic atom/bond/aromatic counts', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0, c.name + ' has parse errors');
        assert.strictEqual(m.atoms.length, c.atoms);
        assert.strictEqual(m.bonds.length, c.bonds);
        assert.strictEqual(aromatic(m), c.arom);
    });
});

// ---------------------------------------------------------------------------
// Heteroaromatic 5- and 6-membered rings (donor vs acceptor N)
// ---------------------------------------------------------------------------
[
    { name: 'pyridazine',    smi: 'c1ccnnc1',     atoms: 6, arom: 6 },
    { name: 'pyrimidine',    smi: 'c1cncnc1',     atoms: 6, arom: 6 },
    { name: 'pyrazine',      smi: 'c1cnccn1',     atoms: 6, arom: 6 },
    { name: 'pyrazole',      smi: 'c1cc[nH]n1',   atoms: 5, arom: 5 },
    { name: 'imidazole',     smi: 'c1c[nH]cn1',   atoms: 5, arom: 5 },
    { name: 'oxazole',       smi: 'c1ocnc1',      atoms: 5, arom: 5 },
    { name: 'thiazole',      smi: 'c1scnc1',      atoms: 5, arom: 5 },
    { name: 'quinoline',     smi: 'c1ccc2ncccc2c1', atoms: 10, arom: 10 },
    { name: 'isoquinoline',  smi: 'c1ccc2cnccc2c1', atoms: 10, arom: 10 },
    { name: 'benzofuran',    smi: 'c1ccc2occc2c1', atoms: 9,  arom: 9 },
    { name: 'benzothiophene', smi: 'c1ccc2sccc2c1', atoms: 9, arom: 9 },
    { name: 'purine',        smi: 'c1[nH]c2ncncc2n1', atoms: 9, arom: 9 },
    { name: 'carbazole',     smi: 'c1ccc2c(c1)[nH]c3ccccc23', atoms: 13, arom: 13 }
].forEach(function(c) {
    test(c.name + ': heteroaromatic atom/aromatic counts', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0, c.name + ' has parse errors');
        assert.strictEqual(m.atoms.length, c.atoms);
        assert.strictEqual(aromatic(m), c.arom);
    });
});

// ---------------------------------------------------------------------------
// Charged species — zwitterions, multi-charge ions, aromatic cations
// ---------------------------------------------------------------------------
[
    { name: 'glycine_zwitterion', smi: '[NH3+]CC(=O)[O-]',          atoms: 5,  netCharge: 0 },
    { name: 'betaine',            smi: 'C[N+](C)(C)CC(=O)[O-]',     atoms: 8,  netCharge: 0 },
    { name: 'alanine_zwitter',    smi: 'C[C@@H]([NH3+])C(=O)[O-]',  atoms: 6,  netCharge: 0 },
    { name: 'azide',              smi: '[N-]=[N+]=[N-]',            atoms: 3,  netCharge: -1 },
    { name: 'nitrate',            smi: '[O-][N+](=O)[O-]',          atoms: 4,  netCharge: -1 },
    { name: 'phosphate_3minus',   smi: '[O-]P(=O)([O-])[O-]',       atoms: 5,  netCharge: -3 },
    { name: 'sulfate',            smi: '[O-]S(=O)(=O)[O-]',         atoms: 5,  netCharge: -2 },
    { name: 'guanidinium',        smi: 'NC(=[NH2+])N',              atoms: 4,  netCharge: 1 },
    { name: 'quaternary_amm',     smi: 'C[N+](C)(C)CCO',            atoms: 7,  netCharge: 1 },
    { name: 'pyridinium_Nmethyl', smi: 'C[n+]1ccccc1',              atoms: 7,  netCharge: 1 },
    { name: 'phenolate',          smi: '[O-]c1ccccc1',              atoms: 7,  netCharge: -1 },
    { name: 'sodium_benzoate',    smi: '[Na+].[O-]C(=O)c1ccccc1',   atoms: 10, netCharge: 0 },
    { name: 'ammonium_acetate',   smi: '[NH4+].CC(=O)[O-]',         atoms: 5,  netCharge: 0 }
].forEach(function(c) {
    test(c.name + ': atom count and net formal charge preserved', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0, c.name + ' has parse errors');
        assert.strictEqual(m.atoms.length, c.atoms);
        assert.strictEqual(totalCharge(m), c.netCharge);
    });
});

// ---------------------------------------------------------------------------
// Bridged / spiro / cage skeletons — ring count and bond count regression
// ---------------------------------------------------------------------------
[
    { name: 'norbornane',   smi: 'C1CC2CCC1C2',          atoms: 7,  bonds: 8 },
    { name: 'bicyclo_221',  smi: 'C1CC2CC1C2',           atoms: 6,  bonds: 7 },
    { name: 'spiro_5_5',    smi: 'C1CCC2(CC1)CCCC2',     atoms: 10, bonds: 11 },
    { name: 'spiro_3_3',    smi: 'C1CC12CC2',            atoms: 5,  bonds: 6 },
    { name: 'adamantane',   smi: 'C1C2CC3CC1CC(C2)C3',   atoms: 10, bonds: 12 },
    { name: 'cubane_skel',  smi: 'C12C3C4C1C5C2C3C45',   atoms: 8,  bonds: 12 }
].forEach(function(c) {
    test(c.name + ': bridged/spiro/cage atom and bond counts', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0, c.name + ' has parse errors');
        assert.strictEqual(m.atoms.length, c.atoms);
        assert.strictEqual(m.bonds.length, c.bonds);
    });
});

// ---------------------------------------------------------------------------
// Macrocycles — 12-membered lactone, 14-membered tetraaza ring
// ---------------------------------------------------------------------------
test('macrolactone (12-membered ring): 15 atoms, 15 bonds', function() {
    var m = P('CC1OC(=O)CCCCCCC(O)CO1');
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 15);
    assert.strictEqual(m.bonds.length, 15);
});

test('cyclam (1,4,8,11-tetraazacyclotetradecane): 12 atoms, 12 bonds', function() {
    var m = P('C1CNCCNCCNCCN1');
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 12);
    assert.strictEqual(m.bonds.length, 12);
});

// ---------------------------------------------------------------------------
// Sugars / nucleotide fragments — sterically demanding stereo + phosphate
// ---------------------------------------------------------------------------
test('glucopyranose (sugar ring): 12 atoms, 12 bonds, no aromatic', function() {
    var m = P('OC1OC(CO)C(O)C(O)C1O');
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 12);
    assert.strictEqual(m.bonds.length, 12);
    assert.strictEqual(aromatic(m), 0);
});

test('disaccharide (glycoside dimer): 23 atoms, 24 bonds', function() {
    var m = P('OC1OC(CO)C(O)C(O)C1OC2OC(CO)C(O)C(O)C2O');
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 23);
    assert.strictEqual(m.bonds.length, 24);
});

test('inositol (cyclohexane hexaol, 6 stereocentres): 12 atoms, 12 bonds', function() {
    var m = P('O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O');
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 12);
    var stereoC = m.atoms.filter(function(a) { return a.symbol === 'C' && a.chirality; });
    assert.strictEqual(stereoC.length, 6, 'expected 6 chiral C atoms');
});

test('ATP-like nucleotide fragment: 27 atoms, 9 aromatic, net charge 0', function() {
    var smi = 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O';
    var m = P(smi);
    assert.strictEqual(m.parseErrors.length, 0);
    assert.strictEqual(m.atoms.length, 27);
    assert.strictEqual(aromatic(m), 9);  // adenine purine = 9 aromatic atoms
    assert.strictEqual(totalCharge(m), 0);
});

// ---------------------------------------------------------------------------
// Pyrimidone / pyridone-fused heterocyclic patterns (uracil, thymine,
// cytosine) — BIME perceives these aromatic on parse via lowercase form
// ---------------------------------------------------------------------------
[
    { name: 'uracil',   smi: 'O=c1cc[nH]c(=O)[nH]1',     atoms: 8, arom: 6 },
    { name: 'thymine',  smi: 'Cc1c[nH]c(=O)[nH]c1=O',    atoms: 9, arom: 6 },
    { name: 'cytosine', smi: 'Nc1ccnc(=O)[nH]1',         atoms: 8, arom: 6 }
].forEach(function(c) {
    test(c.name + ' (pyrimidone): atom count + 6 aromatic ring atoms', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0);
        assert.strictEqual(m.atoms.length, c.atoms);
        assert.strictEqual(aromatic(m), c.arom);
    });
});

// ---------------------------------------------------------------------------
// SMARTS motifs — biochemically relevant functional-group queries
// ---------------------------------------------------------------------------
function nMatches(targetSmi, smartsStr) {
    return SmartsMatch.match(P(targetSmi), SmartsParser.parse(smartsStr)).length;
}

test('amide motif [NX3][CX3](=O) on tripeptide → 2 matches', function() {
    assert.strictEqual(
        nMatches('NCC(=O)N[C@@H](CO)C(=O)NCC(=O)O', '[NX3][CX3](=O)'),
        2
    );
});

test('sulfonamide motif [NX3]S(=O)(=O) on sulfanilamide → 2 matches', function() {
    assert.strictEqual(
        nMatches('NS(=O)(=O)c1ccc(N)cc1', '[NX3]S(=O)(=O)'),
        2
    );
});

test('thiol [SX2H1] on cysteine fragment → 1 match', function() {
    assert.strictEqual(
        nMatches('NCC(=O)N[C@@H](CS)C(=O)NCC(=O)O', '[SX2H1]'),
        1
    );
});

test('carboxylate [CX3](=O)[O-] on glycine zwitterion → 1 match', function() {
    assert.strictEqual(
        nMatches('[NH3+]CC([O-])=O', '[CX3](=O)[O-]'),
        1
    );
});

test('quaternary ammonium [N+] on choline → 1 match', function() {
    assert.strictEqual(
        nMatches('C[N+](C)(C)CCO', '[N+]'),
        1
    );
});

test('aromatic [nH] on tetrazole → 1 match', function() {
    assert.strictEqual(
        nMatches('c1nnn[nH]1', '[nH]'),
        1
    );
});

test('aromatic [n+] on N-methylpyridinium → 1 match', function() {
    assert.strictEqual(
        nMatches('C[n+]1ccccc1', '[n+]'),
        1
    );
});

// ---------------------------------------------------------------------------
// Substructure / MCS — extension scenarios (sub-structure of a larger
// homologue, common scaffolds across drug variants)
// ---------------------------------------------------------------------------
function G(smi) { return new SMSDGraph.SMSDGraph(P(smi)); }
function opts() { return new SMSDGraph.ChemOptions(); }

test('peptide dipeptide substructure of tripeptide', function() {
    var q = G('NCC(=O)NCC(=O)O');
    var t = G('NCC(=O)N[C@@H](CO)C(=O)NCC(=O)O');
    assert.strictEqual(SMSDVF2.isSubstructure(q, t, opts()), true);
});

test('phosphate substructure of diphosphate (chain extension)', function() {
    var q = G('COP(=O)(O)O');
    var t = G('COP(=O)(O)OP(=O)(O)O');
    assert.strictEqual(SMSDVF2.isSubstructure(q, t, opts()), true);
});

test('benzene substructure of disconnected (benzene . acetic acid)', function() {
    var q = G('c1ccccc1');
    var t = G('c1ccccc1.CC(=O)O');
    assert.strictEqual(SMSDVF2.isSubstructure(q, t, opts()), true);
});

test('MCS(glucose, glucoside dimer) covers monomer (12 atoms)', function() {
    var sz = SMSDMCS.findMCSFromSmiles(
        'OC1OC(CO)C(O)C(O)C1O',
        'OC1OC(CO)C(O)C(O)C1OC2OC(CO)C(O)C(O)C2O'
    ).size;
    assert.strictEqual(sz, 12);
});

test('MCS(benzene, naphthalene) covers benzene ring (6 atoms)', function() {
    var sz = SMSDMCS.findMCSFromSmiles('c1ccccc1', 'c1ccc2ccccc2c1').size;
    assert.strictEqual(sz, 6);
});

test('MCS(cyclohexane, decalin) covers cyclohexane ring (6 atoms)', function() {
    var sz = SMSDMCS.findMCSFromSmiles('C1CCCCC1', 'C1CCC2CCCCC2C1').size;
    assert.strictEqual(sz, 6);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
