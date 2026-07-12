/**
 * tests/test_v2_0_19_metabolite_library.js — built-in metabolite repertoire +
 * pathway templates (v2.0.19).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Validates that every drawn-structure metabolite parses in BIME's own
 * SmilesParser with the expected formula, that name/alias lookup works, and
 * that the glycolysis template (mirroring KEGG map00010) is internally
 * consistent (every substrate/product resolves, every arrow type is valid).
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/MetaboliteLibrary.js');

var runner = shim.makeRunner('Metabolite Library (v2.0.19)');
var test = runner.test;
var LIB = MetaboliteLibrary;
var ARROW_TYPES = ['forward', 'reverse', 'reversible', 'resonance'];

console.log('Metabolite Library (v2.0.19)');

function formulaOf(smiles) {
    var m = new Molecule();
    SmilesParser.parse(smiles, m);
    return { atoms: m.atoms.length, formula: Molecule.formulaString(m.elementCounts()) };
}

// ---------------------------------------------------------------------------
// A. library shape
// ---------------------------------------------------------------------------
test('A1 library exposes metabolites, pathways, find()', function() {
    assert.ok(Array.isArray(LIB.metabolites) && LIB.metabolites.length >= 25, 'metabolites array');
    assert.ok(LIB.pathways && LIB.pathways.glycolysis, 'glycolysis template present');
    assert.strictEqual(typeof LIB.find, 'function');
});

test('A2 every entry is well-formed', function() {
    LIB.metabolites.forEach(function(m) {
        assert.ok(m.name && typeof m.name === 'string', 'name: ' + JSON.stringify(m));
        assert.ok(typeof m.smiles === 'string', 'smiles string: ' + m.name);
        assert.ok(m.display === 'structure' || m.display === 'label', 'display enum: ' + m.name);
        assert.ok(Array.isArray(m.aliases), 'aliases array: ' + m.name);
        assert.ok(/^C\d{5}$/.test(m.kegg), 'KEGG id format: ' + m.name + ' = ' + m.kegg);
    });
});

// ---------------------------------------------------------------------------
// B. every drawn structure parses with atoms
// ---------------------------------------------------------------------------
test('B1 all display:structure SMILES parse with >0 atoms', function() {
    LIB.metabolites.filter(function(m) { return m.display === 'structure'; }).forEach(function(m) {
        assert.ok(m.smiles.length > 0, 'structure has SMILES: ' + m.name);
        var r = formulaOf(m.smiles);
        assert.ok(r.atoms > 0, 'parses with atoms: ' + m.name);
    });
});

test('B2 key formulas are correct (chemistry sanity)', function() {
    assert.strictEqual(formulaOf(LIB.find('Glucose').smiles).formula, 'C6H12O6');
    assert.strictEqual(formulaOf(LIB.find('Pyruvate').smiles).formula, 'C3H4O3');
    assert.strictEqual(formulaOf(LIB.find('PEP').smiles).formula, 'C3H5O6P');
    assert.strictEqual(formulaOf(LIB.find('Lactate').smiles).formula, 'C3H6O3');
    // G6P and F6P are isomers — same formula, distinct entries.
    assert.strictEqual(formulaOf(LIB.find('G6P').smiles).formula, 'C6H13O9P');
    assert.strictEqual(formulaOf(LIB.find('F6P').smiles).formula, 'C6H13O9P');
});

// ---------------------------------------------------------------------------
// C. lookup
// ---------------------------------------------------------------------------
test('C1 find() resolves canonical name, alias, and is case-insensitive', function() {
    assert.strictEqual(LIB.find('D-Glucose').kegg, 'C00031');
    assert.strictEqual(LIB.find('Glc').kegg, 'C00031');       // alias
    assert.strictEqual(LIB.find('glucose').kegg, 'C00031');   // case-insensitive alias
    assert.strictEqual(LIB.find('GAP').name, 'Glyceraldehyde 3-phosphate');
    assert.strictEqual(LIB.find('ATP').display, 'label');
    assert.strictEqual(LIB.find('totally-unknown-xyz'), null);
    assert.strictEqual(LIB.find(''), null);
});

test('C2 big cofactors are labelled nodes, not drawn', function() {
    ['ATP', 'ADP', 'NAD+', 'NADH', 'CoA', 'Acetyl-CoA'].forEach(function(n) {
        var m = LIB.find(n);
        assert.ok(m, 'present: ' + n);
        assert.strictEqual(m.display, 'label', n + ' is a labelled node');
    });
});

// ---------------------------------------------------------------------------
// D. glycolysis template (KEGG map00010) consistency
// ---------------------------------------------------------------------------
test('D1 glycolysis has the 11 canonical steps and maps to KEGG map00010', function() {
    var g = LIB.pathways.glycolysis;
    assert.strictEqual(g.kegg, 'map00010');
    assert.strictEqual(g.steps.length, 11);
});

test('D2 every step substrate/product resolves in the library', function() {
    LIB.pathways.glycolysis.steps.forEach(function(s) {
        assert.ok(LIB.find(s.from), 'from resolves: ' + s.from);
        assert.ok(LIB.find(s.to), 'to resolves: ' + s.to);
    });
});

test('D3 every step has a valid arrow type, EC, and enzyme', function() {
    LIB.pathways.glycolysis.steps.forEach(function(s) {
        assert.ok(ARROW_TYPES.indexOf(s.arrowType) >= 0, 'arrowType valid: ' + s.arrowType);
        assert.ok(/^\d+\.\d+\.\d+\.\d+$/.test(s.ec), 'EC format: ' + s.enzyme + ' = ' + s.ec);
        assert.ok(s.enzyme && s.enzyme.length > 0, 'enzyme named');
    });
});

test('D4 the three irreversible commitment steps use forward arrows', function() {
    var byEnzyme = {};
    LIB.pathways.glycolysis.steps.forEach(function(s) { byEnzyme[s.enzyme] = s; });
    assert.strictEqual(byEnzyme['Hexokinase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Phosphofructokinase-1'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Pyruvate kinase'].arrowType, 'forward');
    // a reversible step uses the equilibrium arrow
    assert.strictEqual(byEnzyme['Enolase'].arrowType, 'reversible');
});

// ---------------------------------------------------------------------------
// E. TCA cycle (KEGG map00020) consistency
// ---------------------------------------------------------------------------
test('E1 TCA cycle template is present, shape cycle, and closes', function() {
    var tca = LIB.pathways.tca;
    assert.ok(tca, 'tca template exists');
    assert.strictEqual(tca.shape, 'cycle');
    assert.strictEqual(tca.kegg, 'map00020');
    assert.strictEqual(tca.steps.length, 9);
    // the last step's product must equal the first step's substrate (closure)
    assert.strictEqual(tca.steps[0].from, tca.steps[tca.steps.length - 1].to,
        'cycle closes: L-Malate -> Oxaloacetate -> Citrate');
});

test('E2 every TCA step resolves + has a valid arrow type and EC', function() {
    LIB.pathways.tca.steps.forEach(function(s) {
        assert.ok(LIB.find(s.from), 'from resolves: ' + s.from);
        assert.ok(LIB.find(s.to), 'to resolves: ' + s.to);
        assert.ok(ARROW_TYPES.indexOf(s.arrowType) >= 0, 'arrowType valid: ' + s.arrowType);
        assert.ok(/^\d+\.\d+\.\d+\.\d+$/.test(s.ec), 'EC: ' + s.enzyme + ' = ' + s.ec);
    });
});

test('E3 key TCA-intermediate formulas are correct (chemistry sanity)', function() {
    assert.strictEqual(formulaOf(LIB.find('Citrate').smiles).formula, 'C6H8O7');
    assert.strictEqual(formulaOf(LIB.find('Isocitrate').smiles).formula, 'C6H8O7'); // isomer of citrate
    assert.strictEqual(formulaOf(LIB.find('Succinate').smiles).formula, 'C4H6O4');
    assert.strictEqual(formulaOf(LIB.find('Fumarate').smiles).formula, 'C4H4O4');
    assert.strictEqual(formulaOf(LIB.find('L-Malate').smiles).formula, 'C4H6O5');
    assert.strictEqual(formulaOf(LIB.find('Oxaloacetate').smiles).formula, 'C4H4O5');
    assert.strictEqual(formulaOf(LIB.find('2-Oxoglutarate').smiles).formula, 'C5H6O5');
});

test('E4 the three irreversible TCA steps use forward arrows', function() {
    var byEnzyme = {};
    LIB.pathways.tca.steps.forEach(function(s) { byEnzyme[s.enzyme] = s; });
    // Citrate synthase, isocitrate DH, alpha-KG DH are physiologically irreversible
    assert.strictEqual(byEnzyme['Citrate synthase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Isocitrate dehydrogenase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['2-Oxoglutarate dehydrogenase complex'].arrowType, 'forward');
});

// ---------------------------------------------------------------------------
// F. Urea cycle (Krebs–Henseleit, KEGG map00220) consistency
// ---------------------------------------------------------------------------
test('F1 Urea cycle template is present, shape cycle, 4-step closure', function() {
    var u = LIB.pathways.urea;
    assert.ok(u, 'urea template exists');
    assert.strictEqual(u.shape, 'cycle');
    assert.strictEqual(u.kegg, 'map00220');
    assert.strictEqual(u.steps.length, 4);
    assert.strictEqual(u.steps[0].from, u.steps[u.steps.length - 1].to,
        'cycle closes: Arginine -> Ornithine -> Citrulline');
});

test('F2 every Urea step resolves + valid arrow type and EC', function() {
    LIB.pathways.urea.steps.forEach(function(s) {
        assert.ok(LIB.find(s.from), 'from resolves: ' + s.from);
        assert.ok(LIB.find(s.to), 'to resolves: ' + s.to);
        assert.ok(ARROW_TYPES.indexOf(s.arrowType) >= 0, 'arrowType: ' + s.arrowType);
        assert.ok(/^\d+\.\d+\.\d+\.\d+$/.test(s.ec), 'EC: ' + s.enzyme + ' = ' + s.ec);
    });
});

test('F3 key urea-cycle formulas are correct (chemistry sanity)', function() {
    assert.strictEqual(formulaOf(LIB.find('Ornithine').smiles).formula, 'C5H12N2O2');
    assert.strictEqual(formulaOf(LIB.find('Citrulline').smiles).formula, 'C6H13N3O3');
    assert.strictEqual(formulaOf(LIB.find('Argininosuccinate').smiles).formula, 'C10H18N4O6');
    assert.strictEqual(formulaOf(LIB.find('Arginine').smiles).formula, 'C6H14N4O2');
    assert.strictEqual(formulaOf(LIB.find('Carbamoyl phosphate').smiles).formula, 'CH4NO5P');
    assert.strictEqual(formulaOf(LIB.find('Urea').smiles).formula, 'CH4N2O');
});

test('F4 all four urea-cycle steps are physiologically irreversible (forward)', function() {
    LIB.pathways.urea.steps.forEach(function(s) {
        assert.strictEqual(s.arrowType, 'forward', s.enzyme + ' should be forward');
    });
});

// ---------------------------------------------------------------------------
// G. Pentose Phosphate Pathway (KEGG map00030, branched)
// ---------------------------------------------------------------------------
test('G1 PPP template is present with shape branched + 4-node trunk', function() {
    var p = LIB.pathways.ppp;
    assert.ok(p, 'ppp template exists');
    assert.strictEqual(p.shape, 'branched');
    assert.strictEqual(p.kegg, 'map00030');
    assert.ok(Array.isArray(p.trunk) && p.trunk.length === 4, 'trunk is 4-node oxidative phase');
    assert.strictEqual(p.steps.length, 5, '3 trunk + 2 branch steps');
});

test('G2 every PPP step resolves + valid arrow type + EC', function() {
    LIB.pathways.ppp.steps.forEach(function(s) {
        assert.ok(LIB.find(s.from), 'from resolves: ' + s.from);
        assert.ok(LIB.find(s.to), 'to resolves: ' + s.to);
        assert.ok(ARROW_TYPES.indexOf(s.arrowType) >= 0, 'arrowType: ' + s.arrowType);
        assert.ok(/^\d+\.\d+\.\d+\.\d+$/.test(s.ec), 'EC: ' + s.enzyme + ' = ' + s.ec);
    });
});

test('G3 key PPP-intermediate formulas are correct (chemistry sanity)', function() {
    assert.strictEqual(formulaOf(LIB.find('6PGL').smiles).formula, 'C6H11O9P');
    assert.strictEqual(formulaOf(LIB.find('6PG').smiles).formula, 'C6H13O10P');
    // Ru5P, R5P, Xu5P are isomers — all C5H11O8P
    assert.strictEqual(formulaOf(LIB.find('Ru5P').smiles).formula, 'C5H11O8P');
    assert.strictEqual(formulaOf(LIB.find('R5P').smiles).formula, 'C5H11O8P');
    assert.strictEqual(formulaOf(LIB.find('Xu5P').smiles).formula, 'C5H11O8P');
});

test('G4 the three oxidative-phase steps are forward; non-oxidative isomerase/epimerase reversible', function() {
    var byEnzyme = {};
    LIB.pathways.ppp.steps.forEach(function(s) { byEnzyme[s.enzyme] = s; });
    assert.strictEqual(byEnzyme['Glucose-6-phosphate dehydrogenase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Gluconolactonase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['6-Phosphogluconate dehydrogenase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Ribose-5-phosphate isomerase'].arrowType, 'reversible');
    assert.strictEqual(byEnzyme['Ribulose-5-phosphate 3-epimerase'].arrowType, 'reversible');
});

// ---------------------------------------------------------------------------
// H. Pyruvate fates branch
// ---------------------------------------------------------------------------
test('H1 pyruvate fates template is present with fanout shape + 4 branch reactions', function() {
    var p = LIB.pathways.pyruvate_fates;
    assert.ok(p, 'pyruvate_fates template exists');
    assert.strictEqual(p.shape, 'fanout');
    assert.strictEqual(p.hub, 'Pyruvate');
    assert.strictEqual(p.kegg, 'map00620');
    assert.strictEqual(p.steps.length, 4);
});

test('H2 every pyruvate-fates branch resolves + valid arrow type + EC', function() {
    LIB.pathways.pyruvate_fates.steps.forEach(function(s) {
        assert.ok(LIB.find(s.from), 'from resolves: ' + s.from);
        assert.ok(LIB.find(s.to), 'to resolves: ' + s.to);
        assert.strictEqual(s.from, 'Pyruvate', 'fanout branches start at pyruvate');
        assert.ok(ARROW_TYPES.indexOf(s.arrowType) >= 0, 'arrowType: ' + s.arrowType);
        assert.ok(/^\d+\.\d+\.\d+\.\d+$/.test(s.ec), 'EC: ' + s.enzyme + ' = ' + s.ec);
    });
});

test('H3 pyruvate-fates arrow policy keeps fermentative reduction reversible and oxidative/carboxylating fates forward', function() {
    var byEnzyme = {};
    LIB.pathways.pyruvate_fates.steps.forEach(function(s) { byEnzyme[s.enzyme] = s; });
    assert.strictEqual(byEnzyme['Lactate dehydrogenase'].arrowType, 'reversible');
    assert.strictEqual(byEnzyme['Pyruvate decarboxylase'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Pyruvate dehydrogenase system'].arrowType, 'forward');
    assert.strictEqual(byEnzyme['Pyruvate carboxylase'].arrowType, 'forward');
});

// ---------------------------------------------------------------------------
// I. β-Oxidation spiral (KEGG map00071, loop-iterative)
// ---------------------------------------------------------------------------
test('I1 β-oxidation template is present with shape loop-iterative + 4 reactions', function() {
    var p = LIB.pathways.beta_oxidation;
    assert.ok(p, 'beta_oxidation template exists');
    assert.strictEqual(p.shape, 'loop-iterative');
    assert.strictEqual(p.kegg, 'map00071');
    assert.strictEqual(p.steps.length, 4, '4 steps form one full β-ox turn');
    // closure: step 4's `to` returns to step 1's `from` (the iteration arrow)
    assert.strictEqual(p.steps[3].to, p.steps[0].from, 'spiral closes butyryl → butyryl conceptually');
});

test('I2 every β-oxidation step resolves + valid arrow type + EC', function() {
    LIB.pathways.beta_oxidation.steps.forEach(function(s) {
        assert.ok(LIB.find(s.from), 'from resolves: ' + s.from);
        assert.ok(LIB.find(s.to), 'to resolves: ' + s.to);
        assert.ok(ARROW_TYPES.indexOf(s.arrowType) >= 0, 'arrowType: ' + s.arrowType);
        assert.ok(/^\d+\.\d+\.\d+\.\d+$/.test(s.ec), 'EC: ' + s.enzyme + ' = ' + s.ec);
    });
});

test('I3 key β-oxidation-intermediate formulas are correct (S-methyl thioester abbreviation, chemistry sanity)', function() {
    assert.strictEqual(formulaOf(LIB.find('Butyryl-CoA').smiles).formula, 'C5H10OS');
    assert.strictEqual(formulaOf(LIB.find('trans-Crotonyl-CoA').smiles).formula, 'C5H8OS');
    assert.strictEqual(formulaOf(LIB.find('L-3-Hydroxybutyryl-CoA').smiles).formula, 'C5H10O2S');
    assert.strictEqual(formulaOf(LIB.find('Acetoacetyl-CoA').smiles).formula, 'C5H8O2S');
});

test('I4 β-oxidation arrow policy: oxidation + thiolysis forward, hydration + dehydrogenation reversible', function() {
    var byEnzyme = {};
    LIB.pathways.beta_oxidation.steps.forEach(function(s) { byEnzyme[s.enzyme] = s; });
    assert.strictEqual(byEnzyme['Acyl-CoA dehydrogenase'].arrowType, 'forward', 'FAD oxidation drives irreversibility');
    assert.strictEqual(byEnzyme['Enoyl-CoA hydratase'].arrowType, 'reversible');
    assert.strictEqual(byEnzyme['3-Hydroxyacyl-CoA dehydrogenase'].arrowType, 'reversible');
    assert.strictEqual(byEnzyme['β-Ketothiolase'].arrowType, 'forward', 'thiolytic cleavage with CoA-SH is irreversible');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
