/**
 * tests/test_charge_iso_round_trip.js — Charge / isotope / multi-component
 * round-trip regression coverage for BIME v1.0.2.
 *
 * For every case we assert:
 *   - parse → no errors
 *   - atom count and net formal charge match expectation on parse
 *   - SmilesWriter.write → output re-parses without errors
 *   - re-parsed atom count matches original
 *   - re-parsed net formal charge matches original (charge preserved
 *     through canonical normalisation)
 *   - re-parsed isotope count matches original (isotope labels survive)
 *
 * Cases cover:
 *   - simple cations / anions and high formal charges (+3 perchlorate)
 *   - bare metal ions (Fe³⁺, Cu²⁺, Zn²⁺, Mg²⁺)
 *   - multi-charge anions (oxide, sulfide)
 *   - aromatic charged rings (imidazolium, saccharinate)
 *   - radical-like zwitterion (nitroxide N-oxide)
 *   - explicit-H exotic species (hydronium, methide)
 *   - isotope labels (²H, ¹³C, mixed)
 *   - multi-component salts (MgCl₂, copper acetate)
 *
 * No file I/O, no network, runs in <50 ms.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Charge/Isotope Round-Trip');
var test = runner.test;

console.log('Charge/Isotope Round-Trip');

function P(smi) { return SmilesParser.parse(smi); }
function W(mol) { return SmilesWriter.write(mol); }
function totalCharge(m) { return m.atoms.reduce(function(s, a) { return s + a.charge; }, 0); }
function isotopeCount(m) { return m.atoms.filter(function(a) { return a.isotope > 0; }).length; }

// Each case: name, SMILES, expected atom count, expected net formal charge,
// expected number of atoms carrying an isotope label.
var CASES = [
    // Simple cations / anions
    { name: 'fluoride',         smi: '[F-]',                                 atoms: 1, charge: -1, iso: 0 },
    { name: 'hydronium',        smi: '[OH3+]',                               atoms: 1, charge: 1,  iso: 0 },
    { name: 'cyanide',          smi: '[C-]#N',                               atoms: 2, charge: -1, iso: 0 },
    { name: 'nitronium',        smi: '[N+](=O)=O',                           atoms: 3, charge: 1,  iso: 0 },

    // Bare transition / main-group metal ions
    { name: 'fe3_cation',       smi: '[Fe+3]',                               atoms: 1, charge: 3,  iso: 0 },
    { name: 'cu2_cation',       smi: '[Cu+2]',                               atoms: 1, charge: 2,  iso: 0 },
    { name: 'zn2_cation',       smi: '[Zn+2]',                               atoms: 1, charge: 2,  iso: 0 },

    // Multi-charge anions
    { name: 'oxide_dianion',    smi: '[O-2]',                                atoms: 1, charge: -2, iso: 0 },
    { name: 'sulfide_dianion',  smi: '[S-2]',                                atoms: 1, charge: -2, iso: 0 },

    // High formal charge with negative ligands (perchlorate Cl(VII))
    { name: 'perchlorate',      smi: '[Cl+3]([O-])([O-])([O-])[O-]',         atoms: 5, charge: -1, iso: 0 },

    // Aromatic charged rings
    { name: 'imidazolium',      smi: 'c1[nH+]c[nH]c1',                       atoms: 5, charge: 1,  iso: 0 },
    { name: 'saccharinate',     smi: 'O=C1[N-]S(=O)(=O)c2ccccc12',           atoms: 12, charge: -1, iso: 0 },
    { name: 'phenolate',        smi: '[O-]c1ccccc1',                         atoms: 7, charge: -1, iso: 0 },
    { name: 'tropylium',        smi: '[cH+]1cccccc1',                        atoms: 7, charge: 1,  iso: 0 },

    // Zwitterionic / N-oxide species
    { name: 'trimethylamine_Noxide', smi: 'C[N+]([O-])(C)C',                 atoms: 5, charge: 0,  iso: 0 },
    { name: 'azide',            smi: '[N-]=[N+]=[N-]',                       atoms: 3, charge: -1, iso: 0 },

    // Isotope labels
    { name: 'deuterated_methane', smi: '[2H]C([2H])([2H])[2H]',              atoms: 5, charge: 0,  iso: 4 },
    { name: 'C13_methane',      smi: '[13CH4]',                              atoms: 1, charge: 0,  iso: 1 },
    { name: 'C13_with_D',       smi: '[2H][13CH3]',                          atoms: 2, charge: 0,  iso: 2 },
    { name: 'O18_water',        smi: '[18OH2]',                              atoms: 1, charge: 0,  iso: 1 },

    // Multi-component salts
    { name: 'magnesium_chloride', smi: '[Mg+2].[Cl-].[Cl-]',                 atoms: 3, charge: 0,  iso: 0 },
    { name: 'copper_acetate',   smi: '[Cu+2].CC(=O)[O-].CC(=O)[O-]',         atoms: 9, charge: 0,  iso: 0 },
    { name: 'sodium_phenolate', smi: '[Na+].[O-]c1ccccc1',                   atoms: 8, charge: 0,  iso: 0 }
];

CASES.forEach(function(c) {
    test(c.name + ': parses with expected count + charge + isotopes', function() {
        var m = P(c.smi);
        assert.strictEqual(m.parseErrors.length, 0,
                           'parse errors: ' + m.parseErrors.join('; '));
        assert.strictEqual(m.atoms.length, c.atoms, 'atom count');
        assert.strictEqual(totalCharge(m), c.charge, 'net formal charge');
        assert.strictEqual(isotopeCount(m), c.iso, 'isotope count');
    });

    test(c.name + ': round-trip preserves count + charge + isotopes', function() {
        var m1 = P(c.smi);
        var written = W(m1);
        assert.ok(written && written.length > 0, 'writer produced empty string');
        var m2 = P(written);
        assert.strictEqual(m2.parseErrors.length, 0,
                           're-parse errors on ' + written + ': ' +
                           m2.parseErrors.join('; '));
        assert.strictEqual(m2.atoms.length, c.atoms, 'atom count drift on round-trip');
        assert.strictEqual(totalCharge(m2), c.charge, 'charge drift on round-trip');
        assert.strictEqual(isotopeCount(m2), c.iso, 'isotope drift on round-trip');
    });
});

// ---------------------------------------------------------------------------
// Multi-component salt — fragment count preserved through writer
// ---------------------------------------------------------------------------
test('MgCl2 written form preserves three fragments separated by .', function() {
    var m = P('[Mg+2].[Cl-].[Cl-]');
    var written = W(m);
    var dotCount = (written.match(/\./g) || []).length;
    assert.strictEqual(dotCount, 2, 'expected 2 dots in ' + written);
});

test('copper acetate (3 fragments) preserves three fragments through writer', function() {
    var m = P('[Cu+2].CC(=O)[O-].CC(=O)[O-]');
    var written = W(m);
    var dotCount = (written.match(/\./g) || []).length;
    assert.strictEqual(dotCount, 2, 'expected 2 dots in ' + written);
});

// ---------------------------------------------------------------------------
// Per-atom charge preservation on the multi-charge ion (perchlorate
// must keep one Cl(+3) and four O(-1) atoms — atoms whose individual
// charge is non-zero must round-trip with the same per-atom value).
// ---------------------------------------------------------------------------
test('perchlorate per-atom charges: Cl(+3) and 4×O(-1) preserved', function() {
    var m1 = P('[Cl+3]([O-])([O-])([O-])[O-]');
    var w  = W(m1);
    var m2 = P(w);
    var sortedCharges = m2.atoms.map(function(a) { return a.charge; }).sort();
    assert.deepStrictEqual(sortedCharges, [-1, -1, -1, -1, 3]);
});

// ---------------------------------------------------------------------------
// Specific isotope value preserved exactly (not coerced)
// ---------------------------------------------------------------------------
test('isotope value 13 preserved exactly through round-trip', function() {
    var m1 = P('[13CH4]');
    var w  = W(m1);
    var m2 = P(w);
    assert.strictEqual(m2.atoms[0].isotope, 13);
});

test('isotope value 18 preserved exactly through round-trip', function() {
    var m1 = P('[18OH2]');
    var w  = W(m1);
    var m2 = P(w);
    assert.strictEqual(m2.atoms[0].isotope, 18);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
