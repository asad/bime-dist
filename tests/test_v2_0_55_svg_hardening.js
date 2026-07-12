/**
 * tests/test_v2_0_55_svg_hardening.js — defensive SVG attribute hardening
 * in `AtomTrace.renderInteractiveSvg` (v2.0.55).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Background. The v2.0.54 release-review panel flagged a pre-existing
 * security non-blocker: `AtomTrace.renderInteractiveSvg` interpolates
 * `atom.symbol` into SVG attribute and text content via string
 * concatenation. The shipped path is safe (atom.symbol is canonicalised
 * by `SmilesParser`), but if a future caller hand-builds a Molecule the
 * raw string could end up inside `<text>` and `data-element`. v2.0.55
 * adds a regex assert that rejects anything that isn't a 1-or-2-character
 * Latin element symbol, substituting '?' in its place.
 *
 * Plain Node, no DOM.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/AtomTrace.js');

var runner = shim.makeRunner('SVG hardening (v2.0.55)');
var test = runner.test;

console.log('SVG hardening (v2.0.55)');

// ---------------------------------------------------------------------------
// A. legitimate inputs still render unchanged
// ---------------------------------------------------------------------------
test('A1 a normal SMILES-parsed molecule still renders one <g> per heavy atom', function () {
    var m = new Molecule();
    SmilesParser.parse('CCO', m);
    var svg = AtomTrace.renderInteractiveSvg(m, { width: 100, height: 60 });
    var groups = svg.match(/<g class="bime-trace-atom"/g) || [];
    assert.strictEqual(groups.length, 3, 'three heavy atoms → three <g> groups');
    // Element labels render as plain element symbols.
    assert.ok(svg.indexOf('>O</text>') !== -1, 'oxygen element label rendered');
});

test('A2 two-letter halogen symbols (Cl, Br, etc.) still render correctly', function () {
    var m = new Molecule();
    SmilesParser.parse('CCCl', m);
    var svg = AtomTrace.renderInteractiveSvg(m, { width: 100, height: 60 });
    assert.ok(svg.indexOf('>Cl</text>') !== -1, 'Cl label rendered');
    assert.ok(svg.indexOf('data-element="Cl"') !== -1, 'data-element attribute carries Cl');
});

// ---------------------------------------------------------------------------
// B. tampered / hand-built atoms are sanitised to '?' (no markup smuggling)
// ---------------------------------------------------------------------------
test('B1 atom.symbol containing angle brackets is replaced with "?"', function () {
    // Hand-built Molecule with a tampered atom.symbol.
    var m = new Molecule();
    m.atoms.push({ id: 1, symbol: '<script>x</script>', x: 0, y: 0, charge: 0, isotope: 0 });
    m.bonds = [];
    var svg = AtomTrace.renderInteractiveSvg(m, { width: 100, height: 60 });
    assert.strictEqual(svg.indexOf('<script>'), -1, 'no script tag escapes into the SVG');
    assert.ok(svg.indexOf('data-element="?"') !== -1,
        'data-element substituted to ? for malformed symbol');
    assert.ok(svg.indexOf('>?</text>') !== -1, 'text content substituted to ?');
});

test('B2 atom.symbol that is not a string at all (number, object, null) is replaced with "?"', function () {
    var m = new Molecule();
    m.atoms.push({ id: 1, symbol: 42, x: 0, y: 0, charge: 0, isotope: 0 });
    m.atoms.push({ id: 2, symbol: null, x: 1, y: 0, charge: 0, isotope: 0 });
    m.atoms.push({ id: 3, symbol: { toString: function () { return '<x>'; } },
                   x: 2, y: 0, charge: 0, isotope: 0 });
    m.bonds = [];
    var svg = AtomTrace.renderInteractiveSvg(m, { width: 100, height: 60 });
    // No injected markup from any of the three tampered symbols.
    assert.strictEqual(svg.indexOf('<x>'), -1, 'object.toString() result does not escape');
    // Three '?' substitutions.
    var qm = (svg.match(/data-element="\?"/g) || []).length;
    assert.strictEqual(qm, 3, 'three tampered atoms → three "?" substitutions');
});

test('B3 atom.symbol with 3+ characters or non-Latin chars is replaced with "?"', function () {
    var m = new Molecule();
    m.atoms.push({ id: 1, symbol: 'Carb', x: 0, y: 0, charge: 0, isotope: 0 });
    m.atoms.push({ id: 2, symbol: 'C1', x: 1, y: 0, charge: 0, isotope: 0 });
    m.atoms.push({ id: 3, symbol: 'C"X', x: 2, y: 0, charge: 0, isotope: 0 });
    m.bonds = [];
    var svg = AtomTrace.renderInteractiveSvg(m, { width: 100, height: 60 });
    // None of those should land in the markup.
    assert.strictEqual(svg.indexOf('"Carb"'), -1);
    assert.strictEqual(svg.indexOf('"C1"'), -1);
    assert.strictEqual(svg.indexOf('"C"X"'), -1);
    var qm = (svg.match(/data-element="\?"/g) || []).length;
    assert.strictEqual(qm, 3, 'all three rejected; three "?" substitutions');
});

// ---------------------------------------------------------------------------
// C. AtomTrace.version stamp
// ---------------------------------------------------------------------------
test('C1 AtomTrace.version is at least 2.0.55 (forward-compatible)', function () {
    // The SVG hardening shipped in v2.0.55; later patches keep the assert.
    assert.ok(typeof AtomTrace.version === 'string');
    var parts = AtomTrace.version.split('.').map(Number);
    var atLeast = (parts[0] > 2) || (parts[0] === 2 && (parts[1] > 0 || parts[2] >= 55));
    assert.ok(atLeast, 'AtomTrace.version >= 2.0.55, got ' + AtomTrace.version);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
