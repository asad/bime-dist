/**
 * tests/test_v3_0_0_library_valence_sanity.js — every molecule in the shipped
 * common-molecules library must parse to a CHEMICALLY POSSIBLE structure.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Why this exists
 * ---------------
 * An early "canonicalise every SMILES to aromatic form" pass (commit 984e5a2)
 * over-aromatised NON-aromatic rings — it wrote saturated / quaternary carbons
 * as lowercase aromatic `c`, producing chemically impossible degree-4 aromatic
 * carbons. 38 teaching-relevant molecules (corticosteroids, sex steroids,
 * retinoids, Gabapentin, Sucralose, Paclitaxel, Digoxin, several opioids and
 * antipsychotics, …) parsed to the WRONG molecular formula as a result — a
 * "Never mislabel" liability for a tool used to teach students. v3.0.0 replaces
 * each with a structure verified against its authoritative record (formula +
 * InChIKey connectivity).
 *
 * This test makes the whole class un-regressable: it parses EVERY library entry
 * and fails on any chemically impossible atom — a degree-4 aromatic carbon or a
 * negative implicit-hydrogen count — or any parse error. It is deliberately
 * exhaustive (library-wide), not a per-molecule golden set, so a future bad
 * SMILES anywhere in the library is caught immediately.
 *
 * Real engine, no DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));
shim.loadAll();

var SmilesParser = globalThis.SmilesParser;
var Molecule = globalThis.Molecule;

/* eslint-disable no-eval */
// strict-mode eval scopes the file's `var COMMON_MOLECULES` internally, so
// capture it via a trailing expression (same trick as the v2.4.17 golden test).
var COMMON_MOLECULES = eval(
    fs.readFileSync(path.join(__dirname, '..', 'common-molecules.js'), 'utf8') +
    '\n;COMMON_MOLECULES;');

var runner = shim.makeRunner('Library valence sanity (v3.0.0)');
var test = runner.test;
console.log('Library valence sanity (v3.0.0)');

// A degree-4 aromatic carbon is impossible (an aromatic sp2 carbon has at most 3
// sigma neighbours). A negative implicit-H count means the written bonds already
// exceed the element's valence. Either signals a malformed SMILES.
function impossibleAtoms(m) {
    var probs = [];
    for (var i = 0; i < m.atoms.length; i++) {
        var a = m.atoms[i];
        var el = a.element || a.symbol;
        var deg = m.getBondsOfAtom(a.id).length;
        var h;
        try { h = m.calcHydrogens(a.id); } catch (e) { h = null; }
        if (a.aromatic && el === 'C' && deg >= 4) {
            probs.push('degree-' + deg + ' aromatic C @atom' + a.id);
        }
        if (typeof h === 'number' && h < 0) {
            probs.push('negative implicit H (' + el + '=' + h + ') @atom' + a.id);
        }
    }
    return probs;
}

test('every common-molecules entry parses to a chemically possible structure', function () {
    assert.ok(COMMON_MOLECULES && COMMON_MOLECULES.length > 1000,
        'library loaded (' + (COMMON_MOLECULES && COMMON_MOLECULES.length) + ' entries)');
    var failures = [];
    for (var i = 0; i < COMMON_MOLECULES.length; i++) {
        var e = COMMON_MOLECULES[i];
        var name = Array.isArray(e) ? e[0] : (e.name || e.n);
        var smi = Array.isArray(e) ? e[1] : (e.smiles || e.smi);
        var m;
        try {
            m = SmilesParser.parse(smi);
        } catch (err) {
            failures.push(name + ': PARSE THREW — ' + (err && err.message));
            continue;
        }
        if ((m.parseErrors || []).length) {
            failures.push(name + ': ' + m.parseErrors.length + ' parse error(s) — ' + smi);
            continue;
        }
        var probs = impossibleAtoms(m);
        if (probs.length) {
            failures.push(name + ' [' + smi + '] => ' + probs.join(', '));
        }
    }
    assert.strictEqual(failures.length, 0,
        failures.length + ' library molecule(s) parse to a chemically impossible structure:\n  ' +
        failures.slice(0, 60).join('\n  '));
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
