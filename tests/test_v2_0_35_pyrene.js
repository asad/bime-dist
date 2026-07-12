/**
 * tests/test_v2_0_35_pyrene.js — geometrically-correct peri-condensed
 * pyrene template (v2.0.35).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Pyrene C16H10 — 4 hexagons sharing a central pair (10b, 10c) such that
 * each central atom is in 3 rings. Tests pin atom + bond counts, exact
 * bond lengths, degree distribution, the 4-cycle SSSR, the LOOKUP route
 * registration, and end-to-end parse + layout.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/Layout.js');
require('../editor/Templates.js');

var runner = shim.makeRunner('Pyrene template (v2.0.35)');
var test = runner.test;

console.log('Pyrene template (v2.0.35)');

var BL_EXPECTED = 30;
var BL_TOL = 0.05;

test('A1 16 atoms, 19 bonds, fully aromatic', function() {
    var t = Templates.pyrene();
    assert.strictEqual(t.atoms.length, 16);
    assert.strictEqual(t.bonds.length, 19);
    assert.strictEqual(t.aromatic.length, 16);
});

test('A2 every atom is carbon', function() {
    var t = Templates.pyrene();
    for (var i = 0; i < t.atoms.length; i++) {
        assert.strictEqual(t.atoms[i].symbol, 'C');
    }
});

test('B1 every bond length is exactly BL ± 0.05 px', function() {
    var t = Templates.pyrene();
    for (var i = 0; i < t.bonds.length; i++) {
        var b = t.bonds[i];
        var a = t.atoms[b.a1], c = t.atoms[b.a2];
        var d = Math.sqrt(Math.pow(a.x - c.x, 2) + Math.pow(a.y - c.y, 2));
        assert.ok(Math.abs(d - BL_EXPECTED) < BL_TOL,
            'bond ' + i + ' (atoms ' + b.a1 + '-' + b.a2 + ') length ' + d.toFixed(3));
    }
});

function adjacency(t) {
    var adj = [];
    for (var i = 0; i < t.atoms.length; i++) { adj.push([]); }
    for (var bi = 0; bi < t.bonds.length; bi++) {
        var b = t.bonds[bi];
        adj[b.a1].push(b.a2);
        adj[b.a2].push(b.a1);
    }
    return adj;
}

test('C1 degree distribution: 10 peripheral (deg 2, 1 H each) + 6 internal (deg 3, no H)', function() {
    var t = Templates.pyrene();
    var adj = adjacency(t);
    var counts = { 2: 0, 3: 0 };
    for (var i = 0; i < adj.length; i++) {
        var d = adj[i].length;
        if (counts[d] === undefined) { assert.fail('atom ' + i + ' degree ' + d); }
        counts[d]++;
    }
    assert.strictEqual(counts[2], 10);
    assert.strictEqual(counts[3], 6);
});

function find6Cycles(adj) {
    var n = adj.length;
    var cycles = [];
    var seen = {};
    function key(arr) { return arr.slice().sort(function (a,b) { return a - b; }).join(','); }
    for (var s = 0; s < n; s++) {
        var stack = [{ at: s, path: [s], depth: 0 }];
        while (stack.length) {
            var rec = stack.pop();
            if (rec.depth === 6) {
                if (rec.at === s && rec.path.length === 6) {
                    var k = key(rec.path);
                    if (!seen[k]) { seen[k] = true; cycles.push(rec.path); }
                }
                continue;
            }
            var nbs = adj[rec.at];
            for (var ni = 0; ni < nbs.length; ni++) {
                var nb = nbs[ni];
                if (rec.depth < 5 && rec.path.indexOf(nb) >= 0) continue;
                if (rec.depth === 5 && nb !== s) continue;
                var nextPath = (rec.depth < 5) ? rec.path.concat([nb]) : rec.path;
                stack.push({ at: nb, path: nextPath, depth: rec.depth + 1 });
            }
        }
    }
    return cycles;
}

test('C2 exactly 4 distinct 6-membered cycles (pyrene SSSR)', function() {
    var t = Templates.pyrene();
    var cycles = find6Cycles(adjacency(t));
    assert.strictEqual(cycles.length, 4, 'found ' + cycles.length + ' 6-cycles, expected 4');
});

test('C3 atoms 0 (10b) and 1 (10c) are each in 3 of the 4 rings', function() {
    var t = Templates.pyrene();
    var cycles = find6Cycles(adjacency(t));
    var inRings0 = 0, inRings1 = 0;
    for (var i = 0; i < cycles.length; i++) {
        if (cycles[i].indexOf(0) >= 0) inRings0++;
        if (cycles[i].indexOf(1) >= 0) inRings1++;
    }
    assert.strictEqual(inRings0, 3, 'atom 0 (10b) in 3 rings');
    assert.strictEqual(inRings1, 3, 'atom 1 (10c) in 3 rings');
});

test('D1 TEMPLATE_LOOKUP routes "4:6,6,6,6" through pyrene', function() {
    var fs = require('fs');
    var path = require('path');
    var src = fs.readFileSync(path.join(__dirname, '..', 'editor', 'Layout.js'), 'utf8');
    var m = src.match(/'4:6,6,6,6'\s*:\s*\[([^\]]*)\]/);
    assert.ok(m, 'TEMPLATE_LOOKUP has the 4:6,6,6,6 key');
    assert.ok(m[1].indexOf("'pyrene'") >= 0, 'pyrene is among the candidates');
});

test('E1 BIME parses pyrene SMILES + Layout.layout produces finite coords', function() {
    var mol = new Molecule();
    SmilesParser.parse('c1cc2ccc3cccc4ccc(c1)c2c34', mol);
    assert.strictEqual(mol.atoms.length, 16, 'parser produces 16 atoms');
    Layout.layout(mol);
    assert.strictEqual(mol.atoms.length, 16, 'layout preserves atom count');
    for (var i = 0; i < mol.atoms.length; i++) {
        assert.ok(isFinite(mol.atoms[i].x), 'atom ' + i + ' x finite');
        assert.ok(isFinite(mol.atoms[i].y), 'atom ' + i + ' y finite');
    }
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
