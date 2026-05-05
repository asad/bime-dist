/**
 * tests/test_v1_5_0_features.js — BIME v1.5.0 unified library search.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Covers MolEditor.prototype.searchLibrary() — single dispatch over the
 * COMMON_MOLECULES library for four search modes:
 *
 *   exact         — canonical SMILES equality
 *   substructure  — VF2++ subgraph isomorphism via SmartsMatch
 *   mcs           — Maximum Common Substructure via SMSDMCS
 *   similarity    — Tanimoto over 1024-bit path fingerprints
 *
 * Tests use a small synthetic 5-entry COMMON_MOLECULES fixture so the
 * suite remains fast and deterministic. The library cache is reset
 * before each test so fixtures don't leak.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/SMSDLayout.js');
require('../editor/SMSDBatch.js');
require('../editor/Templates.js');
require('../editor/Renderer.js');
require('../editor/History.js');
require('../editor/Tools.js');
require('../editor/SmartsWriter.js');
require('../editor/ImageExport.js');
require('../editor/ToolbarPrefs.js');
require('../editor/MolEditor.js');

var runner = shim.makeRunner('Library search v1.5.0');
var test = runner.test;

console.log('Library search v1.5.0');

// -------------------------------------------------------------------------
// Tiny synthetic library so tests are fast and deterministic.
// -------------------------------------------------------------------------
var FIXTURE_LIBRARY = [
    { name: 'Methane',     smiles: 'C',                 category: 'basic' },
    { name: 'Ethanol',     smiles: 'CCO',               category: 'basic' },
    { name: 'Benzene',     smiles: 'c1ccccc1',          category: 'aromatic' },
    { name: 'Phenol',      smiles: 'c1ccc(O)cc1',       category: 'aromatic' },
    { name: 'Toluene',     smiles: 'Cc1ccccc1',         category: 'aromatic' },
    { name: 'Acetic acid', smiles: 'CC(=O)O',           category: 'acid' },
    { name: 'Aspirin',     smiles: 'CC(=O)Oc1ccccc1C(=O)O', category: 'drug' }
];

// Setup helpers ----------------------------------------------------------
var savedCM = global.COMMON_MOLECULES;
function withFixture(fn) {
    global.COMMON_MOLECULES = FIXTURE_LIBRARY.slice();
    if (typeof MolEditor !== 'undefined' && MolEditor.resetLibraryCache) {
        MolEditor.resetLibraryCache();
    }
    try { fn(); } finally {
        global.COMMON_MOLECULES = savedCM;
        if (typeof MolEditor !== 'undefined' && MolEditor.resetLibraryCache) {
            MolEditor.resetLibraryCache();
        }
    }
}

// Build a minimal MolEditor that can resolve a SMILES query without any DOM.
// We need: editor.smiles() => string and editor.searchLibrary().
function makeStubEditor(smiles) {
    var mol = new Molecule();
    SmilesParser.parse(smiles, mol);
    var editor = Object.create(MolEditor.prototype);
    editor.molecule = mol;
    editor.smiles = function() { return smiles; };
    editor.render = function() { /* DOM-less test stub */ };
    return editor;
}

// =========================================================================
// (Q) API surface
// =========================================================================

test('Q1. MolEditor.prototype.searchLibrary is a function', function() {
    assert.strictEqual(typeof MolEditor.prototype.searchLibrary, 'function');
});

test('Q2. searchLibrary returns a Promise', function() {
    withFixture(function() {
        var ed = makeStubEditor('CCO');
        var p = ed.searchLibrary('CCO', 'exact');
        assert.ok(p && typeof p.then === 'function', 'expected a Promise');
    });
});

test('Q3. searchLibrary rejects on unknown mode', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('CCO');
        ed.searchLibrary('CCO', 'bogus').then(function() {
            assert.fail('should have rejected');
        }, function(err) {
            assert.ok(/unknown mode/i.test(err.message));
            done && done();
        });
    });
});

// =========================================================================
// (R) Exact match
// =========================================================================

test('R1. Exact match: CCO finds Ethanol', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('CCO');
        ed.searchLibrary('CCO', 'exact').then(function(hits) {
            assert.ok(hits.length >= 1, 'expected at least 1 exact hit');
            var names = hits.map(function(h) { return h.name; });
            assert.ok(names.indexOf('Ethanol') >= 0, 'Ethanol missing: ' + names.join(', '));
            assert.strictEqual(hits[0].mode, 'exact');
            assert.strictEqual(hits[0].score, 1.0);
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('R2. Exact match: c1ccccc1 finds Benzene', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'exact').then(function(hits) {
            var names = hits.map(function(h) { return h.name; });
            assert.ok(names.indexOf('Benzene') >= 0, 'Benzene missing: ' + names.join(', '));
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('R3. Exact match: novel SMILES yields zero hits', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('CCNCCNCCN');
        ed.searchLibrary('CCNCCNCCN', 'exact').then(function(hits) {
            assert.strictEqual(hits.length, 0,
                'unique SMILES should not match anything in the fixture');
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

// =========================================================================
// (S) Similarity (fingerprint Tanimoto)
// =========================================================================

test('S1. Similarity: c1ccccc1 retrieves Benzene/Phenol/Toluene above threshold',
function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'similarity', { threshold: 0.2 }).then(function(hits) {
            var names = hits.map(function(h) { return h.name; });
            // The aromatic six-membered cluster should be in the top hits.
            var aromaticCount = ['Benzene', 'Phenol', 'Toluene'].filter(function(n) {
                return names.indexOf(n) >= 0;
            }).length;
            assert.ok(aromaticCount >= 2,
                'expected >=2 aromatic hits; got ' + aromaticCount + ': ' + names.join(', '));
            assert.ok(hits[0].score >= 0.2, 'top hit should clear threshold');
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('S2. Similarity: respects threshold (high threshold = fewer hits)', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'similarity', { threshold: 0.99 }).then(function(hits) {
            // Only Benzene itself should clear a near-1.0 threshold.
            assert.ok(hits.length <= 2,
                'high threshold should drop most hits; got ' + hits.length);
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('S3. Similarity: respects topN cap', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'similarity', { threshold: 0.0, topN: 2 }).then(function(hits) {
            assert.ok(hits.length <= 2);
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

// =========================================================================
// (T) Substructure (VF2 / SmartsMatch)
// =========================================================================

test('T1. Substructure: c1ccccc1 finds aromatic ring containers', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'substructure').then(function(hits) {
            var names = hits.map(function(h) { return h.name; });
            // Phenol, Toluene, Aspirin, and Benzene all contain a benzene ring.
            var aromaticContainers = ['Benzene', 'Phenol', 'Toluene', 'Aspirin'].filter(function(n) {
                return names.indexOf(n) >= 0;
            }).length;
            assert.ok(aromaticContainers >= 3,
                'expected >=3 ring-containing hits; got ' + aromaticContainers + ': ' + names.join(', '));
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('T2. Substructure: each hit carries matchCount > 0', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'substructure').then(function(hits) {
            for (var i = 0; i < hits.length; i++) {
                assert.ok(hits[i].matchCount > 0,
                    'hit ' + hits[i].name + ' should have matchCount > 0');
                assert.strictEqual(hits[i].mode, 'substructure');
            }
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

// =========================================================================
// (U) MCS
// =========================================================================

test('U1. MCS: c1ccccc1 retrieves benzene-ring-bearing molecules with score >= 0.3',
function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'mcs', { threshold: 0.3, timeoutMs: 5000 }).then(function(hits) {
            assert.ok(hits.length >= 1, 'MCS search should find at least 1 hit');
            for (var i = 0; i < hits.length; i++) {
                assert.strictEqual(hits[i].mode, 'mcs');
                assert.ok(hits[i].score >= 0.3,
                    'hit ' + hits[i].name + ' score ' + hits[i].score + ' should be >= 0.3');
                assert.ok(hits[i].mcsSize > 0);
            }
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('U2. MCS: identity match (Benzene against c1ccccc1) yields score = 1.0',
function(done) {
    withFixture(function() {
        var ed = makeStubEditor('c1ccccc1');
        ed.searchLibrary('c1ccccc1', 'mcs', { threshold: 0.5 }).then(function(hits) {
            // The Benzene entry should appear with a near-1.0 score.
            var benzeneHit = hits.filter(function(h) { return h.name === 'Benzene'; })[0];
            assert.ok(benzeneHit, 'Benzene should be in hits');
            assert.ok(benzeneHit.score >= 0.99,
                'Benzene MCS score should be ~1.0; got ' + benzeneHit.score);
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

// =========================================================================
// (V) Library cache + invalidation
// =========================================================================

test('V1. resetLibraryCache forces re-parse on next call', function(done) {
    withFixture(function() {
        var ed = makeStubEditor('CCO');
        ed.searchLibrary('CCO', 'exact').then(function(hits1) {
            assert.strictEqual(hits1[0].name, 'Ethanol');
            // Mutate the library and reset cache.
            global.COMMON_MOLECULES = [{ name: 'Replaced', smiles: 'CCO', category: 'basic' }];
            MolEditor.resetLibraryCache();
            return ed.searchLibrary('CCO', 'exact');
        }).then(function(hits2) {
            assert.strictEqual(hits2[0].name, 'Replaced',
                'after resetLibraryCache, the new library should be used');
            done && done();
        }).catch(function(e) { assert.fail(e); });
    });
});

test('V2. searchLibrary rejects when COMMON_MOLECULES is empty', function(done) {
    var saved = global.COMMON_MOLECULES;
    global.COMMON_MOLECULES = [];
    if (MolEditor.resetLibraryCache) MolEditor.resetLibraryCache();
    try {
        var ed = makeStubEditor('CCO');
        ed.searchLibrary('CCO', 'exact').then(function() {
            assert.fail('should reject on empty library');
        }, function(err) {
            assert.ok(/library not loaded/i.test(err.message),
                'error should mention library; got ' + err.message);
            done && done();
        });
    } finally {
        global.COMMON_MOLECULES = saved;
        if (MolEditor.resetLibraryCache) MolEditor.resetLibraryCache();
    }
});

test('V3. RDT.version stamp matches current bundle', function() {
    assert.strictEqual(typeof RDT, 'object');
    assert.strictEqual(RDT.version, '1.8.15');
});

// =========================================================================
// (W) User-supplied target library (options.targets)
// =========================================================================

test('W1. options.targets accepts a custom 5-entry list and searches it instead',
function(done) {
    var saved = global.COMMON_MOLECULES;
    global.COMMON_MOLECULES = []; // built-in disabled to prove targets path
    if (MolEditor.resetLibraryCache) MolEditor.resetLibraryCache();
    try {
        var ed = makeStubEditor('CCO');
        var customTargets = [
            { name: 'Custom-Methane', smiles: 'C' },
            { name: 'Custom-Ethanol', smiles: 'CCO' },
            { name: 'Custom-Benzene', smiles: 'c1ccccc1' }
        ];
        ed.searchLibrary('CCO', 'exact', { targets: customTargets }).then(function(hits) {
            assert.strictEqual(hits.length, 1, 'expected 1 exact hit in user lib');
            assert.strictEqual(hits[0].name, 'Custom-Ethanol');
            done && done();
        }).catch(function(e) { assert.fail(e); });
    } finally {
        global.COMMON_MOLECULES = saved;
        if (MolEditor.resetLibraryCache) MolEditor.resetLibraryCache();
    }
});

test('W2. options.targets caps at MolEditor.USER_LIBRARY_LIMIT (100)',
function(done) {
    var ed = makeStubEditor('CCO');
    // Build 150 entries; only the first 100 should be searched.
    var bigTargets = [];
    for (var i = 0; i < 150; i++) {
        bigTargets.push({ name: 'M' + i, smiles: 'CCO' });
    }
    ed.searchLibrary('CCO', 'exact', { targets: bigTargets, topN: 200 }).then(function(hits) {
        assert.ok(hits.length <= 100,
            'user library cap should drop entries beyond the 100th; got ' + hits.length);
        done && done();
    }).catch(function(e) { assert.fail(e); });
});

test('W3. options.targets rejects when the user list is empty', function(done) {
    var ed = makeStubEditor('CCO');
    ed.searchLibrary('CCO', 'exact', { targets: [] }).then(function() {
        assert.fail('should reject on empty user library');
    }, function(err) {
        assert.ok(/empty or unparsable/i.test(err.message),
            'error should mention empty/unparsable; got ' + err.message);
        done && done();
    });
});

test('W4. MolEditor.USER_LIBRARY_LIMIT exposes the cap as a public constant',
function() {
    assert.strictEqual(MolEditor.USER_LIBRARY_LIMIT, 100);
});

// =========================================================================
// (X) Match highlighting on hit click (highlightSearchMatch)
// =========================================================================

test('X1. highlightSearchMatch substructure: ring atoms get atom.highlighted = true',
function() {
    // Load Phenol; query is c1ccccc1 (benzene ring substructure).
    var ed = makeStubEditor('Oc1ccccc1');
    var n = ed.highlightSearchMatch('substructure', 'c1ccccc1');
    assert.ok(n >= 6, 'expected >=6 atoms highlighted (the ring); got ' + n);
    var hi = 0;
    for (var i = 0; i < ed.molecule.atoms.length; i++) {
        if (ed.molecule.atoms[i].highlighted) hi++;
    }
    assert.ok(hi >= 6, 'molecule should carry >=6 highlighted atoms after call');
});

test('X2. highlightSearchMatch exact: every atom gets highlighted', function() {
    var ed = makeStubEditor('CCO');
    var n = ed.highlightSearchMatch('exact', 'CCO');
    assert.strictEqual(n, ed.molecule.atoms.length,
        'exact highlight should mark every atom');
});

test('X3. highlightSearchMatch similarity: zero atoms highlighted (whole-mol score)',
function() {
    var ed = makeStubEditor('CCO');
    var n = ed.highlightSearchMatch('similarity', 'CCO');
    assert.strictEqual(n, 0,
        'similarity highlight is intentionally a no-op (per-atom Tanimoto undefined)');
});

test('X4. highlightSearchMatch mcs: at least 1 matched atom for shared scaffold',
function() {
    // Loaded molecule is Phenol; query is benzene. MCS should be the ring.
    var ed = makeStubEditor('Oc1ccccc1');
    var n = ed.highlightSearchMatch('mcs', 'c1ccccc1');
    assert.ok(n >= 6, 'MCS highlight should mark at least the 6-atom ring; got ' + n);
});

test('X5. highlightSearchMatch resets prior highlights before applying new ones',
function() {
    var ed = makeStubEditor('Oc1ccccc1');
    // Pre-set atom 0 as highlighted with bgColor.
    ed.molecule.atoms[0].highlighted = true;
    ed.molecule.atoms[0].bgColor = '#ff0000';
    // Now run similarity highlight (should be a no-op AND clear the prior highlight).
    ed.highlightSearchMatch('similarity', 'CCO');
    assert.strictEqual(ed.molecule.atoms[0].highlighted, false,
        'prior highlight should be cleared before re-applying');
    assert.strictEqual(ed.molecule.atoms[0].bgColor, null,
        'prior bgColor should be cleared before re-applying');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
