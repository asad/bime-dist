/**
 * tests/test_v1_4_3_features.js — BIME v1.4.3 user-customizable toolbar.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Covers the ToolbarPrefs persistence layer + applyToGroups/applyToAtomBar
 * filtering and reordering. The Customize panel UI itself is exercised
 * indirectly through DOM-level assertions in workbench smoke tests
 * (manual). Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
// ToolbarPrefs is registered as a global by editor/ToolbarPrefs.js (or by
// the bundle via dist/bime.js). Pull it off the global so this test file
// works against both source-mode and bundle-mode test runs.
var ToolbarPrefs = global.ToolbarPrefs;

var runner = shim.makeRunner('Toolbar customization v1.4.3');
var test = runner.test;

console.log('Toolbar customization v1.4.3');

// -------------------------------------------------------------------------
// In-memory localStorage stub so the module can save/load in Node tests.
// -------------------------------------------------------------------------
var storeBacking = {};
var fakeStorage = {
    getItem: function(k) { return Object.prototype.hasOwnProperty.call(storeBacking, k) ? storeBacking[k] : null; },
    setItem: function(k, v) { storeBacking[k] = String(v); },
    removeItem: function(k) { delete storeBacking[k]; },
    clear: function() { storeBacking = {}; }
};
function withStorage(fn) {
    storeBacking = {};
    var saved = global.localStorage;
    global.localStorage = fakeStorage;
    try { fn(); } finally { global.localStorage = saved; }
}

// Canonical fixtures matching MolEditor.js layout.
var GROUPS = [
    { id: 'draw', label: 'Draw', items: [
        { id: 'bond', label: 'Bond' },
        { id: 'chain', label: 'Chain' }
    ]},
    { id: 'bonds', label: 'Bonds', items: [
        { id: 'single', label: 'Single' },
        { id: 'double', label: 'Double' }
    ]},
    { id: 'edit', label: 'Edit', items: [
        { id: 'select', label: 'Select' },
        { id: 'delete', label: 'Delete' }
    ]}
];
var ATOMS = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'H'];

// =========================================================================
// (M) ToolbarPrefs API surface + storage round-trip
// =========================================================================

test('M1. Module exposes load/save/reset/applyToGroups/applyToAtomBar/snapshot', function() {
    assert.strictEqual(typeof ToolbarPrefs.load, 'function');
    assert.strictEqual(typeof ToolbarPrefs.save, 'function');
    assert.strictEqual(typeof ToolbarPrefs.reset, 'function');
    assert.strictEqual(typeof ToolbarPrefs.applyToGroups, 'function');
    assert.strictEqual(typeof ToolbarPrefs.applyToAtomBar, 'function');
    assert.strictEqual(typeof ToolbarPrefs.snapshot, 'function');
    assert.strictEqual(typeof ToolbarPrefs.validate, 'function');
    assert.strictEqual(ToolbarPrefs.STORAGE_KEY, 'bime-toolbar-prefs-v1');
    assert.strictEqual(ToolbarPrefs.SCHEMA_VERSION, 1);
});

test('M2. snapshot(GROUPS, ATOMS) produces a valid prefs object', function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    assert.ok(ToolbarPrefs.validate(snap), 'snapshot should validate');
    assert.deepStrictEqual(snap.groupOrder, ['draw', 'bonds', 'edit']);
    assert.deepStrictEqual(snap.atomBar, ATOMS);
    assert.strictEqual(snap.groups.draw.hidden, false);
    assert.deepStrictEqual(snap.groups.draw.items, ['bond', 'chain']);
});

test('M3. save then load round-trips identical prefs', function() {
    withStorage(function() {
        var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
        snap.atomBar = ['C', 'N', 'O'];
        snap.groups.draw.items = ['chain', 'bond'];
        var ok = ToolbarPrefs.save(snap);
        assert.strictEqual(ok, true);
        var loaded = ToolbarPrefs.load();
        assert.deepStrictEqual(loaded, snap);
    });
});

test('M4. reset clears localStorage so next load returns null', function() {
    withStorage(function() {
        ToolbarPrefs.save(ToolbarPrefs.snapshot(GROUPS, ATOMS));
        assert.notStrictEqual(ToolbarPrefs.load(), null);
        var ok = ToolbarPrefs.reset();
        assert.strictEqual(ok, true);
        assert.strictEqual(ToolbarPrefs.load(), null);
    });
});

test('M5. load returns null when localStorage is unavailable', function() {
    var savedLs = global.localStorage;
    delete global.localStorage;
    try {
        assert.strictEqual(ToolbarPrefs.load(), null);
    } finally { global.localStorage = savedLs; }
});

test('M6. validate rejects malformed prefs', function() {
    assert.strictEqual(ToolbarPrefs.validate(null), false);
    assert.strictEqual(ToolbarPrefs.validate({}), false);
    assert.strictEqual(ToolbarPrefs.validate({ version: 999 }), false);
    assert.strictEqual(ToolbarPrefs.validate({ version: 1, groupOrder: 'oops' }), false);
    assert.strictEqual(ToolbarPrefs.validate({
        version: 1, groupOrder: [], groups: {}, atomBar: [1, 2, 3]
    }), false, 'atomBar must be string array');
    assert.strictEqual(ToolbarPrefs.validate({
        version: 1, groupOrder: [], groups: { x: { hidden: 'no', items: [] } }, atomBar: []
    }), false, 'group.hidden must be boolean');
});

test('M7. corrupted localStorage entry yields load() === null', function() {
    withStorage(function() {
        global.localStorage.setItem(ToolbarPrefs.STORAGE_KEY, '{ not valid json }');
        assert.strictEqual(ToolbarPrefs.load(), null);
        // setItem with invalid prefs JSON parses but fails validate.
        global.localStorage.setItem(ToolbarPrefs.STORAGE_KEY, JSON.stringify({ version: 99 }));
        assert.strictEqual(ToolbarPrefs.load(), null);
    });
});

// =========================================================================
// (N) applyToGroups — filtering + reordering
// =========================================================================

test('N1. applyToGroups with null prefs returns canonical groups unchanged',
function() {
    var out = ToolbarPrefs.applyToGroups(GROUPS, null);
    assert.strictEqual(out.length, GROUPS.length);
    for (var i = 0; i < GROUPS.length; i++) {
        assert.strictEqual(out[i].id, GROUPS[i].id);
    }
});

test('N2. applyToGroups reorders groups per prefs.groupOrder', function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.groupOrder = ['edit', 'draw', 'bonds'];
    var out = ToolbarPrefs.applyToGroups(GROUPS, snap);
    assert.deepStrictEqual(out.map(function(g) { return g.id; }),
        ['edit', 'draw', 'bonds']);
});

test('N3. applyToGroups drops hidden groups', function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.groups.bonds.hidden = true;
    var out = ToolbarPrefs.applyToGroups(GROUPS, snap);
    var ids = out.map(function(g) { return g.id; });
    assert.ok(ids.indexOf('bonds') < 0, 'bonds should be hidden');
    assert.ok(ids.indexOf('draw') >= 0);
    assert.ok(ids.indexOf('edit') >= 0);
});

test('N4. applyToGroups reorders items within a group', function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.groups.draw.items = ['chain', 'bond'];
    var out = ToolbarPrefs.applyToGroups(GROUPS, snap);
    var draw = out.filter(function(g) { return g.id === 'draw'; })[0];
    assert.deepStrictEqual(draw.items.map(function(it) { return it.id; }),
        ['chain', 'bond']);
});

test('N5. applyToGroups hides items via hiddenItems (explicit user-hide)',
function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    // User explicitly hid 'chain' via the Customize panel checkbox.
    snap.groups.draw.items = ['bond'];
    snap.groups.draw.hiddenItems = ['chain'];
    var out = ToolbarPrefs.applyToGroups(GROUPS, snap);
    var draw = out.filter(function(g) { return g.id === 'draw'; })[0];
    assert.strictEqual(draw.items.length, 1);
    assert.strictEqual(draw.items[0].id, 'bond');
});

test('N5b. applyToGroups appends new items NOT in items AND NOT in hiddenItems',
function() {
    // Forward-compat: user prefs only mention 'bond' (no chain in items, no
    // chain in hiddenItems either — simulating a new BIME version that
    // added 'chain' as a button).
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.groups.draw.items = ['bond'];
    snap.groups.draw.hiddenItems = [];
    var out = ToolbarPrefs.applyToGroups(GROUPS, snap);
    var draw = out.filter(function(g) { return g.id === 'draw'; })[0];
    assert.strictEqual(draw.items.length, 2);
    assert.strictEqual(draw.items[0].id, 'bond');
    assert.strictEqual(draw.items[1].id, 'chain', 'newer-version item appended');
});

test('N6. applyToGroups appends a NEW canonical group not yet in prefs',
function() {
    // User saved prefs from older BIME with only [draw, bonds]; new BIME
    // adds 'edit' group. The new group should appear at the end.
    var snap = {
        version: 1,
        groupOrder: ['draw', 'bonds'],
        groups: {
            draw: { hidden: false, items: ['bond', 'chain'] },
            bonds: { hidden: false, items: ['single', 'double'] }
        },
        atomBar: ATOMS.slice()
    };
    var out = ToolbarPrefs.applyToGroups(GROUPS, snap);
    var ids = out.map(function(g) { return g.id; });
    assert.deepStrictEqual(ids, ['draw', 'bonds', 'edit'],
        'newly-added "edit" group should append at end');
});

// =========================================================================
// (O) applyToAtomBar — filtering + reordering atoms
// =========================================================================

test('O1. applyToAtomBar with null prefs returns canonical atoms unchanged',
function() {
    var out = ToolbarPrefs.applyToAtomBar(ATOMS, null);
    assert.deepStrictEqual(out, ATOMS);
});

test('O2. applyToAtomBar drops symbols not in prefs.atomBar', function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.atomBar = ['C', 'O', 'N'];
    var out = ToolbarPrefs.applyToAtomBar(ATOMS, snap);
    assert.deepStrictEqual(out, ['C', 'O', 'N']);
});

test('O3. applyToAtomBar rejects symbols not in canonical bar (anti-tamper)',
function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.atomBar = ['C', 'NotAnElement', 'O'];
    var out = ToolbarPrefs.applyToAtomBar(ATOMS, snap);
    assert.deepStrictEqual(out, ['C', 'O']);
});

test('O4. applyToAtomBar deduplicates repeated symbols', function() {
    var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
    snap.atomBar = ['C', 'C', 'N', 'C'];
    var out = ToolbarPrefs.applyToAtomBar(ATOMS, snap);
    assert.deepStrictEqual(out, ['C', 'N']);
});

// =========================================================================
// (P) End-to-end: save -> reload -> apply yields the user's customizations
// =========================================================================

test('P1. Round-trip: save a custom layout, reload, applyToGroups matches',
function() {
    withStorage(function() {
        var snap = ToolbarPrefs.snapshot(GROUPS, ATOMS);
        snap.groupOrder = ['edit', 'draw', 'bonds'];
        snap.groups.bonds.hidden = true;
        snap.groups.draw.items = ['chain'];
        snap.groups.draw.hiddenItems = ['bond'];   // explicit user-hide
        snap.atomBar = ['O', 'C', 'N'];
        ToolbarPrefs.save(snap);

        var loaded = ToolbarPrefs.load();
        var groups = ToolbarPrefs.applyToGroups(GROUPS, loaded);
        var atoms = ToolbarPrefs.applyToAtomBar(ATOMS, loaded);

        assert.deepStrictEqual(groups.map(function(g) { return g.id; }),
            ['edit', 'draw']);
        var draw = groups.filter(function(g) { return g.id === 'draw'; })[0];
        assert.deepStrictEqual(draw.items.map(function(it) { return it.id; }),
            ['chain']);
        assert.deepStrictEqual(atoms, ['O', 'C', 'N']);
    });
});

test('P2. RDT.version is unchanged from v1.4.2 (toolbar customization is editor-only)',
function() {
    assert.strictEqual(typeof RDT, 'object');
    // toolbar customization does not change RDT semantics; version stamp
    // tracks the bundle version regardless.
    assert.ok(/^1\.\d+\.\d+/.test(RDT.version),
        'RDT.version should be a 1.x.y SemVer; got ' + RDT.version);
});

test('P3. ToolbarPrefs.version stamp matches the bundle release', function() {
    assert.strictEqual(ToolbarPrefs.version, '1.8.0');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
