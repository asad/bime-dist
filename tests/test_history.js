/**
 * tests/test_history.js — Undo/redo ring-buffer + reaction-state regressions.
 *
 * Locks in:
 *   - Cap=100 ring buffer (oldest evicted on overflow)
 *   - Push clears redo stack
 *   - Reaction undo/redo round-trips reactionPlusSigns + arrow type/conditions
 *     (the v1.0.2 Molecule.toJSON / clone() fix)
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('History');
var test = runner.test;

console.log('History (undo/redo)');

function makeMol(seed) {
    var m = new Molecule();
    // Add a single distinct atom so toJSON round-trips uniquely per seed.
    var a = new Molecule.Atom('C', seed, 0);
    m.atoms.push(a);
    m._atomMap[a.id] = a;
    m._adjacency[a.id] = [];
    m.name = 'mol-' + seed;
    return m;
}

// Replicate MolEditor._restoreMolecule (without canvas/renderer plumbing) so
// the test can step through undo/redo with real Molecule instances.
function restore(mol, json) {
    mol.clear();
    if (json.name !== undefined) mol.name = json.name;
    var idMap = {};
    for (var i = 0; i < json.atoms.length; i++) {
        var a = json.atoms[i];
        var atom = mol.addAtom(a.symbol, a.x, a.y);
        atom.charge = a.charge || 0;
        atom.isotope = a.isotope || 0;
        atom.mapNumber = a.mapNumber || 0;
        atom.hydrogens = a.hydrogens !== undefined ? a.hydrogens : -1;
        if (a.aromatic) atom.aromatic = true;
        idMap[a.id] = atom.id;
    }
    for (var i = 0; i < json.bonds.length; i++) {
        var b = json.bonds[i];
        var a1 = (b.atom1 in idMap) ? idMap[b.atom1] : b.atom1;
        var a2 = (b.atom2 in idMap) ? idMap[b.atom2] : b.atom2;
        var bond = mol.addBond(a1, a2, b.type);
        if (bond) bond.stereo = b.stereo || 0;
    }
    if (json.reactionArrow) mol.reactionArrow = json.reactionArrow;
    if (json.reactionPlusSigns) mol.reactionPlusSigns = json.reactionPlusSigns;
    return mol;
}

test('push 100 entries, undo back to start, redo all the way forward', function() {
    var h = new EditorHistory();
    var live = makeMol(0);
    for (var i = 1; i <= 100; i++) {
        h.push(live);
        live = makeMol(i);
    }
    assert.strictEqual(h.canUndo(), true);
    assert.strictEqual(h.canRedo(), false);

    // Undo all 100 — undo() returns the popped JSON snapshot and pushes the
    // current molecule's JSON onto the redo stack.
    var u = 0;
    while (h.canUndo()) {
        var snap = h.undo(live);
        restore(live, snap);
        u++;
    }
    assert.strictEqual(u, 100);
    assert.strictEqual(h.canRedo(), true);

    var r = 0;
    while (h.canRedo()) {
        var snap = h.redo(live);
        restore(live, snap);
        r++;
    }
    assert.strictEqual(r, 100);
});

test('push 200 entries — oldest 100 evicted, length capped at 100', function() {
    var h = new EditorHistory();
    var live = makeMol(0);
    for (var i = 1; i <= 200; i++) {
        h.push(live);
        live = makeMol(i);
    }
    assert.strictEqual(h._undoStack.length, 100);
    // Undo all 100 retained snapshots; each snapshot's name encodes its seed.
    var seenNames = [];
    while (h.canUndo()) {
        var snap = h.undo(live);
        seenNames.push(snap.name);
        restore(live, snap);
    }
    // We pushed snapshots of seed-1..seed-200 (cur snapshot is taken before
    // each `live = makeMol(i)` advance). After cap=100 eviction, the surviving
    // window is the last 100 snapshots = mol-100 .. mol-199. LIFO undo pops
    // newest first.
    assert.strictEqual(seenNames.length, 100);
    assert.strictEqual(seenNames[0], 'mol-199');
    assert.strictEqual(seenNames[99], 'mol-100');
});

test('push, then changed-without-push (push again) clears redo stack', function() {
    var h = new EditorHistory();
    var a = makeMol(1);
    var b = makeMol(2);
    var c = makeMol(3);
    h.push(a);
    h.push(b);
    var undoneSnap = h.undo(c); // moves b's snapshot pop, pushes c to redo
    assert.strictEqual(h.canRedo(), true);
    h.push(c); // new push must clear redo
    assert.strictEqual(h.canRedo(), false);
});

test('reaction undo includes reactionPlusSigns (v1.0.2 fix)', function() {
    var m = new Molecule();
    m.reactionArrow = { x1: 0, y1: 0, x2: 100, y2: 0, type: 'forward', conditions: 'heat' };
    m.reactionPlusSigns = [{ x: 50, y: 10 }, { x: 150, y: 10 }];
    var snap = m.toJSON();
    assert.ok(snap.reactionArrow, 'arrow lost in toJSON');
    assert.strictEqual(snap.reactionArrow.type, 'forward');
    assert.strictEqual(snap.reactionArrow.conditions, 'heat');
    assert.ok(Array.isArray(snap.reactionPlusSigns));
    assert.strictEqual(snap.reactionPlusSigns.length, 2);
    assert.strictEqual(snap.reactionPlusSigns[0].x, 50);

    // clone() must also preserve them.
    var c = m.clone();
    assert.ok(c.reactionPlusSigns);
    assert.strictEqual(c.reactionPlusSigns.length, 2);
    assert.strictEqual(c.reactionArrow.type, 'forward');
    assert.strictEqual(c.reactionArrow.conditions, 'heat');
});

test('clear() empties both undo and redo stacks', function() {
    var h = new EditorHistory();
    h.push(makeMol(1));
    h.push(makeMol(2));
    h.undo(makeMol(3));
    h.clear();
    assert.strictEqual(h.canUndo(), false);
    assert.strictEqual(h.canRedo(), false);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
