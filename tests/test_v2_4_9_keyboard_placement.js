/**
 * tests/test_v2_4_9_keyboard_placement.js — keyboard atom placement / build mode
 * (v2.4.9). Finishes the v2.3.7 accessibility arc (Phase B).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.3.7 made the draw canvas REACHABLE + DESCRIBED (role=img) but explicitly
 * deferred keyboard OPERABILITY for atom placement, noting that a premature
 * role=application — promising key handling the surface did not implement — is
 * worse for a screen-reader user than an honest labelled graphic. v2.4.9
 * delivers that handling: a mode-gated "keyboard build mode" with a
 * directional-growth model. Enter on the focused canvas enters the mode (role
 * flips img->application ONLY then); arrow keys grow a bonded atom from the
 * current atom (or bond to an atom already there — ring closure); letter keys set
 * the element; 1/2/3 set the bond order; Delete removes the current atom; Escape
 * leaves the mode (role reverts to img). Each step is announced through a polite
 * aria-live region. The "current atom" reuses the existing selection highlight.
 *
 * This suite drives the REAL MolEditor._kb* methods against a REAL Molecule (the
 * graph mutation needs no DOM), exercises the role flip through the real
 * enter/exit methods with a recording stub, and pins the keydown routing +
 * mode-gating as a source-shape contract.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'MolEditor.js'));
var MolEditor = globalThis.MolEditor;

var runner = shim.makeRunner('Keyboard atom placement (v2.4.9)');
var test = runner.test;
console.log('Keyboard atom placement (v2.4.9)');

var SRC = fs.readFileSync(path.join(__dirname, '..', 'editor', 'MolEditor.js'), 'utf8');

// A minimal MolEditor-like host: a real Molecule + no-op DOM/editor hooks, with
// the real _kb* prototype methods bound on. Lets us exercise the actual graph
// logic without instantiating the full DOM editor.
function makeHost() {
    var host = {
        molecule: new globalThis.Molecule(),
        saveHistory: function () {}, changed: function () {}, render: function () {},
        canvasDescription: function () { return ''; }, showInfo: function () {},
        _updateCanvasDescription: function () {},
        _kbLive: null, _kbElement: 'C', _kbBondType: 1, _kbCurrentAtomId: null,
        _kbLastDir: { dx: 1, dy: 0 }, _kbBuildMode: false
    };
    ['_kbSetCurrent', '_kbAtomNear', '_kbAreBonded', '_kbAnnounce', '_kbGrow',
     '_kbSetElementForBuild', '_kbDeleteCurrent', '_enterKeyboardBuild', '_exitKeyboardBuild'
    ].forEach(function (m) { host[m] = MolEditor.prototype[m]; });
    return host;
}

// ---------------------------------------------------------------------------
// A. Real graph logic — directional growth, ring closure, retype, delete.
// ---------------------------------------------------------------------------
test('A1 first growth on an empty canvas places one atom and makes it current', function () {
    var h = makeHost();
    h._kbGrow({ dx: 1, dy: 0 });
    assert.strictEqual(h.molecule.atoms.length, 1, 'one atom placed');
    assert.ok(h._kbCurrentAtomId != null, 'the placed atom is current');
    assert.ok(h.molecule.getAtom(h._kbCurrentAtomId).selected, 'the current atom is selected (the visible cursor)');
});

test('A2 arrow growth adds a bonded atom in the given direction', function () {
    var h = makeHost();
    h._kbGrow({ dx: 1, dy: 0 });   // first atom at (0,0)
    var firstId = h._kbCurrentAtomId;
    h._kbGrow({ dx: 1, dy: 0 });   // grow right
    assert.strictEqual(h.molecule.atoms.length, 2, 'a second atom was added');
    assert.strictEqual(h.molecule.bonds.length, 1, 'they are bonded');
    assert.notStrictEqual(h._kbCurrentAtomId, firstId, 'the new atom becomes current');
    var nu = h.molecule.getAtom(h._kbCurrentAtomId);
    assert.ok(Math.abs(nu.x - 30) < 1e-6, 'placed one bond-length to the right');
});

test('A3 growing back into an existing atom closes a ring (bond, no new atom)', function () {
    var h = makeHost();
    h._kbGrow({ dx: 1, dy: 0 });   // a0 (0,0)
    h._kbGrow({ dx: 1, dy: 0 });   // a1 (30,0)
    h._kbGrow({ dx: 0, dy: 1 });   // a2 (30,30)
    h._kbGrow({ dx: -1, dy: 0 });  // a3 (0,30)
    assert.strictEqual(h.molecule.atoms.length, 4, '4 atoms');
    assert.strictEqual(h.molecule.bonds.length, 3, '3 bonds (open chain)');
    h._kbGrow({ dx: 0, dy: -1 });  // up -> lands on a0 (0,0): ring closure
    assert.strictEqual(h.molecule.atoms.length, 4, 'no new atom on ring closure');
    assert.strictEqual(h.molecule.bonds.length, 4, 'a closing bond was added');
});

test('A4 a letter key retypes the current atom and sets the pending element', function () {
    var h = makeHost();
    h._kbGrow({ dx: 1, dy: 0 });
    h._kbSetElementForBuild('O');
    assert.strictEqual(h.molecule.getAtom(h._kbCurrentAtomId).symbol, 'O', 'current atom retyped to O');
    assert.strictEqual(h._kbElement, 'O', 'pending element is O for the next growth');
    h._kbGrow({ dx: 1, dy: 0 });
    assert.strictEqual(h.molecule.getAtom(h._kbCurrentAtomId).symbol, 'O', 'next grown atom uses the pending element');
});

test('A5 the bond order is honoured by growth', function () {
    var h = makeHost();
    h._kbGrow({ dx: 1, dy: 0 });
    h._kbBondType = 2;
    h._kbGrow({ dx: 1, dy: 0 });
    assert.strictEqual(h.molecule.bonds[h.molecule.bonds.length - 1].type, 2, 'the new bond is a double bond');
});

test('A6 delete removes the current atom and re-anchors on a neighbour', function () {
    var h = makeHost();
    h._kbGrow({ dx: 1, dy: 0 });   // a0
    var a0 = h._kbCurrentAtomId;
    h._kbGrow({ dx: 1, dy: 0 });   // a1 (current)
    h._kbDeleteCurrent();
    assert.strictEqual(h.molecule.atoms.length, 1, 'the current atom was removed');
    assert.strictEqual(h._kbCurrentAtomId, a0, 'the surviving neighbour becomes current');
});

// ---------------------------------------------------------------------------
// B. Role flip through the real enter/exit methods (mode-gated application).
// ---------------------------------------------------------------------------
test('B1 entering build mode flips the canvas role to application; exiting reverts to img', function () {
    var roles = [];
    var stubArea = {
        _attrs: {},
        setAttribute: function (k, v) { this._attrs[k] = v; if (k === 'role') roles.push(v); },
        getAttribute: function (k) { return this._attrs[k]; },
        focus: function () {}
    };
    var h = makeHost();
    h._drawArea = stubArea;
    h._enterKeyboardBuild();
    assert.strictEqual(h._kbBuildMode, true, 'build mode is on');
    assert.strictEqual(stubArea.getAttribute('role'), 'application', 'role flipped to application while building');
    h._exitKeyboardBuild();
    assert.strictEqual(h._kbBuildMode, false, 'build mode is off');
    assert.strictEqual(stubArea.getAttribute('role'), 'img', 'role reverts to img on exit');
    assert.deepStrictEqual(roles, ['application', 'img'], 'role flipped application then back to img');
});

// ---------------------------------------------------------------------------
// C. Source-shape: keydown routing, mode-gating, default role, live region.
// ---------------------------------------------------------------------------
test('C1 the canvas defaults to role=img (the v2.3.7 honest-graphic contract is kept)', function () {
    assert.ok(/setAttribute\('role',\s*'img'\)/.test(SRC),
        'the draw area is still constructed as role=img by default (application is mode-gated)');
});

test('C2 the keydown handler routes build keys ONLY when the canvas is focused and no Ctrl/Cmd', function () {
    assert.ok(/e\.target === self\._drawArea && !e\.ctrlKey && !e\.metaKey/.test(SRC),
        'build-mode keys are gated on the draw canvas having focus, excluding Ctrl/Cmd chords');
    assert.ok(/if \(!self\._kbBuildMode\)[\s\S]{0,80}self\._enterKeyboardBuild\(\)/.test(SRC),
        'Enter on the focused canvas (not yet in build mode) enters keyboard build');
    assert.ok(SRC.indexOf('self._kbGrow(KBDIRS[e.key])') !== -1, 'arrow keys route to _kbGrow');
    assert.ok(SRC.indexOf('self._exitKeyboardBuild()') !== -1, 'Escape exits build mode');
    assert.ok(SRC.indexOf('self._kbDeleteCurrent()') !== -1, 'Delete removes the current atom');
});

test('C3 a polite aria-live region is created for the announcements', function () {
    assert.ok(/setAttribute\('aria-live',\s*'polite'\)/.test(SRC), 'a polite aria-live region exists');
    assert.ok(/this\._kbLive\s*=\s*document\.createElement/.test(SRC), 'the live region is a real element (_kbLive)');
});

test('C4 every keyboard mutation goes through saveHistory (so Ctrl+Z undoes it)', function () {
    // _kbGrow / _kbSetElementForBuild / _kbDeleteCurrent each saveHistory before mutating.
    var grow = SRC.slice(SRC.indexOf('prototype._kbGrow'), SRC.indexOf('prototype._kbSetElementForBuild'));
    assert.ok(grow.indexOf('this.saveHistory()') !== -1, '_kbGrow saves an undo frame');
    var del = SRC.slice(SRC.indexOf('prototype._kbDeleteCurrent'), SRC.indexOf('prototype._kbDeleteCurrent') + 700);
    assert.ok(del.indexOf('this.saveHistory()') !== -1, '_kbDeleteCurrent saves an undo frame');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
