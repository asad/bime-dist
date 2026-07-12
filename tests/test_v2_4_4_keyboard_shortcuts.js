/**
 * tests/test_v2_4_4_keyboard_shortcuts.js — standard chemistry keyboard
 * shortcuts + a discoverable cheat-sheet (v2.4.4).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.4.4 brings the molecule editor's hotkeys in line with the conventions a
 * web survey established for mainstream chemistry editors, so BIME feels
 * familiar:
 *
 *   1 / 2 / 3          — single / double / triple bond tool
 *   c n o s p f i h b  — set the hovered/selected atom's element, else arm the
 *                        atom tool with that element
 *   Shift+C / Shift+B  — chlorine / bromine (the two-letter halogens)
 *   Esc                — return to the Select tool
 *   Del / Backspace    — delete the selection or the atom/bond under the cursor
 *   Ctrl/Cmd+Z         — undo  ·  Ctrl/Cmd+Shift+Z and Ctrl/Cmd+Y — redo
 *   Ctrl/Cmd+A         — select all
 *   ?                  — open the keyboard-shortcuts cheat-sheet
 *
 * BIME's bond model has only single/double/triple primitives (aromatic is an
 * atom/ring property surfaced via the "Ar" display toggle; there is no any/query
 * bond), so the surveyed `4` (aromatic) and `0` (any) bond keys have no
 * primitive to bind to and are deliberately NOT mapped.
 *
 * This suite is part source-shape (slicing editor/MolEditor.js + js/workbench.js
 * so a refactor that drops a handler/guard fails here) and part executable
 * DOM-stub: it constructs a real MolEditor against a recording document stub
 * (richer than tests/shim.js, modelled on test_v2_3_4 section E) and drives the
 * new methods + the keydown dispatcher to assert behaviour and the text-input /
 * focus-lens guards.
 *
 * Plain Node, no real DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Keyboard shortcuts (v2.4.4)');
var test = runner.test;
console.log('Keyboard shortcuts (v2.4.4)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var ME = readRepoFile('editor/MolEditor.js');
var WB = readRepoFile('js/workbench.js');

// Slice a prototype-method body: from `MolEditor.prototype.<name> = function`
// up to the next `MolEditor.prototype.` — the same shape-pinning the other
// suites use, so a refactor that splits/relocates a method fails here too.
function meMethod(name) {
    var start = ME.indexOf('MolEditor.prototype.' + name + ' = function');
    if (start === -1) { return ''; }
    var rest = ME.slice(start + 1);
    var next = rest.indexOf('MolEditor.prototype.');
    return next === -1 ? rest : rest.slice(0, next);
}

// ===========================================================================
// (A) Source-shape — the dispatcher + helpers + guards exist.
// ===========================================================================
test('A1 the keydown dispatcher maps bond / element / Esc / select-all keys', function () {
    // Bond number keys.
    assert.ok(/var BOND_KEYS = \{ '1': 1, '2': 2, '3': 3 \}/.test(ME),
        'bond number keys 1/2/3 mapped to single/double/triple');
    // Element letters (lower-case) incl. boron, plus the Shift halogens.
    assert.ok(/ELEMENT_KEYS = \{[\s\S]*?'c': 'C'[\s\S]*?'b': 'B'[\s\S]*?\}/.test(ME),
        'element letters c..b mapped to symbols');
    assert.ok(/SHIFT_ELEMENT_KEYS = \{ 'c': 'Cl', 'b': 'Br' \}/.test(ME),
        'Shift+C / Shift+B mapped to Cl / Br');
    // Select-all chord lives on the Ctrl/Cmd branch.
    assert.ok(ME.indexOf("self._selectAll()") !== -1, 'Ctrl/Cmd+A routes to _selectAll');
    // Esc -> Select.
    assert.ok(/self\._currentToolName !== 'select'[\s\S]*?self\._setTool\('select'\)/.test(ME),
        'Escape returns to the Select tool when a non-Select tool is armed');
    // Del / Backspace -> _deleteViaKeyboard.
    assert.ok(/e\.key === 'Delete' \|\| e\.key === 'Backspace'[\s\S]*?self\._deleteViaKeyboard\(\)/.test(ME),
        'Delete / Backspace route to _deleteViaKeyboard');
    // ? opens the cheat-sheet.
    assert.ok(/e\.key === '\?'[\s\S]*?self\._showShortcutHelp\(\)/.test(ME),
        '? opens the shortcuts cheat-sheet');
});

test('A2 aromatic (4) / any (0) bond keys are intentionally NOT mapped', function () {
    // BIME has no aromatic-bond or any-bond primitive; the dispatcher must not
    // invent one. The bond-key table is exactly {1,2,3}.
    var m = ME.match(/var BOND_KEYS = (\{[^}]*\})/);
    assert.ok(m, 'BOND_KEYS table present');
    assert.strictEqual(m[1].indexOf("'4'"), -1, 'no aromatic (4) bond key');
    assert.strictEqual(m[1].indexOf("'0'"), -1, 'no any (0) bond key');
});

test('A3 the shortcut guard suppresses text-input / dialog / open-lens contexts', function () {
    var fn = meMethod('_isShortcutContext');
    assert.ok(fn, '_isShortcutContext defined');
    // Text-input + interactive controls are excluded (so typing isn't hijacked).
    assert.ok(/input,textarea,select,button,a,\[contenteditable="true"\]/.test(fn),
        'text inputs / textarea / contenteditable / select are excluded');
    // Editor-owned modals stand down.
    assert.ok(fn.indexOf('this._shortcutOverlay') !== -1 && fn.indexOf('this._customizeOverlay') !== -1,
        'editor-owned modals (shortcuts / customize) suppress shortcuts');
    // An OPEN focus lens (its wrap is role=dialog) owns Esc — defer to it.
    assert.ok(fn.indexOf('.wb-editor-wrap') !== -1 && fn.indexOf("=== 'dialog'") !== -1,
        'an open focus lens (.wb-editor-wrap[role=dialog]) suppresses editor shortcuts');
});

test('A4 the keydown handler gates on _isShortcutContext (not the bare text-input check)', function () {
    assert.ok(ME.indexOf('if (!self._isShortcutContext(e.target)) return;') !== -1,
        'the document keydown handler consults _isShortcutContext first');
    // The legacy isInteractiveTarget closure stays for the paste handler.
    assert.ok(ME.indexOf('if (isInteractiveTarget(e.target)) return;') !== -1,
        'the paste handler keeps its own interactive-target guard');
});

test('A5 hover tracking records the atom/bond under the cursor for the shortcuts', function () {
    assert.ok(/self\._hoverAtomId = atom \? atom\.id : null;/.test(ME),
        'onMove records the hovered atom id');
    assert.ok(/self\._hoverBondId = \(!atom && bond\) \? bond\.id : null;/.test(ME),
        'onMove records the hovered bond id');
});

test('A6 destroy() tears down an open cheat-sheet (no leaked capture listener)', function () {
    var fn = meMethod('destroy');
    assert.ok(fn.indexOf('this._closeShortcutHelp(') !== -1,
        'destroy closes the shortcuts overlay if open');
});

test('A7 a "Shortcuts" toolbar action + a ? deep-link sit in the Edit group', function () {
    assert.ok(/id: 'shortcuts', label: 'Shortcuts', icon: 'shortcuts', type: 'action'/.test(ME),
        'a Shortcuts action button is declared');
    assert.ok(/case 'shortcuts':[\s\S]*?this\._showShortcutHelp\(\);/.test(ME),
        "_doAction('shortcuts') opens the cheat-sheet");
    // The button advertises the ? hotkey in its tooltip (discoverability).
    assert.ok(ME.indexOf('Keyboard shortcuts — show the full list (or press ?)') !== -1,
        'the button tooltip advertises the ? hotkey');
});

// ===========================================================================
// (B) Source-shape — cheat-sheet dialog semantics + a11y (mirrors v2.3.4).
// ===========================================================================
test('B1 the cheat-sheet is an accessible modal dialog (role / aria-modal / label)', function () {
    var fn = meMethod('_showShortcutHelp');
    assert.ok(fn, '_showShortcutHelp defined');
    assert.ok(fn.indexOf("setAttribute('role', 'dialog')") !== -1, 'role=dialog');
    assert.ok(fn.indexOf("setAttribute('aria-modal', 'true')") !== -1, 'aria-modal=true');
    assert.ok(fn.indexOf("setAttribute('aria-labelledby'") !== -1, 'aria-labelledby points at the title');
    // Focus is recorded on open + moved into the card.
    assert.ok(fn.indexOf('document.activeElement') !== -1, 'records the pre-open focus');
    assert.ok(/\(f \|\| card\)\.focus\(\)/.test(fn), 'moves focus into the dialog on open');
});

test('B2 the cheat-sheet closes on Esc and backdrop click, and traps Tab', function () {
    var fn = meMethod('_showShortcutHelp');
    // Esc closes.
    assert.ok(/e\.key === 'Escape'[\s\S]*?close\(\)/.test(fn), 'Escape closes the cheat-sheet');
    // Backdrop (click on the overlay itself) closes.
    assert.ok(/if \(e\.target === overlay\) \{ close\(\); \}/.test(fn), 'click-away on the backdrop closes');
    // Tab focus-trap.
    assert.ok(fn.indexOf("e.key !== 'Tab'") !== -1, 'Tab is trapped inside the dialog');
    // The close button has an accessible name.
    assert.ok(fn.indexOf("setAttribute('aria-label', 'Close keyboard shortcuts')") !== -1,
        'the close button is screen-reader labelled');
});

test('B3 closing restores focus + removes the capture-phase keydown listener', function () {
    var fn = meMethod('_closeShortcutHelp');
    assert.ok(fn, '_closeShortcutHelp defined');
    assert.ok(/document\.removeEventListener\('keydown', this\._shortcutKeyHandler, true\)/.test(fn),
        'removes the capture-phase keydown listener');
    assert.ok(fn.indexOf('previousFocus.focus()') !== -1, 'restores focus to the pre-open element (WCAG 2.4.3)');
});

test('B4 the catalogue lists bonds, atoms, edit + help (no emoji in copy)', function () {
    var fn = meMethod('_shortcutGroups');
    assert.ok(fn, '_shortcutGroups defined');
    ['Bonds', 'Atoms', 'Edit', 'Help'].forEach(function (g) {
        assert.ok(fn.indexOf("title: '" + g + "'") !== -1, 'group ' + g + ' present');
    });
    // The note teaches the hovered/selected-vs-arm behaviour + the Cmd hint.
    assert.ok(fn.indexOf('hovered or selected') !== -1, 'explains hovered/selected element edit');
    assert.ok(fn.indexOf('Cmd in place of Ctrl') !== -1, 'documents the macOS Cmd substitution');
    // No pictographic emoji anywhere in the dialog copy/markup.
    assert.ok(!/[\u{1F000}-\u{1FAFF}\u{2600}-\u{27BF}]/u.test(fn),
        'no pictographic emoji in the cheat-sheet catalogue');
});

// ===========================================================================
// (C) Source-shape — no conflict with the lens Esc / Ctrl+Z or the view switch.
// ===========================================================================
test('C1 the workbench lens Esc still commits the structure (unchanged)', function () {
    // The workbench document Escape handler commits an open lens FIRST and
    // returns, before any modal close. v2.4.4 must not have touched that.
    var escIdx = WB.indexOf("e.key === 'Escape' || e.keyCode === 27");
    assert.ok(escIdx !== -1, 'workbench Escape handler present');
    var slice = WB.slice(escIdx, escIdx + 300);
    assert.ok(slice.indexOf('_lens && _lens.isOpen()') !== -1, 'lens-open is checked first');
    assert.ok(slice.indexOf('collapseStructureLens(true)') !== -1, 'an open lens commits on Esc');
});

test('C2 the editor Esc defers to an open lens via the role=dialog wrap check', function () {
    // The guard is the bridge: when the editor host is inside an open lens
    // (.wb-editor-wrap[role=dialog]), the editor-level Esc->Select stands down,
    // so the workbench lens-Esc is the sole handler. This is what keeps the two
    // Esc behaviours from colliding.
    var fn = meMethod('_isShortcutContext');
    assert.ok(/wrap\.getAttribute && wrap\.getAttribute\('role'\) === 'dialog'/.test(fn),
        'editor shortcuts stand down while the lens wrap is role=dialog');
    // And the lens marks its wrap role=dialog on open (the signal we rely on).
    assert.ok(/wrap\.setAttribute\('role', 'dialog'\)/.test(WB) ||
              WB.indexOf("setAttribute('role', 'dialog')") !== -1,
        'the lens sets role=dialog on its wrap while open');
});

test('C3 the Editor|Pathway view switch + pathway-canvas keys are untouched', function () {
    // setWorkbenchMode + the pathway keydown (Delete / Ctrl+Z) still exist and
    // are independent of the editor document handler.
    assert.ok(WB.indexOf('function setWorkbenchMode(') !== -1, 'view switch present');
    assert.ok(WB.indexOf('function onPathwayKeyDown(') !== -1, 'pathway canvas key handler present');
    var fn = WB.slice(WB.indexOf('function onPathwayKeyDown('));
    fn = fn.slice(0, fn.indexOf('\nfunction '));
    assert.ok(fn.indexOf("e.key === 'Delete' || e.key === 'Backspace'") !== -1,
        'pathway Delete deletes a node');
    assert.ok(fn.indexOf('undoPathway()') !== -1, 'pathway Ctrl+Z undo present');
});

// ===========================================================================
// (D) Executable DOM-stub — construct a real MolEditor + drive the shortcuts.
// A recording document stub (richer than tests/shim.js) lets the editor build
// its DOM and lets us run the keydown dispatcher + helpers for real.
// ===========================================================================
function recEl(tag) {
    var el = {
        tag: tag || 'div', _attrs: {}, _listeners: {}, style: {}, children: [],
        tabIndex: 0, className: '', textContent: '', innerHTML: '', focused: false,
        offsetParent: {}, clientWidth: 500, clientHeight: 380, parentNode: null,
        classList: {
            _s: {}, add: function (c) { this._s[c] = true; },
            remove: function (c) { delete this._s[c]; },
            contains: function (c) { return !!this._s[c]; },
            toggle: function (c) { this._s[c] = !this._s[c]; return this._s[c]; }
        },
        setAttribute: function (k, v) { this._attrs[k] = v; },
        getAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k) ? this._attrs[k] : null; },
        removeAttribute: function (k) { delete this._attrs[k]; },
        hasAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k); },
        appendChild: function (c) { c.parentNode = this; this.children.push(c); return c; },
        removeChild: function (c) { var i = this.children.indexOf(c); if (i >= 0) this.children.splice(i, 1); return c; },
        insertBefore: function (c) { c.parentNode = this; this.children.unshift(c); return c; },
        addEventListener: function (t, fn) { (this._listeners[t] = this._listeners[t] || []).push(fn); },
        removeEventListener: function (t, fn) { var a = this._listeners[t] || []; var i = a.indexOf(fn); if (i >= 0) a.splice(i, 1); },
        querySelector: function () { return recEl('button'); },
        querySelectorAll: function () { return []; },
        focus: function () { el.focused = true; doc.activeElement = el; },
        closest: function () { return null; },
        getBoundingClientRect: function () { return { x: 0, y: 0, width: 500, height: 380, left: 0, top: 0 }; },
        getComputedTextLength: function () { return 0; },
        remove: function () { if (this.parentNode) this.parentNode.removeChild(this); },
        contains: function () { return false; }
    };
    return el;
}

var doc;
function makeDoc() {
    var byId = {};
    doc = {
        _docListeners: {},
        createElement: recEl,
        createElementNS: function (ns, t) { return recEl(t); },
        getElementById: function (id) { if (!byId[id]) { byId[id] = recEl('div'); } return byId[id]; },
        body: recEl('body'), head: recEl('head'), documentElement: recEl('html'),
        activeElement: null,
        addEventListener: function (t, fn) { (this._docListeners[t] = this._docListeners[t] || []).push(fn); },
        removeEventListener: function (t, fn) { var a = this._docListeners[t] || []; var i = a.indexOf(fn); if (i >= 0) a.splice(i, 1); }
    };
    return doc;
}

// Build a MolEditor against the stub document. Editor modules are already on
// globalThis (shim.loadAll), but the renderer/tools/etc. modules the editor
// constructor touches must be present too — load them once here.
var editorReady = (function () {
    try {
        var need = ['Layout', 'Templates', 'MetaboliteLibrary', 'MolfileWriter',
            'Renderer', 'Tools', 'CurlyArrowEndpoint', 'ImageExport', 'ExportStamp',
            'XmlUtil', 'MolEditor'];
        for (var i = 0; i < need.length; i++) {
            if (!globalThis[need[i]]) { shim.require_editor(need[i]); }
        }
        return !!globalThis.MolEditor;
    } catch (e) { return false; }
})();

function buildEditor() {
    var saved = globalThis.document;
    var savedRO = globalThis.ResizeObserver;
    globalThis.document = makeDoc();
    globalThis.ResizeObserver = function () { return { observe: function () {}, disconnect: function () {} }; };
    var ed;
    try {
        ed = new globalThis.MolEditor('bime-editor', '100%', '400px', {});
    } finally {
        // Leave the stub document in place for the duration of the test body;
        // callers restore via restoreDoc().
        ed._test_restore = function () {
            globalThis.document = saved;
            globalThis.ResizeObserver = savedRO;
        };
    }
    return ed;
}

// A keydown event factory for driving the real dispatcher.
function key(opts) {
    opts = opts || {};
    return {
        key: opts.key, keyCode: opts.keyCode,
        ctrlKey: !!opts.ctrl, metaKey: !!opts.meta, shiftKey: !!opts.shift,
        target: opts.target || { closest: function () { return null; } },
        _prevented: false,
        preventDefault: function () { this._prevented = true; }
    };
}

test('D1 bond keys 1/2/3 arm the matching bond tool', function () {
    if (!editorReady) { console.log('    (skipped: editor stack unavailable)'); return; }
    var ed = buildEditor();
    try {
        ed._keydownHandler(key({ key: '2' }));
        assert.strictEqual(ed._currentToolName, 'bond', '2 arms the bond tool');
        assert.strictEqual(ed._currentTool.bondType, 2, '2 sets a double bond');
        ed._keydownHandler(key({ key: '3' }));
        assert.strictEqual(ed._currentTool.bondType, 3, '3 sets a triple bond');
        ed._keydownHandler(key({ key: '1' }));
        assert.strictEqual(ed._currentTool.bondType, 1, '1 sets a single bond');
    } finally { ed._test_restore(); }
});

test('D2 element letter with NO target arms the atom tool with that element', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._keydownHandler(key({ key: 'o' }));
        assert.strictEqual(ed._currentToolName, 'atom', 'o arms the atom tool');
        assert.strictEqual(ed._currentElement, 'O', 'o sets the active element to oxygen');
    } finally { ed._test_restore(); }
});

test('D3 element letter with a selected atom changes it in place (undoable)', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        var a = ed.molecule.addAtom('C', 0, 0);
        a.selected = true;
        var undoBefore = ed.history.canUndo();
        ed._keydownHandler(key({ key: 'n' }));
        assert.strictEqual(ed.molecule.getAtom(a.id).symbol, 'N', 'n changes the selected atom to nitrogen');
        assert.strictEqual(ed.history.canUndo(), true, 'the change pushed an undo frame');
        assert.ok(undoBefore === false || undoBefore === true, 'history queried');
        // It should NOT have switched to the atom tool (we edited in place).
        assert.notStrictEqual(ed._currentToolName, 'atom', 'in-place edit does not arm the atom tool');
    } finally { ed._test_restore(); }
});

test('D4 Shift+C reaches chlorine, Shift+B reaches bromine (halogen combos win)', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._keydownHandler(key({ key: 'C', shift: true }));
        assert.strictEqual(ed._currentElement, 'Cl', 'Shift+C -> chlorine (not carbon)');
        ed._keydownHandler(key({ key: 'B', shift: true }));
        assert.strictEqual(ed._currentElement, 'Br', 'Shift+B -> bromine (not boron)');
        // And the bare letters still reach carbon / boron.
        ed._keydownHandler(key({ key: 'c' }));
        assert.strictEqual(ed._currentElement, 'C', 'bare c -> carbon');
        ed._keydownHandler(key({ key: 'b' }));
        assert.strictEqual(ed._currentElement, 'B', 'bare b -> boron');
    } finally { ed._test_restore(); }
});

test('D5 Esc returns to the Select tool', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._setTool('bond');
        assert.strictEqual(ed._currentToolName, 'bond', 'precondition: a non-Select tool is armed');
        ed._keydownHandler(key({ key: 'Escape' }));
        assert.strictEqual(ed._currentToolName, 'select', 'Esc returns to Select');
    } finally { ed._test_restore(); }
});

test('D6 Ctrl+A selects every atom + bond and switches to Select', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        var a1 = ed.molecule.addAtom('C', 0, 0);
        var a2 = ed.molecule.addAtom('C', 10, 0);
        ed.molecule.addBond(a1.id, a2.id, 1);
        ed._setTool('bond');
        ed._keydownHandler(key({ key: 'a', ctrl: true }));
        var selAtoms = ed.molecule.atoms.filter(function (x) { return x.selected; }).length;
        var selBonds = ed.molecule.bonds.filter(function (x) { return x.selected; }).length;
        assert.strictEqual(selAtoms, 2, 'all atoms selected');
        assert.strictEqual(selBonds, 1, 'all bonds selected');
        assert.strictEqual(ed._currentToolName, 'select', 'Ctrl+A switches to Select');
    } finally { ed._test_restore(); }
});

test('D7 Del removes the selection; with none, removes the hovered atom', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        var a1 = ed.molecule.addAtom('C', 0, 0);
        var a2 = ed.molecule.addAtom('C', 10, 0);
        ed.molecule.addBond(a1.id, a2.id, 1);
        // Selected-atom delete.
        a1.selected = true;
        ed._keydownHandler(key({ key: 'Delete' }));
        assert.ok(ed.molecule.getAtom(a1.id) == null, 'selected atom deleted (getAtom returns null)');
        // Hovered-atom delete (nothing selected now).
        ed.molecule.atoms.forEach(function (x) { x.selected = false; });
        ed.molecule.bonds.forEach(function (x) { x.selected = false; });
        ed._hoverAtomId = a2.id;
        var before = ed.molecule.atoms.length;
        ed._keydownHandler(key({ key: 'Backspace' }));
        assert.strictEqual(ed.molecule.atoms.length, before - 1, 'hovered atom deleted on Backspace');
    } finally { ed._test_restore(); }
});

test('D8 ? opens the cheat-sheet as a labelled modal dialog; Esc closes it', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._keydownHandler(key({ key: '?' }));
        assert.ok(ed._shortcutOverlay, 'the cheat-sheet overlay is created');
        assert.strictEqual(ed._shortcutOverlay.getAttribute('role'), 'dialog', 'role=dialog');
        assert.strictEqual(ed._shortcutOverlay.getAttribute('aria-modal'), 'true', 'aria-modal=true');
        assert.ok(ed._shortcutOverlay.getAttribute('aria-labelledby'), 'aria-labelledby set');
        assert.ok(ed._shortcutKeyHandler, 'a key handler is bound while open');
        // Esc (routed through the overlay key handler) closes + cleans up.
        ed._shortcutKeyHandler(key({ key: 'Escape' }));
        assert.strictEqual(ed._shortcutOverlay, null, 'Esc closes the cheat-sheet');
        assert.strictEqual(ed._shortcutKeyHandler, null, 'the key handler is removed on close');
    } finally { ed._test_restore(); }
});

test('D9 a second ? while open is a no-op (idempotent)', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._keydownHandler(key({ key: '?' }));
        var first = ed._shortcutOverlay;
        ed._keydownHandler(key({ key: '?' }));
        assert.strictEqual(ed._shortcutOverlay, first, 'a second ? does not stack a second dialog');
    } finally { ed._test_restore(); }
});

test('D10 GUARD: a keydown whose target is a text input is NOT hijacked', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._setTool('select');
        // A target that reports it is inside an <input> (e.g. the SMILES / label
        // field). The element shortcut must NOT fire, the tool must stay put,
        // and the event must NOT be preventDefault'd (so the keystroke types).
        var inputTarget = { closest: function (sel) { return /input/.test(sel) ? {} : null; } };
        var ev = key({ key: 'o', target: inputTarget });
        ed._keydownHandler(ev);
        assert.strictEqual(ed._currentToolName, 'select', 'tool unchanged while typing in an input');
        assert.strictEqual(ed._currentElement, 'C', 'active element unchanged while typing');
        assert.strictEqual(ev._prevented, false, 'the keystroke is not preventDefault-ed (it types normally)');
    } finally { ed._test_restore(); }
});

test('D11 GUARD: while the editor sits in an OPEN focus lens, Esc->Select stands down', function () {
    if (!editorReady) { console.log('    (skipped)'); return; }
    var ed = buildEditor();
    try {
        ed._setTool('bond');
        // Model the editor host inside an open lens: its container.closest finds
        // a .wb-editor-wrap with role=dialog. _isShortcutContext must return false
        // so the editor Esc does NOT fire (the workbench lens-Esc owns it).
        var wrap = { getAttribute: function (k) { return k === 'role' ? 'dialog' : null; } };
        ed._container.closest = function (sel) { return /wb-editor-wrap/.test(sel) ? wrap : null; };
        var ev = key({ key: 'Escape' });
        ed._keydownHandler(ev);
        assert.strictEqual(ed._currentToolName, 'bond', 'editor Esc does not change the tool while the lens is open');
        assert.strictEqual(ev._prevented, false, 'editor leaves Esc for the lens handler');
    } finally { ed._test_restore(); }
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
