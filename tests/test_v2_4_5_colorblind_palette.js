/**
 * tests/test_v2_4_5_colorblind_palette.js — colour-blind-safe (CVD) palette
 * + opt-in toggle (v2.4.5).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.4.5 adds the COLOUR layer for users who SEE colour but with a
 * colour-vision deficiency (CVD): deuteranopia / protanopia (red-green) and
 * tritanopia. It is purely ADDITIVE on top of the v2.3.7 non-colour cues — the
 * node-type letter badges (M / R / Co / Rxn), the reaction-arrow shapes, and
 * the atom-trace border PATTERNS (solid / dashed / dotted / double) all STAY.
 * This release only re-tints the semantic colours, behind an opt-in toggle,
 * using an established colour-blind-safe palette so adjacent
 * semantics differ in hue AND lightness under simulated CVD.
 *
 * This suite is two halves:
 *
 *   (A–C) Source-shape over css/workbench.css: the `[data-palette="cvd"]`
 *         override EXISTS, remaps the four atom-trace moiety states + the
 *         node-type colours + the quality badges to DISTINCT CVD-safe values,
 *         is opt-in (the default palette is untouched), and composes with dark
 *         mode (it does not key on / fight `[data-theme]`).
 *   (D)   The v2.3.7 non-colour cues are STILL present (additive guarantee):
 *         the border PATTERNS and the M/R/Co/Rxn badge are unchanged.
 *   (E)   Executable DOM-stub: the real js/nav.js is evaluated in a vm against
 *         a recording document / window / localStorage stub, then the
 *         `.palette-toggle` is clicked to assert it SETS, PERSISTS and READS
 *         `data-palette` (localStorage key 'bime-palette'), reflects
 *         aria-pressed, and that a stored preference is applied on load even
 *         before DOMContentLoaded — and survives alongside a dark-mode
 *         data-theme attribute (composition).
 *
 * Plain Node, no real DOM.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var vm = require('vm');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Colour-blind-safe palette + toggle (v2.4.5)');
var test = runner.test;
console.log('Colour-blind-safe palette + toggle (v2.4.5)');

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

var CSS = readRepoFile('css/workbench.css');
var STYLE = readRepoFile('css/style.css');
var NAV = readRepoFile('js/nav.js');
var HTML = readRepoFile('workbench.html');

// The colour-blind-safe palette (hex, upper-case as in the CSS). Used
// to assert the remapped semantics land on this established scheme.
var CVD_SAFE = {
    bluishGreen: '#009E73',
    orange: '#E69F00',
    skyBlue: '#56B4E9',
    vermillion: '#D55E00',
    blue: '#0072B2',
    purple: '#CC79A7'
};

// Extract the value of a single declaration (`prop`) from the FIRST CSS rule
// whose selector text contains `selector`. Returns '' if not found. Good
// enough for these flat, hand-written rules (no nested braces inside).
function declOf(css, selector, prop) {
    var at = css.indexOf(selector);
    if (at === -1) { return ''; }
    var open = css.indexOf('{', at);
    var close = css.indexOf('}', open);
    if (open === -1 || close === -1) { return ''; }
    var body = css.slice(open + 1, close);
    var re = new RegExp('(?:^|;)\\s*' + prop.replace(/[-]/g, '\\-') + '\\s*:\\s*([^;]+)');
    var m = body.match(re);
    return m ? m[1].trim() : '';
}

// ===========================================================================
// A. The [data-palette="cvd"] override exists and is opt-in.
// ===========================================================================
test('A1 a [data-palette="cvd"] override block exists in workbench.css', function () {
    assert.ok(CSS.indexOf('[data-palette="cvd"]') !== -1,
        'workbench.css defines a [data-palette="cvd"] palette override');
    // The CVD design tokens are named once on the attribute root.
    var tokenBlock = declOf(CSS, '[data-palette="cvd"] {', '--cvd-vermillion');
    assert.ok(tokenBlock.toUpperCase().indexOf(CVD_SAFE.vermillion) !== -1,
        'the CVD token block defines --cvd-vermillion as the CVD-safe vermillion');
});

test('A2 the palette is OPT-IN — the default semantics are not redefined', function () {
    // Every CVD override must be SCOPED under [data-palette="cvd"]; there must
    // be no bare re-paint of the semantic rules at the root. Spot-check the
    // four trace-state selectors: the only occurrences after the default rules
    // must be the data-palette-scoped ones.
    var states = ['is-moiety-intact', 'is-moiety-fragmented', 'is-moiety-partial-loss', 'is-moiety-absent'];
    states.forEach(function (st) {
        var re = new RegExp('([^\\n]*)\\.wb-atom-trace-cell\\.' + st + '\\s*\\{', 'g');
        var m;
        var sawDefault = false;
        var sawCvd = false;
        while ((m = re.exec(CSS)) !== null) {
            if (m[1].indexOf('[data-palette="cvd"]') !== -1) { sawCvd = true; }
            else { sawDefault = true; }
        }
        assert.ok(sawDefault, st + ' keeps its DEFAULT rule (unchanged)');
        assert.ok(sawCvd, st + ' gains a [data-palette="cvd"] override');
    });
});

// ===========================================================================
// B. The four atom-trace moiety states remap to DISTINCT CVD-safe colours.
//    (The worst case: default green / amber / orange / red — the whole
//    red-green axis — must become four hues separable under deuteranopia.)
// ===========================================================================
test('B1 each trace state gets a CVD-safe border-color override', function () {
    function cvdBorder(stateClass) {
        return declOf(CSS, '[data-palette="cvd"] .wb-atom-trace-cell.' + stateClass + ' {', 'border-color');
    }
    var intact = cvdBorder('is-moiety-intact');
    var fragmented = cvdBorder('is-moiety-fragmented');
    var partial = cvdBorder('is-moiety-partial-loss');
    var absent = cvdBorder('is-moiety-absent');
    // They resolve (directly or via a --cvd-* token) to the CVD-safe palette colours.
    assert.ok(/--cvd-bluish-green|#009E73/i.test(intact), 'intact => bluish-green');
    assert.ok(/--cvd-orange|#E69F00/i.test(fragmented), 'fragmented => orange');
    assert.ok(/--cvd-sky-blue|#56B4E9/i.test(partial), 'partial-loss => sky-blue');
    assert.ok(/--cvd-vermillion|#D55E00/i.test(absent), 'absent => vermillion');
});

test('B2 the four trace-state CVD colours are pairwise DISTINCT', function () {
    // Resolve each state's border-color token to its concrete hex via the
    // --cvd-* definitions, then assert all four differ (separable under CVD).
    function resolve(stateClass) {
        var raw = declOf(CSS, '[data-palette="cvd"] .wb-atom-trace-cell.' + stateClass + ' {', 'border-color');
        var tok = raw.match(/--([a-z-]+)/);
        if (tok) {
            var hex = declOf(CSS, '[data-palette="cvd"] {', '--' + tok[1]);
            return hex.toUpperCase();
        }
        return raw.toUpperCase();
    }
    var hexes = ['is-moiety-intact', 'is-moiety-fragmented', 'is-moiety-partial-loss', 'is-moiety-absent']
        .map(resolve);
    hexes.forEach(function (h) {
        assert.ok(/^#[0-9A-F]{6}$/.test(h), 'resolved to a concrete hex: ' + h);
    });
    var distinct = {};
    hexes.forEach(function (h) { distinct[h] = true; });
    assert.strictEqual(Object.keys(distinct).length, 4,
        'all four atom-trace states map to DISTINCT CVD-safe hues, got: ' + hexes.join(' '));
    // And they are exactly the four chosen CVD-safe colours.
    assert.deepStrictEqual(hexes,
        [CVD_SAFE.bluishGreen, CVD_SAFE.orange, CVD_SAFE.skyBlue, CVD_SAFE.vermillion],
        'the four states are bluish-green / orange / sky-blue / vermillion');
});

test('B3 the single-atom trace + moiety-source + cofactor-origin marks are re-tinted distinctly', function () {
    var traced = declOf(CSS, '[data-palette="cvd"] .bime-trace-atom.is-traced .bime-trace-atom-mark {', 'stroke');
    var moiety = declOf(CSS, '[data-palette="cvd"] .bime-trace-atom.is-moiety-source .bime-trace-atom-mark {', 'stroke');
    var cofactor = declOf(CSS, '[data-palette="cvd"] .bime-trace-atom.is-cofactor-origin .bime-trace-atom-mark {', 'stroke');
    assert.ok(/--cvd-orange/.test(traced), 'single-trace highlight => CVD orange');
    assert.ok(/--cvd-bluish-green/.test(moiety), 'moiety-source ring => CVD bluish-green');
    assert.ok(/--cvd-purple/.test(cofactor), 'cofactor-origin ring => CVD reddish-purple');
});

// ===========================================================================
// C. Node-type colours + quality badges remap to CVD-safe values.
// ===========================================================================
test('C1 each pathway node TYPE gets a CVD-safe stroke override', function () {
    var residue = declOf(CSS, '[data-palette="cvd"] .wb-pathway-node.is-residue rect {', 'stroke');
    var cofactor = declOf(CSS, '[data-palette="cvd"] .wb-pathway-node.is-cofactor rect {', 'stroke');
    var step = declOf(CSS, '[data-palette="cvd"] .wb-pathway-step rect {', 'stroke');
    var note = declOf(CSS, '[data-palette="cvd"] .wb-pathway-note rect {', 'stroke');
    assert.ok(/--cvd-bluish-green/.test(residue), 'residue stroke => bluish-green');
    assert.ok(/--cvd-orange/.test(cofactor), 'cofactor stroke => orange');
    assert.ok(/--cvd-blue/.test(step), 'reaction step stroke => blue');
    assert.ok(/--cvd-purple/.test(note), 'note stroke => reddish-purple');
});

test('C2 the node-type CVD strokes are pairwise distinct (separable under CVD)', function () {
    function resolveTok(tokName) {
        return declOf(CSS, '[data-palette="cvd"] {', tokName).toUpperCase();
    }
    var strokes = {
        residue: resolveTok('--cvd-bluish-green'),
        cofactor: resolveTok('--cvd-orange'),
        step: resolveTok('--cvd-blue'),
        note: resolveTok('--cvd-purple')
    };
    var vals = Object.keys(strokes).map(function (k) { return strokes[k]; });
    var distinct = {};
    vals.forEach(function (v) {
        assert.ok(/^#[0-9A-F]{6}$/.test(v), 'token resolves to a hex: ' + v);
        distinct[v] = true;
    });
    assert.strictEqual(Object.keys(distinct).length, 4,
        'residue / cofactor / step / note strokes are all distinct under CVD');
});

test('C3 the quality badge good / warn / bad remap to CVD-safe values', function () {
    var good = declOf(CSS, '[data-palette="cvd"] .wb-quality-badge.is-good {', 'border-color');
    var warn = declOf(CSS, '[data-palette="cvd"] .wb-quality-badge.is-warn {', 'border-color');
    var bad = declOf(CSS, '[data-palette="cvd"] .wb-quality-badge.is-bad {', 'border-color');
    assert.ok(/--cvd-bluish-green/.test(good), 'good => bluish-green');
    assert.ok(/--cvd-orange/.test(warn), 'warn => orange');
    assert.ok(/--cvd-vermillion/.test(bad), 'bad => vermillion');
    // good vs bad were green vs red (red-green confusable); now bluish-green vs
    // vermillion — distinct under deuteranopia.
    assert.notStrictEqual(good, bad, 'good and bad no longer collide on the red-green axis');
});

// ===========================================================================
// D. The v2.3.7 NON-COLOUR cues are STILL present (additive guarantee).
// ===========================================================================
test('D1 atom-trace border PATTERNS (solid/dashed/dotted/double) are untouched', function () {
    function defaultBorderStyle(stateClass) {
        // The DEFAULT (non-data-palette) rule still owns the pattern.
        var re = new RegExp('(^|\\n)\\.wb-atom-trace-cell\\.' + stateClass + '\\s*\\{');
        var m = re.exec(CSS);
        assert.ok(m, stateClass + ' default rule present');
        var open = CSS.indexOf('{', m.index);
        var block = CSS.slice(open, CSS.indexOf('}', open));
        var sm = block.match(/border-style:\s*([a-z]+)/);
        assert.ok(sm, stateClass + ' default rule declares a border-style');
        return sm[1];
    }
    assert.strictEqual(defaultBorderStyle('is-moiety-intact'), 'solid', 'intact still solid');
    assert.strictEqual(defaultBorderStyle('is-moiety-fragmented'), 'dashed', 'fragmented still dashed');
    assert.strictEqual(defaultBorderStyle('is-moiety-partial-loss'), 'dotted', 'partial still dotted');
    assert.strictEqual(defaultBorderStyle('is-moiety-absent'), 'double', 'absent still double');
    // The CVD overrides must NOT set border-style (so the pattern is preserved).
    var states = ['is-moiety-intact', 'is-moiety-fragmented', 'is-moiety-partial-loss', 'is-moiety-absent'];
    states.forEach(function (st) {
        var cvd = declOf(CSS, '[data-palette="cvd"] .wb-atom-trace-cell.' + st + ' {', 'border-style');
        assert.strictEqual(cvd, '', st + ' CVD override leaves border-style (the pattern cue) alone');
    });
});

test('D2 the node-type letter badge + the moiety status-tag chip are still styled', function () {
    // The v2.3.7 greyscale badge chip is untouched.
    assert.ok(CSS.indexOf('.wb-pathway-node .wb-pathway-node-badge rect') !== -1,
        'pathway node-type letter badge (M/R/Co/Rxn) chip is still present');
    assert.ok(CSS.indexOf('.wb-atom-trace-moiety-tag') !== -1,
        'atom-trace status-tag (glyph + word) chip is still present');
    // The CVD block must not redefine the badge chip (it stays slate-on-pale —
    // a colour-independent cue), so the cue is not weakened.
    assert.strictEqual(CSS.indexOf('[data-palette="cvd"] .wb-pathway-node-badge'), -1,
        'the CVD palette does not override the (colour-independent) node-type badge');
});

// ===========================================================================
// E. The toggle: a labelled, accessible control that SETS / PERSISTS / READS
//    data-palette via localStorage — driven against the REAL js/nav.js.
// ===========================================================================

// A minimal recording DOM/window for evaluating nav.js. documentElement is the
// <html> the palette attribute lives on; querySelector returns our registered
// stub elements; localStorage is an in-memory map. matchMedia is needed by the
// (unrelated) dark-mode branch so nav.js runs end-to-end.
function makeElementStub(tagClass) {
    return {
        _attrs: {},
        _listeners: {},
        className: tagClass || '',
        classList: {
            _set: {},
            contains: function (c) { return !!this._set[c]; },
            add: function (c) { this._set[c] = true; },
            remove: function (c) { delete this._set[c]; },
            toggle: function (c) { this._set[c] = !this._set[c]; return this._set[c]; }
        },
        innerHTML: '',
        setAttribute: function (k, v) { this._attrs[k] = String(v); },
        getAttribute: function (k) {
            return Object.prototype.hasOwnProperty.call(this._attrs, k) ? this._attrs[k] : null;
        },
        removeAttribute: function (k) { delete this._attrs[k]; },
        hasAttribute: function (k) { return Object.prototype.hasOwnProperty.call(this._attrs, k); },
        addEventListener: function (type, fn) {
            (this._listeners[type] = this._listeners[type] || []).push(fn);
        },
        querySelectorAll: function () { return []; },
        contains: function () { return false; },
        focus: function () {},
        dispatch: function (type, evt) {
            (this._listeners[type] || []).forEach(function (fn) { fn(evt || {}); });
        }
    };
}

// Build a sandbox, register the named query-selector stubs, evaluate nav.js in
// it (which runs applyStoredPalette immediately + registers DOMContentLoaded),
// then fire DOMContentLoaded. Returns the live stubs for assertions.
function bootNav(opts) {
    opts = opts || {};
    var store = opts.store || {};
    var documentElement = makeElementStub('html');
    if (opts.preAttrs) {
        Object.keys(opts.preAttrs).forEach(function (k) {
            documentElement.setAttribute(k, opts.preAttrs[k]);
        });
    }
    var elements = {
        '.palette-toggle': makeElementStub('palette-toggle'),
        '.theme-toggle': makeElementStub('theme-toggle')
        // .nav-hamburger / .nav-links intentionally absent => those branches skip.
    };
    var docListeners = {};
    var documentStub = {
        documentElement: documentElement,
        querySelector: function (sel) {
            return Object.prototype.hasOwnProperty.call(elements, sel) ? elements[sel] : null;
        },
        addEventListener: function (type, fn) {
            (docListeners[type] = docListeners[type] || []).push(fn);
        }
    };
    var sandbox = {
        document: documentStub,
        window: {
            matchMedia: function () { return { matches: false }; }
        },
        localStorage: {
            getItem: function (k) {
                return Object.prototype.hasOwnProperty.call(store, k) ? store[k] : null;
            },
            setItem: function (k, v) { store[k] = String(v); },
            removeItem: function (k) { delete store[k]; }
        },
        console: { log: function () {}, warn: function () {}, error: function () {} }
    };
    vm.createContext(sandbox);
    vm.runInContext(NAV, sandbox, { filename: 'nav.js' });
    // Fire DOMContentLoaded so the toggle wiring runs.
    (docListeners['DOMContentLoaded'] || []).forEach(function (fn) { fn({}); });
    return {
        store: store,
        documentElement: documentElement,
        paletteBtn: elements['.palette-toggle']
    };
}

test('E1 the toggle button exists in the nav with an accessible label + aria-pressed', function () {
    // The control is markup-present, keyboard-reachable (a native <button>),
    // and labelled — next to the dark-mode toggle.
    assert.ok(/<button class="palette-toggle"[^>]*aria-pressed="false"[^>]*aria-label="[^"]+"/.test(HTML),
        'workbench.html has a <button class="palette-toggle"> with aria-pressed + aria-label');
    var navIdx = HTML.indexOf('class="palette-toggle"');
    var themeIdx = HTML.indexOf('class="theme-toggle"');
    assert.ok(navIdx !== -1 && themeIdx !== -1, 'both toggles present in the nav');
    // Styled like the theme toggle (shared bordered icon-button look) + has a
    // visible pressed state for the binary on/off.
    assert.ok(STYLE.indexOf('.palette-toggle') !== -1, 'css/style.css styles .palette-toggle');
    assert.ok(STYLE.indexOf('.palette-toggle[aria-pressed="true"]') !== -1,
        'an aria-pressed=true visual (active) state is styled (not colour-alone)');
});

test('E2 clicking the toggle SETS data-palette="cvd" and PERSISTS it to localStorage', function () {
    var ctx = bootNav();
    // Fresh load, nothing stored => default palette (no attribute).
    assert.strictEqual(ctx.documentElement.getAttribute('data-palette'), null,
        'default load leaves data-palette unset (opt-in)');
    assert.strictEqual(ctx.paletteBtn.getAttribute('aria-pressed'), 'false',
        'toggle starts unpressed');
    // Click ON.
    ctx.paletteBtn.dispatch('click');
    assert.strictEqual(ctx.documentElement.getAttribute('data-palette'), 'cvd',
        'click sets data-palette="cvd" on <html>');
    assert.strictEqual(ctx.store['bime-palette'], 'cvd',
        'the choice is persisted under localStorage["bime-palette"]');
    assert.strictEqual(ctx.paletteBtn.getAttribute('aria-pressed'), 'true',
        'aria-pressed reflects the ON state (WCAG 4.1.2)');
});

test('E3 clicking again toggles OFF and persists the default', function () {
    var ctx = bootNav();
    ctx.paletteBtn.dispatch('click');   // on
    ctx.paletteBtn.dispatch('click');   // off
    assert.strictEqual(ctx.documentElement.getAttribute('data-palette'), null,
        'second click removes the attribute (back to default)');
    assert.strictEqual(ctx.store['bime-palette'], 'default',
        'the OFF choice is persisted (so it overrides any earlier cvd)');
    assert.strictEqual(ctx.paletteBtn.getAttribute('aria-pressed'), 'false',
        'aria-pressed reflects the OFF state');
});

test('E4 a stored "cvd" preference is READ + applied on load, before DOMContentLoaded', function () {
    // applyStoredPalette() runs at script-evaluation time (top of nav.js), so the
    // attribute must be set the moment nav.js is evaluated — we assert it is
    // present immediately after evaluation, independent of the DOMContentLoaded
    // wiring, so a returning CVD user gets no flash of the default palette.
    var store = { 'bime-palette': 'cvd' };
    var documentElement = makeElementStub('html');
    var docListeners = {};
    var sandbox = {
        document: {
            documentElement: documentElement,
            querySelector: function () { return null; },
            addEventListener: function (t, fn) { (docListeners[t] = docListeners[t] || []).push(fn); }
        },
        window: { matchMedia: function () { return { matches: false }; } },
        localStorage: {
            getItem: function (k) { return Object.prototype.hasOwnProperty.call(store, k) ? store[k] : null; },
            setItem: function (k, v) { store[k] = String(v); },
            removeItem: function (k) { delete store[k]; }
        },
        console: { log: function () {} }
    };
    vm.createContext(sandbox);
    vm.runInContext(NAV, sandbox, { filename: 'nav.js' });
    assert.strictEqual(documentElement.getAttribute('data-palette'), 'cvd',
        'a stored cvd preference is applied at nav.js evaluation (pre-DOMContentLoaded)');
    // And the synced button (after DOMContentLoaded) would read pressed — but
    // here there is no button; the point is the early apply does not throw and
    // sets the attribute.
});

test('E5 the CVD palette COMPOSES with dark mode (orthogonal attributes)', function () {
    // Boot with dark mode already on (data-theme="dark") AND a stored cvd
    // preference: both attributes must end up on <html>, neither clobbering the
    // other. Then toggle CVD off: data-theme must be untouched.
    var ctx = bootNav({ store: { 'bime-palette': 'cvd' }, preAttrs: { 'data-theme': 'dark' } });
    assert.strictEqual(ctx.documentElement.getAttribute('data-palette'), 'cvd',
        'stored cvd applied');
    assert.strictEqual(ctx.documentElement.getAttribute('data-theme'), 'dark',
        'dark mode preserved alongside the CVD palette');
    assert.strictEqual(ctx.paletteBtn.getAttribute('aria-pressed'), 'true',
        'the toggle button syncs to pressed for the stored cvd preference on load');
    // Toggle CVD off — dark mode must remain.
    ctx.paletteBtn.dispatch('click');
    assert.strictEqual(ctx.documentElement.getAttribute('data-palette'), null,
        'CVD removed');
    assert.strictEqual(ctx.documentElement.getAttribute('data-theme'), 'dark',
        'dark mode is unaffected by the CVD toggle (composition holds)');
    // The CSS proves the two layers do not fight: the CVD overrides are NOT
    // scoped to a theme, so they apply in light OR dark.
    assert.strictEqual(CSS.indexOf('[data-theme="dark"][data-palette="cvd"]'), -1,
        'the CVD rules are theme-independent (no [data-theme] coupling needed)');
    assert.strictEqual(CSS.indexOf('[data-palette="cvd"][data-theme'), -1,
        'the CVD rules do not require a specific theme to apply');
});

test('E6 a localStorage failure on read/write does not throw (private browsing)', function () {
    // Mirror the dark-mode toggle resilience: a throwing localStorage must not
    // break page init or the toggle click.
    var documentElement = makeElementStub('html');
    var paletteBtn = makeElementStub('palette-toggle');
    var docListeners = {};
    var sandbox = {
        document: {
            documentElement: documentElement,
            querySelector: function (sel) { return sel === '.palette-toggle' ? paletteBtn : null; },
            addEventListener: function (t, fn) { (docListeners[t] = docListeners[t] || []).push(fn); }
        },
        window: { matchMedia: function () { return { matches: false }; } },
        localStorage: {
            getItem: function () { throw new Error('blocked'); },
            setItem: function () { throw new Error('blocked'); },
            removeItem: function () { throw new Error('blocked'); }
        },
        console: { log: function () {} }
    };
    vm.createContext(sandbox);
    assert.doesNotThrow(function () {
        vm.runInContext(NAV, sandbox, { filename: 'nav.js' });
        (docListeners['DOMContentLoaded'] || []).forEach(function (fn) { fn({}); });
    }, 'nav.js survives a throwing localStorage at load');
    assert.doesNotThrow(function () {
        paletteBtn.dispatch('click');
    }, 'the toggle click survives a throwing localStorage');
    // The in-DOM attribute still flips even if persistence failed.
    assert.strictEqual(documentElement.getAttribute('data-palette'), 'cvd',
        'the attribute still applies for the session even when persistence is blocked');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
