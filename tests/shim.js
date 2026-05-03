/**
 * tests/shim.js — Minimal Node.js compatibility shim for BIME editor modules
 *
 * Editor modules use a `(function(global) { ... })(window)` IIFE pattern that
 * attaches to `window`. This shim aliases `window = globalThis` and provides
 * just enough of `document`/`fetch` so accidental DOM/network usage in code
 * paths exercised by tests fails fast instead of silently swallowing errors.
 *
 * Modules are loaded in dependency order:
 *   Molecule -> SmilesParser/SmilesWriter/SmartsParser/SmartsMatch
 *            -> SMSDGraph -> SMSDRings -> SMSDVF2 -> SMSDMCS
 *            -> History
 */
'use strict';

if (typeof globalThis.window === 'undefined') globalThis.window = globalThis;

if (typeof globalThis.document === 'undefined') {
    var stubElement = function() {
        return {
            setAttribute: function() {},
            getAttribute: function() { return null; },
            appendChild: function() {},
            removeChild: function() {},
            style: {},
            textContent: '',
            innerHTML: '',
            children: [],
            classList: { add: function() {}, remove: function() {}, contains: function() { return false; } },
            getComputedTextLength: function() { return 0; },
            getBoundingClientRect: function() { return { x: 0, y: 0, width: 0, height: 0 }; }
        };
    };
    globalThis.document = {
        createElement: stubElement,
        createElementNS: stubElement,
        body: { appendChild: function() {}, removeChild: function() {} },
        documentElement: stubElement(),
        head: stubElement()
    };
}

if (typeof globalThis.fetch === 'undefined') {
    globalThis.fetch = function() {
        throw new Error('fetch not allowed in tests');
    };
}

var path = require('path');

function require_editor(name) {
    return require(path.join(__dirname, '..', 'editor', name + '.js'));
}

// Load editor modules once per process. Returns the populated globalThis so
// tests can pull Molecule, SmilesParser, ... off it.
var loaded = false;
function loadAll() {
    if (loaded) return globalThis;
    require_editor('Molecule');
    require_editor('SmilesParser');
    require_editor('SmilesWriter');
    require_editor('SmartsParser');
    require_editor('SmartsMatch');
    require_editor('CipStereo');
    require_editor('SMSDGraph');
    require_editor('SMSDRings');
    require_editor('SMSDVF2');
    require_editor('SMSDMCS');
    require_editor('History');
    loaded = true;
    return globalThis;
}

// Tiny test runner used by every tests/test_*.js file.
function makeRunner(label) {
    var passed = 0, failed = 0;
    var failures = [];
    var results = [];

    function test(name, fn) {
        try {
            fn();
            passed++;
            results.push({ name: name, pass: true });
            console.log('  ✓ ' + name);
        } catch (e) {
            failed++;
            var reason = (e && e.message) ? e.message : String(e);
            failures.push({ name: name, reason: reason, stack: e && e.stack });
            results.push({ name: name, pass: false, reason: reason });
            console.log('  ✗ ' + name + ' (' + reason + ')');
        }
    }

    function summary() {
        return { label: label, passed: passed, failed: failed, failures: failures, results: results };
    }

    return { test: test, summary: summary };
}

module.exports = {
    require_editor: require_editor,
    loadAll: loadAll,
    makeRunner: makeRunner
};
