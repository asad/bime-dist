/**
 * tests/test_v2_4_10_dark_canvas.js — true dark drawing canvas (v2.4.10).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Dark mode previously kept the editor canvas a light "paper" sheet (hardcoded
 * #333 ink). v2.4.10 makes the drawing surface genuinely dark in dark mode with a
 * light ink, while saturated CPK heteroatom colours keep their hue. The mechanism
 * is theme-toggle-safe with NO re-render:
 *   - Renderer inks bonds / default labels / arrows with `currentColor`, and sets
 *     the drawing svg's inline `color` to `var(--mol-ink, #333)`. Near-black atom
 *     labels (Carbon / Hydrogen / secondary grey) route through currentColor too;
 *     CPK heteroatom colours are emitted explicitly (they read on dark).
 *   - CSS themes the DRAWING svg (svg[role="img"]) via CUSTOM PROPERTIES that apply
 *     regardless of the background cascade: light pins --mol-ink:#333; dark sets
 *     --color-bg:#101c2e + --mol-ink:#cbd5e1. The inline styles consume them, so a
 *     single data-theme swap re-themes the whole depiction live.
 *   - The old dark light-sheet rule no longer forces the editor svg / .me-draw-area
 *     (only the export preview, template thumbs and pathway canvas stay light).
 *
 * Source-shape contract on editor/Renderer.js + css/workbench.css (the level
 * test_v2_4_8 uses), with CSS comments stripped so prose can't spoof a regex.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('True dark drawing canvas (v2.4.10)');
var test = runner.test;
console.log('True dark drawing canvas (v2.4.10)');

function readRepoFile(rel) { return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8'); }
var RJS = readRepoFile('editor/Renderer.js');
var CSS = readRepoFile('css/workbench.css').replace(/\/\*[\s\S]*?\*\//g, '');

// ---------------------------------------------------------------------------
// A. Renderer inks with currentColor + drives the svg colour from --mol-ink.
// ---------------------------------------------------------------------------
test('A1 the drawing svg inks via currentColor and takes its colour from --mol-ink', function () {
    assert.ok(/this\.svg\.style\.color\s*=\s*'var\(--mol-ink, #333\)'/.test(RJS),
        "the svg's inline colour is var(--mol-ink, #333) — themed by the dark-canvas CSS");
    assert.ok(/this\.svg\.style\.background\s*=\s*'var\(--color-bg, white\)'/.test(RJS),
        'the svg background stays var(--color-bg) (dark via the CSS custom property)');
});

test('A2 bonds, the _line/_text defaults and the reaction arrow default to currentColor', function () {
    assert.ok(/bond\.selected \? SELECT_COLOR : 'currentColor'/.test(RJS), 'unselected bonds ink with currentColor');
    assert.ok(/var arrowColor = 'currentColor'/.test(RJS), 'the reaction arrow inks with currentColor');
    assert.ok(/setAttribute\('stroke', stroke \|\| 'currentColor'\)/.test(RJS), '_line stroke defaults to currentColor');
    assert.ok(/setAttribute\('fill', fill \|\| 'currentColor'\)/.test(RJS), '_text fill defaults to currentColor');
    // The old hardcoded #333 ink defaults are gone (only the _NEAR_BLACK lookup keys remain).
    assert.strictEqual(/stroke \|\| '#333'/.test(RJS), false, 'no #333 stroke default remains');
    assert.strictEqual(/fill \|\| '#333'/.test(RJS), false, 'no #333 fill default remains');
});

test('A3 near-black atom labels route through currentColor; CPK heteroatoms keep their hue', function () {
    assert.ok(/function _rLabelColor\(c\)/.test(RJS), '_rLabelColor helper exists');
    assert.ok(/_NEAR_BLACK\[c\.toLowerCase\(\)\]\)\s*\?\s*'currentColor'\s*:\s*c/.test(RJS),
        '_rLabelColor returns currentColor for near-black colours, else the colour itself (CPK)');
    assert.ok(/_rLabelColor\(elem\.color\)/.test(RJS), 'the atom label uses _rLabelColor(elem.color)');
    // Carbon + Hydrogen (#333333) are in the near-black set so they follow the ink.
    assert.ok(/'#333333':\s*1/.test(RJS), 'Carbon/Hydrogen #333333 is treated as near-black ink');
    // The secondary hydrogen label + mol name now follow the ink too.
    assert.ok(/hStr \+ hCountStr, 'currentColor'/.test(RJS), 'the H label inks with currentColor');
    assert.ok(/nameText, 'currentColor'/.test(RJS), 'the molecule-name label inks with currentColor');
});

// ---------------------------------------------------------------------------
// B. CSS themes the drawing svg via custom properties (cascade-safe).
// ---------------------------------------------------------------------------
test('B1 light pins --mol-ink to #333 on the drawing svg', function () {
    assert.ok(/\.wb-editor-wrap svg\[role="img"\]\s*\{\s*--mol-ink:\s*#333;\s*\}/.test(CSS),
        'the base (light) rule pins --mol-ink:#333 on svg[role="img"]');
});

test('B2 dark sets --color-bg dark + --mol-ink light on the drawing svg', function () {
    var i = CSS.indexOf('[data-theme="dark"] .wb-editor-wrap svg[role="img"] {');
    assert.ok(i !== -1, 'a dark rule targets the drawing svg');
    var block = CSS.slice(i, CSS.indexOf('}', i) + 1);
    assert.ok(/--color-bg:\s*#101c2e/.test(block), 'dark sets --color-bg to a dark surface');
    assert.ok(/--mol-ink:\s*#cbd5e1/.test(block), 'dark sets --mol-ink to a light ink');
});

test('B3 neither dark path (explicit toggle OR OS prefers-dark) light-sheets the editor drawing svg', function () {
    // BOTH the `[data-theme="dark"]` rules and the `@media (prefers-color-scheme:
    // dark) :root:not([data-theme="light"])` rules used to list the editor svg as a
    // grouped branch of the #f7fbff light sheet (".wb-editor-wrap svg,"). The OS
    // path's `background:#f7fbff !important` was the real cascade culprit. Both must
    // now leave the editor svg to the custom-prop dark rules. The base rule uses
    // ".wb-editor-wrap svg {" (no comma), so a leftover light-sheet branch is the
    // only thing that produces ".wb-editor-wrap svg,".
    assert.strictEqual(/\.wb-editor-wrap svg\s*,/.test(CSS), false,
        'no light-sheet rule lists ".wb-editor-wrap svg," as a grouped branch (editor svg freed in dark)');
    // The #f7fbff sheet still exists for the non-editor surfaces (pathway canvas).
    assert.ok(/#f7fbff/.test(CSS) && /wb-pathway-canvas/.test(CSS),
        'the #f7fbff sheet still applies to the pathway canvas / export preview');
    // The OS-dark path also carries the custom-prop dark-canvas rule.
    assert.ok(/:root:not\(\[data-theme="light"\]\) \.wb-editor-wrap svg\[role="img"\]/.test(CSS),
        'the @media (prefers-color-scheme: dark) path also themes svg[role="img"] via custom props');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
