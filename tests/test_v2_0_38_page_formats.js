/**
 * tests/test_v2_0_38_page_formats.js — A4/Letter page-format sizing
 * for the pathway canvas (R2; v2.0.38).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/PageFormats.js');

var runner = shim.makeRunner('Page Formats (v2.0.38)');
var test = runner.test;
var PF = PageFormats;

console.log('Page Formats (v2.0.38)');

// ---------------------------------------------------------------------------
// A. format catalog
// ---------------------------------------------------------------------------
test('A1 FORMATS has screen, a4, and letter', function() {
    assert.ok(PF.FORMATS.screen, 'screen present');
    assert.ok(PF.FORMATS.a4,     'a4 present');
    assert.ok(PF.FORMATS.letter, 'letter present');
});

test('A2 every format has portrait + landscape sub-sizes', function() {
    Object.keys(PF.FORMATS).forEach(function (k) {
        var f = PF.FORMATS[k];
        assert.ok(f.portrait  && f.portrait.w  > 0 && f.portrait.h  > 0, k + ' portrait');
        assert.ok(f.landscape && f.landscape.w > 0 && f.landscape.h > 0, k + ' landscape');
    });
});

// ---------------------------------------------------------------------------
// B. exact pixel dimensions match the standards
// ---------------------------------------------------------------------------
test('B1 screen format is the legacy 1200×620 viewBox (backward compat)', function() {
    var s = PF.sizeFor('screen', 'landscape');
    assert.strictEqual(s.w, 1200);
    assert.strictEqual(s.h, 620);
});

test('B2 A4 dimensions match 210×297 mm at 96 DPI', function() {
    // 210 mm × 96 DPI / 25.4 mm/in = 793.7 ≈ 794
    // 297 mm × 96 DPI / 25.4 mm/in = 1122.5 ≈ 1123
    assert.deepStrictEqual(PF.sizeFor('a4', 'portrait'),  { w: 794, h: 1123 });
    assert.deepStrictEqual(PF.sizeFor('a4', 'landscape'), { w: 1123, h: 794 });
});

test('B3 US Letter dimensions match 8.5×11 in at 96 DPI', function() {
    assert.deepStrictEqual(PF.sizeFor('letter', 'portrait'),  { w: 816, h: 1056 });
    assert.deepStrictEqual(PF.sizeFor('letter', 'landscape'), { w: 1056, h: 816 });
});

// ---------------------------------------------------------------------------
// C. viewBoxFor / isKnown
// ---------------------------------------------------------------------------
test('C1 viewBoxFor returns the SVG viewBox string', function() {
    assert.strictEqual(PF.viewBoxFor('a4', 'landscape'), '0 0 1123 794');
});

test('C2 isKnown accepts valid pairs, rejects invalid', function() {
    assert.strictEqual(PF.isKnown('a4', 'portrait'), true);
    assert.strictEqual(PF.isKnown('a4', 'landscape'), true);
    assert.strictEqual(PF.isKnown('letter', 'portrait'), true);
    assert.strictEqual(PF.isKnown('screen', 'landscape'), true);
    assert.strictEqual(PF.isKnown('a3', 'portrait'), false, 'unknown format');
    assert.strictEqual(PF.isKnown('a4', 'diagonal'), false, 'unknown orientation');
    assert.strictEqual(PF.isKnown(undefined, undefined), false);
});

// ---------------------------------------------------------------------------
// D. fallbacks (defensive)
// ---------------------------------------------------------------------------
test('D1 sizeFor falls back to defaults on unknown format', function() {
    var s = PF.sizeFor('quarto', 'landscape');
    // screen landscape = 1200 × 620
    assert.deepStrictEqual(s, { w: 1200, h: 620 });
});

test('D2 sizeFor falls back to landscape on unknown orientation', function() {
    var s = PF.sizeFor('a4', 'banana');
    // a4 landscape
    assert.deepStrictEqual(s, { w: 1123, h: 794 });
});

// ---------------------------------------------------------------------------
// E. formatList
// ---------------------------------------------------------------------------
test('E1 formatList enumerates {key, label} for every format', function() {
    var list = PF.formatList();
    assert.ok(Array.isArray(list));
    assert.ok(list.length >= 3, 'at least screen + a4 + letter');
    var keys = list.map(function (x) { return x.key; }).sort();
    assert.ok(keys.indexOf('a4') >= 0);
    assert.ok(keys.indexOf('letter') >= 0);
    assert.ok(keys.indexOf('screen') >= 0);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
