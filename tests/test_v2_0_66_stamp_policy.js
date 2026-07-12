/**
 * tests/test_v2_0_66_stamp_policy.js — codify the v2.0.65 export-stamp
 * standing rule (v2.0.66).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * v2.0.66 is a policy-documentation release. The user confirmed (with
 * the phrase "always remember 'Encrypted' for public domain
 * InshahAllah") that BIME's notion of "encrypted for public domain"
 * means **tamper-evident signing**, NOT asymmetric encryption. The
 * test below makes that policy text non-removable from `AGENTS.md`,
 * `editor/ExportStamp.js`, and `USER_GUIDE.md` so future releases
 * inherit it. A failing test here means the policy was edited out
 * by accident — restore the original phrasing or update the test
 * deliberately.
 *
 * Plain Node, no DOM. Source-shape only.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require('./shim.js');
shim.loadAll();

var runner = shim.makeRunner('Stamp policy (v2.0.66)');
var test = runner.test;

console.log('Stamp policy (v2.0.66)');

// This suite asserts policy text inside AGENTS.md, an internal-only file
// that is NOT part of the public distribution. When AGENTS.md is absent
// (public distribution), the suite skips cleanly.
if (!fs.existsSync(path.join(__dirname, '..', 'AGENTS.md'))) {
    console.log('  (AGENTS.md not present — internal policy suite skipped)');
    module.exports = function () {
        return {
            label: 'Stamp policy (v2.0.66) [skipped]',
            passed: 0, failed: 0, failures: [], results: []
        };
    };
    return;
}

function readRepoFile(rel) {
    return fs.readFileSync(path.join(__dirname, '..', rel), 'utf8');
}

// ---------------------------------------------------------------------------
// A. AGENTS.md carries the clarification
// ---------------------------------------------------------------------------
test('A1 AGENTS.md has an Encrypted-for-public-domain clarification section', function () {
    var md = readRepoFile('AGENTS.md');
    assert.ok(md.indexOf('Encrypted-for-public-domain') !== -1 ||
              md.indexOf('encrypted for public domain') !== -1 ||
              md.indexOf('encrypted-for-public-domain') !== -1,
        'AGENTS.md mentions the encrypted-for-public-domain phrase explicitly');
});

test('A2 AGENTS.md states that the stamp is tamper-evident, not encrypted', function () {
    var md = readRepoFile('AGENTS.md');
    assert.ok(/tamper.?evident/i.test(md),
        'AGENTS.md describes the stamp as tamper-evident');
    assert.ok(md.indexOf('not asymmetric') !== -1 ||
              md.indexOf('not encrypted') !== -1 ||
              md.indexOf('NOT asymmetric') !== -1,
        'AGENTS.md states explicitly that the stamp is not asymmetric encryption');
});

test('A3 AGENTS.md captures the Apache-2.0-must-stay-readable reasoning', function () {
    var md = readRepoFile('AGENTS.md');
    assert.ok(md.indexOf('Apache-2.0') !== -1,
        'AGENTS.md references the Apache-2.0 licence in the policy');
    assert.ok(/readable by anyone|read,? modify/i.test(md),
        'AGENTS.md states the readability rationale');
});

test('A4 AGENTS.md acknowledges asymmetric signatures as a future roadmap item', function () {
    var md = readRepoFile('AGENTS.md');
    assert.ok(/Ed25519|RSA/i.test(md),
        'AGENTS.md mentions Ed25519/RSA as the asymmetric option');
    assert.ok(/roadmap|future release/i.test(md),
        'AGENTS.md flags asymmetric signing as a future roadmap item');
});

// ---------------------------------------------------------------------------
// B. ExportStamp.js docblock carries the same clarification
// ---------------------------------------------------------------------------
test('B1 ExportStamp.js docblock carries the public-domain clarification', function () {
    var js = readRepoFile('editor/ExportStamp.js');
    // Normalise the docblock's word-wrap so multi-line phrases match.
    // Strip the leading ` * ` of each line, then collapse all whitespace.
    var flat = js.replace(/\n\s*\*\s?/g, ' ').replace(/\s+/g, ' ');
    assert.ok(flat.indexOf('tamper-evident') !== -1 ||
              flat.indexOf('tamper evident') !== -1,
        'ExportStamp.js docblock describes the stamp as tamper-evident');
    assert.ok(/encrypted for public domain/i.test(flat),
        'ExportStamp.js docblock mentions the phrase "encrypted for public domain"');
});

test('B2 ExportStamp.js docblock explicitly says NOT confidentiality', function () {
    var js = readRepoFile('editor/ExportStamp.js');
    var flat = js.replace(/\n\s*\*\s?/g, ' ').replace(/\s+/g, ' ');
    assert.ok(/confidentiality/i.test(flat),
        'ExportStamp.js mentions confidentiality');
    assert.ok(/attribution[\s\S]{1,30}integrity/i.test(flat) ||
              /attribution \+ integrity/i.test(flat),
        'ExportStamp.js promises attribution + integrity');
});

// ---------------------------------------------------------------------------
// C. USER_GUIDE.md surfaces the same explanation
// ---------------------------------------------------------------------------
test('C1 USER_GUIDE.md has a "what that actually means" section on the stamp', function () {
    var md = readRepoFile('USER_GUIDE.md');
    assert.ok(/Encrypted for public domain.{0,40}what that actually means/i.test(md),
        'USER_GUIDE.md carries the explanatory section');
    assert.ok(/tamper.?evident/i.test(md),
        'USER_GUIDE.md surfaces the tamper-evident framing');
});

test('C2 USER_GUIDE.md spells out that confidentiality is NOT promised', function () {
    var md = readRepoFile('USER_GUIDE.md');
    assert.ok(/Confidentiality is not promised/i.test(md) ||
              /does not.{0,30}lock/i.test(md),
        'USER_GUIDE.md explicitly disclaims confidentiality');
});

// ---------------------------------------------------------------------------
// D. functional contract still holds — stamp is non-empty + idempotent
// ---------------------------------------------------------------------------
test('D1 ExportStamp.stampSvg still produces a non-empty stamp', function () {
    require('../editor/ExportStamp.js');
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="60" viewBox="0 0 100 60"></svg>';
    var stamped = ExportStamp.stampSvg(svg, { version: '2.0.66' });
    assert.ok(stamped.indexOf('<metadata id="bime-stamp">') !== -1);
    assert.ok(stamped.indexOf('<bime:version>2.0.66</bime:version>') !== -1);
});

test('D2 ExportStamp.stampText still idempotent under double stamping', function () {
    require('../editor/ExportStamp.js');
    var once  = ExportStamp.stampText('mol', 'mol content');
    var twice = ExportStamp.stampText('mol', once);
    var begins = (twice.match(/# BIME-STAMP-BEGIN/g) || []).length;
    assert.strictEqual(begins, 1, 'double stamping does not nest');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
