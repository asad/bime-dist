/**
 * tests/test_v2_0_x_golden_reactions.js — golden-reactions regression corpus.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Runs BIME's atom-atom-mapping engine (RDT.mapReaction) over a curated
 * set of named reactions and asserts the expected shape of the result
 * (status, mapped-atom count, bond-change presence). The corpus is the
 * ground-truth set the v2.0.10 AAM-quality release was validated
 * against; the design target is 400 entries, growing incrementally
 * with each release.
 *
 * Each reaction is asserted INDIVIDUALLY by id so a failure pinpoints
 * the offending entry. A summary at the end of the run gives a per-class
 * pass-rate breakdown, which makes regressions in a specific reaction
 * class easy to spot at a glance.
 *
 * The corpus lives at tests/data/golden_reactions.json. Each entry
 * declares its class, SMILES, expected status, minimum mapped-atom
 * count, minimum bond-change count, and whether 'formed' / 'cleaved'
 * events are expected. To grow the corpus, add entries to that JSON;
 * the runner picks them up automatically.
 *
 * Adversarial entries (class starts with 'unbalanced') exercise the
 * v2.0.10 strict pre-check — they must return status='unbalanced'
 * with an empty mapping. Identity entries (class 'identity') exercise
 * the canary path where no bond changes are expected.
 */
'use strict';

var assert = require('assert');
var path   = require('path');
var fs     = require('fs');
var shim   = require(path.join(__dirname, 'shim.js'));
shim.loadAll();
require(path.join(__dirname, '..', 'editor', 'RDT.js'));

var CORPUS_PATH = path.join(__dirname, 'data', 'golden_reactions.json');
var CORPUS = JSON.parse(fs.readFileSync(CORPUS_PATH, 'utf8'));

var runner = shim.makeRunner('Golden reactions AAM corpus');
var test   = runner.test;

console.log('Golden reactions AAM corpus (target: ' +
    CORPUS.target_count + ', current: ' + CORPUS.reactions.length + ')');

function parse(smi) {
    var m = SmilesParser.parse(smi);
    if (m.parseErrors.length > 0) {
        throw new Error('parse error: ' + m.parseErrors.join('; '));
    }
    return m;
}

function countEvents(events, type) {
    var n = 0;
    for (var i = 0; i < events.length; i++) {
        if (events[i].type === type) { n++; }
    }
    return n;
}

// Tally per-class pass/fail counts so the end-of-run summary surfaces
// regressions in a single reaction family at a glance.
var perClass = {};
function tally(cls, ok) {
    if (!perClass[cls]) { perClass[cls] = { pass: 0, fail: 0 }; }
    if (ok) { perClass[cls].pass++; }
    else    { perClass[cls].fail++; }
}

CORPUS.reactions.forEach(function (rx) {
    var label = rx.id + ' [' + rx.class + '] — ' + rx.name;
    test(label, function () {
        var ok = true;
        try {
            var m = parse(rx.smiles);
            var res = RDT.mapReaction(m);

            // Status — exact match.
            assert.strictEqual(res.status, rx.expected.status,
                'status mismatch: expected ' + rx.expected.status +
                ', got ' + res.status +
                (res.elementDelta ? ' (elementDelta=' + JSON.stringify(res.elementDelta) + ')' : ''));

            if (rx.expected.status === 'mapped') {
                // Minimum mapped-atom count.
                var mapped = Object.keys(res.mapping || {}).length;
                assert.ok(mapped >= rx.expected.minMapping,
                    'expected at least ' + rx.expected.minMapping +
                    ' mapped atoms; got ' + mapped);

                // Bond-change presence + total count.
                var bc = res.bondChanges || [];
                var totalChanges = countEvents(bc, 'formed') +
                                   countEvents(bc, 'cleaved') +
                                   countEvents(bc, 'orderChange');
                assert.ok(totalChanges >= rx.expected.minBondChanges,
                    'expected at least ' + rx.expected.minBondChanges +
                    ' bond changes; got ' + totalChanges);

                if (rx.expected.hasFormed) {
                    assert.ok(countEvents(bc, 'formed') >= 1,
                        "expected at least one 'formed' event");
                }
                if (rx.expected.hasCleaved) {
                    assert.ok(countEvents(bc, 'cleaved') >= 1,
                        "expected at least one 'cleaved' event");
                }

                // Quality fields surface check — v2.0.10 invariants.
                assert.ok(typeof res.decisiveness === 'number',
                    'decisiveness must be numeric');
                assert.ok(res.decisiveness >= 0 && res.decisiveness <= 1,
                    'decisiveness ' + res.decisiveness + ' must be in [0, 1]');
                assert.ok(Array.isArray(res.candidates),
                    'candidates must be an array');
            } else if (rx.expected.status === 'unbalanced') {
                // The v2.0.10 strict pre-check: empty mapping + elementDelta.
                assert.deepStrictEqual(res.mapping, {},
                    'unbalanced reaction must return empty mapping');
                assert.ok(res.elementDelta && Object.keys(res.elementDelta).length > 0,
                    'unbalanced reaction must report elementDelta');
            }
        } catch (e) {
            ok = false;
            throw e;
        } finally {
            tally(rx.class, ok);
        }
    });
});

// Class summary AFTER all individual reaction tests — surfaced as its
// own pseudo-test so the runner's pass/fail count includes it (i.e. if
// any class drops to 0% pass, this 'test' fails).
test('per-class pass-rate summary', function () {
    var classes = Object.keys(perClass).sort();
    var totalPass = 0, totalFail = 0;
    var lines = [];
    for (var i = 0; i < classes.length; i++) {
        var c = classes[i];
        var stats = perClass[c];
        totalPass += stats.pass;
        totalFail += stats.fail;
        var pct = (stats.pass / (stats.pass + stats.fail) * 100).toFixed(0);
        lines.push('    ' + (pct === '100' ? '✓' : '✗') + ' [' + pct + '%] ' +
            c + ' (' + stats.pass + '/' + (stats.pass + stats.fail) + ')');
    }
    console.log('\n  Per-class pass-rate:');
    for (var j = 0; j < lines.length; j++) { console.log(lines[j]); }
    console.log('\n  Overall: ' + totalPass + '/' + (totalPass + totalFail) +
        ' (' + (totalPass / (totalPass + totalFail) * 100).toFixed(1) + '%)');
    assert.strictEqual(totalFail, 0,
        totalFail + ' golden reactions failed across ' + classes.length + ' classes');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
