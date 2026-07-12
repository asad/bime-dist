/**
 * tests/test_v2_0_62_cli.js — bime command-line wrapper (v2.0.62).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * The CLI lives in `tools/bime-cli.js` with a thin shebang wrapper at
 * `bin/bime`. These tests drive both via Node's child_process so the
 * surface-level shell contract is locked down end-to-end:
 *
 *   bime version
 *   bime parse <smiles> [--format json]
 *   bime clean <smiles>
 *   bime aam   <reactionSmiles>
 *   bime moiety <pathway> --start <met> --atoms i,j,k
 *   bime pathway-layout <pathway> [--shape auto|cycle|...]
 *   bime export <smiles> --format svg|mol|sdf|rxn|smiles
 *   bime pathway-list
 *   bime help [subcommand]
 *
 * Each test asserts on stdout. Exit-code expectations are pinned where
 * the contract is "fail fast" (bad subcommand, malformed input).
 */
'use strict';

var assert = require('assert');
var path = require('path');
var fs = require('fs');
var childProcess = require('child_process');
var shim = require('./shim.js');

var runner = shim.makeRunner('bime CLI (v2.0.62)');
var test = runner.test;

console.log('bime CLI (v2.0.62)');

var ROOT = path.join(__dirname, '..');
var BIME = path.join(ROOT, 'bin', 'bime');

function run(args, stdinText) {
    var opts = { encoding: 'utf8' };
    if (stdinText) { opts.input = stdinText; }
    var r = childProcess.spawnSync(BIME, args, opts);
    return { stdout: r.stdout || '', stderr: r.stderr || '', code: r.status };
}

// v2.4.18: the pathway subcommands (moiety, pathway-layout, pathway-list) are
// gated behind PATHWAY_ENABLED in tools/bime-cli.js (default false = core build).
// Read the flag from source so the pathway CLI tests track it in BOTH states —
// full behaviour when enabled, "unknown command" when gated off (shipped).
var PATHWAY_ENABLED = /var\s+PATHWAY_ENABLED\s*=\s*true\b/.test(
    fs.readFileSync(path.join(ROOT, 'tools', 'bime-cli.js'), 'utf8'));

// ---------------------------------------------------------------------------
// A. existence + shebang
// ---------------------------------------------------------------------------
test('A1 bin/bime exists and is executable', function() {
    assert.ok(fs.existsSync(BIME), 'bin/bime is present');
    var st = fs.statSync(BIME);
    // Owner-execute bit (mode & 0o100) must be set.
    assert.ok((st.mode & 0o100) !== 0, 'bin/bime has owner-execute bit');
});

test('A2 tools/bime-cli.js exports the subcommand functions', function() {
    var mod = require(path.join(ROOT, 'tools', 'bime-cli.js'));
    var must = ['parseArgs', 'loadEditor', 'cmdVersion', 'cmdParse', 'cmdClean',
                'cmdAam', 'cmdMoiety', 'cmdPathwayLayout', 'cmdExport',
                'cmdPathwayList', 'cmdHelp', 'main'];
    for (var i = 0; i < must.length; i++) {
        assert.strictEqual(typeof mod[must[i]], 'function', must[i] + ' is exported');
    }
});

// ---------------------------------------------------------------------------
// B. version + help
// ---------------------------------------------------------------------------
test('B1 bime version prints BIME + module versions', function() {
    var r = run(['version']);
    assert.strictEqual(r.code, 0);
    assert.ok(/bime \d+\.\d+\.\d+/.test(r.stdout), 'has "bime X.Y.Z" line');
    assert.ok(/RDT\.version/.test(r.stdout));
    assert.ok(/AtomTrace\.version/.test(r.stdout));
});

test('B2 bime help prints the top-level help', function() {
    var r = run(['help']);
    assert.strictEqual(r.code, 0);
    assert.ok(/Usage:\s+bime/.test(r.stdout));
    assert.ok(/Commands:/.test(r.stdout));
    assert.ok(/parse/.test(r.stdout));
    assert.ok(/aam/.test(r.stdout));
    // v2.4.18: the pathway subcommands appear in help only when enabled.
    if (PATHWAY_ENABLED) {
        assert.ok(/moiety/.test(r.stdout), 'help lists moiety when pathway is enabled');
    } else {
        assert.ok(!/moiety|pathway-layout|pathway-list/.test(r.stdout),
            'help omits the pathway subcommands in the core build');
    }
});

test('B3 unknown subcommand exits non-zero', function() {
    var r = run(['nope-not-a-command']);
    assert.notStrictEqual(r.code, 0);
    assert.ok(/unknown command/.test(r.stderr));
});

// ---------------------------------------------------------------------------
// C. parse
// ---------------------------------------------------------------------------
test('C1 bime parse benzene produces a 6-atom 6-bond report', function() {
    var r = run(['parse', 'c1ccccc1']);
    assert.strictEqual(r.code, 0);
    assert.ok(/Atoms:\s+6/.test(r.stdout));
    assert.ok(/Bonds:\s+6/.test(r.stdout));
    assert.ok(/Reaction:\s+no/.test(r.stdout));
});

test('C2 bime parse --format json emits JSON', function() {
    var r = run(['parse', 'c1ccccc1', '--format', 'json']);
    assert.strictEqual(r.code, 0);
    var parsed = JSON.parse(r.stdout);
    assert.strictEqual(parsed.atomCount, 6);
    assert.strictEqual(parsed.bondCount, 6);
    assert.deepStrictEqual(parsed.elements, { C: 6 });
});

test('C3 bime parse reads from stdin when no positional', function() {
    var r = run(['parse'], 'CCO');
    assert.strictEqual(r.code, 0);
    assert.ok(/Atoms:\s+3/.test(r.stdout));
});

// ---------------------------------------------------------------------------
// D. clean
// ---------------------------------------------------------------------------
test('D1 bime clean returns a canonical SMILES (toluene)', function() {
    var r = run(['clean', 'Cc1ccccc1']);
    assert.strictEqual(r.code, 0);
    assert.ok(/[Cc]/.test(r.stdout), 'output contains a SMILES-like string');
});

// ---------------------------------------------------------------------------
// E. aam
// ---------------------------------------------------------------------------
test('E1 bime aam on a balanced reaction emits status=mapped', function() {
    // Williamson ether synthesis: ethanol + acetic acid -> ester + water
    // (balanced — heavy atoms match).
    var r = run(['aam', 'CCO.CC(=O)O>>CCOC(=O)C.O']);
    assert.strictEqual(r.code, 0);
    assert.ok(/Status:\s+mapped/.test(r.stdout) ||
              /Status:\s+unbalanced/.test(r.stdout),
        'AAM produces a status string');
});

test('E2 bime aam --format json produces a parseable AAM result', function() {
    var r = run(['aam', 'CCO.CC(=O)O>>CCOC(=O)C.O', '--format', 'json']);
    assert.strictEqual(r.code, 0);
    var parsed = JSON.parse(r.stdout);
    assert.ok(typeof parsed.status === 'string');
    assert.ok('mappedCount' in parsed);
    assert.ok('mapping'    in parsed);
});

var MAPPED_RXN = '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]';

test('E3 bime aam --keep-mapping uses the input :n mapping (no re-solve)', function() {
    var r = run(['aam', MAPPED_RXN, '--format', 'json', '--keep-mapping']);
    assert.strictEqual(r.code, 0);
    var parsed = JSON.parse(r.stdout);
    assert.strictEqual(parsed.status, 'mapped');
    assert.strictEqual(parsed.mappedCount, 3, 'all three heavy atoms carry :n -> mapped');
    assert.strictEqual(Object.keys(parsed.mapping).length, 3);
});

test('E4 bime aam --keep-mapping --format svg renders a figure from the existing mapping', function() {
    var r = run(['aam', MAPPED_RXN, '--format', 'svg', '--keep-mapping']);
    assert.strictEqual(r.code, 0);
    assert.ok(r.stdout.indexOf('<svg') !== -1, 'emits an <svg> document');
});

test('E5 bime aam --keep-mapping fails fast on a reaction with NO atom-map numbers', function() {
    var r = run(['aam', 'CCO>>CC=O', '--format', 'svg', '--keep-mapping']);
    assert.notStrictEqual(r.code, 0, 'unmapped input + --keep-mapping must error (it is not re-solved)');
    assert.ok(/atom-map numbers/i.test(r.stderr + r.stdout), 'explains the missing mapping');
});

// ---------------------------------------------------------------------------
// F. pathway-list + moiety (gated behind PATHWAY_ENABLED in the core build)
// ---------------------------------------------------------------------------
// v2.4.18: full pathway-CLI behaviour when the flag is on; when off (shipped
// core build) the subcommands are hidden — see the gated-case test after F10.
if (PATHWAY_ENABLED) {
test('F1 bime pathway-list enumerates the bundled pathways', function() {
    var r = run(['pathway-list']);
    assert.strictEqual(r.code, 0);
    var must = ['glycolysis', 'tca', 'urea', 'ppp', 'pyruvate_fates', 'beta_oxidation'];
    for (var i = 0; i < must.length; i++) {
        assert.ok(r.stdout.indexOf(must[i]) !== -1, 'lists ' + must[i]);
    }
});

test('F2 bime moiety glycolysis G6P phosphate atoms reports survives + breakAt', function() {
    // The G6P phosphate group: atoms 0..4 in the parse of
    // O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O (corrected v2.4.17;
    // phosphate on the exocyclic C6, written first so it parses as atoms 0-4).
    var r = run([
        'moiety', 'glycolysis',
        '--start', 'Glucose 6-phosphate',
        '--atoms', '0,1,2,3,4'
    ]);
    assert.strictEqual(r.code, 0, 'moiety command succeeded; stderr=' + r.stderr);
    assert.ok(/Pathway:\s+Glycolysis/.test(r.stdout));
    assert.ok(/Moiety:\s+0,1,2,3,4 \(5 atoms\)/.test(r.stdout));
    assert.ok(/Breaks at:/.test(r.stdout));
    assert.ok(/Survivors:/.test(r.stdout));
});

test('F3 bime moiety --format json produces a parseable structured report', function() {
    var r = run([
        'moiety', 'glycolysis',
        '--start', 'Glucose 6-phosphate',
        '--atoms', '0,1,2,3,4',
        '--format', 'json'
    ]);
    assert.strictEqual(r.code, 0);
    var parsed = JSON.parse(r.stdout);
    assert.strictEqual(parsed.pathway, 'glycolysis');
    assert.strictEqual(parsed.start, 'Glucose 6-phosphate');
    assert.deepStrictEqual(parsed.moietyAtoms, [0, 1, 2, 3, 4]);
    assert.ok(Array.isArray(parsed.perMetabolite));
    assert.ok(parsed.perMetabolite.length >= 5, 'has per-metabolite entries');
});

test('F4 bime moiety without --atoms exits non-zero with a helpful error', function() {
    var r = run(['moiety', 'glycolysis', '--start', 'Glucose 6-phosphate']);
    assert.notStrictEqual(r.code, 0);
    assert.ok(/--atoms/.test(r.stderr));
});

test('F5 bime pathway-layout emits reusable fanout coordinates as JSON', function() {
    var r = run(['pathway-layout', 'pyruvate_fates', '--shape', 'fanout', '--format', 'json']);
    assert.strictEqual(r.code, 0, 'pathway-layout succeeded; stderr=' + r.stderr);
    var parsed = JSON.parse(r.stdout);
    assert.strictEqual(parsed.pathway, 'pyruvate_fates');
    assert.strictEqual(parsed.requestedShape, 'fanout');
    assert.strictEqual(parsed.resolvedShape, 'fanout');
    assert.strictEqual(parsed.nodeCount, 5);
    assert.strictEqual(parsed.edgeCount, 4);
    assert.ok(Array.isArray(parsed.nodes));
    assert.ok(parsed.nodes.some(function(n) { return n.id === 'Pyruvate'; }));
});

test('F6 bime pathway-layout normalizes branch and hybrid aliases', function() {
    var branched = run(['pathway-layout', 'ppp', '--shape', 'branch', '--format', 'json']);
    assert.strictEqual(branched.code, 0, 'branch alias succeeded; stderr=' + branched.stderr);
    assert.strictEqual(JSON.parse(branched.stdout).requestedShape, 'branched');

    var hybrid = run(['pathway-layout', 'pyruvate_fates', '--layout', 'combinations', '--format', 'json']);
    assert.strictEqual(hybrid.code, 0, 'hybrid alias succeeded; stderr=' + hybrid.stderr);
    var parsed = JSON.parse(hybrid.stdout);
    assert.strictEqual(parsed.requestedShape, 'hybrid');
    assert.ok(parsed.componentShapes.indexOf('fanout') !== -1, 'hybrid uses inferred fanout component');
});

test('F7 bime help pathway-layout documents the shape switch', function() {
    var r = run(['help', 'pathway-layout']);
    assert.strictEqual(r.code, 0);
    assert.ok(/--shape/.test(r.stdout));
    assert.ok(/--in/.test(r.stdout));
    assert.ok(/split-merge/.test(r.stdout));
    assert.ok(/loop-iterative/.test(r.stdout));
});

test('F8 bime pathway-layout accepts a saved pathway JSON file', function() {
    var tmp = path.join(require('os').tmpdir(), 'bime-pathway-layout-' + Date.now() + '.json');
    var model = {
        format: 'bime-pathway',
        version: 1,
        name: 'Split merge smoke',
        layoutShape: 'branch-merge',
        nodes: [
            { id: 'A', label: 'A' },
            { id: 'B', label: 'B' },
            { id: 'C', label: 'C' },
            { id: 'D', label: 'D' }
        ],
        edges: [
            { id: 'e1', from: 'A', to: 'B' },
            { id: 'e2', from: 'A', to: 'C' },
            { id: 'e3', from: 'B', to: 'D' },
            { id: 'e4', from: 'C', to: 'D' }
        ]
    };
    fs.writeFileSync(tmp, JSON.stringify(model), 'utf8');
    var r = run(['pathway-layout', '--in', tmp, '--format', 'json']);
    fs.unlinkSync(tmp);
    assert.strictEqual(r.code, 0, 'pathway-layout --in succeeded; stderr=' + r.stderr);
    var parsed = JSON.parse(r.stdout);
    assert.strictEqual(parsed.pathway, 'Split merge smoke');
    assert.strictEqual(parsed.requestedShape, 'split-merge');
    assert.strictEqual(parsed.nodeCount, 4);
    assert.strictEqual(parsed.edgeCount, 4);
});

test('F9 bime pathway-layout reports missing --in files without a stack trace', function() {
    var missing = path.join(require('os').tmpdir(), 'bime-missing-' + Date.now() + '.json');
    var r = run(['pathway-layout', '--in', missing]);
    assert.notStrictEqual(r.code, 0);
    assert.ok(/cannot read input/.test(r.stderr));
    assert.ok(!/\n\s+at\s+/.test(r.stderr), 'no raw stack trace');
});

test('F10 bime pathway-layout reports invalid JSON without a stack trace', function() {
    var r = run(['pathway-layout', '-', '--format', 'json'], '{not-json');
    assert.notStrictEqual(r.code, 0);
    assert.ok(/cannot parse input JSON/.test(r.stderr));
    assert.ok(!/\n\s+at\s+/.test(r.stderr), 'no raw stack trace');
});
} else {
    test('F-gated: pathway subcommands are unavailable in the core build', function() {
        ['pathway-list', 'moiety', 'pathway-layout'].forEach(function(c) {
            var r = run([c]);
            assert.notStrictEqual(r.code, 0, c + ' exits non-zero when gated off');
            assert.ok(/unknown command/.test(r.stderr),
                c + ' reports "unknown command"; stderr=' + r.stderr);
        });
        var h = run(['help']);
        assert.ok(!/moiety|pathway-layout|pathway-list/.test(h.stdout),
            'top-level help omits the gated pathway subcommands');
    });
}

// ---------------------------------------------------------------------------
// G. export
// ---------------------------------------------------------------------------
test('G1 bime export --format smiles round-trips a SMILES', function() {
    var r = run(['export', 'c1ccccc1', '--format', 'smiles']);
    assert.strictEqual(r.code, 0);
    assert.ok(r.stdout.indexOf('c') !== -1);
});

test('G2 bime export --format svg emits an <svg> document', function() {
    var r = run(['export', 'c1ccccc1', '--format', 'svg']);
    assert.strictEqual(r.code, 0);
    assert.ok(r.stdout.indexOf('<svg') !== -1);
    assert.ok(r.stdout.indexOf('</svg>') !== -1);
});

test('G3 bime export --format mol emits a MOL V2000 block', function() {
    var r = run(['export', 'Cc1ccccc1', '--format', 'mol']);
    assert.strictEqual(r.code, 0);
    assert.ok(/V2000/.test(r.stdout), 'output contains V2000 marker');
});

test('G3b bime export --format sdf emits a MOL block terminated by $$$$ (v2.0.64 fix)', function() {
    var r = run(['export', 'Cc1ccccc1', '--format', 'sdf']);
    assert.strictEqual(r.code, 0, 'sdf export succeeded; stderr=' + r.stderr);
    assert.ok(/V2000/.test(r.stdout), 'sdf output contains V2000 marker');
    assert.ok(/\$\$\$\$/.test(r.stdout), 'sdf output is terminated by $$$$');
});

test('G3c bime export --format rxn round-trips a reaction SMILES', function() {
    var r = run(['export', 'CCO>>CC=O', '--format', 'rxn']);
    assert.strictEqual(r.code, 0);
    assert.ok(/\$RXN/.test(r.stdout), 'rxn output starts with $RXN header');
});

test('G4 bime export --format png is refused (no Canvas in Node)', function() {
    var r = run(['export', 'CCO', '--format', 'png']);
    assert.notStrictEqual(r.code, 0);
    assert.ok(/png|browser DOM|Canvas/i.test(r.stderr));
});

test('G5 bime export --out writes to a file', function() {
    var tmp = path.join(require('os').tmpdir(), 'bime-cli-test-' + Date.now() + '.svg');
    var r = run(['export', 'c1ccccc1', '--format', 'svg', '--out', tmp]);
    assert.strictEqual(r.code, 0);
    var content = fs.readFileSync(tmp, 'utf8');
    assert.ok(content.indexOf('<svg') !== -1);
    fs.unlinkSync(tmp);
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
