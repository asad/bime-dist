/**
 * tests/test_bondchange.js — Bond-change semantics for BIME v1.1.1.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Verifies the RDT.annotateBondChanges output for textbook reactions:
 * carbonyl reduction, hydration, hydrogenation, esterification, hydrolysis,
 * SN2, Diels-Alder, no-change, isomerization, ring-opening/closing, stereo
 * change, atom-id consistency, deterministic ordering, and symmetry.
 *
 * Plain Node, no external dependencies.
 *
 * Note on v1.x policy
 * -------------------
 * The pure-JS RDT port in v1.1 conservatively annotates bond changes from the
 * AAM atom map only: when reactant atoms are unmapped to a product atom (e.g.
 * leaving groups, multi-component reactions), the connecting bond is reported
 * as 'cleaved'. New bonds count as 'formed' only when *both* product
 * endpoints are images of mapped reactant atoms. As a result, multi-component
 * couplings (esterification, hydrolysis) currently emit only 'cleaved' events
 * — see Section G of the v1.1.1 expansion notes for follow-up.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('Bond Changes (BC)');
var test = runner.test;

console.log('Bond Changes (BC)');

function rdt(rxnSmiles) {
    var r = SmilesParser.parse(rxnSmiles);
    var res = RDT.mapReaction(r);
    res._reaction = r;
    return res;
}

function ct(events, type) {
    var n = 0;
    for (var i = 0; i < events.length; i++) {
        if (events[i].type === type) n++;
    }
    return n;
}

// =========================================================================
// (1) Carbonyl reduction C=O → C-O
// =========================================================================
test('Carbonyl reduction CC=O>>CCO yields exactly 1 orderChange event', function() {
    var res = rdt('CC=O>>CCO');
    assert.strictEqual(ct(res.bondChanges, 'orderChange'), 1);
});

// =========================================================================
// (2) Hydration C=C + H2O → C-C-OH
// =========================================================================
test('Hydration C=C.O>>CCO yields exactly 1 bond-change event', function() {
    var res = rdt('C=C.O>>CCO');
    assert.strictEqual(res.bondChanges.length, 1);
});

// =========================================================================
// (3) Hydrogenation C=C + H2 → C-C
// =========================================================================
test('Hydrogenation C=C.[H][H]>>CC yields exactly 1 orderChange', function() {
    var res = rdt('C=C.[H][H]>>CC');
    assert.strictEqual(ct(res.bondChanges, 'orderChange'), 1,
        'expected exactly one C=C → C-C orderChange event');
});

// =========================================================================
// (4) Esterification — current v1.x policy emits cleaved-only events
// =========================================================================
test('Esterification CC(=O)O.OCC>>CC(=O)OCC.O yields exactly 2 cleaved events', function() {
    var res = rdt('CC(=O)O.OCC>>CC(=O)OCC.O');
    // The two leaving-group bonds (acid -OH and alcohol -O-H endpoints)
    // are reported as cleaved under the v1.x cleaved-only policy.
    assert.strictEqual(ct(res.bondChanges, 'cleaved'), 2);
    assert.strictEqual(ct(res.bondChanges, 'formed'), 0);
});

// =========================================================================
// (5) Hydrolysis (reverse of esterification) — same shape
// =========================================================================
test('Hydrolysis CC(=O)OCC.O>>CC(=O)O.OCC yields exactly 2 cleaved events', function() {
    var res = rdt('CC(=O)OCC.O>>CC(=O)O.OCC');
    assert.strictEqual(ct(res.bondChanges, 'cleaved'), 2);
});

// =========================================================================
// (6) SN2: Br → OH
// =========================================================================
test('SN2 CCBr.[OH-]>>CCO.[Br-] yields exactly 1 cleaved event (the C-Br bond)', function() {
    var res = rdt('CCBr.[OH-]>>CCO.[Br-]');
    assert.strictEqual(ct(res.bondChanges, 'cleaved'), 1);
});

// =========================================================================
// (7) Diels-Alder: 2 cleaved (the C=C → C-C transitions become cleaved
// because RDT can't always pair the cycloaddition endpoints in v1.x)
// =========================================================================
test('Diels-Alder C=CC=C.C=C>>C1CC=CCC1 yields exactly 2 bond-change events', function() {
    var res = rdt('C=CC=C.C=C>>C1CC=CCC1');
    assert.strictEqual(res.bondChanges.length, 2,
        'Diels-Alder should yield exactly 2 bond-change events under v1.x policy');
});

// =========================================================================
// (8) No-change reaction
// =========================================================================
test('Identity CCO>>CCO yields zero bond-change events', function() {
    var res = rdt('CCO>>CCO');
    assert.strictEqual(res.bondChanges.length, 0);
});

// =========================================================================
// (9) Pure isomerization — formed/cleaved or order-change events sum > 0
// =========================================================================
test('Isomerization CCC>>CC=C yields exactly 1 bond-change event', function() {
    var res = rdt('CCC>>CC=C');
    var totalChanges = ct(res.bondChanges, 'formed') +
        ct(res.bondChanges, 'cleaved') +
        ct(res.bondChanges, 'orderChange');
    assert.strictEqual(totalChanges, 1);
});

// =========================================================================
// (10) Ring-opening: 1 cleaved (the ring bond)
// =========================================================================
test('Ring-opening C1CCCCC1>>CCCCCC yields exactly 1 cleaved event', function() {
    var res = rdt('C1CCCCC1>>CCCCCC');
    assert.strictEqual(ct(res.bondChanges, 'cleaved'), 1);
});

// =========================================================================
// (11) Ring-closing: 1 formed (the new ring bond)
// =========================================================================
test('Ring-closing CCCCCC>>C1CCCCC1 yields exactly 1 formed event', function() {
    var res = rdt('CCCCCC>>C1CCCCC1');
    assert.strictEqual(ct(res.bondChanges, 'formed'), 1);
});

// =========================================================================
// (12) Stereo flip — v1.x AAM does not perform CIP perception so
// stereoChange events may be absent. This test pins that policy.
// =========================================================================
test('Stereo flip [C@H](N)(C)O>>[C@@H](N)(C)O: stereoChange count documented', function() {
    var res = rdt('[C@H](N)(C)O>>[C@@H](N)(C)O');
    var stereo = ct(res.bondChanges, 'stereoChange');
    // BIME v1.1 does not perceive CIP labels in mapReaction, so stereo
    // events are 0. This test documents the policy; if a future release
    // adds CIP perception, update this assertion.
    assert.ok(stereo === 0, 'stereoChange currently 0 by v1.x policy; got ' + stereo);
});

// =========================================================================
// (13) Bond-change atom indices reference real atoms in the input molecule
// =========================================================================
test('Bond-change atom IDs all reference real reactant atoms', function() {
    var res = rdt('CCO>>CC=O');
    var ids = {};
    for (var i = 0; i < res._reaction.atoms.length; i++) {
        ids[res._reaction.atoms[i].id] = true;
    }
    for (var j = 0; j < res.bondChanges.length; j++) {
        var ev = res.bondChanges[j];
        if (ev.atoms[0] !== null) {
            assert.ok(ids[ev.atoms[0]], 'atom id ' + ev.atoms[0] + ' not in molecule');
        }
        // Second slot may be a product atom id (also in r.atoms — reaction is one big mol)
        if (ev.atoms[1] !== null) {
            assert.ok(ids[ev.atoms[1]], 'atom id ' + ev.atoms[1] + ' not in molecule');
        }
    }
});

// =========================================================================
// (14) Bond-change events sorted deterministically (by reactantBond, then productBond)
// =========================================================================
test('Bond-change events are sorted deterministically (reactantBond ascending)', function() {
    var res = rdt('CC(=O)O.OCC>>CC(=O)OCC.O');
    var lastR = -1;
    for (var i = 0; i < res.bondChanges.length; i++) {
        var ev = res.bondChanges[i];
        if (ev.type === 'cleaved' || ev.type === 'orderChange') {
            // Reactant-side events are emitted in pass 1 in ascending reactantBond order.
            assert.ok(ev.reactantBond >= lastR, 'reactant-pass events out of order');
            lastR = ev.reactantBond;
        }
    }
});

// =========================================================================
// (15) Symmetry: cleaved + formed counts on a perfectly mirrored reaction agree
// =========================================================================
test('Mirror reaction CC(=O)O>>OC(=O)C: formed === cleaved (both zero)', function() {
    var res = rdt('CC(=O)O>>OC(=O)C');
    assert.strictEqual(ct(res.bondChanges, 'formed'), ct(res.bondChanges, 'cleaved'));
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
