/**
 * tests/test_v1_2_0_features.js — RDT v1.2.0 additive AAM improvements.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Covers the three v1.2.0 additions to RDT.mapReaction:
 *
 *   (A) options.includeHydrogens (default true) — emits 'hydrogenChange'
 *       events when a mapped atom's implicit-H count differs (oxidation,
 *       reduction, dehydration). Implementation cites the same RDT 2016
 *       paper (Rahman SA et al., Bioinformatics 32(13):2065-66) and uses
 *       BIME's existing Molecule.calcHydrogens(atomId).
 *
 *   (B) options.useBipartitePostPass (default false) — refines the four
 *       strategies' greedy component pairing with a max-weight bipartite
 *       matching post-pass (textbook Munkres-Kuhn assignment, citing
 *       Kuhn HW. Naval Research Logistics Quarterly 1955;2:83-97 and
 *       Munkres J. SIAM J 1957;5:32-38). Re-pairs only when total MCS
 *       coverage strictly improves.
 *
 *   (C) options.includeStereo (default true) — runs CIP perception inside
 *       mapReaction so 'stereoChange' events fire when a mapped stereocentre
 *       flips R<->S during the reaction (e.g. SN2 inversion).
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');

var runner = shim.makeRunner('RDT v1.2.0 features');
var test = runner.test;

console.log('RDT v1.2.0 features');

function parse(rxnSmiles) {
    var m = SmilesParser.parse(rxnSmiles);
    assert.strictEqual(m.parseErrors.length, 0,
        'parse error in ' + rxnSmiles + ': ' + m.parseErrors.join('; '));
    return m;
}

function ct(events, type) {
    var n = 0;
    for (var i = 0; i < events.length; i++) { if (events[i].type === type) { n++; } }
    return n;
}

function deltaSum(events, type) {
    var s = 0;
    for (var i = 0; i < events.length; i++) {
        if (events[i].type === type && typeof events[i].deltaH === 'number') {
            s += events[i].deltaH;
        }
    }
    return s;
}

// =========================================================================
// (A) includeHydrogens — implicit-H accounting
// =========================================================================

test('A1. Ethanol oxidation CCO>>CC=O emits exactly 2 hydrogenChange events (-1, -1)', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r);
    assert.strictEqual(ct(res.bondChanges, 'hydrogenChange'), 2,
        'expected 2 hydrogenChange events; got ' + JSON.stringify(res.bondChanges));
    assert.strictEqual(deltaSum(res.bondChanges, 'hydrogenChange'), -2,
        'sum of deltaH should be -2 for ethanol -> acetaldehyde');
});

test('A2. Ethanol oxidation orderChange + hydrogenChange total = 3', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r);
    assert.strictEqual(res.bondChanges.length, 3,
        'expected 1 orderChange + 2 hydrogenChange = 3 events');
});

test('A3. Carbonyl reduction CC=O.[H][H]>>CCO emits 2 hydrogenChange events (+1, +1)', function() {
    var r = parse('CC=O.[H][H]>>CCO');
    var res = RDT.mapReaction(r);
    var hc = ct(res.bondChanges, 'hydrogenChange');
    assert.strictEqual(hc, 2, 'expected 2 hydrogenChange; got ' + hc);
    assert.strictEqual(deltaSum(res.bondChanges, 'hydrogenChange'), 2);
});

test('A4. SN2 CCBr.[OH-]>>CCO.[Br-] emits 0 hydrogenChange (heavy-atom skeleton preserved)', function() {
    var r = parse('CCBr.[OH-]>>CCO.[Br-]');
    var res = RDT.mapReaction(r);
    assert.strictEqual(ct(res.bondChanges, 'hydrogenChange'), 0);
});

test('A5. Aromatic chlorination produces a -1 deltaH on the substituted ring carbon', function() {
    var r = parse('c1ccccc1.ClCl>>c1ccc(Cl)cc1.[ClH]');
    var res = RDT.mapReaction(r);
    var hc = ct(res.bondChanges, 'hydrogenChange');
    assert.ok(hc >= 1, 'expected at least 1 hydrogenChange; got ' + hc);
    var foundMinusOne = false;
    for (var i = 0; i < res.bondChanges.length; i++) {
        if (res.bondChanges[i].type === 'hydrogenChange' && res.bondChanges[i].deltaH === -1) {
            foundMinusOne = true; break;
        }
    }
    assert.ok(foundMinusOne, 'expected a -1 deltaH on the substituted ring carbon');
});

test('A6. options.includeHydrogens=false suppresses all hydrogenChange events', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r, { includeHydrogens: false });
    assert.strictEqual(ct(res.bondChanges, 'hydrogenChange'), 0);
});

test('A7. Identity CCO>>CCO emits 0 hydrogenChange (no atoms change H count)', function() {
    var r = parse('CCO>>CCO');
    var res = RDT.mapReaction(r);
    assert.strictEqual(ct(res.bondChanges, 'hydrogenChange'), 0);
});

test('A8. hydrogenChange events carry deltaH, beforeH, afterH, atom, productAtom', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r);
    var any = null;
    for (var i = 0; i < res.bondChanges.length; i++) {
        if (res.bondChanges[i].type === 'hydrogenChange') { any = res.bondChanges[i]; break; }
    }
    assert.ok(any, 'expected at least one hydrogenChange');
    assert.strictEqual(typeof any.deltaH, 'number');
    assert.strictEqual(typeof any.beforeH, 'number');
    assert.strictEqual(typeof any.afterH, 'number');
    assert.strictEqual(any.afterH - any.beforeH, any.deltaH);
    assert.strictEqual(typeof any.atom, 'number');
    assert.strictEqual(typeof any.productAtom, 'number');
});

test('A9. Hydration C=C.O>>CCO emits 1 orderChange + >=1 hydrogenChange', function() {
    // v1.4.2 rescue: the water oxygen is now mapped, so its H count delta
    // (2 -> 1) is reported alongside the carbon's H gain (2 -> 3). Pre-
    // v1.4.2 only the carbon's hydrogenChange was visible (water O was an
    // orphan reactant component). The orderChange (C=C -> C-C) is the
    // canonical heavy-atom event in either era.
    var r = parse('C=C.O>>CCO');
    var res = RDT.mapReaction(r);
    assert.strictEqual(ct(res.bondChanges, 'orderChange'), 1);
    assert.ok(ct(res.bondChanges, 'hydrogenChange') >= 1,
        'expected at least 1 hydrogenChange (C of C=C gains H, water O loses H)');
});

test('A10. Isomerization CCC>>CC=C emits hydrogenChange entries with summing deltaH = -2', function() {
    var r = parse('CCC>>CC=C');
    var res = RDT.mapReaction(r);
    assert.strictEqual(deltaSum(res.bondChanges, 'hydrogenChange'), -2,
        'CH3-CH2-CH3 -> CH3-CH=CH2: 2 H removed across the new double-bond carbons');
});

// =========================================================================
// (B) useBipartitePostPass — max-weight bipartite component pairing
// =========================================================================

test('B1. _bipartiteComponentPairing exists on the public API', function() {
    assert.strictEqual(typeof RDT._bipartiteComponentPairing, 'function');
});

test('B2. Bipartite pairing on a 1x1 weight matrix returns the single pair', function() {
    var pairs = RDT._bipartiteComponentPairing(1, 1, [[5]]);
    assert.deepStrictEqual(pairs, [{ reactantIdx: 0, productIdx: 0, mcsSize: 5 }]);
});

test('B3. Bipartite pairing prefers max-weight 4+1 over greedy-min 1+3 on a 2x2 matrix', function() {
    // Esterification-shaped matrix: greedy MIN picks the smallest cell first
    // (R[0]xP[1]=1) then is forced into R[1]xP[0]=3 -> total 4. Bipartite
    // matches R[0]xP[0]=4 + R[1]xP[1]=1 -> total 5.
    var weights = [
        [4, 1],
        [3, 1]
    ];
    var pairs = RDT._bipartiteComponentPairing(2, 2, weights);
    var total = 0;
    for (var i = 0; i < pairs.length; i++) { total += pairs[i].mcsSize; }
    assert.strictEqual(total, 5, 'bipartite total should be 5; got ' + total);
});

test('B4. useBipartitePostPass=true on esterification with strategies=[MIN] re-pairs to 5 mapped atoms', function() {
    // v1.4.2: pass useLeftoverRescue:false to test bipartite in isolation
    // — the v1.4.2 leftover-atom rescue would otherwise extend the mapping
    // to 6+ atoms by catching the leaving water O.
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    var res = RDT.mapReaction(r, { useBipartitePostPass: true, strategies: ['MIN'], useLeftoverRescue: false });
    assert.strictEqual(res.strategy, 'MIN');
    assert.strictEqual(Object.keys(res.mapping).length, 5,
        'MIN+bipartite should map 5 atoms (vs. 4 without)');
    assert.strictEqual(res.bipartiteApplied, true,
        'MIN strategy should have triggered the bipartite re-pairing on this reaction');
});

test('B5a. v1.3.0 default (bipartite ON) on esterification with strategies=[MIN] re-pairs to 5 atoms', function() {
    // v1.3.0: useBipartitePostPass defaults to true. Bipartite fires on
    // multi-component esterification because the greedy MIN strategy picks a
    // sub-optimal component pairing that bipartite strictly improves.
    // v1.4.2: pass useLeftoverRescue:false to keep the bipartite-only count.
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    var res = RDT.mapReaction(r, { strategies: ['MIN'], useLeftoverRescue: false });
    assert.strictEqual(Object.keys(res.mapping).length, 5,
        'v1.3.0 default-on bipartite should map 5 atoms (vs 4 for greedy MIN)');
});

test('B5b. Explicit opt-out (v1.2.x legacy) on esterification with strategies=[MIN] keeps the 4-atom greedy', function() {
    // Pass useBipartitePostPass:false to restore the v1.2.x default-off
    // behaviour, pinned so anyone depending on the v1.2.x mapping can
    // deterministically reproduce it. v1.4.2 also adds useLeftoverRescue:false
    // so the rescue pass doesn't extend the mapping past the v1.2.x baseline.
    var r = parse('CC(=O)O.OCC>>CC(=O)OCC.O');
    var res = RDT.mapReaction(r, { strategies: ['MIN'], useBipartitePostPass: false, useLeftoverRescue: false });
    assert.strictEqual(Object.keys(res.mapping).length, 4,
        'with bipartite + rescue both off, MIN-only must give the 4-atom greedy mapping');
});

test('B6. Bipartite is a no-op on a reaction whose greedy is already optimal', function() {
    // Identity: greedy already perfect (3 atoms map). Bipartite mustn't change result.
    var r1 = parse('CCO>>CCO');
    var r2 = parse('CCO>>CCO');
    var withBip = RDT.mapReaction(r1, { useBipartitePostPass: true });
    var withoutBip = RDT.mapReaction(r2);
    assert.strictEqual(Object.keys(withBip.mapping).length, Object.keys(withoutBip.mapping).length);
});

test('B7. Bipartite handles unequal component counts (rN != pN) by zero-padding', function() {
    // 2 reactants, 1 product -> 1 winning pair (the larger MCS)
    var pairs = RDT._bipartiteComponentPairing(2, 1, [[3], [5]]);
    assert.strictEqual(pairs.length, 1);
    assert.strictEqual(pairs[0].reactantIdx, 1);
    assert.strictEqual(pairs[0].productIdx, 0);
    assert.strictEqual(pairs[0].mcsSize, 5);
});

test('B8. Bipartite drops zero-weight pairs (no spurious matches)', function() {
    // 2x2 with one zero column: only one true pair survives.
    var pairs = RDT._bipartiteComponentPairing(2, 2, [[0, 4], [3, 0]]);
    var total = 0;
    for (var i = 0; i < pairs.length; i++) { total += pairs[i].mcsSize; }
    assert.strictEqual(total, 7);
    assert.strictEqual(pairs.length, 2);
});

// =========================================================================
// (C) includeStereo — CIP integration in AAM
// =========================================================================

test('C1. SN2-style direct stereo flip [C@H](N)(C)O>>[C@@H](N)(C)O emits 1 stereoChange (R<->S)', function() {
    var r = parse('[C@H](N)(C)O>>[C@@H](N)(C)O');
    var res = RDT.mapReaction(r);
    assert.strictEqual(ct(res.bondChanges, 'stereoChange'), 1,
        'expected 1 stereoChange; got ' + JSON.stringify(res.bondChanges));
});

test('C2. stereoChange event carries before/after, atom, productAtom, mapNumber fields', function() {
    var r = parse('[C@H](N)(C)O>>[C@@H](N)(C)O');
    var res = RDT.mapReaction(r);
    var ev = null;
    for (var i = 0; i < res.bondChanges.length; i++) {
        if (res.bondChanges[i].type === 'stereoChange') { ev = res.bondChanges[i]; break; }
    }
    assert.ok(ev, 'expected a stereoChange event');
    assert.ok(ev.before === 'R' || ev.before === 'S');
    assert.ok(ev.after === 'R' || ev.after === 'S');
    assert.notStrictEqual(ev.before, ev.after);
    assert.strictEqual(typeof ev.atom, 'number');
    assert.strictEqual(typeof ev.productAtom, 'number');
});

test('C3. options.includeStereo=false suppresses stereoChange events', function() {
    var r = parse('[C@H](N)(C)O>>[C@@H](N)(C)O');
    var res = RDT.mapReaction(r, { includeStereo: false });
    assert.strictEqual(ct(res.bondChanges, 'stereoChange'), 0);
});

test('C4. Achiral reaction CCO>>CC=O emits 0 stereoChange', function() {
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r);
    assert.strictEqual(ct(res.bondChanges, 'stereoChange'), 0);
});

test('C5. Mapping does not crash if CipStereo is unavailable (best-effort behaviour)', function() {
    // Temporarily hide CipStereo; mapReaction must still complete without throw.
    var saved = global.CipStereo;
    try {
        global.CipStereo = undefined;
        var r = parse('[C@H](N)(C)O>>[C@@H](N)(C)O');
        var res = RDT.mapReaction(r);
        assert.ok(res, 'mapReaction returned a result');
        assert.ok(Array.isArray(res.bondChanges));
        // With CIP unavailable, no stereoChange events fire (cipLabel stays empty).
        assert.strictEqual(ct(res.bondChanges, 'stereoChange'), 0);
    } finally {
        global.CipStereo = saved;
    }
});

test('C6. CIP perception inside mapReaction does not regress mapping count', function() {
    // CCO>>CC=O still maps all 3 heavy atoms with stereo perception on.
    var r = parse('CCO>>CC=O');
    var res = RDT.mapReaction(r);
    assert.strictEqual(Object.keys(res.mapping).length, 3);
});

test('C7. RDT.version matches current bundle version', function() {
    // RDT.version tracks the live bundle. The v1.2.0 features tested in
    // this file ship in any version >= 1.2.0 — but the version stamp must
    // match versions.json / SMSDVersion.bimeVersion exactly so consumers
    // can trust runtime introspection.
    assert.strictEqual(RDT.version, '1.8.15');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
