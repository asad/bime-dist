/**
 * tests/test_v1_4_0_features.js — RDT v1.4.0 mol-mol pair highlight features.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Covers the v1.4.0 additions to RDT.mapReaction and Renderer:
 *
 *   (A) result.componentPairs — flat list of
 *       { reactantCompIdx, productCompIdx, mcsSize, paletteIndex } describing
 *       which reactant component pairs with which product component, and
 *       which palette colour the renderer should use (paletteIndex 0..N-1
 *       for paired components, -1 for spectator / unpaired components from
 *       stoichiometric mismatch).
 *
 *   (B) result.confidence — a normalised mapping confidence in [0, 1]:
 *           confidence = mappedAtoms / totalHeavyAtoms ×
 *                        (1 / (1 + 0.1 × bondChanges))
 *       Identity reactions hit 1.0; textbook reactions land in [0.5, 1.0];
 *       wholly disjoint reactions stay below 0.5.
 *
 *   (C) Renderer per-atom halo overlay — Renderer.componentPairColors() helper,
 *       Renderer.colorAtoms toggle, halo SVG <circle> elements with
 *       opacity = 0.5 and the correct palette fill. Halos render BEFORE bonds
 *       (in bgLayer) so atoms still display on top.
 *
 * Plain Node, no external dependencies.
 */
'use strict';

var assert = require('assert');
var shim = require('./shim.js');
shim.loadAll();
require('../editor/RDT.js');
require('../editor/Layout.js');
require('../editor/Renderer.js');

var runner = shim.makeRunner('RDT v1.4.0 mol-mol features');
var test = runner.test;

console.log('RDT v1.4.0 mol-mol features');

function parse(rxnSmiles) {
    var m = SmilesParser.parse(rxnSmiles);
    assert.strictEqual(m.parseErrors.length, 0,
        'parse error in ' + rxnSmiles + ': ' + m.parseErrors.join('; '));
    return m;
}

function pairedCount(componentPairs) {
    if (!componentPairs) return 0;
    var n = 0;
    for (var i = 0; i < componentPairs.length; i++) {
        if (componentPairs[i].paletteIndex >= 0) { n++; }
    }
    return n;
}

// =========================================================================
// (A) componentPairs population
// =========================================================================

test('A1. Single-component reaction CCO>>CC=O yields 1 component pair', function() {
    var r = RDT.mapReaction(parse('CCO>>CC=O'));
    assert.ok(Array.isArray(r.componentPairs), 'componentPairs should be an array');
    assert.strictEqual(pairedCount(r.componentPairs), 1);
    assert.strictEqual(r.componentPairs[0].reactantCompIdx, 0);
    assert.strictEqual(r.componentPairs[0].productCompIdx, 0);
    assert.strictEqual(r.componentPairs[0].paletteIndex, 0);
});

test('A2. Two-reactant + two-product esterification yields >=2 pairs', function() {
    // v1.4.2: leftover-atom rescue may add a 3rd paired entry when an
    // anchored unmapped atom (the leaving OH oxygen) gets matched to the
    // water product O. The rescue is opt-out via useLeftoverRescue:false;
    // here we test the default-on behaviour.
    var r = RDT.mapReaction(parse('CC(=O)O.OCC>>CC(=O)OCC.O'));
    assert.ok(pairedCount(r.componentPairs) >= 2,
        'expected >=2 paired entries (rescue may add one); got ' +
        pairedCount(r.componentPairs));
    // Distinct paletteIndex starting from 0; 0..N-1 should all appear.
    var idxs = r.componentPairs
        .map(function(p) { return p.paletteIndex; })
        .filter(function(p) { return p >= 0; })
        .sort(function(a, b) { return a - b; });
    assert.ok(idxs.length >= 2);
    assert.strictEqual(idxs[0], 0);
});

test('A3. Spectator [Na+] participates in pairing — >=3 paired', function() {
    // v1.4.2: leftover-atom rescue may add a 4th paired component when a
    // bipartite-paired component has unmapped atoms anchored to a mapped
    // neighbour (e.g. carbonyl O of acetic acid that the acyl-side MCS
    // didn't catch). Three is the minimum; rescue may exceed it.
    var r = RDT.mapReaction(parse('CC(=O)O.OCC.[Na+]>>CC(=O)OCC.O.[Na+]'));
    assert.ok(pairedCount(r.componentPairs) >= 3,
        'at least 3 components should be paired (rescue may add more); got ' +
        pairedCount(r.componentPairs));
});

test('A4. Stoichiometric mismatch (3 reactants -> 2 products) emits paletteIndex=-1', function() {
    // v1.4.2: rescue may pair more components than the bipartite step alone
    // when anchored atoms exist. The structural invariant we still preserve:
    // there is at least one paletteIndex = -1 entry per truly-unpaired
    // reactant component (atoms with no mapped neighbour anywhere).
    var r = RDT.mapReaction(parse('CC.OO.NN>>CO.CO'));
    assert.ok(pairedCount(r.componentPairs) >= 2,
        'expected at least 2 paired components; got ' + pairedCount(r.componentPairs));
    var unpaired = r.componentPairs.filter(function(p) { return p.paletteIndex === -1; });
    assert.ok(unpaired.length >= 1,
        'expected at least one unpaired component with paletteIndex=-1; got ' + unpaired.length);
});

test('A5. Reaction with no reactant side returns componentPairs = []', function() {
    var r = RDT.mapReaction(parse('>>CCO'));
    assert.ok(Array.isArray(r.componentPairs));
    assert.strictEqual(r.componentPairs.length, 0);
});

// =========================================================================
// (B) confidence
// =========================================================================

test('B1. Identity reaction CCO>>CCO has confidence === 1.0', function() {
    var r = RDT.mapReaction(parse('CCO>>CCO'));
    assert.strictEqual(typeof r.confidence, 'number');
    assert.ok(r.confidence > 0.999, 'expected confidence ~1.0; got ' + r.confidence);
});

test('B2. Disjoint reaction CC.OO.NN>>CO.CO has confidence < 0.5', function() {
    var r = RDT.mapReaction(parse('CC.OO.NN>>CO.CO'));
    assert.ok(r.confidence < 0.5,
        'disjoint reaction should score < 0.5; got ' + r.confidence);
});

test('B3. Esterification CC(=O)O.OCC>>CC(=O)OCC.O confidence in [0.5, 1.0]', function() {
    var r = RDT.mapReaction(parse('CC(=O)O.OCC>>CC(=O)OCC.O'));
    assert.ok(r.confidence >= 0.5 && r.confidence <= 1.0,
        'esterification should score in [0.5, 1.0]; got ' + r.confidence);
});

// =========================================================================
// (C) Palette assignment
// =========================================================================

test('C1. First N pairs get paletteIndex 0..N-1', function() {
    var r = RDT.mapReaction(parse('CC.OO.NN>>CC.OO.NN'));
    var paired = r.componentPairs.filter(function(p) { return p.paletteIndex >= 0; });
    var idxs = paired.map(function(p) { return p.paletteIndex; }).sort(function(a, b) { return a - b; });
    assert.deepStrictEqual(idxs, [0, 1, 2]);
});

test('C2. Palette index can exceed length (renderer wraps modulo)', function() {
    // 12 distinct components on each side -> paletteIndex up to 11. The
    // renderer's halo lookup wraps modulo COMPONENT_PAIR_PALETTE.length.
    var rxn = 'C.O.N.S.[F].[Cl].[Br].[I].P.CC.CO.CN>>C.O.N.S.[F].[Cl].[Br].[I].P.CC.CO.CN';
    var r = RDT.mapReaction(parse(rxn));
    var maxIdx = -1;
    for (var i = 0; i < r.componentPairs.length; i++) {
        var p = r.componentPairs[i];
        if (p.paletteIndex > maxIdx) { maxIdx = p.paletteIndex; }
    }
    assert.ok(maxIdx >= Renderer.COMPONENT_PAIR_PALETTE.length,
        'expected at least one paletteIndex >= ' + Renderer.COMPONENT_PAIR_PALETTE.length +
        '; got max=' + maxIdx);
    // Compute renderer halo colour for the highest paletteIndex — must wrap.
    var wrapped = Renderer.COMPONENT_PAIR_PALETTE[maxIdx % Renderer.COMPONENT_PAIR_PALETTE.length];
    assert.ok(typeof wrapped === 'string' && wrapped.charAt(0) === '#');
});

test('C3. Unpaired components get paletteIndex === -1', function() {
    var r = RDT.mapReaction(parse('CC.OO.NN>>CO.CO'));
    var unpaired = r.componentPairs.filter(function(p) { return p.paletteIndex === -1; });
    assert.strictEqual(unpaired.length, 1);
    assert.strictEqual(unpaired[0].mcsSize, 0);
});

// =========================================================================
// (D) Renderer integration
// =========================================================================

test('D1. Renderer.componentPairColors(componentPairs) returns N+1 entries', function() {
    var pairs = [
        { paletteIndex: 0 },
        { paletteIndex: 1 },
        { paletteIndex: -1 }
    ];
    var colors = Renderer.componentPairColors(pairs);
    assert.strictEqual(colors.length, pairs.length + 1,
        'expected N+1 colours (pair colours + trailing neutral grey)');
    // Pair colours come from the palette.
    assert.strictEqual(colors[0], Renderer.COMPONENT_PAIR_PALETTE[0]);
    assert.strictEqual(colors[1], Renderer.COMPONENT_PAIR_PALETTE[1]);
    // -1 yields the neutral grey halo.
    assert.strictEqual(colors[2], Renderer.COMPONENT_PAIR_NEUTRAL);
    // Trailing entry is always the neutral grey for index -1 lookups.
    assert.strictEqual(colors[3], Renderer.COMPONENT_PAIR_NEUTRAL);
});

// The remaining D-tests need a real Renderer with a DOM — set up a minimal
// container shim and exercise the halo paths. Re-uses tests/shim.js's
// document stub.
function makeContainer() {
    var children = [];
    return {
        appendChild: function(c) { children.push(c); return c; },
        removeChild: function(c) {
            var i = children.indexOf(c);
            if (i >= 0) children.splice(i, 1);
            return c;
        },
        children: children
    };
}

// Build a minimal recording stub for SVG nodes so we can count halo circles.
// The shim's existing document stub returns shared dummy objects, so we
// install a richer stub that tracks tag, attributes, and parent/child links.
function installRecordingDom() {
    var saved = global.document;
    function makeNode(tag) {
        return {
            tagName: tag,
            children: [],
            _attrs: {},
            style: {},
            textContent: '',
            classList: { add: function() {}, remove: function() {}, contains: function() { return false; } },
            setAttribute: function(k, v) { this._attrs[k] = v; },
            getAttribute: function(k) { return this._attrs[k] !== undefined ? this._attrs[k] : null; },
            appendChild: function(c) { c._parent = this; this.children.push(c); return c; },
            removeChild: function(c) {
                var i = this.children.indexOf(c);
                if (i >= 0) this.children.splice(i, 1);
                return c;
            },
            insertBefore: function(c, ref) { c._parent = this; this.children.push(c); return c; },
            getComputedTextLength: function() { return 0; },
            getBoundingClientRect: function() { return { x: 0, y: 0, width: 600, height: 400, left: 0, top: 0 }; },
            firstChild: null,
            get _firstChild() { return this.children[0] || null; }
        };
    }
    Object.defineProperty(makeNode({}), 'firstChild', {});
    global.document = {
        createElement: function(tag) { return makeNode(tag); },
        createElementNS: function(_ns, tag) { return makeNode(tag); },
        body: { appendChild: function() {}, removeChild: function() {} },
        documentElement: makeNode('html'),
        head: makeNode('head')
    };
    return saved;
}

function restoreDom(saved) { global.document = saved; }

// Walk an SVG tree and collect all <circle> elements (so we can find halos).
function collectCircles(root) {
    var out = [];
    function walk(node) {
        if (!node || !node.children) return;
        for (var i = 0; i < node.children.length; i++) {
            var c = node.children[i];
            if (c && c.tagName === 'circle') { out.push(c); }
            walk(c);
        }
    }
    walk(root);
    return out;
}

function buildRendererForReaction(rxnSmiles) {
    var savedDom = installRecordingDom();
    try {
        var rxn = parse(rxnSmiles);
        // Run AAM — also assigns mapNumbers and computes componentPairs.
        var result = RDT.mapReaction(rxn);
        // Translate componentCompIdx -> reactantAtomIds/productAtomIds the
        // renderer expects. Mirrors MolEditor._runRdtAutoMap.
        var pal = Renderer.COMPONENT_PAIR_PALETTE;
        var atomIdPairs = RDT.deriveComponentPairs(result, pal);
        var compPairs = result.componentPairs || [];
        var rendererPairs = atomIdPairs.slice();
        for (var ci = 0; ci < rendererPairs.length; ci++) {
            var entry = rendererPairs[ci];
            var match = null;
            for (var cj = 0; cj < compPairs.length; cj++) {
                if (compPairs[cj].reactantCompIdx === entry.reactantComponentIdx &&
                    compPairs[cj].productCompIdx === entry.productComponentIdx) {
                    match = compPairs[cj]; break;
                }
            }
            entry.paletteIndex = match ? match.paletteIndex : ci;
        }

        var container = makeContainer();
        var renderer = new Renderer(container, 600, 400);
        renderer.setMolecule(rxn);
        renderer.componentPairs = rendererPairs;
        renderer.showComponentPairs = true;
        return { renderer: renderer, savedDom: savedDom, atomCount: rxn.atoms.length };
    } catch (e) {
        restoreDom(savedDom);
        throw e;
    }
}

test('D2. Renderer.colorAtoms = false suppresses halo circles', function() {
    var ctx = buildRendererForReaction('CCO>>CC=O');
    try {
        ctx.renderer.colorAtoms = false;
        ctx.renderer.render();
        var circles = collectCircles(ctx.renderer.bgLayer);
        var halos = circles.filter(function(c) {
            return c._attrs && c._attrs['class'] === 'bime-pair-halo';
        });
        assert.strictEqual(halos.length, 0, 'colorAtoms=false should emit zero halos');
    } finally {
        restoreDom(ctx.savedDom);
    }
});

test('D3. Renderer halos appear with correct fill colour when colorAtoms = true', function() {
    var ctx = buildRendererForReaction('CCO>>CC=O');
    try {
        ctx.renderer.colorAtoms = true;
        ctx.renderer.render();
        var circles = collectCircles(ctx.renderer.bgLayer);
        var halos = circles.filter(function(c) {
            return c._attrs && c._attrs['class'] === 'bime-pair-halo';
        });
        assert.ok(halos.length > 0, 'expected at least one halo circle');
        // Single component pair (paletteIndex 0) -> first palette tint.
        assert.strictEqual(halos[0]._attrs.fill, Renderer.COMPONENT_PAIR_PALETTE[0],
            'halo fill should match the first palette colour');
    } finally {
        restoreDom(ctx.savedDom);
    }
});

test('D4. Halo opacity is 0.5', function() {
    var ctx = buildRendererForReaction('CCO>>CC=O');
    try {
        ctx.renderer.colorAtoms = true;
        ctx.renderer.render();
        var circles = collectCircles(ctx.renderer.bgLayer);
        var halos = circles.filter(function(c) {
            return c._attrs && c._attrs['class'] === 'bime-pair-halo';
        });
        assert.ok(halos.length > 0, 'expected halo circles');
        assert.strictEqual(halos[0]._attrs.opacity, '0.5',
            'halo opacity must be exactly 0.5');
    } finally {
        restoreDom(ctx.savedDom);
    }
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
