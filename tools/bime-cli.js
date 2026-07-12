/**
 * tools/bime-cli.js — command-line wrapper around BIME's editor engine.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Headless Node access to the same engine the workbench uses in the
 * browser. Subcommands:
 *
 *   bime version
 *   bime parse    <smiles>                                 — SMILES diagnostic
 *   bime clean    <smiles> [--out file]                    — Layout.layout
 *   bime aam      <smiles>  [--format text|json|svg] [--out f] — atom-atom mapping
 *                              (svg = publication-quality mapped-reaction figure:
 *                               colour mol + map numbers + bond changes + trace)
 *   bime moiety   <pathwayId> --start <met> --atoms i,j,k  — moiety trace
 *                              [--criterion strict|present]
 *                              [--format text|json]
 *   bime pathway-layout <pathwayId|--in file> [--shape auto|cycle|...] — layout
 *                              [--format text|json] [--out f]
 *   bime export   <smiles>  --format svg|mol|sdf|smiles    — file export
 *                            [--width N] [--height N] [--out f]
 *   bime pathway-list                                       — list bundled pathways
 *   bime help     [subcommand]                              — help text
 *
 * Conventions:
 *   - Input is the first positional arg, OR `--in <file>`, OR stdin (`-`).
 *   - Output goes to stdout, OR `--out <file>`.
 *   - Non-zero exit on user error or engine failure; useful from shell scripts.
 *   - No external npm deps; only Node built-ins + the editor source files.
 *
 * The CLI loads the editor modules via the same shim the test battery
 * uses (`tests/shim.js`), so any behaviour visible to the workbench is
 * also visible from the CLI.
 */
'use strict';

var fs = require('fs');
var path = require('path');

var ROOT = path.join(__dirname, '..');
var SHIM_PATH = path.join(ROOT, 'tests', 'shim.js');

// v2.4.18: the pathway / metabolic-map features ship as a SEPARATE release, so
// the public core CLI is molecule + reactions only. PATHWAY_ENABLED gates the
// pathway subcommands (moiety, pathway-layout, pathway-list) out of both the
// dispatcher and the help text from a single source of truth. The engine and the
// cmd* functions are preserved (still exported), so flipping this to true
// restores the full CLI with no other change. Mirrors js/workbench.js.
var PATHWAY_ENABLED = false;

function loadEditor() {
    // Pre-bundled binaries (bun --compile) statically load the whole editor
    // before main() runs, so the engine is already on globalThis. Skip the
    // filesystem-relative dynamic requires below (a static bundler cannot
    // follow them). Under Node/pkg the engine is not yet present, so this is
    // a no-op on first call and the normal load path runs.
    if (globalThis.RDT && globalThis.SmilesParser && globalThis.Layout) { return globalThis; }
    var shim = require(SHIM_PATH);
    shim.loadAll();
    require(path.join(ROOT, 'editor', 'RDT.js'));
    require(path.join(ROOT, 'editor', 'AtomTrace.js'));
    require(path.join(ROOT, 'editor', 'MetaboliteLibrary.js'));
    try { require(path.join(ROOT, 'editor', 'PathwayLayout.js')); } catch (e) { /* pathway engine ships separately */ }
    try { require(path.join(ROOT, 'editor', 'PathwayIO.js')); } catch (e) { /* pathway engine ships separately */ }
    require(path.join(ROOT, 'editor', 'Layout.js'));
    require(path.join(ROOT, 'editor', 'ImageExport.js'));
    require(path.join(ROOT, 'editor', 'MolfileWriter.js'));
    require(path.join(ROOT, 'editor', 'ExportStamp.js'));
    return globalThis;
}

// ---------- arg parsing -----------------------------------------------------

function parseArgs(argv) {
    // argv = process.argv.slice(2); first non-flag positional is the subcommand,
    // second non-flag positional (if any) is the input.
    var args = { _: [], flags: {} };
    var i = 0;
    while (i < argv.length) {
        var a = argv[i];
        if (a.indexOf('--') === 0) {
            var key = a.slice(2);
            var next = argv[i + 1];
            if (next && next.indexOf('--') !== 0) {
                args.flags[key] = next;
                i += 2;
            } else {
                args.flags[key] = true;
                i++;
            }
        } else {
            args._.push(a);
            i++;
        }
    }
    return args;
}

function readInput(positional, inFlag) {
    // A positional '-' means stdin (matches the documented `--in -` contract);
    // without this it was treated as a literal SMILES and silently mis-parsed.
    if (positional === '-') { return fs.readFileSync(0, 'utf8').trim(); }
    if (positional) { return positional; }
    if (inFlag) {
        if (inFlag === '-') {
            return fs.readFileSync(0, 'utf8').trim();   // stdin
        }
        try {
            return fs.readFileSync(inFlag, 'utf8').trim();
        } catch (e) {
            die('cannot read input "' + inFlag + '": ' + ((e && e.message) || String(e)));
        }
    }
    if (!process.stdin.isTTY) {
        return fs.readFileSync(0, 'utf8').trim();
    }
    return null;
}

function writeOutput(text, outFlag) {
    if (outFlag) {
        fs.writeFileSync(outFlag, text);
        return;
    }
    process.stdout.write(text);
    if (text.length && text.charAt(text.length - 1) !== '\n') { process.stdout.write('\n'); }
}

function die(message, code) {
    process.stderr.write('bime: ' + message + '\n');
    process.exit(code || 1);
}

// ---------- subcommands -----------------------------------------------------

function cmdVersion() {
    var versions = require(path.join(ROOT, 'versions.json'));
    loadEditor();
    var out = [
        'bime ' + versions.bime,
        'RDT.version         ' + (globalThis.RDT && RDT.version),
        'AtomTrace.version   ' + (globalThis.AtomTrace && AtomTrace.version),
        'ToolbarPrefs.version ' + (globalThis.ToolbarPrefs && ToolbarPrefs.version)
    ].join('\n');
    writeOutput(out);
}

function cmdParse(args) {
    loadEditor();
    var smi = readInput(args._[1], args.flags['in']);
    if (!smi) { die('parse: needs a SMILES (positional, --in, or stdin)'); }
    var rxn = SmilesParser.parse(smi);
    if (rxn.parseErrors && rxn.parseErrors.length) {
        die('parse failed: ' + rxn.parseErrors.join('; '));
    }
    var report = {
        smiles: smi,
        atomCount: rxn.atoms.length,
        bondCount: rxn.bonds.length,
        elements: {},
        charges: 0,
        stereo: 0,
        hasReactionArrow: !!rxn.reactionArrow
    };
    for (var i = 0; i < rxn.atoms.length; i++) {
        var sym = rxn.atoms[i].symbol;
        report.elements[sym] = (report.elements[sym] || 0) + 1;
        if ((rxn.atoms[i].charge || 0) !== 0) { report.charges++; }
        if (rxn.atoms[i].stereo) { report.stereo++; }
    }
    if (args.flags['format'] === 'json') {
        writeOutput(JSON.stringify(report, null, 2), args.flags['out']);
    } else {
        var lines = [
            'SMILES:        ' + report.smiles,
            'Atoms:         ' + report.atomCount,
            'Bonds:         ' + report.bondCount,
            'Charges:       ' + report.charges,
            'Stereo atoms:  ' + report.stereo,
            'Reaction:      ' + (report.hasReactionArrow ? 'yes' : 'no'),
            'Elements:'
        ];
        var keys = Object.keys(report.elements).sort();
        for (var k = 0; k < keys.length; k++) {
            lines.push('  ' + keys[k] + ': ' + report.elements[keys[k]]);
        }
        writeOutput(lines.join('\n'), args.flags['out']);
    }
}

function cmdClean(args) {
    loadEditor();
    var smi = readInput(args._[1], args.flags['in']);
    if (!smi) { die('clean: needs a SMILES'); }
    var mol = SmilesParser.parse(smi);
    if (mol.parseErrors && mol.parseErrors.length) {
        die('clean: parse failed: ' + mol.parseErrors.join('; '));
    }
    if (typeof Layout !== 'undefined' && Layout.layout) {
        for (var i = 0; i < mol.atoms.length; i++) {
            mol.atoms[i].x = (i % 6) * 0.5 - 1.25;
            mol.atoms[i].y = Math.floor(i / 6) * 0.5 - 1.25;
        }
        Layout.layout(mol);
    }
    if (args.flags['format'] === 'json') {
        var coords = mol.atoms.map(function(a) {
            return { id: a.id, symbol: a.symbol, x: a.x, y: a.y };
        });
        writeOutput(JSON.stringify({
            smiles: SmilesWriter.write(mol),
            coords: coords
        }, null, 2), args.flags['out']);
    } else {
        writeOutput(SmilesWriter.write(mol), args.flags['out']);
    }
}

// Build a mapReaction-shaped result from a reaction whose atoms ALREADY carry
// :n atom-map numbers — WITHOUT re-running the solver. Lets `aam --keep-mapping`
// render a figure (or report) from a curated / externally-supplied mapping
// instead of remapping it. The mapping is taken verbatim from the input; bond
// changes and the colour grouping are derived from it (not re-solved).
function resultFromExistingMapping(rxn) {
    var sides = RDT._splitReactionSides(rxn);
    var prodByMap = {};
    sides.products.forEach(function (m) {
        m.atoms.forEach(function (a) { if (a.mapNumber > 0) { prodByMap[a.mapNumber] = a.id; } });
    });
    var mapping = {}, mapped = 0;
    sides.reactants.forEach(function (m) {
        m.atoms.forEach(function (a) {
            if (a.mapNumber > 0 && prodByMap[a.mapNumber] != null) {
                mapping[a.id] = prodByMap[a.mapNumber]; mapped++;
            }
        });
    });
    if (mapped === 0) {
        die('aam --keep-mapping: the reaction carries no atom-map numbers ' +
            '(e.g. "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]"). ' +
            'Drop --keep-mapping to compute the mapping automatically.');
    }
    var res = {
        status: 'mapped', mapping: mapping, sides: sides, mappedCount: mapped,
        bondChanges: RDT.annotateBondChanges(sides.reactants, sides.products, mapping, {}),
        timedOut: false, warnings: [], keptMapping: true
    };
    try { res.componentPairs = RDT.deriveComponentPairs(res); } catch (e) { /* optional */ }
    try { res.confidence = RDT.deriveConfidence(res); } catch (e) { /* optional */ }
    try { var g = RDT.gradeMapping(res); res.quality = (g && g.quality) || g; } catch (e) { /* optional */ }
    return res;
}

function cmdAam(args) {
    loadEditor();
    var smi = readInput(args._[1], args.flags['in']);
    if (!smi) { die('aam: needs a reaction SMILES (A>>B form)'); }
    var rxn = SmilesParser.parse(smi);
    if (rxn.parseErrors && rxn.parseErrors.length) {
        die('aam: parse failed: ' + rxn.parseErrors.join('; '));
    }
    var opts = {};
    if (args.flags['no-chemistry-aware']) { opts.chemistryAware = false; }
    if (args.flags['timeout-ms']) { opts.timeoutMs = parseInt(args.flags['timeout-ms'], 10); }
    // --keep-mapping: use the input :n atom-map numbers as-is (no re-solve);
    // for rendering a figure from an existing/curated mapping.
    var res = args.flags['keep-mapping']
        ? resultFromExistingMapping(rxn)
        : RDT.mapReaction(rxn, opts);
    var format = args.flags['format'] || 'text';
    if (format === 'svg') {
        // v2.4.13: render the mapped reaction as a publication-quality figure —
        // colour mol + atom-atom map numbers + bond changes (red cleaved /
        // green formed / amber order) + sub-fragment trace colouring. (No
        // confidence caption.) Refine geometry first; Layout preserves the
        // reactant -> arrow -> product scheme.
        if (typeof ImageExport === 'undefined' || !ImageExport.toReactionMapSVG) {
            die('aam: SVG export (ImageExport.toReactionMapSVG) not available in this build');
        }
        // toReactionMapSVG owns the layout (refine geometry + place the arrow);
        // do NOT pre-run Layout.layout here — it would scramble the scheme.
        var rmOpts = {};
        if (args.flags['screen']) { rmOpts.publication = false; }    // default = publication; --screen for the on-screen look
        if (args.flags['transparent']) { rmOpts.background = 'transparent'; }
        if (args.flags['reaction-center']) { rmOpts.showReactionCenter = true; }
        if (args.flags['no-map-numbers']) { rmOpts.showMapNumbers = false; }
        if (args.flags['no-labels']) { rmOpts.labels = false; }   // suppress compound name/id captions
        if (args.flags['ids']) { rmOpts.showIds = true; }         // caption with KEGG id instead of name
        var _rmW = parseInt(args.flags['width'], 10);
        var _rmH = parseInt(args.flags['height'], 10);
        if (isFinite(_rmW) && _rmW > 0) { rmOpts.width = _rmW; }   // ignore non-numeric --width (fall back to auto-fit)
        if (isFinite(_rmH) && _rmH > 0) { rmOpts.height = _rmH; }
        var rmSvg = ImageExport.toReactionMapSVG(rxn, res, rmOpts);
        if (args.flags['no-stamp'] !== true && typeof ExportStamp !== 'undefined') {
            rmSvg = ExportStamp.stampSvg(rmSvg);
        }
        writeOutput(rmSvg, args.flags['out']);
    } else if (format === 'json') {
        var slim = {
            status: res.status,
            score: res.score,
            confidence: res.confidence,
            decisiveness: res.decisiveness,
            quality: res.quality,
            timedOut: !!res.timedOut,
            mappedCount: Object.keys(res.mapping || {}).length,
            mapping: res.mapping,
            bondChanges: res.bondChanges,
            reactionClass: res.reactionClass,
            classReason: res.classReason,
            diagnostic: res.diagnostic
        };
        writeOutput(JSON.stringify(slim, null, 2), args.flags['out']);
    } else {
        var lines = [
            'Status:       ' + res.status,
            'Quality:      ' + (res.quality || '-'),
            'Score:        ' + res.score,
            'Mapped atoms: ' + Object.keys(res.mapping || {}).length,
            'Bond changes: ' + (res.bondChanges ? res.bondChanges.length : 0),
            'Reaction:     ' + (res.reactionClass || '(unknown)'),
            'Confidence:   ' + (res.confidence || 0).toFixed(3),
            'Decisiveness: ' + (res.decisiveness || 0).toFixed(3),
            'Timed out:    ' + (!!res.timedOut)
        ];
        writeOutput(lines.join('\n'), args.flags['out']);
    }
}

function cmdPathwayList() {
    loadEditor();
    if (!MetaboliteLibrary || !MetaboliteLibrary.pathways) {
        die('pathway-list: MetaboliteLibrary not loaded');
    }
    var ids = Object.keys(MetaboliteLibrary.pathways);
    var lines = ['Bundled pathways:'];
    for (var i = 0; i < ids.length; i++) {
        var pw = MetaboliteLibrary.pathways[ids[i]];
        lines.push('  ' + ids[i] + '  —  ' + pw.name + ' (' + pw.steps.length + ' steps, ' + (pw.shape || 'ladder') + ')');
    }
    writeOutput(lines.join('\n'));
}

function pathwayOrderAndEdges(pw) {
    if (pw && Array.isArray(pw.nodes) && Array.isArray(pw.edges)) {
        var nodeSeen = {};
        var nodeOrder = [];
        pw.nodes.forEach(function(n) {
            var id = (n && n.id !== undefined) ? String(n.id) : '';
            if (id && !nodeSeen[id]) {
                nodeSeen[id] = true;
                nodeOrder.push(id);
            }
        });
        var nodeSet = {};
        nodeOrder.forEach(function(id) { nodeSet[id] = true; });
        var modelEdges = [];
        pw.edges.forEach(function(e) {
            if (!e || e.from === undefined || e.to === undefined) { return; }
            var from = String(e.from);
            var to = String(e.to);
            if (nodeSet[from] && nodeSet[to] && from !== to) {
                modelEdges.push({ from: from, to: to });
            }
        });
        return {
            order: nodeOrder,
            nodes: nodeOrder.map(function(id) { return { id: id }; }),
            edges: modelEdges
        };
    }
    var order = [];
    var seen = {};
    pw.steps.forEach(function(s) {
        [s.from, s.to].forEach(function(name) {
            if (!seen[name]) {
                seen[name] = true;
                order.push(name);
            }
        });
    });
    return {
        order: order,
        nodes: order.map(function(name) { return { id: name }; }),
        edges: pw.steps.map(function(s) { return { from: s.from, to: s.to }; })
    };
}

function formatNumber(n) {
    if (typeof n !== 'number' || !isFinite(n)) { return n; }
    return Math.round(n * 1000) / 1000;
}

function loadPathwayModelForLayout(args) {
    var pathwayId = args._[1];
    if (pathwayId && pathwayId !== '-') {
        var bundled = MetaboliteLibrary.pathways[pathwayId];
        if (!bundled) { die('pathway-layout: unknown pathway "' + pathwayId + '"'); }
        return { id: pathwayId, name: bundled.name, model: bundled, defaultShape: bundled.shape || 'auto' };
    }

    var raw = (pathwayId === '-') ? fs.readFileSync(0, 'utf8').trim() : readInput(null, args.flags['in']);
    if (!raw) { die('pathway-layout: needs a pathway id or --in <bime-pathway.json>'); }
    var model;
    try {
        if (PathwayIO && PathwayIO.parseJSON) {
            model = PathwayIO.parseJSON(raw);
        } else {
            model = JSON.parse(raw);
        }
    } catch (e) {
        try {
            model = JSON.parse(raw);
        } catch (jsonError) {
            die('pathway-layout: cannot parse input JSON: ' + ((jsonError && jsonError.message) || String(jsonError)));
        }
    }
    return {
        id: model.name || 'custom',
        name: model.name || 'Custom pathway',
        model: model,
        defaultShape: model.layoutShape || 'auto'
    };
}

function cmdPathwayLayout(args) {
    loadEditor();
    if (!MetaboliteLibrary || !MetaboliteLibrary.pathways) {
        die('pathway-layout: MetaboliteLibrary not loaded');
    }
    if (!PathwayLayout || !PathwayLayout.compute) {
        die('pathway-layout: PathwayLayout not loaded');
    }
    var source = loadPathwayModelForLayout(args);
    var pw = source.model;

    var requested = args.flags['shape'] || args.flags['layout'] || args.flags['layout-shape'] || source.defaultShape;
    var shape = PathwayLayout.normalizeLayoutShape(requested);
    var graph = pathwayOrderAndEdges(pw);
    if (!graph.nodes.length) { die('pathway-layout: input has no layout nodes'); }
    var laid = PathwayLayout.compute(graph.nodes, graph.edges, {
        layoutShape: shape,
        orientation: args.flags['orientation'] || 'horizontal',
        order: graph.order,
        hub: pw.hub || null,
        trunk: pw.trunk || null
    });
    var out = {
        pathway: source.id,
        name: source.name,
        requestedShape: shape,
        resolvedShape: laid.stats && laid.stats.layoutShape || shape,
        componentShapes: laid.stats && laid.stats.layoutShapes || [],
        nodeCount: graph.nodes.length,
        edgeCount: graph.edges.length,
        componentCount: laid.componentCount || 0,
        layerCount: laid.layerCount || 0,
        stats: laid.stats || {},
        nodes: graph.order.map(function(name) {
            var p = laid.positions[name] || { x: 0, y: 0 };
            return { id: name, x: formatNumber(p.x), y: formatNumber(p.y) };
        }),
        edges: graph.edges
    };

    if ((args.flags['format'] || 'text') === 'json') {
        writeOutput(JSON.stringify(out, null, 2), args.flags['out']);
        return;
    }

    var lines = [
        'Pathway:     ' + out.name + ' (' + out.pathway + ')',
        'Shape:       ' + out.requestedShape,
        'Resolved:    ' + out.resolvedShape + (out.componentShapes.length ? ' [' + out.componentShapes.join(', ') + ']' : ''),
        'Nodes:       ' + out.nodeCount,
        'Edges:       ' + out.edgeCount,
        'Components:  ' + out.componentCount,
        'Layer count: ' + out.layerCount,
        '',
        'Positions:'
    ];
    for (var i = 0; i < out.nodes.length; i++) {
        var n = out.nodes[i];
        lines.push('  ' + (n.id + '                              ').slice(0, 30) +
            ' x=' + n.x + ' y=' + n.y);
    }
    writeOutput(lines.join('\n'), args.flags['out']);
}

function cmdMoiety(args) {
    loadEditor();
    var pathwayId = args._[1];
    if (!pathwayId) { die('moiety: needs a pathway id (try: bime pathway-list)'); }
    var pw = MetaboliteLibrary.pathways[pathwayId];
    if (!pw) { die('moiety: unknown pathway "' + pathwayId + '"'); }
    var startName = args.flags['start'] || pw.steps[0].from;
    var atomsArg = args.flags['atoms'];
    if (!atomsArg) { die('moiety: needs --atoms <comma-separated-indices> (atom indices in the start metabolite)'); }
    var moietyAtoms = String(atomsArg).split(',').map(function(s) { return parseInt(s.trim(), 10); }).filter(function(n) { return !isNaN(n); });
    if (!moietyAtoms.length) { die('moiety: --atoms is empty'); }

    // Build the metabolite order + DAG edges (same as the workbench does).
    var order = [], seen = {};
    pw.steps.forEach(function(s) { [s.from, s.to].forEach(function(n) { if (!seen[n]) { seen[n] = true; order.push(n); } }); });
    var nameToIdx = {};
    order.forEach(function(n, i) { nameToIdx[n.toLowerCase()] = i; });
    var startIdx = nameToIdx[startName.toLowerCase()];
    if (startIdx === undefined) {
        die('moiety: start metabolite "' + startName + '" not found in pathway "' + pathwayId + '"');
    }
    var smilesAll = order.map(function(n) { return (MetaboliteLibrary.find(n) || {}).smiles || ''; });
    var edges = pw.steps.map(function(s) {
        return { from: nameToIdx[s.from.toLowerCase()], to: nameToIdx[s.to.toLowerCase()] };
    });
    // Slice from startIdx onwards (re-index edges 0-based).
    var slice = smilesAll.slice(startIdx);
    var sliceEdges = [];
    edges.forEach(function(e) {
        if (e.from >= startIdx && e.to >= startIdx) {
            sliceEdges.push({ from: e.from - startIdx, to: e.to - startIdx });
        }
    });
    var tr = AtomTrace.tracePathway(slice, { edges: sliceEdges });
    if (!tr) { die('moiety: tracePathway returned null'); }
    var criterion = args.flags['criterion'] || 'strict';
    var mt = AtomTrace.deriveMoietyTrace(tr, moietyAtoms, {
        metaboliteSmiles: slice, criterion: criterion
    });
    if (mt.error) { die('moiety: ' + mt.error); }

    var format = args.flags['format'] || 'text';
    if (format === 'json') {
        writeOutput(JSON.stringify({
            pathway: pathwayId,
            start: order[startIdx],
            moietyAtoms: mt.moietyAtoms,
            criterion: mt.criterion,
            breakAt: (mt.breakAt !== null) ? order[startIdx + mt.breakAt] : null,
            survivors: mt.survivors.map(function(i) { return order[startIdx + i]; }),
            perMetabolite: mt.perMetabolite.map(function(p) {
                return {
                    name: order[startIdx + p.metabolite],
                    status: p.status,
                    alive: p.alive,
                    present: p.present.length,
                    missing: p.missing.length,
                    preservedBonds: p.preservedBonds.length,
                    brokenBonds: p.brokenBonds.length
                };
            })
        }, null, 2), args.flags['out']);
    } else {
        var lines = [
            'Pathway:   ' + pw.name + ' (' + pathwayId + ')',
            'Start:     ' + order[startIdx],
            'Moiety:    ' + mt.moietyAtoms.join(',') + ' (' + mt.moietySize + ' atoms)',
            'Criterion: ' + mt.criterion,
            'Breaks at: ' + (mt.breakAt !== null ? order[startIdx + mt.breakAt] : '(never)'),
            'Survivors: ' + mt.survivors.length + ' / ' + mt.perMetabolite.length,
            '',
            'Per metabolite:'
        ];
        for (var i = 0; i < mt.perMetabolite.length; i++) {
            var p = mt.perMetabolite[i];
            lines.push('  ' + (order[startIdx + i] + '                          ').slice(0, 30) +
                ' status=' + p.status +
                ' atoms=' + p.present.length + '/' + mt.moietySize +
                ' bonds=' + p.preservedBonds.length + '/' + mt.startBonds.length);
        }
        writeOutput(lines.join('\n'), args.flags['out']);
    }
}

function cmdExport(args) {
    loadEditor();
    var smi = readInput(args._[1], args.flags['in']);
    if (!smi) { die('export: needs a SMILES'); }
    var format = (args.flags['format'] || 'smiles').toLowerCase();
    var mol = SmilesParser.parse(smi);
    if (mol.parseErrors && mol.parseErrors.length) {
        die('export: parse failed: ' + mol.parseErrors.join('; '));
    }
    if (typeof Layout !== 'undefined' && Layout.layout && mol.atoms.length) {
        for (var i = 0; i < mol.atoms.length; i++) {
            mol.atoms[i].x = (i % 6) * 0.5 - 1.25;
            mol.atoms[i].y = Math.floor(i / 6) * 0.5 - 1.25;
        }
        Layout.layout(mol);
    }
    // v2.0.65: default-on provenance stamping. The CLI runs every
    // export through ExportStamp before emit so the artefact carries
    // BIME version + SHA-384 fingerprint + Apache-2.0 attribution.
    // Users can opt out with `--no-stamp` for downstream tools that
    // refuse to parse comment-prefixed input.
    var stamping = (args.flags['no-stamp'] !== true);
    var out;
    if (format === 'smiles') {
        out = SmilesWriter.write(mol);
    } else if (format === 'svg') {
        if (typeof ImageExport === 'undefined' || !ImageExport.toSVG) {
            die('export: ImageExport.toSVG not available in this build');
        }
        out = ImageExport.toSVG(mol, {
            width:  parseInt(args.flags['width'],  10) || 320,
            height: parseInt(args.flags['height'], 10) || 240,
            padding: 12,
            background: 'transparent',
            showAtomLabels: true
        });
    } else if (format === 'mol') {
        if (typeof MolfileWriter === 'undefined' || !MolfileWriter.write) {
            die('export: MolfileWriter.write not available in this build');
        }
        out = MolfileWriter.write(mol);
    } else if (format === 'sdf') {
        if (typeof MolfileWriter === 'undefined' || !MolfileWriter.writeSDF) {
            die('export: MolfileWriter.writeSDF not available in this build');
        }
        // v2.0.64: writeSDF takes a single Molecule (not an array) — fixed
        // from the v2.0.62 first cut which guessed at the signature.
        out = MolfileWriter.writeSDF(mol);
    } else if (format === 'rxn') {
        if (typeof MolfileWriter === 'undefined' || !MolfileWriter.writeRXN) {
            die('export: MolfileWriter.writeRXN not available in this build');
        }
        if (!mol.reactionArrow) {
            die('export: --format rxn requires a reaction SMILES (A>>B)');
        }
        out = MolfileWriter.writeRXN(mol);
    } else if (format === 'png' || format === 'pdf') {
        die('export: ' + format + ' needs a browser DOM (PNG/PDF go through Canvas); use --format svg/mol/sdf/smiles from the CLI');
    } else {
        die('export: unknown --format "' + format + '" (use svg, mol, sdf, smiles)');
    }
    // v2.0.65: provenance stamp (default on, opt out with --no-stamp).
    if (stamping && typeof ExportStamp !== 'undefined') {
        if (format === 'svg') {
            out = ExportStamp.stampSvg(out);
        } else {
            out = ExportStamp.stampText(format, out);
        }
    }
    writeOutput(out, args.flags['out']);
}

// ---------- help text -------------------------------------------------------

var HELP_TOP = [
    'bime — command-line wrapper around the BIME engine.',
    '',
    'Usage:  bime <command> [args] [flags]',
    '',
    'Commands:',
    '  version                          Print bundle + module versions.',
    '  parse <smiles>                   SMILES diagnostic (atoms, bonds, elements).',
    '  clean <smiles>                   Run Layout.layout; print canonical SMILES.',
    '  aam   <reactionSmiles>           Atom-atom mapping (RDT.mapReaction).'
].concat(PATHWAY_ENABLED ? [
    '  moiety <pathwayId> --start <met> --atoms i,j,k',
    '                                   Trace a connected-subgraph moiety through a',
    '                                   bundled pathway. See `pathway-list` for ids.',
    '  pathway-layout <pathwayId|--in>  Compute a reusable pathway layout.'
] : []).concat([
    '  export <smiles> --format <fmt>   Export to svg | mol | sdf | smiles.'
]).concat(PATHWAY_ENABLED ? [
    '  pathway-list                     List the bundled pathway templates.'
] : []).concat([
    '  help [command]                   Show this help, or per-command help.',
    '',
    'Global flags:',
    '  --in   <file|-> | --out <file>   I/O paths (default: stdin/stdout).',
    '  --format <text|json>             Output format where applicable.',
    '',
    'Run `bime help <command>` for per-command flags.',
    ''
]).join('\n');

var HELP_BY_CMD = {
    aam: 'bime aam <reactionSmiles>\n' +
         '  Flags:\n' +
         '    --format text|json|svg   (default: text)\n' +
         '    --no-chemistry-aware Disable the v2.0.57 aromatic-ring + stereo penalty.\n' +
         '    --timeout-ms N       Per-call timeout for the AAM solver (default: 2000).\n' +
         '    --keep-mapping       Use the input :n atom-map numbers as-is (no re-solve);\n' +
         '                         render/report a curated mapping instead of remapping it.\n' +
         '  SVG (v2.4.13) — a publication-quality mapped-reaction figure: colour mol\n' +
         '  + atom-atom map numbers + bond changes (red cleaved / green formed /\n' +
         '  amber order) + sub-fragment trace colouring. Auto-sized; white background.\n' +
         '    --screen             Use the on-screen look instead of the print preset.\n' +
         '    --transparent        Transparent background (drop onto any page).\n' +
         '    --reaction-center    Add dashed reaction-centre rings.\n' +
         '    --no-map-numbers     Omit the atom-atom map numbers.\n' +
         '    --ids                Caption components with the KEGG id, not the name.\n' +
         '    --no-labels          Omit the per-component compound-name captions.\n' +
         '    --width N --height N Pin the canvas size (default: auto-fit to content).\n' +
         '    --out file           Write the SVG to a file (default: stdout).\n' +
         '    --no-stamp           Opt out of the BIME provenance stamp.\n',
    moiety: 'bime moiety <pathwayId> --start <metabolite> --atoms i,j,k\n' +
            '  Flags:\n' +
            '    --start <name>          Start metabolite (default: pathway step[0].from).\n' +
            '    --atoms i,j,k          Comma-separated atom indices in the start metabolite.\n' +
            '    --criterion strict|present  Strict (default) requires bond preservation.\n' +
            '    --format text|json\n',
    'pathway-layout': 'bime pathway-layout <pathwayId|--in file> [--shape auto|hierarchical|ladder|cycle|branched|fanout|fanin|split-merge|hybrid|loop-iterative]\n' +
                      '  Flags:\n' +
                      '    --shape <name>          Preferred layout shape (default: pathway default).\n' +
                      '    --layout <name>         Alias for --shape.\n' +
                      '    --in file               Layout a saved BIME pathway JSON file.\n' +
                      '    --format text|json      Structured JSON includes nodes, edges, stats.\n' +
                      '    --out file              Write the report to a file.\n',
    "export": 'bime export <smiles> --format svg|mol|sdf|rxn|smiles\n' +
              '  Flags:\n' +
              '    --format svg|mol|sdf|rxn|smiles\n' +
              '    --width N --height N     SVG canvas size (svg only).\n' +
              '    --out file               Write to file (default: stdout).\n' +
              '    --no-stamp               v2.0.65: opt out of the BIME provenance\n' +
              '                              stamp (default ON — embeds BIME version,\n' +
              '                              SHA-384 fingerprint, Apache-2.0\n' +
              '                              attribution into every artefact).\n',
    parse: 'bime parse <smiles> [--format text|json]\n',
    clean: 'bime clean <smiles> [--format text|json]\n',
    'pathway-list': 'bime pathway-list — list ids + names of the bundled pathway templates.\n',
    version: 'bime version — print BIME + RDT + AtomTrace + ToolbarPrefs versions.\n'
};

// v2.4.18: drop the gated pathway subcommands from per-command help too, so
// `bime help moiety` in the core build falls back to the top-level help.
if (!PATHWAY_ENABLED) {
    delete HELP_BY_CMD.moiety;
    delete HELP_BY_CMD['pathway-layout'];
    delete HELP_BY_CMD['pathway-list'];
}

function cmdHelp(args) {
    var sub = args._[1];
    if (sub && HELP_BY_CMD[sub]) {
        writeOutput(HELP_BY_CMD[sub]);
        return;
    }
    writeOutput(HELP_TOP);
}

// ---------- main ------------------------------------------------------------

function main() {
    var args = parseArgs(process.argv.slice(2));
    var cmd = args._[0];
    if (!cmd || args.flags['help'] === true) {
        cmdHelp(args);
        process.exit(cmd ? 0 : (args.flags['help'] ? 0 : 1));
        return;
    }
    switch (cmd) {
        case 'version':      return cmdVersion();
        case 'parse':        return cmdParse(args);
        case 'clean':        return cmdClean(args);
        case 'aam':          return cmdAam(args);
        case 'moiety':       return PATHWAY_ENABLED ? cmdMoiety(args)
                                 : die('unknown command "' + cmd + '". Run `bime help` for a list.');
        case 'pathway-layout': return PATHWAY_ENABLED ? cmdPathwayLayout(args)
                                 : die('unknown command "' + cmd + '". Run `bime help` for a list.');
        case 'export':       return cmdExport(args);
        case 'pathway-list': return PATHWAY_ENABLED ? cmdPathwayList()
                                 : die('unknown command "' + cmd + '". Run `bime help` for a list.');
        case 'help':         return cmdHelp(args);
        default:
            die('unknown command "' + cmd + '". Run `bime help` for a list.');
    }
}

if (require.main === module) {
    try { main(); }
    catch (e) { die((e && e.message) || String(e)); }
}

module.exports = {
    parseArgs: parseArgs,
    loadEditor: loadEditor,
    cmdVersion: cmdVersion,
    cmdParse: cmdParse,
    cmdClean: cmdClean,
    cmdAam: cmdAam,
    cmdMoiety: cmdMoiety,
    cmdPathwayLayout: cmdPathwayLayout,
    cmdExport: cmdExport,
    cmdPathwayList: cmdPathwayList,
    cmdHelp: cmdHelp,
    main: main
};
