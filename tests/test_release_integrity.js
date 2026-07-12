/**
 * tests/test_release_integrity.js
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Guards the release wiring that tends to fragment: browser source script
 * order, export string safety, and SDF metadata delimiters.
 */
'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var shim = require(path.join(__dirname, 'shim.js'));

shim.loadAll();
shim.require_editor('Renderer');
shim.require_editor('MolfileWriter');
shim.require_editor('ImageExport');
shim.require_editor('Templates');
shim.require_editor('MolEditor');

var runner = shim.makeRunner('Release integrity and export hardening');
var test = runner.test;

var ROOT = path.join(__dirname, '..');
var EDITOR_DIR = path.join(ROOT, 'editor');
var editorFiles = require(path.join(ROOT, 'tools', 'editor-files.js'));
var Molecule = globalThis.Molecule;
var ImageExport = globalThis.ImageExport;
var MolfileWriter = globalThis.MolfileWriter;
var Templates = globalThis.Templates;
var Renderer = globalThis.Renderer;
var MolEditor = globalThis.MolEditor;

function editorScriptsForPage(page) {
    var html = fs.readFileSync(path.join(ROOT, page), 'utf8');
    var scripts = [];
    html.replace(/<script\s+src="editor\/([^"?]+\.js)(?:\?v=[^"]*)?"/g,
        function (_m, file) {
            scripts.push(file);
            return _m;
        });
    return scripts;
}

function readWorkbenchHtml() {
    return fs.readFileSync(path.join(ROOT, 'workbench.html'), 'utf8');
}

function readWorkbenchJs() {
    return fs.readFileSync(path.join(ROOT, 'js', 'workbench.js'), 'utf8');
}

function readText(file) {
    // Tolerate files that are intentionally absent from the public
    // distribution (the public sync drops internal-only docs). A missing
    // file has no URL to validate, so an empty body is the correct
    // no-op for the URL-hygiene checks below.
    try { return fs.readFileSync(path.join(ROOT, file), 'utf8'); }
    catch (e) { return ''; }
}

function distance(a, b) {
    var dx = a.x - b.x;
    var dy = a.y - b.y;
    return Math.sqrt(dx * dx + dy * dy);
}

function bondKey(a, b) {
    return a < b ? a + ':' + b : b + ':' + a;
}

function templateGeometryHardIssues(name, tmpl) {
    var issues = [];
    var atoms = tmpl && tmpl.atoms ? tmpl.atoms : [];
    var bonds = tmpl && tmpl.bonds ? tmpl.bonds : [];
    var bondLength = Molecule.BOND_LENGTH || 30;
    var minBond = bondLength * 0.45;
    var maxBond = bondLength * 2.2;
    var severeClash = bondLength * 0.33;
    var bonded = {};
    var seenBonds = {};

    if (!atoms.length) {
        issues.push('no atoms');
    }

    for (var i = 0; i < atoms.length; i++) {
        var atom = atoms[i];
        if (!atom || typeof atom.symbol !== 'string' || !atom.symbol) {
            issues.push('atom ' + i + ' missing symbol');
        }
        if (!atom || !isFinite(atom.x) || !isFinite(atom.y)) {
            issues.push('atom ' + i + ' has non-finite coordinate');
        }
    }

    for (var b = 0; b < bonds.length; b++) {
        var bond = bonds[b];
        if (!bond || bond.a1 < 0 || bond.a2 < 0 || bond.a1 >= atoms.length || bond.a2 >= atoms.length) {
            issues.push('bond ' + b + ' has invalid atom index');
            continue;
        }
        if (bond.a1 === bond.a2) {
            issues.push('bond ' + b + ' is a self-bond');
            continue;
        }

        var key = bondKey(bond.a1, bond.a2);
        if (seenBonds[key]) {
            issues.push('bond ' + b + ' duplicates ' + key);
        }
        seenBonds[key] = true;
        bonded[key] = true;

        var len = distance(atoms[bond.a1], atoms[bond.a2]);
        if (len < minBond || len > maxBond) {
            issues.push(name + ' bond ' + b + ' length ' + len.toFixed(2));
        }
    }

    for (var j = 0; j < atoms.length; j++) {
        for (var k = j + 1; k < atoms.length; k++) {
            var d = distance(atoms[j], atoms[k]);
            if (d < 0.25) {
                issues.push('atoms ' + j + '/' + k + ' duplicate coordinates');
            } else if (!bonded[bondKey(j, k)] && d < severeClash) {
                issues.push('atoms ' + j + '/' + k + ' severe clash ' + d.toFixed(2));
            }
        }
    }

    return issues;
}

function templateSignature(tmpl) {
    var atoms = tmpl.atoms.map(function(atom) {
        return [
            atom.symbol,
            Math.round(atom.x * 1000) / 1000,
            Math.round(atom.y * 1000) / 1000
        ].join(':');
    });
    var bonds = tmpl.bonds.map(function(bond) {
        return [
            Math.min(bond.a1, bond.a2),
            Math.max(bond.a1, bond.a2),
            bond.type
        ].join(':');
    }).sort();
    return atoms.join('|') + '//' + bonds.join('|') + '//' + JSON.stringify(tmpl.aromatic || []);
}

test('canonical editor manifest has no duplicates and every file exists', function() {
    var seen = {};
    editorFiles.FILES.forEach(function(file) {
        assert.ok(!seen[file], 'duplicate editor file: ' + file);
        seen[file] = true;
        assert.ok(fs.existsSync(path.join(EDITOR_DIR, file)), 'missing editor file: ' + file);
    });
});

test('source pages are owned by the canonical editor manifest', function() {
    ['index.html', 'examples.html', 'test.html'].forEach(function(page) {
        var html = fs.readFileSync(path.join(ROOT, page), 'utf8');
        assert.ok(html.indexOf('<!-- BIME_SOURCE_SCRIPTS:start -->') >= 0, page + ' missing source marker start');
        assert.ok(html.indexOf('<!-- BIME_SOURCE_SCRIPTS:end -->') >= 0, page + ' missing source marker end');
        assert.deepStrictEqual(editorScriptsForPage(page), editorFiles.FILES, page + ' editor scripts drifted');
    });
});

test('public release surfaces point to the public distribution repository', function() {
    var publicFiles = [
        'README.md',
        'USER_GUIDE.md',
        'CHANGES.md',
        'CONTRIBUTING.md',
        'SUPPORT.md',
        'CITATION.cff',
        'LICENSE.txt',
        'index.html',
        'workbench.html',
        'examples.html',
        'docs.html',
        'screenshots.html',
        'test.html',
        'user-guide.html',
        '.github/ISSUE_TEMPLATE/bug_report.md',
        '.github/ISSUE_TEMPLATE/config.yml',
        '.github/REPO_DESCRIPTION.md',
        'docs/EDUCATORS.md',
        'docs/EMBED.md',
        'docs/HOSTING.md',
        'docs/RESEARCHERS.md',
        'docs/STUDENTS.md',
        'docs/USAGE.md'
    ];
    var privatePublicUrl = /(?:https?:\/\/)?(?:www\.)?(?:github\.com\/asad\/bime|asad\.github\.io\/bime)(?!-dist)(?:\.git|\/|$|[\s"'<>)])/;
    var malformedDistUrl = /\bbime-dist(?:-dist)+\b/;
    var failures = [];

    publicFiles.forEach(function(file) {
        var body = readText(file);
        if (privatePublicUrl.test(body)) {
            failures.push(file + ': private repository or Pages URL');
        }
        if (malformedDistUrl.test(body)) {
            failures.push(file + ': malformed bime-dist URL');
        }
    });

    assert.deepStrictEqual(failures, []);
});

test('SVG export escapes metadata and rejects unsafe background paint', function() {
    var mol = new Molecule();
    mol.addAtom('C', 0, 0);
    mol.name = 'bad\'"><script>alert(1)</script>';

    var svg = ImageExport.toSVG(mol, {
        background: '" onload="alert(1)',
        width: 120,
        height: 90
    });

    assert.strictEqual(svg.indexOf('<script>'), -1);
    assert.strictEqual(svg.indexOf('onload='), -1);
    assert.ok(svg.indexOf('&lt;script&gt;alert(1)&lt;/script&gt;') >= 0);
    assert.ok(svg.indexOf('bad&#39;&quot;&gt;') >= 0);
    assert.ok(svg.indexOf('fill="#ffffff"') >= 0);
});

test('SDF metadata strips record delimiters and normalises field newlines', function() {
    var mol = new Molecule();
    mol.addAtom('C', 0, 0);
    mol.name = 'A\n$$$$\nB';
    mol.comment = 'C\r\nD';

    var sdf = MolfileWriter.writeSDF(mol);
    var delimiterLines = sdf.split(/\r?\n/).filter(function(line) {
        return line === '$$$$';
    });

    assert.strictEqual(delimiterLines.length, 1);
    assert.ok(sdf.indexOf('> <NAME>\nA B\n') >= 0);
    assert.ok(sdf.indexOf('> <COMMENT>\nC D\n') >= 0);
});

test('live renderer honours depiction-only wedge and dash stereo', function() {
    var mol = new Molecule();
    var substituent = mol.addAtom('Cl', 30, 0);
    var centre = mol.addAtom('C', 0, 0);
    var bond = mol.addBond(substituent.id, centre.id, Molecule.BOND_SINGLE);
    bond.depictStereo = Molecule.STEREO_WEDGE;
    bond.depictStereoFromAtom = centre.id;

    var renderer = new Renderer(document.createElement('div'), 120, 80);
    renderer.setMolecule(mol);
    var captured = null;
    renderer._drawWedgeBond = function(x1, y1, x2, y2) {
        captured = { x1: x1, y1: y1, x2: x2, y2: y2 };
    };
    renderer._drawDashBond = function() {
        throw new Error('dash renderer called for wedge bond');
    };

    renderer._drawBond(bond);
    assert.ok(captured, 'wedge renderer was not called for depictStereo');
    assert.ok(captured.x1 < captured.x2, 'wedge should be drawn from recorded stereo centre');

    bond.depictStereo = Molecule.STEREO_DASH;
    var dashCalled = false;
    renderer._drawWedgeBond = function() {
        throw new Error('wedge renderer called for dash bond');
    };
    renderer._drawDashBond = function() { dashCalled = true; };
    renderer._drawBond(bond);
    assert.ok(dashCalled, 'dash renderer was not called for depictStereo');
});

test('local molecular file import rejects oversized and unsupported files', function() {
    assert.strictEqual(MolEditor.validateImportFile({
        name: 'library.sdf',
        size: 1024,
        type: 'chemical/x-mdl-sdfile'
    }), '');
    assert.ok(/too large/.test(MolEditor.validateImportFile({
        name: 'huge.sdf',
        size: MolEditor.MAX_IMPORT_BYTES + 1,
        type: 'chemical/x-mdl-sdfile'
    })));
    assert.ok(/Unsupported molecular file type/.test(MolEditor.validateImportFile({
        name: 'binary.cdx',
        size: 128,
        type: 'application/octet-stream'
    })));
});

test('adamantane template has the canonical 10-atom 12-bond cage graph', function() {
    var tmpl = Templates.adamantane();
    assert.strictEqual(tmpl.atoms.length, 10);
    assert.strictEqual(tmpl.bonds.length, 12);

    var degrees = tmpl.atoms.map(function() { return 0; });
    tmpl.bonds.forEach(function(bond) {
        degrees[bond.a1]++;
        degrees[bond.a2]++;
    });
    degrees.sort(function(a, b) { return a - b; });
    assert.deepStrictEqual(degrees, [2, 2, 2, 2, 2, 2, 3, 3, 3, 3]);
});

test('stored scaffold templates pass hard geometry checks and SVG rendering', function() {
    var names = Templates.list();
    var failures = [];
    assert.ok(names.length >= 50, 'unexpectedly small scaffold/template list');

    names.forEach(function(name) {
        var tmpl = Templates[name]();
        var issues = templateGeometryHardIssues(name, tmpl);
        if (issues.length) {
            failures.push(name + ': ' + issues.slice(0, 4).join('; '));
        }

        var mol = new Molecule();
        var atomIds = Templates.apply(mol, tmpl, 0, 0);
        if (atomIds.length !== tmpl.atoms.length) {
            failures.push(name + ': apply returned ' + atomIds.length + ' atoms for ' + tmpl.atoms.length + ' template atoms');
            return;
        }

        var svg;
        try {
            svg = ImageExport.toSVG(mol, {
                background: '#ffffff',
                width: 320,
                height: 220
            });
        } catch (err) {
            failures.push(name + ': SVG export threw ' + err.message);
            return;
        }

        if (svg.indexOf('<svg') < 0) {
            failures.push(name + ': SVG export did not return SVG markup');
        }
    });

    assert.deepStrictEqual(failures, []);
});

test('public scaffold template list hides exact duplicate aliases', function() {
    var names = Templates.list();
    var seenNames = {};
    var seenSignatures = {};
    var duplicates = [];

    names.forEach(function(name) {
        assert.strictEqual(seenNames[name], undefined, 'duplicate template name ' + name);
        seenNames[name] = true;
        var sig = templateSignature(Templates[name]());
        if (seenSignatures[sig]) {
            duplicates.push(seenSignatures[sig] + ' / ' + name);
        } else {
            seenSignatures[sig] = name;
        }
    });

    assert.deepStrictEqual(duplicates, []);

    var aliases = Templates.aliases || {};
    Object.keys(aliases).forEach(function(alias) {
        assert.strictEqual(names.indexOf(alias), -1, alias + ' should be hidden from the public drawer list');
        assert.strictEqual(typeof Templates[alias], 'function', alias + ' alias should remain callable');
        assert.strictEqual(typeof Templates[aliases[alias]], 'function', aliases[alias] + ' canonical template should exist');
    });
});

test('workbench setSafeSvg sanitizes root and child SVG markup', function() {
    var js = readWorkbenchJs();
    var start = js.indexOf('function setSafeSvg');
    var end = js.indexOf('/* ---- v1.8.1: structure-preview tooltip', start);
    assert.ok(start >= 0 && end > start, 'setSafeSvg function not found in workbench.html');
    var fn = js.slice(start, end);
    assert.ok(fn.indexOf('[root].concat(Array.prototype.slice.call(root.querySelectorAll') >= 0,
        'root <svg> node must be included in attribute sanitization');
    assert.ok(/script,foreignObject,iframe,object,embed/.test(fn),
        'active SVG container elements must be removed');
    assert.ok(/attrName\.indexOf\('on'\)\s*===\s*0/.test(fn),
        'event-handler attributes must be stripped');
    assert.ok(/javascript:/.test(fn),
        'javascript: URLs must be stripped');
});

test('workbench compound lookup and reaction labels are guarded', function() {
    var html = readWorkbenchHtml();
    var js = readWorkbenchJs();
    assert.ok(/connect-src 'self' https:\/\/pubchem\.ncbi\.nlm\.nih\.gov/.test(html),
        'CSP must scope optional network lookup to the configured compound host only');
    assert.ok(html.indexOf('id="bime-lookup-input"') >= 0,
        'compound lookup input should exist');
    assert.ok(html.indexOf('Name, CAS, IUPAC, formula, CID, or InChI') >= 0 &&
            html.indexOf('Name / CAS / IUPAC') >= 0,
        'compound lookup UI should advertise name, CAS, IUPAC, formula, CID, and InChI modes');
    assert.ok(js.indexOf('function compoundLookupUrl') >= 0,
        'compound lookup URL builder should exist');
    assert.ok(/encodeURIComponent\(q\)/.test(js) && /encodeURIComponent\(formula\)/.test(js),
        'lookup URL builder must encode user query strings');
    assert.ok(js.indexOf('fastformula/') >= 0 && js.indexOf('/cids/JSON?MaxRecords=10') >= 0,
        'formula lookup should resolve CIDs before requesting properties');
    assert.ok(js.indexOf('inchikey/') >= 0 && js.indexOf('inchi/cids/JSON?inchi=') >= 0,
        'InChIKey and full InChI lookup should be supported');
    assert.ok(js.indexOf('function compoundPropertyUrlForCids') >= 0,
        'CID property URL builder should exist for formula/InChI lookup');
    var blockedProductHosts = [
        ['free', 'chem', 'draw.com'].join(''),
        ['rev', 'vity', 'signals.com'].join('')
    ];
    blockedProductHosts.forEach(function (host) {
        assert.ok(html.indexOf(host) < 0 && js.indexOf(host) < 0,
            'workbench must not call external product pages for lookup');
    });
    assert.ok(js.indexOf('editor.readGenericMolecularInput(hit.smiles)') >= 0,
        'lookup hits must load through the normal molecular input parser');
    assert.ok(html.indexOf('id="wb-reaction-conditions"') >= 0 &&
            html.indexOf('id="wb-reaction-yield"') >= 0 &&
            html.indexOf('id="wb-reaction-note"') >= 0,
        'reaction condition, yield, and note controls should exist');
    assert.ok(js.indexOf('function sanitizeReactionLabelPart') >= 0 &&
            js.indexOf('function composeReactionLabel') >= 0,
        'reaction labels should be sanitised and composed in one path');
    assert.ok(js.indexOf('arrow.conditions = composeReactionLabel(clean)') >= 0,
        'reaction labels should render through the existing arrow conditions field');
});

test('workbench keeps scripts external and actions declarative', function() {
    var html = readWorkbenchHtml();
    var js = readWorkbenchJs();
    var uncommented = html.replace(/<!--[\s\S]*?-->/g, '');
    assert.ok(!/<script(?![^>]*\bsrc=)[^>]*>/.test(uncommented),
        'workbench should not use inline script blocks');
    assert.ok(html.indexOf('onclick=') < 0, 'workbench controls should not use inline onclick handlers');
    assert.ok(/script-src 'self'/.test(html) && !/script-src 'self' 'unsafe-inline'/.test(html),
        'workbench CSP should not require inline scripts');
    assert.ok(html.indexOf('css/workbench.css') >= 0, 'workbench stylesheet should be external');
    assert.ok(html.indexOf('js/workbench.js') >= 0, 'workbench controller should be external');
    assert.ok(js.indexOf('function bindWorkbenchActions') >= 0 &&
            js.indexOf("closest('[data-wb-action]')") >= 0,
        'workbench action dispatch should be delegated through data-wb-action');
    assert.ok(html.indexOf('wb-drawer-strip') >= 0 &&
            html.indexOf('data-wb-action="toggle-workbench-section"') >= 0 &&
            html.indexOf('id="bime-search-section"') >= 0 &&
            /id="bime-search-section"[^>]*hidden/.test(html),
        'secondary lookup/search/import surfaces should be collapsed behind explicit drawers');
    assert.ok(html.indexOf('id="wb-pathway-section"') >= 0 &&
            html.indexOf('id="wb-pathway-canvas"') >= 0 &&
            /id="wb-pathway-section"[^>]*hidden/.test(html) &&
            html.indexOf('data-wb-action="set-workbench-mode"') >= 0 &&
            html.indexOf('data-mode="pathway"') >= 0,
        'pathway canvas should be collapsed by default and available through the editor mode toggle');
    assert.ok(js.indexOf('function initPathwayCanvas') >= 0 &&
            js.indexOf('function renderPathwayCanvas') >= 0 &&
            js.indexOf('function exportPathwaySvg') >= 0 &&
            js.indexOf('function renderMechanismArrow') >= 0 &&
            js.indexOf('function layoutPathwayCanvas') >= 0,
        'pathway canvas should have centralized SVG rendering and export wiring');
    var initWorkbenchBlock = js.slice(js.indexOf('function initWorkbenchUX'), js.indexOf('function initWorkbenchDisclosures'));
    var openedBlock = js.slice(js.indexOf('function onWorkbenchSectionOpened'), js.indexOf('function bindWorkbenchActions'));
    assert.ok(initWorkbenchBlock.indexOf('initPathwayCanvas()') < 0 &&
            openedBlock.indexOf("id === 'wb-pathway-section'") >= 0 &&
            openedBlock.indexOf('initPathwayCanvas()') >= 0,
        'pathway canvas should lazy-initialize only when its collapsed drawer opens');
    assert.ok(js.indexOf('document.createElementNS(PATHWAY_NS') >= 0 &&
            js.indexOf('function sanitizePathwayText') >= 0 &&
            html.indexOf('data-wb-action="set-pathway-tool"') >= 0 &&
            html.indexOf('data-tool="curly"') >= 0 &&
            html.indexOf('id="wb-pathway-arrow-kind"') >= 0,
        'pathway canvas should use SVG DOM creation, sanitized text, declarative controls, and mechanism arrows');
    assert.ok(html.indexOf('data-tool="residue"') >= 0 &&
            html.indexOf('id="wb-pathway-shortcut"') >= 0 &&
            html.indexOf('data-wb-action="update-pathway-selection"') >= 0 &&
            js.indexOf('PATHWAY_SHORTCUTS') >= 0 &&
            js.indexOf('function updatePathwaySelectionFromControls') >= 0 &&
            js.indexOf("return value === 'custom' ? '' : value") >= 0 &&
            js.indexOf("kind === 'cofactor'") >= 0,
        'pathway canvas should support editable amino-acid and cofactor shortcuts without creating custom-labeled nodes');
    assert.ok(html.indexOf('id="wb-pathway-bg-input"') >= 0 &&
            html.indexOf('data-wb-action="import-pathway-background"') >= 0 &&
            js.indexOf('function importPathwayBackgroundFile') >= 0 &&
            js.indexOf('window.BIME_REACTION_IMAGE_RECOGNIZER') >= 0 &&
            js.indexOf('function applyPathwayDiagramResult') >= 0,
        'pathway canvas should support safe image/PDF reference import and private local recognition hooks');
    assert.ok(html.indexOf('data-wb-panel="templates" data-wb-collapsed="true"') >= 0,
        'template drawer should be collapsed by default so template SVGs are not generated on first paint');
    assert.ok(js.indexOf('function initWorkbenchDisclosures') >= 0 &&
            js.indexOf('function toggleWorkbenchSection') >= 0,
        'collapsed workbench drawers should have centralized accessible toggle wiring');
    assert.ok(js.indexOf('waitForBuiltInLibraryAndSearch') >= 0,
        'built-in search should wait for the async library before running');
    assert.ok(js.indexOf('parseSdfBlockToLibraryEntry') >= 0,
        'user-library SDF parsing should use the MOL parser path');
    assert.ok(js.indexOf('function moleculePerfSignature') >= 0 &&
            js.indexOf('_insightsCacheKey') >= 0,
        'workbench should cache repeated molecule insight calculations in one frame');
    assert.ok(js.indexOf('function moleculeChemSignature') >= 0,
        'SMILES output should use a connectivity cache key so coordinate-only drags avoid writer hotpaths');
    assert.ok(html.indexOf('id="wb-structure-name"') >= 0 &&
            html.indexOf('id="wb-structure-comment"') >= 0 &&
            html.indexOf('data-wb-action="apply-structure-properties"') >= 0 &&
            js.indexOf('function applyStructureProperties') >= 0 &&
            js.indexOf('function sanitizeStructureProperty') >= 0,
        'Properties tab should expose safe editable molecule/reaction identity metadata');
    assert.ok(js.indexOf("['Name', editor && editor.molecule && editor.molecule.name") >= 0 &&
            js.indexOf("var parts = [mol.atoms.length, mol.bonds ? mol.bonds.length : 0, mol.name || '', mol.comment || ''];") >= 0,
        'Properties and metadata-sensitive exports should reflect edited molecule names');
    assert.ok(html.indexOf('id="wb-reaction-property-editor"') >= 0 &&
            html.indexOf('data-wb-action="apply-property-reaction-labels"') >= 0 &&
            js.indexOf('function applyPropertyReactionLabels') >= 0 &&
            js.indexOf('function populateReactionPropertyControls') >= 0,
        'Properties tab should expose reaction label metadata when a reaction arrow is present');
    assert.ok(fs.readFileSync(path.join(ROOT, 'editor', 'Renderer.js'), 'utf8').indexOf('_screenRectCacheFrame') >= 0,
        'renderer coordinate conversion should cache SVG bounds per frame for pointer-move hotpaths');
    assert.ok(js.indexOf('function queueSearchThumb') >= 0 &&
            js.indexOf('function queueMolThumb') >= 0 &&
            js.indexOf('_searchThumbQueue') >= 0,
        'search result previews should render lazily off the critical search path');
    assert.ok(js.indexOf('function queueTemplateThumb') >= 0 &&
            js.indexOf('_templateThumbQueue') >= 0,
        'template thumbnails should render lazily after the collapsed template drawer opens');
    assert.ok(js.indexOf('function templateSearchText') >= 0 &&
            js.indexOf('Templates.aliases') >= 0,
        'template search should match compatibility aliases hidden from the public drawer list');
    assert.ok(js.indexOf('function scheduleLibraryWarmup') >= 0 &&
            js.indexOf("load('c1ccccc1');\n        warmLibraryCacheWhenReady") < 0,
        'library-cache warming should not compete with first editor paint');
    assert.ok(!/src="common-molecules\.js[^"]*"/.test(html) &&
            js.indexOf('function ensureBuiltInLibraryScript') >= 0 &&
            js.indexOf('common-molecules.js') >= 0 &&
            js.indexOf('function scheduleLibraryWarmup') >= 0,
        'common molecule library should be lazy-loaded only when search opens or runs');
    assert.ok(js.indexOf('opts.yieldEvery') >= 0,
        'workbench search should request chunked library scans for responsiveness');
});

module.exports = runner.summary;
if (require.main === module) {
    var s = runner.summary();
    console.log('\n' + s.passed + ' passed, ' + s.failed + ' failed');
    process.exit(s.failed > 0 ? 1 : 0);
}
