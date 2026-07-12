/* BIME Workbench controller. Deliberately external to keep CSP tight and button wiring testable. */
var editor;
// Focus-lens (merge Stage 2): structure-node inline editor. Guarded so an older
// bundle without FocusLens cleanly disables the lens (legacy bridge remains).
var _lens = (typeof FocusLens !== 'undefined' && FocusLens.create) ? FocusLens.create() : null;
// v2.3.4 (lens hardening): the element focused BEFORE the lens opened, restored
// on collapse (WCAG 2.4.3) — mirrors `_modalLastFocus` for the footer modal. And
// the active Tab-trap keydown handler, bound on the overlay while open and torn
// down on collapse so focus can't escape behind the dialog (WCAG 4.1.2). Both
// are safe as module globals because the focus lens is a singleton — at most one
// is open at a time (a hop closes the current lens before opening the next).
var _lensLastFocus = null;
var _lensTrapHandler = null;
// Hybrid Stage A (v2.4.0): which standalone VIEW the workbench is showing —
// 'pathway' (the canvas, the default boot view) or 'editor' (the full molecule/
// reaction editor as the in-flow hero). This is orthogonal to the focus lens:
// the lens is an in-place OVERLAY that can open over a node while the Pathway
// view is showing, so `_wbView` only tracks the standalone surface, never the
// lens. Default 'pathway' protects ?pathway= deep-links, the empty-canvas CTA,
// and first-run. setWorkbenchMode flips it and stamps the matching
// `.wb-view-*` class on `.wb-main`; the CSS owns the actual show/hide.
//
// v2.4.18: the pathway / metabolic-map editor ships as a SEPARATE release; the
// public core build is the molecule + reaction editor only. PATHWAY_ENABLED is
// the single source of truth that gates every pathway entry point — the default
// boot view, the ?pathway= deep-link, the affordances that live outside the
// pathway section, and the command palette. The pathway engine code, markup, and
// CSS are all preserved, so flipping this one constant to true restores the full
// pathway workspace with no other change.
var PATHWAY_ENABLED = false;
var _wbView = PATHWAY_ENABLED ? 'pathway' : 'editor';
var _outputFrame = 0;
var _searchSeq = 0;
var _activeOutputTab = 'smiles';
var _activeOutputFormat = 'smiles';
var _templateSvgCache = {};
var _templateThumbQueue = [];
var _templateThumbFrame = 0;
var _templateDrawerTimer = 0;
var _templateDrawerSeq = 0;
var _templateRecent = [];
var _commandItems = [];
var _workbenchFrame = 0;
var _lookupSeq = 0;
var _insightsCacheKey = '';
var _insightsCacheValue = null;
var _imageOutputCacheKey = '';
var _imageOutputCacheSvg = null;
var _outputTextCacheKey = '';
var _outputTextCache = {};
var _searchThumbQueue = [];
var _searchThumbFrame = 0;
var _libraryWarmupQueued = false;
var _libraryScriptPromise = null;
var _lastWarningCount = 0;
var _pathwayRenderFrame = 0;
var COMPOUND_LOOKUP_PROPS = 'MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,IUPACName,Title,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity';
var PATHWAY_NS = 'http://www.w3.org/2000/svg';
var PATHWAY_XLINK_NS = 'http://www.w3.org/1999/xlink';
var PATHWAY_STRUCTURE_TEXT_LIMIT = 512 * 1024;
var PATHWAY_SHORTCUTS = {
    Gly: { name: 'Glycine', kind: 'residue' },
    Ala: { name: 'Alanine', kind: 'residue' },
    Val: { name: 'Valine', kind: 'residue' },
    Leu: { name: 'Leucine', kind: 'residue' },
    Ile: { name: 'Isoleucine', kind: 'residue' },
    Ser: { name: 'Serine', kind: 'residue' },
    Thr: { name: 'Threonine', kind: 'residue' },
    Cys: { name: 'Cysteine', kind: 'residue' },
    Met: { name: 'Methionine', kind: 'residue' },
    Asp: { name: 'Aspartate', kind: 'residue' },
    Glu: { name: 'Glutamate', kind: 'residue' },
    Asn: { name: 'Asparagine', kind: 'residue' },
    Gln: { name: 'Glutamine', kind: 'residue' },
    Lys: { name: 'Lysine', kind: 'residue' },
    Arg: { name: 'Arginine', kind: 'residue' },
    His: { name: 'Histidine', kind: 'residue' },
    Phe: { name: 'Phenylalanine', kind: 'residue' },
    Tyr: { name: 'Tyrosine', kind: 'residue' },
    Trp: { name: 'Tryptophan', kind: 'residue' },
    Pro: { name: 'Proline', kind: 'residue' },
    ATP: { name: 'ATP', kind: 'cofactor' },
    ADP: { name: 'ADP', kind: 'cofactor' },
    'NAD+': { name: 'NAD+', kind: 'cofactor' },
    NADH: { name: 'NADH', kind: 'cofactor' },
    'NADP+': { name: 'NADP+', kind: 'cofactor' },
    NADPH: { name: 'NADPH', kind: 'cofactor' },
    CoA: { name: 'CoA', kind: 'cofactor' },
    Pi: { name: 'Inorganic phosphate', kind: 'cofactor' },
    PPi: { name: 'Pyrophosphate', kind: 'cofactor' },
    'H+': { name: 'Proton', kind: 'cofactor' },
    H2O: { name: 'Water', kind: 'cofactor' },
    CO2: { name: 'Carbon dioxide', kind: 'cofactor' }
};
var _pathway = {
    tool: 'select',
    scale: 1,
    offsetX: 0,
    offsetY: 0,
    nextId: 1,
    backgrounds: [],
    nodes: [],
    edges: [],
    notes: [],
    steps: [],
    mechanismArrows: [],
    compartments: [],
    selectedType: '',
    selectedId: null,
    pendingArrow: null,
    pendingMechanismArrow: null,
    drag: null,
    initialized: false,
    // v2.0.38: print page-format awareness. 'screen' = legacy 1200×620
    // viewBox; 'a4' / 'letter' switch to the corresponding paper size.
    pageFormat: 'screen',
    pageOrientation: 'landscape',
    // Reusable topology-aware layout strategy for Clean Up. `auto` lets the
    // layout engine choose cycle, branch, merge, or layered flow from graph shape.
    layoutShape: 'auto',
    // v2.0.54: moiety-trace selection state. While `moietySet` has ≥ 2
    // entries the atom-trace strip runs `AtomTrace.deriveMoietyTrace` and
    // colours each cell by status (intact / fragmented / partial-loss /
    // absent). Reset when the user clicks a non-start-metabolite atom or
    // explicitly clears the moiety.
    moietyStartMetIdx: null,
    moietySet: []
};
var PATHWAY_BACKGROUND_IMAGE_LIMIT = 12 * 1024 * 1024;
var PATHWAY_BACKGROUND_PDF_LIMIT = 24 * 1024 * 1024;

function afterPaint(fn) {
    // Run `fn` after first paint, but NEVER depend solely on
    // requestAnimationFrame: rAF is throttled/paused in background or hidden
    // tabs, which would otherwise leave the whole workbench un-booted (a blank
    // editor) until the tab gains focus. Race a double-rAF against a short
    // timeout — whichever fires first wins — guarded so `fn` runs exactly once.
    var ran = false;
    function run() { if (ran) { return; } ran = true; fn(); }
    if (typeof requestAnimationFrame === 'function') {
        requestAnimationFrame(function() { requestAnimationFrame(run); });
    }
    setTimeout(run, 250);
}

function runIdle(fn) {
    if (typeof requestIdleCallback === 'function') {
        requestIdleCallback(fn, { timeout: 2500 });
    } else {
        setTimeout(fn, 250);
    }
}

document.addEventListener('DOMContentLoaded', function() {
    editor = new MolEditor('bime-editor', '100%', '100%', { options: 'hydrogens' });

    editor.setCallBack('AfterStructureModified', function() {
        queueOutputUpdate();
        scheduleWorkbenchUpdate();
    });

    initWorkbenchUX();

    afterPaint(function() {
        // Merge Stage 3 (v2.3.1): the pathway canvas is ONE continuous surface —
        // the permanent workspace. Mount it on boot (no Molecule/Pathway mode
        // swap exists anymore); structure editing happens through the focus
        // lens. The molecule editor still loads any deep-linked structure so it
        // remains the capture source-of-truth, but it only ever surfaces as the
        // lens overlay over a node.
        var _bootMain = document.querySelector('.wb-main');
        if (PATHWAY_ENABLED) {
            // Hybrid Stage A (v2.4.0): the DEFAULT boot view is Pathway — it
            // protects ?pathway= deep-links, the v2.3.5 empty-canvas CTA, and
            // first-run. Mount the canvas and stamp `.wb-view-pathway` on
            // `.wb-main` from first paint (so the editor stays out of flow via
            // the base rule) and select the Pathway tab. We do NOT route boot
            // through setWorkbenchMode('pathway') so the deep-link load below
            // runs exactly as before; this only sets the standalone-view class.
            openWorkbenchSection('wb-pathway-section');
            if (_bootMain) {
                _bootMain.classList.add('wb-view-pathway');
                _bootMain.classList.remove('wb-view-editor');
            }
            syncWorkbenchModeButtons(_wbView);
        } else {
            // v2.4.18: boot straight into
            // the molecule/reaction editor as the standalone hero and hide every
            // pathway affordance that lives OUTSIDE the (now-unmounted) pathway
            // section (the section itself is hidden by the Editor-view CSS).
            hidePathwayAffordances();
            setWorkbenchMode('editor');
        }
        var initial = getInitialWorkbenchLoad();
        if (PATHWAY_ENABLED && initial && initial.pathway && typeof loadExamplePathway === 'function') {
            loadExamplePathway(initial.pathway);
        } else if (initial && initial.input) {
            loadNamed(initial.input, initial.name);
        } else if (PATHWAY_ENABLED) {
            // v2.3.5 (P0-5): no deep-link content, so neither load path above ran
            // and nothing else paints the canvas on boot. Render the empty
            // pathway explicitly so the onboarding CTA (drawEmptyCanvasCta) shows
            // on FIRST paint, not only after the first user action.
            renderPathwayCanvas();
        }
        // v2.3.5 (P0-4): no DEFAULT boot molecule. The editor is display:none
        // until a focus lens opens, so loading a structure here put SMILES +
        // properties into Editor Output for a molecule the user could see
        // NOWHERE — the canvas stayed empty while Output showed benzene's
        // SMILES, so the two panels disagreed and read as a bug. On a
        // no-deep-link boot we leave the editor empty: Editor Output shows its
        // empty-state placeholder and the canvas shows the onboarding CTA
        // (drawEmptyCanvasCta). Deep links (?smiles= / ?input= / ?pathway=)
        // still load above.
    });
});

/* ---- Core functions ---- */
// v2.4.18: hide the pathway affordances that live OUTSIDE #wb-pathway-section
// (the section itself is hidden by the Editor-view CSS). Inline display:none
// beats the affordances' class-level display rules, so the hide is bulletproof;
// the markup stays in the DOM.
function hidePathwayAffordances() {
    if (typeof document === 'undefined') return;
    var selectors = [
        '.wb-modebar',                                    // Editor|Pathway view switch (only Editor remains)
        '[data-wb-action="reaction-to-pathway"]',         // "To Pathway" (Reaction Scheme panel)
        '[data-wb-action="editor-add-to-pathway"]',       // "Add to pathway" (Editor Output row)
        '.wb-mobile-dock [data-wb-action="pathway-fit"]', // mobile dock Fit (pathway canvas)
        '#wb-atom-trace-section'                           // pathway-only atom-trace strip
    ];
    for (var i = 0; i < selectors.length; i++) {
        var nodes = document.querySelectorAll(selectors[i]);
        for (var j = 0; j < nodes.length; j++) {
            nodes[j].hidden = true;
            if (nodes[j].style) { nodes[j].style.setProperty('display', 'none', 'important'); }
        }
    }
}

function getInitialWorkbenchLoad() {
    if (typeof URLSearchParams === 'undefined') return null;
    try {
        var qs = new URLSearchParams(window.location.search || '');
        // v2.0.21: ?pathway=<id> opens a built-in example pathway directly.
        var pathway = (qs.get('pathway') || '').trim().toLowerCase();
        if (pathway) {
            if (typeof MetaboliteLibrary !== 'undefined' && MetaboliteLibrary.pathways &&
                MetaboliteLibrary.pathways[pathway]) {
                return { pathway: pathway };
            }
            return null;
        }
        var input = (qs.get('smiles') || qs.get('input') || '').trim();
        if (!input || input.length > 10000) return null;
        var name = (qs.get('name') || '').trim();
        if (name.length > 80) name = name.slice(0, 80);
        return { input: input, name: name };
    } catch (e) {
        return null;
    }
}

function queueOutputUpdate() {
    if (_outputFrame) return;
    var raf = typeof requestAnimationFrame === 'function'
        ? requestAnimationFrame
        : function(cb) { return setTimeout(cb, 16); };
    _outputFrame = raf(function() {
        _outputFrame = 0;
        updateOutputNow();
    });
}

function cancelOutputUpdate() {
    if (!_outputFrame) return;
    if (typeof cancelAnimationFrame === 'function') {
        cancelAnimationFrame(_outputFrame);
    } else {
        clearTimeout(_outputFrame);
    }
    _outputFrame = 0;
}

function updateOutputNow() {
    if (!editor) return;
    var out = document.getElementById('smiles-out');
    if (out) out.value = outputTextFor(_activeOutputFormat || _activeOutputTab || 'smiles');
    updateOutputPanels();
    updateStats();
    updateReactionEquationReadout();
    // v2.3.8 (P1-4): mirror the SMILES + warning count into the open lens's live
    // readout strip (a no-op when no lens is open). Reuses the cached signals
    // this function already computed — no extra recompute.
    updateLensReadout();
}

function load(s) {
    if (!editor) return;
    editor.readGenericMolecularInput(s);
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
}

function loadReactionScheme(s, labelParts) {
    load(s);
    if (!editor || !editor.molecule || !editor.molecule.reactionArrow) return;
    setReactionLabelParts(labelParts || {});
    editor.render();
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
}

function loadNamed(s, name) {
    if (!editor) return;
    editor.readGenericMolecularInput(s);
    if (name && editor.molecule) {
        editor.molecule.name = name;
        editor.renderer.showMolName = true;
        editor.render();
    }
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
}

function doExport(f) {
    if (!editor) return;
    cancelOutputUpdate();
    var o = document.getElementById('smiles-out');
    if (!o) return;
    _activeOutputFormat = f || 'smiles';
    var text = outputTextFor(_activeOutputFormat);
    if (_activeOutputFormat === 'rxn' && !text && editor.showInfo) {
        editor.showInfo('RXN export requires a reaction arrow.');
    }
    o.value = text;
    setOutputTab(f === 'smiles' ? 'smiles' : 'mol', true);
    updateOutputPanels();
}

function copyOut() {
    var out = document.getElementById('smiles-out');
    var t = out ? out.value : '';
    if (!t || !navigator.clipboard || !navigator.clipboard.writeText) return;
    navigator.clipboard.writeText(t).then(function() {
        if (editor && editor.showInfo) editor.showInfo('Copied!');
    });
}
function updateStats() {
    if (!editor) return;
    var el = document.getElementById('stats');
    if (!el) return;
    var insights = getMoleculeInsights();
    el.textContent = insights.atoms + ' atoms \u00b7 ' + insights.bonds + ' bonds \u00b7 ' +
        insights.rings + ' rings \u00b7 ' + insights.formula;
}

function outputTextFor(format) {
    if (!editor) return '';
    var layoutSensitive = format === 'mol' || format === 'mol3' || format === 'sdf' || format === 'rxn';
    var cacheKey = (layoutSensitive ? moleculePerfSignature : moleculeChemSignature)(editor.molecule);
    if (_outputTextCacheKey !== cacheKey) {
        _outputTextCacheKey = cacheKey;
        _outputTextCache = {};
    }
    if (_outputTextCache[format] !== undefined) return _outputTextCache[format];
    var value = '';
    try {
        if (format === 'mol') value = editor.molFile(false);
        else if (format === 'mol3') value = editor.molFile(true);
        if (format === 'sdf' && typeof MolfileWriter !== 'undefined' && MolfileWriter.writeSDF) {
            value = MolfileWriter.writeSDF(editor.molecule);
        }
        else if (format === 'rxn') {
            if (editor.rxnFile) value = editor.rxnFile(false);
            else if (typeof MolfileWriter !== 'undefined' && MolfileWriter.writeRXN) {
                value = MolfileWriter.writeRXN(editor.molecule);
            }
        } else if (format !== 'mol' && format !== 'mol3') {
            value = editor.smiles();
        }
    } catch (e) {
        value = '';
    }
    _outputTextCache[format] = value || '';
    return _outputTextCache[format];
}

function setOutputTab(tab, keepFormat) {
    _activeOutputTab = tab || 'smiles';
    if (!keepFormat) {
        _activeOutputFormat = _activeOutputTab === 'mol' ? 'mol' : 'smiles';
    }
    var tabs = ['smiles', 'mol', 'image', 'props', 'warnings'];
    for (var i = 0; i < tabs.length; i++) {
        var id = tabs[i];
        var tabEl = document.getElementById('wb-tab-' + id);
        if (tabEl) tabEl.setAttribute('aria-selected', id === _activeOutputTab ? 'true' : 'false');
    }
    var textPanel = document.getElementById('wb-panel-text');
    var imagePanel = document.getElementById('wb-panel-image');
    var propsPanel = document.getElementById('wb-panel-props');
    var warningsPanel = document.getElementById('wb-panel-warnings');
    if (textPanel) {
        textPanel.classList.toggle('is-active', _activeOutputTab === 'smiles' || _activeOutputTab === 'mol');
        textPanel.setAttribute('aria-labelledby', _activeOutputTab === 'mol' ? 'wb-tab-mol' : 'wb-tab-smiles');
    }
    if (imagePanel) imagePanel.classList.toggle('is-active', _activeOutputTab === 'image');
    if (propsPanel) propsPanel.classList.toggle('is-active', _activeOutputTab === 'props');
    if (warningsPanel) warningsPanel.classList.toggle('is-active', _activeOutputTab === 'warnings');
    var out = document.getElementById('smiles-out');
    if (out && (_activeOutputTab === 'smiles' || _activeOutputTab === 'mol')) {
        out.value = outputTextFor(_activeOutputFormat);
    }
    updateOutputPanels();
}

function updateOutputPanels() {
    if (_activeOutputTab === 'image') renderImageOutput();
    if (_activeOutputTab === 'props') renderPropertiesOutput();
    if (_activeOutputTab === 'warnings') renderWarningsOutput();
}

function renderImageOutput() {
    var el = document.getElementById('wb-output-image');
    if (!el) return;
    if (_activeOutputTab !== 'image') return;
    if (!editor || !editor.molecule || editor.molecule.isEmpty()) {
        _imageOutputCacheKey = '';
        _imageOutputCacheSvg = null;
        setSafeSvg(el, '', 'Draw or load a structure to preview SVG output.');
        return;
    }
    if (typeof ImageExport === 'undefined' || !ImageExport.toSVG) {
        setSafeSvg(el, '', 'Image export is unavailable in this bundle.');
        return;
    }
    try {
        var key = moleculePerfSignature(editor.molecule) + '|760x300|image';
        var svg = _imageOutputCacheKey === key ? _imageOutputCacheSvg : null;
        if (!svg) {
            svg = ImageExport.toSVG(editor.molecule, {
                width: 760,
                height: 300,
                padding: 16,
                background: 'transparent',
                showAtomLabels: true
            });
            _imageOutputCacheKey = key;
            _imageOutputCacheSvg = svg || null;
        }
        setSafeSvg(el, svg, 'Preview unavailable.');
    } catch (e) {
        setSafeSvg(el, '', 'Preview unavailable.');
    }
}

function renderPropertiesOutput() {
    var el = document.getElementById('wb-properties-output');
    if (!el) return;
    if (_activeOutputTab !== 'props') return;
    populateStructurePropertyControls();
    populateReactionPropertyControls();
    el.textContent = '';
    var insights = getMoleculeInsights();
    var rows = [
        ['Type', structureKindLabel()],
        ['Name', editor && editor.molecule && editor.molecule.name ? editor.molecule.name : '—'],
        ['Comment', editor && editor.molecule && editor.molecule.comment ? editor.molecule.comment : '—'],
        ['Formula', insights.formula],
        ['Mol. Weight', insights.molWeightLabel],
        ['Atoms', insights.atoms],
        ['Bonds', insights.bonds],
        ['Rings', insights.rings],
        ['Net Charge', insights.charge],
        ['HBD / HBA', insights.hbd + ' / ' + insights.hba],
        ['Rotatable Bonds', insights.rotatableBonds],
        ['Quality', insights.qualityLabel]
    ];
    var lookup = currentCompoundLookupData();
    if (lookup) {
        rows.push(['Compound CID', lookup.cid || '—']);
        rows.push(['IUPAC', lookup.iupacName || '—']);
        rows.push(['XLogP', formatProperty(lookup.xlogp)]);
        rows.push(['TPSA', formatProperty(lookup.tpsa)]);
        rows.push(['Exact Mass', formatProperty(lookup.exactMass)]);
    }
    var smilesNow = editor && editor.smiles ? editor.smiles() : '';
    var rxnEq = reactionFormulaFromSmiles(smilesNow);
    if (rxnEq) {
        rows.push(['Equation', rxnEq.equation]);
    }
    var rxn = reactionInsightsFromSmiles(smilesNow);
    if (rxn) {
        rows.push(['Reactants', rxn.reactants.join(' + ') || '—']);
        rows.push(['Products', rxn.products.join(' + ') || '—']);
        rows.push(['Arrow Label', reactionLabelForArrow(editor.molecule && editor.molecule.reactionArrow) || '—']);
        rows.push(['Atom Balance', rxn.balanced ? 'Balanced' : 'Check']);
    }
    rows.forEach(function(row) {
        var box = document.createElement('div');
        box.className = 'wb-prop';
        var k = document.createElement('span');
        k.className = 'wb-prop-k';
        k.textContent = row[0];
        var v = document.createElement('span');
        v.className = 'wb-prop-v';
        v.textContent = String(row[1]);
        box.appendChild(k);
        box.appendChild(v);
        el.appendChild(box);
    });
}

function structureKindLabel() {
    return editor && editor.molecule && editor.molecule.reactionArrow ? 'Reaction' : 'Molecule';
}

function sanitizeStructureProperty(value, maxLen) {
    return String(value || '')
        .replace(/[\r\n\t]+/g, ' ')
        .replace(/\s+/g, ' ')
        .trim()
        .slice(0, maxLen || 120);
}

function populateStructurePropertyControls() {
    var mol = editor && editor.molecule;
    var kind = structureKindLabel();
    setText('wb-property-kind', kind);
    var active = document.activeElement;
    var editing = active && (
        active.id === 'wb-structure-name' ||
        active.id === 'wb-structure-comment' ||
        active.id === 'wb-structure-show-name'
    );
    if (editing) return;
    setInputValue('wb-structure-name', mol ? (mol.name || '') : '');
    setInputValue('wb-structure-comment', mol ? (mol.comment || '') : '');
    var show = document.getElementById('wb-structure-show-name');
    if (show) {
        show.checked = !!(editor && editor.renderer && editor.renderer.showMolName && mol && mol.name);
        show.disabled = !mol || !mol.atoms || mol.atoms.length === 0;
    }
}

function invalidateWorkbenchOutputCaches() {
    _outputTextCacheKey = '';
    _outputTextCache = {};
    _imageOutputCacheKey = '';
    _imageOutputCacheSvg = null;
    _insightsCacheKey = '';
    _insightsCacheValue = null;
}

function applyStructureProperties() {
    if (!editor || !editor.molecule) return;
    var nameInput = document.getElementById('wb-structure-name');
    var commentInput = document.getElementById('wb-structure-comment');
    var showInput = document.getElementById('wb-structure-show-name');
    var name = sanitizeStructureProperty(nameInput && nameInput.value, 80);
    var comment = sanitizeStructureProperty(commentInput && commentInput.value, 180);
    var showName = !!(showInput && showInput.checked && name);
    if (showInput && !name) showInput.checked = false;
    var mol = editor.molecule;
    var changed = (mol.name || '') !== name ||
        (mol.comment || '') !== comment ||
        !!(editor.renderer && editor.renderer.showMolName) !== showName;
    if (!changed) {
        setText('wb-structure-property-status', 'Properties already up to date.');
        return;
    }
    if (editor.saveHistory) editor.saveHistory();
    mol.name = name;
    mol.comment = comment;
    setInputValue('wb-structure-name', name);
    setInputValue('wb-structure-comment', comment);
    if (editor.renderer) editor.renderer.showMolName = showName;
    if (typeof editor._updateNameStatus === 'function') editor._updateNameStatus();
    if (editor.render) editor.render();
    invalidateWorkbenchOutputCaches();
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
    setText('wb-structure-property-status', (structureKindLabel() === 'Reaction' ? 'Reaction' : 'Molecule') + ' properties applied.');
}

function clearStructureName() {
    setInputValue('wb-structure-name', '');
    var show = document.getElementById('wb-structure-show-name');
    if (show) show.checked = false;
    applyStructureProperties();
}

function renderWarningsOutput() {
    var target = document.getElementById('wb-warnings-output');
    if (!target) return;
    if (_activeOutputTab !== 'warnings') return;
    renderWarningList(target, getMoleculeInsights().warnings);
}

function warmLibraryCacheWhenReady() {
    if (!editor || window.BIME_LIBRARY_WARMED) return;
    if (!builtInLibraryReady()) {
        ensureBuiltInLibraryScript().then(warmLibraryCacheWhenReady).catch(function() {});
        return;
    }
    window.BIME_LIBRARY_WARMED = true;
    runIdle(function() {
        if (!editor || typeof editor.searchLibrary !== 'function') return;
        editor.searchLibrary('C', 'exact', { topN: 1 }).catch(function() {});
    });
}

function scheduleLibraryWarmup() {
    if (_libraryWarmupQueued || window.BIME_LIBRARY_WARMED) return;
    _libraryWarmupQueued = true;
    runIdle(function() {
        _libraryWarmupQueued = false;
        warmLibraryCacheWhenReady();
    });
}

function currentAssetQuery() {
    var scripts = document.getElementsByTagName ? document.getElementsByTagName('script') : [];
    for (var i = 0; i < scripts.length; i++) {
        var src = scripts[i].getAttribute && scripts[i].getAttribute('src');
        if (src && src.indexOf('js/workbench.js') >= 0) {
            var q = src.indexOf('?');
            return q >= 0 ? src.slice(q) : '';
        }
    }
    return '';
}

function ensureBuiltInLibraryScript() {
    if (builtInLibraryReady()) return Promise.resolve(COMMON_MOLECULES);
    if (_libraryScriptPromise) return _libraryScriptPromise;
    _libraryScriptPromise = new Promise(function(resolve, reject) {
        var existing = document.querySelector && document.querySelector('script[data-bime-common-molecules]');
        if (!existing && document.getElementsByTagName) {
            var scripts = document.getElementsByTagName('script');
            for (var i = 0; i < scripts.length; i++) {
                var src = scripts[i].getAttribute && scripts[i].getAttribute('src');
                if (src && src.indexOf('common-molecules.js') >= 0) {
                    existing = scripts[i];
                    break;
                }
            }
        }
        var script = existing || document.createElement('script');
        script.setAttribute('data-bime-common-molecules', 'true');
        script.async = true;
        script.onload = function() {
            if (builtInLibraryReady()) resolve(COMMON_MOLECULES);
            else reject(new Error('common molecule library loaded empty'));
        };
        script.onerror = function() {
            reject(new Error('could not load common molecule library'));
        };
        if (!existing) {
            script.src = 'common-molecules.js' + currentAssetQuery();
            (document.body || document.head || document.documentElement).appendChild(script);
        }
    });
    return _libraryScriptPromise;
}

/* ---- Workbench UX: inspector, templates, output tabs, command surface ---- */
function scheduleWorkbenchUpdate() {
    if (_workbenchFrame) return;
    var raf = typeof requestAnimationFrame === 'function'
        ? requestAnimationFrame
        : function(cb) { return setTimeout(cb, 16); };
    _workbenchFrame = raf(function() {
        _workbenchFrame = 0;
        updateInspector();
        if (!_outputFrame) updateOutputPanels();
    });
}

function initWorkbenchUX() {
    try {
        var stored = localStorage.getItem('bime-template-recent');
        _templateRecent = stored ? JSON.parse(stored) : [];
        if (!_templateRecent || !_templateRecent.slice) _templateRecent = [];
    } catch (e) {
        _templateRecent = [];
    }
    initWorkbenchDisclosures();
    bindWorkbenchActions();
    initTemplateDrawer();
    initCommandPalette();
    initCompoundLookup();
    initPropertyControls();
    initReactionLabelControls();
    initEditorSurfaceEvents();
    setOutputTab('smiles', true);
    scheduleWorkbenchUpdate();
}

function initWorkbenchDisclosures() {
    var panels = document.querySelectorAll('.wb-panel[data-wb-collapsible]');
    for (var i = 0; i < panels.length; i++) {
        enhanceWorkbenchPanel(panels[i], i);
    }
    syncWorkbenchDrawerButtons();
}

function enhanceWorkbenchPanel(panel, index) {
    var head = panel.querySelector('.wb-panel-head');
    var body = panel.querySelector('.wb-panel-body');
    if (!head || !body) return;
    var panelName = panel.getAttribute('data-wb-panel') || ('panel-' + index);
    var bodyId = body.id || ('wb-panel-body-' + panelName);
    body.id = bodyId;
    head.setAttribute('role', 'button');
    head.setAttribute('tabindex', '0');
    head.setAttribute('aria-controls', bodyId);
    var collapsed = panel.getAttribute('data-wb-collapsed') === 'true';
    setWorkbenchPanelCollapsed(panel, collapsed);
    head.addEventListener('click', function(ev) {
        if (ev.target && ev.target.closest && ev.target.closest('button,input,select,textarea,a')) return;
        setWorkbenchPanelCollapsed(panel, !panel.classList.contains('is-collapsed'));
    });
    head.addEventListener('keydown', function(ev) {
        if (ev.key === 'Enter' || ev.key === ' ') {
            ev.preventDefault();
            setWorkbenchPanelCollapsed(panel, !panel.classList.contains('is-collapsed'));
        }
    });
}

function setWorkbenchPanelCollapsed(panel, collapsed) {
    var body = panel && panel.querySelector ? panel.querySelector('.wb-panel-body') : null;
    var head = panel && panel.querySelector ? panel.querySelector('.wb-panel-head') : null;
    if (!panel || !body || !head) return;
    panel.classList.toggle('is-collapsed', !!collapsed);
    body.hidden = !!collapsed;
    head.setAttribute('aria-expanded', collapsed ? 'false' : 'true');
    if (!collapsed) {
        var name = panel.getAttribute('data-wb-panel') || '';
        if (name === 'templates') renderTemplateDrawer();
        if (name === 'command') renderCommandResults();
        if (name === 'warnings') {
            var warnPanel = document.getElementById('wb-warning-panel');
            if (warnPanel) renderWarningList(warnPanel, getMoleculeInsights().warnings);
        }
    }
}

function openWorkbenchPanel(name) {
    var panel = document.querySelector('.wb-panel[data-wb-panel="' + name + '"]');
    if (panel) setWorkbenchPanelCollapsed(panel, false);
}

function isWorkbenchPanelExpanded(name) {
    var panel = document.querySelector('.wb-panel[data-wb-panel="' + name + '"]');
    return !panel || !panel.classList.contains('is-collapsed');
}

function toggleWorkbenchSection(button) {
    if (!button) return;
    var targetId = button.getAttribute('data-target');
    var group = button.getAttribute('data-target-group');
    var open = button.getAttribute('aria-expanded') !== 'true';
    if (targetId) {
        setWorkbenchSectionOpen(targetId, open);
    } else if (group) {
        setWorkbenchSectionGroupOpen(group, open);
    }
}

function setWorkbenchSectionOpen(id, open) {
    var section = document.getElementById(id);
    if (!section) return;
    section.hidden = !open;
    syncWorkbenchDrawerButtons();
    syncWorkbenchLayoutMode();
    if (open) onWorkbenchSectionOpened(id);
}

function setWorkbenchSectionGroupOpen(group, open) {
    var sections = document.querySelectorAll('[data-wb-section-group="' + group + '"]');
    for (var i = 0; i < sections.length; i++) sections[i].hidden = !open;
    syncWorkbenchDrawerButtons();
    syncWorkbenchLayoutMode();
}

function openWorkbenchSection(id) {
    setWorkbenchSectionOpen(id, true);
}

// Hybrid Stage A (v2.4.0): setWorkbenchMode now drives the explicit
// "Editor | Pathway" VIEW SWITCH, layered on the v2.3.x focus-lens continuum.
// It KEEPS its older compatibility contract (the source-shape tests pin both
// the literal `openWorkbenchSection('wb-pathway-section')` here and the
// `setWorkbenchMode('molecule')` callers elsewhere) and ADDS real switching:
//   - the pathway canvas stays MOUNTED in every view (we always
//     `openWorkbenchSection('wb-pathway-section')` — never hide it via
//     setWorkbenchSectionOpen(...,false), which v2.3.1 A3 forbids); the
//     standalone surface is chosen by the `.wb-view-*` class on `.wb-main`,
//     and CSS owns the show/hide;
//   - 'editor' shows the full molecule/reaction editor as the in-flow hero
//     and re-fits it (it was display:none until now, so it needs a measure);
//   - any other value (incl. the legacy 'molecule'/'pathway' callers) lands on
//     the Pathway view — the default, deep-link-safe surface.
// The focus lens is ORTHOGONAL: it is an in-place overlay keyed on
// `.wb-lens-open` (with `!important`), so it still opens over a node while the
// Pathway view is showing. We re-fit the editor whenever it becomes visible —
// as the hero (Editor view) or as the lens overlay — so a just-shown,
// previously-unsized #bime-editor can't render under-fitted.
function setWorkbenchMode(mode) {
    // v2.4.18: with the pathway editor hidden the only standalone view is the
    // molecule/reaction editor, so force it — every legacy setWorkbenchMode(...)
    // caller (reaction->pathway, add-to-pathway, imports) then lands on Editor.
    if (!PATHWAY_ENABLED) { mode = 'editor'; }
    // The canvas is the permanent workspace — keep it mounted in every view.
    openWorkbenchSection('wb-pathway-section');
    _wbView = (mode === 'editor') ? 'editor' : 'pathway';
    var main = document.querySelector('.wb-main');
    if (main) {
        main.classList.toggle('wb-view-editor', _wbView === 'editor');
        main.classList.toggle('wb-view-pathway', _wbView !== 'editor');
    }
    syncWorkbenchModeButtons(_wbView);
    // Hybrid Stage D (v2.4.3): on mobile the dock (Draw/Fit/Inspect) is a
    // PATHWAY-canvas affordance — Draw opens the lens, Fit fits the canvas. In
    // the standalone Editor view the editor itself is the drawing surface and
    // "Add to pathway" lives in the Output row, so the dock's Draw/Fit are
    // redundant/confusing there. Hide the dock while the Editor view owns the
    // screen (its own flag, orthogonal to the lens-sheet `.wb-dock-hidden`), and
    // restore it in the Pathway view where Draw still opens the lens. Off-mobile
    // the dock is already display:none, so this is a harmless no-op.
    setMobileDockEditorHidden(_wbView === 'editor');
    if (_wbView === 'editor') {
        // Hybrid Stage B (v2.4.1): entering the Editor view clears any pathway
        // node selection so a subsequent "Add to pathway" is unambiguous — it
        // always stamps a FRESH node and can never edit a stale selection.
        if (typeof selectPathwayItem === 'function') { selectPathwayItem('', null); }
        // The editor is the hero now (it was out of flow): re-measure + redraw.
        refitMoleculeEditorSoon();
    } else if (_lens && _lens.isOpen()) {
        // Back on the canvas with a lens open over a node — keep it fitted.
        refitMoleculeEditorSoon();
    }
}

// v2.0.17: re-measure the molecule editor container and resize/redraw once the
// browser has laid it out. Used when the single-view toggle re-shows the editor.
function refitMoleculeEditorSoon() {
    if (!editor) return;
    var attempts = 0;
    var apply = function() {
        var host = document.getElementById('bime-editor');
        if (host && editor.setSize && host.clientWidth > 0) {
            editor.setSize(host.clientWidth, host.clientHeight);
            if (editor.render) editor.render();
            return;
        }
        // Container not laid out yet (the editor was just un-hidden): retry on
        // the next few frames so a slow layout pass can't leave the canvas
        // under-fitted, then render regardless as a final fallback.
        if (attempts++ < 3 && typeof requestAnimationFrame === 'function') {
            requestAnimationFrame(apply);
            return;
        }
        if (editor.render) editor.render();
    };
    if (typeof requestAnimationFrame === 'function') { requestAnimationFrame(apply); }
    else { setTimeout(apply, 16); }
}

// v2.0.17: dependency-free footer modal (Acknowledgments / License / Contact).
// Backdrop click, the close button, and Escape all dismiss it; focus moves to
// the dialog and is restored to the trigger on close.
var _modalLastFocus = null;
function openWorkbenchModal(name) {
    var root = document.getElementById('wb-modal-root');
    if (!root || !name) return;
    var panes = root.querySelectorAll('[data-modal-pane]');
    var found = false;
    for (var i = 0; i < panes.length; i++) {
        var match = panes[i].getAttribute('data-modal-pane') === name;
        panes[i].hidden = !match;
        if (match) found = true;
    }
    if (!found) return;
    var titleEl = document.getElementById('wb-modal-title');
    var titles = { ack: 'Acknowledgments', license: 'License', contact: 'Contact' };
    if (titleEl) titleEl.textContent = titles[name] || 'About BIME';
    _modalLastFocus = document.activeElement;
    root.hidden = false;
    if (document.body) document.body.classList.add('wb-modal-open');
    var closeBtn = root.querySelector('.wb-modal-close');
    if (closeBtn && closeBtn.focus) { try { closeBtn.focus(); } catch (e) {} }
}
function closeWorkbenchModal() {
    var root = document.getElementById('wb-modal-root');
    if (!root || root.hidden) return;
    root.hidden = true;
    if (document.body) document.body.classList.remove('wb-modal-open');
    if (_modalLastFocus && _modalLastFocus.focus) { try { _modalLastFocus.focus(); } catch (e) {} }
    _modalLastFocus = null;
}

function syncWorkbenchDrawerButtons() {
    var buttons = document.querySelectorAll('[data-wb-action="toggle-workbench-section"]');
    for (var i = 0; i < buttons.length; i++) {
        var btn = buttons[i];
        var targetId = btn.getAttribute('data-target');
        var group = btn.getAttribute('data-target-group');
        var expanded = false;
        if (targetId) {
            var target = document.getElementById(targetId);
            expanded = !!(target && !target.hidden);
        } else if (group) {
            var sections = document.querySelectorAll('[data-wb-section-group="' + group + '"]');
            for (var j = 0; j < sections.length; j++) {
                if (!sections[j].hidden) { expanded = true; break; }
            }
        }
        btn.setAttribute('aria-expanded', expanded ? 'true' : 'false');
    }
    syncWorkbenchLayoutMode();
}

// Merge Stage 3 (v2.3.1): COMPATIBILITY SHIM. The full-view mode swap is gone,
// so this no longer flips the layout between an editor hero and a pathway hero.
// The canvas is the permanent workspace; the editor surfaces only as the lens.
// We keep the function symbol and its callers intact, and still stamp the
// legacy `.is-pathway-open` marker on `.wb-main` (now meaning "the continuum
// workspace is mounted") so any CSS / test keyed on it keeps resolving.
function syncWorkbenchLayoutMode() {
    var main = document.querySelector('.wb-main');
    var pathway = document.getElementById('wb-pathway-section');
    var pathwayOpen = !!(pathway && !pathway.hidden);
    if (!main) return;
    main.classList.toggle('is-pathway-open', pathwayOpen);
    // Hybrid Stage A (v2.4.0): the active VIEW is owned by `_wbView` (the
    // Editor | Pathway switch), not by whether the canvas section is mounted
    // (it always is now). Reflect the current view on the tabs so opening other
    // drawers can't steal the Editor tab's selected state.
    syncWorkbenchModeButtons(_wbView);
}

// Hybrid Stage A (v2.4.0): reflect the active standalone view on the
// Editor | Pathway tabs. They are role="tab", so `aria-selected` is the primary
// state; we also keep `aria-pressed` in sync for any legacy CSS/tests that key
// on it. `mode` is the view name ('editor' / 'pathway').
function syncWorkbenchModeButtons(mode) {
    var buttons = document.querySelectorAll('.wb-mode-btn[data-mode]');
    for (var i = 0; i < buttons.length; i++) {
        var on = buttons[i].getAttribute('data-mode') === mode ? 'true' : 'false';
        buttons[i].setAttribute('aria-selected', on);
        buttons[i].setAttribute('aria-pressed', on);
    }
}

function onWorkbenchSectionOpened(id) {
    if (id === 'bime-search-section') scheduleLibraryWarmup();
    if (id === 'wb-pathway-section') {
        initPathwayCanvas();
        renderPathwayCanvas();
    }
}

function bindWorkbenchActions() {
    if (document.body && document.body.getAttribute('data-wb-actions-bound') === 'true') return;
    if (document.body) document.body.setAttribute('data-wb-actions-bound', 'true');
    // v2.0.17: Escape closes the footer modal.
    // Merge Stage 2: when the focus lens is open, Escape commits the edited
    // structure back onto the node and collapses the overlay — handled first
    // (and short-circuited) so it takes precedence over the modal close.
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Escape' || e.keyCode === 27) {
            if (_lens && _lens.isOpen()) {
                e.preventDefault();
                collapseStructureLens(true);
                return;
            }
            // v2.3.6 (P0-7): Escape also dismisses the mobile Inspector bottom
            // sheet (after the lens, before the footer modal).
            var insp = document.getElementById('wb-inspector-panel');
            if (insp && insp.classList && insp.classList.contains('wb-inspector-open')) {
                e.preventDefault();
                setMobileInspectorOpen(false);
                return;
            }
            var root = document.getElementById('wb-modal-root');
            if (root && !root.hidden) { e.preventDefault(); closeWorkbenchModal(); }
        }
    });
    // v2.3.6 (P0-6): re-evaluate the lens layout when the viewport crosses the
    // mobile breakpoint while a lens is open (e.g. a phone rotates, or a desktop
    // window is resized narrow). positionStructureLens re-decides sheet-vs-anchored
    // each call, so re-running it here is all that's needed.
    if (typeof window !== 'undefined' && window.addEventListener) {
        window.addEventListener('resize', function() {
            if (_lens && _lens.isOpen()) { positionStructureLens(); }
        });
    }
    // v2.0.58: marquee selection mousedown/move/up.
    document.addEventListener('mousedown', handleAtomTraceMouseDown);
    document.addEventListener('mousemove', handleAtomTraceMouseMove);
    document.addEventListener('mouseup', handleAtomTraceMouseUp);
    // v2.0.61: right-click context menu in pathway mode + global
    // dismiss-on-outside-click / dismiss-on-Escape.
    document.addEventListener('contextmenu', handlePathwayCanvasContextMenu);
    document.addEventListener('dblclick', handlePathwayCanvasDoubleClick);
    document.addEventListener('mousedown', handleLensClickAway, true);
    document.addEventListener('mousedown', function(e) {
        if (!_pathwayContextMenu) { return; }
        if (e.target && _pathwayContextMenu.contains(e.target)) { return; }
        hidePathwayContextMenu();
    }, true);
    document.addEventListener('keydown', function(e) {
        if (_pathwayContextMenu && (e.key === 'Escape' || e.keyCode === 27)) {
            hidePathwayContextMenu();
        }
    });
    document.addEventListener('click', function(e) {
        // v2.0.58: a marquee drag just finalised — swallow the trailing click
        // so we don't accidentally clear or re-trigger the selection.
        if (_suppressNextAtomTraceClick) {
            _suppressNextAtomTraceClick = false;
            e.preventDefault();
            return;
        }
        // v2.0.22: atom-trace strip click — handled before the wb-action dispatch.
        if (e.target && e.target.closest && e.target.closest('.bime-trace-atom')) {
            e.preventDefault();
            // v2.0.54: shift-click on a start-metabolite atom toggles the
            // atom into the moiety selection. Plain click resets the
            // moiety and runs single-atom trace as before.
            handleTraceAtomClick(e.target, !!e.shiftKey);
            return;
        }
        // v2.0.54: "Clear moiety" link in the atom-trace status bar.
        if (e.target && e.target.closest &&
                e.target.closest('[data-wb-action="clear-moiety-trace"]')) {
            e.preventDefault();
            clearMoietyTrace();
            return;
        }
        var target = e.target && e.target.closest ? e.target.closest('[data-wb-action]') : null;
        if (!target) return;
        var action = target.getAttribute('data-wb-action');
        if (!action) return;
        e.preventDefault();
        if (action === 'clean2d') clean2DWorkbench();
        else if (action === 'best-fit') bestFitWorkbench();
        else if (action === 'set-workbench-mode') setWorkbenchMode(target.getAttribute('data-mode'));
        else if (action === 'toggle-rs') toggleRSWorkbench();
        else if (action === 'toggle-ez') toggleEZWorkbench();
        else if (action === 'apply-structure-properties') applyStructureProperties();
        else if (action === 'clear-structure-name') clearStructureName();
        else if (action === 'apply-reaction-labels') applyReactionLabels();
        else if (action === 'clear-reaction-labels') clearReactionLabels();
        else if (action === 'apply-property-reaction-labels') applyPropertyReactionLabels();
        else if (action === 'clear-property-reaction-labels') clearPropertyReactionLabels();
        else if (action === 'toggle-workbench-section') toggleWorkbenchSection(target);
        else if (action === 'set-pathway-tool') setPathwayTool(target.getAttribute('data-tool'));
        else if (action === 'add-pathway-current') addCurrentMoleculeToPathway();
        else if (action === 'reaction-to-pathway') convertCurrentReactionToPathway();
        // Hybrid Stage B (v2.4.1): from the Editor view, drop the current drawing
        // onto the pathway (molecule -> node; reaction -> reactant/step/product),
        // then switch to the Pathway view.
        else if (action === 'editor-add-to-pathway') addEditorDrawingToPathway();
        // v2.0.60: load the selected pathway node's structure back into
        // the molecule editor for in-place revision.
        else if (action === 'edit-pathway-structure') editPathwayNodeStructure();
        // v2.0.61: context-menu action handlers. Each dispatches into a
        // small bridge that re-reads the data-* attributes the menu item
        // carries so the action has full context regardless of the
        // global selection state.
        else if (action === 'context-edit-node-structure') {
            hidePathwayContextMenu();
            contextEditNodeStructure(target.getAttribute('data-context-id'));
        }
        else if (action === 'context-update-node-labels') {
            hidePathwayContextMenu();
            contextUpdateNodeLabels(target.getAttribute('data-context-type'),
                                    target.getAttribute('data-context-id'));
        }
        else if (action === 'context-delete-item') {
            hidePathwayContextMenu();
            contextDeleteItem(target.getAttribute('data-context-type'),
                              target.getAttribute('data-context-id'));
        }
        else if (action === 'context-copy-smiles') {
            hidePathwayContextMenu();
            contextCopySmiles(target.getAttribute('data-context-id'));
        }
        else if (action === 'context-duplicate-node') {
            hidePathwayContextMenu();
            contextDuplicateNode(target.getAttribute('data-context-id'));
        }
        else if (action === 'context-reverse-edge') {
            hidePathwayContextMenu();
            contextReverseEdge(target.getAttribute('data-context-id'));
        }
        else if (action === 'context-cycle-arrow-type') {
            hidePathwayContextMenu();
            contextCycleArrowType(target.getAttribute('data-context-id'));
        }
        else if (action === 'context-add-here') {
            hidePathwayContextMenu();
            contextAddHere(target.getAttribute('data-context-kind'),
                           target.getAttribute('data-context-x'),
                           target.getAttribute('data-context-y'));
        }
        else if (action === 'update-pathway-selection') updatePathwaySelectionFromControls();
        else if (action === 'import-pathway-background') openPathwayBackgroundPicker();
        else if (action === 'recognize-pathway-background') recognizePathwayBackground();
        else if (action === 'remove-pathway-background') removePathwayBackground();
        else if (action === 'pathway-zoom') zoomPathway(target.getAttribute('data-dir'));
        else if (action === 'pathway-fit') fitPathwayCanvas();
        else if (action === 'pathway-layout') layoutPathwayCanvas();
        else if (action === 'set-pathway-layout-shape') setPathwayLayoutShape(target.getAttribute('data-shape'));
        else if (action === 'export-pathway-svg') exportPathwaySvg();
        // v2.0.30: pathway JSON round-trip + KEGG KGML import
        else if (action === 'export-pathway-json') exportPathwayJson();
        else if (action === 'import-pathway-json') openPathwayJsonPicker();
        else if (action === 'import-pathway-kgml') openPathwayKgmlPicker();
        // v2.0.32: SBGN-ML + BioPAX export
        else if (action === 'export-pathway-sbgnml') exportPathwaySbgnml();
        else if (action === 'export-pathway-biopax') exportPathwayBiopax();
        // v2.0.37: SBGN-ML + BioPAX import (round-trip closure)
        else if (action === 'import-pathway-sbgnml') openPathwaySbgnmlPicker();
        else if (action === 'import-pathway-biopax') openPathwayBiopaxPicker();
        // v2.0.38: page-format toggle (A4 / Letter / screen)
        else if (action === 'set-pathway-page-format') setPathwayPageFormat(target.getAttribute('data-format'), target.getAttribute('data-orientation'));
        else if (action === 'clear-pathway') clearPathwayCanvas();
        else if (action === 'clear-curly-arrows') {
            // v2.0.29: drop all IUPAC mechanism curly arrows on the current molecule.
            if (editor && editor.molecule && editor.molecule.clearCurlyArrows) {
                editor.saveHistory && editor.saveHistory();
                editor.molecule.clearCurlyArrows();
                editor.render && editor.render();
                if (typeof editor.changed === 'function') { editor.changed(); }
            }
        }
        else if (action === 'load-example-pathway') loadExamplePathway(target.getAttribute('data-id'));
        else if (action === 'open-modal') openWorkbenchModal(target.getAttribute('data-modal'));
        else if (action === 'close-modal') closeWorkbenchModal();
        else if (action === 'export-preset') runExportPreset(target.getAttribute('data-kind'));
        else if (action === 'export') doExport(target.getAttribute('data-format'));
        else if (action === 'output-tab') setOutputTab(target.getAttribute('data-tab'));
        else if (action === 'copy-output') copyOut();
        else if (action === 'validate') doValidate();
        else if (action === 'clear-editor') {
            if (editor && editor.reset) {
                editor.reset();
                queueOutputUpdate();
                scheduleWorkbenchUpdate();
            }
        } else if (action === 'load-named') {
            loadNamed(target.getAttribute('data-smiles') || '', target.getAttribute('data-name') || '');
        } else if (action === 'load') {
            load(target.getAttribute('data-smiles') || '');
        } else if (action === 'load-reaction') {
            loadReactionScheme(target.getAttribute('data-smiles') || '', {
                conditions: target.getAttribute('data-conditions') || '',
                yield: target.getAttribute('data-yield') || '',
                note: target.getAttribute('data-note') || ''
            });
        } else if (action === 'load-smiles-input') {
            loadFromSmilesInput();
        } else if (action === 'paste-clipboard') {
            pasteFromClipboard();
        } else if (action === 'set-tool') {
            if (editor && editor._setTool) editor._setTool(target.getAttribute('data-tool'));
        } else if (action === 'focus-command') {
            focusCommandPalette();
        } else if (action === 'focus-output') {
            setOutputTab('smiles');
            var out = document.getElementById('smiles-out');
            if (out) out.focus();
        } else if (action === 'dock-draw') {
            // v2.3.6 (P0-7): the mobile dock's primary "Draw" action. Routes to the
            // SAME continuum entry point as the Tools-row / command-grid "Draw
            // Structure" (addCurrentMoleculeToPathway -> newPathwayNodeInLens drops a
            // fresh node and opens it in the lens; on a phone that lens is the
            // full-screen sheet). A dedicated dock action keeps add-pathway-current at
            // its two pinned placements while giving the dock the right behaviour.
            addCurrentMoleculeToPathway();
        } else if (action === 'toggle-inspector') {
            // v2.3.6 (P0-7): the mobile dock's "Inspect" toggle reveals/hides the
            // Inspector (Selection / Warnings / Command palette) as a bottom sheet.
            toggleMobileInspector();
        } else if (action === 'close-inspector') {
            // v2.3.6 (P0-7): the X inside the inspector bottom sheet.
            setMobileInspectorOpen(false);
        }
    });
}

// v2.3.6 (P0-6/P0-7): the shared mobile breakpoint. MUST match the CSS
// `@media screen and (max-width: 720px)` block — the lens-sheet layout and the
// inspector bottom sheet are CSS-owned at that width, and the JS only toggles
// the classes / skips the desktop anchoring when the viewport is this narrow.
var WB_MOBILE_BREAKPOINT = 720;

// v2.3.6 (P0-6): true on a narrow (phone) viewport. Prefers matchMedia (so it
// tracks the SAME query the CSS uses) and falls back to window.innerWidth.
// Guarded for the test/DOM-stub environment (no window -> desktop).
function isNarrowViewport() {
  if (typeof window === 'undefined') { return false; }
  if (window.matchMedia) {
    try { return window.matchMedia('(max-width: ' + WB_MOBILE_BREAKPOINT + 'px)').matches; }
    catch (e) {}
  }
  return typeof window.innerWidth === 'number' && window.innerWidth <= WB_MOBILE_BREAKPOINT;
}

// v2.3.6 (P0-7): show/hide the Inspector as a mobile bottom sheet. Toggling
// `.wb-inspector-open` overrides the `display:none` the <=720px block sets, so
// the Selection details, Warnings, and Command palette get a touch entry point.
// Keeps the dock toggle's aria-expanded in sync (WCAG 4.1.2) and, when opening,
// expands the collapsed Command panel + moves focus to the inspector so the
// sheet is reachable; when closing, returns focus to the dock toggle.
function setMobileInspectorOpen(open) {
  var inspector = document.getElementById('wb-inspector-panel');
  if (!inspector || !inspector.classList) { return; }
  if (open) { inspector.classList.add('wb-inspector-open'); }
  else { inspector.classList.remove('wb-inspector-open'); }
  var toggle = document.querySelector('.wb-mobile-dock [data-wb-action="toggle-inspector"]');
  if (toggle && toggle.setAttribute) {
    toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
  }
  if (open) {
    // Surface the Command palette (it boots collapsed) so the essentials are
    // visible the moment the sheet opens, then move focus into the sheet.
    openWorkbenchPanel('command');
    if (inspector.focus) {
      if (!inspector.getAttribute('tabindex')) { inspector.setAttribute('tabindex', '-1'); }
      try { inspector.focus(); } catch (eF) {}
    }
  } else if (toggle && toggle.focus) {
    try { toggle.focus(); } catch (eT) {}
  }
}

// v2.3.6 (P0-7): flip the inspector bottom sheet open/closed from the dock.
function toggleMobileInspector() {
  var inspector = document.getElementById('wb-inspector-panel');
  var isOpen = !!(inspector && inspector.classList && inspector.classList.contains('wb-inspector-open'));
  setMobileInspectorOpen(!isOpen);
}

// v2.3.6 (P0-7): show/hide the mobile dock. Used to get the dock out of the way
// while the full-screen lens sheet owns the screen (the sheet's own sticky
// Done/Cancel are the exits). `.wb-dock-hidden` is keyed in the <=720px block;
// off-mobile the dock is already display:none, so this is a harmless no-op.
function setMobileDockHidden(hidden) {
  var dock = document.querySelector('.wb-mobile-dock');
  if (!dock || !dock.classList) { return; }
  if (hidden) { dock.classList.add('wb-dock-hidden'); }
  else { dock.classList.remove('wb-dock-hidden'); }
}

// v2.4.3 (Hybrid Stage D): hide/show the mobile dock for the standalone Editor
// VIEW. The dock's Draw/Fit/Inspect are Pathway-canvas affordances; in the
// Editor view the editor IS the drawing surface and "Add to pathway" lives in
// the Output row, so the dock is redundant/confusing there. This uses its OWN
// `.wb-dock-editor-hidden` flag, orthogonal to the lens-sheet `.wb-dock-hidden`,
// so the two never fight (collapsing a lens can't reveal the dock mid-Editor-
// view, and leaving the Editor view can't reveal it mid-lens-sheet). Keyed in
// the <=720px block; off-mobile the dock is already display:none — a no-op.
function setMobileDockEditorHidden(hidden) {
  var dock = document.querySelector('.wb-mobile-dock');
  if (!dock || !dock.classList) { return; }
  if (hidden) { dock.classList.add('wb-dock-editor-hidden'); }
  else { dock.classList.remove('wb-dock-editor-hidden'); }
}

function initEditorSurfaceEvents() {
    if (!editor || !editor.renderer || !editor.renderer.svg) return;
    var svg = editor.renderer.svg;
    ['click', 'mouseup', 'keyup', 'touchend'].forEach(function(evName) {
        svg.addEventListener(evName, function() { scheduleWorkbenchUpdate(); });
    });
    svg.addEventListener('dragover', function(e) {
        if (hasTemplateDrag(e)) {
            e.preventDefault();
            e.dataTransfer.dropEffect = 'copy';
        }
    });
    svg.addEventListener('drop', function(e) {
        if (!hasTemplateDrag(e)) return;
        e.preventDefault();
        var name = e.dataTransfer.getData('application/x-bime-template');
        if (!name) return;
        var pos = editor.renderer.screenToMol(e.clientX, e.clientY);
        stampTemplate(name, pos.x, pos.y);
    });
}

function hasTemplateDrag(e) {
    var dt = e && e.dataTransfer;
    if (!dt || !dt.types) return false;
    for (var i = 0; i < dt.types.length; i++) {
        if (dt.types[i] === 'application/x-bime-template') return true;
    }
    return false;
}

/* ---- Pathway canvas: lightweight SVG workspace for annotated metabolic maps ---- */
function initPathwayCanvas() {
    var svg = document.getElementById('wb-pathway-canvas');
    if (!svg || _pathway.initialized) return;
    _pathway.initialized = true;
    svg.addEventListener('pointerdown', onPathwayPointerDown);
    document.addEventListener('pointermove', onPathwayPointerMove);
    document.addEventListener('pointerup', endPathwayDrag);
    svg.addEventListener('keydown', onPathwayKeyDown);
    svg.addEventListener('wheel', onPathwayWheel, { passive: false });
    var bgInput = document.getElementById('wb-pathway-bg-input');
    if (bgInput) {
        bgInput.addEventListener('change', function(e) {
            importPathwayBackgroundFile(e.target && e.target.files ? e.target.files[0] : null);
            bgInput.value = '';
        });
    }
    var shortcut = document.getElementById('wb-pathway-shortcut');
    if (shortcut) {
        shortcut.addEventListener('change', function() {
            applyPathwayShortcutToLabel(shortcut.value);
        });
    }
    setPathwayTool(_pathway.tool || 'select', true);
    renderPathwayCanvas();
}

function setPathwayTool(tool, silent) {
    if (!/^(select|metabolite|residue|arrow|curly|step|note|compartment)$/.test(tool || '')) tool = 'select';
    _pathway.tool = tool;
    if (tool !== 'arrow') _pathway.pendingArrow = null;
    if (tool !== 'curly') _pathway.pendingMechanismArrow = null;
    var buttons = document.querySelectorAll('.wb-pathway-tool');
    for (var i = 0; i < buttons.length; i++) {
        buttons[i].setAttribute('aria-pressed',
            buttons[i].getAttribute('data-tool') === tool ? 'true' : 'false');
    }
    if (!silent) {
        if (tool === 'select') {
            pathwayStatus('Select and drag items. Hold Shift while dragging the background to pan.');
        } else if (tool === 'arrow') {
            pathwayStatus('Click a metabolite node, then click a second node to add an arrow.');
        } else if (tool === 'curly') {
            pathwayStatus('Click a structure or step, then click the electron-flow destination.');
        } else if (tool === 'residue') {
            pathwayStatus('Choose a shortcut or edit the label, then click the canvas.');
        } else if (tool === 'step') {
            pathwayStatus('Click the canvas to add a mechanism step.');
        } else if (tool === 'note') {
            pathwayStatus('Click the canvas to add a comment.');
        } else {
            pathwayStatus('Click the canvas to add a ' + tool + '.');
        }
    }
}

function focusPathwayCanvas() {
    openWorkbenchSection('wb-pathway-section');
    var svg = document.getElementById('wb-pathway-canvas');
    if (svg) {
        svg.focus();
        svg.scrollIntoView({ block: 'center', behavior: 'smooth' });
    }
}

function pathwayStatus(text) {
    var el = document.getElementById('wb-pathway-status');
    if (el) el.textContent = text || '';
}

function sanitizePathwayText(value, maxLen) {
    var max = maxLen || 80;
    return String(value || '')
        .replace(/[\u0000-\u001f\u007f<>]/g, ' ')
        .replace(/\s+/g, ' ')
        .trim()
        .slice(0, max);
}

function pathwayInputLabel(fallback, maxLen) {
    var input = document.getElementById('wb-pathway-label');
    var value = sanitizePathwayText(input && input.value, maxLen || 80);
    if (!value || value === 'Metabolite') value = sanitizePathwayText(fallback, maxLen || 80);
    return value || 'Metabolite';
}

function pathwayShortcutValue() {
    var select = document.getElementById('wb-pathway-shortcut');
    var value = sanitizePathwayText(select && select.value, 18);
    return value === 'custom' ? '' : value;
}

function pathwayShortcutMeta(value) {
    value = sanitizePathwayText(value, 18);
    return Object.prototype.hasOwnProperty.call(PATHWAY_SHORTCUTS, value) ? PATHWAY_SHORTCUTS[value] : null;
}

function applyPathwayShortcutToLabel(value) {
    var meta = pathwayShortcutMeta(value);
    var input = document.getElementById('wb-pathway-label');
    if (!meta || !input) return;
    input.value = value;
}

function pathwayNodeKindForLabel(label, fallback) {
    var meta = pathwayShortcutMeta(label);
    if (meta && meta.kind) return meta.kind;
    return fallback || 'metabolite';
}

function pathwayNodeSubtitleForLabel(label, fallback) {
    var meta = pathwayShortcutMeta(label);
    if (meta && meta.name && meta.name !== label) return meta.name;
    return fallback || '';
}

function updatePathwayControlsFromSelection() {
    var item = findPathwayItem(_pathway.selectedType, _pathway.selectedId);
    if (!item) return;
    var input = document.getElementById('wb-pathway-label');
    if (input && item.label) input.value = item.label;
    var shortcut = document.getElementById('wb-pathway-shortcut');
    if (shortcut) shortcut.value = pathwayShortcutMeta(item.label) ? item.label : 'custom';
    // v2.0.18: reflect the selected edge's reaction-arrow type.
    var edgeArrow = document.getElementById('wb-pathway-edge-arrow');
    if (edgeArrow && _pathway.selectedType === 'edge') {
        var t = item.arrowType;
        edgeArrow.value = (PATHWAY_ARROW_TYPES.indexOf(t) >= 0) ? t : 'forward';
    }
}

function pathwayVisibleCenter() {
    // v2.0.42: routed through CanvasView so the screen→world arithmetic
    // matches the molecule editor's screen→world arithmetic byte-for-byte.
    // The fallback inline form is preserved for runs where the bundle
    // doesn't include CanvasView (e.g. source-mode dev loads).
    if (typeof CanvasView !== 'undefined' && CanvasView.screenToWorld) {
        var size = pathwayPageSize();
        return CanvasView.screenToWorld(_pathway, size.w / 2, size.h / 2);
    }
    return {
        x: (600 - _pathway.offsetX) / _pathway.scale,
        y: (310 - _pathway.offsetY) / _pathway.scale
    };
}

function newPathwayId(prefix) {
    var id = prefix + _pathway.nextId;
    _pathway.nextId += 1;
    return id;
}

function pathwayEventPoint(e) {
    var svg = document.getElementById('wb-pathway-canvas');
    var rect = svg && svg.getBoundingClientRect ? svg.getBoundingClientRect() : null;
    if (!rect || !rect.width || !rect.height) return { x: 0, y: 0 };
    // v2.0.42: viewport-aware screen→world. The DOM rect is converted to
    // the active page-format viewBox coords, then run through
    // CanvasView.screenToWorld for the world coords.
    var size = pathwayPageSize();
    var sx = ((e.clientX - rect.left) / rect.width) * size.w;
    var sy = ((e.clientY - rect.top) / rect.height) * size.h;
    if (typeof CanvasView !== 'undefined' && CanvasView.screenToWorld) {
        return CanvasView.screenToWorld(_pathway, sx, sy);
    }
    return {
        x: (sx - _pathway.offsetX) / _pathway.scale,
        y: (sy - _pathway.offsetY) / _pathway.scale
    };
}

function pathwayScreenDelta(e, drag) {
    var svg = document.getElementById('wb-pathway-canvas');
    var rect = svg && svg.getBoundingClientRect ? svg.getBoundingClientRect() : null;
    if (!rect || !rect.width || !rect.height) return { x: 0, y: 0 };
    // Convert screen-pixel drag deltas to world units using the LIVE page
    // size (matches pathwayEventPoint). Hardcoding 1200/620 was the 'screen'
    // format only, so background pan drifted (wrong rate + wrong aspect) on
    // A4 / Letter, whose viewBox is 1123x794 / 1056x816 etc.
    var size = pathwayPageSize();
    return {
        x: (e.clientX - drag.clientX) * (size.w / rect.width),
        y: (e.clientY - drag.clientY) * (size.h / rect.height)
    };
}

function pathwayEventItem(target) {
    var el = target && target.closest ? target.closest('[data-pathway-id]') : null;
    if (!el) return null;
    return {
        type: el.getAttribute('data-pathway-type'),
        id: el.getAttribute('data-pathway-id')
    };
}

function pathwayCollection(type) {
    if (type === 'background') return _pathway.backgrounds;
    if (type === 'node') return _pathway.nodes;
    if (type === 'note') return _pathway.notes;
    if (type === 'step') return _pathway.steps;
    if (type === 'compartment') return _pathway.compartments;
    if (type === 'edge') return _pathway.edges;
    if (type === 'mechanism') return _pathway.mechanismArrows;
    return [];
}

function findPathwayItem(type, id) {
    var items = pathwayCollection(type);
    for (var i = 0; i < items.length; i++) {
        if (items[i].id === id) return items[i];
    }
    return null;
}

function selectPathwayItem(type, id) {
    _pathway.selectedType = type || '';
    _pathway.selectedId = id || null;
    updatePathwayControlsFromSelection();
}

function onPathwayPointerDown(e) {
    var svg = document.getElementById('wb-pathway-canvas');
    if (!svg || (e.button && e.button !== 0 && e.button !== 1)) return;
    svg.focus();
    var item = pathwayEventItem(e.target);
    var point = pathwayEventPoint(e);
    if (item) {
        if (_pathway.tool === 'arrow' && item.type === 'node') {
            handlePathwayArrowTarget(item.id);
            e.preventDefault();
            return;
        }
        if (_pathway.tool === 'curly' && canUseMechanismAnchor(item.type)) {
            handleMechanismArrowTarget(item.type, item.id);
            e.preventDefault();
            return;
        }
        selectPathwayItem(item.type, item.id);
        if (_pathway.tool === 'select' &&
                canDragPathwayItem(item.type)) {
            var selected = findPathwayItem(item.type, item.id);
            if (selected) {
                // v2.0.47: snapshot the pre-drag state so Ctrl+Z returns
                // the dragged item to its original position.
                if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
                _pathway.drag = {
                    type: 'item',
                    itemType: item.type,
                    id: item.id,
                    startX: point.x,
                    startY: point.y,
                    origX: selected.x,
                    origY: selected.y
                };
                svg.classList.add('is-dragging');
            }
        }
        renderPathwayCanvas();
        e.preventDefault();
        return;
    }
    if (_pathway.tool === 'metabolite') {
        var metaboliteLabel = pathwayInputLabel(pathwayShortcutValue() || 'Metabolite');
        addPathwayNode(point.x, point.y, {
            label: metaboliteLabel,
            kind: pathwayNodeKindForLabel(metaboliteLabel, 'metabolite'),
            subtitle: pathwayNodeSubtitleForLabel(metaboliteLabel, '')
        });
    } else if (_pathway.tool === 'residue') {
        var residueLabel = pathwayInputLabel(pathwayShortcutValue() || 'Gly');
        addPathwayNode(point.x, point.y, {
            label: residueLabel,
            kind: pathwayNodeKindForLabel(residueLabel, 'residue'),
            subtitle: pathwayNodeSubtitleForLabel(residueLabel, '')
        });
    } else if (_pathway.tool === 'step') {
        addPathwayStep(point.x, point.y, pathwayInputLabel('Mechanism step'));
    } else if (_pathway.tool === 'note') {
        addPathwayNote(point.x, point.y, pathwayInputLabel('Comment'));
    } else if (_pathway.tool === 'compartment') {
        addPathwayCompartment(point.x, point.y, pathwayInputLabel('Compartment'));
    } else if (_pathway.tool === 'select' && (e.shiftKey || e.button === 1)) {
        _pathway.drag = {
            type: 'pan',
            clientX: e.clientX,
            clientY: e.clientY,
            offsetX: _pathway.offsetX,
            offsetY: _pathway.offsetY
        };
        svg.classList.add('is-dragging');
    } else {
        selectPathwayItem('', null);
        if (_pathway.tool === 'arrow') {
            _pathway.pendingArrow = null;
            pathwayStatus('Click a metabolite node to start an arrow.');
        } else if (_pathway.tool === 'curly') {
            _pathway.pendingMechanismArrow = null;
            pathwayStatus('Click a structure or step to start an electron-flow arrow.');
        }
        renderPathwayCanvas();
    }
    e.preventDefault();
}

function onPathwayPointerMove(e) {
    var drag = _pathway.drag;
    if (!drag) return;
    if (drag.type === 'pan') {
        var delta = pathwayScreenDelta(e, drag);
        _pathway.offsetX = drag.offsetX + delta.x;
        _pathway.offsetY = drag.offsetY + delta.y;
        schedulePathwayRender();
        e.preventDefault();
        return;
    }
    var item = findPathwayItem(drag.itemType, drag.id);
    if (!item) return;
    var point = pathwayEventPoint(e);
    item.x = drag.origX + (point.x - drag.startX);
    item.y = drag.origY + (point.y - drag.startY);
    schedulePathwayRender();
    e.preventDefault();
}

function endPathwayDrag() {
    if (!_pathway.drag) return;
    _pathway.drag = null;
    var svg = document.getElementById('wb-pathway-canvas');
    if (svg) svg.classList.remove('is-dragging');
}

// v2.0.43: PathwayHistory instance. Lazy-instantiated on first edit so
// the cost is zero for users who never touch the pathway canvas. Cleared
// on clearPathwayCanvas / loadExamplePathway / file-import.
var _pathwayHistory = null;
function pathwayHistory() {
    if (!_pathwayHistory) {
        if (typeof PathwayHistory === 'undefined') { return null; }
        _pathwayHistory = new PathwayHistory({ maxDepth: 100 });
    }
    return _pathwayHistory;
}

function pushPathwayHistory() {
    var h = pathwayHistory();
    if (h) { h.push(_pathway); }
}

function clearPathwayHistory() {
    if (_pathwayHistory) { _pathwayHistory.clear(); }
}

function undoPathway() {
    var h = pathwayHistory();
    if (!h || !h.canUndo()) {
        pathwayStatus('Nothing to undo.');
        return;
    }
    var prev = h.undo();
    if (prev) { applyPathwayStateFromHistory(prev); }
    pathwayStatus('Undo (' + h.depth() + ' snapshot' + (h.depth() === 1 ? '' : 's') + ' remaining).');
}

function redoPathway() {
    var h = pathwayHistory();
    if (!h || !h.canRedo()) {
        pathwayStatus('Nothing to redo.');
        return;
    }
    var next = h.redo();
    if (next) { applyPathwayStateFromHistory(next); }
    pathwayStatus('Redo.');
}

// Apply a snapshot to the live _pathway without disturbing the page-
// format toggle or the history itself. Mirrors
// applyPathwayStateFromImport but does NOT push a new snapshot (that
// would corrupt the undo stack).
function applyPathwayStateFromHistory(snapshot) {
    if (!snapshot) { return; }
    _pathway.nodes = (snapshot.nodes || []).slice();
    _pathway.edges = (snapshot.edges || []).slice();
    _pathway.notes = (snapshot.notes || []).slice();
    _pathway.steps = (snapshot.steps || []).slice();
    _pathway.compartments = (snapshot.compartments || []).slice();
    _pathway.mechanismArrows = (snapshot.mechanismArrows || []).slice();
    if (typeof snapshot.scale === 'number')   { _pathway.scale = snapshot.scale; }
    if (typeof snapshot.offsetX === 'number') { _pathway.offsetX = snapshot.offsetX; }
    if (typeof snapshot.offsetY === 'number') { _pathway.offsetY = snapshot.offsetY; }
    _pathway.selectedType = '';
    _pathway.selectedId = null;
    _pathway.pendingArrow = null;
    _pathway.pendingMechanismArrow = null;
    renderPathwayCanvas();
}

function onPathwayKeyDown(e) {
    // v2.0.43: Ctrl+Z (undo), Ctrl+Shift+Z / Ctrl+Y (redo) drive
    // PathwayHistory when the pathway canvas has focus.
    //
    // Merge Stage 4 (v2.3.2): focus-scoped undo/redo. While a focus lens is
    // OPEN the molecule editor owns Ctrl+Z/Y at atom/bond granularity (its
    // ephemeral focus-session history); the pathway timeline is paused so the
    // per-atom edits do NOT leak into PathwayHistory. The whole session lands
    // as exactly one pathway snapshot on commit (collapseStructureLens).
    //
    // Arbitration vs MolEditor's OWN document-level keydown (editor/MolEditor.js
    // ~1632, gated by isInteractiveTarget):
    //   - This handler is bound on #wb-pathway-canvas (the SVG); MolEditor's is
    //     bound on `document`. The editor overlay (.wb-editor-wrap > #bime-editor)
    //     is a SEPARATE DOM subtree, NOT a descendant of the pathway SVG.
    //   - (a) Focus inside the editor overlay: the event never passes through the
    //     pathway SVG, so THIS handler does not fire; MolEditor's document handler
    //     runs _doAction once. Single application.
    //   - (b) Focus on the pathway SVG (lens open): THIS handler fires first
    //     (deeper node, bubble phase). We route to the editor AND stopPropagation
    //     so MolEditor's document handler never sees the event. Single application.
    // Net rule: exactly one undo per Ctrl+Z in every case; pathway undo is
    // suppressed entirely while the lens is open.
    var ctrlOrMeta = e.ctrlKey || e.metaKey;
    if (ctrlOrMeta && (e.key === 'z' || e.key === 'Z' || e.key === 'y' || e.key === 'Y')) {
        if (_lens && _lens.isOpen()) {
            if (editor && editor._doAction) {
                var redoZ = (e.key === 'z' || e.key === 'Z') ? e.shiftKey : true;
                editor._doAction(redoZ ? 'redo' : 'undo');
            }
            e.preventDefault();
            e.stopPropagation();
            return;
        }
    }
    if (ctrlOrMeta && (e.key === 'z' || e.key === 'Z')) {
        if (e.shiftKey) { redoPathway(); } else { undoPathway(); }
        e.preventDefault();
        return;
    }
    if (ctrlOrMeta && (e.key === 'y' || e.key === 'Y')) {
        redoPathway();
        e.preventDefault();
        return;
    }
    // Merge Stage 5 (v2.3.3): Enter-to-focus parity. With NO lens open and a
    // single structure NODE selected, Enter opens the focus lens on it — the
    // keyboard twin of double-click (handlePathwayCanvasDoubleClick). Guarded so
    // it never fires while the lens is already open (the Ctrl+Z/Y lens block
    // above has already handled the focused-session chords).
    if (e.key === 'Enter' && !(_lens && _lens.isOpen())
            && _pathway.selectedType === 'node' && _pathway.selectedId) {
        if (openStructureLens(_pathway.selectedId)) {
            e.preventDefault();
            return;
        }
    }
    if (e.key === 'Delete' || e.key === 'Backspace') {
        pushPathwayHistory();
        deletePathwaySelection();
        e.preventDefault();
    } else if (e.key === 'Escape') {
        _pathway.pendingArrow = null;
        _pathway.pendingMechanismArrow = null;
        selectPathwayItem('', null);
        renderPathwayCanvas();
        pathwayStatus('Selection cleared.');
        e.preventDefault();
    } else if (e.key === '+' || e.key === '=') {
        zoomPathway('in');
        e.preventDefault();
    } else if (e.key === '-' || e.key === '_') {
        zoomPathway('out');
        e.preventDefault();
    } else if (/^Arrow/.test(e.key || '')) {
        pushPathwayHistory();
        nudgePathwaySelection(e.key, e.shiftKey ? 40 : 10);
        e.preventDefault();
    }
}

function onPathwayWheel(e) {
    if (!e.ctrlKey && !e.metaKey) return;
    e.preventDefault();
    zoomPathway(e.deltaY < 0 ? 'in' : 'out');
}

function addPathwayNode(x, y, opts) {
    opts = opts || {};
    // v2.0.43: snapshot the pre-edit state so Ctrl+Z can restore it.
    // We push BEFORE the mutation so undo returns to the state without
    // this node, not the state including it. Silent batch operations push
    // one shared history entry before inserting several related nodes.
    if (!opts.silent && !opts.batch && typeof pushPathwayHistory === 'function') {
        pushPathwayHistory();
    }
    var label = sanitizePathwayText(opts.label || 'Metabolite', 80);
    var kind = sanitizePathwayText(opts.kind || pathwayNodeKindForLabel(label, 'metabolite'), 24) || 'metabolite';
    var subtitle = sanitizePathwayText(opts.subtitle || pathwayNodeSubtitleForLabel(label, ''), 40);
    var structure = normalizePathwayStructurePayload(opts.structure);
    if (structure && structure.type === 'reaction') { kind = 'reaction'; }
    var compact = !opts.imageHref && (kind === 'residue' || kind === 'cofactor');
    var imageW = kind === 'reaction' ? 240 : 180;
    var imageH = kind === 'reaction' ? 124 : 112;
    var node = {
        id: newPathwayId('node-'),
        x: Math.round(x - (opts.imageHref ? imageW / 2 : (compact ? 52 : 82))),
        y: Math.round(y - (opts.imageHref ? imageH / 2 : (compact ? 34 : 46))),
        w: opts.imageHref ? imageW : (compact ? 104 : 164),
        h: opts.imageHref ? imageH : (compact ? 68 : 92),
        label: label,
        kind: kind,
        subtitle: subtitle,
        imageHref: opts.imageHref || '',
        structure: structure
    };
    _pathway.nodes.push(node);
    if (!opts.silent) {
        selectPathwayItem('node', node.id);
        renderPathwayCanvas();
        pathwayStatus((kind === 'reaction' ? 'Reaction' : (kind === 'residue' ? 'Residue' : (kind === 'cofactor' ? 'Cofactor' : 'Metabolite'))) + ' added. Drag it to position it.');
    }
    return node;
}

function updatePathwaySelectionFromControls() {
    if (!_pathway.selectedType || !_pathway.selectedId) {
        pathwayStatus('Select a pathway item to update.');
        return;
    }
    var item = findPathwayItem(_pathway.selectedType, _pathway.selectedId);
    if (!item) return;
    // v2.0.47: snapshot the pre-edit state so the property-panel
    // "Apply" button is undoable.
    if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var label = pathwayInputLabel(item.label || 'Metabolite', 90);
    if (item.label !== undefined) item.label = label;
    if (_pathway.selectedType === 'node') {
        item.kind = pathwayNodeKindForLabel(label, item.kind || 'metabolite');
        if (!item.imageHref) {
            item.subtitle = pathwayNodeSubtitleForLabel(label, item.subtitle || '');
            if (item.kind === 'residue' || item.kind === 'cofactor') {
                item.w = Math.max(92, Math.min(150, 56 + label.length * 8));
                item.h = 68;
            }
        }
    } else if (_pathway.selectedType === 'step') {
        item.comment = sanitizePathwayText(item.comment || 'Residues, cofactors, evidence, or transition-state note', 90);
    } else if (_pathway.selectedType === 'edge' || _pathway.selectedType === 'mechanism') {
        item.label = label === 'Metabolite' ? '' : label;
        if (_pathway.selectedType === 'edge') {
            // v2.0.18: apply the chosen reaction-arrow type to the edge.
            var edgeArrow = document.getElementById('wb-pathway-edge-arrow');
            if (edgeArrow && PATHWAY_ARROW_TYPES.indexOf(edgeArrow.value) >= 0) {
                item.arrowType = edgeArrow.value;
            }
        }
    }
    renderPathwayCanvas();
    pathwayStatus('Selected pathway item updated.');
}

function addPathwayNote(x, y, label, opts) {
    opts = opts || {};
    // v2.0.47: snapshot pre-edit state for Ctrl+Z (skip on silent batch
    // imports so loadExamplePathway / file-import don't fragment history).
    if (!opts.silent && typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var note = {
        id: newPathwayId('note-'),
        x: Math.round(x - 86),
        y: Math.round(y - 34),
        w: 172,
        h: 68,
        label: sanitizePathwayText(label || 'Comment', 90)
    };
    _pathway.notes.push(note);
    if (!opts.silent) {
        selectPathwayItem('note', note.id);
        renderPathwayCanvas();
        pathwayStatus('Comment added.');
    }
    return note;
}

function addPathwayStep(x, y, label, opts) {
    opts = opts || {};
    if (!opts.silent && typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var step = {
        id: newPathwayId('step-'),
        number: _pathway.steps.length + 1,
        x: Math.round(x - 120),
        y: Math.round(y - 45),
        w: 240,
        h: 90,
        label: sanitizePathwayText(label || 'Mechanism step', 88),
        comment: sanitizePathwayText('Residues, cofactors, evidence, or transition-state note', 90)
    };
    _pathway.steps.push(step);
    if (!opts.silent) {
        selectPathwayItem('step', step.id);
        renderPathwayCanvas();
        pathwayStatus('Mechanism step added.');
    }
    return step;
}

function addPathwayCompartment(x, y, label, opts) {
    opts = opts || {};
    if (!opts.silent && typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var compartment = {
        id: newPathwayId('compartment-'),
        x: Math.round(x - 190),
        y: Math.round(y - 105),
        w: 380,
        h: 210,
        label: sanitizePathwayText(label || 'Compartment', 70)
    };
    _pathway.compartments.push(compartment);
    if (!opts.silent) {
        selectPathwayItem('compartment', compartment.id);
        renderPathwayCanvas();
        pathwayStatus('Compartment added.');
    }
    return compartment;
}

function openPathwayBackgroundPicker() {
    openWorkbenchSection('wb-pathway-section');
    var input = document.getElementById('wb-pathway-bg-input');
    if (input) input.click();
}

function importPathwayBackgroundFile(file) {
    if (!file) return;
    var name = sanitizePathwayText(file.name || 'reference', 80);
    var type = String(file.type || '').toLowerCase();
    var isPdf = type === 'application/pdf' || /\.pdf$/i.test(name);
    var isImage = /^image\/(png|jpe?g|webp|gif)$/.test(type) || /\.(png|jpe?g|webp|gif)$/i.test(name);
    if (!isPdf && !isImage) {
        pathwayStatus('Unsupported background file. Use PNG, JPG, WEBP, GIF, or PDF.');
        return;
    }
    var limit = isPdf ? PATHWAY_BACKGROUND_PDF_LIMIT : PATHWAY_BACKGROUND_IMAGE_LIMIT;
    if (file.size > limit) {
        pathwayStatus('Background file is too large for the local workbench import.');
        return;
    }
    if (isPdf) {
        addPathwayPdfReference(name, file.size);
        return;
    }
    var reader = new FileReader();
    reader.onload = function(ev) {
        var href = String(ev.target && ev.target.result || '');
        if (!/^data:image\/(png|jpe?g|webp|gif);base64,/i.test(href)) {
            pathwayStatus('Could not read image background safely.');
            return;
        }
        addPathwayImageBackground(name, href);
    };
    reader.onerror = function() {
        pathwayStatus('Could not read image background.');
    };
    reader.readAsDataURL(file);
}

function addPathwayImageBackground(name, href) {
    var img = new Image();
    img.onload = function() {
        var maxW = 980;
        var maxH = 520;
        var w = img.naturalWidth || 900;
        var h = img.naturalHeight || 500;
        var scale = Math.min(maxW / Math.max(1, w), maxH / Math.max(1, h), 1);
        w = Math.max(180, Math.round(w * scale));
        h = Math.max(120, Math.round(h * scale));
        var bg = {
            id: newPathwayId('bg-'),
            kind: 'image',
            name: name,
            href: href,
            x: Math.round(600 - w / 2),
            y: Math.round(310 - h / 2),
            w: w,
            h: h,
            meta: 'Traceable reaction/pathway image'
        };
        _pathway.backgrounds.push(bg);
        selectPathwayItem('background', bg.id);
        renderPathwayCanvas();
        pathwayStatus('Image background imported. Trace, annotate, or run a local recognizer.');
    };
    img.onerror = function() {
        pathwayStatus('Could not decode image background.');
    };
    img.src = href;
}

function addPathwayPdfReference(name, size) {
    var bg = {
        id: newPathwayId('bg-'),
        kind: 'pdf',
        name: name,
        href: '',
        x: 330,
        y: 210,
        w: 540,
        h: 190,
        meta: 'PDF reference imported locally (' + formatFileSize(size) + ')'
    };
    _pathway.backgrounds.push(bg);
    selectPathwayItem('background', bg.id);
    renderPathwayCanvas();
    pathwayStatus('PDF reference added. Convert a page to image for visual tracing.');
}

function formatFileSize(bytes) {
    bytes = Math.max(0, bytes || 0);
    if (bytes >= 1024 * 1024) return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
    if (bytes >= 1024) return Math.round(bytes / 1024) + ' KB';
    return bytes + ' B';
}

function selectedPathwayBackground() {
    if (_pathway.selectedType === 'background') {
        return findPathwayItem('background', _pathway.selectedId);
    }
    return _pathway.backgrounds.length ? _pathway.backgrounds[_pathway.backgrounds.length - 1] : null;
}

function removePathwayBackground() {
    var bg = selectedPathwayBackground();
    if (!bg) {
        pathwayStatus('No image or PDF background to remove.');
        return;
    }
    for (var i = _pathway.backgrounds.length - 1; i >= 0; i--) {
        if (_pathway.backgrounds[i].id === bg.id) _pathway.backgrounds.splice(i, 1);
    }
    if (_pathway.selectedType === 'background' && _pathway.selectedId === bg.id) {
        selectPathwayItem('', null);
    }
    renderPathwayCanvas();
    pathwayStatus('Background removed.');
}

function recognizePathwayBackground() {
    var bg = selectedPathwayBackground();
    if (!bg) {
        pathwayStatus('Import an image background before recognition.');
        return;
    }
    if (bg.kind !== 'image') {
        pathwayStatus('Recognition expects an image background. Convert the PDF page to an image first.');
        return;
    }
    var recognizer = window.BIME_REACTION_IMAGE_RECOGNIZER;
    if (typeof recognizer !== 'function') {
        pathwayStatus('No local reaction-image recognizer is registered. The image is ready for manual tracing.');
        return;
    }
    pathwayStatus('Recognizing reaction image locally...');
    setTimeout(function() {
        Promise.resolve().then(function() {
            return recognizer({
                name: bg.name,
                href: bg.href,
                width: bg.w,
                height: bg.h
            });
        }).then(function(result) {
            applyPathwayRecognitionResult(result);
        }).catch(function(err) {
            pathwayStatus('Recognition failed: ' + sanitizePathwayText(err && err.message ? err.message : err, 90));
        });
    }, 0);
}

function applyPathwayRecognitionResult(result) {
    if (result && result.pathway) {
        applyPathwayDiagramResult(result.pathway);
        return;
    }
    var text = '';
    if (typeof result === 'string') text = result;
    else if (result && result.rxn) text = result.rxn;
    else if (result && result.mol) text = result.mol;
    else if (result && result.smiles) text = result.smiles;
    text = String(text || '').trim();
    if (!text) {
        pathwayStatus('Recognizer returned no MOL, RXN, or SMILES text.');
        return;
    }
    try {
        editor.readGenericMolecularInput(text);
        addCurrentMoleculeToPathway();
        pathwayStatus('Recognized structure loaded into the editor and pathway.');
    } catch (e) {
        pathwayStatus('Recognizer output could not be loaded.');
    }
}

function applyPathwayDiagramResult(diagram) {
    if (!diagram || typeof diagram !== 'object') {
        pathwayStatus('Recognizer returned no pathway diagram.');
        return;
    }
    var nodeMap = {};
    var itemMap = {};
    var importedNodes = diagram.nodes || diagram.metabolites || [];
    for (var i = 0; i < importedNodes.length; i++) {
        var src = importedNodes[i] || {};
        var label = sanitizePathwayText(src.label || src.name || src.shortName || 'Metabolite', 80);
        var node = addPathwayNode(
            safePathwayNumber(src.x, 160 + i * 190),
            safePathwayNumber(src.y, 250),
            {
                label: label,
                subtitle: sanitizePathwayText(src.subtitle || src.fullName || pathwayNodeSubtitleForLabel(label, ''), 40),
                kind: sanitizePathwayText(src.kind || pathwayNodeKindForLabel(label, 'metabolite'), 24),
                silent: true
            });
        nodeMap[src.id || src.key || label] = node.id;
        itemMap[src.id || src.key || label] = { type: 'node', id: node.id };
    }
    var importedCompartments = diagram.compartments || [];
    for (var c = 0; c < importedCompartments.length; c++) {
        var compSrc = importedCompartments[c] || {};
        var comp = addPathwayCompartment(
            safePathwayNumber(compSrc.x, 260 + c * 420),
            safePathwayNumber(compSrc.y, 260),
            sanitizePathwayText(compSrc.label || compSrc.name || 'Compartment', 70),
            { silent: true });
        if (compSrc.w) comp.w = Math.max(180, safePathwayNumber(compSrc.w, comp.w));
        if (compSrc.h) comp.h = Math.max(120, safePathwayNumber(compSrc.h, comp.h));
        itemMap[compSrc.id || compSrc.key || comp.label] = { type: 'compartment', id: comp.id };
    }
    var importedSteps = diagram.steps || [];
    for (var s = 0; s < importedSteps.length; s++) {
        var stepSrc = importedSteps[s] || {};
        var step = addPathwayStep(
            safePathwayNumber(stepSrc.x, 220 + s * 260),
            safePathwayNumber(stepSrc.y, 100),
            sanitizePathwayText(stepSrc.label || stepSrc.name || 'Mechanism step', 88),
            { silent: true });
        if (stepSrc.comment || stepSrc.note) step.comment = sanitizePathwayText(stepSrc.comment || stepSrc.note, 90);
        itemMap[stepSrc.id || stepSrc.key || step.label] = { type: 'step', id: step.id };
    }
    var importedNotes = diagram.notes || diagram.comments || [];
    for (var n = 0; n < importedNotes.length; n++) {
        var noteSrc = importedNotes[n] || {};
        var note = addPathwayNote(
            safePathwayNumber(noteSrc.x, 140 + n * 240),
            safePathwayNumber(noteSrc.y, 455),
            sanitizePathwayText(noteSrc.label || noteSrc.text || 'Comment', 90),
            { silent: true });
        itemMap[noteSrc.id || noteSrc.key || note.label] = { type: 'note', id: note.id };
    }
    var importedEdges = diagram.edges || diagram.reactions || [];
    for (var e = 0; e < importedEdges.length; e++) {
        var edgeSrc = importedEdges[e] || {};
        var fromId = nodeMap[edgeSrc.from] || nodeMap[edgeSrc.source];
        var toId = nodeMap[edgeSrc.to] || nodeMap[edgeSrc.target];
        if (fromId && toId && fromId !== toId) {
            _pathway.edges.push({
                id: newPathwayId('edge-'),
                from: fromId,
                to: toId,
                label: sanitizePathwayText(edgeSrc.label || edgeSrc.name || '', 48),
                // v2.0.18: carry the reaction-arrow type through import; a
                // `reversible` flag (as KGML/SBML reactions provide) maps to the
                // equilibrium arrow, otherwise honour an explicit arrowType.
                arrowType: (PATHWAY_ARROW_TYPES.indexOf(edgeSrc.arrowType) >= 0)
                    ? edgeSrc.arrowType
                    : (edgeSrc.reversible === true ? 'reversible' : 'forward')
            });
        }
    }
    var importedMechanisms = diagram.mechanismArrows || diagram.mechanisms || diagram.curlyArrows || [];
    for (var m = 0; m < importedMechanisms.length; m++) {
        var mechSrc = importedMechanisms[m] || {};
        var fromItem = itemMap[mechSrc.from] || itemMap[mechSrc.source];
        var toItem = itemMap[mechSrc.to] || itemMap[mechSrc.target];
        if (fromItem && toItem &&
                canUseMechanismAnchor(fromItem.type) &&
                canUseMechanismAnchor(toItem.type) &&
                !(fromItem.type === toItem.type && fromItem.id === toItem.id)) {
            _pathway.mechanismArrows.push({
                id: newPathwayId('mech-'),
                fromType: fromItem.type,
                fromId: fromItem.id,
                toType: toItem.type,
                toId: toItem.id,
                kind: mechSrc.kind === 'single' ? 'single' : 'pair',
                label: sanitizePathwayText(mechSrc.label || mechSrc.name || '', 48)
            });
        }
    }
    renderPathwayCanvas();
    fitPathwayCanvas();
    pathwayStatus('Recognized pathway diagram imported for editing.');
}

function safePathwayNumber(value, fallback) {
    value = Number(value);
    return isFinite(value) ? value : fallback;
}

function pathwayStructureText(value) {
    value = String(value || '');
    if (value.length > PATHWAY_STRUCTURE_TEXT_LIMIT) {
        return value.slice(0, PATHWAY_STRUCTURE_TEXT_LIMIT);
    }
    return value;
}

function normalizePathwayStructurePayload(value) {
    value = value || {};
    var out = {
        type: value.type === 'reaction' ? 'reaction' : 'molecule',
        smiles: pathwayStructureText(value.smiles),
        mol: pathwayStructureText(value.mol),
        rxn: pathwayStructureText(value.rxn),
        title: sanitizePathwayText(value.title || '', 90),
        comment: sanitizePathwayText(value.comment || '', 180)
    };
    if (!out.smiles && !out.mol && !out.rxn && !out.title && !out.comment) return null;
    return out;
}

function pathwayNodeStructurePayload(node) {
    if (!node) return null;
    if (node.structure && typeof node.structure === 'object') {
        return normalizePathwayStructurePayload(node.structure);
    }
    return normalizePathwayStructurePayload({
        type: node.structureType,
        smiles: node.structureSmiles,
        mol: node.structureMolfile,
        rxn: node.structureRxnfile,
        title: node.structureTitle,
        comment: node.structureComment
    });
}

function pathwayNodeStructureInput(node) {
    var structure = pathwayNodeStructurePayload(node);
    if (structure) {
        if (structure.type === 'reaction' && structure.rxn) return structure.rxn;
        if (structure.mol) return structure.mol;
        if (structure.smiles) return structure.smiles;
        if (structure.rxn) return structure.rxn;
    }
    return '';
}

function genericPathwayStructureLabel(value) {
    value = sanitizePathwayText(value || '', 90);
    return !value || value === 'Metabolite' || value === 'Structure' || value === 'Reaction';
}

function captureEditorStructureForPathway() {
    if (!editor || !editor.molecule || !editor.molecule.atoms || editor.molecule.atoms.length === 0) {
        return null;
    }
    var mol = editor.molecule;
    var isReaction = !!mol.reactionArrow;
    // Label fallback chain. A reaction keeps its named-or-'Reaction' fallback
    // unchanged. For a PLAIN MOLECULE with no editor/mol name, derive the
    // molecular formula (e.g. an unnamed ethanol lands as 'C2H6O' rather than
    // the generic 'Structure') via the existing moleculeFormulaLabel helper,
    // and keep 'Structure' only as the final fallback if no formula is derivable.
    // v2.4.3 (Hybrid Stage D): folds the v2.4.1 unnamed-molecule auto-tag nit.
    var name = (editor.getMolName && editor.getMolName()) || mol.name || '';
    var title;
    if (isReaction) {
        title = name || 'Reaction';
    } else {
        title = name || moleculeFormulaLabel(mol, '') || 'Structure';
    }
    var label = pathwayInputLabel(title, 90);
    var smiles = '';
    var molfile = '';
    var rxnfile = '';
    try { smiles = editor.smiles ? editor.smiles() : ''; } catch (e) { smiles = ''; }
    try { molfile = editor.molFile ? editor.molFile(false) : ''; } catch (e2) { molfile = ''; }
    if (isReaction) {
        try { rxnfile = editor.rxnFile ? editor.rxnFile(false) : ''; } catch (e3) { rxnfile = ''; }
    }
    var imageHref = '';
    if (typeof ImageExport !== 'undefined' && ImageExport.toSVG) {
        try {
            var svg = ImageExport.toSVG(mol, {
                width: isReaction ? 240 : 180,
                height: isReaction ? 88 : 82,
                padding: 8,
                background: 'transparent',
                showAtomLabels: true
            });
            imageHref = pathwaySvgDataUri(svg);
        } catch (e4) {
            imageHref = '';
        }
    }
    var structure = normalizePathwayStructurePayload({
        type: isReaction ? 'reaction' : 'molecule',
        smiles: smiles,
        mol: molfile,
        rxn: rxnfile,
        title: title,
        comment: mol.comment || ''
    });
    return {
        label: label || (isReaction ? 'Reaction' : 'Structure'),
        subtitle: sanitizePathwayText(smiles || (isReaction ? 'Reaction' : ''), 42),
        kind: isReaction ? 'reaction' : pathwayNodeKindForLabel(label, 'metabolite'),
        imageHref: imageHref,
        structure: structure
    };
}

function moleculeFormulaLabel(mol, fallback) {
    if (!mol || typeof Molecule === 'undefined' || !Molecule.formulaString) return fallback || '';
    try {
        var counts = mol.elementCounts ? mol.elementCounts() : null;
        var formula = counts ? Molecule.formulaString(counts) : '';
        return formula || fallback || '';
    } catch (e) {
        return fallback || '';
    }
}

function pathwayImageHrefForMolecule(mol, width, height) {
    if (!mol || typeof ImageExport === 'undefined' || !ImageExport.toSVG) return '';
    try {
        var svg = ImageExport.toSVG(mol, {
            width: width || 180,
            height: height || 96,
            padding: 8,
            background: 'transparent',
            showAtomLabels: true
        });
        return pathwaySvgDataUri(svg);
    } catch (e) {
        return '';
    }
}

function capturePathwayComponentPayload(mol, role, index) {
    var title = role + ' ' + index;
    var formula = moleculeFormulaLabel(mol, title);
    var smiles = '';
    var molfile = '';
    try {
        smiles = (typeof SmilesWriter !== 'undefined' && SmilesWriter.write) ?
            SmilesWriter.write(mol) : '';
    } catch (e) { smiles = ''; }
    try {
        molfile = (typeof MolfileWriter !== 'undefined' && MolfileWriter.write) ?
            MolfileWriter.write(mol, false) : '';
    } catch (e2) { molfile = ''; }
    return {
        label: sanitizePathwayText(formula || title, 80),
        subtitle: sanitizePathwayText(smiles || title, 40),
        kind: 'metabolite',
        imageHref: pathwayImageHrefForMolecule(mol, 172, 96),
        structure: normalizePathwayStructurePayload({
            type: 'molecule',
            smiles: smiles,
            mol: molfile,
            rxn: '',
            title: title,
            comment: role + ' component from drawn reaction'
        })
    };
}

function reactionPathwayArrowType(arrow) {
    if (!arrow) return 'forward';
    return (arrow.type === 'reverse' || arrow.type === 'retro') ? 'reverse' : 'forward';
}

// v2.4.2 (Hybrid Stage C): version stamped INSIDE every edge.aam we write, so a
// future mapping-engine revision can detect + recompute stale linkages without
// guessing. Bumped only when the stored mapping's MEANING changes (not on every
// release). Distinct from the module release stamp.
var AAM_LINK_VERSION = '2.4.2';

function addPathwayEdgeRecord(fromId, toId, label, arrowType, aam) {
    if (!fromId || !toId || fromId === toId) return null;
    var edge = {
        id: newPathwayId('edge-'),
        from: fromId,
        to: toId,
        label: sanitizePathwayText(label || '', 48),
        arrowType: PATHWAY_ARROW_TYPES.indexOf(arrowType) >= 0 ? arrowType : 'forward'
    };
    // v2.4.2 (Hybrid Stage C): optional atom-atom-mapping payload. ADDITIVE —
    // existing callers pass nothing, so the edge is byte-identical to v2.4.1.
    // When present, `aam` carries { mapping, score, status, version } where
    // `mapping` is a fromToIndex-style POSITIONAL per-SMILES index map (the
    // same space AtomTrace.tracePathway consumes), NOT RDT internal atom-ids.
    // linkImportedStepAAM is the writer; PathwayIO round-trips it.
    if (aam && typeof aam === 'object' && aam.mapping && typeof aam.mapping === 'object') {
        edge.aam = {
            mapping: aam.mapping,
            score: (typeof aam.score === 'number') ? aam.score : null,
            status: (typeof aam.status === 'string') ? aam.status : '',
            version: (typeof aam.version === 'string') ? aam.version : AAM_LINK_VERSION
        };
    }
    _pathway.edges.push(edge);
    return edge;
}

function convertCurrentReactionToPathway() {
    if (!editor || !editor.molecule || !editor.molecule.reactionArrow) {
        if (typeof setWorkbenchMode === 'function') { setWorkbenchMode('molecule'); }
        pathwayStatus('Draw or load a reaction arrow first, then convert it to a pathway step.');
        return;
    }
    if (typeof RDT === 'undefined' || !RDT._splitReactionSides) {
        pathwayStatus('Reaction splitting is unavailable in this build.');
        return;
    }
    var split = RDT._splitReactionSides(editor.molecule);
    var reactants = split.reactants || [];
    var products = split.products || [];
    if (!reactants.length || !products.length) {
        pathwayStatus('Reaction conversion needs at least one reactant and one product.');
        return;
    }
    var reactionPayload = captureEditorStructureForPathway();
    if (!reactionPayload || !reactionPayload.structure) {
        pathwayStatus('Could not capture the current reaction for the pathway canvas.');
        return;
    }
    var arrow = editor.molecule.reactionArrow;
    var arrowLabel = reactionLabelForArrow(arrow);
    var reactionTitle = sanitizePathwayText(
        (editor.getMolName && editor.getMolName()) || arrowLabel || 'Reaction step', 90);
    reactionPayload.label = genericPathwayStructureLabel(reactionPayload.label) ? reactionTitle : reactionPayload.label;
    reactionPayload.structure.title = reactionTitle;
    reactionPayload.structure.comment = reactionPayload.structure.comment || 'Full drawn reaction payload';

    openWorkbenchSection('wb-pathway-section');
    if (typeof setWorkbenchMode === 'function') { setWorkbenchMode('pathway'); }
    if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }

    var center = pathwayVisibleCenter();
    var reactionNode = addPathwayNode(center.x, center.y, {
        label: reactionPayload.label,
        subtitle: reactionPayload.subtitle,
        kind: 'reaction',
        imageHref: reactionPayload.imageHref,
        structure: reactionPayload.structure,
        silent: true,
        batch: true
    });
    var step = addPathwayStep(center.x, center.y - 155, reactionTitle, { silent: true, batch: true });
    if (step) {
        step.comment = sanitizePathwayText(arrowLabel || 'Drawn reaction step', 90);
    }

    var yGap = 122;
    var leftX = center.x - 360;
    var rightX = center.x + 360;
    var baseYReact = center.y - ((reactants.length - 1) * yGap) / 2;
    var baseYProd = center.y - ((products.length - 1) * yGap) / 2;
    var i, component, node;
    for (i = 0; i < reactants.length; i++) {
        component = capturePathwayComponentPayload(reactants[i], 'Reactant', i + 1);
        node = addPathwayNode(leftX, baseYReact + i * yGap, {
            label: component.label,
            subtitle: component.subtitle,
            kind: component.kind,
            imageHref: component.imageHref,
            structure: component.structure,
            silent: true,
            batch: true
        });
        if (node && reactionNode) {
            addPathwayEdgeRecord(node.id, reactionNode.id, '', 'forward');
        }
    }
    var productArrow = reactionPathwayArrowType(arrow);
    for (i = 0; i < products.length; i++) {
        component = capturePathwayComponentPayload(products[i], 'Product', i + 1);
        node = addPathwayNode(rightX, baseYProd + i * yGap, {
            label: component.label,
            subtitle: component.subtitle,
            kind: component.kind,
            imageHref: component.imageHref,
            structure: component.structure,
            silent: true,
            batch: true
        });
        if (node && reactionNode) {
            addPathwayEdgeRecord(reactionNode.id, node.id, '', productArrow);
        }
    }
    selectPathwayItem('node', reactionNode.id);
    renderPathwayCanvas();
    fitPathwayCanvas();
    pathwayStatus('Converted reaction to pathway step: ' + reactants.length +
        ' reactant' + (reactants.length === 1 ? '' : 's') + ' -> ' +
        products.length + ' product' + (products.length === 1 ? '' : 's') + '.');
}

// v2.4.2 (Hybrid Stage C): make an IMPORTED reaction step atom-atom-mapped and
// link it so cross-pathway atom-tracing chains automatically across user-built
// pathways. Called AFTER convertCurrentReactionToPathway has placed the
// reactant -> reaction-step -> product nodes (and wired the edges) — it does
// NOT re-implement any of that. Given the reaction-step node id (defaults to
// the freshly-selected node), it:
//   (a) collects the boundary reaction's reactant SMILES (nodes whose edge
//       points INTO the step) and product SMILES (nodes the step points OUT to)
//       from each node's structure.smiles;
//   (b) runs the SAME id->positional-index path AtomTrace itself uses
//       (AtomTrace.mapStepIndexMap — the shared helper mapStep delegates to), so
//       the stored mapping is in the EXACT fromToIndex positional space the
//       trace chainer consumes — NEVER RDT internal atom-ids;
//   (c) attaches { mapping, score, status, version } onto the connecting
//       edge(s) via the edge record's aam field.
// Graceful no-op (like every other optional engine path) when RDT / AtomTrace /
// SMILES / a usable mapping is unavailable, so a build without the mapping
// engine still imports reactions exactly as in v2.4.1.
function linkImportedStepAAM(reactionNodeId) {
    if (typeof AtomTrace === 'undefined' || !AtomTrace.mapStepIndexMap ||
            typeof RDT === 'undefined' || !RDT.mapReaction) {
        return null;
    }
    if (!_pathway || !Array.isArray(_pathway.edges) || !Array.isArray(_pathway.nodes)) {
        return null;
    }
    var stepId = reactionNodeId || (_pathway.selectedType === 'node' ? _pathway.selectedId : null);
    if (!stepId) { return null; }
    var stepNode = findPathwayItem('node', stepId);
    if (!stepNode || stepNode.kind !== 'reaction') { return null; }

    // Partition the step's edges into reactant-side (into the step) and
    // product-side (out of the step). Keep the edge records so we can stamp the
    // aam back onto the actual connecting edges.
    var inEdges = [], outEdges = [];
    for (var e = 0; e < _pathway.edges.length; e++) {
        var ed = _pathway.edges[e];
        if (ed.to === stepId) { inEdges.push(ed); }
        else if (ed.from === stepId) { outEdges.push(ed); }
    }
    if (!inEdges.length || !outEdges.length) { return null; }

    var reactantSmiles = pathwayStepSideSmiles(inEdges, 'from');
    var productSmiles = pathwayStepSideSmiles(outEdges, 'to');
    if (!reactantSmiles || !productSmiles) { return null; }

    var result;
    try {
        result = AtomTrace.mapStepIndexMap(reactantSmiles, productSmiles);
    } catch (err) {
        return null;
    }
    if (!result || !result.fromToIndex) { return null; }
    var mapped = 0;
    for (var k in result.fromToIndex) {
        if (result.fromToIndex.hasOwnProperty(k)) { mapped++; }
    }
    if (mapped === 0) { return null; }

    // The boundary mapping describes the WHOLE reactant-side >> product-side
    // reaction in fromToIndex positional space. Stamp the SAME payload onto
    // every connecting edge so the live-pathway trace (which collapses
    // reactant -> step -> product into a direct reactant -> product trace edge)
    // can pick it up from either side.
    var aam = {
        mapping: result.fromToIndex,
        score: (typeof result.score === 'number') ? result.score : null,
        status: result.status || '',
        version: AAM_LINK_VERSION
    };
    var i;
    for (i = 0; i < inEdges.length; i++) { inEdges[i].aam = cloneAamPayload(aam); }
    for (i = 0; i < outEdges.length; i++) { outEdges[i].aam = cloneAamPayload(aam); }
    return aam;
}

// v2.4.2: collect the `.`-joined SMILES of the structure nodes on one side of a
// reaction-step node. `endKey` is 'from' for reactant-side edges (the structure
// node is the edge's `from`) or 'to' for product-side edges. Returns '' when any
// participating structure node lacks a usable SMILES (caller treats '' as a
// no-op signal — we never fabricate a partial reaction).
function pathwayStepSideSmiles(edges, endKey) {
    var parts = [];
    for (var i = 0; i < edges.length; i++) {
        var nodeId = edges[i][endKey];
        var node = findPathwayItem('node', nodeId);
        if (!node || node.kind === 'reaction') { continue; }
        var smi = node.structure && node.structure.smiles ? String(node.structure.smiles).trim() : '';
        if (!smi) { return ''; }
        parts.push(smi);
    }
    if (!parts.length) { return ''; }
    return parts.join('.');
}

// v2.4.2: defensive deep-ish copy of an aam payload so each edge owns its own
// record (the mapping object is shared-by-value otherwise, which would let a
// later in-place edit on one edge silently mutate the sibling edge's mapping).
function cloneAamPayload(aam) {
    var mapping = {};
    for (var k in aam.mapping) {
        if (aam.mapping.hasOwnProperty(k)) { mapping[k] = aam.mapping[k]; }
    }
    return { mapping: mapping, score: aam.score, status: aam.status, version: aam.version };
}

function applyPathwayStructurePayload(node, payload) {
    if (!node || !payload) return;
    var nextLabel = sanitizePathwayText(payload.label || '', 90);
    if (nextLabel && (!node.label || genericPathwayStructureLabel(node.label) || !genericPathwayStructureLabel(nextLabel))) {
        node.label = nextLabel;
    }
    node.subtitle = sanitizePathwayText(payload.subtitle || '', 42);
    node.kind = payload.kind === 'reaction' ? 'reaction' : pathwayNodeKindForLabel(node.label, node.kind || 'metabolite');
    node.imageHref = payload.imageHref || '';
    node.structure = normalizePathwayStructurePayload(payload.structure);
    var wantsReactionSize = node.kind === 'reaction' || (node.structure && node.structure.type === 'reaction');
    if (node.imageHref) {
        var nextW = wantsReactionSize ? 240 : 180;
        var nextH = wantsReactionSize ? 124 : 112;
        if (node.w !== nextW || node.h !== nextH) {
            var cx = node.x + node.w / 2;
            var cy = node.y + node.h / 2;
            node.w = nextW;
            node.h = nextH;
            node.x = Math.round(cx - node.w / 2);
            node.y = Math.round(cy - node.h / 2);
        }
    }
}

function addCurrentMoleculeToPathway() {
    if (!editor || !editor.molecule || !editor.molecule.atoms || editor.molecule.atoms.length === 0) {
        // Merge Stage 3 (v2.3.1): the continuum has no separate Molecule view to
        // switch back to. When "Draw Structure" is clicked with an empty editor,
        // drop a FRESH structure node on the canvas and open the focus lens on it
        // (newPathwayNodeInLens) — the user draws directly in the lens and reads
        // its SMILES/MOL/RXN while focused. This preserves the single-molecule
        // workflow: draw one structure -> read its output, all on one surface.
        // Falls back to the legacy setWorkbenchMode('molecule') shim only when
        // the lens is unavailable (older bundle without FocusLens).
        openWorkbenchSection('wb-pathway-section');
        if (_lens && newPathwayNodeInLens()) { return; }
        var emptyHint = 'Editor is empty — use "Draw Structure" again to start a new structure once the lens is available.';
        if (typeof setWorkbenchMode === 'function') {
            setWorkbenchMode('molecule');
        }
        pathwayStatus(emptyHint);
        return;
    }
    openWorkbenchSection('wb-pathway-section');
    var payload = captureEditorStructureForPathway();
    if (!payload) {
        pathwayStatus('No molecule or reaction is ready to add to the pathway canvas.');
        return;
    }
    // v2.0.60: edit-in-place. If a pathway node is currently selected,
    // OVERWRITE its structure (imageHref + subtitle) instead of stamping
    // a new node. This is the natural fix for nodes that landed with the
    // wrong / missing structure — select the node, draw the correction
    // in the editor, click "Draw Structure" and the node updates in
    // place. To force-create a fresh node when a node is selected,
    // deselect first (Select tool + click empty canvas).
    if (_pathway.selectedType === 'node' && _pathway.selectedId) {
        var node = findPathwayItem('node', _pathway.selectedId);
        if (node) {
            applyPathwayStructurePayload(node, payload);
            renderPathwayCanvas();
            pathwayStatus('Updated ' + (payload.kind === 'reaction' ? 'reaction' : 'structure') + ' on "' + node.label + '" (' + node.id + ').');
            return;
        }
    }
    var center = pathwayVisibleCenter();
    addPathwayNode(center.x, center.y, {
        label: payload.label,
        subtitle: payload.subtitle,
        kind: payload.kind,
        imageHref: payload.imageHref,
        structure: payload.structure
    });
}

// Hybrid Stage B (v2.4.1): the "Add to pathway" entry point from the EDITOR
// view. It is a thin ENTRY POINT + auto-tag glue layered on the existing
// import spine — it never re-implements the capture/split/stamp logic:
//   - a drawn REACTION (editor.molecule.reactionArrow truthy) ->
//     convertCurrentReactionToPathway(), which splits the reaction via
//     RDT._splitReactionSides and drops reactant nodes + one kind:'reaction'
//     step node + product nodes (wired with addPathwayEdgeRecord);
//   - a plain non-empty MOLECULE -> addCurrentMoleculeToPathway(), which stamps
//     one node. We CLEAR the pathway selection first so it always stamps a FRESH
//     node and can never overwrite a stale selection (addCurrentMoleculeToPathway
//     keeps its own edit-in-place branch for the Pathway-view "Draw Structure"
//     workflow; here we deliberately bypass it).
// The empty-editor case is guarded with a friendly status and drops no node.
// After a successful import we switch to the Pathway view + fit the canvas so the
// user sees the result, and surface a unified `Added "<label>" to the pathway.`
// status. The <label> reuses the existing namer (captureEditorStructureForPathway
// -> pathwayInputLabel / moleculeFormulaLabel / pathwayNodeKindForLabel) — no new
// label logic here.
//
// v2.4.2 (Hybrid Stage C): after the REACTION branch imports the
// reactant/step/product subgraph, linkImportedStepAAM atom-atom-maps the
// boundary step and stores the mapping on the connecting edges, so cross-pathway
// atom tracing chains automatically across user-built pathways. Only the
// reaction path needs this — a molecule lands as a single node with no step to
// map. The call is a graceful no-op if the mapping engine / SMILES are
// unavailable, so the import itself is unaffected.
function addEditorDrawingToPathway() {
    if (!editor || !editor.molecule || !editor.molecule.atoms || editor.molecule.atoms.length === 0) {
        pathwayStatus('Draw a molecule or reaction in the editor first, then "Add to pathway".');
        return;
    }
    // Capture the auto-tag label up front (the editor contents do not change as
    // we import, so this is the label the new node/step will carry). Reuses the
    // shared namer; falls back to a sensible default if capture is unavailable.
    var payload = captureEditorStructureForPathway();
    var isReaction = !!editor.molecule.reactionArrow;
    var label = (payload && payload.label) || (isReaction ? 'Reaction' : 'Structure');
    if (isReaction) {
        // Reaction -> reactant/step/product subgraph. convertCurrentReactionToPathway
        // owns the split + node/edge stamping (and its own guards for one-sided
        // reactions / an unavailable splitter), so we delegate wholesale.
        convertCurrentReactionToPathway();
        // v2.4.2: AAM-link the freshly imported step (convert selects the new
        // reaction-step node, so linkImportedStepAAM defaults to it). No-op when
        // the engine / SMILES are unavailable.
        linkImportedStepAAM(_pathway && _pathway.selectedType === 'node' ? _pathway.selectedId : null);
    } else {
        // Plain molecule -> a single FRESH node. Clear any pathway selection first
        // so addCurrentMoleculeToPathway's edit-in-place branch never fires here.
        selectPathwayItem('', null);
        addCurrentMoleculeToPathway();
    }
    // Show the result: land the user on the Pathway view and fit the canvas.
    if (typeof setWorkbenchMode === 'function') { setWorkbenchMode('pathway'); }
    if (typeof fitPathwayCanvas === 'function') { fitPathwayCanvas(); }
    // v2.4.2: refresh the atom-trace strip for the live pathway. It surfaces
    // once the user has built >=2 structure nodes (e.g. two imported reaction
    // steps that share a metabolite), so a single first import stays hidden but
    // the second un-hides it. No-op for an example pathway (that lane is driven
    // by loadExamplePathway).
    if (typeof showAtomTraceStrip === 'function' && !(_pathway && _pathway.loadedExampleId)) {
        showAtomTraceStrip();
    }
    pathwayStatus('Added "' + label + '" to the pathway.');
}

// v2.0.60: load the SELECTED pathway node's structure into the molecule
// editor for revision.
//
// Merge Stage 3 (v2.3.1): in the continuum this routes to the Stage-2 focus
// lens — the editor floats IN PLACE over the selected node (openStructureLens),
// Esc / click-away commits it back. The legacy body below (load(source) into
// the main editor + setWorkbenchMode('molecule')) is retained as the fallback
// for bundles without FocusLens, and keeps the source-shape-pinned substrings
// (pathwayNodeStructureInput / load(source) / setWorkbenchMode('molecule')).
function editPathwayNodeStructure() {
    if (_pathway.selectedType !== 'node' || !_pathway.selectedId) {
        pathwayStatus('Select a pathway node first, then double-click it (or click "Edit Structure") to edit its structure in the focus lens.');
        return;
    }
    var node = findPathwayItem('node', _pathway.selectedId);
    if (!node) {
        pathwayStatus('Pathway node not found — selection may be stale.');
        return;
    }
    // Continuum: prefer the in-place focus lens (overlay editor over the node).
    if (_lens && openStructureLens(_pathway.selectedId)) { return; }
    var source = pathwayNodeStructureInput(node);
    var loaded = false;
    if (source && editor && editor.readGenericMolecularInput) {
        try {
            load(source);
            loaded = !!(editor.molecule && editor.molecule.atoms && editor.molecule.atoms.length);
            var structure = pathwayNodeStructurePayload(node);
            if (loaded && structure) {
                if (structure.title && editor.setMolName) { editor.setMolName(structure.title); }
                if (editor.molecule) { editor.molecule.comment = structure.comment || ''; }
                if (editor.render) { editor.render(); }
                queueOutputUpdate();
                scheduleWorkbenchUpdate();
            }
        } catch (e) { loaded = false; }
    }
    if (typeof setWorkbenchMode === 'function') {
        setWorkbenchMode('molecule');
    }
    if (loaded) {
        pathwayStatus('Loaded "' + node.label + '" into the editor. Edit it, then click "Draw Structure" to save back to the node.');
    } else {
        pathwayStatus('Could not parse the saved structure for "' + node.label + '". Draw a fresh molecule or reaction and click "Draw Structure" to save it onto the node.');
    }
}

// v2.0.19: render a metabolite SMILES to a compact SVG thumbnail (data URI) for
// use as a pathway-node image. Parses, lays out, and exports — independent of
// the molecule editor's current contents.
function metaboliteImageHref(smiles) {
    if (!smiles || typeof SmilesParser === 'undefined' || typeof Molecule === 'undefined') return '';
    if (typeof ImageExport === 'undefined' || !ImageExport.toSVG) return '';
    try {
        var m = new Molecule();
        SmilesParser.parse(smiles, m);
        if (!m.atoms.length) return '';
        if (typeof Layout !== 'undefined' && Layout.layout) {
            // Seed atoms on a small deterministic grid so the layout starts from
            // a non-degenerate (not all-coincident) state, then lay out.
            for (var i = 0; i < m.atoms.length; i++) {
                m.atoms[i].x = (i % 6) * 0.5 - 1.25;
                m.atoms[i].y = Math.floor(i / 6) * 0.5 - 1.25;
            }
            Layout.layout(m);
        }
        var svg = ImageExport.toSVG(m, {
            width: 172, height: 96, padding: 8, background: 'transparent', showAtomLabels: true
        });
        return pathwaySvgDataUri(svg);
    } catch (e) {
        return '';
    }
}

function normalizePathwayLayoutShape(shape) {
    if (typeof PathwayLayout !== 'undefined' && PathwayLayout.normalizeLayoutShape) {
        return PathwayLayout.normalizeLayoutShape(shape);
    }
    var s = String(shape || 'auto').toLowerCase().trim();
    s = s.replace(/\s+/g, '-').replace(/_/g, '-');
    if (s === 'linear' || s === 'line') { return 'ladder'; }
    if (s === 'branch') { return 'branched'; }
    if (s === 'fan-out' || s === 'hub' || s === 'hub-and-spoke') { return 'fanout'; }
    if (s === 'fan-in' || s === 'convergent' || s === 'converge' || s === 'merge' || s === 'sink') { return 'fanin'; }
    if (s === 'splitmerge' || s === 'branch-merge' || s === 'branchmerge' ||
            s === 'fork-join' || s === 'forkjoin' || s === 'diamond') { return 'split-merge'; }
    if (s === 'mixed' || s === 'combo' || s === 'combined' ||
            s === 'combination' || s === 'combinations') { return 'hybrid'; }
    if (s === 'loop' || s === 'iterative') { return 'loop-iterative'; }
    if (s === 'auto' || s === 'hierarchical' || s === 'ladder' || s === 'cycle' ||
            s === 'fanout' || s === 'fanin' || s === 'branched' ||
            s === 'split-merge' || s === 'hybrid' || s === 'loop-iterative') {
        return s;
    }
    return 'auto';
}

// v2.0.19: build a canonical pathway from the bundled MetaboliteLibrary,
// drawing each small metabolite as a real structure (big cofactors as labelled
// nodes) and each reaction with its enzyme/EC label + arrow type. Template
// shapes give common pathway maps a useful starting geometry before the
// general Clean Up pass is applied.
function loadExamplePathway(id) {
    if (typeof MetaboliteLibrary === 'undefined' || !MetaboliteLibrary.pathways) {
        pathwayStatus('Metabolite library is unavailable in this build.');
        return;
    }
    var pw = MetaboliteLibrary.pathways[id];
    if (!pw) { pathwayStatus('Unknown example pathway: ' + id); return; }
    setWorkbenchMode('pathway');
    clearPathwayCanvas();

    // unique metabolites in first-appearance order
    var order = [], seen = {};
    pw.steps.forEach(function(s) {
        [s.from, s.to].forEach(function(n) { if (!seen[n]) { seen[n] = true; order.push(n); } });
    });

    var nodeByName = {};
    var shape = pw.shape || 'ladder';
    _pathway.layoutShape = normalizePathwayLayoutShape(shape);

    var layoutNodes = [];
    var layoutEdges = [];
    for (var ln = 0; ln < order.length; ln++) { layoutNodes.push({ id: order[ln] }); }
    for (var le = 0; le < pw.steps.length; le++) {
        layoutEdges.push({ from: pw.steps[le].from, to: pw.steps[le].to });
    }
    var positionsByName = {};
    if (typeof PathwayLayout !== 'undefined' && PathwayLayout.compute) {
        var laid = PathwayLayout.compute(layoutNodes, layoutEdges, {
            layoutShape: _pathway.layoutShape,
            order: order,
            hub: pw.hub || null,
            trunk: pw.trunk || null
        });
        positionsByName = (laid && laid.positions) || {};
    }
    var N = order.length;
    for (var i = 0; i < N; i++) {
        var name = order[i];
        var meta = MetaboliteLibrary.find(name) || { name: name, display: 'label', smiles: '', kegg: '', aliases: [] };
        var label = (meta.aliases && meta.aliases.length) ? meta.aliases[0] : meta.name;
        var img = (meta.display === 'structure') ? metaboliteImageHref(meta.smiles) : '';
        var structure = meta.smiles ? normalizePathwayStructurePayload({
            type: 'molecule',
            smiles: meta.smiles,
            mol: '',
            rxn: '',
            title: meta.name || label,
            comment: meta.kegg ? ('Template metabolite ' + meta.kegg) : ''
        }) : null;
        var p = positionsByName[name] || { x: 300, y: 120 + i * 130 };
        var node = addPathwayNode(p.x, p.y, {
            label: label,
            subtitle: meta.name + (meta.kegg ? ' · ' + meta.kegg : ''),
            kind: (meta.display === 'structure') ? 'metabolite' : 'cofactor',
            imageHref: img,
            structure: structure
        });
        nodeByName[name] = node.id;
    }

    // one edge per reaction, enzyme + EC + cofactors as the label, arrow type set
    for (var j = 0; j < pw.steps.length; j++) {
        var s = pw.steps[j];
        var fromId = nodeByName[s.from], toId = nodeByName[s.to];
        if (!fromId || !toId) continue;
        var lbl = s.enzyme + (s.ec ? ' (' + s.ec + ')' : '');
        if (s.cofactors) { lbl += '  ·  ' + s.cofactors; }
        _pathway.edges.push({
            id: newPathwayId('edge-'),
            from: fromId,
            to: toId,
            label: sanitizePathwayText(lbl, 64),
            arrowType: (PATHWAY_ARROW_TYPES.indexOf(s.arrowType) >= 0) ? s.arrowType : 'forward'
        });
    }

    fitPathwayCanvas();
    renderPathwayCanvas();
    pathwayStatus('Loaded ' + pw.name + ' — ' + order.length + ' metabolites, ' + pw.steps.length + ' reactions.');
    // v2.0.22: surface the atom-trace strip for the loaded example.
    _pathway.loadedExampleId = id;
    showAtomTraceStrip(null, null);
}

// v2.0.22: build the metabolite list for an example pathway in first-appearance
// order (same order loadExamplePathway uses), each entry carrying the parsed +
// laid-out Molecule and its repertoire metadata.
function buildExamplePathwayMetabolites(id) {
    if (typeof MetaboliteLibrary === 'undefined' || !MetaboliteLibrary.pathways[id]) return null;
    var pw = MetaboliteLibrary.pathways[id];
    var order = [], seen = {};
    pw.steps.forEach(function(s) { [s.from, s.to].forEach(function(n) { if (!seen[n]) { seen[n] = true; order.push(n); } }); });
    var out = [];
    for (var i = 0; i < order.length; i++) {
        var meta = MetaboliteLibrary.find(order[i]) || { name: order[i], smiles: '', aliases: [], kegg: '' };
        var mol = new Molecule();
        if (meta.smiles && typeof SmilesParser !== 'undefined') {
            try { SmilesParser.parse(meta.smiles, mol); } catch (e) { /* skip */ }
        }
        if (mol.atoms.length && typeof Layout !== 'undefined' && Layout.layout) {
            for (var ai = 0; ai < mol.atoms.length; ai++) {
                mol.atoms[ai].x = (ai % 6) * 0.5 - 1.25;
                mol.atoms[ai].y = Math.floor(ai / 6) * 0.5 - 1.25;
            }
            Layout.layout(mol);
        }
        out.push({ name: order[i], meta: meta, mol: mol });
    }
    return out;
}

// v2.3.7 (P1-5c): short non-colour status tag for a moiety-trace verdict. Pairs
// the distinct CSS border PATTERN (solid/dashed/dotted/double) with a glyph +
// word so the verdict is readable in greyscale and announced by a screen
// reader, not carried by hue alone. Mirrors AtomTrace.deriveMoietyTrace's
// status vocabulary (intact / fragmented / partial-loss / absent).
var ATOM_TRACE_MOIETY_TAGS = {
    'intact': '✓ intact',
    'fragmented': '⚠ fragmented',
    'partial-loss': '◐ partial',
    'absent': '✕ absent'
};

function atomTraceMoietyTag(status) {
    return Object.prototype.hasOwnProperty.call(ATOM_TRACE_MOIETY_TAGS, status) ?
        ATOM_TRACE_MOIETY_TAGS[status] : '';
}

// v2.0.22: render the interactive atom-trace strip. Each metabolite becomes a
// live inline SVG (AtomTrace.renderInteractiveSvg) whose atoms are clickable
// and CSS-highlightable. If startMetaboliteIdx/startAtomIdx are provided, run
// AtomTrace.tracePathway over metabolites[start..] and highlight the traced
// atom in each cell + show a status line.
function showAtomTraceStrip(startMetaboliteIdx, startAtomIdx) {
    var section = document.getElementById('wb-atom-trace-section');
    var grid = document.getElementById('wb-atom-trace-grid');
    var status = document.getElementById('wb-atom-trace-status');
    if (!section || !grid) return;
    if (typeof AtomTrace === 'undefined' || !AtomTrace.renderInteractiveSvg) {
        section.hidden = true;
        return;
    }
    var id = _pathway && _pathway.loadedExampleId;
    // v2.0.22 / v2.4.2: the metabolite array driving the strip. For a loaded
    // EXAMPLE pathway this is the cached fast lane (parse + Layout once, reuse
    // across click-retraces). For a USER-built pathway (loadedExampleId falsy)
    // we derive it on the fly from the live structure nodes — this is the
    // Stage-C generalisation that brings the trace strip to user pathways.
    var mets;
    if (id) {
        mets = _pathway.loadedExamplePathwayMetabolites;
        if (!mets || !mets.length) {
            mets = buildExamplePathwayMetabolites(id);
            if (!mets || !mets.length) { section.hidden = true; return; }
            _pathway.loadedExamplePathwayMetabolites = mets;
        }
    } else {
        // Live user pathway: need at least two structure nodes (with SMILES) to
        // have anything to chain a trace through. Otherwise stay hidden.
        mets = buildLivePathwayMetabolites();
        if (!mets || mets.length < 2) { section.hidden = true; return; }
    }

    // optional trace over the sub-pathway starting at startMetaboliteIdx
    var labels = null, traceInfo = null;
    // v2.0.54: optional moiety-trace overlay state. Populated only when
    // the user has built a connected moiety set on this start metabolite.
    var moietyTraceInfo = null;
    var moietyConnectError = false;
    // v2.0.59: cofactor lineage labels lifted into the full-strip frame.
    // cofactorLabelsByCell[j] is null OR array of per-atom cofactor-origin
    // records [{cofactor, atomIdx, step}, …]. Populated by v2.0.56's
    // tracePathway when the example pathway carries cofactor strings.
    var cofactorLabelsByCell = null;
    var hasTrace = (startMetaboliteIdx !== null && startMetaboliteIdx !== undefined &&
                    startAtomIdx !== null && startAtomIdx !== undefined);
    if (hasTrace && startMetaboliteIdx < mets.length - 1) {
        // v2.4.2: the slice / DAG-edges / cofactor-spec / cached-mapping
        // derivation now lives in buildPathwayTraceInputs, which handles BOTH
        // the example fast lane (DAG + cofactors from the MetaboliteLibrary
        // step list, as before — calls AtomTrace.parseCofactorString per step)
        // AND a live user pathway (edges from _pathway.edges re-indexed to the
        // slice, preferring a stored edge.aam mapping). showAtomTraceStrip just
        // assembles traceOpts and runs the trace + render.
        var inputs = buildPathwayTraceInputs(mets, startMetaboliteIdx);
        var slice = inputs.slice;
        var traceEdges = inputs.traceEdges;
        var traceCofactors = inputs.traceCofactors;
        var traceOpts = {};
        if (traceEdges.length > 0) {
            traceOpts.edges = traceEdges;
            // v2.0.59: hand the cofactor spec to AtomTrace so it builds
            // balanced rxns (from + cofIn >> to + cofOut) at each edge.
            // cofactorLabels[m][k] on the result records atoms whose
            // ancestry includes a cofactor source.
            if (traceCofactors && traceCofactors.length === traceEdges.length) {
                traceOpts.cofactors = traceCofactors;
            }
            // v2.4.2: when the live pathway stored a per-edge mapping (edge.aam)
            // prefer it — AtomTrace.tracePathway uses it verbatim (same
            // positional space) instead of recomputing via mapReaction.
            if (inputs.edgeMaps && inputs.edgeMaps.length === traceEdges.length) {
                traceOpts.edgeMaps = inputs.edgeMaps;
            }
        }
        var result = (slice.length >= 2 && AtomTrace.tracePathway) ?
            AtomTrace.tracePathway(slice, traceOpts) : null;
        if (result && result.traces && result.traces[startAtomIdx]) {
            traceInfo = result.traces[startAtomIdx];
            labels = [];
            for (var p = 0; p < mets.length; p++) { labels[p] = []; }
            for (var li = 0; li < result.labels.length; li++) {
                labels[startMetaboliteIdx + li] = result.labels[li];
            }
            // v2.0.59: lift cofactorLabels into the full-strip frame so the
            // render loop can stamp a per-atom is-cofactor-origin halo.
            if (Array.isArray(result.cofactorLabels)) {
                cofactorLabelsByCell = [];
                for (var cp = 0; cp < mets.length; cp++) { cofactorLabelsByCell[cp] = null; }
                for (var cli = 0; cli < result.cofactorLabels.length; cli++) {
                    cofactorLabelsByCell[startMetaboliteIdx + cli] = result.cofactorLabels[cli];
                }
            }
        }
        // v2.0.54: if the user has built a moiety on this start metabolite,
        // run the post-processor and surface per-cell aliveness + bond
        // preservation. The result is keyed against the slice (0..n-1) and
        // we lift it back into the full-strip frame by adding startMetaboliteIdx.
        if (result && Array.isArray(_pathway.moietySet) &&
                _pathway.moietySet.length >= 2 &&
                _pathway.moietyStartMetIdx === startMetaboliteIdx &&
                typeof AtomTrace.deriveMoietyTrace === 'function') {
            var moietyResult = AtomTrace.deriveMoietyTrace(
                result, _pathway.moietySet, { metaboliteSmiles: slice });
            if (moietyResult && !moietyResult.error && moietyResult.perMetabolite) {
                moietyTraceInfo = {
                    base: moietyResult,
                    startOffset: startMetaboliteIdx,
                    moietyAtoms: _pathway.moietySet.slice()
                };
            } else if (moietyResult && moietyResult.error === 'not-connected') {
                moietyConnectError = true;
            }
        }
    }

    grid.textContent = '';
    // v2.0.54: derive per-cell moiety status before the render loop so we
    // can stamp the right border class on every cell. moietyByCell[j] is
    // either null (no moiety in flight) or { status, alive, present, ... }.
    var moietyByCell = null;
    if (moietyTraceInfo && moietyTraceInfo.base) {
        moietyByCell = [];
        var base = moietyTraceInfo.base;
        for (var mc = 0; mc < mets.length; mc++) { moietyByCell.push(null); }
        for (var pm = 0; pm < base.perMetabolite.length; pm++) {
            moietyByCell[moietyTraceInfo.startOffset + pm] = base.perMetabolite[pm];
        }
    }
    for (var j = 0; j < mets.length; j++) {
        var cell = document.createElement('div');
        cell.className = 'wb-atom-trace-cell';
        cell.setAttribute('data-metabolite-index', String(j));
        // v2.0.54: per-cell moiety status border.
        // v2.3.7 (P1-5c): the border PATTERN (CSS) now distinguishes the four
        // verdicts without colour; pair it with a short status glyph + word in
        // the title so the verdict also reads in greyscale and to a screen
        // reader (not encoded by hue alone).
        var moietyTag = '';
        if (moietyByCell && moietyByCell[j]) {
            var st = moietyByCell[j].status;
            if (st && st !== 'unreachable') {
                cell.classList.add('is-moiety-' + st);
                moietyTag = atomTraceMoietyTag(st);
            }
        }
        var title = document.createElement('div');
        title.className = 'wb-atom-trace-cell-title';
        var alias = (mets[j].meta && mets[j].meta.aliases && mets[j].meta.aliases.length) ?
            mets[j].meta.aliases[0] : mets[j].name;
        if (moietyTag) {
            var tagEl = document.createElement('span');
            tagEl.className = 'wb-atom-trace-moiety-tag';
            tagEl.textContent = moietyTag;
            tagEl.setAttribute('title', 'Moiety ' + moietyByCell[j].status);
            title.appendChild(tagEl);
            title.appendChild(document.createTextNode(alias));
        } else {
            title.textContent = alias;
        }
        cell.appendChild(title);
        var holder = document.createElement('div');
        holder.className = 'wb-atom-trace-cell-svg';
        // BIME-generated SVG only (no user input goes through here) — innerHTML
        // is the simplest way to embed inline SVG so atoms remain addressable.
        holder.innerHTML = AtomTrace.renderInteractiveSvg(mets[j].mol, { width: 168, height: 100 });
        // v2.0.54: highlight every moiety image atom with a secondary halo
        // (independent of the primary single-atom highlight below). For the
        // start metabolite, the source moiety atoms get the same halo so the
        // user sees the selection.
        if (moietyByCell && moietyByCell[j]) {
            var imgList = moietyByCell[j].imageAtoms || [];
            for (var ia = 0; ia < imgList.length; ia++) {
                var imgG = holder.querySelector('.bime-trace-atom[data-atom-index="' + imgList[ia] + '"]');
                if (imgG) { imgG.classList.add('is-moiety-source'); }
            }
        }
        // v2.0.59: cofactor-origin halo. Every atom whose ancestry includes
        // a cofactor IN (parsed from the pathway-step `cofactors` string)
        // gets a distinct class + an SVG <title> tooltip listing the
        // cofactor sources. Renders alongside (not in place of) the primary
        // single-atom and moiety highlights.
        if (cofactorLabelsByCell && cofactorLabelsByCell[j]) {
            var cofBag = cofactorLabelsByCell[j];
            for (var ck = 0; ck < cofBag.length; ck++) {
                var cofRecs = cofBag[ck];
                if (!cofRecs || !cofRecs.length) { continue; }
                var cofG = holder.querySelector('.bime-trace-atom[data-atom-index="' + ck + '"]');
                if (!cofG) { continue; }
                cofG.classList.add('is-cofactor-origin');
                // De-dup cofactor names for the tooltip — atom may have
                // multiple records pointing at the same cofactor.
                var nameSet = {};
                var names = [];
                for (var cri = 0; cri < cofRecs.length; cri++) {
                    var cname = cofRecs[cri].cofactor;
                    if (!nameSet[cname]) { nameSet[cname] = true; names.push(cname); }
                }
                var titleEl = document.createElementNS('http://www.w3.org/2000/svg', 'title');
                titleEl.textContent = 'Came from cofactor: ' + names.join(', ');
                cofG.appendChild(titleEl);
            }
        }
        if (labels && labels[j] && labels[j].length) {
            for (var ki = 0; ki < labels[j].length; ki++) {
                // v2.0.53: labels[j][ki] is now a multi-valued array (a single
                // downstream atom can have >1 ancestor under convergence —
                // e.g. F1,6BP → DHAP → GAP + F1,6BP → GAP direct both
                // converge on GAP). Highlight if the active startAtomIdx
                // is anywhere in the bag.
                var bag = labels[j][ki];
                var matched = (typeof bag === 'number')
                    ? (bag === startAtomIdx)
                    : (Array.isArray(bag) && bag.indexOf(startAtomIdx) !== -1);
                if (matched) {
                    var atomG = holder.querySelector('.bime-trace-atom[data-atom-index="' + ki + '"]');
                    if (atomG) {
                        atomG.classList.add('is-traced');
                        // v2.0.24: chemist-style lineage label — show the
                        // start atom's 1-indexed number on every downstream
                        // copy (matches the [C1]/[C2] convention used in
                        // textbook carbon-fate diagrams).
                        var hit = atomG.querySelector('.bime-trace-atom-hit');
                        if (hit) {
                            var num = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                            num.setAttribute('class', 'bime-trace-number');
                            num.setAttribute('pointer-events', 'none');
                            num.setAttribute('x', String(parseFloat(hit.getAttribute('cx')) + 9));
                            num.setAttribute('y', String(parseFloat(hit.getAttribute('cy')) - 6));
                            num.textContent = String(startAtomIdx + 1);
                            atomG.appendChild(num);
                        }
                    }
                }
            }
        }
        cell.appendChild(holder);
        grid.appendChild(cell);
    }

    if (status) {
        // v2.0.54: moiety status wins over single-atom status when active.
        if (moietyConnectError) {
            status.textContent = 'Moiety must be connected — add the bridging atom (shift-click).';
        } else if (moietyTraceInfo && moietyTraceInfo.base) {
            var mb = moietyTraceInfo.base;
            var sizeTxt = mb.moietySize + '-atom moiety';
            var msg2;
            if (mb.breakAt === null) {
                var lastIdx = (mb.perMetabolite.length - 1) + moietyTraceInfo.startOffset;
                msg2 = sizeTxt + ' (' + mb.moietyAtoms.map(function(a) { return '#' + (a + 1); }).join(',') +
                    ') of ' + mets[moietyTraceInfo.startOffset].name +
                    ' → survives intact to ' + mets[lastIdx].name;
            } else {
                var breakName = mets[moietyTraceInfo.startOffset + mb.breakAt].name;
                var breakEntry = mb.perMetabolite[mb.breakAt];
                msg2 = sizeTxt + ' of ' + mets[moietyTraceInfo.startOffset].name +
                    ' → ' + breakEntry.status + ' at ' + breakName +
                    ' (preservedBonds ' + breakEntry.preservedBonds.length + '/' + mb.startBonds.length +
                    ', atoms ' + breakEntry.present.length + '/' + mb.moietySize + ')';
            }
            status.innerHTML = '';
            status.appendChild(document.createTextNode(msg2 + ' '));
            var clearLink = document.createElement('button');
            clearLink.type = 'button';
            clearLink.className = 'wb-link';
            clearLink.setAttribute('data-wb-action', 'clear-moiety-trace');
            clearLink.textContent = 'Clear moiety';
            status.appendChild(clearLink);
            status.appendChild(document.createTextNode(
                ' (shift-click to extend/remove atoms; click any atom for single-atom trace).'));
        } else if (traceInfo) {
            var endIdx = startMetaboliteIdx + (traceInfo.path.length - 1);
            var endName = mets[Math.min(endIdx, mets.length - 1)].name;
            var msg = 'Tracing ' + traceInfo.element + ' #' + (startAtomIdx + 1) + ' of ' +
                mets[startMetaboliteIdx].name + ' → ' +
                (traceInfo.lostAtStep === null ? 'reaches ' + endName :
                    'lost in the step into ' + mets[startMetaboliteIdx + traceInfo.lostAtStep].name);
            status.textContent = msg + ' (click any atom to retrace, shift-click to start a moiety).';
        } else {
            status.textContent = 'Click an atom in any structure to trace its fate forward through the pathway (shift-click to start a moiety).';
        }
    }
    section.hidden = false;
}

// v2.4.2 (Hybrid Stage C): build the trace metabolite array for a LIVE user
// pathway (no loadedExampleId). Each STRUCTURE node (kind !== 'reaction') that
// carries a usable structure.smiles becomes one strip cell, in _pathway.nodes
// order. Mirrors the shape buildExamplePathwayMetabolites returns
// ({name, meta:{smiles, aliases}, mol}) so the render loop + trace setup are
// lane-agnostic, plus a `nodeId` so buildPathwayTraceInputs can re-index the
// live edges onto the slice. Reaction-step nodes and node-less / SMILES-less
// nodes are skipped (they are not metabolites to trace through). Returns [] when
// the engine is unavailable.
function buildLivePathwayMetabolites() {
    if (typeof Molecule === 'undefined' || typeof SmilesParser === 'undefined') { return []; }
    if (!_pathway || !Array.isArray(_pathway.nodes)) { return []; }
    var out = [];
    for (var i = 0; i < _pathway.nodes.length; i++) {
        var node = _pathway.nodes[i];
        if (!node || node.kind === 'reaction') { continue; }
        var smi = node.structure && node.structure.smiles ? String(node.structure.smiles).trim() : '';
        if (!smi) { continue; }
        var mol = new Molecule();
        try { SmilesParser.parse(smi, mol); } catch (e) { continue; }
        if (!mol.atoms.length) { continue; }
        // Lay out atom coordinates the same way the example lane does so the
        // interactive SVG renders sensibly.
        if (typeof Layout !== 'undefined' && Layout.layout) {
            for (var ai = 0; ai < mol.atoms.length; ai++) {
                mol.atoms[ai].x = (ai % 6) * 0.5 - 1.25;
                mol.atoms[ai].y = Math.floor(ai / 6) * 0.5 - 1.25;
            }
            Layout.layout(mol);
        }
        var name = node.label || node.subtitle || ('Structure ' + (out.length + 1));
        out.push({
            name: name,
            nodeId: node.id,
            meta: { smiles: smi, aliases: [], kegg: '' },
            mol: mol
        });
    }
    return out;
}

// v2.4.2 (Hybrid Stage C): derive the AtomTrace.tracePathway inputs (slice +
// DAG edges + optional cofactor specs + optional cached per-edge maps) for the
// metabolite array `mets`, sub-sliced from startMetaboliteIdx. Two lanes:
//
//   EXAMPLE (loadedExampleId truthy) — slice from mets[].meta.smiles, DAG edges
//   from the MetaboliteLibrary step list (re-indexed into the slice's 0-based
//   namespace), and a parallel cofactor spec parsed per step via
//   AtomTrace.parseCofactorString. Byte-for-byte the v2.0.59 behaviour (the
//   cached fast lane), now just lifted into a helper.
//
//   LIVE user pathway — slice from mets[].meta.smiles, DAG edges from
//   _pathway.edges COLLAPSED through reaction-step nodes (reactant -> step ->
//   product becomes a direct reactant -> product trace edge) and re-indexed onto
//   the slice. When a connecting edge carries a stored edge.aam (written by
//   linkImportedStepAAM, in the SAME fromToIndex positional space), its mapping
//   is surfaced as edgeMaps[e] so tracePathway uses it verbatim; otherwise
//   tracePathway recomputes via its own mapReaction. No cofactor strings on the
//   live lane (user-drawn reactions don't carry them), so traceCofactors stays
//   all-empty.
//
// Returns { slice, traceEdges, traceCofactors, edgeMaps }. traceEdges / edgeMaps
// are re-indexed to the slice (0 == startMetaboliteIdx).
function buildPathwayTraceInputs(mets, startMetaboliteIdx) {
    var slice = [];
    for (var s = startMetaboliteIdx; s < mets.length; s++) {
        slice.push((mets[s].meta && mets[s].meta.smiles) || '');
    }
    var traceEdges = [];
    var traceCofactors = [];
    var edgeMaps = [];
    var id = _pathway && _pathway.loadedExampleId;

    if (id && typeof MetaboliteLibrary !== 'undefined' &&
            MetaboliteLibrary.pathways && MetaboliteLibrary.pathways[id]) {
        // ---- Example fast lane (unchanged DAG + cofactor derivation) ----
        var pwSteps = MetaboliteLibrary.pathways[id].steps || [];
        var nameToIdxRel = {};
        for (var ni = 0; ni < mets.length; ni++) {
            var canonName = mets[ni].name;
            var aliases = (mets[ni].meta && mets[ni].meta.aliases) || [];
            nameToIdxRel[String(canonName).toLowerCase()] = ni;
            for (var na = 0; na < aliases.length; na++) {
                var aliasKey = String(aliases[na]).toLowerCase();
                if (nameToIdxRel[aliasKey] === undefined) {
                    nameToIdxRel[aliasKey] = ni;
                }
            }
        }
        for (var psi = 0; psi < pwSteps.length; psi++) {
            var st = pwSteps[psi];
            var fIdx = nameToIdxRel[String(st.from).toLowerCase()];
            var tIdx = nameToIdxRel[String(st.to).toLowerCase()];
            if (fIdx === undefined || tIdx === undefined) { continue; }
            if (fIdx < startMetaboliteIdx || tIdx < startMetaboliteIdx) { continue; }
            traceEdges.push({ from: fIdx - startMetaboliteIdx, to: tIdx - startMetaboliteIdx });
            // Parse the step's cofactor string into a structured spec.
            if (typeof AtomTrace.parseCofactorString === 'function') {
                traceCofactors.push(AtomTrace.parseCofactorString(st.cofactors || ''));
            } else {
                traceCofactors.push({ in: [], out: [] });
            }
        }
        // Example lane never injects cached maps — it recomputes as before.
        edgeMaps = [];
    } else {
        // ---- Live user pathway lane ----
        var nodeIdToRel = {};
        for (var mi = 0; mi < mets.length; mi++) {
            if (mets[mi].nodeId) { nodeIdToRel[mets[mi].nodeId] = mi; }
        }
        var collapsed = collapseLivePathwayEdges();
        for (var ce = 0; ce < collapsed.length; ce++) {
            var c = collapsed[ce];
            var fRel = nodeIdToRel[c.from];
            var tRel = nodeIdToRel[c.to];
            if (fRel === undefined || tRel === undefined) { continue; }
            if (fRel < startMetaboliteIdx || tRel < startMetaboliteIdx) { continue; }
            traceEdges.push({ from: fRel - startMetaboliteIdx, to: tRel - startMetaboliteIdx });
            traceCofactors.push({ in: [], out: [] });
            // Prefer a stored mapping (edge.aam.mapping) — same positional space.
            edgeMaps.push((c.aam && c.aam.mapping && typeof c.aam.mapping === 'object') ?
                c.aam.mapping : null);
        }
        // If no edge carried a cached map, drop edgeMaps entirely so the trace
        // recomputes everything (and the all-null array doesn't masquerade as
        // an override set).
        var anyCached = false;
        for (var ek = 0; ek < edgeMaps.length; ek++) { if (edgeMaps[ek]) { anyCached = true; break; } }
        if (!anyCached) { edgeMaps = []; }
    }
    return { slice: slice, traceEdges: traceEdges, traceCofactors: traceCofactors, edgeMaps: edgeMaps };
}

// v2.4.2: collapse the live pathway's edges into structure-node -> structure-node
// trace edges. A direct edge between two structure nodes passes through as-is
// (carrying any edge.aam). An edge into a reaction-step node is paired with each
// edge out of that step to form reactant -> product trace edges, carrying the
// step's stored aam (preferring the OUT edge's aam, then the IN edge's). Returns
// [{ from: nodeId, to: nodeId, aam }]. Reaction-step nodes themselves are never
// trace endpoints.
function collapseLivePathwayEdges() {
    var out = [];
    if (!_pathway || !Array.isArray(_pathway.edges) || !Array.isArray(_pathway.nodes)) { return out; }
    var kindById = {};
    for (var n = 0; n < _pathway.nodes.length; n++) {
        if (_pathway.nodes[n]) { kindById[_pathway.nodes[n].id] = _pathway.nodes[n].kind; }
    }
    var edges = _pathway.edges;
    var e, ed;
    // Index step in/out edges.
    var stepIn = {}, stepOut = {};
    for (e = 0; e < edges.length; e++) {
        ed = edges[e];
        if (kindById[ed.to] === 'reaction') {
            if (!stepIn[ed.to]) { stepIn[ed.to] = []; }
            stepIn[ed.to].push(ed);
        }
        if (kindById[ed.from] === 'reaction') {
            if (!stepOut[ed.from]) { stepOut[ed.from] = []; }
            stepOut[ed.from].push(ed);
        }
    }
    // Direct structure -> structure edges.
    for (e = 0; e < edges.length; e++) {
        ed = edges[e];
        if (kindById[ed.from] !== 'reaction' && kindById[ed.to] !== 'reaction') {
            out.push({ from: ed.from, to: ed.to, aam: ed.aam || null });
        }
    }
    // Collapse reactant -> step -> product into reactant -> product.
    for (var stepId in stepIn) {
        if (!stepIn.hasOwnProperty(stepId)) { continue; }
        var ins = stepIn[stepId];
        var outs = stepOut[stepId] || [];
        for (var ii = 0; ii < ins.length; ii++) {
            for (var oo = 0; oo < outs.length; oo++) {
                var aam = outs[oo].aam || ins[ii].aam || null;
                out.push({ from: ins[ii].from, to: outs[oo].to, aam: aam });
            }
        }
    }
    return out;
}

function clearAtomTraceStrip() {
    var section = document.getElementById('wb-atom-trace-section');
    var grid = document.getElementById('wb-atom-trace-grid');
    if (grid) grid.textContent = '';
    if (section) section.hidden = true;
    if (typeof _pathway !== 'undefined' && _pathway) {
        _pathway.loadedExampleId = null;
        _pathway.loadedExamplePathwayMetabolites = null;
    }
}

// v2.0.22: react to an atom click anywhere in the trace strip.
// v2.0.54: when `shift` is true and the click lands on the START metabolite
// of the active trace, toggle the atom into the moiety selection instead
// of starting a new single-atom trace.
function handleTraceAtomClick(target, shift) {
    if (!target || !target.closest) return;
    var atomG = target.closest('.bime-trace-atom');
    if (!atomG) return;
    var cell = atomG.closest('[data-metabolite-index]');
    if (!cell) return;
    var atomIdx = parseInt(atomG.getAttribute('data-atom-index'), 10);
    var metIdx = parseInt(cell.getAttribute('data-metabolite-index'), 10);
    if (!isFinite(atomIdx) || !isFinite(metIdx)) return;
    // Shift-click on the active start metabolite toggles moiety membership.
    if (shift && _pathway.moietyStartMetIdx === metIdx) {
        toggleMoietyAtom(atomIdx);
        return;
    }
    // Shift-click on a non-start metabolite (or first click on a fresh
    // start) initialises moiety mode rooted at that metabolite.
    if (shift) {
        _pathway.moietyStartMetIdx = metIdx;
        _pathway.moietySet = [atomIdx];
        showAtomTraceStrip(metIdx, atomIdx);
        return;
    }
    // Plain click resets any moiety selection and runs single-atom trace.
    _pathway.moietyStartMetIdx = null;
    _pathway.moietySet = [];
    showAtomTraceStrip(metIdx, atomIdx);
}

// v2.0.54: toggle an atom in/out of the moiety set, then re-render.
function toggleMoietyAtom(atomIdx) {
    var set = _pathway.moietySet || [];
    var i = set.indexOf(atomIdx);
    if (i >= 0) {
        set.splice(i, 1);
    } else {
        set.push(atomIdx);
        set.sort(function(a, b) { return a - b; });
    }
    _pathway.moietySet = set;
    if (set.length === 0) {
        _pathway.moietyStartMetIdx = null;
        // Fall back to whichever atom was last single-traced — but if the
        // user cleared every moiety atom we just re-render the strip with
        // no active trace.
        showAtomTraceStrip();
    } else {
        // Single-atom trace stays anchored on the first moiety atom so the
        // UI has SOMETHING highlighted as the "primary".
        showAtomTraceStrip(_pathway.moietyStartMetIdx, set[0]);
    }
}

// v2.0.61: right-click context-menu state. Single floating <div> per
// workbench, repopulated + repositioned on each contextmenu event.
// `target` carries the {type, id} that was clicked (or {type:'canvas'}
// for an empty-canvas right-click).
var _pathwayContextMenu = null;

// Build + show the context menu at (clientX, clientY) for the given target.
// `items` is an array of { label, action, data, divider, disabled }.
function showPathwayContextMenu(items, clientX, clientY) {
    hidePathwayContextMenu();
    if (!items || !items.length) { return; }
    var menu = document.createElement('div');
    menu.id = 'wb-pathway-context-menu';
    menu.className = 'wb-context-menu';
    menu.setAttribute('role', 'menu');
    for (var i = 0; i < items.length; i++) {
        var item = items[i];
        if (item.divider) {
            var sep = document.createElement('div');
            sep.className = 'wb-context-menu-divider';
            sep.setAttribute('aria-hidden', 'true');
            menu.appendChild(sep);
            continue;
        }
        var btn = document.createElement('button');
        btn.type = 'button';
        btn.className = 'wb-context-menu-item';
        if (item.disabled) { btn.classList.add('is-disabled'); btn.disabled = true; }
        btn.setAttribute('role', 'menuitem');
        btn.setAttribute('data-wb-action', item.action);
        if (item.data) {
            for (var k in item.data) {
                if (item.data.hasOwnProperty(k)) {
                    btn.setAttribute('data-' + k, String(item.data[k]));
                }
            }
        }
        btn.textContent = item.label;
        menu.appendChild(btn);
    }
    // Provisional placement; clamp to viewport after measure.
    menu.style.left = clientX + 'px';
    menu.style.top = clientY + 'px';
    document.body.appendChild(menu);
    var rect = menu.getBoundingClientRect();
    var maxX = (window.innerWidth || rect.right) - rect.width - 6;
    var maxY = (window.innerHeight || rect.bottom) - rect.height - 6;
    if (clientX > maxX) { menu.style.left = Math.max(6, maxX) + 'px'; }
    if (clientY > maxY) { menu.style.top  = Math.max(6, maxY) + 'px'; }
    _pathwayContextMenu = menu;
}

function hidePathwayContextMenu() {
    if (_pathwayContextMenu && _pathwayContextMenu.parentNode) {
        _pathwayContextMenu.parentNode.removeChild(_pathwayContextMenu);
    }
    _pathwayContextMenu = null;
}

// Resolve the right-clicked pathway target from the event.
function resolvePathwayContextTarget(e) {
    var node = e.target;
    while (node && node !== document.body) {
        if (node.getAttribute && node.getAttribute('data-pathway-type')) {
            return {
                type: node.getAttribute('data-pathway-type'),
                id:   node.getAttribute('data-pathway-id')
            };
        }
        node = node.parentNode;
    }
    // Empty canvas — anchor the menu at the click point in world space.
    return { type: 'canvas' };
}

// Build the menu items list for a given target.
function pathwayContextMenuItems(target, worldX, worldY) {
    var items = [];
    if (target.type === 'node') {
        items.push({ label: 'Edit Structure',           action: 'context-edit-node-structure', data: { 'context-id': target.id } });
        items.push({ label: 'Update Labels from form', action: 'context-update-node-labels', data: { 'context-id': target.id } });
        items.push({ label: 'Duplicate node',           action: 'context-duplicate-node',     data: { 'context-id': target.id } });
        items.push({ label: 'Copy SMILES',              action: 'context-copy-smiles',        data: { 'context-id': target.id } });
        items.push({ divider: true });
        items.push({ label: 'Delete',                   action: 'context-delete-item',        data: { 'context-type': 'node', 'context-id': target.id } });
    } else if (target.type === 'edge') {
        items.push({ label: 'Reverse arrow',            action: 'context-reverse-edge',       data: { 'context-id': target.id } });
        items.push({ label: 'Cycle arrow type',         action: 'context-cycle-arrow-type',   data: { 'context-id': target.id } });
        items.push({ divider: true });
        items.push({ label: 'Delete',                   action: 'context-delete-item',        data: { 'context-type': 'edge', 'context-id': target.id } });
    } else if (target.type === 'step' || target.type === 'note' || target.type === 'compartment') {
        items.push({ label: 'Update Labels from form', action: 'context-update-node-labels', data: { 'context-type': target.type, 'context-id': target.id } });
        items.push({ divider: true });
        items.push({ label: 'Delete',                   action: 'context-delete-item',        data: { 'context-type': target.type, 'context-id': target.id } });
    } else {
        // Empty canvas — drop-here actions + canvas-level commands.
        items.push({ label: 'Add Metabolite here',     action: 'context-add-here', data: { 'context-kind': 'metabolite', 'context-x': worldX, 'context-y': worldY } });
        items.push({ label: 'Add Residue here',         action: 'context-add-here', data: { 'context-kind': 'residue',    'context-x': worldX, 'context-y': worldY } });
        items.push({ label: 'Add Step here',            action: 'context-add-here', data: { 'context-kind': 'step',       'context-x': worldX, 'context-y': worldY } });
        items.push({ label: 'Add Note here',            action: 'context-add-here', data: { 'context-kind': 'note',       'context-x': worldX, 'context-y': worldY } });
        items.push({ divider: true });
        items.push({ label: 'Clean Up Layout',          action: 'pathway-layout' });
        items.push({ label: 'Fit',                       action: 'pathway-fit' });
    }
    return items;
}

function handlePathwayCanvasContextMenu(e) {
    var svg = document.getElementById('wb-pathway-canvas');
    if (!svg || !svg.contains(e.target)) { return; }
    e.preventDefault();
    var target = resolvePathwayContextTarget(e);
    if (target.type !== 'canvas') {
        // Auto-select the right-clicked item so subsequent context
        // actions (Edit Structure, Update Labels, …) act on it.
        selectPathwayItem(target.type, target.id);
        renderPathwayCanvas();
    }
    var world = pathwayEventPoint(e);
    var items = pathwayContextMenuItems(target, world.x, world.y);
    showPathwayContextMenu(items, e.clientX, e.clientY);
}

function handlePathwayCanvasDoubleClick(e) {
    var svg = document.getElementById('wb-pathway-canvas');
    if (!svg || !svg.contains(e.target)) { return; }
    var target = resolvePathwayContextTarget(e);
    if (target.type !== 'node' || !target.id) { return; }
    e.preventDefault();
    selectPathwayItem('node', target.id);
    renderPathwayCanvas();
    // Merge Stage 2: prefer the in-place focus lens (overlay editor anchored
    // over the node). Fall back to the legacy bridge (load into the main editor)
    // when the lens is disabled — e.g. an older bundle without FocusLens.
    if (!(_lens && openStructureLens(target.id))) {
        editPathwayNodeStructure();
    }
}

// Action handlers for the context menu — most delegate to existing
// editing primitives.
function contextEditNodeStructure(id) {
    if (id) { selectPathwayItem('node', id); }
    editPathwayNodeStructure();
}
function contextUpdateNodeLabels(type, id) {
    if (type && id) { selectPathwayItem(type, id); }
    updatePathwaySelectionFromControls();
}
function contextDeleteItem(type, id) {
    if (type && id) { selectPathwayItem(type, id); }
    deletePathwaySelection();
}
function contextCopySmiles(id) {
    var node = findPathwayItem('node', id);
    var structure = pathwayNodeStructurePayload(node);
    var smiles = (structure && structure.smiles) || (node && node.subtitle) || '';
    if (!node || !smiles) {
        pathwayStatus('No SMILES on this node — draw a structure first.');
        return;
    }
    if (navigator && navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(smiles).then(function() {
            pathwayStatus('Copied SMILES to clipboard: ' + smiles);
        }, function() {
            pathwayStatus('Could not copy SMILES (clipboard refused). SMILES: ' + smiles);
        });
    } else {
        pathwayStatus('Clipboard unavailable. SMILES: ' + smiles);
    }
}
function contextDuplicateNode(id) {
    var node = findPathwayItem('node', id);
    if (!node) { return; }
    pushPathwayHistory && pushPathwayHistory();
    var copy = addPathwayNode(node.x + 40, node.y + 40, {
        label: node.label,
        subtitle: node.subtitle,
        kind: node.kind,
        imageHref: node.imageHref,
        structure: pathwayNodeStructurePayload(node)
    });
    if (copy) { pathwayStatus('Duplicated "' + node.label + '" → ' + copy.id + '.'); }
}
function contextReverseEdge(id) {
    var edge = findPathwayItem('edge', id);
    if (!edge) { return; }
    pushPathwayHistory && pushPathwayHistory();
    var tmp = edge.from; edge.from = edge.to; edge.to = tmp;
    renderPathwayCanvas();
    pathwayStatus('Reversed arrow ' + id + '.');
}
function contextCycleArrowType(id) {
    var edge = findPathwayItem('edge', id);
    if (!edge) { return; }
    pushPathwayHistory && pushPathwayHistory();
    var order = (typeof PATHWAY_ARROW_TYPES !== 'undefined') ? PATHWAY_ARROW_TYPES :
        ['forward', 'reverse', 'reversible', 'resonance'];
    var idx = order.indexOf(edge.arrowType || 'forward');
    edge.arrowType = order[(idx + 1) % order.length];
    renderPathwayCanvas();
    pathwayStatus('Arrow type → ' + edge.arrowType + '.');
}
function contextAddHere(kind, x, y) {
    pushPathwayHistory && pushPathwayHistory();
    var nx = +x, ny = +y;
    if (kind === 'metabolite') {
        addPathwayNode(nx, ny, { label: pathwayInputLabel('Metabolite'), kind: 'metabolite' });
    } else if (kind === 'residue') {
        addPathwayNode(nx, ny, { label: pathwayInputLabel('Residue'), kind: 'residue' });
    } else if (kind === 'step') {
        addPathwayStep(nx, ny, pathwayInputLabel('Mechanism step'));
    } else if (kind === 'note') {
        addPathwayNote(nx, ny, pathwayInputLabel('Comment'));
    }
}

// ---------------------------------------------------------------------------
// Focus-lens glue (merge Stage 2). Double-clicking a structure node expands it
// IN PLACE: the always-mounted molecule editor is repositioned as a floating
// overlay panel over the node, pre-loaded with the node's structure. Escape or
// click-away captures the edited structure back onto the node (re-rasterising
// its thumbnail) and collapses the overlay. We never leave the canvas and never
// create/destroy the editor — this reuses the existing capture/apply bridge.
//
// _lens (module global, near `var editor`) owns ONLY which node is focused
// (editor/FocusLens.js, DOM-free). These functions are the DOM glue: they load
// the node into the editor, position the overlay, and commit on collapse.
//
// v2.3.4 (focus-lens trust + accessibility) hardens this interaction: the lens
// is now an explicit, accessible DIALOG. The small helpers below build its
// Done/Cancel header (P0-2), manage the Tab focus-trap (P0-3), and flash the
// node on save (P0-2). openStructureLens / collapseStructureLens delegate to
// them so those two functions stay single, un-split, source-shape-stable.
// ---------------------------------------------------------------------------

// P0-2: build the lens overlay's header bar (node label + Done + Cancel) and
// prepend it inside `.wb-editor-wrap`, above the editor chrome. Done commits
// (collapseStructureLens(true)); Cancel discards (collapseStructureLens(false)).
// Any stale header is removed first so re-opening / hopping never doubles it.
function buildStructureLensHeader(node) {
    var wrap = document.querySelector('.wb-editor-wrap');
    if (!wrap) { return null; }
    removeStructureLensHeader();
    var label = (node && node.label) ? String(node.label) : 'structure';
    var header = document.createElement('div');
    header.className = 'wb-lens-header';
    var title = document.createElement('span');
    title.className = 'wb-lens-title';
    title.textContent = 'Edit structure: ' + label;
    var actions = document.createElement('div');
    actions.className = 'wb-lens-actions';
    var cancel = document.createElement('button');
    cancel.type = 'button';
    cancel.className = 'wb-lens-btn wb-lens-cancel';
    cancel.textContent = 'Cancel';
    cancel.addEventListener('click', function() { collapseStructureLens(false); });
    var done = document.createElement('button');
    done.type = 'button';
    done.className = 'wb-lens-btn wb-lens-done';
    done.textContent = 'Done';
    done.addEventListener('click', function() { collapseStructureLens(true); });
    actions.appendChild(cancel);
    actions.appendChild(done);
    header.appendChild(title);
    header.appendChild(actions);
    if (wrap.firstChild) { wrap.insertBefore(header, wrap.firstChild); }
    else { wrap.appendChild(header); }
    return header;
}

// P0-2: tear down the lens header on collapse (idempotent).
function removeStructureLensHeader() {
    var existing = document.querySelector('.wb-editor-wrap .wb-lens-header');
    if (existing && existing.parentNode) { existing.parentNode.removeChild(existing); }
}

// v2.3.8 (P1-4): a COMPACT live readout strip inside the lens. While a focus
// lens is open the Inspector / Editor Output dock BELOW the canvas — off-screen
// on a phone (where the lens is a full-screen sheet) and easy to miss on
// desktop — so a user drawing in the lens could not see the SMILES or warnings
// their edit produced without closing the lens (which commits). This footer bar
// (appended at the BOTTOM of `.wb-editor-wrap`, below the editor host) shows the
// current SMILES (truncated with an ellipsis when long) and a warning-COUNT
// badge, both refreshed from the SAME change->output path that already drives
// Editor Output (updateLensReadout, called at the tail of updateOutputNow — no
// extra recompute). The whole strip is aria-live="polite" so a screen-reader
// user editing in the canvas hears the SMILES / warning updates. Built here in
// openStructureLens; torn down in collapseStructureLens (mirrors the v2.3.4
// header). Any stale strip is removed first so re-opening / hopping never
// doubles it.
function buildStructureLensReadout(node) {
  var wrap = document.querySelector('.wb-editor-wrap');
  if (!wrap) { return null; }
  removeStructureLensReadout();
  var bar = document.createElement('div');
  bar.className = 'wb-lens-readout';
  // Polite live region: announce SMILES / warning changes without interrupting
  // the user's drawing keystrokes (WCAG 4.1.3 status messages).
  bar.setAttribute('aria-live', 'polite');
  bar.setAttribute('aria-atomic', 'true');
  // SMILES line. A button so it is keyboard-reachable + announces as actionable;
  // clicking copies the full (untruncated) SMILES via the existing copy path.
  var smiles = document.createElement('button');
  smiles.type = 'button';
  smiles.className = 'wb-lens-readout-smiles';
  smiles.title = 'Copy SMILES';
  smiles.textContent = '—';
  smiles.addEventListener('click', function() { copyLensReadoutSmiles(); });
  // Warning-count badge. Text glyph + count (no colour-only meaning, no emoji).
  var warn = document.createElement('span');
  warn.className = 'wb-lens-readout-warn';
  warn.textContent = '0 warnings';
  bar.appendChild(smiles);
  bar.appendChild(warn);
  wrap.appendChild(bar);
  // Fill it immediately so the strip is correct the instant the lens opens
  // (the just-loaded node structure), not only after the next edit.
  updateLensReadout();
  return bar;
}

// v2.3.8 (P1-4): tear down the live readout strip on collapse (idempotent).
function removeStructureLensReadout() {
  var existing = document.querySelector('.wb-editor-wrap .wb-lens-readout');
  if (existing && existing.parentNode) { existing.parentNode.removeChild(existing); }
}

// v2.3.8 (P1-4): copy the CURRENT (full, untruncated) SMILES the readout is
// showing. Reuses the editor's own smiles() + the existing clipboard pattern
// (copyOut copies #smiles-out; here the source is the live molecule so the copy
// works even with the SMILES tab not selected). A no-op without clipboard.
function copyLensReadoutSmiles() {
  if (!editor || !editor.molecule) { return; }
  var text = '';
  try { text = outputTextFor('smiles') || ''; } catch (e) { text = ''; }
  if (!text) { return; }
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(function() {
      if (editor && editor.showInfo) { editor.showInfo('Copied SMILES'); }
    });
  }
}

// v2.3.8 (P1-4): refresh the lens readout from the SAME cached signals Editor
// Output already computed this frame — outputTextFor('smiles') (cached per chem
// signature) and getMoleculeInsights().warnings (cached per perf signature).
// No extra parse / layout work: this is a pure DOM write. A no-op unless a lens
// is open and the strip exists. Called at the tail of updateOutputNow (the one
// RAF-debounced per-change update path) so it tracks every draw automatically.
function updateLensReadout() {
  if (!_lens || !_lens.isOpen()) { return; }
  var bar = document.querySelector('.wb-editor-wrap .wb-lens-readout');
  if (!bar) { return; }
  var smilesEl = bar.querySelector('.wb-lens-readout-smiles');
  var warnEl = bar.querySelector('.wb-lens-readout-warn');
  var mol = editor && editor.molecule;
  var empty = !mol || !mol.atoms || mol.atoms.length === 0;
  var smiles = '';
  if (!empty) {
    try { smiles = outputTextFor('smiles') || ''; } catch (e) { smiles = ''; }
  }
  if (smilesEl) {
    if (empty) {
      smilesEl.textContent = 'No structure yet';
      smilesEl.setAttribute('disabled', 'disabled');
      smilesEl.removeAttribute('title');
    } else {
      smilesEl.textContent = smiles ? truncateMiddle(smiles, 48) : '(no SMILES)';
      smilesEl.removeAttribute('disabled');
      smilesEl.title = smiles ? ('Copy SMILES: ' + smiles) : 'Copy SMILES';
    }
  }
  if (warnEl) {
    var warnings = [];
    try { warnings = (getMoleculeInsights() || {}).warnings || []; } catch (e2) { warnings = []; }
    var n = warnings.length;
    // Non-colour cue: a triangular alert glyph + the literal count + the word
    // "warning(s)" so the badge reads in greyscale and to a screen reader.
    warnEl.textContent = n === 0 ? '0 warnings'
      : ('⚠ ' + n + ' warning' + (n === 1 ? '' : 's'));
    warnEl.classList.toggle('has-warnings', n > 0);
    warnEl.setAttribute('title', n === 0 ? 'No structure warnings'
      : (n + ' structure warning' + (n === 1 ? '' : 's') + ' — see Editor Output > Warnings'));
  }
}

// v2.3.8 (P1-4): elide the MIDDLE of a long SMILES so both the head (atoms) and
// tail stay visible in the narrow strip — more legible than a trailing cut for
// a chemist eyeballing a structure. Short strings pass through unchanged.
function truncateMiddle(s, max) {
  s = String(s == null ? '' : s);
  if (s.length <= max) { return s; }
  var keep = max - 1;
  var head = Math.ceil(keep / 2);
  var tail = Math.floor(keep / 2);
  return s.slice(0, head) + '…' + s.slice(s.length - tail);
}

// P0-3: the overlay's focusable elements, in DOM order — the Done/Cancel buttons
// plus any focusable editor chrome — used for the initial focus move and the Tab
// trap. Defensive against a missing wrap (returns []).
function lensFocusableElements() {
    var wrap = document.querySelector('.wb-editor-wrap');
    if (!wrap || !wrap.querySelectorAll) { return []; }
    var nodes = wrap.querySelectorAll(
        'button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])');
    var out = [];
    for (var i = 0; i < nodes.length; i++) {
        var el = nodes[i];
        if (el.disabled) { continue; }
        if (el.getAttribute && el.getAttribute('tabindex') === '-1') { continue; }
        out.push(el);
    }
    return out;
}

// P0-3: Tab/Shift+Tab cycle WITHIN the overlay while the lens is open, so focus
// never escapes behind the modal dialog (WCAG 2.4.3 / 4.1.2). Bound on the wrap
// in openStructureLens (kept in `_lensTrapHandler`) and removed on collapse.
function onLensTrapKeydown(e) {
    if (e.key !== 'Tab' && e.keyCode !== 9) { return; }
    var focusables = lensFocusableElements();
    if (!focusables.length) { return; }
    var first = focusables[0];
    var last = focusables[focusables.length - 1];
    var active = (typeof document !== 'undefined') ? document.activeElement : null;
    if (e.shiftKey) {
        if (active === first || focusables.indexOf(active) === -1) {
            if (last && last.focus) { last.focus(); }
            e.preventDefault();
        }
    } else {
        if (active === last || focusables.indexOf(active) === -1) {
            if (first && first.focus) { first.focus(); }
            e.preventDefault();
        }
    }
}

// P0-2: brief save-confirmation flash on the node group once an edit commits.
// Adds `.wb-node-saved-flash` for ~900ms (the CSS animation honours
// prefers-reduced-motion). Guarded so a missing node element is a no-op.
function flashSavedNode(nodeId) {
    var el = document.querySelector('[data-pathway-id="' + nodeId + '"]');
    if (!el || !el.classList) { return; }
    el.classList.add('wb-node-saved-flash');
    if (typeof setTimeout === 'function') {
        setTimeout(function() {
            if (el.classList) { el.classList.remove('wb-node-saved-flash'); }
        }, 900);
    }
}

// P1-1 (v2.3.4): de-eager the click-away. The old handler committed + collapsed
// on ANY mousedown outside the overlay/focused node — so a stray click on the
// inspector or output panel silently destroyed the focus session. Now, while a
// lens is open, a mousedown is routed by WHERE it lands:
//   - on ANOTHER editable structure node  -> commit this lens, then HOP straight
//     into that node's lens (resolved via resolvePathwayContextTarget, the same
//     way double-click / context menu resolve a node). A direct hop, not a
//     dismiss.
//   - inside the overlay (drawing atoms / the lens header), the focused node, OR
//     other known app chrome (the Inspector, the Editor Output panel) -> IGNORE.
//   - only a TRUE click on the bare canvas / page background commits + collapses
//     (the original click-away, preserved).
function handleLensClickAway(e) {
    if (!_lens || !_lens.isOpen()) { return; }
    var hasClosest = e.target && e.target.closest;
    // Mousedowns inside the overlay editor (which includes the lens header) are
    // the user working in the lens — never dismiss.
    if (hasClosest && e.target.closest('.wb-editor-wrap')) { return; }
    var focusedId = _lens.focused();
    // A mousedown on the focused node itself is not a dismiss either.
    if (hasClosest && e.target.closest('[data-pathway-id="' + focusedId + '"]')) { return; }
    // A mousedown on ANOTHER editable structure node hops the lens there.
    var target = resolvePathwayContextTarget(e);
    if (target && target.type === 'node' && target.id && target.id !== focusedId) {
        collapseStructureLens(true);
        openStructureLens(target.id);
        return;
    }
    // Mousedowns inside other known app chrome (inspector / output panel) must
    // NOT dismiss — only the bare canvas / page background does.
    if (hasClosest && e.target.closest('.wb-inspector, .wb-output-section')) { return; }
    collapseStructureLens(true);
}

function openStructureLens(nodeId) {
    if (!_lens) { return false; }
    var node = findPathwayItem('node', nodeId);
    if (!node) { return false; }
    selectPathwayItem('node', nodeId);
    renderPathwayCanvas();
    // Load the node's saved chemistry into the editor (RXN -> MOL -> SMILES),
    // preserving its title/comment, exactly like editPathwayNodeStructure.
    var source = pathwayNodeStructureInput(node);
    if (source) {
        try {
            load(source);
            var st = pathwayNodeStructurePayload(node);
            if (st) {
                if (st.title && editor.setMolName) { editor.setMolName(st.title); }
                if (editor.molecule) { editor.molecule.comment = st.comment || ''; }
            }
            if (editor.render) { editor.render(); }
        } catch (e) {}
    }
    // Merge Stage 4 (v2.3.2): give the focus session a CLEAN editor-history
    // baseline. load() -> readGenericMolecularInput() calls saveHistory() before
    // it clears, so the editor undo stack would otherwise carry the PREVIOUS
    // molecule (or a prior focus session). Clearing here means session Ctrl+Z
    // bottoms out at the just-loaded structure and can't escape into another
    // node's chemistry. (The blank-node path from newPathwayNodeInLens has no
    // source to load; clearing still drops any stale prior-molecule history.)
    if (editor && editor.history && editor.history.clear) { editor.history.clear(); }
    _lens.open(nodeId);
    var main = document.querySelector('.wb-main');
    if (main && main.classList) { main.classList.add('wb-lens-open'); }
    var nodeEl = document.querySelector('[data-pathway-id="' + nodeId + '"]');
    if (nodeEl && nodeEl.classList) { nodeEl.classList.add('is-lens-focused'); }
    // v2.3.4 (P0-2): give the overlay an explicit Done/Cancel header so the lens
    // has a visible commit + a real cancel (no longer destructive-by-default with
    // no escape hatch). Built before positioning so the panel sizes around it.
    buildStructureLensHeader(node);
    // v2.3.8 (P1-4): a compact live SMILES + warning-count strip at the BOTTOM
    // of the overlay so the user reads what their drawing produces WITHOUT
    // closing the lens (the Inspector / Editor Output is docked off-screen below
    // the canvas). Built after the header so the panel sizes around both; fed by
    // the existing change->output path (updateLensReadout in updateOutputNow).
    buildStructureLensReadout(node);
    var wrap = document.querySelector('.wb-editor-wrap');
    // v2.3.4 (P0-3): the overlay is now an accessible modal dialog (WCAG 4.1.2).
    // Record the pre-lens focus to restore on collapse (WCAG 2.4.3 — mirrors the
    // footer modal's _modalLastFocus), label the dialog, then move focus INTO it
    // (the Done button) and trap Tab so focus can't slip behind the overlay.
    _lensLastFocus = (typeof document !== 'undefined') ? document.activeElement : null;
    if (wrap && wrap.setAttribute) {
        wrap.setAttribute('role', 'dialog');
        wrap.setAttribute('aria-modal', 'true');
        wrap.setAttribute('aria-label', 'Edit structure: ' + node.label);
    }
    positionStructureLens();
    var doneBtn = wrap ? wrap.querySelector('.wb-lens-done') : null;
    if (doneBtn && doneBtn.focus) { try { doneBtn.focus(); } catch (eF) {} }
    if (wrap && wrap.addEventListener) {
        _lensTrapHandler = onLensTrapKeydown;
        wrap.addEventListener('keydown', _lensTrapHandler, true);
    }
    pathwayStatus('Editing "' + node.label + '" — Done (or Esc) saves; Cancel discards.');
    return true;
}

// Merge Stage 3 (v2.3.1): the single-molecule entry point for the continuum.
// "Draw Structure" with an empty editor drops a FRESH blank structure node onto
// the canvas and opens the focus lens on it, so the user draws one standalone
// structure and reads its SMILES/MOL/RXN in the lens — the editor output stays
// available while focused, exactly like the retired standalone Molecule mode.
// Reuses addPathwayNode + openStructureLens; the empty editor is cleared first
// so the lens opens onto a blank drawing surface. Returns true when the lens
// opened (so the caller can short-circuit its legacy fallback).
function newPathwayNodeInLens() {
    if (!_lens) { return false; }
    if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    // Start from a clean editor so the new node's lens opens onto a blank
    // drawing surface (the empty node carries no structure to load). The caller
    // only reaches here when the editor is already empty; clearing is defensive
    // so the helper is correct if reused.
    if (editor && editor.molecule && editor.molecule.clear) {
        try {
            editor.molecule.clear();
            if (editor.render) { editor.render(); }
            queueOutputUpdate();
        } catch (e) {}
    }
    var center = pathwayVisibleCenter();
    var node = addPathwayNode(center.x, center.y, {
        label: 'Structure',
        kind: 'metabolite',
        silent: true
    });
    if (!node) { return false; }
    if (!openStructureLens(node.id)) { return false; }
    pathwayStatus('New structure node — draw it in the lens; its SMILES/MOL is in Editor Output. Esc or click away to save it onto the node.');
    return true;
}

function collapseStructureLens(commit) {
    if (!_lens || !_lens.isOpen()) { return; }
    var id = _lens.focused();
    var node = findPathwayItem('node', id);
    // v2.3.4 (P0-2): commit applies + confirms; CANCEL (commit === false) DISCARDS
    // — the whole pushPathwayHistory + applyPathwayStructurePayload block is
    // guarded on `commit`, so a cancelled session leaves the node's stored
    // structure (and the pathway history) completely untouched.
    var committed = false;
    if (commit && node) {
        var payload = captureEditorStructureForPathway();
        if (payload) {
            pushPathwayHistory();
            applyPathwayStructurePayload(node, payload);
            committed = true;
        }
    }
    _lens.close();
    var main = document.querySelector('.wb-main');
    if (main && main.classList) { main.classList.remove('wb-lens-open'); }
    // v2.3.6 (P0-6): drop the full-screen sheet layout class (if a mobile lens was
    // open) and restore the mobile dock the sheet had hidden.
    if (main && main.classList) { main.classList.remove('wb-lens-sheet'); }
    setMobileDockHidden(false);
    var nodeEl = document.querySelector('[data-pathway-id="' + id + '"]');
    if (nodeEl && nodeEl.classList) { nodeEl.classList.remove('is-lens-focused'); }
    // v2.3.4 (P0-2): tear down the Done/Cancel header. (P0-3): drop the dialog
    // ARIA + the Tab-trap and restore focus to whatever held it before the lens
    // opened (WCAG 2.4.3), guarded in case that element is gone.
    removeStructureLensHeader();
    // v2.3.8 (P1-4): tear down the live SMILES/warning readout strip too.
    removeStructureLensReadout();
    // Clear the inline positioning so the editor returns to its normal flow
    // (and stays display:none under .is-pathway-open) once the lens collapses.
    var wrap = document.querySelector('.wb-editor-wrap');
    if (wrap && wrap.style) {
        wrap.style.position = '';
        wrap.style.left = '';
        wrap.style.top = '';
        wrap.style.width = '';
        wrap.style.height = '';
        wrap.style.zIndex = '';
    }
    if (wrap && wrap.removeAttribute) {
        wrap.removeAttribute('role');
        wrap.removeAttribute('aria-modal');
        wrap.removeAttribute('aria-label');
    }
    if (wrap && _lensTrapHandler && wrap.removeEventListener) {
        wrap.removeEventListener('keydown', _lensTrapHandler, true);
    }
    _lensTrapHandler = null;
    if (_lensLastFocus && _lensLastFocus.focus) {
        try { _lensLastFocus.focus(); } catch (eR) {}
    }
    _lensLastFocus = null;
    renderPathwayCanvas();
    // v2.3.4 (P0-2): confirm the save AFTER the re-render (the flash class is set
    // on the freshly-drawn node group). A discard stays silent.
    if (committed) {
        flashSavedNode(id);
        pathwayStatus('Saved "' + node.label + '".');
    }
}

// v2.3.6 (P0-6): toggle the full-screen lens SHEET layout. positionStructureLens
// delegates here so the pinned function stays readable. When `on` is true we let
// CSS own the sizing (`.wb-lens-sheet` on `.wb-main`, sized in the <=720px block),
// strip any inline left/top/width/height the desktop anchor path wrote (so they
// can't fight the sheet), hide the mobile dock, and refit the editor draw surface.
// When `on` is false we simply drop the class + restore the dock; the caller's
// desktop anchor math then writes fresh inline coords. Idempotent.
function applyLensSheetLayout(wrap, on) {
  var main = document.querySelector('.wb-main');
  if (on) {
    if (main && main.classList) { main.classList.add('wb-lens-sheet'); }
    // Hand sizing to CSS: clear the inline coords a prior desktop layout set.
    wrap.style.position = '';
    wrap.style.left = '';
    wrap.style.top = '';
    wrap.style.width = '';
    wrap.style.height = '';
    wrap.style.zIndex = '';
    // The sheet covers the screen; keep the dock out of its way.
    setMobileDockHidden(true);
  } else {
    if (main && main.classList) { main.classList.remove('wb-lens-sheet'); }
    setMobileDockHidden(false);
  }
  // Re-fit the editor's draw surface to the new panel size (retry across a few
  // frames so a slow layout pass can't leave the canvas under-fitted), mirroring
  // the desktop refit at the tail of positionStructureLens.
  if (editor && editor.setSize) {
    var host = document.getElementById('bime-editor');
    var attempts = 0;
    var apply = function() {
      if (host && host.clientWidth > 0) {
        editor.setSize(host.clientWidth, host.clientHeight);
        if (editor.render) { editor.render(); }
        return;
      }
      if (attempts++ < 3 && typeof requestAnimationFrame === 'function') {
        requestAnimationFrame(apply);
        return;
      }
      if (editor.render) { editor.render(); }
    };
    if (typeof requestAnimationFrame === 'function') { requestAnimationFrame(apply); }
    else { apply(); }
  }
}

function positionStructureLens() {
    if (!_lens || !_lens.isOpen()) { return; }
    var wrap = document.querySelector('.wb-editor-wrap');
    if (!wrap || !wrap.style) { return; }
    var id = _lens.focused();
    var node = findPathwayItem('node', id);
    // v2.3.6 (P0-6): on a narrow (phone) viewport, present the editor as a
    // FULL-SCREEN SHEET instead of the node-anchored floating panel. The anchored
    // panel clamped to a 380x320 min, which on a ~360px phone covered the node,
    // overflowed, and offered no on-screen exit (Esc is unavailable on a touch
    // keyboard). Here we DETECT the narrow viewport, hand sizing to CSS via the
    // `.wb-lens-sheet` class on `.wb-main`, clear any inline desktop coords from a
    // prior wider layout, refit the editor, hide the dock so it can't fight the
    // sheet, and SKIP the node-anchoring/flip math below. The v2.3.4 lens contract
    // (sticky Done/Cancel header + role=dialog/aria-modal/focus-trap set in
    // openStructureLens) still applies to the sheet. Desktop behaviour below is
    // unchanged above the breakpoint.
    if (isNarrowViewport()) {
        applyLensSheetLayout(wrap, true);
        return;
    }
    // Wider viewport: ensure the sheet layout is off (e.g. the user resized a
    // phone-width window back to desktop while the lens stayed open) and the
    // dock is restored, then fall through to the desktop anchor/flip math.
    applyLensSheetLayout(wrap, false);
    // Anchor rect in CSS pixels (fixed/viewport coords). Prefer the rendered
    // node element's own bounding box; fall back to the FocusLens model rect
    // mapped through the canvas SVG's viewBox->client ratio when the DOM rect
    // is unavailable (e.g. node off-screen / not yet painted).
    var anchor = null;
    var nodeEl = document.querySelector('[data-pathway-id="' + id + '"]');
    if (nodeEl && nodeEl.getBoundingClientRect) {
        var nb = nodeEl.getBoundingClientRect();
        if (nb && (nb.width > 0 || nb.height > 0)) {
            anchor = { left: nb.left, top: nb.top, width: nb.width, height: nb.height };
        }
    }
    if (!anchor && node) {
        var svg = document.getElementById('wb-pathway-canvas');
        var svgRect = svg && svg.getBoundingClientRect ? svg.getBoundingClientRect() : null;
        var model = (typeof FocusLens !== 'undefined' && FocusLens.lensRectForNode)
            ? FocusLens.lensRectForNode(_pathway, node)
            : { x: node.x, y: node.y, w: node.w, h: node.h };
        if (svgRect && model) {
            // viewBox is "0 0 1200 620" (see workbench.html); map viewBox px to
            // client px by the rendered/viewBox width ratio.
            var vbW = 1200, vbH = 620;
            var sx = svgRect.width / vbW;
            var sy = svgRect.height / vbH;
            anchor = {
                left: svgRect.left + model.x * sx,
                top: svgRect.top + model.y * sy,
                width: model.w * sx,
                height: model.h * sy
            };
        }
    }
    var vw = (typeof window !== 'undefined' && window.innerWidth) ? window.innerWidth : 1280;
    var vh = (typeof window !== 'undefined' && window.innerHeight) ? window.innerHeight : 800;
    if (!anchor) {
        // Last resort: centre a default-sized panel in the viewport.
        anchor = { left: vw / 2 - 240, top: vh / 2 - 200, width: 480, height: 400 };
    }
    // Comfortable editing size: at least 380x320, growing toward the node's
    // scaled on-screen size, but never larger than the viewport (with margin).
    var MARGIN = 12;
    var MIN_W = 380, MIN_H = 320;
    var maxW = Math.max(MIN_W, vw - 2 * MARGIN);
    var maxH = Math.max(MIN_H, vh - 2 * MARGIN);
    var w = Math.min(maxW, Math.max(MIN_W, anchor.width));
    var h = Math.min(maxH, Math.max(MIN_H, anchor.height));
    // Anchor the panel over/near the node: start aligned to the node's top-left.
    // Merge Stage 5 (v2.3.3): if that overflows the viewport edge, FLIP the panel
    // to the opposite side of the node instead of letting the clamp drag it back
    // OVER the node — a node near the right edge gets the panel to its LEFT
    // (anchor.left - w - GAP), a node near the bottom gets it ABOVE
    // (anchor.top - h - GAP). The min-size + viewport clamp below is still the
    // final safety (e.g. a node so wide the flip also overflows).
    var GAP = 8;
    var left = anchor.left;
    var top = anchor.top;
    if (left + w > vw - MARGIN) { left = anchor.left - w - GAP; }
    if (top + h > vh - MARGIN) { top = anchor.top - h - GAP; }
    // Final safety: keep the whole panel inside the viewport.
    if (left + w > vw - MARGIN) { left = vw - MARGIN - w; }
    if (top + h > vh - MARGIN) { top = vh - MARGIN - h; }
    if (left < MARGIN) { left = MARGIN; }
    if (top < MARGIN) { top = MARGIN; }
    wrap.style.position = 'fixed';
    wrap.style.left = Math.round(left) + 'px';
    wrap.style.top = Math.round(top) + 'px';
    wrap.style.width = Math.round(w) + 'px';
    wrap.style.height = Math.round(h) + 'px';
    wrap.style.zIndex = '50';
    // Re-fit the editor's draw surface to the panel's inner size, mirroring
    // refitMoleculeEditorSoon (retry across a few frames so a slow layout pass
    // can't leave the canvas under-fitted), then render.
    if (editor && editor.setSize) {
        var host = document.getElementById('bime-editor');
        var attempts = 0;
        var apply = function() {
            if (host && host.clientWidth > 0) {
                editor.setSize(host.clientWidth, host.clientHeight);
                if (editor.render) { editor.render(); }
                return;
            }
            if (attempts++ < 3 && typeof requestAnimationFrame === 'function') {
                requestAnimationFrame(apply);
                return;
            }
            if (editor.render) { editor.render(); }
        };
        if (typeof requestAnimationFrame === 'function') { requestAnimationFrame(apply); }
        else { apply(); }
    } else if (editor && editor.render) {
        editor.render();
    }
}

// v2.0.58: marquee (drag-rectangle) moiety selection state.
// _atomTraceDrag is non-null only between mousedown and mouseup over an
// atom-trace cell SVG. When the user releases the mouse after moving more
// than a few pixels, every atom whose hit-target falls inside the rect is
// added to the active moiety (or, with Shift, EXTENDED into it).
var _atomTraceDrag = null;
// Suppress the click that fires immediately after a marquee drag so we
// don't accidentally re-trigger handleTraceAtomClick on top of the
// drag-finalised selection.
var _suppressNextAtomTraceClick = false;
var DRAG_THRESHOLD_PX = 4;

function handleAtomTraceMouseDown(e) {
    // Only react to primary-button mousedowns inside an atom-trace cell SVG.
    if (e.button !== 0) { return; }
    var holder = e.target && e.target.closest ? e.target.closest('.wb-atom-trace-cell-svg') : null;
    if (!holder) { return; }
    var svg = holder.querySelector('svg.bime-trace-svg');
    var cell = holder.closest('[data-metabolite-index]');
    if (!svg || !cell) { return; }
    var cellIdx = parseInt(cell.getAttribute('data-metabolite-index'), 10);
    if (!isFinite(cellIdx)) { return; }
    var rect = svg.getBoundingClientRect();
    _atomTraceDrag = {
        cellIdx: cellIdx,
        svg: svg,
        startX: e.clientX - rect.left,
        startY: e.clientY - rect.top,
        currentX: e.clientX - rect.left,
        currentY: e.clientY - rect.top,
        shiftKey: !!e.shiftKey,
        rectEl: null,
        moved: false
    };
}

function handleAtomTraceMouseMove(e) {
    if (!_atomTraceDrag) { return; }
    var svg = _atomTraceDrag.svg;
    var rect = svg.getBoundingClientRect();
    var x = e.clientX - rect.left;
    var y = e.clientY - rect.top;
    var dx = x - _atomTraceDrag.startX;
    var dy = y - _atomTraceDrag.startY;
    if (Math.sqrt(dx * dx + dy * dy) > DRAG_THRESHOLD_PX) {
        _atomTraceDrag.moved = true;
    }
    _atomTraceDrag.currentX = x;
    _atomTraceDrag.currentY = y;
    if (_atomTraceDrag.moved) {
        renderAtomTraceMarquee(svg, _atomTraceDrag);
        if (e.preventDefault) { e.preventDefault(); }
    }
}

function handleAtomTraceMouseUp(e) {
    if (!_atomTraceDrag) { return; }
    var drag = _atomTraceDrag;
    _atomTraceDrag = null;
    if (drag.rectEl && drag.rectEl.parentNode) {
        drag.rectEl.parentNode.removeChild(drag.rectEl);
    }
    if (!drag.moved) {
        // Genuine click — let the existing click dispatcher handle it.
        return;
    }
    // Marquee selection: compute box + hit-test atoms.
    var x0 = Math.min(drag.startX, drag.currentX);
    var x1 = Math.max(drag.startX, drag.currentX);
    var y0 = Math.min(drag.startY, drag.currentY);
    var y1 = Math.max(drag.startY, drag.currentY);
    var atoms = drag.svg.querySelectorAll('.bime-trace-atom');
    var selected = [];
    for (var i = 0; i < atoms.length; i++) {
        var hit = atoms[i].querySelector('.bime-trace-atom-hit');
        if (!hit) { continue; }
        var cx = parseFloat(hit.getAttribute('cx'));
        var cy = parseFloat(hit.getAttribute('cy'));
        if (cx >= x0 && cx <= x1 && cy >= y0 && cy <= y1) {
            var idx = parseInt(atoms[i].getAttribute('data-atom-index'), 10);
            if (isFinite(idx)) { selected.push(idx); }
        }
    }
    _suppressNextAtomTraceClick = true;
    if (selected.length === 0) {
        // Empty drag — leave any existing moiety alone.
        return;
    }
    if (drag.shiftKey && _pathway.moietyStartMetIdx === drag.cellIdx) {
        // Extend the existing moiety.
        var existing = _pathway.moietySet || [];
        for (var s = 0; s < selected.length; s++) {
            if (existing.indexOf(selected[s]) === -1) { existing.push(selected[s]); }
        }
        existing.sort(function(a, b) { return a - b; });
        _pathway.moietySet = existing;
    } else {
        // Fresh selection rooted at this cell.
        _pathway.moietyStartMetIdx = drag.cellIdx;
        _pathway.moietySet = selected.slice().sort(function(a, b) { return a - b; });
    }
    showAtomTraceStrip(drag.cellIdx, _pathway.moietySet[0]);
    if (e.preventDefault) { e.preventDefault(); }
}

function renderAtomTraceMarquee(svg, drag) {
    var x0 = Math.min(drag.startX, drag.currentX);
    var x1 = Math.max(drag.startX, drag.currentX);
    var y0 = Math.min(drag.startY, drag.currentY);
    var y1 = Math.max(drag.startY, drag.currentY);
    if (!drag.rectEl) {
        drag.rectEl = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        drag.rectEl.setAttribute('class', 'bime-trace-marquee');
        drag.rectEl.setAttribute('fill', 'rgba(13,148,136,0.10)');
        drag.rectEl.setAttribute('stroke', '#0d9488');
        drag.rectEl.setAttribute('stroke-width', '1.5');
        drag.rectEl.setAttribute('stroke-dasharray', '4 2');
        drag.rectEl.setAttribute('pointer-events', 'none');
        svg.appendChild(drag.rectEl);
    }
    drag.rectEl.setAttribute('x', x0);
    drag.rectEl.setAttribute('y', y0);
    drag.rectEl.setAttribute('width', x1 - x0);
    drag.rectEl.setAttribute('height', y1 - y0);
}

// v2.0.54: explicit "Clear moiety" affordance — public so the click
// dispatcher can call it.
function clearMoietyTrace() {
    var firstAtom = (_pathway.moietySet && _pathway.moietySet[0]) || 0;
    var startIdx = _pathway.moietyStartMetIdx;
    _pathway.moietyStartMetIdx = null;
    _pathway.moietySet = [];
    if (startIdx !== null && startIdx !== undefined) {
        showAtomTraceStrip(startIdx, firstAtom);
    } else {
        showAtomTraceStrip();
    }
}

function pathwaySvgDataUri(svg) {
    svg = String(svg || '');
    if (!svg || /<script/i.test(svg) || /<foreignObject/i.test(svg) || /\son[a-z]+\s*=/i.test(svg)) return '';
    return 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
}

function handlePathwayArrowTarget(nodeId) {
    if (!_pathway.pendingArrow) {
        _pathway.pendingArrow = nodeId;
        selectPathwayItem('node', nodeId);
        renderPathwayCanvas();
        pathwayStatus('Choose a second metabolite node to finish the arrow.');
        return;
    }
    if (_pathway.pendingArrow === nodeId) {
        _pathway.pendingArrow = null;
        pathwayStatus('Arrow cancelled.');
        renderPathwayCanvas();
        return;
    }
    var edgeLabel = pathwayInputLabel('', 48);
    if (edgeLabel === 'Metabolite') edgeLabel = '';
    // v2.0.47: snapshot the pre-commit state so Ctrl+Z removes the
    // newly-committed arrow rather than just clearing the selection.
    if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var edge = {
        id: newPathwayId('edge-'),
        from: _pathway.pendingArrow,
        to: nodeId,
        label: edgeLabel
    };
    _pathway.edges.push(edge);
    _pathway.pendingArrow = null;
    selectPathwayItem('edge', edge.id);
    renderPathwayCanvas();
    pathwayStatus('Arrow added.');
}

function canDragPathwayItem(type) {
    return type === 'background' || type === 'node' || type === 'note' || type === 'step' || type === 'compartment';
}

function canUseMechanismAnchor(type) {
    return type === 'node' || type === 'step' || type === 'note';
}

function mechanismArrowKind() {
    var select = document.getElementById('wb-pathway-arrow-kind');
    return select && select.value === 'single' ? 'single' : 'pair';
}

function handleMechanismArrowTarget(type, id) {
    if (!_pathway.pendingMechanismArrow) {
        _pathway.pendingMechanismArrow = { type: type, id: id };
        selectPathwayItem(type, id);
        renderPathwayCanvas();
        pathwayStatus('Choose a second structure or step to finish the electron-flow arrow.');
        return;
    }
    if (_pathway.pendingMechanismArrow.type === type && _pathway.pendingMechanismArrow.id === id) {
        _pathway.pendingMechanismArrow = null;
        pathwayStatus('Electron-flow arrow cancelled.');
        renderPathwayCanvas();
        return;
    }
    var label = pathwayInputLabel('', 48);
    if (label === 'Metabolite') label = '';
    // v2.0.47: snapshot pre-commit state so Ctrl+Z removes the
    // newly-committed mechanism arrow.
    if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var arrow = {
        id: newPathwayId('mech-'),
        fromType: _pathway.pendingMechanismArrow.type,
        fromId: _pathway.pendingMechanismArrow.id,
        toType: type,
        toId: id,
        kind: mechanismArrowKind(),
        label: label
    };
    _pathway.mechanismArrows.push(arrow);
    _pathway.pendingMechanismArrow = null;
    selectPathwayItem('mechanism', arrow.id);
    renderPathwayCanvas();
    pathwayStatus((arrow.kind === 'single' ? 'Single-electron' : 'Electron-pair') + ' curly arrow added.');
}

function nodeCenter(node) {
    return {
        x: node.x + node.w / 2,
        y: node.y + node.h / 2
    };
}

function pathwayItemCenter(type, item) {
    if (!item) return null;
    return {
        x: item.x + (item.w || 1) / 2,
        y: item.y + (item.h || 1) / 2
    };
}

function edgeEndpoints(edge) {
    var from = findPathwayItem('node', edge.from);
    var to = findPathwayItem('node', edge.to);
    if (!from || !to) return null;
    var a = nodeCenter(from);
    var b = nodeCenter(to);
    var dx = b.x - a.x;
    var dy = b.y - a.y;
    var len = Math.sqrt(dx * dx + dy * dy) || 1;
    var fromPad = Math.min(from.w, from.h) * 0.45;
    var toPad = Math.min(to.w, to.h) * 0.45;
    return {
        x1: a.x + dx / len * fromPad,
        y1: a.y + dy / len * fromPad,
        x2: b.x - dx / len * toPad,
        y2: b.y - dy / len * toPad,
        mx: (a.x + b.x) / 2,
        my: (a.y + b.y) / 2
    };
}

function renderPathwayCanvas() {
    var viewport = document.getElementById('wb-pathway-viewport');
    if (!viewport) return;
    _pathwayRenderFrame = 0;
    viewport.textContent = '';
    // v2.0.75: emit the viewport transform via the shared CanvasSurface helper
    // (byte-identical string); fall back to the inline literal when the bundle
    // predates CanvasSurface.
    if (typeof CanvasSurface !== 'undefined' && CanvasSurface.applyTransform) {
        CanvasSurface.applyTransform(viewport, _pathway);
    } else {
        viewport.setAttribute('transform', 'translate(' + _pathway.offsetX + ' ' + _pathway.offsetY + ') scale(' + _pathway.scale + ')');
    }
    // v2.0.75: rebuild the geometric hit-test registry as we paint. The DOM
    // `closest('[data-pathway-id]')` path in pathwayEventItem stays the primary
    // picker; this registry is the world-coordinate parallel the Stage-2 focus
    // lens needs (pick a node from a world point). Populated in paint order so
    // pick() returns the topmost item.
    rebuildPathwayHitRegistry();
    for (var b = 0; b < _pathway.backgrounds.length; b++) renderPathwayBackground(viewport, _pathway.backgrounds[b]);
    for (var c = 0; c < _pathway.compartments.length; c++) renderPathwayCompartment(viewport, _pathway.compartments[c]);
    for (var e = 0; e < _pathway.edges.length; e++) renderPathwayEdge(viewport, _pathway.edges[e]);
    for (var n = 0; n < _pathway.nodes.length; n++) renderPathwayNode(viewport, _pathway.nodes[n]);
    for (var s = 0; s < _pathway.steps.length; s++) renderPathwayStep(viewport, _pathway.steps[s]);
    for (var t = 0; t < _pathway.notes.length; t++) renderPathwayNote(viewport, _pathway.notes[t]);
    for (var m = 0; m < _pathway.mechanismArrows.length; m++) renderMechanismArrow(viewport, _pathway.mechanismArrows[m]);
    // v2.3.5 (P0-5): a fresh, empty canvas teaches nothing on first run — a
    // blank crosshair grid hides the focus lens (the product's core value).
    // Paint a centered, muted onboarding CTA when NOTHING has been placed yet.
    // It is drawn ONLY in the empty branch, so the very next render after a
    // node/edge/etc. is added drops it. pointer-events:none keeps it from
    // swallowing the double-click that starts drawing (the user can click
    // "through" the hint).
    if (isPathwayEmpty()) { drawEmptyCanvasCta(viewport); }
    // Merge Stage 2: keep the focus-lens overlay tracking the node across any
    // re-render (pan / zoom / wheel / layout). Guarded so a closed lens (or an
    // older bundle without FocusLens) costs nothing. open/collapse re-render
    // with the lens already toggled, so this never double-positions or recurses.
    if (_lens && _lens.isOpen()) { positionStructureLens(); }
}

// v2.3.5 (P0-5): the canvas is "empty" only when EVERY pickable + paintable
// layer is empty — the same all-collections-clear state clearPathwayCanvas()
// resets to. Kept as a tiny predicate so renderPathwayCanvas reads cleanly and
// the onboarding CTA branch can't drift from the real notion of "empty".
function isPathwayEmpty() {
    return _pathway.nodes.length === 0 &&
        _pathway.edges.length === 0 &&
        _pathway.steps.length === 0 &&
        _pathway.notes.length === 0 &&
        _pathway.mechanismArrows.length === 0 &&
        _pathway.compartments.length === 0 &&
        _pathway.backgrounds.length === 0;
}

// v2.3.5 (P0-5): paint the first-run onboarding prompt into the SVG viewport.
// Centered on the active page, muted (--color tokens), terse, on-brand. The
// whole group is pointer-events:none + aria-hidden so it neither swallows the
// double-click that starts drawing nor narrates to screen readers (the canvas
// SVG already carries role/aria-label). No animation — a subtle static chevron
// points down toward the Example controls; honoring prefers-reduced-motion is
// trivial because nothing moves. Drawn only from renderPathwayCanvas's empty
// branch, so any added content removes it on the next render.
function drawEmptyCanvasCta(viewport) {
    var size = pathwayPageSize();
    var cx = size.w / 2;
    var cy = size.h / 2;
    var g = pathwaySvgEl('g', {
        'class': 'wb-pathway-cta',
        'aria-hidden': 'true',
        'pointer-events': 'none'
    });
    g.appendChild(pathwaySvgEl('text', {
        'class': 'wb-pathway-cta-title',
        x: cx,
        y: cy - 18,
        'text-anchor': 'middle',
        fill: 'var(--color-text-muted)',
        text: 'Your canvas is empty'
    }));
    g.appendChild(pathwaySvgEl('text', {
        'class': 'wb-pathway-cta-line',
        x: cx,
        y: cy + 12,
        'text-anchor': 'middle',
        fill: 'var(--color-text-muted)',
        text: 'Double-click anywhere to draw a structure — or load an Example below'
    }));
    // v2.4.3 (Hybrid Stage D): a second, equally terse line surfaces the
    // standalone Editor view as a third way in — draw a molecule or reaction
    // there, then "Add to pathway" drops it onto this canvas. This is how a
    // new user discovers the Editor tab + the import flow; it stays muted
    // (--color tokens), pointer-events:none, aria-hidden like the rest of the CTA.
    g.appendChild(pathwaySvgEl('text', {
        'class': 'wb-pathway-cta-line',
        x: cx,
        y: cy + 34,
        'text-anchor': 'middle',
        fill: 'var(--color-text-muted)',
        text: 'Or open the Editor tab to draw, then Add to pathway'
    }));
    // A subtle downward chevron toward the Example controls beneath the canvas.
    g.appendChild(pathwaySvgEl('path', {
        'class': 'wb-pathway-cta-chevron',
        d: 'M ' + (cx - 12) + ' ' + (cy + 60) + ' L ' + cx + ' ' + (cy + 72) + ' L ' + (cx + 12) + ' ' + (cy + 60),
        fill: 'none',
        stroke: 'var(--color-text-muted)',
        'stroke-width': 2,
        'stroke-linecap': 'round',
        'stroke-linejoin': 'round'
    }));
    viewport.appendChild(g);
    return g;
}

// v2.0.75: module-level world-coordinate hit registry for the pathway canvas.
// Rebuilt on every render from node/compartment/note bounding boxes. The Stage-2
// focus lens picks against this; today it is a tested parallel to the DOM picker.
var _pathwayHitRegistry = (typeof CanvasSurface !== 'undefined' && CanvasSurface.HitRegistry)
    ? CanvasSurface.HitRegistry()
    : null;

function rebuildPathwayHitRegistry() {
    if (!_pathwayHitRegistry) { return null; }
    _pathwayHitRegistry.clear();
    var i;
    // Register the PICKABLE layers in their relative paint order:
    // compartments (back) -> nodes -> steps -> notes (front). This is a SUBSET
    // of the full SVG paint order — backgrounds, edges, and mechanism arrows
    // are intentionally NOT registered (the focus lens never picks them), so
    // topmost-wins is correct for the registered set. Stage 2 decides whether
    // edges need registration before pathwayHitPick becomes the sole picker.
    for (i = 0; i < _pathway.compartments.length; i++) {
        var cm = _pathway.compartments[i];
        _pathwayHitRegistry.add(cm.id, 'compartment', cm.x, cm.y, cm.w, cm.h);
    }
    for (i = 0; i < _pathway.nodes.length; i++) {
        var nd = _pathway.nodes[i];
        _pathwayHitRegistry.add(nd.id, 'node', nd.x, nd.y, nd.w, nd.h);
    }
    for (i = 0; i < _pathway.steps.length; i++) {
        var st = _pathway.steps[i];
        if (typeof st.w === 'number' && typeof st.h === 'number') {
            _pathwayHitRegistry.add(st.id, 'step', st.x, st.y, st.w, st.h);
        }
    }
    for (i = 0; i < _pathway.notes.length; i++) {
        var nt = _pathway.notes[i];
        if (typeof nt.w === 'number' && typeof nt.h === 'number') {
            _pathwayHitRegistry.add(nt.id, 'note', nt.x, nt.y, nt.w, nt.h);
        }
    }
    return _pathwayHitRegistry;
}

// v2.0.75: world-coordinate pick (Stage-2 focus lens). Returns {id,type} or null.
function pathwayHitPick(worldX, worldY) {
    return _pathwayHitRegistry ? _pathwayHitRegistry.pick(worldX, worldY) : null;
}

function schedulePathwayRender() {
    if (_pathwayRenderFrame) return;
    var raf = typeof requestAnimationFrame === 'function'
        ? requestAnimationFrame
        : function(cb) { return setTimeout(cb, 16); };
    _pathwayRenderFrame = raf(renderPathwayCanvas);
}

function pathwaySvgEl(tag, attrs) {
    var el = document.createElementNS(PATHWAY_NS, tag);
    attrs = attrs || {};
    for (var key in attrs) {
        if (!Object.prototype.hasOwnProperty.call(attrs, key)) continue;
        if (key === 'text') {
            el.textContent = attrs[key];
        } else if (key === 'href') {
            el.setAttribute('href', attrs[key]);
            el.setAttributeNS(PATHWAY_XLINK_NS, 'xlink:href', attrs[key]);
        } else {
            el.setAttribute(key, attrs[key]);
        }
    }
    return el;
}

function pathwayGroup(type, id, label) {
    var g = pathwaySvgEl('g', {
        'class': 'wb-pathway-' + type + (selectedPathwayItem(type, id) ? ' is-selected' : ''),
        'data-pathway-type': type,
        'data-pathway-id': id,
        'role': 'button',
        'aria-label': label || type
    });
    return g;
}

function selectedPathwayItem(type, id) {
    return _pathway.selectedType === type && _pathway.selectedId === id;
}

function pendingMechanismAnchor(type, id) {
    return !!(_pathway.pendingMechanismArrow &&
        _pathway.pendingMechanismArrow.type === type &&
        _pathway.pendingMechanismArrow.id === id);
}

function renderPathwayCompartment(viewport, item) {
    var g = pathwayGroup('compartment', item.id, item.label);
    g.appendChild(pathwaySvgEl('rect', {
        x: item.x,
        y: item.y,
        width: item.w,
        height: item.h,
        rx: 8,
        fill: 'rgba(241,245,249,0.62)',
        stroke: selectedPathwayItem('compartment', item.id) ? '#0d9488' : '#94a3b8',
        'stroke-width': selectedPathwayItem('compartment', item.id) ? 3 : 2,
        'stroke-dasharray': '8 8'
    }));
    g.appendChild(pathwaySvgEl('text', {
        x: item.x + 16,
        y: item.y + 28,
        fill: '#475569',
        'font-size': 18,
        'font-weight': 700,
        text: truncatePathwayText(item.label, 34)
    }));
    viewport.appendChild(g);
}

function renderPathwayBackground(viewport, item) {
    var g = pathwayGroup('background', item.id, item.name);
    if (item.kind === 'image') {
        g.appendChild(pathwaySvgEl('image', {
            x: item.x,
            y: item.y,
            width: item.w,
            height: item.h,
            href: item.href,
            preserveAspectRatio: 'xMidYMid meet',
            opacity: 0.72
        }));
        g.appendChild(pathwaySvgEl('rect', {
            x: item.x,
            y: item.y,
            width: item.w,
            height: item.h,
            fill: 'none',
            stroke: selectedPathwayItem('background', item.id) ? '#0d9488' : '#94a3b8',
            'stroke-width': selectedPathwayItem('background', item.id) ? 3 : 1.5,
            'stroke-dasharray': selectedPathwayItem('background', item.id) ? '0' : '8 8'
        }));
        g.appendChild(pathwaySvgEl('text', {
            x: item.x + 12,
            y: item.y + 24,
            fill: '#475569',
            'font-size': 15,
            'font-weight': 700,
            text: truncatePathwayText(item.name, 44)
        }));
    } else {
        g.appendChild(pathwaySvgEl('rect', {
            x: item.x,
            y: item.y,
            width: item.w,
            height: item.h,
            rx: 8,
            fill: 'rgba(255,255,255,0.86)',
            stroke: selectedPathwayItem('background', item.id) ? '#0d9488' : '#94a3b8',
            'stroke-width': selectedPathwayItem('background', item.id) ? 3 : 2
        }));
        g.appendChild(pathwaySvgEl('text', {
            x: item.x + 22,
            y: item.y + 48,
            fill: '#475569',
            'font-size': 15,
            'font-weight': 700,
            text: truncatePathwayText(item.name, 46)
        }));
        g.appendChild(pathwaySvgEl('text', {
            x: item.x + 22,
            y: item.y + 82,
            fill: '#64748b',
            'font-size': 12,
            'font-weight': 650,
            'class': 'wb-pathway-bg-meta',
            text: truncatePathwayText(item.meta, 64)
        }));
        g.appendChild(pathwaySvgEl('text', {
            x: item.x + 22,
            y: item.y + 116,
            fill: '#64748b',
            'font-size': 12,
            'font-weight': 650,
            'class': 'wb-pathway-bg-meta',
            text: 'Use as a reference, or import an exported page image for tracing.'
        }));
    }
    viewport.appendChild(g);
}

// v2.0.18: chemistry-standard reaction-arrow types for pathway edges. Drawn
// explicitly in SVG (font glyphs are unreliable for harpoons). Semantics per
// IUPAC convention:
//   forward    →   single full head (default)
//   reverse    ←   full head at the start
//   reversible ⇌   two offset half-harpoons (dynamic EQUILIBRIUM, not ↔)
//   resonance  ↔   one shaft, full heads at both ends (resonance, NOT equilibrium)
var PATHWAY_ARROW_TYPES = ['forward', 'reverse', 'reversible', 'resonance'];
var PATHWAY_ARROW_LABELS = {
    forward: 'Forward →',
    reverse: 'Reverse ←',
    reversible: 'Reversible ⇌',
    resonance: 'Resonance ↔'
};

function pathwayArrowGeom(x1, y1, x2, y2) {
    var dx = x2 - x1, dy = y2 - y1;
    var len = Math.sqrt(dx * dx + dy * dy) || 1;
    var ux = dx / len, uy = dy / len;
    return { ux: ux, uy: uy, px: -uy, py: ux, len: len };
}

// Draw a reaction arrow of the given type into group `g`, between the two
// trimmed endpoints, using stroke colour/width `color`/`w`.
function drawPathwayEdgeArrow(g, type, pts, geo, color, w) {
    var HEAD = 11, WIDE = 5;
    function line(ax, ay, bx, by, dash) {
        var attrs = {
            d: 'M ' + ax + ' ' + ay + ' L ' + bx + ' ' + by,
            fill: 'none', stroke: color, 'stroke-width': w, 'stroke-linecap': 'round'
        };
        // v2.3.7 (P1-5b): resonance (↔, heads at BOTH ends) otherwise shares an
        // identical solid shaft with forward/reverse, so the only difference from
        // a forward arrow was a second arrowhead. A dashed shaft makes resonance
        // unmistakable in greyscale and reinforces that it is NOT a reaction /
        // equilibrium arrow. forward / reverse / reversible stay solid.
        if (dash) { attrs['stroke-dasharray'] = '7 5'; }
        return pathwaySvgEl('path', attrs);
    }
    function fullHead(tx, ty, ux, uy) {
        var bx = tx - ux * HEAD, by = ty - uy * HEAD;
        return pathwaySvgEl('path', {
            d: 'M ' + tx + ' ' + ty +
               ' L ' + (bx + geo.px * WIDE) + ' ' + (by + geo.py * WIDE) +
               ' L ' + (bx - geo.px * WIDE) + ' ' + (by - geo.py * WIDE) + ' Z',
            fill: color, stroke: 'none'
        });
    }
    // single-barb harpoon: from the tip back to ONE side only, `side` = +1/-1
    function harpoon(tx, ty, ux, uy, side) {
        var bx = tx - ux * (HEAD + 1), by = ty - uy * (HEAD + 1);
        return pathwaySvgEl('path', {
            d: 'M ' + tx + ' ' + ty + ' L ' + (bx + geo.px * (WIDE + 1) * side) + ' ' + (by + geo.py * (WIDE + 1) * side),
            fill: 'none', stroke: color, 'stroke-width': w, 'stroke-linecap': 'round'
        });
    }
    if (type === 'reversible') {
        var off = 4;
        // top shaft left→right with a centre-side harpoon at the right end
        g.appendChild(line(pts.x1 + geo.px * off, pts.y1 + geo.py * off,
                           pts.x2 + geo.px * off, pts.y2 + geo.py * off));
        g.appendChild(harpoon(pts.x2 + geo.px * off, pts.y2 + geo.py * off, geo.ux, geo.uy, -1));
        // bottom shaft right→left with a centre-side harpoon at the left end
        g.appendChild(line(pts.x2 - geo.px * off, pts.y2 - geo.py * off,
                           pts.x1 - geo.px * off, pts.y1 - geo.py * off));
        g.appendChild(harpoon(pts.x1 - geo.px * off, pts.y1 - geo.py * off, -geo.ux, -geo.uy, 1));
        return;
    }
    g.appendChild(line(pts.x1, pts.y1, pts.x2, pts.y2, type === 'resonance'));
    if (type === 'forward' || type === 'resonance') { g.appendChild(fullHead(pts.x2, pts.y2, geo.ux, geo.uy)); }
    if (type === 'reverse' || type === 'resonance') { g.appendChild(fullHead(pts.x1, pts.y1, -geo.ux, -geo.uy)); }
}

function renderPathwayEdge(viewport, edge) {
    var pts = edgeEndpoints(edge);
    if (!pts) return;
    var sel = selectedPathwayItem('edge', edge.id);
    var color = sel ? '#0d9488' : '#334155';
    var w = sel ? 4 : 3;
    var g = pathwayGroup('edge', edge.id, edge.label || 'pathway arrow');
    var type = edge.arrowType;
    if (PATHWAY_ARROW_TYPES.indexOf(type) < 0) { type = 'forward'; }
    drawPathwayEdgeArrow(g, type, pts, pathwayArrowGeom(pts.x1, pts.y1, pts.x2, pts.y2), color, w);
    if (edge.label) {
        g.appendChild(pathwaySvgEl('text', {
            x: pts.mx,
            y: pts.my - 10,
            fill: '#334155',
            'font-size': 15,
            'font-weight': 650,
            'text-anchor': 'middle',
            text: truncatePathwayText(edge.label, 28)
        }));
    }
    viewport.appendChild(g);
}

function mechanismArrowPoints(arrow) {
    var from = findPathwayItem(arrow.fromType, arrow.fromId);
    var to = findPathwayItem(arrow.toType, arrow.toId);
    var a = pathwayItemCenter(arrow.fromType, from);
    var b = pathwayItemCenter(arrow.toType, to);
    if (!a || !b) return null;
    var dx = b.x - a.x;
    var dy = b.y - a.y;
    var len = Math.sqrt(dx * dx + dy * dy) || 1;
    var fromPad = Math.min(from.w || 90, from.h || 60) * 0.38;
    var toPad = Math.min(to.w || 90, to.h || 60) * 0.38;
    var x1 = a.x + dx / len * fromPad;
    var y1 = a.y + dy / len * fromPad;
    var x2 = b.x - dx / len * toPad;
    var y2 = b.y - dy / len * toPad;
    var curve = Math.min(125, Math.max(42, len * 0.22));
    var nx = -dy / len;
    var ny = dx / len;
    return {
        x1: x1,
        y1: y1,
        x2: x2,
        y2: y2,
        cx: (x1 + x2) / 2 + nx * curve,
        cy: (y1 + y2) / 2 + ny * curve
    };
}

function renderMechanismArrow(viewport, arrow) {
    var pts = mechanismArrowPoints(arrow);
    if (!pts) return;
    var g = pathwayGroup('mechanism', arrow.id, arrow.label || 'electron-flow arrow');
    var cls = g.getAttribute('class') || '';
    if (arrow.kind === 'single') cls += ' is-single';
    g.setAttribute('class', cls);
    g.appendChild(pathwaySvgEl('path', {
        d: 'M ' + pts.x1 + ' ' + pts.y1 + ' Q ' + pts.cx + ' ' + pts.cy + ' ' + pts.x2 + ' ' + pts.y2,
        fill: 'none',
        stroke: selectedPathwayItem('mechanism', arrow.id) ? '#0d9488' : '#b45309',
        'stroke-width': selectedPathwayItem('mechanism', arrow.id) ? 4 : 3,
        'stroke-linecap': 'round',
        'stroke-dasharray': arrow.kind === 'single' ? '7 5' : '',
        'marker-end': arrow.kind === 'single' ? 'url(#wb-pathway-curly-head-single)' : 'url(#wb-pathway-curly-head-pair)'
    }));
    if (arrow.label) {
        g.appendChild(pathwaySvgEl('text', {
            x: pts.cx,
            y: pts.cy - 8,
            fill: '#92400e',
            'font-size': 14,
            'font-weight': 700,
            'text-anchor': 'middle',
            text: truncatePathwayText(arrow.label, 28)
        }));
    }
    viewport.appendChild(g);
}

function renderPathwayNode(viewport, node) {
    var g = pathwayGroup('node', node.id, node.label);
    var selected = selectedPathwayItem('node', node.id) ||
        _pathway.pendingArrow === node.id ||
        pendingMechanismAnchor('node', node.id);
    var cls = g.getAttribute('class') || '';
    if (node.kind) cls += ' is-' + sanitizePathwayText(node.kind, 24);
    // v2.3.4 (P0-1): every pathway node opens the focus lens (double-click or
    // Enter — see handlePathwayCanvasDoubleClick / onPathwayKeyDown), so flag
    // the group `is-structure` for the hover/selected editable affordance and
    // give it a native <title> tooltip that doubles as the SVG accessible name
    // advertising the edit gesture. Non-node items (steps/notes/edges) never get
    // this — they have their own renderers and are not lens-openable.
    cls += ' is-structure';
    g.setAttribute('class', cls);
    g.appendChild(pathwaySvgEl('title', { text: 'Double-click or press Enter to edit structure' }));
    var isResidue = node.kind === 'residue';
    var isCofactor = node.kind === 'cofactor';
    var isReaction = node.kind === 'reaction' || (node.structure && node.structure.type === 'reaction');
    var fill = isReaction ? '#eff6ff' : (isResidue ? '#ecfdf5' : (isCofactor ? '#fff7ed' : '#ffffff'));
    var stroke = selected ? '#0d9488' : (isReaction ? '#2563eb' : (isResidue ? '#059669' : (isCofactor ? '#f97316' : '#334155')));
    g.appendChild(pathwaySvgEl('rect', {
        x: node.x,
        y: node.y,
        width: node.w,
        height: node.h,
        rx: isReaction ? 10 : (isResidue || isCofactor ? 18 : 8),
        fill: fill,
        stroke: stroke,
        'stroke-width': selected ? 3 : 2
    }));
    if (node.imageHref) {
        g.appendChild(pathwaySvgEl('image', {
            x: node.x + 10,
            y: node.y + 8,
            width: node.w - 20,
            height: 62,
            href: node.imageHref,
            preserveAspectRatio: 'xMidYMid meet'
        }));
        g.appendChild(pathwaySvgEl('text', {
            x: node.x + node.w / 2,
            y: node.y + node.h - 28,
            fill: '#111827',
            'font-size': 16,
            'font-weight': 700,
            'text-anchor': 'middle',
            'dominant-baseline': 'middle',
            text: truncatePathwayText(node.label, 22)
        }));
        if (node.subtitle) {
            g.appendChild(pathwaySvgEl('text', {
                x: node.x + node.w / 2,
                y: node.y + node.h - 11,
                fill: '#64748b',
                'font-size': 12,
                'font-weight': 600,
                'text-anchor': 'middle',
                'dominant-baseline': 'middle',
                'class': 'wb-pathway-node-subtitle',
                text: truncatePathwayText(node.subtitle, 26)
            }));
        }
    } else {
        g.appendChild(pathwaySvgEl('text', {
            x: node.x + node.w / 2,
            y: node.y + (node.subtitle ? node.h / 2 - 8 : node.h / 2),
            fill: '#111827',
            'font-size': isResidue || isCofactor ? 18 : 16,
            'font-weight': 800,
            'text-anchor': 'middle',
            'dominant-baseline': 'middle',
            text: truncatePathwayText(node.label, 24)
        }));
        if (node.subtitle) {
            g.appendChild(pathwaySvgEl('text', {
                x: node.x + node.w / 2,
                y: node.y + node.h / 2 + 14,
                fill: '#64748b',
                'font-size': 11,
                'font-weight': 650,
                'text-anchor': 'middle',
                'dominant-baseline': 'middle',
                'class': 'wb-pathway-node-subtitle',
                text: truncatePathwayText(node.subtitle, 18)
            }));
        }
    }
    // v2.3.7 (P1-5a): node TYPE is otherwise carried only by fill/stroke colour
    // (metabolite white, residue green, cofactor orange, reaction blue), which a
    // greyscale or colour-blind viewer cannot read. Stamp a small uppercase
    // letter badge in the node's top-left corner so the type is perceivable
    // without colour. Drawn LAST so it sits above the rect + any structure image.
    appendPathwayNodeTypeBadge(g, node);
    viewport.appendChild(g);
}

// v2.3.7 (P1-5a): redundant non-colour cue for pathway node type. Returns the
// short uppercase tag for a node's kind (the same semantic the fill/stroke
// colour encodes), or '' for an untyped/structure-only node where no colour
// distinction is drawn either. Kept tiny + uppercase so it reads as a type
// badge, not a label.
var PATHWAY_NODE_TYPE_BADGE = {
    metabolite: 'M',
    residue: 'R',
    cofactor: 'Co',
    reaction: 'Rxn'
};

function pathwayNodeTypeBadgeText(node) {
    if (!node) { return ''; }
    var kind = node.kind;
    if (!kind && node.structure && node.structure.type === 'reaction') { kind = 'reaction'; }
    if (kind === 'reaction' || (node.structure && node.structure.type === 'reaction')) { return 'Rxn'; }
    if (Object.prototype.hasOwnProperty.call(PATHWAY_NODE_TYPE_BADGE, kind)) {
        return PATHWAY_NODE_TYPE_BADGE[kind];
    }
    // Default structure node (no kind / 'metabolite') — the white node still
    // wants a cue so it is told apart from the coloured ones in greyscale.
    return 'M';
}

function appendPathwayNodeTypeBadge(g, node) {
    var text = pathwayNodeTypeBadgeText(node);
    if (!text) { return; }
    var bw = text.length > 1 ? 22 : 16;
    var badge = pathwaySvgEl('g', { 'class': 'wb-pathway-node-badge' });
    badge.appendChild(pathwaySvgEl('rect', {
        x: node.x + 4,
        y: node.y + 4,
        width: bw,
        height: 15,
        rx: 4
    }));
    badge.appendChild(pathwaySvgEl('text', {
        x: node.x + 4 + bw / 2,
        y: node.y + 4 + 8,
        'text-anchor': 'middle',
        'dominant-baseline': 'middle',
        text: text
    }));
    g.appendChild(badge);
}

function renderPathwayStep(viewport, step) {
    var g = pathwayGroup('step', step.id, step.label);
    var selected = selectedPathwayItem('step', step.id) || pendingMechanismAnchor('step', step.id);
    g.appendChild(pathwaySvgEl('rect', {
        x: step.x,
        y: step.y,
        width: step.w,
        height: step.h,
        rx: 8,
        fill: '#eef2ff',
        stroke: selected ? '#0d9488' : '#4f46e5',
        'stroke-width': selected ? 3 : 2
    }));
    g.appendChild(pathwaySvgEl('text', {
        x: step.x + 14,
        y: step.y + 24,
        fill: '#312e81',
        'font-size': 15,
        'font-weight': 700,
        text: 'Step ' + step.number + ': ' + truncatePathwayText(step.label, 24)
    }));
    g.appendChild(pathwaySvgEl('text', {
        x: step.x + 14,
        y: step.y + 50,
        fill: '#64748b',
        'font-size': 12,
        'font-weight': 650,
        'class': 'wb-pathway-step-meta',
        text: truncatePathwayText(step.comment, 34)
    }));
    g.appendChild(pathwaySvgEl('text', {
        x: step.x + 14,
        y: step.y + 70,
        fill: '#64748b',
        'font-size': 12,
        'font-weight': 650,
        'class': 'wb-pathway-step-meta',
        text: 'curly arrows: ' + mechanismCountForStep(step.id)
    }));
    viewport.appendChild(g);
}

function mechanismCountForStep(stepId) {
    var count = 0;
    for (var i = 0; i < _pathway.mechanismArrows.length; i++) {
        var arrow = _pathway.mechanismArrows[i];
        if ((arrow.fromType === 'step' && arrow.fromId === stepId) ||
                (arrow.toType === 'step' && arrow.toId === stepId)) {
            count += 1;
        }
    }
    return count;
}

function renderPathwayNote(viewport, note) {
    var g = pathwayGroup('note', note.id, note.label);
    var selected = selectedPathwayItem('note', note.id) || pendingMechanismAnchor('note', note.id);
    g.appendChild(pathwaySvgEl('rect', {
        x: note.x,
        y: note.y,
        width: note.w,
        height: note.h,
        rx: 8,
        fill: '#fff7ed',
        stroke: selected ? '#0d9488' : '#f59e0b',
        'stroke-width': selected ? 3 : 2
    }));
    g.appendChild(pathwaySvgEl('text', {
        x: note.x + 12,
        y: note.y + 28,
        fill: '#78350f',
        'font-size': 16,
        'font-weight': 650,
        text: truncatePathwayText(note.label, 30)
    }));
    viewport.appendChild(g);
}

function truncatePathwayText(value, maxLen) {
    value = sanitizePathwayText(value, maxLen + 8);
    if (value.length <= maxLen) return value;
    return value.slice(0, Math.max(0, maxLen - 3)) + '...';
}

function deletePathwaySelection() {
    if (!_pathway.selectedType || !_pathway.selectedId) return;
    // v2.0.47: snapshot before delete so Ctrl+Z restores the item.
    // (The keyboard-driven Delete handler in onPathwayKeyDown already
    // pushes before calling this — that push is a no-op duplicate
    // because PathwayHistory dedupes consecutive identical snapshots.)
    if (typeof pushPathwayHistory === 'function') { pushPathwayHistory(); }
    var type = _pathway.selectedType;
    var id = _pathway.selectedId;
    var collection = pathwayCollection(type);
    for (var i = collection.length - 1; i >= 0; i--) {
        if (collection[i].id === id) collection.splice(i, 1);
    }
    if (type === 'node') {
        for (var e = _pathway.edges.length - 1; e >= 0; e--) {
            if (_pathway.edges[e].from === id || _pathway.edges[e].to === id) _pathway.edges.splice(e, 1);
        }
    }
    if (type === 'node' || type === 'step' || type === 'note') {
        for (var m = _pathway.mechanismArrows.length - 1; m >= 0; m--) {
            if ((_pathway.mechanismArrows[m].fromType === type && _pathway.mechanismArrows[m].fromId === id) ||
                    (_pathway.mechanismArrows[m].toType === type && _pathway.mechanismArrows[m].toId === id)) {
                _pathway.mechanismArrows.splice(m, 1);
            }
        }
    }
    if (_pathway.pendingArrow === id) _pathway.pendingArrow = null;
    if (_pathway.pendingMechanismArrow &&
            _pathway.pendingMechanismArrow.type === type &&
            _pathway.pendingMechanismArrow.id === id) {
        _pathway.pendingMechanismArrow = null;
    }
    selectPathwayItem('', null);
    renderPathwayCanvas();
    pathwayStatus('Selected item deleted.');
}

function nudgePathwaySelection(key, amount) {
    if (!_pathway.selectedType || !_pathway.selectedId) return;
    var item = findPathwayItem(_pathway.selectedType, _pathway.selectedId);
    if (!item || item.x === undefined || item.y === undefined) return;
    if (key === 'ArrowLeft') item.x -= amount;
    if (key === 'ArrowRight') item.x += amount;
    if (key === 'ArrowUp') item.y -= amount;
    if (key === 'ArrowDown') item.y += amount;
    renderPathwayCanvas();
}

function zoomPathway(dir) {
    // v2.0.42: route through CanvasView.zoom so the min/max clamping
    // and factor semantics live in one place. Offset is unaffected
    // (no pivot — CanvasView leaves it alone), matching pre-v2.0.42.
    var factor = (dir === 'out') ? 0.85 : 1.18;
    if (typeof CanvasView !== 'undefined' && CanvasView.zoom) {
        var zoomed = CanvasView.zoom(_pathway, factor, { minScale: 0.35, maxScale: 3.4 });
        _pathway.scale = zoomed.scale;
        _pathway.offsetX = zoomed.offsetX;
        _pathway.offsetY = zoomed.offsetY;
    } else {
        var next = _pathway.scale * factor;
        _pathway.scale = Math.max(0.35, Math.min(3.4, next));
    }
    renderPathwayCanvas();
    pathwayStatus('Pathway zoom ' + Math.round(_pathway.scale * 100) + '%.');
}

function setPathwayLayoutShape(shape) {
    _pathway.layoutShape = normalizePathwayLayoutShape(shape);
    layoutPathwayCanvas();
}

function pathwayBounds() {
    var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    function addRect(item) {
        minX = Math.min(minX, item.x);
        minY = Math.min(minY, item.y);
        maxX = Math.max(maxX, item.x + item.w);
        maxY = Math.max(maxY, item.y + item.h);
    }
    for (var i = 0; i < _pathway.compartments.length; i++) addRect(_pathway.compartments[i]);
    for (var j = 0; j < _pathway.nodes.length; j++) addRect(_pathway.nodes[j]);
    for (var k = 0; k < _pathway.notes.length; k++) addRect(_pathway.notes[k]);
    for (var s = 0; s < _pathway.steps.length; s++) addRect(_pathway.steps[s]);
    for (var b = 0; b < _pathway.backgrounds.length; b++) addRect(_pathway.backgrounds[b]);
    if (!isFinite(minX)) return null;
    return { x: minX, y: minY, w: maxX - minX, h: maxY - minY };
}

// v2.0.28: Clean Up. When the pathway has any edges, runs the Sugiyama-style
// hierarchical layout in editor/PathwayLayout.js (connected components →
// longest-path layering → barycentric crossing reduction). Falls back to the
// pre-v2.0.28 linear arrangement when there are no edges (or PathwayLayout
// is unavailable, e.g. running the workbench against the unbundled editor
// during development without re-running the build).
function layoutPathwayCanvas() {
    var i;
    var nodeCount = _pathway.nodes.length;
    var stepCount = _pathway.steps.length;
    var edgeCount = _pathway.edges ? _pathway.edges.length : 0;
    var usedSugiyama = false;
    var lastStats = null;

    if (edgeCount > 0 && typeof PathwayLayout !== 'undefined' && PathwayLayout.compute) {
        // Edge-driven hierarchical layout. Pass {id} pairs so the algorithm
        // doesn't see any extra state, then offset by node half-width/height
        // (the algorithm returns NODE CENTRES; _pathway.nodes store top-left).
        var algoNodes = [];
        for (i = 0; i < _pathway.nodes.length; i++) { algoNodes.push({ id: _pathway.nodes[i].id }); }
        var algoEdges = [];
        for (i = 0; i < _pathway.edges.length; i++) {
            algoEdges.push({ from: _pathway.edges[i].from, to: _pathway.edges[i].to });
        }
        var requestedShape = normalizePathwayLayoutShape(_pathway.layoutShape || 'auto');
        var laid = PathwayLayout.compute(algoNodes, algoEdges, {
            orientation: 'horizontal',
            layoutShape: requestedShape
        });
        if (laid && laid.positions) {
            for (i = 0; i < _pathway.nodes.length; i++) {
                var nd = _pathway.nodes[i];
                var p = laid.positions[nd.id];
                if (p) {
                    nd.x = Math.round(p.x - nd.w / 2);
                    nd.y = Math.round(p.y - nd.h / 2);
                }
            }
            usedSugiyama = true;
            lastStats = laid.stats;
            if (requestedShape !== 'auto' && lastStats && lastStats.layoutShape) {
                _pathway.layoutShape = lastStats.layoutShape;
            }
        }
    }

    // Mechanism step boxes — for the moment still arranged in a single header
    // row above the graph. (Curly-arrow mechanism support is queued for v2.0.29.)
    var maxCount = Math.max(nodeCount, stepCount, 1);
    var gap = Math.max(190, Math.min(300, 940 / Math.max(1, maxCount - 1 || 1)));
    var startX = 130;
    if (maxCount === 1) startX = 510;
    for (i = 0; i < _pathway.steps.length; i++) {
        _pathway.steps[i].number = i + 1;
        _pathway.steps[i].x = Math.round(startX + i * gap - _pathway.steps[i].w / 2);
        _pathway.steps[i].y = 74;
    }
    // Fall-through linear arrangement for nodes when there are no edges.
    if (!usedSugiyama) {
        var nodeY = stepCount ? 245 : 160;
        for (i = 0; i < _pathway.nodes.length; i++) {
            _pathway.nodes[i].x = Math.round(startX + i * gap - _pathway.nodes[i].w / 2);
            _pathway.nodes[i].y = nodeY;
        }
    }
    for (i = 0; i < _pathway.notes.length; i++) {
        _pathway.notes[i].x = Math.round(120 + (i % 4) * 250);
        _pathway.notes[i].y = 425 + Math.floor(i / 4) * 92;
    }
    var contentBounds = pathwayContentBounds();
    for (i = 0; i < _pathway.compartments.length; i++) {
        if (contentBounds) {
            _pathway.compartments[i].x = Math.max(20, Math.round(contentBounds.x - 54 + i * 18));
            _pathway.compartments[i].y = Math.max(20, Math.round(contentBounds.y - 42 + i * 18));
            _pathway.compartments[i].w = Math.min(1120, Math.round(contentBounds.w + 108));
            _pathway.compartments[i].h = Math.min(540, Math.round(contentBounds.h + 84));
        }
    }
    fitPathwayCanvas();
    if (usedSugiyama && lastStats) {
        var shapeLabel = lastStats.layoutShape || 'hierarchical';
        pathwayStatus('Clean Up: ' + nodeCount + ' nodes, ' + shapeLabel + ' layout, ' + lastStats.crossings + ' crossings.');
    } else {
        pathwayStatus('Pathway and mechanism steps arranged.');
    }
}

function pathwayContentBounds() {
    var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    function addRect(item) {
        minX = Math.min(minX, item.x);
        minY = Math.min(minY, item.y);
        maxX = Math.max(maxX, item.x + item.w);
        maxY = Math.max(maxY, item.y + item.h);
    }
    for (var i = 0; i < _pathway.nodes.length; i++) addRect(_pathway.nodes[i]);
    for (var j = 0; j < _pathway.notes.length; j++) addRect(_pathway.notes[j]);
    for (var k = 0; k < _pathway.steps.length; k++) addRect(_pathway.steps[k]);
    if (!isFinite(minX)) return null;
    return { x: minX, y: minY, w: maxX - minX, h: maxY - minY };
}

function fitPathwayCanvas() {
    var bounds = pathwayBounds();
    if (!bounds) {
        _pathway.scale = 1;
        _pathway.offsetX = 0;
        _pathway.offsetY = 0;
        renderPathwayCanvas();
        pathwayStatus('Pathway canvas reset.');
        return;
    }
    // v2.0.42: routed through CanvasView.fit so the centre-and-pad math
    // is shared with future Renderer.fit calls. Page-format aware: uses
    // the active pageFormat / pageOrientation viewBox as the viewport.
    // The pre-v2.0.42 implementation used a fixed 54-pixel pad against
    // the 1200×620 viewBox; CanvasView's padding is a fraction (0.045
    // ≈ 54/1200), so the resulting scale is byte-equivalent.
    var size = pathwayPageSize();
    if (typeof CanvasView !== 'undefined' && CanvasView.fit) {
        var v = CanvasView.fit(bounds,
            { w: size.w, h: size.h },
            { padding: 54 / size.w });
        // Clamp the resulting scale to the legacy [0.35, 2.2] range.
        var s = Math.max(0.35, Math.min(2.2, v.scale));
        if (s !== v.scale) {
            // Re-derive offset to keep content centred under the clamp.
            v = CanvasView.fit(bounds, { w: size.w, h: size.h }, { padding: 54 / size.w });
            _pathway.scale = s;
            _pathway.offsetX = size.w / 2 - (bounds.x + bounds.w / 2) * s;
            _pathway.offsetY = size.h / 2 - (bounds.y + bounds.h / 2) * s;
        } else {
            _pathway.scale = v.scale;
            _pathway.offsetX = v.offsetX;
            _pathway.offsetY = v.offsetY;
        }
    } else {
        var pad = 54;
        var scaleX = (size.w - pad * 2) / Math.max(1, bounds.w);
        var scaleY = (size.h - pad * 2) / Math.max(1, bounds.h);
        _pathway.scale = Math.max(0.35, Math.min(2.2, Math.min(scaleX, scaleY)));
        _pathway.offsetX = size.w / 2 - (bounds.x + bounds.w / 2) * _pathway.scale;
        _pathway.offsetY = size.h / 2 - (bounds.y + bounds.h / 2) * _pathway.scale;
    }
    renderPathwayCanvas();
    pathwayStatus('Pathway fit to canvas.');
}

function exportPathwaySvg() {
    var svg = document.getElementById('wb-pathway-canvas');
    if (!svg) return;
    // v2.0.38: use the active page-format size for the exported SVG.
    var size = pathwayPageSize();
    var clone = svg.cloneNode(true);
    clone.setAttribute('xmlns', PATHWAY_NS);
    clone.setAttribute('width', String(size.w));
    clone.setAttribute('height', String(size.h));
    clone.setAttribute('viewBox', '0 0 ' + size.w + ' ' + size.h);
    var serialized = new XMLSerializer().serializeToString(clone);
    if (typeof ImageExport !== 'undefined' && ImageExport.downloadSVG) {
        ImageExport.downloadSVG(serialized, 'bime-pathway-canvas.svg');
    } else {
        downloadTextFile(serialized, 'bime-pathway-canvas.svg', 'image/svg+xml;charset=utf-8');
    }
    pathwayStatus('Pathway SVG exported (' + size.w + ' × ' + size.h + ').');
}

// v2.0.38: page-format setter. Updates the live SVG viewBox + width/height
// so the on-screen canvas matches print dimensions, then re-renders.
function setPathwayPageFormat(format, orientation) {
    if (typeof PageFormats === 'undefined') {
        pathwayStatus('PageFormats module not loaded.');
        return;
    }
    if (!PageFormats.isKnown(format, orientation)) {
        pathwayStatus('Unknown page format: ' + format + ' / ' + orientation);
        return;
    }
    _pathway.pageFormat = format;
    _pathway.pageOrientation = orientation;
    var svg = document.getElementById('wb-pathway-canvas');
    if (svg) {
        var vb = PageFormats.viewBoxFor(format, orientation);
        svg.setAttribute('viewBox', vb);
    }
    if (typeof renderPathwayCanvas === 'function') { renderPathwayCanvas(); }
    if (typeof fitPathwayCanvas === 'function') { fitPathwayCanvas(); }
    var pageLabel = (PageFormats.FORMATS[format] && PageFormats.FORMATS[format].label) || format;
    pathwayStatus('Page format: ' + pageLabel + ' ' + orientation + '.');
}

// Helper: current page size (w, h). Falls back to legacy 1200×620 if
// PageFormats isn't loaded.
function pathwayPageSize() {
    if (typeof PageFormats !== 'undefined' && PageFormats.sizeFor) {
        return PageFormats.sizeFor(_pathway.pageFormat || 'screen', _pathway.pageOrientation || 'landscape');
    }
    return { w: 1200, h: 620 };
}

// v2.0.32: SBGN-ML export (PD subset)
function exportPathwaySbgnml() {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.serializeSBGNML) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    var xml = PathwayIO.serializeSBGNML(_pathway);
    downloadTextFile(xml, 'bime-pathway.sbgn', 'application/xml;charset=utf-8');
    pathwayStatus('SBGN-ML exported (' + _pathway.nodes.length + ' nodes, ' + _pathway.edges.length + ' edges).');
}

// v2.0.32: BioPAX Level 3 export (OWL/XML)
function exportPathwayBiopax() {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.serializeBioPAX) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    var owl = PathwayIO.serializeBioPAX(_pathway);
    downloadTextFile(owl, 'bime-pathway.owl', 'application/rdf+xml;charset=utf-8');
    pathwayStatus('BioPAX Level 3 exported (' + _pathway.nodes.length + ' entities, ' + _pathway.edges.length + ' reactions).');
}

// v2.0.30: pathway round-trip — native BIME JSON
function exportPathwayJson() {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.serializeJSON) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    var payload = PathwayIO.serializeJSON(_pathway);
    var json = JSON.stringify(payload, null, 2);
    downloadTextFile(json, 'bime-pathway.json', 'application/json;charset=utf-8');
    pathwayStatus('Pathway JSON exported (' + payload.nodes.length + ' nodes, ' + payload.edges.length + ' edges).');
}

// v2.0.30: pathway round-trip — restore from native BIME JSON
function openPathwayJsonPicker() {
    var input = document.getElementById('wb-pathway-json-input');
    if (!input) {
        input = document.createElement('input');
        input.type = 'file';
        input.id = 'wb-pathway-json-input';
        input.accept = 'application/json,.json';
        input.style.display = 'none';
        input.setAttribute('aria-hidden', 'true');
        input.setAttribute('tabindex', '-1');
        input.addEventListener('change', function(e) {
            var file = e.target && e.target.files ? e.target.files[0] : null;
            if (file) { readPathwayTextFile(file, applyImportedJsonText); }
            try { input.value = ''; } catch (err) {}
        });
        document.body.appendChild(input);
    }
    input.click();
}

// v2.0.30: KEGG KGML import — XML subset (entry + reaction)
function openPathwayKgmlPicker() {
    var input = document.getElementById('wb-pathway-kgml-input');
    if (!input) {
        input = document.createElement('input');
        input.type = 'file';
        input.id = 'wb-pathway-kgml-input';
        input.accept = 'application/xml,text/xml,.xml,.kgml';
        input.style.display = 'none';
        input.setAttribute('aria-hidden', 'true');
        input.setAttribute('tabindex', '-1');
        input.addEventListener('change', function(e) {
            var file = e.target && e.target.files ? e.target.files[0] : null;
            if (file) { readPathwayTextFile(file, applyImportedKgmlText); }
            try { input.value = ''; } catch (err) {}
        });
        document.body.appendChild(input);
    }
    input.click();
}

// v2.0.37: SBGN-ML import (PD subset)
function openPathwaySbgnmlPicker() {
    var input = document.getElementById('wb-pathway-sbgnml-input');
    if (!input) {
        input = document.createElement('input');
        input.type = 'file';
        input.id = 'wb-pathway-sbgnml-input';
        input.accept = 'application/xml,text/xml,.xml,.sbgn';
        input.style.display = 'none';
        input.setAttribute('aria-hidden', 'true');
        input.setAttribute('tabindex', '-1');
        input.addEventListener('change', function(e) {
            var file = e.target && e.target.files ? e.target.files[0] : null;
            if (file) { readPathwayTextFile(file, applyImportedSbgnmlText); }
            try { input.value = ''; } catch (err) {}
        });
        document.body.appendChild(input);
    }
    input.click();
}

// v2.0.37: BioPAX Level 3 import
function openPathwayBiopaxPicker() {
    var input = document.getElementById('wb-pathway-biopax-input');
    if (!input) {
        input = document.createElement('input');
        input.type = 'file';
        input.id = 'wb-pathway-biopax-input';
        input.accept = 'application/rdf+xml,application/xml,.owl,.xml,.biopax';
        input.style.display = 'none';
        input.setAttribute('aria-hidden', 'true');
        input.setAttribute('tabindex', '-1');
        input.addEventListener('change', function(e) {
            var file = e.target && e.target.files ? e.target.files[0] : null;
            if (file) { readPathwayTextFile(file, applyImportedBiopaxText); }
            try { input.value = ''; } catch (err) {}
        });
        document.body.appendChild(input);
    }
    input.click();
}

function applyImportedSbgnmlText(text, fileName) {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.parseSBGNML) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    try {
        var loaded = PathwayIO.parseSBGNML(text);
        applyPathwayStateFromImport(loaded, fileName || 'SBGN-ML');
        pathwayStatus('SBGN-ML imported (' + loaded.nodes.length + ' nodes, ' +
                      loaded.edges.length + ' edges).');
    } catch (err) {
        pathwayStatus('SBGN-ML import failed: ' + (err && err.message ? err.message : 'malformed'));
    }
}

function applyImportedBiopaxText(text, fileName) {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.parseBioPAX) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    try {
        var loaded = PathwayIO.parseBioPAX(text);
        applyPathwayStateFromImport(loaded, fileName || 'BioPAX');
        pathwayStatus('BioPAX imported (' + loaded.nodes.length + ' entities, ' +
                      loaded.edges.length + ' reactions). Run Clean Up to lay out.');
    } catch (err) {
        pathwayStatus('BioPAX import failed: ' + (err && err.message ? err.message : 'malformed'));
    }
}

function readPathwayTextFile(file, then) {
    if (file.size > 4 * 1024 * 1024) {
        pathwayStatus('File too large (>4 MB).');
        return;
    }
    var reader = new FileReader();
    reader.onload = function() { then(String(reader.result || ''), file.name); };
    reader.onerror = function() { pathwayStatus('Read failed: ' + (reader.error && reader.error.message || 'unknown')); };
    reader.readAsText(file);
}

function applyImportedJsonText(text, fileName) {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.parseJSON) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    try {
        var loaded = PathwayIO.parseJSON(text);
        applyPathwayStateFromImport(loaded, fileName || 'JSON');
        pathwayStatus('Imported ' + (fileName || 'pathway.json') + ' (' +
                      loaded.nodes.length + ' nodes, ' + loaded.edges.length + ' edges).');
    } catch (err) {
        pathwayStatus('JSON import failed: ' + (err && err.message ? err.message : 'malformed'));
    }
}

function applyImportedKgmlText(text, fileName) {
    if (typeof PathwayIO === 'undefined' || !PathwayIO.parseKGML) {
        pathwayStatus('PathwayIO module not loaded.');
        return;
    }
    try {
        var loaded = PathwayIO.parseKGML(text);
        applyPathwayStateFromImport(loaded, fileName || 'KGML');
        pathwayStatus('KGML imported (' + loaded.nodes.length + ' compound nodes, ' +
                      loaded.edges.length + ' reaction edges).');
    } catch (err) {
        pathwayStatus('KGML import failed: ' + (err && err.message ? err.message : 'malformed'));
    }
}

function applyPathwayStateFromImport(loaded, label) {
    setWorkbenchMode('pathway');
    clearPathwayCanvas();
    if (loaded.scale)   { _pathway.scale = loaded.scale; }
    if (loaded.offsetX !== undefined) { _pathway.offsetX = loaded.offsetX; }
    if (loaded.offsetY !== undefined) { _pathway.offsetY = loaded.offsetY; }
    // v2.0.44: restore page format from the JSON envelope. Falls back to
    // screen / landscape if the imported file is a v2.0.30..43 export
    // that didn't carry the new fields.
    if (typeof loaded.pageFormat === 'string' && typeof loaded.pageOrientation === 'string') {
        if (typeof setPathwayPageFormat === 'function') {
            setPathwayPageFormat(loaded.pageFormat, loaded.pageOrientation);
        }
    }
    _pathway.layoutShape = normalizePathwayLayoutShape(loaded.layoutShape || 'auto');
    _pathway.nodes = loaded.nodes.slice();
    _pathway.edges = loaded.edges.slice();
    _pathway.notes = (loaded.notes || []).slice();
    _pathway.steps = (loaded.steps || []).slice();
    _pathway.compartments = (loaded.compartments || []).slice();
    _pathway.mechanismArrows = (loaded.mechanismArrows || []).slice();
    // Reset next-id counter past anything we just imported so freshly added
    // shapes can't collide with imported ids.
    var maxN = 0;
    function bump(items) {
        for (var i = 0; i < items.length; i++) {
            var m = /(\d+)$/.exec(items[i].id || '');
            if (m) { var v = parseInt(m[1], 10); if (v > maxN) { maxN = v; } }
        }
    }
    bump(_pathway.nodes); bump(_pathway.edges); bump(_pathway.notes);
    bump(_pathway.steps); bump(_pathway.compartments);
    _pathway.nextId = maxN + 1;
    renderPathwayCanvas();
    if (typeof fitPathwayCanvas === 'function') { fitPathwayCanvas(); }
}

function downloadTextFile(text, filename, type) {
    var blob = new Blob([text], { type: type || 'text/plain;charset=utf-8' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(url); }, 0);
}

function clearPathwayCanvas() {
    _pathway.nodes = [];
    _pathway.edges = [];
    _pathway.notes = [];
    _pathway.steps = [];
    _pathway.mechanismArrows = [];
    _pathway.compartments = [];
    _pathway.backgrounds = [];
    _pathway.nextId = 1;
    _pathway.selectedType = '';
    _pathway.selectedId = null;
    _pathway.pendingArrow = null;
    _pathway.pendingMechanismArrow = null;
    _pathway.drag = null;
    _pathway.layoutShape = 'auto';
    renderPathwayCanvas();
    pathwayStatus('Pathway canvas cleared.');
    // v2.0.22: also tear down the atom-trace strip — it is bound to whatever
    // example was loaded, and that example no longer exists.
    if (typeof clearAtomTraceStrip === 'function') { clearAtomTraceStrip(); }
    // v2.0.43: clearing the canvas is a fresh-start signal — drop the
    // undo stack so Ctrl+Z doesn't time-travel back into the previous
    // example. The very next edit pushes the new empty-state as the
    // baseline snapshot.
    if (typeof clearPathwayHistory === 'function') { clearPathwayHistory(); }
}

function getSelectedAtoms() {
    if (!editor || !editor.molecule) return [];
    return editor.molecule.atoms.filter(function(a) { return !!a.selected; });
}

function getSelectedBonds() {
    if (!editor || !editor.molecule) return [];
    return editor.molecule.bonds.filter(function(b) { return !!b.selected; });
}

function getMoleculeInsights() {
    var mol = editor && editor.molecule;
    var empty = !mol || !mol.atoms || mol.atoms.length === 0;
    var key = moleculePerfSignature(mol);
    if (_insightsCacheValue && _insightsCacheKey === key) {
        return _insightsCacheValue;
    }
    var q = evaluateLayoutQuality(mol);
    var warnings = empty ? [] : buildStructureWarnings(mol, q);
    var counts = empty ? {} : formulaCountsForMol(mol);
    var rings = 0;
    if (!empty && mol.findRings) {
        try { rings = (mol.findRings(20) || []).length; } catch (e) { rings = 0; }
    }
    _insightsCacheKey = key;
    _insightsCacheValue = {
        atoms: empty ? 0 : mol.atoms.length,
        bonds: empty ? 0 : mol.bonds.length,
        rings: rings,
        charge: empty ? 0 : netCharge(mol),
        formula: empty ? 'No structure' : formulaFromCounts(counts),
        molWeight: empty ? 0 : molecularWeightFromCounts(counts),
        molWeightLabel: empty ? '—' : formatMass(molecularWeightFromCounts(counts)),
        hbd: empty ? 0 : donorCountForMol(mol),
        hba: empty ? 0 : acceptorCountForMol(mol),
        rotatableBonds: empty ? 0 : rotatableBondCountForMol(mol),
        warnings: warnings,
        quality: q,
        qualityLabel: qualityLabel(q, warnings, empty)
    };
    return _insightsCacheValue;
}

function moleculePerfSignature(mol) {
    if (!mol || !mol.atoms) return 'empty';
    var parts = [mol.atoms.length, mol.bonds ? mol.bonds.length : 0, mol.name || '', mol.comment || ''];
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        parts.push(a.id, a.symbol, a.charge || 0, a.isotope || 0,
            a.hydrogens === undefined ? '' : a.hydrogens,
            a.aromatic ? 1 : 0,
            Math.round((a.x || 0) * 100),
            Math.round((a.y || 0) * 100));
    }
    if (mol.bonds) {
        for (var j = 0; j < mol.bonds.length; j++) {
            var b = mol.bonds[j];
            parts.push(b.id, b.atom1, b.atom2, b.type, b.stereo || 0);
        }
    }
    if (mol.reactionArrow) {
        parts.push('rxn',
            Math.round((mol.reactionArrow.x1 || 0) * 100),
            Math.round((mol.reactionArrow.y1 || 0) * 100),
            Math.round((mol.reactionArrow.x2 || 0) * 100),
            Math.round((mol.reactionArrow.y2 || 0) * 100),
            mol.reactionArrow.conditions || '');
    }
    return parts.join('|');
}

function moleculeChemSignature(mol) {
    if (!mol || !mol.atoms) return 'empty';
    var parts = [mol.atoms.length, mol.bonds ? mol.bonds.length : 0, mol.name || ''];
    for (var i = 0; i < mol.atoms.length; i++) {
        var a = mol.atoms[i];
        parts.push(a.id, a.symbol, a.charge || 0, a.isotope || 0,
            a.hydrogens === undefined ? '' : a.hydrogens,
            a.aromatic ? 1 : 0, a.mapNumber || 0);
    }
    if (mol.bonds) {
        for (var j = 0; j < mol.bonds.length; j++) {
            var b = mol.bonds[j];
            parts.push(b.id, b.atom1, b.atom2, b.type, b.stereo || 0);
        }
    }
    if (mol.reactionArrow) {
        parts.push('rxn', mol.reactionArrow.conditions || '');
    }
    return parts.join('|');
}

function evaluateLayoutQuality(mol) {
    if (!mol || !mol.atoms || mol.atoms.length === 0) return null;
    try {
        if (typeof SDG !== 'undefined' && SDG.LayoutQuality && SDG.LayoutQuality.evaluate) {
            return SDG.LayoutQuality.evaluate(mol);
        }
    } catch (e) {}
    return null;
}

function qualityLabel(q, warnings, empty) {
    if (empty) return 'No structure';
    if (!q) return warnings.length ? 'Review' : 'Ready';
    if (q.hardFailures > 0 || q.severeClashes > 0 || q.duplicateAtoms > 0) return 'Needs cleanup';
    if (q.crossings > 0 || q.bondWarnings > 1 || warnings.length > 0 || q.penalty > 12) return 'Review';
    return 'Clean';
}

function buildStructureWarnings(mol, q) {
    var warnings = [];
    if (!mol) return warnings;
    var ringAtoms = ringAtomIdSet(mol);
    if (q) {
        if (q.nonFiniteAtoms) warnings.push(q.nonFiniteAtoms + ' atom coordinate issue' + plural(q.nonFiniteAtoms) + '.');
        if (q.invalidBonds) warnings.push(q.invalidBonds + ' invalid bond reference' + plural(q.invalidBonds) + '.');
        if (q.duplicateAtoms) warnings.push(q.duplicateAtoms + ' overlapping atom pair' + plural(q.duplicateAtoms) + '.');
        if (q.severeClashes) warnings.push(q.severeClashes + ' severe atom clash' + plural(q.severeClashes) + '.');
        if (q.closePairs) warnings.push(q.closePairs + ' close nonbonded atom pair' + plural(q.closePairs) + '.');
        if (q.shortBonds) warnings.push(q.shortBonds + ' very short bond' + plural(q.shortBonds) + '.');
        if (q.longBonds) warnings.push(q.longBonds + ' stretched bond' + plural(q.longBonds) + '.');
        if (q.crossings) warnings.push(q.crossings + ' bond crossing' + plural(q.crossings) + '.');
        if (q.acuteAngles) warnings.push(q.acuteAngles + ' very acute bond angle' + plural(q.acuteAngles) + '.');
    }
    for (var i = 0; i < mol.bonds.length; i++) {
        var b = mol.bonds[i];
        if (!mol.getAtom(b.atom1) || !mol.getAtom(b.atom2)) {
            warnings.push('Bond ' + b.id + ' points to a missing atom.');
        }
        if (b.type < 1 || b.type > 3) {
            warnings.push('Bond ' + b.id + ' has unsupported order ' + b.type + '.');
        }
    }
    for (var a = 0; a < mol.atoms.length; a++) {
        var atom = mol.atoms[a];
        var bondSum = mol.bondOrderSum ? mol.bondOrderSum(atom.id) : 0;
        var charge = atom.charge || 0;
        var hardLimit = atom.symbol === 'C' ? 4 :
            atom.symbol === 'N' ? 5 :
            atom.symbol === 'O' ? 3 :
            atom.symbol === 'H' ? 1 : 0;
        if (hardLimit && bondSum + Math.max(0, -charge) > hardLimit) {
            if (!isFusedAromaticLikeValence(mol, atom, bondSum, ringAtoms)) {
                warnings.push(atom.symbol + atom.id + ' may exceed normal valence.');
            }
        }
        if (charge > 4 || charge < -4) {
            warnings.push(atom.symbol + atom.id + ' has an unusually large charge.');
        }
    }
    return warnings.slice(0, 12);
}

function plural(n) { return n === 1 ? '' : 's'; }

function isFusedAromaticLikeValence(mol, atom, bondSum, ringAtoms) {
    if (!atom || !ringAtoms || !ringAtoms[atom.id]) return false;
    if (atom.symbol !== 'C' && atom.symbol !== 'N') return false;
    var bonds = mol.getBondsOfAtom ? mol.getBondsOfAtom(atom.id) : [];
    if (bonds.length > 3 || bondSum > 5) return false;
    for (var i = 0; i < bonds.length; i++) {
        var other = mol.getAtom(bonds[i].otherAtom(atom.id));
        if (!other || !ringAtoms[other.id]) return false;
    }
    return true;
}

function formulaCountsForMol(mol) {
    var counts = {};
    if (!mol || !mol.atoms) return counts;
    for (var i = 0; i < mol.atoms.length; i++) {
        var atom = mol.atoms[i];
        counts[atom.symbol] = (counts[atom.symbol] || 0) + 1;
        var h = 0;
        try { h = mol.calcHydrogens ? mol.calcHydrogens(atom.id) : 0; } catch (e) { h = 0; }
        if (h > 0) counts.H = (counts.H || 0) + h;
    }
    return counts;
}

function formulaFromCounts(counts) {
    var keys = Object.keys(counts).sort();
    if (counts.C) keys = ['C'].concat(counts.H ? ['H'] : [], keys.filter(function(k) { return k !== 'C' && k !== 'H'; }));
    return keys.map(function(k) { return k + (counts[k] === 1 ? '' : counts[k]); }).join('');
}

var BIME_ATOMIC_WEIGHTS = {
    H: 1.008, C: 12.011, N: 14.007, O: 15.999, F: 18.998, P: 30.974,
    S: 32.06, Cl: 35.45, Br: 79.904, I: 126.904, B: 10.81, Si: 28.085,
    Se: 78.971, Na: 22.990, K: 39.098, Ca: 40.078, Mg: 24.305,
    Fe: 55.845, Zn: 65.38, Cu: 63.546, Li: 6.94, Al: 26.982, Sn: 118.71
};

function molecularWeightFromCounts(counts) {
    var total = 0;
    Object.keys(counts || {}).forEach(function(symbol) {
        var weight = BIME_ATOMIC_WEIGHTS[symbol] ||
            (typeof Molecule !== 'undefined' && Molecule.ELEMENTS && Molecule.ELEMENTS[symbol] && Molecule.ELEMENTS[symbol].mass) || 0;
        total += weight * counts[symbol];
    });
    return total;
}

function formatMass(value) {
    return value > 0 ? (Math.round(value * 1000) / 1000).toFixed(3) : '—';
}

function formatProperty(value) {
    if (value === undefined || value === null || value === '') return '—';
    var n = Number(value);
    if (isFinite(n)) return String(Math.round(n * 1000) / 1000);
    return String(value);
}

function donorCountForMol(mol) {
    var count = 0;
    for (var i = 0; i < mol.atoms.length; i++) {
        var atom = mol.atoms[i];
        if (!/^(N|O|S)$/.test(atom.symbol)) continue;
        var h = 0;
        try { h = mol.calcHydrogens ? mol.calcHydrogens(atom.id) : 0; } catch (e) { h = 0; }
        if (h > 0) count++;
    }
    return count;
}

function acceptorCountForMol(mol) {
    var count = 0;
    for (var i = 0; i < mol.atoms.length; i++) {
        var atom = mol.atoms[i];
        if (!/^(N|O|S|P)$/.test(atom.symbol)) continue;
        if ((atom.charge || 0) > 0) continue;
        var bondSum = 0;
        try { bondSum = mol.bondOrderSum ? mol.bondOrderSum(atom.id) : 0; } catch (e) { bondSum = 0; }
        if (atom.symbol === 'N' && bondSum >= 4) continue;
        count++;
    }
    return count;
}

function rotatableBondCountForMol(mol) {
    if (!mol || !mol.bonds) return 0;
    var ringBonds = ringBondIdSet(mol);
    var total = 0;
    for (var i = 0; i < mol.bonds.length; i++) {
        var bond = mol.bonds[i];
        if (bond.type !== 1 || bond.stereo) continue;
        if (ringBonds[bond.id]) continue;
        var a1 = mol.getAtom(bond.atom1);
        var a2 = mol.getAtom(bond.atom2);
        if (!a1 || !a2 || a1.symbol === 'H' || a2.symbol === 'H') continue;
        if (heavyDegree(mol, a1.id) > 1 && heavyDegree(mol, a2.id) > 1) total++;
    }
    return total;
}

function heavyDegree(mol, atomId) {
    var bonds = mol.getBondsOfAtom ? mol.getBondsOfAtom(atomId) : [];
    var n = 0;
    for (var i = 0; i < bonds.length; i++) {
        var other = mol.getAtom(bonds[i].otherAtom(atomId));
        if (other && other.symbol !== 'H') n++;
    }
    return n;
}

function ringBondIdSet(mol) {
    var set = {};
    if (!mol || !mol.findRings) return set;
    var rings = [];
    try { rings = mol.findRings(20) || []; } catch (e) { return set; }
    rings.forEach(function(ring) {
        var atoms = ring.atoms || ring;
        for (var i = 0; i < atoms.length; i++) {
            var a = atoms[i];
            var b = atoms[(i + 1) % atoms.length];
            var bond = mol.getBondBetween ? mol.getBondBetween(a, b) : null;
            if (bond) set[bond.id] = true;
        }
    });
    return set;
}

function ringAtomIdSet(mol) {
    var set = {};
    if (!mol || !mol.findRings) return set;
    var rings = [];
    try { rings = mol.findRings(20) || []; } catch (e) { return set; }
    rings.forEach(function(ring) {
        var atoms = ring.atoms || ring;
        for (var i = 0; i < atoms.length; i++) set[atoms[i]] = true;
    });
    return set;
}

function netCharge(mol) {
    var sum = 0;
    for (var i = 0; i < mol.atoms.length; i++) sum += mol.atoms[i].charge || 0;
    return sum;
}

function updateInspector() {
    var insights = getMoleculeInsights();
    setText('wb-atoms-chip', insights.atoms + ' atom' + plural(insights.atoms));
    setText('wb-bonds-chip', insights.bonds + ' bond' + plural(insights.bonds));
    setText('wb-rings-chip', insights.rings + ' ring' + plural(insights.rings));
    setText('wb-warning-count', String(insights.warnings.length));
    setText('wb-quality-detail', insights.quality ? ('Penalty ' + Math.round(insights.quality.penalty * 10) / 10) : 'Ready');

    var badge = document.getElementById('wb-quality-badge');
    if (badge) {
        badge.className = 'wb-quality-badge';
        var cls = insights.qualityLabel === 'Clean' ? 'is-good' :
            insights.qualityLabel === 'Needs cleanup' ? 'is-bad' :
            insights.qualityLabel === 'No structure' ? '' : 'is-warn';
        if (cls) badge.classList.add(cls);
        badge.textContent = insights.qualityLabel;
    }
    var warnPanel = document.getElementById('wb-warning-panel');
    if (insights.warnings.length > 0 && _lastWarningCount === 0) openWorkbenchPanel('warnings');
    _lastWarningCount = insights.warnings.length;
    if (warnPanel && isWorkbenchPanelExpanded('warnings')) renderWarningList(warnPanel, insights.warnings);
    populateStructurePropertyControls();
    populateReactionLabelControls();
    populateReactionPropertyControls();
    renderSelectionInspector();
}

function setText(id, text) {
    var el = document.getElementById(id);
    if (el) el.textContent = text;
}

function renderWarningList(target, warnings) {
    target.textContent = '';
    if (!warnings || warnings.length === 0) {
        var p = document.createElement('p');
        p.className = 'wb-output-muted';
        p.style.margin = '0';
        p.textContent = 'No warnings yet.';
        target.appendChild(p);
        return;
    }
    var ul = document.createElement('ul');
    ul.className = 'wb-warn-list';
    warnings.forEach(function(w) {
        var li = document.createElement('li');
        li.textContent = w;
        ul.appendChild(li);
    });
    target.appendChild(ul);
}

function sanitizeReactionLabelPart(value) {
    return String(value || '').replace(/[\r\n\t]+/g, ' ').replace(/\s+/g, ' ').trim().slice(0, 80);
}

function composeReactionLabel(parts) {
    parts = parts || {};
    var items = [];
    var conditions = sanitizeReactionLabelPart(parts.conditions);
    var yieldText = sanitizeReactionLabelPart(parts.yield);
    var note = sanitizeReactionLabelPart(parts.note);
    if (conditions) items.push(conditions);
    if (yieldText) items.push('yield ' + yieldText);
    if (note) items.push(note);
    return items.join('; ');
}

function setReactionLabelParts(parts) {
    var arrow = editor && editor.molecule && editor.molecule.reactionArrow;
    if (!arrow) return false;
    var clean = {
        conditions: sanitizeReactionLabelPart(parts.conditions),
        yield: sanitizeReactionLabelPart(parts.yield),
        note: sanitizeReactionLabelPart(parts.note)
    };
    arrow.labelParts = clean;
    arrow.conditions = composeReactionLabel(clean);
    return true;
}

function reactionLabelPartsForArrow(arrow) {
    if (!arrow) return { conditions: '', yield: '', note: '' };
    if (arrow.labelParts) {
        return {
            conditions: sanitizeReactionLabelPart(arrow.labelParts.conditions),
            yield: sanitizeReactionLabelPart(arrow.labelParts.yield),
            note: sanitizeReactionLabelPart(arrow.labelParts.note)
        };
    }
    return { conditions: sanitizeReactionLabelPart(arrow.conditions), yield: '', note: '' };
}

function reactionLabelForArrow(arrow) {
    return arrow ? (arrow.conditions || composeReactionLabel(reactionLabelPartsForArrow(arrow))) : '';
}

function populateReactionLabelControls() {
    var arrow = editor && editor.molecule && editor.molecule.reactionArrow;
    var status = document.getElementById('wb-reaction-label-status');
    var help = document.getElementById('wb-reaction-label-help');
    var ids = [
        'wb-reaction-conditions',
        'wb-reaction-yield',
        'wb-reaction-note',
        'wb-prop-reaction-conditions',
        'wb-prop-reaction-yield',
        'wb-prop-reaction-note'
    ];
    var active = document.activeElement;
    var editing = ids.indexOf(active && active.id) >= 0;
    if (status) status.textContent = arrow ? (reactionLabelForArrow(arrow) ? 'Labeled' : 'Ready') : 'No arrow';
    if (help && !editing) help.textContent = arrow ?
        'Labels render above the reaction arrow and export in SVG/PNG.' :
        'Load or draw a reaction arrow before applying a label.';
    if (editing) return;
    var parts = reactionLabelPartsForArrow(arrow);
    setInputValue('wb-reaction-conditions', parts.conditions);
    setInputValue('wb-reaction-yield', parts.yield);
    setInputValue('wb-reaction-note', parts.note);
    setInputValue('wb-prop-reaction-conditions', parts.conditions);
    setInputValue('wb-prop-reaction-yield', parts.yield);
    setInputValue('wb-prop-reaction-note', parts.note);
}

function setInputValue(id, value) {
    var el = document.getElementById(id);
    if (el && el.value !== value) el.value = value;
}

function applyReactionLabels() {
    applyReactionLabelsFromIds({
        conditions: 'wb-reaction-conditions',
        yield: 'wb-reaction-yield',
        note: 'wb-reaction-note'
    }, 'wb-reaction-label-help');
}

function applyPropertyReactionLabels() {
    applyReactionLabelsFromIds({
        conditions: 'wb-prop-reaction-conditions',
        yield: 'wb-prop-reaction-yield',
        note: 'wb-prop-reaction-note'
    }, 'wb-reaction-property-status');
}

function applyReactionLabelsFromIds(ids, statusId) {
    if (!editor || !editor.molecule || !editor.molecule.reactionArrow) {
        setText(statusId || 'wb-reaction-label-help', 'Load or draw a reaction arrow first.');
        return;
    }
    if (editor.saveHistory) editor.saveHistory();
    setReactionLabelParts({
        conditions: document.getElementById(ids.conditions) && document.getElementById(ids.conditions).value,
        yield: document.getElementById(ids.yield) && document.getElementById(ids.yield).value,
        note: document.getElementById(ids.note) && document.getElementById(ids.note).value
    });
    editor.render();
    invalidateWorkbenchOutputCaches();
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
    var msg = reactionLabelForArrow(editor.molecule.reactionArrow) ? 'Reaction label applied.' : 'Reaction label cleared.';
    setText(statusId || 'wb-reaction-label-help', msg);
    if (statusId !== 'wb-reaction-label-help') setText('wb-reaction-label-help', msg);
    if (statusId !== 'wb-reaction-property-status') setText('wb-reaction-property-status', msg);
}

function clearReactionLabels() {
    setInputValue('wb-reaction-conditions', '');
    setInputValue('wb-reaction-yield', '');
    setInputValue('wb-reaction-note', '');
    applyReactionLabels();
}

function clearPropertyReactionLabels() {
    setInputValue('wb-prop-reaction-conditions', '');
    setInputValue('wb-prop-reaction-yield', '');
    setInputValue('wb-prop-reaction-note', '');
    applyPropertyReactionLabels();
}

function populateReactionPropertyControls() {
    var arrow = editor && editor.molecule && editor.molecule.reactionArrow;
    var section = document.getElementById('wb-reaction-property-editor');
    if (section) section.hidden = !arrow;
    var status = document.getElementById('wb-reaction-property-status');
    if (status) status.textContent = arrow ? (reactionLabelForArrow(arrow) ? 'Labeled' : 'Ready') : 'No arrow';
    var active = document.activeElement;
    var editing = active && active.id && active.id.indexOf('wb-prop-reaction-') === 0;
    if (editing) return;
    var parts = reactionLabelPartsForArrow(arrow);
    setInputValue('wb-prop-reaction-conditions', parts.conditions);
    setInputValue('wb-prop-reaction-yield', parts.yield);
    setInputValue('wb-prop-reaction-note', parts.note);
}

function renderSelectionInspector() {
    var panel = document.getElementById('wb-selection-panel');
    if (!panel) return;
    var atoms = getSelectedAtoms();
    var bonds = getSelectedBonds();
    setText('wb-selection-count', atoms.length + bonds.length ? (atoms.length + 'A / ' + bonds.length + 'B') : 'None');
    panel.textContent = '';
    if (!atoms.length && !bonds.length) {
        var p = document.createElement('p');
        p.className = 'wb-output-muted';
        p.style.margin = '0';
        p.textContent = 'Select an atom or bond to edit chemistry details.';
        panel.appendChild(p);
        return;
    }
    if (atoms.length) panel.appendChild(atomInspector(atoms[0], atoms.length));
    if (bonds.length) panel.appendChild(bondInspector(bonds[0], bonds.length));
}

function atomInspector(atom, count) {
    var wrap = document.createElement('div');
    wrap.className = 'wb-mini-grid';
    wrap.appendChild(fieldSelect('Atom', 'wb-inspector-atom-symbol',
        ['C','N','O','S','P','F','Cl','Br','I','H','B','Si','Se','R'], atom.symbol));
    wrap.appendChild(fieldNumber('Charge', 'wb-inspector-atom-charge', atom.charge || 0, -4, 4));
    wrap.appendChild(fieldNumber('Isotope', 'wb-inspector-atom-isotope', atom.isotope || 0, 0, 300));
    wrap.appendChild(fieldNumber('Map', 'wb-inspector-atom-map', atom.mapNumber || 0, 0, 999));
    var btn = document.createElement('button');
    btn.className = 'wb-btn wb-btn-solid';
    btn.type = 'button';
    btn.style.gridColumn = '1 / -1';
    btn.textContent = count > 1 ? 'Apply to first selected atom' : 'Apply atom changes';
    btn.addEventListener('click', function() { applyInspectorAtom(atom.id); });
    wrap.appendChild(btn);
    return wrap;
}

function bondInspector(bond, count) {
    var wrap = document.createElement('div');
    wrap.className = 'wb-mini-grid';
    wrap.style.marginTop = '0.55rem';
    wrap.appendChild(fieldSelect('Order', 'wb-inspector-bond-order', ['1','2','3'], String(bond.type || 1)));
    wrap.appendChild(fieldSelect('Stereo', 'wb-inspector-bond-stereo', ['0','1','6'], String(bond.stereo || 0)));
    var btn = document.createElement('button');
    btn.className = 'wb-btn wb-btn-solid';
    btn.type = 'button';
    btn.style.gridColumn = '1 / -1';
    btn.textContent = count > 1 ? 'Apply to first selected bond' : 'Apply bond changes';
    btn.addEventListener('click', function() { applyInspectorBond(bond.id); });
    wrap.appendChild(btn);
    return wrap;
}

function fieldSelect(labelText, id, values, selected) {
    var wrap = fieldWrap(labelText, id);
    var select = document.createElement('select');
    select.id = id;
    select.className = 'wb-select';
    values.forEach(function(value) {
        var opt = document.createElement('option');
        opt.value = value;
        opt.textContent = value === '0' ? 'None' : value === '1' && id.indexOf('stereo') >= 0 ? 'Wedge' : value === '6' ? 'Dash' : value;
        if (value === selected) opt.selected = true;
        select.appendChild(opt);
    });
    wrap.appendChild(select);
    return wrap;
}

function fieldNumber(labelText, id, value, min, max) {
    var wrap = fieldWrap(labelText, id);
    var input = document.createElement('input');
    input.id = id;
    input.className = 'wb-input';
    input.type = 'number';
    input.min = String(min);
    input.max = String(max);
    input.value = String(value);
    wrap.appendChild(input);
    return wrap;
}

function fieldWrap(labelText, id) {
    var wrap = document.createElement('div');
    wrap.className = 'wb-field';
    var label = document.createElement('label');
    if (id) label.htmlFor = id;
    label.textContent = labelText;
    wrap.appendChild(label);
    return wrap;
}

function applyInspectorAtom(atomId) {
    if (!editor || !editor.molecule) return;
    var atom = editor.molecule.getAtom(atomId);
    if (!atom) return;
    editor.saveHistory();
    atom.symbol = valueOf('wb-inspector-atom-symbol', atom.symbol);
    atom.charge = intValue('wb-inspector-atom-charge', atom.charge || 0);
    atom.isotope = intValue('wb-inspector-atom-isotope', atom.isotope || 0);
    atom.mapNumber = intValue('wb-inspector-atom-map', atom.mapNumber || 0);
    editor.changed();
    scheduleWorkbenchUpdate();
}

function applyInspectorBond(bondId) {
    if (!editor || !editor.molecule) return;
    var bond = editor.molecule.getBond(bondId);
    if (!bond) return;
    editor.saveHistory();
    bond.type = Math.max(1, Math.min(3, intValue('wb-inspector-bond-order', bond.type || 1)));
    bond.stereo = intValue('wb-inspector-bond-stereo', bond.stereo || 0);
    editor.changed();
    scheduleWorkbenchUpdate();
}

function valueOf(id, fallback) {
    var el = document.getElementById(id);
    return el ? el.value : fallback;
}

function intValue(id, fallback) {
    var n = parseInt(valueOf(id, fallback), 10);
    return isFinite(n) ? n : fallback;
}

function clean2DWorkbench() {
    if (!editor) return;
    if (editor._autoLayout) editor._autoLayout();
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
}

function bestFitWorkbench() {
    if (!editor) return;
    if (editor._zoomToFit) editor._zoomToFit();
    else if (editor.renderer && editor.renderer.centerMolecule) { editor.renderer.centerMolecule(); editor.render(); }
    scheduleWorkbenchUpdate();
}

function toggleRSWorkbench() {
    if (editor && editor._doAction) editor._doAction('ciprs');
    scheduleWorkbenchUpdate();
}

function toggleEZWorkbench() {
    if (editor && editor._doAction) editor._doAction('cipez');
    scheduleWorkbenchUpdate();
}

function initTemplateDrawer() {
    var search = document.getElementById('wb-template-search');
    var category = document.getElementById('wb-template-category');
    var grid = document.getElementById('wb-template-grid');
    if (search) search.addEventListener('input', scheduleTemplateDrawerRender);
    if (category) category.addEventListener('change', renderTemplateDrawer);
    if (grid && grid.getAttribute('data-wb-template-bound') !== 'true') {
        grid.setAttribute('data-wb-template-bound', 'true');
        grid.addEventListener('click', function(ev) {
            var card = ev.target && ev.target.closest ? ev.target.closest('.wb-template-card') : null;
            if (card && card.getAttribute('data-template')) applyTemplateFromDrawer(card.getAttribute('data-template'));
        });
        grid.addEventListener('dragstart', function(ev) {
            var card = ev.target && ev.target.closest ? ev.target.closest('.wb-template-card') : null;
            if (!card || !card.getAttribute('data-template') || !ev.dataTransfer) return;
            ev.dataTransfer.setData('application/x-bime-template', card.getAttribute('data-template'));
            ev.dataTransfer.effectAllowed = 'copy';
        });
    }
    if (typeof Templates !== 'undefined' && Templates.list) {
        try { setText('wb-template-count', String(Templates.list().length)); } catch (e) {}
    }
    if (isWorkbenchPanelExpanded('templates')) renderTemplateDrawer();
}

function scheduleTemplateDrawerRender() {
    if (_templateDrawerTimer) clearTimeout(_templateDrawerTimer);
    _templateDrawerTimer = setTimeout(function() {
        _templateDrawerTimer = 0;
        renderTemplateDrawer();
    }, 80);
}

function renderTemplateDrawer() {
    var grid = document.getElementById('wb-template-grid');
    if (!grid || typeof Templates === 'undefined' || !Templates.list) return;
    _templateThumbQueue = [];
    _templateDrawerSeq++;
    var termEl = document.getElementById('wb-template-search');
    var catEl = document.getElementById('wb-template-category');
    var term = (termEl && termEl.value || '').trim().toLowerCase();
    var cat = (catEl && catEl.value) || 'all';
    var names = Templates.list().sort();
    if (cat === 'recent') {
        names = _templateRecent.filter(function(n) { return names.indexOf(n) >= 0; });
    } else if (cat !== 'all') {
        names = names.filter(function(n) { return templateCategory(n) === cat; });
    }
    if (term) {
        names = names.filter(function(n) { return templateSearchText(n).indexOf(term) >= 0; });
    }
    setText('wb-template-count', String(names.length));
    grid.textContent = '';
    if (!names.length) {
        var p = document.createElement('p');
        p.className = 'wb-output-muted';
        p.textContent = 'No templates match.';
        grid.appendChild(p);
        return;
    }
    var frag = document.createDocumentFragment();
    names.forEach(function(name) { frag.appendChild(templateCard(name, _templateDrawerSeq)); });
    grid.appendChild(frag);
}

function templateCard(name, seq) {
    var btn = document.createElement('button');
    btn.className = 'wb-template-card';
    btn.type = 'button';
    btn.draggable = true;
    btn.setAttribute('data-template', name);
    btn.setAttribute('aria-label', 'Add template ' + prettyTemplateName(name));
    var thumb = document.createElement('div');
    thumb.className = 'wb-template-thumb';
    queueTemplateThumb(thumb, name, seq);
    var label = document.createElement('div');
    label.className = 'wb-template-name';
    var text = document.createElement('span');
    text.textContent = prettyTemplateName(name);
    var meta = document.createElement('span');
    meta.className = 'wb-template-meta';
    meta.textContent = templateCategory(name);
    label.appendChild(text);
    label.appendChild(meta);
    btn.appendChild(thumb);
    btn.appendChild(label);
    return btn;
}

function queueTemplateThumb(container, name, seq) {
    if (!container) return;
    if (_templateSvgCache[name] !== undefined) {
        setSafeSvg(container, _templateSvgCache[name], 'Preview');
        return;
    }
    container.textContent = 'Preview';
    _templateThumbQueue.push({ container: container, name: name, seq: seq });
    pumpTemplateThumbQueue();
}

function pumpTemplateThumbQueue() {
    if (_templateThumbFrame) return;
    var schedule = typeof requestIdleCallback === 'function'
        ? function(cb) { return requestIdleCallback(cb, { timeout: 160 }); }
        : (typeof requestAnimationFrame === 'function'
            ? function(cb) { return requestAnimationFrame(function() { cb({ timeRemaining: function() { return 8; } }); }); }
            : function(cb) { return setTimeout(function() { cb({ timeRemaining: function() { return 8; } }); }, 16); });
    _templateThumbFrame = schedule(function(deadline) {
        _templateThumbFrame = 0;
        var processed = 0;
        while (_templateThumbQueue.length > 0) {
            if (deadline && deadline.timeRemaining && deadline.timeRemaining() < 4 && processed > 0) break;
            if (processed >= 3) break;
            var job = _templateThumbQueue.shift();
            processed++;
            if (!job || job.seq !== _templateDrawerSeq || !job.container ||
                    !document.documentElement.contains(job.container)) continue;
            setSafeSvg(job.container, renderTemplateThumb(job.name), 'Preview');
        }
        if (_templateThumbQueue.length > 0) pumpTemplateThumbQueue();
    });
}

function renderTemplateThumb(name) {
    if (_templateSvgCache[name] !== undefined) return _templateSvgCache[name];
    if (typeof Molecule === 'undefined' || typeof Templates === 'undefined' ||
        typeof ImageExport === 'undefined' || !ImageExport.toSVG) {
        _templateSvgCache[name] = null;
        return null;
    }
    try {
        var mol = new Molecule();
        Templates.apply(mol, name, 0, 0);
        var svg = ImageExport.toSVG(mol, {
            width: 144,
            height: 88,
            padding: 8,
            background: 'transparent',
            showAtomLabels: true
        });
        _templateSvgCache[name] = svg || null;
    } catch (e) {
        _templateSvgCache[name] = null;
    }
    return _templateSvgCache[name];
}

function applyTemplateFromDrawer(name) {
    if (!editor || !editor.molecule) return;
    var mol = editor.molecule;
    var x = 0, y = 0;
    if (!mol.isEmpty()) {
        var b = mol.getBounds();
        x = b.x + b.w + 75;
        y = b.y + b.h / 2;
    }
    stampTemplate(name, x, y);
}

function stampTemplate(name, x, y) {
    if (!editor || !editor.molecule || typeof Templates === 'undefined' || !Templates.apply) return;
    if (typeof Templates[name] !== 'function') {
        if (editor.showInfo) editor.showInfo('Template unavailable: ' + prettyTemplateName(name));
        return;
    }
    var mol = editor.molecule;
    editor.saveHistory();
    var ids = [];
    try {
        ids = Templates.apply(mol, name, x || 0, y || 0);
    } catch (e) {
        if (editor.showInfo) editor.showInfo('Could not add ' + prettyTemplateName(name));
        return;
    }
    if (!ids || !ids.length) return;
    mol.atoms.forEach(function(a) { a.selected = ids.indexOf(a.id) >= 0; });
    mol.bonds.forEach(function(b) { b.selected = false; });
    rememberTemplate(name);
    editor.render();
    if (mol.atoms.length === ids.length && editor.renderer && editor.renderer.centerMolecule) {
        editor.renderer.centerMolecule();
        editor.render();
    }
    if (editor.showInfo) editor.showInfo('Added ' + prettyTemplateName(name));
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
}

function rememberTemplate(name) {
    _templateRecent = [name].concat(_templateRecent.filter(function(n) { return n !== name; })).slice(0, 12);
    try { localStorage.setItem('bime-template-recent', JSON.stringify(_templateRecent)); } catch (e) {}
}

function prettyTemplateName(name) {
    return String(name || '').replace(/_/g, ' ').replace(/([a-z])([A-Z])/g, '$1 $2');
}

function templateAliasesFor(name) {
    var aliases = [];
    var map = (typeof Templates !== 'undefined' && Templates.aliases) ? Templates.aliases : null;
    if (!map) return aliases;
    for (var alias in map) {
        if (Object.prototype.hasOwnProperty.call(map, alias) && map[alias] === name) {
            aliases.push(alias);
        }
    }
    return aliases;
}

function templateSearchText(name) {
    var parts = [String(name || ''), prettyTemplateName(name)];
    templateAliasesFor(name).forEach(function(alias) {
        parts.push(alias);
        parts.push(prettyTemplateName(alias));
    });
    return parts.join(' ').toLowerCase();
}

function templateCategory(name) {
    var n = String(name || '').toLowerCase();
    if (/^(benzene|cyclopropane|cyclobutane|cyclopentane|cyclohexane|cycloheptane)$/.test(n)) return 'ring';
    if (/amino|purine|sugar|hexose|pentose|xanthine|pyranose|furanose|peptide/.test(n)) return 'bio';
    if (/naph|anthra|pyrene|phenanth|fluorene|fused|indole|quin|benz|chrom|carbaz|pteridine|acrid/.test(n)) return 'fused';
    if (/pyr|imid|azole|azine|thioph|furan|morph|pip|lactam|tetraz|triaz|oxaz|thiaz/.test(n)) return 'hetero';
    if (/steroid|penam|cepham|carbapenem|tetracycline|adamantane|tropane|aporphine|quinuclidine|morphinan/.test(n)) return 'scaffold';
    if (/cycl|ring/.test(n)) return 'ring';
    return 'scaffold';
}

function initCommandPalette() {
    _commandItems = [
        { label: 'Clean 2D layout', keys: 'clean layout sdg', run: clean2DWorkbench },
        { label: 'Best fit structure', keys: 'fit zoom center', run: bestFitWorkbench },
        { label: 'Compound lookup', keys: 'compound cas name formula cid inchi inchikey', run: focusCompoundLookup },
        { label: 'Open pathway canvas', keys: 'pathway canvas metabolic network annotate map', run: focusPathwayCanvas },
        { label: 'Layout pathway canvas', keys: 'pathway mechanism steps auto layout', run: function() { openWorkbenchSection('wb-pathway-section'); layoutPathwayCanvas(); } },
        { label: 'Convert reaction to pathway step', keys: 'reaction pathway step transform convert', run: convertCurrentReactionToPathway },
        { label: 'Reaction labels', keys: 'conditions yield solvent catalyst arrow note', run: focusReactionLabels },
        { label: 'Add benzene', keys: 'template ring benzene', run: function() { applyTemplateFromDrawer('benzene'); } },
        { label: 'Add cyclohexane', keys: 'template ring cyclohexane', run: function() { applyTemplateFromDrawer('cyclohexane'); } },
        { label: 'Add naphthalene', keys: 'template fused naphthalene', run: function() { applyTemplateFromDrawer('naphthalene'); } },
        { label: 'Show R/S labels', keys: 'stereo chirality rs', run: toggleRSWorkbench },
        { label: 'Show E/Z labels', keys: 'stereo double ez', run: toggleEZWorkbench },
        { label: 'Export publication SVG', keys: 'export svg publication print journal figure', run: function() { runExportPreset('pub-svg'); } },
        { label: 'Export transparent SVG', keys: 'export svg publication transparent journal figure', run: function() { runExportPreset('pub-svg-transparent'); } },
        { label: 'Export transparent PNG', keys: 'export png image transparent publication', run: function() { runExportPreset('png'); } },
        { label: 'Copy current output', keys: 'copy clipboard smiles mol', run: copyOut },
        { label: 'Search library', keys: 'search library mcs substructure', run: runLibrarySearch },
        { label: 'Warnings tab', keys: 'warnings validate valence', run: function() { setOutputTab('warnings'); doValidate(); } }
    ];
    // v2.4.18: drop the pathway command-palette entries in the core build (the
    // pathway editor is hidden behind PATHWAY_ENABLED). The entries stay in the
    // source literal above so the pathway release re-enables them by flipping it.
    if (!PATHWAY_ENABLED) {
        _commandItems = _commandItems.filter(function(it) {
            return !/pathway/i.test(it.keys || '') && !/pathway/i.test(it.label || '');
        });
    }
    var input = document.getElementById('wb-command-input');
    if (input) {
        input.addEventListener('focus', renderCommandResults);
        input.addEventListener('input', renderCommandResults);
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') {
                e.preventDefault();
                var first = document.querySelector('#wb-command-results .wb-command-item');
                if (first) first.click();
            } else if (e.key === 'Escape') {
                input.value = '';
                renderCommandResults();
            }
        });
    }
    document.addEventListener('keydown', function(e) {
        var key = String(e.key || '').toLowerCase();
        if ((e.ctrlKey || e.metaKey) && key === 'k') {
            e.preventDefault();
            focusCommandPalette();
        }
    });
    if (isWorkbenchPanelExpanded('command')) renderCommandResults();
}

function initCompoundLookup() {
    var input = document.getElementById('bime-lookup-input');
    var btn = document.getElementById('bime-lookup-run');
    if (btn) btn.addEventListener('click', runCompoundLookup);
    if (input) {
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') {
                e.preventDefault();
                runCompoundLookup();
            }
        });
    }
}

function initPropertyControls() {
    ['wb-structure-name', 'wb-structure-comment'].forEach(function(id) {
        var input = document.getElementById(id);
        if (!input) return;
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter' && (id !== 'wb-structure-comment' || e.ctrlKey || e.metaKey)) {
                e.preventDefault();
                applyStructureProperties();
            }
        });
    });
    var show = document.getElementById('wb-structure-show-name');
    if (show) {
        show.addEventListener('change', applyStructureProperties);
    }
}

function initReactionLabelControls() {
    [
        'wb-reaction-conditions',
        'wb-reaction-yield',
        'wb-reaction-note',
        'wb-prop-reaction-conditions',
        'wb-prop-reaction-yield',
        'wb-prop-reaction-note'
    ].forEach(function(id) {
        var input = document.getElementById(id);
        if (!input) return;
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') {
                e.preventDefault();
                if (id.indexOf('wb-prop-reaction-') === 0) applyPropertyReactionLabels();
                else applyReactionLabels();
            }
        });
    });
}

function focusCompoundLookup() {
    openWorkbenchSection('bime-compound-lookup');
    var input = document.getElementById('bime-lookup-input');
    if (!input) return;
    input.focus();
    input.select();
    input.scrollIntoView({ block: 'center', behavior: 'smooth' });
}

function focusReactionLabels() {
    openWorkbenchPanel('reaction');
    var input = document.getElementById('wb-reaction-conditions');
    if (!input) return;
    input.focus();
    input.select();
    input.scrollIntoView({ block: 'center', behavior: 'smooth' });
}

function focusCommandPalette() {
    openWorkbenchPanel('command');
    var input = document.getElementById('wb-command-input');
    if (input) {
        input.focus();
        input.select();
        input.scrollIntoView({ block: 'center', behavior: 'smooth' });
    }
}

function renderCommandResults() {
    var list = document.getElementById('wb-command-results');
    if (!list) return;
    var input = document.getElementById('wb-command-input');
    var term = (input && input.value || '').trim().toLowerCase();
    var hits = _commandItems.filter(function(item) {
        return !term || item.label.toLowerCase().indexOf(term) >= 0 || item.keys.indexOf(term) >= 0;
    }).slice(0, 6);
    list.textContent = '';
    hits.forEach(function(item) {
        var btn = document.createElement('button');
        btn.className = 'wb-command-item';
        btn.type = 'button';
        btn.textContent = item.label;
        btn.addEventListener('click', function() {
            item.run();
            if (input) input.value = '';
            renderCommandResults();
        });
        list.appendChild(btn);
    });
}

// v2.4.7: shared base options for the publication export presets -- pulls the
// current editor size / title / H-display, then flags the journal-standard
// publication preset (Helvetica labels, Kekule double bonds, CMYK-safe colours,
// heavier bonds, no watermark). Callers add only the background (white vs
// transparent) and, for PNG, the resolution scale.
function publicationExportOptions(extra) {
    var base = (editor && editor._exportOptions) ? editor._exportOptions() : {};
    var opts = {
        width: base.width,
        height: base.height,
        title: base.title,
        showHydrogens: base.showHydrogens,
        publication: true
    };
    if (extra) {
        for (var k in extra) { if (extra.hasOwnProperty(k)) opts[k] = extra[k]; }
    }
    return opts;
}

// v3.0.3: browser SVG exports carry the SAME BIME provenance as the CLI —
// bime:version metadata, Apache-2.0 attribution, SHA fingerprint, and the
// bottom-right wordmark/glyph. Applied HERE at the workbench layer via
// ExportStamp.stampSvg, deliberately NOT inside ImageExport, so the CLI output,
// the regression tests, and the on-screen preview stay unstamped. stampSvg is
// idempotent, so re-exporting the same figure never double-stamps.
function stampExportSvg(svg) {
    if (svg && typeof ExportStamp !== 'undefined' && ExportStamp.stampSvg) {
        try { return ExportStamp.stampSvg(svg); } catch (e) { return svg; }
    }
    return svg;
}

function runExportPreset(kind) {
    if (!editor || !editor.molecule || editor.molecule.isEmpty()) return;
    var haveIE = (typeof ImageExport !== 'undefined');

    // v2.4.7: journal-standard publication SVG -- white (opaque figure) or
    // transparent (drop onto any page / slide).
    if ((kind === 'pub-svg' || kind === 'pub-svg-transparent') &&
            haveIE && ImageExport.toPublicationSVG && ImageExport.downloadSVG) {
        var transparent = (kind === 'pub-svg-transparent');
        var svg = ImageExport.toPublicationSVG(editor.molecule,
            publicationExportOptions({ background: transparent ? 'transparent' : '#ffffff' }));
        ImageExport.downloadSVG(stampExportSvg(svg), (editor._exportFilename ? editor._exportFilename('svg') : 'molecule.svg'));
        if (editor.showInfo) editor.showInfo('Downloaded publication SVG' + (transparent ? ' (transparent)' : ''));
        return;
    }
    if (kind === 'pub-svg' && editor._exportPrintSVG) {   // graceful fallback
        editor._exportPrintSVG();
        return;
    }

    // v2.4.7: publication-quality transparent PNG (4x; >=1024px enforced by toPNG).
    if (kind === 'png' && haveIE && ImageExport.toPNG && ImageExport.downloadPNG) {
        if (editor.showInfo) editor.showInfo('Exporting transparent PNG...');
        ImageExport.toPNG(editor.molecule,
            publicationExportOptions({ background: 'transparent', scale: 4 }))
            .then(function(blob) {
                ImageExport.downloadPNG(blob, (editor._exportFilename ? editor._exportFilename('png') : 'molecule.png'));
                if (editor.showInfo) editor.showInfo('Downloaded transparent PNG (4x)');
            })
            .catch(function(err) {
                if (editor.showInfo) editor.showInfo('PNG export failed: ' + (err && err.message ? err.message : err));
            });
        return;
    }
    if (kind === 'png' && editor._exportPNGWithImageExport) {   // graceful fallback
        editor._exportPNGWithImageExport();
        return;
    }

    if (kind === 'compact-svg' && haveIE && ImageExport.toSVG) {
        var compact = ImageExport.toSVG(editor.molecule, {
            width: 420,
            height: 260,
            padding: 12,
            background: 'transparent',
            showAtomLabels: true
        });
        var out = document.getElementById('smiles-out');
        if (out) out.value = compact || '';
        setOutputTab('image', true);
        var image = document.getElementById('wb-output-image');
        if (image) setSafeSvg(image, compact, 'Preview unavailable.');
        if (editor.showInfo) editor.showInfo('Compact SVG ready in output');
    }
}

function setSafeSvg(container, svg, fallbackText) {
    if (!container) return;
    container.textContent = '';
    if (!svg) {
        var empty = document.createElement('span');
        empty.className = 'wb-output-muted';
        empty.textContent = fallbackText || '(no preview available)';
        container.appendChild(empty);
        return;
    }
    if (typeof DOMParser === 'undefined') {
        container.textContent = fallbackText || '(preview unavailable)';
        return;
    }
    try {
        var doc = new DOMParser().parseFromString(svg, 'image/svg+xml');
        var root = doc && doc.documentElement;
        if (!root || root.nodeName.toLowerCase() !== 'svg') {
            container.textContent = fallbackText || '(preview unavailable)';
            return;
        }
        var unsafe = root.querySelectorAll('script,foreignObject,iframe,object,embed');
        for (var i = unsafe.length - 1; i >= 0; i--) {
            unsafe[i].parentNode.removeChild(unsafe[i]);
        }
        var nodes = [root].concat(Array.prototype.slice.call(root.querySelectorAll('*')));
        for (var n = 0; n < nodes.length; n++) {
            var attrs = Array.prototype.slice.call(nodes[n].attributes || []);
            attrs.forEach(function(attr) {
                var attrName = attr.name.toLowerCase();
                var value = String(attr.value || '').trim().toLowerCase();
                if (attrName.indexOf('on') === 0 || value.indexOf('javascript:') === 0) {
                    nodes[n].removeAttribute(attr.name);
                }
            });
        }
        container.appendChild(document.importNode(root, true));
    } catch (e) {
        container.textContent = fallbackText || '(preview unavailable)';
    }
}

/* ---- v1.8.1: structure-preview tooltip on search-result hover ---- */
var _molSvgCache = {};
var _molSvgCacheKeys = [];

function rememberMolSvgCache(smiles, svg) {
    if (_molSvgCache[smiles] === undefined) _molSvgCacheKeys.push(smiles);
    _molSvgCache[smiles] = svg || null;
    while (_molSvgCacheKeys.length > 160) {
        var old = _molSvgCacheKeys.shift();
        if (old !== undefined) delete _molSvgCache[old];
    }
}

function _renderMolSvg(smiles, profile) {
    profile = profile || 'tooltip';
    var cacheKey = profile + '\n' + smiles;
    if (_molSvgCache[cacheKey] !== undefined) return _molSvgCache[cacheKey];
    if (typeof Molecule === 'undefined' || typeof SmilesParser === 'undefined' ||
        typeof Layout === 'undefined' || typeof ImageExport === 'undefined') {
        rememberMolSvgCache(cacheKey, null);
        return null;
    }
    try {
        var m = new Molecule();
        SmilesParser.parse(smiles, m);
        Layout.layout(m);
        var opts = profile === 'thumb'
            ? {
                width: 156,
                height: 104,
                padding: 10,
                background: 'transparent',
                showAtomLabels: true,
                showHydrogens: false,
                bondWidth: 2.2,
                doubleBondGap: 3.6,
                tripleBondGap: 4.4,
                fontSize: 11
            }
            : {
                width: 240,
                height: 174,
                padding: 10,
                background: 'transparent',
                showAtomLabels: true,
                showHydrogens: false,
                bondWidth: 2,
                fontSize: 12
            };
        var svg = ImageExport.toSVG(m, opts);
        rememberMolSvgCache(cacheKey, svg || null);
        return _molSvgCache[cacheKey];
    } catch (e) {
        rememberMolSvgCache(cacheKey, null);
        return null;
    }
}

function _renderMolSvgForThumb(smiles) {
    return _renderMolSvg(smiles, 'thumb');
}

function queueSearchThumb(container, smiles, seq) {
    queueMolThumb(container, smiles, function() { return seq === _searchSeq; });
}

function queueMolThumb(container, smiles, guard) {
    if (!container) return;
    container.textContent = 'Preview';
    _searchThumbQueue.push({ container: container, smiles: smiles, guard: guard });
    pumpSearchThumbQueue();
}

function pumpSearchThumbQueue() {
    if (_searchThumbFrame) return;
    var schedule = typeof requestIdleCallback === 'function'
        ? function(cb) { return requestIdleCallback(cb, { timeout: 120 }); }
        : (typeof requestAnimationFrame === 'function'
            ? function(cb) { return requestAnimationFrame(function() { cb({ timeRemaining: function() { return 8; } }); }); }
            : function(cb) { return setTimeout(function() { cb({ timeRemaining: function() { return 8; } }); }, 16); });
    _searchThumbFrame = schedule(function(deadline) {
        _searchThumbFrame = 0;
        var start = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
        var processed = 0;
        while (_searchThumbQueue.length > 0) {
            if (deadline && deadline.timeRemaining && deadline.timeRemaining() < 4 && processed > 0) break;
            if (processed >= 2) break;
            var job = _searchThumbQueue.shift();
            processed++;
            if (!job || (job.guard && !job.guard()) || !job.container ||
                    !document.documentElement.contains(job.container)) continue;
            setSafeSvg(job.container, _renderMolSvgForThumb(job.smiles), '(no preview)');
            var now = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
            if (now - start > 12) break;
        }
        if (_searchThumbQueue.length > 0) pumpSearchThumbQueue();
    });
}

function setMolTooltipPreview(svgEl, svg) {
    setSafeSvg(svgEl, svg, '(no preview available)');
}

function showMolTooltip(name, smiles, x, y) {
    var tt = document.getElementById('wb-mol-tooltip');
    if (!tt) return;
    var nameEl = tt.querySelector('.wb-mol-tooltip-name');
    var svgEl  = tt.querySelector('.wb-mol-tooltip-svg');
    var smiEl  = tt.querySelector('.wb-mol-tooltip-smi');
    if (nameEl) nameEl.textContent = name || smiles;
    if (smiEl)  smiEl.textContent  = smiles || '';
    if (svgEl) {
        if (_molSvgCache[smiles] !== undefined) {
            setMolTooltipPreview(svgEl, _molSvgCache[smiles]);
        } else {
            setSafeSvg(svgEl, null, 'Preview loading...');
            queueMolThumb(svgEl, smiles, function() {
                return tt.classList.contains('is-visible') && tt.getAttribute('data-smiles') === smiles;
            });
        }
    }
    tt.setAttribute('data-smiles', smiles || '');
    tt.classList.add('is-visible');
    tt.setAttribute('aria-hidden', 'false');
    moveMolTooltip(x, y);
}

function moveMolTooltip(x, y) {
    var tt = document.getElementById('wb-mol-tooltip');
    if (!tt || !tt.classList.contains('is-visible')) return;
    // Position to the right of the cursor by default; flip left if it
    // would overflow the viewport.
    var w = tt.offsetWidth || 240;
    var h = tt.offsetHeight || 200;
    var left = x + 14;
    var top  = y + 14;
    if (left + w > window.innerWidth - 8) left = x - w - 14;
    if (top + h > window.innerHeight - 8) top = window.innerHeight - h - 8;
    if (top < 8) top = 8;
    tt.style.left = Math.max(8, left) + 'px';
    tt.style.top  = top + 'px';
}

function hideMolTooltip() {
    var tt = document.getElementById('wb-mol-tooltip');
    if (!tt) return;
    tt.classList.remove('is-visible');
    tt.setAttribute('aria-hidden', 'true');
    tt.removeAttribute('data-smiles');
}

function runCompoundLookup() {
    var input = document.getElementById('bime-lookup-input');
    var modeEl = document.getElementById('bime-lookup-mode');
    var status = document.getElementById('bime-lookup-status');
    var list = document.getElementById('bime-lookup-results');
    if (!input || !status || !list) return;
    var query = input.value.trim();
    if (!query) {
        status.textContent = 'Enter a compound identifier';
        list.style.display = 'none';
        return;
    }
    var mode = (modeEl && modeEl.value) || 'auto';
    mode = mode === 'auto' ? detectCompoundLookupMode(query) : mode;
    var url = compoundLookupUrl(query, mode);
    if (!url) {
        status.textContent = 'Unsupported lookup query';
        list.style.display = 'none';
        return;
    }
    var seq = ++_lookupSeq;
    var t0 = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    status.style.color = '';
    status.textContent = 'Searching compounds...';
    list.innerHTML = '';
    list.style.display = 'none';
    fetchCompoundResults(url, mode, query)
        .then(function(hits) {
            if (seq !== _lookupSeq) return;
            hits = hits.slice(0, mode === 'formula' ? 10 : 6);
            var elapsed = Math.round(((typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now()) - t0);
            if (!hits.length) {
                status.textContent = 'No compound matches (' + elapsed + ' ms)';
                return;
            }
            status.textContent = hits.length + ' compound match' + (hits.length === 1 ? '' : 'es') +
                ' · ' + elapsed + ' ms';
            renderCompoundLookupResults(hits);
        })
        .catch(function(err) {
            if (seq !== _lookupSeq) return;
            status.style.color = '#dc2626';
            status.textContent = 'Compound lookup failed: ' + (err && err.message ? err.message : err);
        });
}

function fetchCompoundResults(url, mode, query) {
    return fetchJsonWithTimeout(url, 12000).then(function(data) {
        var cids = data && data.IdentifierList && data.IdentifierList.CID;
        if (Array.isArray(cids) && cids.length) {
            return fetchJsonWithTimeout(compoundPropertyUrlForCids(cids.slice(0, mode === 'formula' ? 10 : 6)), 12000)
                .then(function(propertyData) {
                    return normalizeCompoundResults(propertyData, query);
                });
        }
        return normalizeCompoundResults(data, query);
    });
}

function detectCompoundLookupMode(query) {
    var q = String(query || '').trim();
    if (/^(cid\s*)?\d+$/i.test(q)) return 'cid';
    if (/^InChI=/i.test(q) || /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/i.test(q)) return 'inchi';
    if (/^[A-Z][A-Za-z0-9]*$/.test(q) && /^([A-Z][a-z]?\d*)+$/.test(q)) return 'formula';
    return 'name';
}

function compoundLookupUrl(query, mode) {
    var base = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/';
    var q = String(query || '').trim();
    if (mode === 'cid') {
        var cid = q.replace(/^cid\s*/i, '').replace(/[^\d]/g, '');
        if (!cid) return '';
        return base + 'cid/' + encodeURIComponent(cid) + '/property/' + COMPOUND_LOOKUP_PROPS + '/JSON';
    }
    if (mode === 'formula') {
        var formula = q.replace(/\s+/g, '');
        return base + 'fastformula/' + encodeURIComponent(formula) + '/cids/JSON?MaxRecords=10&MaxSeconds=5';
    }
    if (mode === 'inchi') {
        if (/^InChI=/i.test(q)) {
            return base + 'inchi/cids/JSON?inchi=' + encodeURIComponent(q);
        }
        return base + 'inchikey/' + encodeURIComponent(q) + '/property/' + COMPOUND_LOOKUP_PROPS + '/JSON';
    }
    return base + 'name/' + encodeURIComponent(q) + '/property/' + COMPOUND_LOOKUP_PROPS + '/JSON';
}

function compoundPropertyUrlForCids(cids) {
    var cidList = (cids || []).map(function(cid) {
        return String(cid).replace(/[^\d]/g, '');
    }).filter(Boolean).join(',');
    if (!cidList) return '';
    return 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' +
        encodeURIComponent(cidList) + '/property/' + COMPOUND_LOOKUP_PROPS + '/JSON';
}

function fetchJsonWithTimeout(url, timeoutMs) {
    var controller = typeof AbortController !== 'undefined' ? new AbortController() : null;
    var timer = controller ? setTimeout(function() { controller.abort(); }, timeoutMs || 12000) : 0;
    return fetch(url, {
        headers: { Accept: 'application/json' },
        signal: controller ? controller.signal : undefined
    }).then(function(resp) {
        if (timer) clearTimeout(timer);
        if (!resp.ok) throw new Error('HTTP ' + resp.status);
        return resp.json();
    }, function(err) {
        if (timer) clearTimeout(timer);
        throw err;
    });
}

function normalizeCompoundResults(data, fallbackName) {
    var props = data && data.PropertyTable && data.PropertyTable.Properties;
    if (!Array.isArray(props)) return [];
    var seen = {};
    return props.map(function(row) {
        var smiles = row.SMILES || row.IsomericSMILES || row.CanonicalSMILES || row.ConnectivitySMILES || '';
        return {
            cid: row.CID || '',
            name: row.Title || row.IUPACName || fallbackName || 'Compound',
            smiles: smiles,
            formula: row.MolecularFormula || '',
            molecularWeight: row.MolecularWeight || '',
            iupacName: row.IUPACName || '',
            xlogp: row.XLogP,
            exactMass: row.ExactMass,
            monoisotopicMass: row.MonoisotopicMass,
            tpsa: row.TPSA,
            complexity: row.Complexity
        };
    }).filter(function(hit) {
        if (!hit.smiles || seen[hit.cid || hit.smiles]) return false;
        seen[hit.cid || hit.smiles] = true;
        return true;
    });
}

function renderCompoundLookupResults(hits) {
    var list = document.getElementById('bime-lookup-results');
    if (!list) return;
    var seq = _lookupSeq;
    list.innerHTML = '';
    list.style.display = 'block';
    hits.forEach(function(hit) {
        var li = document.createElement('li');
        li.className = 'wb-search-hit';
        li.setAttribute('role', 'button');
        li.setAttribute('tabindex', '0');
        li.setAttribute('aria-label', 'Load ' + hit.name + ', compound CID ' + hit.cid);

        var thumb = document.createElement('div');
        thumb.className = 'wb-search-hit-thumb';
        thumb.setAttribute('aria-hidden', 'true');
        queueMolThumb(thumb, hit.smiles, function() { return seq === _lookupSeq; });

        var body = document.createElement('div');
        body.className = 'wb-search-hit-body';
        var title = document.createElement('div');
        title.className = 'wb-search-hit-title';
        title.textContent = hit.name;
        var meta = document.createElement('div');
        meta.className = 'wb-search-hit-meta';
        meta.textContent = [
            hit.cid ? 'CID ' + hit.cid : '',
            hit.formula || '',
            hit.molecularWeight ? 'MW ' + formatProperty(hit.molecularWeight) : '',
            hit.xlogp !== undefined ? 'XLogP ' + formatProperty(hit.xlogp) : ''
        ].filter(Boolean).join(' · ');
        var smi = document.createElement('div');
        smi.className = 'wb-search-hit-smiles';
        smi.textContent = hit.smiles;
        body.appendChild(title);
        body.appendChild(meta);
        body.appendChild(smi);
        li.appendChild(thumb);
        li.appendChild(body);
        li.addEventListener('click', function() { loadCompoundHit(hit); });
        li.addEventListener('keydown', function(ev) {
            if (ev.key === 'Enter' || ev.key === ' ') {
                ev.preventDefault();
                loadCompoundHit(hit);
            }
        });
        li.addEventListener('mouseenter', function(ev) { showMolTooltip(hit.name, hit.smiles, ev.clientX, ev.clientY); });
        li.addEventListener('mousemove', function(ev) { moveMolTooltip(ev.clientX, ev.clientY); });
        li.addEventListener('mouseleave', hideMolTooltip);
        li.addEventListener('focus', function() {
            var rect = li.getBoundingClientRect();
            showMolTooltip(hit.name, hit.smiles, rect.right, rect.top);
        });
        li.addEventListener('blur', hideMolTooltip);
        list.appendChild(li);
    });
}

function loadCompoundHit(hit) {
    if (!editor || !hit || !hit.smiles) return;
    editor.readGenericMolecularInput(hit.smiles);
    if (editor.molecule) {
        editor.molecule.name = hit.name || hit.iupacName || '';
        editor.molecule.data = editor.molecule.data || {};
        editor.molecule.data.compoundLookup = hit;
        editor.renderer.showMolName = !!editor.molecule.name;
        editor.render();
    }
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
    var status = document.getElementById('bime-lookup-status');
    if (status) status.textContent = 'Loaded ' + (hit.name || ('CID ' + hit.cid));
}

function currentCompoundLookupData() {
    return editor && editor.molecule && editor.molecule.data && editor.molecule.data.compoundLookup;
}

function reactionInsightsFromSmiles(smiles) {
    if (!smiles || smiles.indexOf('>>') < 0 || typeof SmilesParser === 'undefined' || typeof Molecule === 'undefined') return null;
    var sides = smiles.split('>>');
    var reactants = (sides[0] || '').split('.').filter(Boolean);
    var products = (sides[1] || '').split('.').filter(Boolean);
    var rCounts = {};
    var pCounts = {};
    var rLabels = reactants.map(function(smi) { return componentFormulaLabel(smi, rCounts); }).filter(Boolean);
    var pLabels = products.map(function(smi) { return componentFormulaLabel(smi, pCounts); }).filter(Boolean);
    return {
        reactants: rLabels,
        products: pLabels,
        balanced: sameFormulaCounts(rCounts, pCounts)
    };
}

// v2.0.16: coefficient-collapsed reaction equation for the Editor Output.
// Splits a reaction SMILES into per-side component molecules, then defers the
// grouping/balance/formatting to the bundled, unit-tested
// Molecule.reactionEquation. Identical species are grouped by canonical SMILES
// so "[H][H].[H][H].O=O>>O.O" renders as "2 H2 + O2 -> 2 H2O".
function reactionFormulaFromSmiles(smiles) {
    if (!smiles || smiles.indexOf('>>') < 0) return null;
    if (typeof SmilesParser === 'undefined' || typeof Molecule === 'undefined' ||
        typeof Molecule.reactionEquation !== 'function') return null;
    var sides = smiles.split('>>');
    var rSmis = (sides[0] || '').split('.').filter(Boolean);
    var pSmis = (sides[sides.length - 1] || '').split('.').filter(Boolean);
    if (!rSmis.length && !pSmis.length) return null;
    function parseAll(list) {
        var out = [];
        for (var i = 0; i < list.length; i++) {
            try {
                var m = new Molecule();
                SmilesParser.parse(list[i], m);
                if (m.atoms.length) out.push(m);
            } catch (e) { /* skip unparseable component */ }
        }
        return out;
    }
    function canonKey(m) {
        try {
            return (typeof SmilesWriter !== 'undefined' && SmilesWriter.write)
                ? SmilesWriter.write(m)
                : Molecule.formulaString(m.elementCounts());
        } catch (e) {
            return Molecule.formulaString(m.elementCounts());
        }
    }
    var rMols = parseAll(rSmis), pMols = parseAll(pSmis);
    if (!rMols.length && !pMols.length) return null;
    return Molecule.reactionEquation(rMols, pMols, canonKey);
}

// v2.0.16: paint the always-visible reaction-equation readout under the Editor
// Output. Shows the coefficient equation + a balance badge when a reaction
// arrow is present, and a discoverability hint when several disconnected
// fragments exist but no arrow has been drawn yet.
function updateReactionEquationReadout() {
    var el = document.getElementById('wb-reaction-equation');
    if (!el) return;
    var mol = editor && editor.molecule;
    var hasArrow = !!(mol && mol.reactionArrow);
    var fragments = (mol && mol.getComponents) ? mol.getComponents().length : 0;
    el.textContent = '';
    if (hasArrow) {
        var info = reactionFormulaFromSmiles(editor && editor.smiles ? editor.smiles() : '');
        if (info && (info.reactants.length || info.products.length)) {
            el.hidden = false;
            var eq = document.createElement('span');
            eq.className = 'wb-reaction-equation-text';
            eq.textContent = info.equation;
            var badge = document.createElement('span');
            badge.className = 'wb-reaction-equation-badge ' +
                (info.balanced ? 'is-balanced' : 'is-unbalanced');
            badge.textContent = info.balanced ? 'Balanced' : 'Unbalanced';
            el.appendChild(eq);
            el.appendChild(badge);
            return;
        }
    }
    if (!hasArrow && fragments >= 2) {
        el.hidden = false;
        var hint = document.createElement('span');
        hint.className = 'wb-reaction-equation-hint';
        hint.textContent = 'Add a → reaction arrow (left rail ▸ Reaction ▸ Arrow) ' +
            'to turn these ' + fragments + ' fragments into a stoichiometric equation.';
        el.appendChild(hint);
        return;
    }
    el.hidden = true;
}

function componentFormulaLabel(smiles, accum) {
    try {
        var mol = new Molecule();
        SmilesParser.parse(smiles, mol);
        var counts = formulaCountsForMol(mol);
        Object.keys(counts).forEach(function(k) { accum[k] = (accum[k] || 0) + counts[k]; });
        return formulaFromCounts(counts);
    } catch (e) {
        return '';
    }
}

function sameFormulaCounts(a, b) {
    var keys = {};
    Object.keys(a || {}).forEach(function(k) { keys[k] = true; });
    Object.keys(b || {}).forEach(function(k) { keys[k] = true; });
    return Object.keys(keys).every(function(k) { return (a[k] || 0) === (b[k] || 0); });
}

/**
 * v1.5.1: Validate now reads from the live canvas (smiles-out / editor.smiles())
 * rather than a removed Import-Structure textarea.
 */
function doValidate() {
    var el = document.getElementById('validate-result');
    var v = ((editor && editor.smiles && editor.smiles()) || '').trim();
    if (!v) { if (el) el.textContent = 'Draw or load a molecule first'; return; }
    if (!editor || typeof editor.validateSmiles !== 'function') { return; }
    var r = editor.validateSmiles(v);
    if (!el) return;
    if (r.valid) {
        el.style.color = r.warnings && r.warnings.length ? '#ca8a04' : '#059669';
        el.textContent = '\u2713 Valid \u2014 ' + r.atoms + ' atoms' + (r.warnings && r.warnings.length ? ' \u26a0' : '');
    } else {
        el.style.color = '#dc2626';
        el.textContent = '\u2717 ' + (r.error || 'Invalid');
    }
}


/* ---- v1.5.0: Search-type radio + source toggle + Search button wired to
       editor.searchLibrary(). Radio cards pick the algorithm; the source
       toggle picks the target set (built-in 1181 OR user-supplied max-100);
       clicking a hit loads it AND highlights the matched part on the canvas
       via editor.highlightSearchMatch(). */
window.BIME_SEARCH_TYPE = 'substructure';
window.BIME_SEARCH_SOURCE = 'builtin';
window.BIME_USER_LIBRARY = [];
window.BIME_LAST_QUERY = '';

document.addEventListener('DOMContentLoaded', function() {
    // v1.5.1: search type now lives in a single dropdown below the editor.
    var modeSelect = document.getElementById('bime-search-mode');
    if (modeSelect) {
        window.BIME_SEARCH_TYPE = modeSelect.value || 'substructure';
        modeSelect.addEventListener('change', function() {
            window.BIME_SEARCH_TYPE = modeSelect.value || 'substructure';
        });
    }
    var srcRadios = document.querySelectorAll('input[name="bime-search-source"]');
    for (var j = 0; j < srcRadios.length; j++) {
        srcRadios[j].addEventListener('change', function(e) {
            if (!e.target.checked) return;
            window.BIME_SEARCH_SOURCE = e.target.value;
            var zone = document.getElementById('bime-user-lib-zone');
            if (zone) { zone.style.display = (e.target.value === 'user') ? 'block' : 'none'; }
        });
    }
    var btn = document.getElementById('bime-search-run');
    if (btn) { btn.addEventListener('click', runLibrarySearch); }

    var pickBtn = document.getElementById('bime-user-lib-pick');
    var fileIn = document.getElementById('bime-user-lib-input');
    if (pickBtn && fileIn) {
        pickBtn.addEventListener('click', function() { fileIn.click(); });
        pickBtn.addEventListener('dragover', function(e) { e.preventDefault(); pickBtn.style.background = 'rgba(13,148,136,0.18)'; });
        pickBtn.addEventListener('dragleave', function() { pickBtn.style.background = 'rgba(13,148,136,0.04)'; });
        pickBtn.addEventListener('drop', function(e) {
            e.preventDefault();
            pickBtn.style.background = 'rgba(13,148,136,0.04)';
            var f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
            if (f) { loadUserLibraryFile(f); }
        });
        fileIn.addEventListener('change', function() {
            if (fileIn.files && fileIn.files[0]) { loadUserLibraryFile(fileIn.files[0]); }
        });
    }
});

/**
 * Parse a multi-record .smi / .sdf / .mol / .txt file into an array of
 * {name, smiles}. Caps at MolEditor.USER_LIBRARY_LIMIT (100). Stores into
 * window.BIME_USER_LIBRARY and updates the status line.
 *
 * Format heuristics:
 *   - .sdf with $$$$ separators: split on $$$$, parse each block as MOL,
 *     extract a SMILES via SmilesWriter, use the first non-blank line as
 *     the name.
 *   - line-per-molecule (.smi / .smiles / .txt): each line is "SMILES name"
 *     (whitespace-separated; rest of line = name; if no name, use line
 *     index).
 *   - single .mol: one entry.
 */
function validateWorkbenchImportFile(file) {
    if (typeof MolEditor !== 'undefined' && MolEditor.validateImportFile) {
        return MolEditor.validateImportFile(file);
    }
    if (!file) return 'No file selected';
    var maxBytes = 5 * 1024 * 1024;
    if (typeof file.size === 'number' && file.size > maxBytes) {
        return 'File is too large for local import (max 5 MB)';
    }
    var name = String(file.name || '').toLowerCase();
    var ext = name.indexOf('.') >= 0 ? name.split('.').pop() : '';
    if (ext && !/^(mol|sdf|rxn|smi|smiles|txt)$/.test(ext)) {
        return 'Unsupported molecular file type .' + ext;
    }
    return '';
}

function readWorkbenchMolecularTextFile(file, onload, onerror) {
    if (typeof MolEditor !== 'undefined' && MolEditor.readMolecularTextFile) {
        return MolEditor.readMolecularTextFile(file, onload, onerror);
    }
    var err = validateWorkbenchImportFile(file);
    if (err) {
        if (onerror) onerror(new Error(err));
        return false;
    }
    var reader = new FileReader();
    reader.onload = function(ev) {
        var text = ev && ev.target ? ev.target.result : '';
        onload(typeof text === 'string' ? text : String(text || ''));
    };
    reader.onerror = function() {
        if (onerror) onerror(new Error('Could not read file'));
    };
    reader.readAsText(file);
    return true;
}

function loadUserLibraryFile(file) {
    var statusEl = document.getElementById('bime-user-lib-status');
    readWorkbenchMolecularTextFile(file, function(text) {
        var entries = parseMultiMolFile(text, file.name);
        var cap = (typeof MolEditor !== 'undefined' && MolEditor.USER_LIBRARY_LIMIT) || 100;
        var dropped = Math.max(0, entries.length - cap);
        entries = entries.slice(0, cap);
        window.BIME_USER_LIBRARY = entries;
        if (statusEl) {
            statusEl.textContent = entries.length + ' molecule' + (entries.length === 1 ? '' : 's') +
                ' loaded' + (dropped > 0 ? ' (' + dropped + ' over the 100-cap dropped)' : '');
        }
    }, function(err) {
        if (statusEl) { statusEl.textContent = err && err.message ? err.message : 'Could not read file'; }
    });
}

function parseMultiMolFile(text, filename) {
    var out = [];
    if (!text) return out;
    // Apply the 100-cap during parsing so a pathological 100k-line file
    // doesn't allocate a transient ~5-10MB intermediate array. The
    // _buildUserLibraryCache call downstream re-applies the same cap as a
    // belt-and-braces guard; this just stops parsing early when we have
    // enough.
    var cap = (typeof MolEditor !== 'undefined' && MolEditor.USER_LIBRARY_LIMIT)
        ? MolEditor.USER_LIBRARY_LIMIT : 100;
    var lower = (filename || '').toLowerCase();
    if (/\.sdf$/.test(lower) || /\$\$\$\$/.test(text)) {
        // SDF with $$$$ separators. Prefer the real MOL parser + SMILES
        // writer path; fall back to scanning for SMILES-bearing property
        // lines for lightweight vendor SDF exports.
        var blocks = text.split(/\n\$\$\$\$\s*\n?/);
        for (var i = 0; i < blocks.length && out.length < cap; i++) {
            var block = blocks[i].trim();
            if (!block) continue;
            var parsed = parseSdfBlockToLibraryEntry(block, i + 1);
            if (parsed) {
                out.push(parsed);
                continue;
            }
            var firstLine = firstNonEmptyLine(block) || ('record-' + (i + 1));
            var lines = block.split(/\n/);
            for (var li = 0; li < lines.length; li++) {
                var ln = lines[li].trim();
                if (!ln || /^M\s+END/.test(ln) || /^>/.test(ln)) continue;
                if (looksLikeSmiles(ln)) {
                    out.push({ name: firstLine, smiles: ln });
                    break;
                }
            }
        }
        return out;
    }
    // Otherwise treat as line-per-molecule .smi/.smiles/.txt.
    var rows = text.split(/\r?\n/);
    var idx = 0;
    for (var r = 0; r < rows.length && out.length < cap; r++) {
        var row = rows[r].trim();
        if (!row || row.charAt(0) === '#') continue;
        var parts = row.split(/\s+/);
        var smi = parts[0];
        if (!looksLikeSmiles(smi)) continue;
        idx++;
        var nameRest = parts.slice(1).join(' ').trim() || ('mol-' + idx);
        out.push({ name: nameRest, smiles: smi });
    }
    return out;
}

function firstNonEmptyLine(text) {
    var lines = String(text || '').split(/\r?\n/);
    for (var i = 0; i < lines.length; i++) {
        var line = lines[i].trim();
        if (line) return line;
    }
    return '';
}

function sdfFieldValue(block, fieldName) {
    var re = new RegExp('^>\\\\s*<'+ fieldName +'>\\\\s*(?:\\r\\n|\\r|\\n)([^\\r\\n]*)', 'mi');
    var m = String(block || '').match(re);
    return m && m[1] ? m[1].trim() : '';
}

function parseSdfBlockToLibraryEntry(block, index) {
    if (!editor || typeof Molecule === 'undefined' || typeof SmilesWriter === 'undefined') return null;
    if (!editor._parseMolV2000 || !editor._parseMolV3000) return null;
    var molEnd = block.indexOf('M  END');
    if (molEnd < 0) return null;
    var molPart = block.substring(0, molEnd + 6);
    var mol = new Molecule();
    try {
        if (molPart.indexOf('V3000') >= 0) {
            editor._parseMolV3000(molPart, mol);
        } else {
            editor._parseMolV2000(molPart, mol);
        }
        if (!mol.atoms || mol.atoms.length === 0) return null;
        var smiles = SmilesWriter.write(mol);
        if (!smiles) return null;
        return {
            name: sdfFieldValue(block, 'NAME') || mol.name || firstNonEmptyLine(block) || ('record-' + index),
            smiles: smiles
        };
    } catch (e) {
        return null;
    }
}

function looksLikeSmiles(s) {
    if (!s) return false;
    if (/[A-Za-z0-9@+\-\[\]\(\)=#\.\\\/]/.test(s) === false) return false;
    // Reject obvious non-SMILES (tab-separated header lines, MOL header).
    if (/\s/.test(s)) return false;
    if (/^M\s+/.test(s)) return false;
    return s.length >= 1 && s.length < 1000;
}

function runLibrarySearch() {
    openWorkbenchSection('bime-search-section');
    scheduleLibraryWarmup();
    var statusEl = document.getElementById('bime-search-status');
    var listEl = document.getElementById('bime-search-results');
    if (!editor || typeof editor.searchLibrary !== 'function') {
        if (statusEl) { statusEl.textContent = 'searchLibrary unavailable (older bundle)'; }
        return;
    }
    var query = (editor.smiles && editor.smiles()) || '';
    if (!query) {
        statusEl.textContent = 'Draw a molecule first or paste a SMILES';
        listEl.style.display = 'none';
        return;
    }
    var mode = window.BIME_SEARCH_TYPE || 'substructure';
    var source = window.BIME_SEARCH_SOURCE || 'builtin';
    if (source === 'builtin' && !builtInLibraryReady()) {
        statusEl.style.color = '';
        statusEl.textContent = 'Built-in library is still loading... search will run automatically.';
        listEl.style.display = 'none';
        listEl.innerHTML = '';
        waitForBuiltInLibraryAndSearch(0);
        return;
    }
    var opts = { topN: 10, threshold: 0.3, timeoutMs: 10000 };
    if (source === 'user') {
        if (!window.BIME_USER_LIBRARY || window.BIME_USER_LIBRARY.length === 0) {
            statusEl.textContent = 'Drop a .smi or .sdf file first (max 100 molecules)';
            return;
        }
        opts.targets = window.BIME_USER_LIBRARY;
    }
    // v1.8.5: ring-match-only toggle. Read once per Search click so the
    // user can flip it between runs. Only meaningful for substructure
    // and MCS modes; quietly ignored for exact and similarity.
    var rmroEl = document.getElementById('bime-rmro');
    if (rmroEl && rmroEl.checked) {
        opts.ringMatchesRingOnly = true;
    }
    window.BIME_LAST_QUERY = query;
    var label = (source === 'user')
        ? (window.BIME_USER_LIBRARY.length + ' user molecules')
        : '1181 molecules';
    statusEl.style.color = '';
    statusEl.textContent = 'Searching ' + label + ' in ' + mode + ' mode...';
    listEl.style.display = 'none';
    listEl.innerHTML = '';
    _searchThumbQueue = [];
    var t0 = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    var seq = ++_searchSeq;
    opts.yieldEvery = mode === 'mcs' ? 6 : 64;
    opts.onProgress = function(done, total) {
        if (seq !== _searchSeq || !statusEl || !total) return;
        statusEl.textContent = 'Searching ' + label + ' in ' + mode + ' mode... ' +
            Math.min(total, done) + '/' + total;
    };
    afterPaint(function() {
        if (seq !== _searchSeq) return;
        editor.searchLibrary(query, mode, opts)
        .then(function(hits) {
            if (seq !== _searchSeq) return;
            var now = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
            var elapsed = Math.round(now - t0);
            if (!hits || hits.length === 0) {
                statusEl.textContent = 'No matches in ' + mode + ' mode (' + elapsed + ' ms)';
                listEl.style.display = 'none';
                return;
            }
            statusEl.textContent = hits.length + ' hit' + (hits.length === 1 ? '' : 's') +
                ' · ' + mode + ' · ' + elapsed + ' ms · click a hit to load + highlight';
            listEl.style.display = 'block';
            listEl.innerHTML = '';
            var frag = document.createDocumentFragment();
            hits.forEach(function(h) {
                var li = document.createElement('li');
                li.className = 'wb-search-hit';
                // v1.8.1: native title removed — replaced by hover preview
                // tooltip below. Keep as data-attr for accessibility.
                li.setAttribute('data-smiles', h.smiles);
                li.setAttribute('data-name', h.name);
                li.setAttribute('role', 'button');
                li.setAttribute('tabindex', '0');
                li.setAttribute('aria-label',
                    h.name + ', SMILES ' + h.smiles + ', score ' + h.score.toFixed(3));
                var lbl = h.name + ' (' + h.score.toFixed(3) + ')';
                if (h.mcsSize) { lbl += ' · MCS=' + h.mcsSize; }
                if (h.matchCount) { lbl += ' · matches=' + h.matchCount; }

                var thumb = document.createElement('div');
                thumb.className = 'wb-search-hit-thumb';
                thumb.setAttribute('aria-hidden', 'true');
                queueSearchThumb(thumb, h.smiles, seq);

                var body = document.createElement('div');
                body.className = 'wb-search-hit-body';
                var title = document.createElement('div');
                title.className = 'wb-search-hit-title';
                title.textContent = h.name;
                var meta = document.createElement('div');
                meta.className = 'wb-search-hit-meta';
                meta.textContent = lbl.replace(h.name + ' ', '');
                var smi = document.createElement('div');
                smi.className = 'wb-search-hit-smiles';
                smi.textContent = h.smiles;
                body.appendChild(title);
                body.appendChild(meta);
                body.appendChild(smi);
                li.appendChild(thumb);
                li.appendChild(body);
                li.addEventListener('click', function() {
                    editor.readGenericMolecularInput(h.smiles);
                    if (typeof queueOutputUpdate === 'function') { queueOutputUpdate(); }
                    if (typeof scheduleWorkbenchUpdate === 'function') { scheduleWorkbenchUpdate(); }
                    if (typeof editor.highlightSearchMatch === 'function') {
                        var highlightOpts = {};
                        var rmroHighlightEl = document.getElementById('bime-rmro');
                        if (rmroHighlightEl && rmroHighlightEl.checked) {
                            highlightOpts.ringMatchesRingOnly = true;
                        }
                        var n = editor.highlightSearchMatch(mode, window.BIME_LAST_QUERY, highlightOpts);
                        if (n > 0) {
                            statusEl.textContent = 'Loaded ' + h.name + ' · ' + n +
                                ' matched atoms highlighted';
                        }
                    }
                });
                li.addEventListener('keydown', function(ev) {
                    if (ev.key === 'Enter' || ev.key === ' ') {
                        ev.preventDefault();
                        li.click();
                    }
                });
                // v1.8.1: structure-preview hover.
                li.addEventListener('mouseenter', function (ev) {
                    showMolTooltip(h.name, h.smiles, ev.clientX, ev.clientY);
                });
                li.addEventListener('mousemove', function (ev) {
                    moveMolTooltip(ev.clientX, ev.clientY);
                });
                li.addEventListener('mouseleave', function () {
                    hideMolTooltip();
                });
                li.addEventListener('focus', function () {
                    var rect = li.getBoundingClientRect();
                    showMolTooltip(h.name, h.smiles, rect.right, rect.top);
                });
                li.addEventListener('blur', function () {
                    hideMolTooltip();
                });
                frag.appendChild(li);
            });
            listEl.appendChild(frag);
        })
        .catch(function(err) {
            if (seq !== _searchSeq) return;
            statusEl.style.color = '#dc2626';
            statusEl.textContent = 'Search failed: ' + (err && err.message ? err.message : err);
            listEl.style.display = 'none';
        });
    });
}

function builtInLibraryReady() {
    return typeof COMMON_MOLECULES !== 'undefined' &&
        COMMON_MOLECULES &&
        COMMON_MOLECULES.length > 0;
}

function waitForBuiltInLibraryAndSearch(attempt) {
    if (builtInLibraryReady()) {
        runLibrarySearch();
        return;
    }
    ensureBuiltInLibraryScript().then(function() {
        runLibrarySearch();
    }).catch(function(err) {
        var statusEl = document.getElementById('bime-search-status');
        if (statusEl) {
            statusEl.style.color = '#dc2626';
            statusEl.textContent = (err && err.message) ? err.message : 'Built-in library did not finish loading.';
        }
    });
}

/* ---- SMILES input + Load button ---- */
function loadFromSmilesInput() {
    var inp = document.getElementById('smiles-input');
    if (!inp) return;
    var v = inp.value.trim();
    if (!v || !editor) return;
    editor.readGenericMolecularInput(v);
    queueOutputUpdate();
    scheduleWorkbenchUpdate();
}

/* ---- Paste from Clipboard ---- */
function pasteFromClipboard() {
    if (!editor) return;
    if (!navigator.clipboard || !navigator.clipboard.readText) {
        if (editor.showInfo) editor.showInfo('Clipboard API unavailable');
        return;
    }
    navigator.clipboard.readText().then(function(text) {
        var t = (text || '').trim();
        if (!t) return;
        editor.readGenericMolecularInput(t);
        queueOutputUpdate();
        scheduleWorkbenchUpdate();
    }).catch(function() {
        if (editor.showInfo) editor.showInfo('Clipboard read failed');
    });
}

/* ---- Drag-drop / file-input wiring ---- */
document.addEventListener('DOMContentLoaded', function() {
    var dz = document.getElementById('bime-dropzone');
    var fi = document.getElementById('bime-file-input');
    if (!dz || !fi) return;

    function handleFile(file) {
        if (!file || !editor) return;
        readWorkbenchMolecularTextFile(file, function(text) {
            editor.readGenericMolecularInput(text);
            queueOutputUpdate();
            scheduleWorkbenchUpdate();
        }, function(err) {
            if (editor.showInfo) {
                editor.showInfo(err && err.message ? err.message : 'Could not read file');
            }
        });
    }

    dz.addEventListener('click', function() { fi.click(); });
    dz.addEventListener('keydown', function(e) {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); fi.click(); }
    });
    fi.addEventListener('change', function(e) {
        var f = e.target.files && e.target.files[0];
        handleFile(f);
        fi.value = '';
    });

    ['dragenter', 'dragover'].forEach(function(ev) {
        dz.addEventListener(ev, function(e) {
            e.preventDefault(); e.stopPropagation();
            dz.classList.add('is-dragging');
        });
    });
    ['dragleave', 'drop'].forEach(function(ev) {
        dz.addEventListener(ev, function(e) {
            e.preventDefault(); e.stopPropagation();
            if (ev === 'dragleave' && e.target !== dz) return;
            dz.classList.remove('is-dragging');
        });
    });
    dz.addEventListener('drop', function(e) {
        var f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
        handleFile(f);
    });

    // Prevent the browser from navigating away on stray drops outside the zone.
    window.addEventListener('dragover', function(e) { e.preventDefault(); });
    window.addEventListener('drop', function(e) {
        if (!dz.contains(e.target)) e.preventDefault();
    });
});
