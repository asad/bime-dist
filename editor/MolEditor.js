/**
 * BIME — Modern SVG Molecule Editor
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Main editor class. Drop-in API-compatible replacement for JSApplet.JSME.
 * Creates an SVG-based drawing canvas with toolbar, tools, and chemical I/O.
 */
(function(global) {
    'use strict';

    var BRAND_COLOR = '#0d9488';

    /** Escape user-controlled strings before inserting into innerHTML. */
    function escapeHTML(s) {
        return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
    }

    // =========================================================================
    // Toolbar layout — grouped by chemistry category
    // v1.5.2: Each group has a `placement` of 'top' or 'left'. Top-placed
    // groups render in the horizontal toolbar above the canvas (most-used
    // bond/edit/view actions). Left-placed groups render in a vertical
    // icon-only rail beside the canvas (mode-switching tools, reaction,
    // analysis, export). Splitting the bar this way matches the layout of
    // ChemDraw / Marvin / Ketcher and keeps the top row uncluttered. The
    // group order, ids, and items are unchanged so that ToolbarPrefs (the
    // Customize panel) and the (i) info popovers continue to work.
    // =========================================================================
    var TOOLBAR_GROUPS = [
        {
            id: 'draw', label: 'Draw', placement: 'left',
            items: [
                { id: 'bond',   label: 'Bond',  icon: 'bond',  type: 'tool' },
                { id: 'chain',  label: 'Chain', icon: 'chain', type: 'tool' },
            ]
        },
        {
            id: 'bonds', label: 'Bonds', placement: 'top',
            items: [
                { id: 'single', label: 'Single', icon: 'single', type: 'bond', bondType: 1 },
                { id: 'double', label: 'Double', icon: 'double', type: 'bond', bondType: 2 },
                { id: 'triple', label: 'Triple', icon: 'triple', type: 'bond', bondType: 3 },
                { id: 'stereo', label: 'Wedge',  icon: 'stereo', type: 'tool' },
            ]
        },
        {
            id: 'rings', label: 'Rings', collapsed: true, placement: 'left',
            items: [
                { id: 'ring3', label: '3',  icon: 'ring3', type: 'ring', size: 3 },
                { id: 'ring4', label: '4',  icon: 'ring4', type: 'ring', size: 4 },
                { id: 'ring5', label: '5',  icon: 'ring5', type: 'ring', size: 5 },
                { id: 'ring6', label: '6',  icon: 'ring6', type: 'ring', size: 6 },
                { id: 'ring7', label: '7',  icon: 'ring7', type: 'ring', size: 7 },
            ]
        },
        {
            id: 'edit', label: 'Edit', placement: 'top',
            items: [
                { id: 'select', label: 'Select', icon: 'select', type: 'tool' },
                { id: 'delete', label: 'Delete', icon: 'del',   type: 'tool' },
                { id: 'move',   label: 'Move',   icon: 'movea', type: 'tool' },
                { id: 'undo',   label: 'Undo',   icon: 'undo',  type: 'action' },
                { id: 'redo',   label: 'Redo',   icon: 'redo',  type: 'action' },
                { id: 'clear',  label: 'Clear',  icon: 'clean', type: 'action' },
            ]
        },
        {
            id: 'view', label: 'View', placement: 'top',
            items: [
                { id: 'autolayout', label: 'Layout',  icon: 'autolayout', type: 'action' },
                { id: 'zoomin',     label: 'Zoom +',  icon: 'zoomin',     type: 'action' },
                { id: 'zoomout',    label: 'Zoom \u2013', icon: 'zoomout',   type: 'action' },
                { id: 'zoomfit',    label: 'Fit',     icon: 'zoomfit',    type: 'action' },
                { id: 'ciprs',      label: 'R/S',     icon: 'ciprs',      type: 'action' },
                { id: 'cipez',      label: 'E/Z',     icon: 'cipez',      type: 'action' },
                { id: 'toggleh',    label: 'H',       icon: 'toggleh',    type: 'action' },
                { id: 'togglearo',  label: 'Ar',      icon: 'togglearo',  type: 'action' },
            ]
        },
        {
            id: 'reaction', label: 'Reaction', placement: 'left',
            items: [
                { id: 'reaction', label: 'Arrow', icon: 'react', type: 'tool' },
                { id: 'mapping',  label: 'Map',   icon: '123',   type: 'tool' },
                { id: 'togglemap', label: 'Map #', icon: '123', type: 'action' },
                { id: 'rdtautomap', label: 'Auto-map', icon: 'rdtautomap', type: 'action' },
                { id: 'togglecolors', label: 'Colors', icon: 'togglecolors', type: 'action' },
                { id: 'togglepairs', label: 'Pairs', icon: 'togglepairs', type: 'action' },
            ]
        },
        {
            id: 'actions', label: '', placement: 'left',
            items: [
                { id: 'validate', label: 'Validate', icon: 'validate', type: 'action' },
                { id: 'smarts',   label: 'SMARTS',   icon: 'smarts',   type: 'action' },
                { id: 'mcs',      label: 'MCS',      icon: 'smsd',     type: 'action' },
                { id: 'simsearch', label: 'Sim',     icon: 'simsearch', type: 'action' },
            ]
        },
        {
            id: 'export', label: 'Export', placement: 'left',
            items: [
                { id: 'name',      label: 'Name',  icon: 'name',    type: 'action' },
                { id: 'exportsvg', label: 'SVG',   icon: 'svg',     type: 'action' },
                { id: 'exportpng', label: 'PNG',   icon: 'png',     type: 'action' },
                { id: 'exportprint', label: 'Print', icon: 'print', type: 'action' },
                { id: 'copyclip', label: 'Copy',  icon: 'copy',    type: 'action' },
            ]
        },
        // v1.4.3: Customize toolbar — opens a Customize panel that lets the
        // user show/hide groups + buttons, reorder them, and pick which
        // atom-bar elements appear. Preferences persist to localStorage.
        {
            id: 'custom', label: 'Custom', placement: 'top',
            items: [
                { id: 'customize', label: 'Customize', icon: 'customize', type: 'action' },
            ]
        }
    ];

    var ATOM_BAR = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'H'];

    // Periodic table layout (18-column standard)
    var PTABLE_ROWS = [
        ['H','','','','','','','','','','','','','','','','','He'],
        ['Li','Be','','','','','','','','','','','B','C','N','O','F','Ne'],
        ['Na','Mg','','','','','','','','','','','Al','Si','P','S','Cl','Ar'],
        ['K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'],
        ['Rb','Sr','','','','','','','','','','','In','Sn','Sb','Te','I','Xe'],
        ['Cs','Ba','','','','','','','','','','','','Pb','Bi','','','Rn']
    ];

    // Element colours for the periodic table popup (chemistry-standard CPK-like)
    var PTABLE_COLORS = {
        'H':'#333333','He':'#d9ffff','Li':'#cc80ff','Be':'#c2ff00',
        'B':'#f59e0b','C':'#333333','N':'#2563eb','O':'#dc2626',
        'F':'#16a34a','Ne':'#b3e3f5','Na':'#ab5cf2','Mg':'#8aff00',
        'Al':'#bfa6a6','Si':'#64748b','P':'#ea580c','S':'#ca8a04',
        'Cl':'#16a34a','Ar':'#80d1e3','K':'#8f40d4','Ca':'#3dff00',
        'Sc':'#e6e6e6','Ti':'#bfc2c7','V':'#a6a6ab','Cr':'#8a99c7',
        'Mn':'#9c7ac7','Fe':'#e06633','Co':'#f09050','Ni':'#50d050',
        'Cu':'#c88033','Zn':'#7d80b0','Ga':'#c28f8f','Ge':'#668f8f',
        'As':'#bd93f9','Se':'#ca8a04','Br':'#9333ea','Kr':'#5cb8d1',
        'Rb':'#702eb0','Sr':'#00ff00','In':'#a67573','Sn':'#668080',
        'Sb':'#9e63b5','Te':'#ca8a04','I':'#7c3aed','Xe':'#429eb0',
        'Cs':'#57178f','Ba':'#00c900','Pb':'#575961','Bi':'#9e4fb5',
        'Rn':'#428296'
    };

    // Fragment abbreviations for the Fg popup
    var FRAGMENT_LIST = [
        { label: 'Me', smiles: 'C', name: 'Methyl' },
        { label: 'Et', smiles: 'CC', name: 'Ethyl' },
        { label: 'iPr', smiles: 'CC(C)', name: 'Isopropyl' },
        { label: 'tBu', smiles: 'CC(C)(C)', name: 'tert-Butyl' },
        { label: 'Ph', smiles: 'c1ccccc1', name: 'Phenyl' },
        { label: 'Bn', smiles: 'Cc1ccccc1', name: 'Benzyl' },
        { label: 'Ac', smiles: 'CC(=O)', name: 'Acetyl' },
        { label: 'Boc', smiles: 'CC(C)(C)OC(=O)', name: 'tert-Butoxycarbonyl' },
        { label: 'Cbz', smiles: 'O=C(O)Cc1ccccc1', name: 'Benzyloxycarbonyl' },
        { label: 'Ts', smiles: 'Cc1ccc(S(=O)(=O))cc1', name: 'Tosyl' },
        { label: 'Ms', smiles: 'CS(=O)(=O)', name: 'Mesyl' },
        { label: 'Tf', smiles: 'FC(F)(F)S(=O)(=O)', name: 'Triflyl' },
        { label: 'OMe', smiles: 'CO', name: 'Methoxy' },
        { label: 'OEt', smiles: 'CCO', name: 'Ethoxy' },
        { label: 'NH2', smiles: 'N', name: 'Amino' },
        { label: 'NO2', smiles: '[N+](=O)[O-]', name: 'Nitro' },
        { label: 'COOH', smiles: 'C(=O)O', name: 'Carboxyl' },
        { label: 'CHO', smiles: 'C=O', name: 'Aldehyde' },
        { label: 'CN', smiles: 'C#N', name: 'Cyano' },
        { label: 'SO3H', smiles: 'S(=O)(=O)O', name: 'Sulfonic acid' }
    ];

    // =========================================================================
    // MolEditor constructor
    // =========================================================================
    function MolEditor(containerId, width, height, params) {
        params = params || {};
        this._containerId = containerId;
        this._container = typeof containerId === 'string' ? document.getElementById(containerId) : containerId;
        if (!this._container) throw new Error('BIME: container "' + containerId + '" not found');

        this._width = width || '100%';
        this._height = height || '400px';
        this._options = {};
        this._callbacks = {};
        this._currentTool = null;
        this._currentElement = 'C';
        this._stickyMapNumber = 0;
        this._nameRequestSeq = 0;

        // Data
        this.molecule = new Molecule();
        this.history = new EditorHistory();

        // Parse params
        if (params.options) this._parseOptions(params.options);

        // Build UI
        this._buildUI();

        // Add right-click/context menu on the Name toolbar button for editing
        var nameBtn = this._toolbarButtons['name'];
        if (nameBtn) {
            var self2 = this;
            nameBtn.addEventListener('contextmenu', function(e) {
                e.preventDefault();
                self2._editMolNameInline();
            });
        }

        // Create renderer
        var drawW = this._drawArea.clientWidth || 500;
        var drawH = this._drawArea.clientHeight || 380;
        this.renderer = new Renderer(this._drawArea, drawW, drawH);
        this.renderer.setMolecule(this.molecule);

        // Center the view
        this.renderer.offsetX = drawW / 2;
        this.renderer.offsetY = drawH / 2;

        // Set default tool
        this._setTool('bond');

        // Bind mouse events
        this._bindEvents();

        // Load initial structure if provided
        if (params.smiles) {
            this.readGenericMolecularInput(params.smiles);
        } else if (params.mol) {
            this.readGenericMolecularInput(params.mol);
        } else if (params.jme) {
            this.readGenericMolecularInput(params.jme);
        }

        // Apply options
        if (params.options) this.options(params.options);

        this.render();
    }

    // =========================================================================
    // UI Construction
    // =========================================================================
    MolEditor.prototype._buildUI = function() {
        this._container.innerHTML = '';
        this._container.style.position = 'relative';
        this._container.style.width = this._width;
        this._container.style.height = this._height;
        this._container.style.display = 'flex';
        this._container.style.flexDirection = 'column';
        this._container.style.border = '1px solid var(--color-border, #e2e8f0)';
        this._container.style.borderRadius = '8px';
        this._container.style.overflow = 'hidden';
        this._container.style.fontFamily = '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
        this._container.style.touchAction = 'none';
        this._container.style.webkitTouchCallout = 'none';
        this._container.style.webkitUserSelect = 'none';
        this._container.style.maxWidth = '100%';

        // Inject scoped CSS for responsive behaviour and toolbar styling
        this._injectStyles();

        // v1.5.2: Two-axis toolbar — top horizontal row for high-frequency
        // bond/edit/view actions, plus a left vertical rail for mode-y
        // tools (Draw, Rings, Reaction, Analysis, Export). The atom bar
        // continues to live as a horizontal strip beside the top toolbar
        // so popups (periodic table, fragments) keep their anchor point.

        // ---------- Top toolbar area: top groups + atom bar ----------
        var toolbarArea = document.createElement('div');
        toolbarArea.className = 'bime-toolbar-area';
        toolbarArea.style.cssText = 'display:flex;flex-direction:row;flex-wrap:wrap;align-items:stretch;background:var(--color-surface,#f8fafc);border-bottom:1px solid var(--color-border,#e2e8f0);flex-shrink:0;';

        // Main top tool row (only top-placed groups are routed here)
        this._toolbar = document.createElement('div');
        this._toolbar.className = 'bime-toolbar';
        this._toolbar.setAttribute('role', 'toolbar');
        this._toolbar.setAttribute('aria-label', 'Top toolbar — bonds, edit, view');
        this._toolbar.style.cssText = 'display:flex;align-items:stretch;gap:0;padding:0;flex-wrap:wrap;min-height:32px;flex:1 1 auto;';
        toolbarArea.appendChild(this._toolbar);

        // Atom bar (horizontal, sits inline with toolbar on wide screens; wraps below on narrow)
        this._atomBar = document.createElement('div');
        this._atomBar.className = 'bime-atom-bar';
        this._atomBar.style.cssText = 'display:flex;align-items:center;gap:2px;padding:2px 6px;background:var(--color-bg,white);flex-wrap:wrap;flex:0 1 auto;border-left:1px solid var(--color-border,#e2e8f0);';
        this._buildAtomBar();
        toolbarArea.appendChild(this._atomBar);

        this._container.appendChild(toolbarArea);

        // ---------- Main row: left rail + drawing canvas ----------
        var mainRow = document.createElement('div');
        mainRow.className = 'bime-main-row';
        mainRow.style.cssText = 'display:flex;flex-direction:row;flex:1;min-height:0;';

        // Left vertical rail (icon-only, ~48px wide). Mode-y tools live here.
        this._leftRail = document.createElement('div');
        this._leftRail.className = 'bime-left-rail';
        this._leftRail.setAttribute('role', 'toolbar');
        this._leftRail.setAttribute('aria-orientation', 'vertical');
        this._leftRail.setAttribute('aria-label', 'Side toolbar — drawing tools, rings, reaction, analysis, export');
        this._leftRail.style.cssText = 'display:flex;flex-direction:column;align-items:stretch;width:48px;flex-shrink:0;background:var(--color-surface,#f8fafc);border-right:1px solid var(--color-border,#e2e8f0);overflow-y:auto;overflow-x:visible;padding:2px 0;gap:0;';
        mainRow.appendChild(this._leftRail);

        // Drawing area
        this._drawArea = document.createElement('div');
        this._drawArea.style.cssText = 'flex:1;position:relative;overflow:hidden;background:var(--color-bg,white);min-width:0;';
        mainRow.appendChild(this._drawArea);

        this._container.appendChild(mainRow);

        // Build top + left toolbars (after both target containers exist)
        this._buildToolbar();

        // Status bar
        this._statusBar = document.createElement('div');
        this._statusBar.style.cssText = 'padding:2px 6px;font-size:11px;color:var(--color-text-muted,#64748b);background:var(--color-surface,#f8fafc);border-top:1px solid var(--color-border,#e2e8f0);display:flex;justify-content:space-between;align-items:center;flex-shrink:0;';
        this._statusBar.innerHTML = '<span id="bime-info"></span><span style="color:' + BRAND_COLOR + ';font-weight:600;font-size:10px">BIME \u2014 BioInception Molecular Editor</span>';
        this._container.appendChild(this._statusBar);
    };

    /**
     * Inject scoped styles for responsive toolbar and active states.
     */
    MolEditor.prototype._injectStyles = function() {
        if (document.getElementById('bime-toolbar-styles')) return;
        var style = document.createElement('style');
        style.id = 'bime-toolbar-styles';
        style.textContent = [
            /* Group container */
            '.bime-group{display:flex;align-items:stretch;position:relative}',
            '.bime-toolbar .bime-group+.bime-group{border-left:1px solid var(--color-border,#e2e8f0)}',

            /* Group label — vertical 8px badge that hugs the group */
            '.bime-group-label{display:flex;align-items:center;justify-content:center;writing-mode:vertical-rl;transform:rotate(180deg);padding:2px 0;margin:2px 1px 2px 2px;font-size:8px;font-weight:700;letter-spacing:0.06em;text-transform:uppercase;color:var(--color-text-muted,#94a3b8);pointer-events:none;white-space:nowrap;line-height:1}',

            /* Tool button — compact 30px height */
            '.bime-btn{display:inline-flex;align-items:center;gap:4px;height:30px;min-width:34px;padding:0 6px;border:none;border-radius:0;background:none;cursor:pointer;font-size:11px;font-weight:500;color:var(--color-text,#334155);transition:background 0.12s,color 0.12s;white-space:nowrap;position:relative;outline:none}',
            '.bime-btn:hover{background:var(--color-primary-light,#f0fdfa)}',
            '.bime-btn:focus-visible{box-shadow:inset 0 0 0 2px ' + BRAND_COLOR + '}',
            '.bime-btn.active{background:' + BRAND_COLOR + ';color:white}',
            '.bime-btn.active svg{stroke:white}',
            '.bime-btn svg{flex-shrink:0}',

            /* Label text — hidden on narrow screens */
            '.bime-btn .bime-lbl{display:inline}',

            /* Action buttons (validate, etc) — outlined style, slightly shorter */
            '.bime-btn-action{border:1px solid var(--color-border,#e2e8f0);border-radius:5px;margin:3px 2px;height:24px;font-weight:600;font-size:10px;color:' + BRAND_COLOR + '}',
            '.bime-btn-action:hover{background:' + BRAND_COLOR + ';color:white;border-color:' + BRAND_COLOR + '}',
            '.bime-btn-action:hover svg{stroke:white}',

            /* Atom button — slightly compressed */
            '.bime-atom-btn{display:inline-flex;align-items:center;justify-content:center;min-width:28px;height:24px;padding:0 5px;border:1px solid var(--color-border,#e2e8f0);border-radius:5px;background:var(--color-bg,white);cursor:pointer;font-size:12px;font-weight:700;transition:all 0.12s;outline:none}',
            '.bime-atom-btn:hover{border-color:' + BRAND_COLOR + ';background:var(--color-primary-light,#f0fdfa)}',
            '.bime-atom-btn:focus-visible{box-shadow:0 0 0 2px ' + BRAND_COLOR + '}',
            '.bime-atom-btn.active{background:' + BRAND_COLOR + ';color:white!important;border-color:' + BRAND_COLOR + '}',

            /* Rings flyout (top toolbar — flyout opens DOWN) */
            '.bime-rings-toggle{position:relative}',
            '.bime-rings-flyout{display:none;position:absolute;top:100%;left:0;z-index:100;background:var(--color-surface,#f8fafc);border:1px solid var(--color-border,#e2e8f0);border-radius:6px;box-shadow:0 4px 16px rgba(0,0,0,0.10);padding:4px;flex-direction:row;gap:2px;min-width:max-content}',
            '.bime-rings-toggle.open .bime-rings-flyout{display:flex}',

            /* Group info (i) button — small icon-only at right edge */
            '.bime-group-info{display:inline-flex;align-items:center;justify-content:center;width:14px;height:14px;align-self:center;margin:0 3px 0 1px;padding:0;border:1px solid ' + BRAND_COLOR + ';border-radius:50%;background:transparent;cursor:pointer;font-size:9px;font-weight:700;font-family:inherit;line-height:1;color:' + BRAND_COLOR + ';opacity:0.55;transition:opacity 0.12s,background 0.12s,color 0.12s;outline:none;flex-shrink:0}',
            '.bime-group-info:hover,.bime-group-info[aria-expanded="true"]{opacity:1;background:' + BRAND_COLOR + ';color:white}',
            '.bime-group-info:focus-visible{opacity:1;box-shadow:0 0 0 2px ' + BRAND_COLOR + '}',

            /* Group info popover — teal-bordered light card */
            '.bime-group-popover{position:absolute;top:100%;right:0;z-index:200;margin-top:2px;max-width:min(320px,calc(100vw - 16px));padding:8px 10px;border:1px solid ' + BRAND_COLOR + ';border-radius:6px;background:var(--color-bg,white);box-shadow:0 6px 18px rgba(13,148,136,0.18);font-size:11px;line-height:1.45;color:var(--color-text,#334155);font-weight:400;white-space:normal}',
            '.bime-group-popover::before{content:"";position:absolute;top:-5px;right:6px;width:8px;height:8px;background:var(--color-bg,white);border-left:1px solid ' + BRAND_COLOR + ';border-top:1px solid ' + BRAND_COLOR + ';transform:rotate(45deg)}',

            /* v1.5.2: Left vertical rail — icon-only, ~48px wide, mode-y tools live here */
            '.bime-left-rail .bime-group{flex-direction:column;align-items:stretch;border-left:none;border-top:1px solid transparent;padding:2px 0}',
            '.bime-left-rail .bime-group+.bime-group{border-top:1px solid var(--color-border,#e2e8f0);margin-top:2px;padding-top:4px}',
            '.bime-left-rail .bime-group-label{writing-mode:horizontal-tb;transform:none;font-size:7px;text-align:center;margin:1px 2px 2px;padding:0;color:var(--color-text-muted,#94a3b8);letter-spacing:0.05em}',
            '.bime-left-rail .bime-btn{justify-content:center;padding:0;min-width:0;width:42px;height:34px;margin:1px auto;border-radius:5px;gap:0}',
            '.bime-left-rail .bime-btn .bime-lbl{display:none}',
            '.bime-left-rail .bime-btn-action{width:42px;height:30px;margin:2px auto;padding:0;border-radius:5px;gap:0}',
            /* (i) info button at the bottom of a left-rail group */
            '.bime-left-rail .bime-group-info{align-self:center;margin:2px auto 0}',
            /* Popover from a left-rail group should appear to the right, not below */
            '.bime-left-rail .bime-group-popover{top:auto;right:auto;bottom:0;left:100%;margin-top:0;margin-left:6px;max-width:min(280px,calc(100vw - 80px))}',
            '.bime-left-rail .bime-group-popover::before{top:auto;right:auto;bottom:8px;left:-5px;border-left:1px solid ' + BRAND_COLOR + ';border-top:1px solid ' + BRAND_COLOR + ';border-right:none;border-bottom:none;transform:rotate(-45deg)}',
            /* Rings flyout in left rail — opens to the RIGHT */
            '.bime-left-rail .bime-rings-flyout{top:0;left:100%;margin-left:4px;flex-direction:row}',

            /* Responsive: hide labels below 640px; atom bar wraps to its own row below 900px */
            '@media(max-width:900px){.bime-toolbar{flex:1 1 100%}.bime-atom-bar{flex:1 1 100%;border-left:none!important;border-top:1px solid var(--color-border,#e2e8f0)}}',
            /* On very narrow screens collapse the left rail to icons only (already icon-only)
               and shrink the top toolbar labels. */
            '@media(max-width:640px){.bime-btn .bime-lbl{display:none}.bime-toolbar .bime-group-label{display:none}.bime-btn{padding:0 4px;min-width:30px}.bime-toolbar .bime-group-info{display:none}.bime-left-rail{width:42px}.bime-left-rail .bime-btn,.bime-left-rail .bime-btn-action{width:38px}}',
        ].join('\n');
        document.head.appendChild(style);
    };

    /**
     * Help text shown by the (i) info button on each labelled toolbar group.
     * Keys match TOOLBAR_GROUPS[i].id. Unlabelled "actions" group is skipped.
     */
    var GROUP_INFO = {
        draw:     'Drawing tools — select atoms, draw bonds, lasso, eraser. Click an atom in the periodic table to set the active element.',
        bonds:    'Bond drawing — single, double, triple, aromatic, wedge (up), hash (down). Click an existing bond to cycle order.',
        rings:    'Ring templates — click to drop a 3-/4-/5-/6-/7-membered or aromatic benzene ring at the cursor.',
        edit:     'Editing — Select, Delete, Move atoms or fragments; Undo/Redo recent edits; Clear empties the canvas.',
        view:     'View controls — auto-layout, zoom in/out/fit, toggle CIP R/S and E/Z stereodescriptors, explicit hydrogens, and aromatic display.',
        reaction: 'Reaction tools — Arrow draws a reaction arrow; Map sets atom-atom map numbers manually; Map # toggles their display; Auto-map runs RDT (Rahman 2016) atom-atom mapping; Colors toggles per-atom RDT-style halos; Pairs toggles the molecule-pair overlay.',
        export:   'Export — Name (toggle/edit molecule name), SVG/PNG/Print/Copy save the canvas as image, vector or to clipboard.'
    };

    /**
     * Append a small (i) info button to a labelled toolbar group. Clicking
     * opens a teal-bordered tooltip popover with help text for that group.
     * Escape, outside-click, or clicking another (i) closes it.
     */
    MolEditor.prototype._appendGroupInfo = function(groupEl, group) {
        var self = this;
        var helpText = GROUP_INFO[group.id];
        if (!helpText) return;

        var popoverId = 'bime-group-info-' + group.id;
        var infoBtn = document.createElement('button');
        infoBtn.className = 'bime-group-info';
        infoBtn.type = 'button';
        infoBtn.textContent = 'i';
        infoBtn.setAttribute('aria-label', group.label + ' group help');
        infoBtn.setAttribute('aria-expanded', 'false');
        infoBtn.setAttribute('aria-describedby', popoverId);
        infoBtn.title = group.label + ' help';

        infoBtn.addEventListener('click', function(e) {
            e.stopPropagation();
            // Toggle: close if already open
            var existing = groupEl.querySelector('.bime-group-popover');
            if (existing) {
                self._closeGroupPopover(infoBtn, existing);
                return;
            }
            // Close any other open group popover first
            self._closeAllGroupPopovers();

            var pop = document.createElement('div');
            pop.className = 'bime-group-popover';
            pop.id = popoverId;
            pop.setAttribute('role', 'tooltip');
            pop.textContent = helpText;
            groupEl.appendChild(pop);
            infoBtn.setAttribute('aria-expanded', 'true');

            // Outside-click closes
            var outside = function(ev) {
                if (!pop.contains(ev.target) && ev.target !== infoBtn) {
                    self._closeGroupPopover(infoBtn, pop);
                    document.removeEventListener('click', outside, true);
                    document.removeEventListener('keydown', escClose, true);
                }
            };
            // Escape closes
            var escClose = function(ev) {
                if (ev.key === 'Escape') {
                    self._closeGroupPopover(infoBtn, pop);
                    infoBtn.focus();
                    document.removeEventListener('click', outside, true);
                    document.removeEventListener('keydown', escClose, true);
                }
            };
            setTimeout(function() {
                document.addEventListener('click', outside, true);
                document.addEventListener('keydown', escClose, true);
            }, 0);
        });

        groupEl.appendChild(infoBtn);
    };

    MolEditor.prototype._closeGroupPopover = function(btn, pop) {
        if (pop && pop.parentNode) pop.parentNode.removeChild(pop);
        if (btn) btn.setAttribute('aria-expanded', 'false');
    };

    MolEditor.prototype._closeAllGroupPopovers = function() {
        // v1.5.2: popovers can live in either the top toolbar or the left rail
        var rails = [this._toolbar, this._leftRail];
        for (var r = 0; r < rails.length; r++) {
            if (!rails[r]) continue;
            var pops = rails[r].querySelectorAll('.bime-group-popover');
            for (var i = 0; i < pops.length; i++) {
                var p = pops[i];
                var grp = p.parentNode;
                var btn = grp ? grp.querySelector('.bime-group-info') : null;
                this._closeGroupPopover(btn, p);
            }
        }
    };

    MolEditor.prototype._buildToolbar = function() {
        var self = this;
        this._toolbarButtons = {};
        this._atomBarButtons = {};

        // v1.4.3: apply user-customized toolbar layout from localStorage
        // before building the DOM. ToolbarPrefs.applyToGroups returns a
        // filtered + reordered copy of TOOLBAR_GROUPS; if no prefs exist
        // (or they're invalid) the canonical default is used unchanged.
        var prefs = (typeof ToolbarPrefs !== 'undefined') ? ToolbarPrefs.load() : null;
        var groupsToBuild = (typeof ToolbarPrefs !== 'undefined' && ToolbarPrefs.applyToGroups)
            ? ToolbarPrefs.applyToGroups(TOOLBAR_GROUPS, prefs)
            : TOOLBAR_GROUPS;

        // v1.5.2: route each group to either the top toolbar or the left
        // vertical rail based on its `placement`. Default is 'top' so any
        // user-customised group missing the new property still renders.
        // Look up the canonical placement (prefs strip the property).
        var canonPlacement = {};
        for (var i = 0; i < TOOLBAR_GROUPS.length; i++) {
            canonPlacement[TOOLBAR_GROUPS[i].id] = TOOLBAR_GROUPS[i].placement || 'top';
        }

        groupsToBuild.forEach(function(group) {
            var groupEl = document.createElement('div');
            groupEl.className = 'bime-group';

            var placement = group.placement || canonPlacement[group.id] || 'top';
            groupEl.setAttribute('data-placement', placement);

            // Group label (skip for unlabelled action group)
            if (group.label) {
                var lbl = document.createElement('span');
                lbl.className = 'bime-group-label';
                lbl.textContent = group.label;
                groupEl.appendChild(lbl);
            }

            // Rings group gets a toggle/flyout
            if (group.id === 'rings') {
                self._buildRingsGroup(groupEl, group);
            } else {
                group.items.forEach(function(item) {
                    var isAction = (group.id === 'actions');
                    var btn = self._createToolButton(item, isAction);
                    groupEl.appendChild(btn);
                });
            }

            // (i) info popover at the right edge — only for labelled groups
            if (group.label) {
                self._appendGroupInfo(groupEl, group);
            }

            // Append to the chosen rail. Fall back to top toolbar if the
            // left rail isn't available (e.g. tests building isolated DOM).
            var target = (placement === 'left' && self._leftRail) ? self._leftRail : self._toolbar;
            target.appendChild(groupEl);
        });
    };

    /**
     * Build the Rings group with a toggle button that opens a flyout.
     */
    MolEditor.prototype._buildRingsGroup = function(groupEl, group) {
        var self = this;
        groupEl.classList.add('bime-rings-toggle');

        // Toggle button shows Ring 6 by default
        var toggleBtn = document.createElement('button');
        toggleBtn.className = 'bime-btn';
        toggleBtn.title = 'Rings';
        toggleBtn.innerHTML = self._getToolIcon({ icon: 'ring6' }) +
            '<span class="bime-lbl">Rings</span>' +
            '<svg width="10" height="10" viewBox="0 0 10 10" fill="none" stroke="currentColor" stroke-width="1.5"><polyline points="2,3.5 5,6.5 8,3.5"/></svg>';
        toggleBtn.addEventListener('click', function(e) {
            e.stopPropagation();
            var isOpen = groupEl.classList.toggle('open');
            if (isOpen) {
                // Close on outside click
                var closer = function() { groupEl.classList.remove('open'); document.removeEventListener('click', closer); };
                setTimeout(function() { document.addEventListener('click', closer); }, 0);
            }
        });
        groupEl.appendChild(toggleBtn);
        self._toolbarButtons['_rings_toggle'] = toggleBtn;

        // Flyout panel
        var flyout = document.createElement('div');
        flyout.className = 'bime-rings-flyout';
        group.items.forEach(function(item) {
            var btn = self._createToolButton(item, false);
            btn.addEventListener('click', function() {
                // Update toggle icon to match selected ring
                var oldSvg = toggleBtn.querySelector('svg');
                var temp = document.createElement('div');
                temp.innerHTML = self._getToolIcon(item);
                var newSvg = temp.firstChild;
                if (oldSvg && newSvg) oldSvg.parentNode.replaceChild(newSvg, oldSvg);
                groupEl.classList.remove('open');
            });
            flyout.appendChild(btn);
        });
        groupEl.appendChild(flyout);
    };

    /**
     * Create a single toolbar button (tool, bond, ring, or action).
     */
    MolEditor.prototype._createToolButton = function(item, isActionStyle) {
        var self = this;
        var btn = document.createElement('button');
        btn.className = isActionStyle ? 'bime-btn bime-btn-action' : 'bime-btn';
        btn.type = 'button';
        btn.title = item.label;
        // The label text is also placed in a visible <span class="bime-lbl">,
        // but on small toolbars it can collapse to icon-only. Provide an
        // explicit accessible name and a pressed state for tool buttons.
        btn.setAttribute('aria-label', item.label);
        if (item.type === 'tool' || item.type === 'bond' || item.type === 'ring') {
            btn.setAttribute('aria-pressed', 'false');
        }
        btn.innerHTML = self._getToolIcon(item) + '<span class="bime-lbl">' + item.label + '</span>';

        btn.addEventListener('click', function() {
            if (item.type === 'action') {
                self._doAction(item.id);
            } else if (item.type === 'tool') {
                if (self._currentToolName === item.id && item.id !== 'bond') {
                    self._setTool('bond');
                } else {
                    self._setTool(item.id);
                }
            } else if (item.type === 'bond') {
                self._setTool('bond');
                self._currentTool.bondType = item.bondType;
                self._highlightToolbarButton(item.id);
            } else if (item.type === 'ring') {
                if (self._currentToolName === 'ring' && self._currentRingSize === item.size) {
                    self._setTool('bond');
                } else {
                    self._setTool('ring', item.size);
                    self._currentTool.ringSize = item.size;
                    self._highlightToolbarButton(item.id);
                }
            }
        });

        if (item.disabled) {
            btn.style.opacity = '0.35';
            btn.style.cursor = 'not-allowed';
            btn.title = (btn.title || item.label) + ' (coming soon)';
        }
        self._toolbarButtons[item.id] = btn;
        return btn;
    };

    MolEditor.prototype._buildAtomBar = function() {
        var self = this;
        var colors = Molecule.ELEMENTS;

        // v1.6.1 fix: ensure _atomBarButtons is initialised before
        // _buildAtomBar populates it. Previously _buildToolbar (which sets
        // this map) was called AFTER _buildAtomBar in _buildUI, so the
        // first atom button assignment crashed with "Cannot set properties
        // of undefined". Defensive init here makes the function callable
        // independently of construction order.
        if (!self._atomBarButtons) self._atomBarButtons = {};

        // Leading label
        var lbl = document.createElement('span');
        lbl.style.cssText = 'font-size:9px;font-weight:700;letter-spacing:0.04em;text-transform:uppercase;color:var(--color-text-muted,#94a3b8);margin-right:4px;white-space:nowrap;';
        lbl.textContent = 'Atoms';
        this._atomBar.appendChild(lbl);

        // v1.4.3: apply user-customized atom-bar from localStorage. Falls
        // back to the canonical ATOM_BAR if no prefs are stored.
        var prefs = (typeof ToolbarPrefs !== 'undefined') ? ToolbarPrefs.load() : null;
        var atomsToBuild = (typeof ToolbarPrefs !== 'undefined' && ToolbarPrefs.applyToAtomBar)
            ? ToolbarPrefs.applyToAtomBar(ATOM_BAR, prefs)
            : ATOM_BAR;

        atomsToBuild.forEach(function(sym) {
            var btn = document.createElement('button');
            var elem = colors[sym] || { color: '#333' };
            btn.className = 'bime-atom-btn';
            btn.type = 'button';
            btn.textContent = sym;
            btn.title = (elem.name || sym);
            btn.setAttribute('aria-label', (elem.name ? elem.name + ' (' + sym + ')' : sym));
            btn.setAttribute('aria-pressed', 'false');
            btn.style.color = elem.color;

            btn.addEventListener('click', function() {
                self._currentElement = sym;
                self._setTool('atom');
                self._currentTool.symbol = sym;
                self._highlightAtomButton(sym);
            });

            self._atomBarButtons[sym] = btn;
            self._atomBar.appendChild(btn);
        });

        // Separator
        var sep = document.createElement('span');
        sep.style.cssText = 'width:1px;height:20px;background:var(--color-border,#e2e8f0);margin:0 4px;flex-shrink:0;';
        this._atomBar.appendChild(sep);

        // Periodic Table button
        var ptBtn = document.createElement('button');
        ptBtn.className = 'bime-atom-btn';
        ptBtn.type = 'button';
        ptBtn.style.cssText = 'width:28px;height:28px;display:flex;align-items:center;justify-content:center;';
        ptBtn.title = 'Periodic Table — select any element';
        ptBtn.setAttribute('aria-label', 'Open periodic table to select any element');
        ptBtn.setAttribute('aria-haspopup', 'dialog');
        ptBtn.innerHTML = '<svg aria-hidden="true" focusable="false" viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="' + BRAND_COLOR + '" stroke-width="1.5"><rect x="2" y="2" width="20" height="20" rx="3"/><rect x="5" y="5" width="4" height="4" rx="1" fill="' + BRAND_COLOR + '" stroke="none"/><rect x="10" y="5" width="4" height="4" rx="1" fill="' + BRAND_COLOR + '" stroke="none" opacity="0.6"/><rect x="15" y="5" width="4" height="4" rx="1" fill="' + BRAND_COLOR + '" stroke="none" opacity="0.3"/><rect x="5" y="10" width="4" height="4" rx="1" fill="' + BRAND_COLOR + '" stroke="none" opacity="0.6"/><rect x="5" y="15" width="4" height="4" rx="1" fill="' + BRAND_COLOR + '" stroke="none" opacity="0.3"/></svg>';
        ptBtn.addEventListener('click', function() { self._showPeriodicTable(); });
        this._atomBar.appendChild(ptBtn);

        // Fragment Abbreviations button
        var fgBtn = document.createElement('button');
        fgBtn.className = 'bime-atom-btn';
        fgBtn.type = 'button';
        fgBtn.style.cssText = 'width:28px;height:28px;display:flex;align-items:center;justify-content:center;';
        fgBtn.title = 'Fragment abbreviations — Ph, Me, Et, Boc, etc.';
        fgBtn.setAttribute('aria-label', 'Open fragment abbreviations menu (Ph, Me, Et, Boc and others)');
        fgBtn.setAttribute('aria-haspopup', 'menu');
        fgBtn.innerHTML = '<svg aria-hidden="true" focusable="false" viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="' + BRAND_COLOR + '" stroke-width="1.5"><circle cx="8" cy="8" r="4"/><circle cx="16" cy="8" r="4"/><line x1="12" y1="8" x2="12" y2="8"/><circle cx="12" cy="16" r="4"/><line x1="8" y1="12" x2="12" y2="12"/><line x1="16" y1="12" x2="12" y2="16" stroke-dasharray="2,1"/></svg>';
        fgBtn.addEventListener('click', function() { self._showFragmentMenu(); });
        this._atomBar.appendChild(fgBtn);
    };

    /**
     * Highlight the active atom button, clear others.
     */
    MolEditor.prototype._highlightAtomButton = function(activeSym) {
        for (var sym in this._atomBarButtons) {
            var btn = this._atomBarButtons[sym];
            var isActive = (sym === activeSym && this._currentToolName === 'atom');
            if (isActive) {
                btn.classList.add('active');
            } else {
                btn.classList.remove('active');
            }
            // Mirror visual state to assistive tech (WCAG 4.1.2 Name/Role/Value)
            btn.setAttribute('aria-pressed', isActive ? 'true' : 'false');
        }
    };

    // =========================================================================
    // Periodic Table Popup
    // =========================================================================

    /**
     * Remove any open popup (periodic table or fragment menu).
     */
    MolEditor.prototype._closePopup = function(id) {
        var old = this._container.querySelector('#' + id);
        if (old) old.remove();
    };

    /**
     * Create shared popup CSS. Injected once.
     */
    MolEditor.prototype._injectPopupStyles = function() {
        if (document.getElementById('bime-popup-styles')) return;
        var style = document.createElement('style');
        style.id = 'bime-popup-styles';
        style.textContent = [
            '.bime-popup{position:absolute;z-index:200;background:var(--color-bg,white);border:1px solid var(--color-border,#e2e8f0);border-radius:8px;box-shadow:0 8px 24px rgba(0,0,0,0.15);padding:8px;max-height:360px;overflow-y:auto}',
            '.bime-popup-title{font-size:10px;font-weight:700;text-transform:uppercase;letter-spacing:0.04em;color:var(--color-text-muted,#94a3b8);margin-bottom:6px;padding:0 2px}',
            '.bime-pt-grid{display:grid;grid-template-columns:repeat(18,1fr);gap:1px}',
            '.bime-pt-cell{display:flex;align-items:center;justify-content:center;width:26px;height:24px;border:1px solid transparent;border-radius:3px;cursor:pointer;font-size:10px;font-weight:700;transition:all 0.1s;background:var(--color-surface,#f8fafc)}',
            '.bime-pt-cell:hover{border-color:' + BRAND_COLOR + ';background:var(--color-primary-light,#f0fdfa);transform:scale(1.15)}',
            '.bime-pt-cell.empty{visibility:hidden;cursor:default}',
            '.bime-pt-cell.empty:hover{transform:none;border-color:transparent}',
            '.bime-frag-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:3px}',
            '.bime-frag-btn{display:flex;flex-direction:column;align-items:center;justify-content:center;padding:6px 4px;border:1px solid var(--color-border,#e2e8f0);border-radius:5px;cursor:pointer;transition:all 0.1s;background:var(--color-surface,#f8fafc)}',
            '.bime-frag-btn:hover{border-color:' + BRAND_COLOR + ';background:var(--color-primary-light,#f0fdfa)}',
            '.bime-frag-btn .frag-abbr{font-size:12px;font-weight:700;color:' + BRAND_COLOR + '}',
            '.bime-frag-btn .frag-name{font-size:8px;color:var(--color-text-muted,#94a3b8);margin-top:1px}'
        ].join('\n');
        document.head.appendChild(style);
    };

    /**
     * Show a periodic table popup near the atom bar.
     * Clicking an element sets the current drawing tool to AtomTool with that element.
     */
    MolEditor.prototype._showPeriodicTable = function() {
        var self = this;
        this._injectPopupStyles();

        // Toggle off if already open
        var old = this._container.querySelector('#bime-ptable');
        if (old) { old.remove(); return; }

        // Close any other popup
        this._closePopup('bime-frags');

        var popup = document.createElement('div');
        popup.id = 'bime-ptable';
        popup.className = 'bime-popup';
        popup.style.left = '8px';
        popup.style.top = (this._atomBar.offsetTop + this._atomBar.offsetHeight + 2) + 'px';

        var title = document.createElement('div');
        title.className = 'bime-popup-title';
        title.textContent = 'Periodic Table';
        popup.appendChild(title);

        var grid = document.createElement('div');
        grid.className = 'bime-pt-grid';

        for (var r = 0; r < PTABLE_ROWS.length; r++) {
            var row = PTABLE_ROWS[r];
            for (var c = 0; c < row.length; c++) {
                var sym = row[c];
                var cell = document.createElement('div');
                cell.className = 'bime-pt-cell';
                if (!sym) {
                    cell.classList.add('empty');
                } else {
                    cell.textContent = sym;
                    cell.title = sym;
                    cell.style.color = PTABLE_COLORS[sym] || '#333';
                    (function(s) {
                        cell.addEventListener('click', function(e) {
                            e.stopPropagation();
                            self._currentElement = s;
                            self._setTool('atom');
                            if (self._currentTool) self._currentTool.symbol = s;
                            self._highlightAtomButton(s);
                            popup.remove();
                            self.showInfo('Element: ' + s);
                        });
                    })(sym);
                }
                grid.appendChild(cell);
            }
        }

        popup.appendChild(grid);
        this._container.appendChild(popup);

        // Close on click outside or Escape
        var closeHandler = function(e) {
            if (!popup.contains(e.target)) {
                popup.remove();
                document.removeEventListener('mousedown', closeHandler);
                document.removeEventListener('keydown', escHandler);
            }
        };
        var escHandler = function(e) {
            if (e.key === 'Escape') {
                popup.remove();
                document.removeEventListener('mousedown', closeHandler);
                document.removeEventListener('keydown', escHandler);
            }
        };
        setTimeout(function() {
            document.addEventListener('mousedown', closeHandler);
            document.addEventListener('keydown', escHandler);
        }, 0);
    };

    // =========================================================================
    // Fragment Abbreviations Popup
    // =========================================================================

    /**
     * Show a popup with common functional group fragments.
     * Clicking a fragment loads its SMILES into the editor.
     */
    MolEditor.prototype._showFragmentMenu = function() {
        var self = this;
        this._injectPopupStyles();

        // Toggle off if already open
        var old = this._container.querySelector('#bime-frags');
        if (old) { old.remove(); return; }

        // Close any other popup
        this._closePopup('bime-ptable');

        var popup = document.createElement('div');
        popup.id = 'bime-frags';
        popup.className = 'bime-popup';
        popup.style.left = '8px';
        popup.style.top = (this._atomBar.offsetTop + this._atomBar.offsetHeight + 2) + 'px';
        popup.style.minWidth = '260px';

        var title = document.createElement('div');
        title.className = 'bime-popup-title';
        title.textContent = 'Fragment Abbreviations';
        popup.appendChild(title);

        var grid = document.createElement('div');
        grid.className = 'bime-frag-grid';

        for (var i = 0; i < FRAGMENT_LIST.length; i++) {
            var frag = FRAGMENT_LIST[i];
            var btn = document.createElement('div');
            btn.className = 'bime-frag-btn';
            btn.title = frag.name + ' (' + frag.smiles + ')';

            var abbr = document.createElement('span');
            abbr.className = 'frag-abbr';
            abbr.textContent = frag.label;
            btn.appendChild(abbr);

            var name = document.createElement('span');
            name.className = 'frag-name';
            name.textContent = frag.name;
            btn.appendChild(name);

            (function(f) {
                btn.addEventListener('click', function(e) {
                    e.stopPropagation();
                    self.saveHistory();
                    self.readGenericMolecularInput(f.smiles);
                    self.showInfo('Fragment: ' + f.name);
                    popup.remove();
                });
            })(frag);

            grid.appendChild(btn);
        }

        popup.appendChild(grid);
        this._container.appendChild(popup);

        // Close on click outside or Escape
        var closeHandler = function(e) {
            if (!popup.contains(e.target)) {
                popup.remove();
                document.removeEventListener('mousedown', closeHandler);
                document.removeEventListener('keydown', escHandler);
            }
        };
        var escHandler = function(e) {
            if (e.key === 'Escape') {
                popup.remove();
                document.removeEventListener('mousedown', closeHandler);
                document.removeEventListener('keydown', escHandler);
            }
        };
        setTimeout(function() {
            document.addEventListener('mousedown', closeHandler);
            document.addEventListener('keydown', escHandler);
        }, 0);
    };

    MolEditor.prototype._getToolIcon = function(item) {
        var icons = {
            bond:   '<line x1="6" y1="18" x2="18" y2="6"/>',
            single: '<line x1="4" y1="12" x2="20" y2="12"/>',
            double: '<line x1="4" y1="10" x2="20" y2="10"/><line x1="4" y1="14" x2="20" y2="14"/>',
            triple: '<line x1="4" y1="8" x2="20" y2="8"/><line x1="4" y1="12" x2="20" y2="12"/><line x1="4" y1="16" x2="20" y2="16"/>',
            ring3:  '<polygon points="12,5 4,19 20,19" fill="none"/>',
            ring4:  '<rect x="5" y="5" width="14" height="14" fill="none"/>',
            ring5:  '<polygon points="12,4 3,11 6,21 18,21 21,11" fill="none"/>',
            ring6:  '<polygon points="12,4 4,8 4,16 12,20 20,16 20,8" fill="none"/>',
            ring7:  '<polygon points="12,3 5,6 3,13 7,20 17,20 21,13 19,6" fill="none"/>',
            chain:  '<polyline points="4,18 8,8 14,18 20,8"/><circle cx="4" cy="18" r="1.5" fill="currentColor" stroke="none"/><circle cx="20" cy="8" r="1.5" fill="currentColor" stroke="none"/>',
            stereo: '<polygon points="4,12 20,6 20,18" fill="' + BRAND_COLOR + '" fill-opacity="0.25" stroke="' + BRAND_COLOR + '"/>',
            select: '<rect x="4" y="4" width="16" height="16" rx="2" fill="none" stroke-dasharray="3,2"/><polyline points="7,12 10,16 17,7" stroke-width="2"/>',
            del:    '<polyline points="3,6 5,6 21,6"/><path d="M19 6v14a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2V6m3 0V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v2"/>',
            movea:  '<path d="M5 9l4-4 4 4"/><path d="M9 5v14"/><path d="M19 15l-4 4-4-4"/><path d="M15 19V5"/>',
            undo:   '<polyline points="1 4 1 10 7 10"/><path d="M3.51 15a9 9 0 1 0 2.13-9.36L1 10"/>',
            redo:   '<polyline points="23 4 23 10 17 10"/><path d="M20.49 15a9 9 0 1 1-2.12-9.36L23 10"/>',
            clean:  '<rect x="3" y="3" width="18" height="18" rx="2"/><line x1="9" y1="9" x2="15" y2="15"/><line x1="15" y1="9" x2="9" y2="15"/>',
            react:  '<line x1="5" y1="12" x2="19" y2="12"/><polyline points="15 8 19 12 15 16"/>',
            '123':  '<path d="M4 15s1-1 4-1 5 2 8 2 4-1 4-1V3s-1 1-4 1-5-2-8-2-4 1-4 1z"/><line x1="4" y1="22" x2="4" y2="15"/>',
            validate: '<path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/><polyline points="22 4 12 14.01 9 11.01"/>',
            smarts:   '<circle cx="12" cy="12" r="8" stroke-dasharray="3,2"/><text x="12" y="15" font-size="8" fill="currentColor" stroke="none" text-anchor="middle" font-weight="bold">S</text>',
            mcs:      '<circle cx="8" cy="10" r="6"/><circle cx="16" cy="14" r="6"/><path d="M12 7v10" stroke-dasharray="2,1" opacity="0.5"/><circle cx="12" cy="12" r="3" fill="currentColor" opacity="0.25" stroke="none"/>',
            simsearch:'<circle cx="10.5" cy="10.5" r="7" stroke-width="1.6"/><line x1="15.8" y1="15.8" x2="21" y2="21" stroke-width="2.2"/><line x1="8" y1="7" x2="12" y2="7" stroke-width="1.2"/><line x1="12" y1="7" x2="13.5" y2="10.5" stroke-width="1.2"/><line x1="13.5" y1="10.5" x2="10.5" y2="13" stroke-width="1.2"/><line x1="10.5" y1="13" x2="7" y2="10.5" stroke-width="1.2"/><line x1="7" y1="10.5" x2="8" y2="7" stroke-width="1.2"/><circle cx="8" cy="7" r="1.2" fill="currentColor" stroke="none"/><circle cx="12" cy="7" r="1.2" fill="currentColor" stroke="none"/><circle cx="13.5" cy="10.5" r="1.2" fill="currentColor" stroke="none"/><circle cx="10.5" cy="13" r="1.2" fill="currentColor" stroke="none"/><circle cx="7" cy="10.5" r="1.2" fill="currentColor" stroke="none"/>',
            autolayout: '<rect x="3" y="3" width="18" height="18" rx="2"/><line x1="7" y1="8" x2="17" y2="8"/><line x1="7" y1="12" x2="14" y2="12"/><line x1="7" y1="16" x2="17" y2="16"/>',
            zoomin:  '<circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/><line x1="11" y1="8" x2="11" y2="14"/><line x1="8" y1="11" x2="14" y2="11"/>',
            zoomout: '<circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/><line x1="8" y1="11" x2="14" y2="11"/>',
            zoomfit: '<polyline points="15 3 21 3 21 9"/><polyline points="9 21 3 21 3 15"/><polyline points="21 15 21 21 15 21"/><polyline points="3 9 3 3 9 3"/>',
            ciprs: '<text x="3" y="14" font-size="9" fill="currentColor" font-weight="bold" font-style="italic">R</text><text x="11" y="14" font-size="7" fill="currentColor">/</text><text x="14" y="14" font-size="9" fill="currentColor" font-weight="bold" font-style="italic">S</text>',
            cipez: '<text x="3" y="14" font-size="9" fill="currentColor" font-weight="bold" font-style="italic">E</text><text x="11" y="14" font-size="7" fill="currentColor">/</text><text x="14" y="14" font-size="9" fill="currentColor" font-weight="bold" font-style="italic">Z</text>',
            toggleh: '<text x="7" y="16" font-size="14" fill="currentColor" font-weight="bold">H</text>',
            togglearo: '<circle cx="12" cy="12" r="6" stroke-dasharray="3,2"/><circle cx="12" cy="12" r="2" fill="currentColor" stroke="none"/>',
            name: '<text x="12" y="16" font-size="11" fill="currentColor" stroke="none" text-anchor="middle" font-weight="bold">Aa</text>',
            svg: '<rect x="3" y="3" width="18" height="18" rx="2"/><text x="12" y="15" font-size="8" fill="currentColor" stroke="none" text-anchor="middle" font-weight="bold">SVG</text>',
            png: '<rect x="3" y="3" width="18" height="18" rx="2"/><text x="12" y="15" font-size="8" fill="currentColor" stroke="none" text-anchor="middle" font-weight="bold">PNG</text>',
            print: '<polyline points="6 9 6 2 18 2 18 9"/><path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"/><rect x="6" y="14" width="12" height="8"/>',
            copy: '<rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/>',
            ptable: '<text x="4" y="10" font-size="7" fill="currentColor" font-weight="bold" stroke="none">Pt</text><rect x="2" y="12" width="20" height="10" rx="2" stroke="currentColor" fill="none" stroke-width="1"/>',
            fg: '<text x="12" y="10" font-size="7" fill="currentColor" font-weight="bold" stroke="none" text-anchor="middle">Fg</text><rect x="3" y="12" width="8" height="6" rx="1" stroke="currentColor" fill="none" stroke-width="0.8"/><rect x="13" y="12" width="8" height="6" rx="1" stroke="currentColor" fill="none" stroke-width="0.8"/><rect x="3" y="19" width="8" height="3" rx="1" stroke="currentColor" fill="none" stroke-width="0.8"/>',
            rdtautomap: '<line x1="4" y1="12" x2="20" y2="12"/><polyline points="16 8 20 12 16 16"/><circle cx="6" cy="6" r="1.6" fill="currentColor" stroke="none"/><circle cx="18" cy="18" r="1.6" fill="currentColor" stroke="none"/><path d="M6 6 L9 9" stroke-dasharray="2,1.4" opacity="0.7"/><path d="M18 18 L15 15" stroke-dasharray="2,1.4" opacity="0.7"/>',
            togglecolors: '<circle cx="7" cy="9" r="3.5" fill="#a5d8d2" stroke="none"/><circle cx="14" cy="9" r="3.5" fill="#fbd5a5" stroke="none"/><circle cx="10.5" cy="15" r="3.5" fill="#c7e8b9" stroke="none"/>',
            togglepairs: '<rect x="3" y="6" width="7" height="12" rx="1" fill="none"/><rect x="14" y="6" width="7" height="12" rx="1" fill="none"/><line x1="10" y1="12" x2="14" y2="12" stroke-dasharray="2,1.5"/>',
            // v1.4.3: customize icon — sliders / preferences glyph
            customize: '<line x1="4" y1="6" x2="14" y2="6"/><circle cx="16" cy="6" r="2" fill="white"/><line x1="4" y1="12" x2="9" y2="12"/><circle cx="11" cy="12" r="2" fill="white"/><line x1="4" y1="18" x2="16" y2="18"/><circle cx="18" cy="18" r="2" fill="white"/>'
        };
        var svg = icons[item.icon] || icons[item.id] || '';
        return '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="' + BRAND_COLOR + '" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">' + svg + '</svg>';
    };

    MolEditor.prototype._highlightToolbarButton = function(activeId) {
        // Clear atom bar selection when a non-atom tool is chosen
        if (activeId !== 'atom') {
            for (var sym in this._atomBarButtons) {
                var ab = this._atomBarButtons[sym];
                ab.classList.remove('active');
                if (ab.hasAttribute('aria-pressed')) ab.setAttribute('aria-pressed', 'false');
            }
        }

        for (var id in this._toolbarButtons) {
            var btn = this._toolbarButtons[id];
            if (id === '_rings_toggle') continue; // managed separately
            var pressed = (id === activeId);
            if (pressed) {
                btn.classList.add('active');
            } else {
                btn.classList.remove('active');
            }
            // Reflect state for assistive tech on toggle buttons
            if (btn.hasAttribute('aria-pressed')) {
                btn.setAttribute('aria-pressed', pressed ? 'true' : 'false');
            }
        }

        // If a ring was selected, also highlight the rings toggle
        var ringsToggle = this._toolbarButtons['_rings_toggle'];
        if (ringsToggle) {
            var isRing = (activeId && activeId.indexOf('ring') === 0);
            if (isRing) {
                ringsToggle.classList.add('active');
            } else {
                ringsToggle.classList.remove('active');
            }
        }
    };

    // =========================================================================
    // Tool management
    // =========================================================================
    MolEditor.prototype._setTool = function(toolName, arg) {
        if (this._currentTool) this._currentTool.deactivate();
        this._currentToolName = toolName;
        if (toolName === 'ring') { this._currentRingSize = arg || 6; }

        switch (toolName) {
            case 'atom':    this._currentTool = new EditorTools.AtomTool(this, this._currentElement); break;
            case 'bond':    this._currentTool = new EditorTools.BondTool(this); break;
            case 'select':  this._currentTool = new EditorTools.SelectTool(this); break;
            case 'delete':  this._currentTool = new EditorTools.DeleteTool(this); break;
            case 'move':    this._currentTool = new EditorTools.MoveTool(this); break;
            case 'ring':    this._currentTool = new EditorTools.RingTool(this, arg || 6); break;
            case 'chain':   this._currentTool = new EditorTools.ChainTool(this); break;
            case 'stereo':  this._currentTool = new EditorTools.StereoBondTool(this); break;
            case 'reaction':this._currentTool = new EditorTools.ReactionTool(this); break;
            case 'mapping': this._currentTool = new EditorTools.MapTool(this); break;
            default:        this._currentTool = new EditorTools.BondTool(this); break;
        }
        this._currentTool.activate();
        this._highlightToolbarButton(toolName);
    };

    MolEditor.prototype._doAction = function(action) {
        switch (action) {
            case 'undo':
                var state = this.history.undo(this.molecule);
                if (state) { this._restoreMolecule(state); this.render(); this._updateNameStatus(); this._fireCallback('AfterStructureModified'); }
                break;
            case 'redo':
                state = this.history.redo(this.molecule);
                if (state) { this._restoreMolecule(state); this.render(); this._updateNameStatus(); this._fireCallback('AfterStructureModified'); }
                break;
            case 'clear':
                this.saveHistory();
                this.molecule.clear();
                // v1.4.0: dropping the molecule also drops any pair-highlight
                // boxes the renderer may be holding from a previous Auto-map.
                this.clearComponentPairs();
                this.changed();
                break;
            case 'validate':
                this._promptValidateSmiles();
                break;
            case 'smarts':
                this._promptSmartsSearch();
                break;
            case 'mcs':
                this._promptMCS();
                break;
            case 'togglemap':
                var on = this.toggleMapNumbers();
                var tmBtn = this._toolbarButtons['togglemap'];
                if (tmBtn) tmBtn.style.opacity = on ? '1' : '0.4';
                this.showInfo(on ? 'Mapping numbers shown' : 'Mapping numbers hidden');
                break;
            case 'rdtautomap':
                this._runRdtAutoMap();
                break;
            case 'togglepairs':
                var ppOn = this.toggleComponentPairs();
                var ppBtn = this._toolbarButtons['togglepairs'];
                if (ppBtn) ppBtn.style.opacity = ppOn ? '1' : '0.4';
                this.showInfo(ppOn ? 'Pair highlight shown' : 'Pair highlight hidden');
                break;
            case 'togglecolors':
                var caOn = this.toggleColorAtoms();
                var caBtn = this._toolbarButtons['togglecolors'];
                if (caBtn) caBtn.style.opacity = caOn ? '1' : '0.4';
                this.showInfo(caOn ? 'Atom colours shown' : 'Atom colours hidden');
                break;
            case 'simsearch':
                this._promptSimSearch();
                break;
            case 'ciprs':
                this.renderer.showCipRS = !this.renderer.showCipRS;
                var rsBtn = this._toolbarButtons['ciprs'];
                if (rsBtn) rsBtn.style.opacity = this.renderer.showCipRS ? '1' : '0.4';
                this.render();
                this.showInfo(this.renderer.showCipRS ? 'R/S stereo labels shown' : 'R/S labels hidden');
                break;
            case 'cipez':
                this.renderer.showCipEZ = !this.renderer.showCipEZ;
                var ezBtn = this._toolbarButtons['cipez'];
                if (ezBtn) ezBtn.style.opacity = this.renderer.showCipEZ ? '1' : '0.4';
                this.render();
                this.showInfo(this.renderer.showCipEZ ? 'E/Z stereo labels shown' : 'E/Z labels hidden');
                break;
            case 'toggleh':
                this.renderer.showHydrogens = !this.renderer.showHydrogens;
                var hBtn = this._toolbarButtons['toggleh'];
                if (hBtn) hBtn.style.opacity = this.renderer.showHydrogens ? '1' : '0.4';
                this.render();
                this.showInfo(this.renderer.showHydrogens ? 'Hydrogens shown' : 'Hydrogens hidden');
                break;
            case 'togglearo':
                this._aromaticStyle = this._aromaticStyle === 'circle' ? 'kekule' : 'circle';
                this.renderer.aromaticStyle = this._aromaticStyle;
                var aroBtn = this._toolbarButtons['togglearo'];
                if (aroBtn) aroBtn.style.opacity = this._aromaticStyle === 'circle' ? '1' : '0.4';
                this.render();
                this.showInfo(this._aromaticStyle === 'circle' ? 'Aromatic circles shown' : 'Kekul\u00e9 bonds shown');
                break;
            case 'autolayout':
                this._autoLayout();
                break;
            case 'zoomin':
                this.renderer.scale = Math.min(5, this.renderer.scale * 1.25);
                this.render();
                break;
            case 'zoomout':
                this.renderer.scale = Math.max(0.2, this.renderer.scale * 0.8);
                this.render();
                break;
            case 'zoomfit':
                this._zoomToFit();
                break;
            case 'name':
                this._toggleOrEditMolName();
                break;
            case 'exportsvg':
                this._exportSVGWithImageExport();
                break;
            case 'exportpng':
                this._exportPNGWithImageExport();
                break;
            case 'exportprint':
                this._exportPrintSVG();
                break;
            case 'copyclip':
                this._copyToClipboard();
                break;
            case 'customize':
                this._openCustomizePanel();
                break;
        }
    };

    // =========================================================================
    // Mouse event handling
    // =========================================================================
    MolEditor.prototype._bindEvents = function() {
        var self = this;
        var svg = this.renderer.svg;

        // Unified pointer handler — works with mouse and touch
        function getPos(e) {
            var cx, cy;
            if (e.touches && e.touches.length > 0) {
                cx = e.touches[0].clientX;
                cy = e.touches[0].clientY;
            } else if (e.changedTouches && e.changedTouches.length > 0) {
                cx = e.changedTouches[0].clientX;
                cy = e.changedTouches[0].clientY;
            } else {
                cx = e.clientX;
                cy = e.clientY;
            }
            return self.renderer.screenToMol(cx, cy);
        }

        function onDown(e) {
            if (self._options.depict) return;
            if (e.touches) e.preventDefault(); // prevent scroll on touch
            var pos = getPos(e);
            if (self._currentTool) self._currentTool.onMouseDown(e, self.molecule, pos);
        }

        function onMove(e) {
            if (self._options.depict) return;
            var pos = getPos(e);

            // Hover highlighting (mouse only — skip on touch to avoid jank)
            if (!e.touches) {
                // Preserve sticky map highlight state
                var stickyMap = self._stickyMapNumber || 0;
                self._clearHighlights();

                // Restore sticky map highlights if set
                if (stickyMap > 0) {
                    self._stickyMapNumber = stickyMap;
                    self._highlightMapPartners(stickyMap);
                }

                var atom = self.molecule.getAtomAt(pos.x, pos.y);
                if (atom) {
                    atom.highlighted = true;
                    self._fireCallback('AtomHighlight', { atom: atom.id });

                    // Cross-highlight mapped partners on hover
                    if (atom.mapNumber > 0 && atom.mapNumber !== stickyMap) {
                        var ids = self._highlightMapPartners(atom.mapNumber);
                        self._fireCallback('MapHighlight', { mapNumber: atom.mapNumber, atoms: ids });
                        self.showInfo('Map :' + atom.mapNumber);
                    } else if (atom.mapNumber > 0 && atom.mapNumber === stickyMap) {
                        self.showInfo('Map :' + atom.mapNumber + ' (pinned)');
                    }
                } else {
                    var bond = self.molecule.getBondAt(pos.x, pos.y);
                    if (bond) {
                        bond.highlighted = true;
                        self._fireCallback('BondHighlight', { bond: bond.id });
                    }
                    if (!stickyMap) self.showInfo('');
                }
                if (!self._renderPending) {
                    self._renderPending = true;
                    requestAnimationFrame(function() {
                        self._renderPending = false;
                        self.render();
                    });
                }
            }

            if (self._currentTool) self._currentTool.onMouseMove(e, self.molecule, pos);
        }

        function onUp(e) {
            if (self._options.depict) return;
            var pos = getPos(e);
            if (self._currentTool) self._currentTool.onMouseUp(e, self.molecule, pos);
        }

        // v1.4.1 (bug-fix #4): every SVG handler is now named and stored on
        // self._svgHandlers so destroy() can remove all of them. Previous
        // versions only tracked onDown/onMove/onUp, leaving 8 anonymous
        // listeners (click, dblclick, two touchstarts, two touchmoves,
        // touchend, wheel) bound across the SVG and one ResizeObserver
        // callback. Long-running apps that rebuild editors leak each one.
        var onClick = function(e) {
            if (self._options.depict) return;
            var pos = getPos(e);
            var atom = self.molecule.getAtomAt(pos.x, pos.y);
            if (atom && atom.mapNumber > 0) {
                if (self._stickyMapNumber === atom.mapNumber) {
                    // Same mapping clicked again — unpin
                    self._clearMapHighlights();
                    self.showInfo('');
                } else {
                    // Pin this mapping
                    self._clearMapHighlights();
                    self._stickyMapNumber = atom.mapNumber;
                    var ids = self._highlightMapPartners(atom.mapNumber);
                    self._fireCallback('MapHighlight', { mapNumber: atom.mapNumber, atoms: ids, sticky: true });
                    self.showInfo('Map :' + atom.mapNumber + ' (pinned)');
                }
                self.render();
            } else if (!atom) {
                // Clicked empty space — clear sticky map highlight
                if (self._stickyMapNumber) {
                    self._clearMapHighlights();
                    self.showInfo('');
                    self.render();
                }
            }
        };
        var onDblClick = function(e) {
            if (self._options.depict) return;
            var nameBounds = self.renderer.getMolNameBounds();
            if (!nameBounds) return;
            var rect = svg.getBoundingClientRect();
            var vx = (e.clientX - rect.left) / rect.width * self.renderer.width;
            var vy = (e.clientY - rect.top) / rect.height * self.renderer.height;
            if (vx >= nameBounds.x && vx <= nameBounds.x + nameBounds.w &&
                vy >= nameBounds.y && vy <= nameBounds.y + nameBounds.h) {
                e.preventDefault();
                self._editMolNameInline();
            }
        };
        var lastPinchDist = 0;
        var onPinchStart = function(e) {
            if (e.touches.length === 2) {
                var dx = e.touches[0].clientX - e.touches[1].clientX;
                var dy = e.touches[0].clientY - e.touches[1].clientY;
                lastPinchDist = Math.sqrt(dx * dx + dy * dy);
            }
        };
        var onPinchMove = function(e) {
            if (e.touches.length === 2) {
                e.preventDefault();
                var dx = e.touches[0].clientX - e.touches[1].clientX;
                var dy = e.touches[0].clientY - e.touches[1].clientY;
                var dist = Math.sqrt(dx * dx + dy * dy);
                if (lastPinchDist > 0) {
                    var factor = dist / lastPinchDist;
                    self.renderer.scale = Math.max(0.3, Math.min(5, self.renderer.scale * factor));
                    self.render();
                }
                lastPinchDist = dist;
            }
        };
        var onWheel = function(e) {
            e.preventDefault();
            var factor = e.deltaY > 0 ? 0.9 : 1.1;
            self.renderer.scale = Math.max(0.3, Math.min(5, self.renderer.scale * factor));
            self.render();
        };

        // Store references for destroy(). All handlers go in one bag so the
        // teardown path is data-driven and additions stay symmetric.
        self._svgHandlers = {
            onDown: onDown, onMove: onMove, onUp: onUp,
            onClick: onClick, onDblClick: onDblClick,
            onPinchStart: onPinchStart, onPinchMove: onPinchMove,
            onWheel: onWheel
        };

        // Mouse events
        svg.addEventListener('mousedown', onDown);
        svg.addEventListener('mousemove', onMove);
        svg.addEventListener('mouseup', onUp);

        svg.addEventListener('click', onClick);
        svg.addEventListener('dblclick', onDblClick);

        // Touch events (mobile/tablet) — the SAME onDown/onMove/onUp bag also
        // serves as the touch handler. Attaching twice is intentional (one
        // listener per event type), but destroy() removes both attachments.
        svg.addEventListener('touchstart', onDown, { passive: false });
        svg.addEventListener('touchmove', onMove, { passive: false });
        svg.addEventListener('touchend', onUp);
        svg.addEventListener('touchstart', onPinchStart, { passive: false });
        svg.addEventListener('touchmove', onPinchMove, { passive: false });
        svg.addEventListener('wheel', onWheel, { passive: false });

        // Keyboard shortcuts
        // FIX: store reference so it can be removed if editor is destroyed (prevents memory leak)
        self._keydownHandler = function(e) {
            if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;
            if (e.ctrlKey || e.metaKey) {
                if (e.key === 'z') { e.preventDefault(); self._doAction(e.shiftKey ? 'redo' : 'undo'); }
                else if (e.key === 'y') { e.preventDefault(); self._doAction('redo'); }
            } else if (e.key === 'Delete' || e.key === 'Backspace') {
                if (self._currentToolName === 'select') {
                    // Delete selected atoms
                    var toRemove = [];
                    for (var i = 0; i < self.molecule.atoms.length; i++) {
                        if (self.molecule.atoms[i].selected) toRemove.push(self.molecule.atoms[i].id);
                    }
                    if (toRemove.length > 0) {
                        self.saveHistory();
                        for (var j = 0; j < toRemove.length; j++) self.molecule.removeAtom(toRemove[j]);
                        self.render();
                        self.showInfo(toRemove.length + ' atoms deleted');
                    }
                } else {
                    self._setTool('delete');
                }
            }
        };
        document.addEventListener('keydown', self._keydownHandler);

        // Ctrl+V paste SMILES / MOL text into editor
        self._pasteHandler = function(e) {
            if (self._destroyed) return;
            if (self._options.depict) return;
            if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;
            var text = (e.clipboardData || window.clipboardData).getData('text');
            if (text && text.trim().length > 0 && text.trim().length < 10000) {
                e.preventDefault();
                self.readGenericMolecularInput(text.trim());
                self.showInfo('Pasted: ' + text.trim().substring(0, 30) + (text.length > 30 ? '...' : ''));
            }
        };
        document.addEventListener('paste', self._pasteHandler);

        // Drag-and-drop MOL/SDF files onto the drawing area. v1.4.1 (bug-fix
        // #4): named handlers stored on self._dropHandlers so destroy() can
        // remove them.
        var container = self._drawArea;
        self._dropTarget = container;
        self._dropHandlers = {
            onDragOver: function(e) { e.preventDefault(); e.dataTransfer.dropEffect = 'copy'; },
            onDrop: function(e) {
                e.preventDefault();
                if (self._destroyed) return;
                if (self._options.depict) return;
                var file = e.dataTransfer.files[0];
                if (!file) return;
                var reader = new FileReader();
                reader.onload = function(ev) {
                    self.readGenericMolecularInput(ev.target.result);
                    self.showInfo('Loaded: ' + file.name);
                };
                reader.readAsText(file);
            }
        };
        container.addEventListener('dragover', self._dropHandlers.onDragOver);
        container.addEventListener('drop', self._dropHandlers.onDrop);

        // Responsive resize observer
        if (typeof ResizeObserver !== 'undefined') {
            self._resizeObserver = new ResizeObserver(function() {
                var w = self._drawArea.clientWidth;
                var h = self._drawArea.clientHeight;
                if (w > 0 && h > 0) {
                    self.renderer.setSize(w, h);
                    self.render();
                }
            });
            self._resizeObserver.observe(self._drawArea);
        }
    };

    /**
     * Tear down the editor: remove global listeners, timers, and observers.
     * v1.4.1 (bug-fix #4): now removes ALL event listeners — was previously
     * leaking 9 (click, dblclick, two touchstarts, two touchmoves, touchend,
     * wheel, dragover/drop on the drawing area).
     */
    MolEditor.prototype.destroy = function() {
        if (this._destroyed) return;
        this._destroyed = true;
        document.removeEventListener('keydown', this._keydownHandler);
        document.removeEventListener('paste', this._pasteHandler);
        if (this._resizeObserver) this._resizeObserver.disconnect();
        // Remove SVG event listeners
        var svg = this.renderer && this.renderer.svg;
        var h = this._svgHandlers;
        if (svg && h) {
            svg.removeEventListener('mousedown',  h.onDown);
            svg.removeEventListener('mousemove',  h.onMove);
            svg.removeEventListener('mouseup',    h.onUp);
            svg.removeEventListener('click',      h.onClick);
            svg.removeEventListener('dblclick',   h.onDblClick);
            svg.removeEventListener('touchstart', h.onDown);
            svg.removeEventListener('touchmove',  h.onMove);
            svg.removeEventListener('touchend',   h.onUp);
            svg.removeEventListener('touchstart', h.onPinchStart);
            svg.removeEventListener('touchmove',  h.onPinchMove);
            svg.removeEventListener('wheel',      h.onWheel);
        }
        // Remove drop-zone listeners on the drawing-area container.
        if (this._dropTarget && this._dropHandlers) {
            this._dropTarget.removeEventListener('dragover', this._dropHandlers.onDragOver);
            this._dropTarget.removeEventListener('drop',     this._dropHandlers.onDrop);
        }
        this._svgHandlers = null;
        this._dropHandlers = null;
        this._dropTarget = null;
    };

    MolEditor.prototype._clearHighlights = function() {
        // FIX: also clear bgColor halos. Without this, SMARTS-match teal
        // backgrounds linger after a generic highlight reset (e.g. ESC,
        // selection clear) until the next match call.
        this.molecule.atoms.forEach(function(a) {
            a.highlighted = false;
            a.mapHighlighted = false;
            a.bgColor = null;
        });
        this.molecule.bonds.forEach(function(b) {
            b.highlighted = false;
            b.bgColor = null;
        });
    };

    // -----------------------------------------------------------------
    // Atom-atom mapping cross-highlight
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // -----------------------------------------------------------------

    /**
     * Highlight all atoms sharing the given mapNumber across the molecule.
     * Sets atom.mapHighlighted = true on every atom with that map number.
     * Returns the array of matched atom ids.
     */
    MolEditor.prototype._highlightMapPartners = function(mapNumber) {
        var ids = [];
        if (!mapNumber || mapNumber <= 0) return ids;
        this.molecule.atoms.forEach(function(a) {
            if (a.mapNumber === mapNumber) {
                a.highlighted = true;
                a.mapHighlighted = true;
                ids.push(a.id);
            }
        });
        return ids;
    };

    /**
     * Clear only the map-highlight sticky state (leaves normal highlights alone).
     */
    MolEditor.prototype._clearMapHighlights = function() {
        this.molecule.atoms.forEach(function(a) { a.mapHighlighted = false; });
        this._stickyMapNumber = 0;
    };

    // =========================================================================
    // State management
    // =========================================================================
    /** Aromatic display style: 'circle' (inner dashed circle) or 'kekule' (alternating bonds). */
    MolEditor.prototype._aromaticStyle = 'circle';

    MolEditor.prototype.saveHistory = function() {
        this.history.push(this.molecule);
    };

    MolEditor.prototype.changed = function() {
        this.render();
        this._fireCallback('AfterStructureModified');
    };

    MolEditor.prototype.render = function() {
        this.renderer.render();
    };

    MolEditor.prototype._restoreMolecule = function(json) {
        this.molecule.clear();
        if (json.name !== undefined) this.molecule.name = json.name;
        // FIX: also restore MOL header program / comment lines so undo/redo
        // round-trips the full header (added when MOL load saved them).
        if (json.program !== undefined) this.molecule.program = json.program;
        if (json.comment !== undefined) this.molecule.comment = json.comment;
        var idMap = {};
        for (var i = 0; i < json.atoms.length; i++) {
            var a = json.atoms[i];
            var atom = this.molecule.addAtom(a.symbol, a.x, a.y);
            atom.charge = a.charge || 0;
            atom.isotope = a.isotope || 0;
            atom.mapNumber = a.mapNumber || 0;
            atom.hydrogens = a.hydrogens !== undefined ? a.hydrogens : -1;
            // FIX: restore aromatic flag so post-undo renders preserve aromatic styling.
            if (a.aromatic) atom.aromatic = true;
            // FIX: restore chirality (otherwise SMILES @ / @@ is silently dropped
            // on undo), CIP label, AAM cross-highlight pin and SMARTS bgColor halo.
            // Without these, undo/redo silently strips stereo and any active
            // map-partner highlight, contradicting the I/O round-trip contract.
            if (a.chirality) atom.chirality = a.chirality;
            if (a.cipLabel)  atom.cipLabel  = a.cipLabel;
            // FIX: restore radical multiplicity. Radicals from MOL `M  RAD` or
            // V3000 `RAD=` are stored on atom.radical (1=singlet, 2=doublet,
            // 3=triplet); without this, the first undo silently drops them.
            if (a.radical) atom.radical = a.radical;
            if (a.mapHighlighted) {
                atom.mapHighlighted = true;
                atom.highlighted    = true; // renderer pairs the two flags
            }
            if (a.bgColor)   atom.bgColor   = a.bgColor;
            idMap[a.id] = atom.id;
        }
        for (var i = 0; i < json.bonds.length; i++) {
            var b = json.bonds[i];
            // FIX: use 'in' check; || would fail if the mapped id happened to be 0 (falsy)
            var a1 = (b.atom1 in idMap) ? idMap[b.atom1] : b.atom1;
            var a2 = (b.atom2 in idMap) ? idMap[b.atom2] : b.atom2;
            var bond = this.molecule.addBond(a1, a2, b.type);
            if (bond) {
                bond.stereo = b.stereo || 0;
                // FIX: restore bond CIP label (E/Z) and bgColor highlight halo.
                if (b.cipLabel) bond.cipLabel = b.cipLabel;
                if (b.bgColor)  bond.bgColor  = b.bgColor;
            }
        }
        if (json.reactionArrow) this.molecule.reactionArrow = json.reactionArrow;
        // FIX: also restore reactionPlusSigns so undo/redo of reaction edits
        // doesn't drop the plus signs between components.
        if (json.reactionPlusSigns) this.molecule.reactionPlusSigns = json.reactionPlusSigns;
        // v1.4.1 (bug-fix #2): clear stale RDT halo state. componentPairs and
        // subFragments reference atom IDs from the pre-undo molecule;
        // _restoreMolecule reissues IDs through addAtom's monotonic counter
        // (idMap above), so the cached halo atom-id arrays are invalid and
        // would either render no halos (best case) or splash colour onto
        // unrelated atoms whose IDs happened to overlap. Force a clean slate.
        this.renderer.componentPairs = [];
        this.molecule.subFragments = [];
        this.molecule.componentPairs = [];
        this.renderer.setMolecule(this.molecule);
    };

    // =========================================================================
    // Callbacks
    // =========================================================================
    MolEditor.prototype.setCallBack = function(name, fn) {
        this._callbacks[name] = fn;
    };

    MolEditor.prototype.getCallBack = function(name) {
        return this._callbacks[name] || null;
    };

    MolEditor.prototype._fireCallback = function(name, extra) {
        var fn = this._callbacks[name];
        if (!fn) return;
        var event = {
            src: this,
            action: name,
            atom: (extra && extra.atom) || 0,
            bond: (extra && extra.bond) || 0,
            molecule: 0,
            argument: null
        };
        fn(event);
    };

    // =========================================================================
    // Public API — JSME-compatible
    // =========================================================================
    MolEditor.prototype.smiles = function() {
        if (typeof SmilesWriter !== 'undefined') return SmilesWriter.write(this.molecule);
        return ''; // SmilesWriter not loaded yet
    };

    MolEditor.prototype.molFile = function(v3000) {
        if (typeof MolfileWriter !== 'undefined') return MolfileWriter.write(this.molecule, v3000);
        return '';
    };

    MolEditor.prototype.jmeFile = function() {
        return ''; // JME format not implemented
    };

    MolEditor.prototype.readGenericMolecularInput = function(input) {
        if (!input || !input.trim()) return;
        input = input.trim();
        // Save history so the load is undoable via Ctrl+Z
        this.saveHistory();
        // Invalidate any in-flight async name lookups
        this._nameRequestSeq = (this._nameRequestSeq || 0) + 1;
        // Only show name if the molecule already has one (e.g. from MOL file header)
        // Detect format
        if (input.indexOf('$RXN') >= 0) {
            // RXN file
            this.molecule.clear();
            this._parseRxnFile(input);
            this.renderer.centerMolecule();
            this.changed();
            this._updateNameStatus();
        } else if (input.indexOf('V2000') >= 0 || input.indexOf('V3000') >= 0 || input.indexOf('M  END') >= 0) {
            // MOL or SDF format — parse the molblock up to M  END, then check for SDF > <NAME> field
            this.molecule.clear();
            var molEnd = input.indexOf('M  END');
            var molPart = molEnd >= 0 ? input.substring(0, molEnd + 6) : input;
            if (molPart.indexOf('V3000') >= 0) {
                this._parseMolV3000(molPart, this.molecule);
            } else {
                this._parseMolV2000(molPart, this.molecule);
            }
            // If this is an SDF block, extract data fields after M  END.
            // FIX: accept classic-Mac \r as a line separator too.
            var nameMatch = input.match(/^> <NAME>(?:\r\n|\r|\n)([^\r\n]*)/m);
            if (nameMatch && nameMatch[1].trim()) {
                this.molecule.name = nameMatch[1].trim();
            }
            // FIX: also extract > <COMMENT> if present (round-trips with writer).
            var commentMatch = input.match(/^> <COMMENT>(?:\r\n|\r|\n)([^\r\n]*)/m);
            if (commentMatch && commentMatch[1].trim()) {
                this.molecule.comment = commentMatch[1].trim();
            }
            this.renderer.centerMolecule();
            this.changed();
            this._updateNameStatus();
        } else {
            // Assume SMILES
            if (typeof SmilesParser !== 'undefined') {
                this.molecule.clear();
                SmilesParser.parse(input, this.molecule);
                // If SmilesParser extracted a name from the input (e.g. "CCO ethanol"),
                // keep it; otherwise try the common-molecules database and NCI Cactus
                if (!this.molecule.name) {
                    this.molecule.name = this._lookupCommonName(input) || '';
                }
                // External name lookup removed in v1.0.1 for privacy.
                this.renderer.centerMolecule();
                this.changed();
                this._updateNameStatus();
            }
        }
    };

    // =========================================================================
    // MOL V2000 Parser
    // =========================================================================

    /**
     * Parse a MOL V2000 format string into a Molecule object.
     * Reads the atom block (including atom-atom mapping column), bond block,
     * and M  CHG / M  ISO property lines.
     *
     * @param {string} molText  The full MOL file text
     * @param {Molecule} mol    Target molecule to populate
     */
    MolEditor.prototype._parseMolV2000 = function(molText, mol) {
        // FIX: accept classic Mac (\r), Unix (\n), and Windows (\r\n) line
        // endings.  The previous /\r?\n/ split silently treated a \r-only
        // file (rare but legal) as a single line and produced an empty mol.
        var lines = molText.split(/\r\n|\r|\n/);
        if (lines.length < 4) return;

        // Header: line 0 = name, line 1 = program/timestamp stamp, line 2 = comment.
        // FIX: round-trip program (line 1) and comment (line 2) so users do not
        //      lose external authoring metadata when the file is re-saved.
        mol.name    = (lines[0] || '').replace(/\s+$/, '');
        mol.program = (lines[1] || '').replace(/\s+$/, '');
        mol.comment = (lines[2] || '').replace(/\s+$/, '');

        // Counts line (line 3): aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
        var countsLine = lines[3] || '';
        var nAtoms = parseInt(countsLine.substring(0, 3).trim(), 10) || 0;
        var nBonds = parseInt(countsLine.substring(3, 6).trim(), 10) || 0;

        // Atom block: lines 4 .. 4+nAtoms-1
        var atomIdxMap = {}; // 1-based index -> atom id
        var lineIdx = 4;
        for (var i = 0; i < nAtoms; i++) {
            if (lineIdx + i >= lines.length) break;
            var aLine = lines[lineIdx + i];
            if (aLine.length < 34) continue;

            var x = parseFloat(aLine.substring(0, 10).trim()) * 30;
            var y = -parseFloat(aLine.substring(10, 20).trim()) * 30;
            // z = aLine.substring(20, 30) — ignored for 2D
            var sym = aLine.substring(31, 34).trim();

            var atom = mol.addAtom(sym, x, y);
            atomIdxMap[i + 1] = atom.id;

            // FIX: parse mass-difference field (cols 35-36) as an isotope
            // fallback.  Spec value is delta from element nominal mass.
            // M ISO will override later if also present.
            if (aLine.length >= 36) {
                var massDiff = parseInt(aLine.substring(34, 36).trim(), 10) || 0;
                if (massDiff !== 0) {
                    var elem = Molecule.ELEMENTS[atom.symbol];
                    if (elem && elem.mass) atom.isotope = elem.mass + massDiff;
                }
            }

            // Charge field (old-style, columns 37-39)
            if (aLine.length >= 39) {
                var chgCode = parseInt(aLine.substring(36, 39).trim(), 10) || 0;
                // FIX: code 4 is doublet radical (no charge); set radical=2.
                if (chgCode === 4) {
                    atom.radical = 2;
                } else {
                    atom.charge = _molChargeFromCode(chgCode);
                }
            }

            // Atom-atom mapping (columns 60-63, the mmm field)
            if (aLine.length >= 63) {
                var aam = parseInt(aLine.substring(60, 63).trim(), 10) || 0;
                if (aam > 0) atom.mapNumber = aam;
            }
        }

        // Bond block: lines 4+nAtoms .. 4+nAtoms+nBonds-1
        lineIdx = 4 + nAtoms;
        for (var i = 0; i < nBonds; i++) {
            if (lineIdx + i >= lines.length) break;
            var bLine = lines[lineIdx + i];
            if (bLine.length < 12) continue;

            var a1 = parseInt(bLine.substring(0, 3).trim(), 10);
            var a2 = parseInt(bLine.substring(3, 6).trim(), 10);
            var bt = parseInt(bLine.substring(6, 9).trim(), 10) || 1;
            var stereoCode = parseInt(bLine.substring(9, 12).trim(), 10) || 0;

            // Bond type: 1=single, 2=double, 3=triple, 4=aromatic
            var bondType = (bt >= 1 && bt <= 3) ? bt : Molecule.BOND_SINGLE;
            var id1 = atomIdxMap[a1];
            var id2 = atomIdxMap[a2];
            if (id1 !== undefined && id2 !== undefined) {
                var bond = mol.addBond(id1, id2, bondType);
                if (bond) {
                    if (stereoCode === 1) bond.stereo = Molecule.STEREO_WEDGE;
                    else if (stereoCode === 6) bond.stereo = Molecule.STEREO_DASH;
                }
            }
        }

        // Property lines (M  CHG, M  ISO)
        lineIdx = 4 + nAtoms + nBonds;
        for (var i = lineIdx; i < lines.length; i++) {
            var pLine = lines[i];
            if (pLine.indexOf('M  END') >= 0) break;

            if (pLine.indexOf('M  CHG') === 0) {
                var count = parseInt(pLine.substring(6, 9).trim(), 10) || 0;
                for (var ci = 0; ci < count; ci++) {
                    var off = 9 + ci * 8;
                    var idx = parseInt(pLine.substring(off, off + 4).trim(), 10);
                    var chg = parseInt(pLine.substring(off + 4, off + 8).trim(), 10);
                    var aid = atomIdxMap[idx];
                    if (aid !== undefined) {
                        var a = mol.getAtom(aid);
                        if (a) a.charge = chg;
                    }
                }
            } else if (pLine.indexOf('M  ISO') === 0) {
                var count = parseInt(pLine.substring(6, 9).trim(), 10) || 0;
                for (var ci = 0; ci < count; ci++) {
                    var off = 9 + ci * 8;
                    var idx = parseInt(pLine.substring(off, off + 4).trim(), 10);
                    var iso = parseInt(pLine.substring(off + 4, off + 8).trim(), 10);
                    var aid = atomIdxMap[idx];
                    if (aid !== undefined) {
                        var a = mol.getAtom(aid);
                        if (a) a.isotope = iso;
                    }
                }
            } else if (pLine.indexOf('M  RAD') === 0) {
                // FIX: parse radical block.  Format: M  RAD nn8 (idx,rad) pairs.
                // rad: 1=singlet, 2=doublet, 3=triplet.
                var count = parseInt(pLine.substring(6, 9).trim(), 10) || 0;
                for (var ci = 0; ci < count; ci++) {
                    var off = 9 + ci * 8;
                    var idx = parseInt(pLine.substring(off, off + 4).trim(), 10);
                    var rad = parseInt(pLine.substring(off + 4, off + 8).trim(), 10) || 0;
                    var aid = atomIdxMap[idx];
                    if (aid !== undefined) {
                        var a = mol.getAtom(aid);
                        if (a) a.radical = rad;
                    }
                }
            }
            // Unknown M  XXX blocks (STY, SAL, ALS, RGP, ...) fall through
            // and are silently skipped by the loop, per CTfile spec.
        }
    };

    /**
     * Convert old-style MOL charge code to actual charge value.
     * 0=0, 1=+3, 2=+2, 3=+1, 4=doublet radical, 5=-1, 6=-2, 7=-3
     */
    function _molChargeFromCode(code) {
        switch (code) {
            case 1: return 3;
            case 2: return 2;
            case 3: return 1;
            case 5: return -1;
            case 6: return -2;
            case 7: return -3;
            default: return 0;
        }
    }

    // =========================================================================
    // MOL V3000 Parser
    // =========================================================================

    /**
     * Parse a MOL V3000 format string into a Molecule object.
     * Reads CTAB atom/bond blocks with optional properties (CHG, MASS, AAM).
     *
     * @param {string} molText  The full MOL file text
     * @param {Molecule} mol    Target molecule to populate
     */
    MolEditor.prototype._parseMolV3000 = function(molText, mol) {
        // FIX: accept classic Mac (\r) line endings as well as \n and \r\n.
        var rawLines = molText.split(/\r\n|\r|\n/);

        // FIX: V3000 supports line continuation with a trailing '-'.  Stitch
        //      continuation fragments together so long atom/bond records (with
        //      many keyword props) parse correctly.  Only V30 records use this.
        var lines = [];
        for (var li = 0; li < rawLines.length; li++) {
            var ln = rawLines[li];
            while (ln.indexOf('M  V30 ') === 0 && /-\s*$/.test(ln) && li + 1 < rawLines.length) {
                ln = ln.replace(/-\s*$/, '');
                var next = rawLines[li + 1] || '';
                if (next.indexOf('M  V30 ') === 0) next = next.substring(7);
                ln = ln + next;
                li++;
            }
            lines.push(ln);
        }

        // FIX: round-trip program (line 1) and comment (line 2) like V2000.
        mol.name    = (lines[0] || '').replace(/\s+$/, '');
        mol.program = (lines[1] || '').replace(/\s+$/, '');
        mol.comment = (lines[2] || '').replace(/\s+$/, '');

        var inAtomBlock = false;
        var inBondBlock = false;
        var atomIdxMap = {}; // V3000 idx -> atom id

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];

            if (line.indexOf('M  V30 BEGIN ATOM') >= 0) { inAtomBlock = true; continue; }
            if (line.indexOf('M  V30 END ATOM') >= 0)   { inAtomBlock = false; continue; }
            if (line.indexOf('M  V30 BEGIN BOND') >= 0) { inBondBlock = true; continue; }
            if (line.indexOf('M  V30 END BOND') >= 0)   { inBondBlock = false; continue; }

            if (inAtomBlock && line.indexOf('M  V30 ') === 0) {
                // M  V30 idx type x y z aamap [keyword=value ...]
                var parts = line.substring(7).trim().split(/\s+/);
                if (parts.length < 5) continue;
                var idx = parseInt(parts[0], 10);
                var sym = parts[1];
                var x = parseFloat(parts[2]) * 30;
                var y = -parseFloat(parts[3]) * 30;

                var atom = mol.addAtom(sym, x, y);
                atomIdxMap[idx] = atom.id;

                // FIX: parts[5] is the positional atom-atom map number per
                //      V3000 spec ('aamap' field).  Only treat it as AAM if
                //      it is a positive integer with no '=' sign.  Previously
                //      this field was silently skipped, so AAMs in spec-
                //      compliant V3000 RXN files were lost.
                var keywordStart = 5;
                if (parts.length > 5 && parts[5].indexOf('=') < 0) {
                    var aamPos = parseInt(parts[5], 10);
                    if (!isNaN(aamPos) && aamPos > 0) atom.mapNumber = aamPos;
                    keywordStart = 6;
                }

                // Parse optional keyword properties: CHG=n, MASS=n, RAD=n,
                //                                    AAM=n (BIME extension), CFG=n
                for (var pi = keywordStart; pi < parts.length; pi++) {
                    var prop = parts[pi];
                    if (prop.indexOf('CHG=') === 0) {
                        atom.charge = parseInt(prop.substring(4), 10) || 0;
                    } else if (prop.indexOf('MASS=') === 0) {
                        atom.isotope = parseInt(prop.substring(5), 10) || 0;
                    } else if (prop.indexOf('AAM=') === 0) {
                        atom.mapNumber = parseInt(prop.substring(4), 10) || 0;
                    } else if (prop.indexOf('RAD=') === 0) {
                        // FIX: V3000 RAD= per-atom radical multiplicity.
                        // 1=singlet, 2=doublet, 3=triplet (CTfile spec).
                        atom.radical = parseInt(prop.substring(4), 10) || 0;
                    } else if (prop.indexOf('ATTCHPT=') === 0) {
                        // Legacy: treat ATTCHPT as AAM for backward compatibility
                        atom.mapNumber = parseInt(prop.substring(8), 10) || 0;
                    }
                    // Unknown V3000 atom keywords are ignored silently per spec.
                }
            }

            if (inBondBlock && line.indexOf('M  V30 ') === 0) {
                // M  V30 idx type atom1 atom2 [props...]
                var parts = line.substring(7).trim().split(/\s+/);
                if (parts.length < 4) continue;
                var bt = parseInt(parts[1], 10) || 1;
                var a1 = parseInt(parts[2], 10);
                var a2 = parseInt(parts[3], 10);

                var bondType = (bt >= 1 && bt <= 3) ? bt : Molecule.BOND_SINGLE;
                var id1 = atomIdxMap[a1];
                var id2 = atomIdxMap[a2];
                if (id1 !== undefined && id2 !== undefined) {
                    var bond = mol.addBond(id1, id2, bondType);
                    // Parse optional bond properties
                    for (var pi = 4; pi < parts.length; pi++) {
                        var prop = parts[pi];
                        if (prop.indexOf('CFG=') === 0) {
                            var cfg = parseInt(prop.substring(4), 10) || 0;
                            if (bond) {
                                if (cfg === 1) bond.stereo = Molecule.STEREO_WEDGE;
                                else if (cfg === 3) bond.stereo = Molecule.STEREO_DASH;
                            }
                        }
                    }
                }
            }
        }
    };

    // =========================================================================
    // RXN File Parser
    // =========================================================================

    /**
     * Parse a $RXN file into the molecule, setting up reaction arrow and
     * preserving atom-atom mapping numbers from the MOL blocks.
     *
     * RXN format:
     *   $RXN
     *   [name line]
     *   [program line]
     *   [comment line]
     *   rrrppp        (counts: reactants, products)
     *   $MOL
     *   [MOL block 1 - reactant]
     *   ...
     *   $MOL
     *   [MOL block N - product]
     *
     * @param {string} rxnText  The full RXN file text
     */
    MolEditor.prototype._parseRxnFile = function(rxnText) {
        var mol = this.molecule;
        // FIX: accept all three line-ending conventions.
        var lines = rxnText.split(/\r\n|\r|\n/);

        // Find $RXN header.  FIX: also detect 'V3000' marker on the same line
        //                   for V3000 RXN files ("$RXN V3000").
        var headerIdx = -1;
        var rxnIsV3000 = false;
        for (var i = 0; i < lines.length; i++) {
            var trimmed = lines[i].trim();
            if (trimmed === '$RXN' || trimmed.indexOf('$RXN') === 0) {
                headerIdx = i;
                if (trimmed.indexOf('V3000') >= 0) rxnIsV3000 = true;
                break;
            }
        }
        if (headerIdx < 0) return;

        // RXN header lines after $RXN: name, program, comment, counts.
        // FIX: preserve program (line 2) and comment (line 3) onto the
        //      molecule so the metadata round-trips on save.
        mol.name    = (lines[headerIdx + 1] || '').replace(/\s+$/, '');
        mol.program = (lines[headerIdx + 2] || '').replace(/\s+$/, '');
        mol.comment = (lines[headerIdx + 3] || '').replace(/\s+$/, '');

        var countsIdx = headerIdx + 4;
        if (countsIdx >= lines.length) return;
        var countsLine = lines[countsIdx].trim();
        var nReactants, nProducts;
        if (rxnIsV3000) {
            // FIX: V3000 RXN uses a free-format counts line that begins with
            //      'M  V30 COUNTS nR nP [nA]' — fall back to whitespace split.
            var m3 = countsLine.match(/COUNTS\s+(\d+)\s+(\d+)/);
            if (m3) {
                nReactants = parseInt(m3[1], 10) || 0;
                nProducts  = parseInt(m3[2], 10) || 0;
            } else {
                // Try whitespace-separated counts on a bare line
                var partsC = countsLine.split(/\s+/);
                nReactants = parseInt(partsC[0], 10) || 0;
                nProducts  = parseInt(partsC[1], 10) || 0;
            }
        } else {
            nReactants = parseInt(countsLine.substring(0, 3).trim(), 10) || 0;
            nProducts  = parseInt(countsLine.substring(3, 6).trim(), 10) || 0;
            // Note: the optional nAgents (cols 6-9) is intentionally ignored.
        }

        // Extract $MOL blocks.  FIX: $MOL marker must start at column 0 to
        //                       avoid matching the literal inside comments.
        var molBlocks = [];
        var blockStart = -1;
        for (var i = countsIdx + 1; i < lines.length; i++) {
            if (lines[i].indexOf('$MOL') === 0) {
                if (blockStart >= 0) {
                    molBlocks.push(lines.slice(blockStart + 1, i).join('\n'));
                }
                blockStart = i;
            }
        }
        // Last block
        if (blockStart >= 0) {
            molBlocks.push(lines.slice(blockStart + 1).join('\n'));
        }

        // Parse each MOL block into temporary molecules
        var BOND_LENGTH = Molecule.BOND_LENGTH;
        var FRAG_GAP  = BOND_LENGTH * 1.5;
        var ARROW_PAD = BOND_LENGTH * 1.2;
        var ARROW_LEN = BOND_LENGTH * 3;
        var PLUS_WIDTH = BOND_LENGTH * 0.6;

        var reactantMols = [];
        var productMols = [];

        for (var mi = 0; mi < molBlocks.length; mi++) {
            var tmpMol = new Molecule();
            if (molBlocks[mi].indexOf('V3000') >= 0) {
                this._parseMolV3000(molBlocks[mi], tmpMol);
            } else {
                this._parseMolV2000(molBlocks[mi], tmpMol);
            }
            if (mi < nReactants) {
                reactantMols.push(tmpMol);
            } else {
                productMols.push(tmpMol);
            }
        }

        // Layout reactants and products with arrow
        var globalMinY = Infinity, globalMaxY = -Infinity;
        var allMols = reactantMols.concat(productMols);
        for (var mi = 0; mi < allMols.length; mi++) {
            var b = allMols[mi].getBounds();
            if (b.y < globalMinY) globalMinY = b.y;
            if (b.y + b.h > globalMaxY) globalMaxY = b.y + b.h;
        }
        var globalCentreY = (globalMinY + globalMaxY) / 2;

        var cursorX = 0;
        var plusPositions = [];

        // Place reactants
        for (var ri = 0; ri < reactantMols.length; ri++) {
            var rBounds = reactantMols[ri].getBounds();
            var fragCentreY = rBounds.y + rBounds.h / 2;
            var shiftX = cursorX - rBounds.x;
            var shiftY = globalCentreY - fragCentreY;
            this._copyMolInto(mol, reactantMols[ri], shiftX, shiftY);
            cursorX += rBounds.w + FRAG_GAP;
            if (ri < reactantMols.length - 1) {
                plusPositions.push({ x: cursorX - FRAG_GAP / 2, y: globalCentreY });
                cursorX += PLUS_WIDTH;
            }
        }

        // Arrow
        var arrowX1 = cursorX + ARROW_PAD;
        var arrowX2 = arrowX1 + ARROW_LEN;
        mol.reactionArrow = {
            x1: arrowX1, y1: globalCentreY,
            x2: arrowX2, y2: globalCentreY,
            type: 'forward',
            conditions: ''
        };

        // Place products
        cursorX = arrowX2 + ARROW_PAD;
        for (var pi = 0; pi < productMols.length; pi++) {
            var pBounds = productMols[pi].getBounds();
            var fragCentreY = pBounds.y + pBounds.h / 2;
            var shiftX = cursorX - pBounds.x;
            var shiftY = globalCentreY - fragCentreY;
            this._copyMolInto(mol, productMols[pi], shiftX, shiftY);
            cursorX += pBounds.w + FRAG_GAP;
            if (pi < productMols.length - 1) {
                plusPositions.push({ x: cursorX - FRAG_GAP / 2, y: globalCentreY });
                cursorX += PLUS_WIDTH;
            }
        }

        mol.reactionPlusSigns = plusPositions;
    };

    /**
     * Copy all atoms and bonds from srcMol into destMol with coordinate offsets.
     * Preserves all atom properties including mapNumber.
     */
    MolEditor.prototype._copyMolInto = function(destMol, srcMol, dx, dy) {
        var idMap = {};
        for (var i = 0; i < srcMol.atoms.length; i++) {
            var a = srcMol.atoms[i];
            var newAtom = destMol.addAtom(a.symbol, a.x + dx, a.y + dy);
            newAtom.charge    = a.charge;
            newAtom.isotope   = a.isotope;
            newAtom.mapNumber = a.mapNumber;
            newAtom.hydrogens = a.hydrogens;
            newAtom.aromatic  = a.aromatic;
            newAtom.chirality = a.chirality;
            idMap[a.id] = newAtom.id;
        }
        for (var i = 0; i < srcMol.bonds.length; i++) {
            var b = srcMol.bonds[i];
            var bond = destMol.addBond(idMap[b.atom1], idMap[b.atom2], b.type);
            if (bond) bond.stereo = b.stereo;
        }
    };

    MolEditor.prototype.options = function(optString) {
        this._parseOptions(optString);
        this.renderer.depictMode = !!this._options.depict;
        this.renderer.showHydrogens = this._options.hydrogens !== false;
        this.renderer.showMapNumbers = true;
        if (this._options.depict) {
            this.renderer.svg.style.cursor = 'default';
        } else {
            this.renderer.svg.style.cursor = 'crosshair';
        }
        this.render();
    };

    MolEditor.prototype._parseOptions = function(optString) {
        if (!optString) return;
        var parts = optString.split(/[,\s]+/);
        for (var i = 0; i < parts.length; i++) {
            var opt = parts[i].trim().toLowerCase();
            if (opt.indexOf('no') === 0) {
                this._options[opt.substring(2)] = false;
            } else {
                this._options[opt] = true;
            }
        }
    };

    MolEditor.prototype.reset = function() {
        this.saveHistory();
        this.molecule.clear();
        this.changed();
    };

    MolEditor.prototype.clear = function() { this.reset(); };

    MolEditor.prototype.setSize = function(w, h) {
        this._container.style.width = w;
        this._container.style.height = h;
        var drawW = this._drawArea.clientWidth || 500;
        var drawH = this._drawArea.clientHeight || 380;
        this.renderer.setSize(drawW, drawH);
        this.render();
    };

    MolEditor.prototype.repaint = function() { this.render(); };

    /**
     * Re-run the 2D layout algorithm on the current molecule.
     * Useful after manual edits that create overlapping atoms.
     */
    MolEditor.prototype._layoutAttempt = 0;
    MolEditor.prototype._autoLayout = function() {
        if (this.molecule.isEmpty()) return;
        this.saveHistory();

        if (typeof Layout === 'undefined' || !Layout.layout) {
            this.renderer.centerMolecule();
            this.render();
            return;
        }

        var mol = this.molecule;
        var self = this;
        var attempts = 3;
        var bestCoords = null;
        var bestCrossings = 999999;

        // Try multiple layout seeds, keep the one with fewest crossings
        for (var attempt = 0; attempt < attempts; attempt++) {
            // Reset coords with slight random perturbation for variety
            var seed = ++self._layoutAttempt;
            for (var i = 0; i < mol.atoms.length; i++) {
                // Deterministic pseudo-random offset based on seed + atom index
                var hash = ((seed * 2654435761 + i * 40503) >>> 0) & 0xFFFF;
                mol.atoms[i].x = (hash & 0xFF) * 0.01 - 1.28;
                mol.atoms[i].y = ((hash >> 8) & 0xFF) * 0.01 - 1.28;
            }

            Layout.layout(mol);

            // Count crossings
            var crossings = 0;
            for (var bi = 0; bi < mol.bonds.length; bi++) {
                var b1 = mol.bonds[bi];
                var a1 = mol.getAtom(b1.atom1), a2 = mol.getAtom(b1.atom2);
                if (!a1 || !a2) continue;
                for (var bj = bi + 1; bj < mol.bonds.length; bj++) {
                    var b2 = mol.bonds[bj];
                    if (b2.atom1 === b1.atom1 || b2.atom1 === b1.atom2 ||
                        b2.atom2 === b1.atom1 || b2.atom2 === b1.atom2) continue;
                    var a3 = mol.getAtom(b2.atom1), a4 = mol.getAtom(b2.atom2);
                    if (!a3 || !a4) continue;
                    var d1 = (a4.x-a3.x)*(a1.y-a3.y) - (a4.y-a3.y)*(a1.x-a3.x);
                    var d2 = (a4.x-a3.x)*(a2.y-a3.y) - (a4.y-a3.y)*(a2.x-a3.x);
                    var d3 = (a2.x-a1.x)*(a3.y-a1.y) - (a2.y-a1.y)*(a3.x-a1.x);
                    var d4 = (a2.x-a1.x)*(a4.y-a1.y) - (a2.y-a1.y)*(a4.x-a1.x);
                    if (((d1>0&&d2<0)||(d1<0&&d2>0))&&((d3>0&&d4<0)||(d3<0&&d4>0))) crossings++;
                }
            }

            if (crossings < bestCrossings) {
                bestCrossings = crossings;
                bestCoords = [];
                for (var ai = 0; ai < mol.atoms.length; ai++) {
                    bestCoords.push({ x: mol.atoms[ai].x, y: mol.atoms[ai].y });
                }
            }

            if (bestCrossings === 0) break; // perfect layout found
        }

        // Apply the best layout
        if (bestCoords) {
            for (var ai = 0; ai < mol.atoms.length; ai++) {
                mol.atoms[ai].x = bestCoords[ai].x;
                mol.atoms[ai].y = bestCoords[ai].y;
            }
        }

        this.renderer.centerMolecule();
        this.render();
        this.showInfo('Layout applied' + (bestCrossings > 0 ? ' (' + bestCrossings + ' crossings)' : '') + ' — Ctrl+Z to undo');
        this._fireCallback('AfterStructureModified');
    };

    /**
     * Zoom to fit the entire molecule in the viewport with padding.
     */
    MolEditor.prototype._zoomToFit = function() {
        if (this.molecule.isEmpty()) return;
        var bounds = this.molecule.getBounds();
        if (bounds.w < 1 && bounds.h < 1) return;
        var padX = 60, padY = 60;
        var scaleX = (this.renderer.width - padX) / (bounds.w || 1);
        var scaleY = (this.renderer.height - padY) / (bounds.h || 1);
        this.renderer.scale = Math.min(scaleX, scaleY, 3); // cap at 3x
        this.renderer.scale = Math.max(this.renderer.scale, 0.2); // min 0.2x
        this.renderer.offsetX = (this.renderer.width / this.renderer.scale - bounds.w) / 2 - bounds.x;
        this.renderer.offsetY = (this.renderer.height / this.renderer.scale - bounds.h) / 2 - bounds.y;
        this.render();
        this.showInfo('Zoom: ' + Math.round(this.renderer.scale * 100) + '%');
    };

    MolEditor.prototype.totalNumberOfAtoms = function() { return this.molecule.atoms.length; };
    MolEditor.prototype.totalNumberOfBonds = function() { return this.molecule.bonds.length; };

    MolEditor.prototype.getAtom = function(molIdx, atomIdx) {
        var atom = this.molecule.atoms[atomIdx - 1]; // 1-indexed
        if (!atom) return null;
        return {
            x: atom.x, y: atom.y, atom: atom.symbol, charge: atom.charge,
            isotope: atom.isotope, mapNumber: atom.mapNumber,
            numExplicitHydrogens: this.molecule.calcHydrogens(atom.id),
            markingNumber: 0, backGroundColor: 0
        };
    };

    MolEditor.prototype.getBond = function(molIdx, bondIdx) {
        var bond = this.molecule.bonds[bondIdx - 1]; // 1-indexed
        if (!bond) return null;
        return {
            atom1: bond.atom1, atom2: bond.atom2, bondType: bond.type,
            stereo: bond.stereo, reactionComponent: 0, markingNumber: 0, backGroundColor: 0
        };
    };

    MolEditor.prototype.showInfo = function(msg) {
        var info = this._container.querySelector('#bime-info');
        if (info) info.textContent = msg || '';
    };

    /**
     * Validate a SMILES string: parse it and report atom/bond count or error.
     * Returns { valid: bool, atoms: n, bonds: n, smiles: canonicalSmiles, error: string|null }
     */
    MolEditor.prototype.validateSmiles = function(smilesStr) {
        if (!smilesStr || !smilesStr.trim()) return { valid: false, atoms: 0, bonds: 0, smiles: '', error: 'Empty input', warnings: [] };
        try {
            var testMol = new Molecule();
            testMol.parseErrors = [];
            testMol.parseWarnings = [];
            SmilesParser.parse(smilesStr.trim(), testMol);

            // Check for parse errors (fatal)
            if (testMol.parseErrors && testMol.parseErrors.length > 0) {
                return { valid: false, atoms: 0, bonds: 0, smiles: '', error: testMol.parseErrors[0], warnings: testMol.parseWarnings || [] };
            }
            if (testMol.atoms.length === 0) {
                return { valid: false, atoms: 0, bonds: 0, smiles: '', error: 'No atoms parsed', warnings: [] };
            }

            var canonical = SmilesWriter.write(testMol);
            var warnings = testMol.parseWarnings || [];
            return { valid: true, atoms: testMol.atoms.length, bonds: testMol.bonds.length, smiles: canonical, error: null, warnings: warnings };
        } catch (e) {
            return { valid: false, atoms: 0, bonds: 0, smiles: '', error: e.message || 'Parse error', warnings: [] };
        }
    };

    /**
     * Show a SMILES validation prompt inline in the status bar.
     */
    MolEditor.prototype._promptValidateSmiles = function() {
        var self = this;
        var bar = this._statusBar;

        // Build inline input
        bar.innerHTML = '<div style="display:flex;align-items:center;gap:0.4rem;flex:1">' +
            '<input id="bime-validate-input" type="text" placeholder="Paste SMILES to validate..." ' +
            'style="flex:1;padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;font-size:11px;font-family:var(--font-mono,monospace);background:var(--color-bg,white);color:var(--color-text,#333)">' +
            '<button id="bime-validate-btn" style="padding:2px 8px;border:1px solid ' + BRAND_COLOR + ';border-radius:4px;background:' + BRAND_COLOR + ';color:white;font-size:11px;font-weight:600;cursor:pointer">Check</button>' +
            '<span id="bime-validate-result" style="font-size:11px;font-weight:600"></span>' +
            '<button id="bime-validate-close" style="padding:2px 6px;border:none;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:13px" title="Close">&times;</button>' +
            '</div>';

        var input = bar.querySelector('#bime-validate-input');
        var btn = bar.querySelector('#bime-validate-btn');
        var result = bar.querySelector('#bime-validate-result');
        var closeBtn = bar.querySelector('#bime-validate-close');

        function doValidate() {
            var val = input.value.trim();
            var res = self.validateSmiles(val);
            if (res.valid) {
                var msg = '\u2713 Valid \u2014 ' + res.atoms + ' atoms, ' + res.bonds + ' bonds';
                if (res.warnings && res.warnings.length > 0) {
                    result.style.color = '#ca8a04';
                    msg += ' \u26a0 ' + res.warnings[0];
                } else {
                    result.style.color = '#059669';
                }
                result.textContent = msg;
            } else {
                result.style.color = '#dc2626';
                result.textContent = '\u2717 ' + (res.error || 'Invalid');
            }
        }

        function doClose() {
            self._resetStatusBar();
        }

        function doLoadAndClose() {
            var val = input.value.trim();
            var res = self.validateSmiles(val);
            if (res.valid) {
                self.readGenericMolecularInput(val);
                doClose();
            }
        }

        btn.addEventListener('click', doValidate);
        closeBtn.addEventListener('click', doClose);
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') {
                if (e.shiftKey) doLoadAndClose();
                else doValidate();
            }
            if (e.key === 'Escape') doClose();
        });
        input.focus();
    };

    // =========================================================================
    // Status bar utilities
    // =========================================================================

    /** Reset the status bar to its default state. */
    MolEditor.prototype._resetStatusBar = function() {
        this._statusBar.innerHTML = '<span id="bime-info"></span><span style="color:' + BRAND_COLOR + ';font-weight:600;font-size:10px">BIME \u2014 BioInception Molecular Editor</span>';
    };

    // =========================================================================
    // Atom-atom mapping (AAM) \u2014 removed in v1.0.3
    //
    // The Java RDT (ReactionDecoder) backend was removed in v1.0.3. BIME now
    // runs entirely in the browser with zero external dependencies and zero
    // external network calls. A pure-JS AAM replacement is planned for
    // v1.1.0; the AAM toolbar button is hidden until then.
    // =========================================================================

    // =========================================================================
    // SMARTS substructure search
    // =========================================================================

    /**
     * Show a SMARTS search prompt inline in the status bar.
     */
    MolEditor.prototype._promptSmartsSearch = function() {
        var self = this;
        var bar = this._statusBar;

        bar.innerHTML = '<div style="display:flex;align-items:center;gap:0.4rem;flex:1">' +
            '<input id="bime-smarts-input" type="text" placeholder="Paste SMARTS pattern (e.g. [#6]~[#7]) ..." ' +
            'style="flex:1;padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;font-size:11px;font-family:var(--font-mono,monospace);background:var(--color-bg,white);color:var(--color-text,#333)">' +
            '<button id="bime-smarts-btn" style="padding:2px 8px;border:1px solid ' + BRAND_COLOR + ';border-radius:4px;background:' + BRAND_COLOR + ';color:white;font-size:11px;font-weight:600;cursor:pointer">Search</button>' +
            '<span id="bime-smarts-result" style="font-size:11px;font-weight:600"></span>' +
            '<button id="bime-smarts-clear" style="padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:10px;font-weight:600" title="Clear highlights">Clear</button>' +
            '<button id="bime-smarts-close" style="padding:2px 6px;border:none;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:13px" title="Close">&times;</button>' +
            '</div>';

        var input = bar.querySelector('#bime-smarts-input');
        var btn = bar.querySelector('#bime-smarts-btn');
        var result = bar.querySelector('#bime-smarts-result');
        var clearBtn = bar.querySelector('#bime-smarts-clear');
        var closeBtn = bar.querySelector('#bime-smarts-close');

        function doSearch() {
            var val = input.value.trim();
            if (!val) { result.textContent = ''; return; }
            var count = self.matchSmarts(val);
            if (count > 0) {
                result.style.color = '#059669';
                result.textContent = '\u2713 ' + count + ' match' + (count > 1 ? 'es' : '');
            } else {
                result.style.color = '#dc2626';
                result.textContent = '\u2717 No match';
            }
            self.render();
        }

        function doClear() {
            self._clearSmartsHighlights();
            self.render();
            result.textContent = '';
        }

        function doClose() {
            self._clearSmartsHighlights();
            self.render();
            self._resetStatusBar();
        }

        btn.addEventListener('click', doSearch);
        clearBtn.addEventListener('click', doClear);
        closeBtn.addEventListener('click', doClose);
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') doSearch();
            if (e.key === 'Escape') doClose();
        });
        input.focus();
    };

    /**
     * Clear SMARTS highlighting from all atoms and bonds.
     */
    MolEditor.prototype._clearSmartsHighlights = function() {
        // FIX: preserve sticky AAM map-partner highlights. SMARTS uses
        // (highlighted + bgColor) for matches, while AAM uses
        // (highlighted + mapHighlighted) for sticky pins. Clearing
        // `highlighted` for atoms whose mapHighlighted flag is set
        // would silently drop the user's pinned map-partner view.
        for (var i = 0; i < this.molecule.atoms.length; i++) {
            var a = this.molecule.atoms[i];
            if (!a.mapHighlighted) a.highlighted = false;
            a.bgColor = null;
        }
        for (var i = 0; i < this.molecule.bonds.length; i++) {
            this.molecule.bonds[i].highlighted = false;
            this.molecule.bonds[i].bgColor = null;
        }
    };

    /**
     * Write SMARTS for the current molecule (if it has query constraints).
     * @returns {string} SMARTS pattern
     */
    MolEditor.prototype.smarts = function() {
        if (typeof SmartsWriter !== 'undefined') return SmartsWriter.write(this.molecule);
        return '';
    };

    /**
     * Match a SMARTS string against the current molecule.
     * Highlights matching atoms in teal.
     * @param {string} smartsString
     * @returns {number} number of matches found
     */
    MolEditor.prototype.matchSmarts = function(smartsString) {
        if (typeof SmartsMatch === 'undefined' || typeof SmartsParser === 'undefined') return 0;
        return SmartsMatch.highlightSmarts(this.molecule, smartsString);
    };

    /**
     * Load a SMARTS pattern as a query molecule into the editor.
     * @param {string} smartsString
     */
    MolEditor.prototype.readSmarts = function(smartsString) {
        if (typeof SmartsParser === 'undefined') return;
        var queryMol = SmartsParser.parse(smartsString);
        if (!queryMol) return;
        this.saveHistory();
        this.molecule.clear();
        // Copy query molecule into editor molecule
        var idMap = {};
        for (var i = 0; i < queryMol.atoms.length; i++) {
            var qa = queryMol.atoms[i];
            var a = this.molecule.addAtom(qa.symbol, qa.x, qa.y);
            a.aromatic = qa.aromatic;
            a.queryConstraints = qa.queryConstraints;
            a.isQuery = true;
            if (qa.hydrogens >= 0) a.hydrogens = qa.hydrogens;
            idMap[qa.id] = a.id;
        }
        for (var i = 0; i < queryMol.bonds.length; i++) {
            var qb = queryMol.bonds[i];
            var bond = this.molecule.addBond(idMap[qb.atom1], idMap[qb.atom2], qb.type);
            if (bond) {
                bond.queryType = qb.queryType;
                bond.isQuery = true;
            }
        }
        this.molecule.isQuery = true;
        this.renderer.centerMolecule();
        this.changed();
    };

    MolEditor.prototype.getMolecularAreaGraphicsString = function() {
        return this.renderer.getSVGString();
    };

    MolEditor.prototype.hasMolecule = function() { return !this.molecule.isEmpty(); };
    MolEditor.prototype.numberOfMolecules = function() { return this.molecule.getComponents().length; };
    MolEditor.prototype.isDepictMode = function() { return !!this._options.depict; };

    MolEditor.prototype.getSupportedFileFormats = function() {
        return ['SMILES', 'MOL V2000', 'MOL V3000', 'RXN'];
    };

    // =========================================================================
    // SMSD (Small Molecule Subgraph Detector) integration — in-browser only
    // MCS, substructure search, and Tanimoto similarity via the bundled
    // SMSDMCS / SMSDGraph / SMSDBatch JS modules. The Java SMSD server was
    // removed in v1.0.3; everything runs entirely in the browser.
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // @see https://github.com/asad/SMSD
    // =========================================================================

    /**
     * Find the Maximum Common Substructure between the current molecule and
     * a second molecule given as SMILES.
     * @param {string} smiles2 - SMILES of the second molecule
     * @returns {Promise<Object>} {mcs, atoms1, atoms2, tanimoto}
     */
    MolEditor.prototype.findMCS = function(smiles2, chemOpts, mcsOpts) {
        var self = this;
        var smiles1 = this.smiles();
        if (!smiles1) return Promise.reject(new Error('No molecule in editor'));
        if (!smiles2)  return Promise.reject(new Error('Second SMILES required'));

        // Built-in JavaScript MCS — runs entirely in the browser, no server.
        if (typeof SMSDMCS !== 'undefined' && typeof SMSDGraph !== 'undefined' && SMSDMCS.findMCS) {
            try {
                var mol2 = new Molecule();
                SmilesParser.parse(smiles2, mol2);
                var SG = SMSDGraph.SMSDGraph || SMSDGraph;
                var g1 = new SG(self.molecule);
                var g2 = new SG(mol2);
                var mcsResult = SMSDMCS.findMCS(g1, g2, chemOpts, mcsOpts);
                var rawMapping = mcsResult.mapping || {};
                // v1.4.1 (bug-fix #15): SMSDMCS returns graph-index mappings;
                // translate to atom-IDs so downstream `_applyMCSHighlights`
                // (which keys into `mol.getAtom(id)`) is correct on
                // molecules with non-monotonic atom IDs (after deletes).
                var atomIdMapping = (typeof SMSDMCS.translateToAtomIds === 'function')
                    ? SMSDMCS.translateToAtomIds(rawMapping, g1, g2)
                    : rawMapping;
                var atoms1 = Object.keys(atomIdMapping).map(Number);
                var atoms2 = [];
                for (var ki = 0; ki < atoms1.length; ki++) {
                    atoms2.push(+atomIdMapping[atoms1[ki]]);
                }
                var localData = {
                    mcs: '',
                    atoms1: atoms1,
                    atoms2: atoms2,
                    mapping: atomIdMapping,
                    tanimoto: mcsResult.tanimoto || 0,
                    size: mcsResult.size || 0
                };
                self.clearHighlights();
                self._lastMCSResult = localData;
                self._lastMCSQuery = smiles2;
                self._applyMCSHighlights(localData);
                self.render();
                return Promise.resolve(localData);
            } catch(e) {
                return Promise.reject(new Error('MCS search failed: ' + (e.message || e)));
            }
        }

        return Promise.reject(new Error('MCS search not available (SMSDMCS module missing)'));
    };

    /**
     * Apply MCS highlighting to the current molecule.
     */
    MolEditor.prototype._applyMCSHighlights = function(data) {
        var mcsSmiles = data.mcs || '';
        // Clear previous MCS map numbers
        for (var c = 0; c < this.molecule.atoms.length; c++) {
            this.molecule.atoms[c].mcsMapNumber = 0;
        }

        var matchedAtoms = 0;
        var matchedBonds = 0;
        var mapCounter = 1;
        var mapping = data.mapping || {};

        if (mcsSmiles && typeof SmartsMatch !== 'undefined' && typeof SmartsParser !== 'undefined') {
            var matches = SmartsMatch.matchSmarts(this.molecule, mcsSmiles);
            if (matches && matches.length > 0) {
                var smartsMapping = matches[0];
                var atomSet = {};
                for (var qId in smartsMapping) {
                    if (!smartsMapping.hasOwnProperty(qId)) continue;
                    var targetAtomId = smartsMapping[qId];
                    var atom = this.molecule.getAtom(targetAtomId);
                    if (atom) {
                        atom.bgColor = '#99f6e4';
                        atom.mcsMapNumber = mapCounter++;
                        atomSet[targetAtomId] = true;
                        matchedAtoms++;
                    }
                }
                for (var b = 0; b < this.molecule.bonds.length; b++) {
                    var bond = this.molecule.bonds[b];
                    if (atomSet[bond.atom1] && atomSet[bond.atom2]) {
                        bond.bgColor = '#99f6e4';
                        matchedBonds++;
                    }
                }
            }
        } else if (data.atoms1 && data.atoms1.length > 0) {
            var atomSet2 = {};
            for (var i = 0; i < data.atoms1.length; i++) {
                var idx = data.atoms1[i];
                if (idx < this.molecule.atoms.length) {
                    this.molecule.atoms[idx].bgColor = '#99f6e4';
                    this.molecule.atoms[idx].mcsMapNumber = mapCounter++;
                    atomSet2[this.molecule.atoms[idx].id] = true;
                    matchedAtoms++;
                }
            }
            for (var b = 0; b < this.molecule.bonds.length; b++) {
                var bond = this.molecule.bonds[b];
                if (atomSet2[bond.atom1] && atomSet2[bond.atom2]) {
                    bond.bgColor = '#99f6e4';
                    matchedBonds++;
                }
            }
        }
        data._matchedAtoms = matchedAtoms;
        data._matchedBonds = matchedBonds;
    };

    /**
     * Fingerprint-based similarity search against common molecules library.
     */
    MolEditor.prototype._doSimilaritySearch = function(smiles1, smiles2, resultEl, exportBtn) {
        var self = this;
        var SG = (typeof SMSDGraph !== 'undefined') ? (SMSDGraph.SMSDGraph || SMSDGraph) : null;
        var Batch = (typeof SMSDBatch !== 'undefined') ? SMSDBatch : null;
        if (!SG || !Batch || !Batch.generateFingerprint || !Batch.tanimotoFingerprint) {
            resultEl.style.color = '#dc2626';
            resultEl.textContent = 'Fingerprint search not available';
            return;
        }

        // If smiles2 given, do pairwise similarity
        if (smiles2) {
            var mol1 = this.molecule;
            var mol2 = new Molecule();
            SmilesParser.parse(smiles2, mol2);
            var g1 = new SG(mol1);
            var g2 = new SG(mol2);
            var fp1 = Batch.generateFingerprint(g1);
            var fp2 = Batch.generateFingerprint(g2);
            var tani = Batch.tanimotoFingerprint(fp1, fp2);
            resultEl.style.color = '#059669';
            resultEl.textContent = 'Tanimoto: ' + tani.toFixed(3);
            resultEl.title = 'Path fingerprint similarity (1024-bit)';
            return;
        }

        // Search against common molecules library
        if (typeof COMMON_MOLECULES === 'undefined' || !COMMON_MOLECULES.length) {
            resultEl.style.color = '#dc2626';
            resultEl.textContent = 'No molecule library loaded';
            return;
        }

        var queryMol = this.molecule;
        var gQ = new SG(queryMol);
        var fpQ = Batch.generateFingerprint(gQ);
        var hits = [];
        for (var i = 0; i < COMMON_MOLECULES.length; i++) {
            var entry = COMMON_MOLECULES[i];
            try {
                var m = new Molecule();
                SmilesParser.parse(entry.smiles, m);
                var g = new SG(m);
                var fp = Batch.generateFingerprint(g);
                var score = Batch.tanimotoFingerprint(fpQ, fp);
                if (score > 0.3) {
                    hits.push({ name: entry.name, smiles: entry.smiles, score: score });
                }
            } catch(e) { /* skip invalid */ }
        }
        hits.sort(function(a, b) { return b.score - a.score; });
        var top = hits.slice(0, 10);
        if (top.length === 0) {
            resultEl.style.color = '#ca8a04';
            resultEl.textContent = 'No similar molecules found (threshold 0.3)';
            return;
        }
        resultEl.style.color = '#059669';
        var topNames = top.map(function(h) { return h.name + ' (' + h.score.toFixed(2) + ')'; });
        resultEl.textContent = 'Top: ' + topNames.slice(0, 5).join(', ');
        resultEl.title = top.map(function(h) { return h.name + ': ' + h.score.toFixed(3) + ' [' + h.smiles + ']'; }).join('\n');

        // Make clicking a result load it
        resultEl.style.cursor = 'pointer';
        resultEl.onclick = function() {
            if (top.length > 0) {
                self.readGenericMolecularInput(top[0].smiles);
                self.showInfo('Loaded: ' + top[0].name);
            }
        };
    };

    /**
     * Show similarity search prompt in the status bar.
     */
    MolEditor.prototype._promptSimSearch = function() {
        var self = this;
        var bar = this._statusBar;
        bar.innerHTML = '<div style="display:flex;align-items:center;gap:0.4rem;flex:1">' +
            '<input id="bime-sim-input" type="text" placeholder="SMILES to compare (leave empty to search library)..." ' +
            'style="flex:1;padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;font-size:11px;font-family:var(--font-mono,monospace);background:var(--color-bg,white);color:var(--color-text,#333)">' +
            '<button id="bime-sim-btn" style="padding:2px 8px;border:1px solid ' + BRAND_COLOR + ';border-radius:4px;background:' + BRAND_COLOR + ';color:white;font-size:11px;font-weight:600;cursor:pointer;white-space:nowrap">Search</button>' +
            '<span id="bime-sim-result" style="font-size:11px;font-weight:600"></span>' +
            '<button id="bime-sim-clear" style="padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:10px;font-weight:600">Clear</button>' +
            '<button id="bime-sim-close" style="padding:2px 6px;border:none;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:13px">&times;</button>' +
            '</div>';

        var input = bar.querySelector('#bime-sim-input');
        var btn = bar.querySelector('#bime-sim-btn');
        var result = bar.querySelector('#bime-sim-result');
        var clearBtn = bar.querySelector('#bime-sim-clear');
        var closeBtn = bar.querySelector('#bime-sim-close');

        function doSearch() {
            var smiles1 = self.smiles();
            if (!smiles1) {
                result.style.color = '#dc2626';
                result.textContent = 'Draw a molecule first';
                return;
            }
            result.style.color = '#ca8a04';
            result.textContent = 'Searching...';
            try {
                self._doSimilaritySearch(smiles1, input.value.trim(), result);
            } catch(e) {
                result.style.color = '#dc2626';
                result.textContent = 'Error: ' + (e.message || 'unknown');
            }
        }

        btn.addEventListener('click', doSearch);
        input.addEventListener('keydown', function(e) { if (e.key === 'Enter') doSearch(); });
        clearBtn.addEventListener('click', function() {
            self.clearHighlights();
            self.render();
            result.textContent = '';
        });
        closeBtn.addEventListener('click', function() {
            self.clearHighlights();
            self.render();
            bar.innerHTML = '';
        });
        input.focus();
    };

    /**
     * Toggle display of atom-atom mapping numbers on all atoms.
     */
    MolEditor.prototype.toggleMapNumbers = function() {
        this.renderer.showMapNumbers = !this.renderer.showMapNumbers;
        this.render();
        return this.renderer.showMapNumbers;
    };

    /**
     * Run the Reaction Decoder Tool (RDT) atom-atom mapping on the current
     * reaction. Requires a reactionArrow on the molecule; otherwise shows
     * an info message and returns. Mutates atom.mapNumber in place and
     * re-renders. The bond-change report is exposed via the status line.
     */
    MolEditor.prototype._runRdtAutoMap = function() {
        if (!this.molecule || !this.molecule.reactionArrow) {
            this.showInfo('Auto-map: draw a reaction arrow first');
            return;
        }
        if (typeof RDT === 'undefined' || !RDT.mapReaction) {
            this.showInfo('Auto-map: RDT module not loaded');
            return;
        }
        var t0 = Date.now();
        var result;
        try {
            result = RDT.mapReaction(this.molecule, { timeoutMs: 10000 });
        } catch (e) {
            // v1.4.1 (bug-fix #1): bail BEFORE saveHistory so a failed AAM
            // doesn't leave a no-op snapshot atop the undo stack.
            this.showInfo('Auto-map failed: ' + (e.message || e));
            return;
        }
        // Only commit a history entry once we have a real mapping result —
        // moves the undo cost off the failure path.
        this.saveHistory();
        var elapsed = Date.now() - t0;
        // v1.4.1 (bug-fix #3): clear stale halo caches BEFORE writing new
        // ones, so a fallback path doesn't leave the previous run's state
        // alongside the current one. Both `subFragments` and `componentPairs`
        // are caches consumers may read independently of the renderer; they
        // must not drift relative to each other or to the live molecule.
        this.molecule.subFragments = [];
        this.molecule.componentPairs = [];
        // Make sure the mapping is visible in the renderer (Number atoms ON).
        this.renderer.showMapNumbers = true;
        // Sync the "Map #" toolbar toggle button state.
        var tmBtn0 = this._toolbarButtons && this._toolbarButtons['togglemap'];
        if (tmBtn0) { tmBtn0.style.opacity = '1'; }
        // v1.4.1 (bug-fix #23): respect a user who has explicitly toggled
        // halos OFF — only force colorAtoms on if it has not been touched
        // (undefined) or is currently true. If the user pressed the Colors
        // button to disable halos, keep their preference.
        if (this.renderer.colorAtoms !== false) { this.renderer.colorAtoms = true; }
        // v1.4.1: populate the renderer's per-atom halo overlay from the AAM
        // result using per-MCS-sub-fragment groups (preserved-bond union-find).
        // Each "sub-fragment" is a maximal connected sub-graph of mapped atoms
        // whose mapped product atoms are also bonded with the same skeleton —
        // that is, a rigid scaffold piece that survived intact through the
        // reaction. Bond-change events break sub-fragments. This produces the
        // RDT-style three-colour visualisation in classic enzyme diagrams
        // (e.g. blue scaffold + green scaffold + orange transferred phosphate).
        var rendererPairs = [];
        var nSubFrags = 0;
        // v1.4.1 (bug-fix #22): allow callers to set MolEditor.options
        // 'minsubfragsize=N' to filter out single-atom sub-fragments and
        // similar clutter. Default 1 (keep all).
        var minSize = 1;
        if (this._options && typeof this._options.minsubfragsize === 'number') {
            minSize = this._options.minsubfragsize;
        }
        // v1.4.1 (bug-fix #18): track whether the sub-fragment derivation
        // succeeded, even when it returned []. A correctly-empty result (zero
        // mapped atoms) should NOT trigger the v1.4.0 fallback — that would
        // paint molecule-pair halos for a reaction the AAM rejected.
        var subFragsSucceeded = false;
        if (typeof RDT.deriveSubFragments === 'function') {
            try {
                var subFrags = RDT.deriveSubFragments(result, { minSize: minSize });
                rendererPairs = subFrags.slice();
                nSubFrags = subFrags.length;
                this.molecule.subFragments = subFrags;
                subFragsSucceeded = true;
            } catch (sfe) {
                // v1.4.1 (bug-fix #17): surface the failure rather than
                // silently swallowing it. A surfaced warning lets users know
                // the renderer has fallen back to component-pair halos.
                if (typeof console !== 'undefined' && console.warn) {
                    console.warn('BIME: RDT.deriveSubFragments failed', sfe);
                }
                rendererPairs = [];
                subFragsSucceeded = false;
            }
        }
        // v1.4.1 fallback: ONLY when sub-fragments truly errored (not when
        // they correctly returned []) does the v1.4.0 component-pair painter
        // kick in — otherwise zero-mapped reactions get spurious halos.
        if (!subFragsSucceeded && typeof RDT.deriveComponentPairs === 'function') {
            try {
                var pal = (typeof Renderer !== 'undefined' && Renderer.COMPONENT_PAIR_PALETTE)
                    ? Renderer.COMPONENT_PAIR_PALETTE : null;
                var atomIdPairs = RDT.deriveComponentPairs(result, pal);
                var compPairs = (result && result.componentPairs) ? result.componentPairs : [];
                rendererPairs = atomIdPairs.slice();
                for (var ci = 0; ci < rendererPairs.length; ci++) {
                    var entry = rendererPairs[ci];
                    var match = null;
                    for (var cj = 0; cj < compPairs.length; cj++) {
                        if (compPairs[cj].reactantCompIdx === entry.reactantComponentIdx &&
                            compPairs[cj].productCompIdx === entry.productComponentIdx) {
                            match = compPairs[cj];
                            break;
                        }
                    }
                    if (match) { entry.paletteIndex = match.paletteIndex; }
                    else { entry.paletteIndex = ci; }
                }
                this.molecule.componentPairs = compPairs;
            } catch (cpe) {
                if (typeof console !== 'undefined' && console.warn) {
                    console.warn('BIME: RDT.deriveComponentPairs failed', cpe);
                }
                rendererPairs = [];
            }
        }
        this.renderer.componentPairs = rendererPairs;
        this.renderer.showComponentPairs = true;
        // Sync the toolbar pair-toggle and colors-toggle button states.
        var ppBtn = this._toolbarButtons && this._toolbarButtons['togglepairs'];
        if (ppBtn) { ppBtn.style.opacity = '1'; }
        var caBtn0 = this._toolbarButtons && this._toolbarButtons['togglecolors'];
        if (caBtn0) { caBtn0.style.opacity = this.renderer.colorAtoms ? '1' : '0.4'; }

        this.render();
        var nMapped = 0;
        for (var i = 0; i < this.molecule.atoms.length; i++) {
            if (this.molecule.atoms[i].mapNumber > 0) { nMapped++; }
        }
        var nChanges = (result && result.bondChanges) ? result.bondChanges.length : 0;
        // v1.4.1 (bug-fix #26): the status line distinguishes
        //   "components paired"  — how many reactant molecules pair with how
        //                          many product molecules (whole-molecule
        //                          mapping coverage),
        // from
        //   "rigid sub-fragments" — how many connected mapped sub-graphs
        //                           survived through the reaction without
        //                           any bond change (the colour groups the
        //                           user actually sees on screen).
        // The two numbers usually differ — and disagreeing is the *signal*
        // that the chemistry has bond-change events that subdivided a
        // component pair into more than one rigid piece.
        var nPairs = (result && result.componentPairs)
            ? result.componentPairs.filter(function(p) { return p.paletteIndex >= 0; }).length
            : 0;
        // v1.4.0: surface the normalised confidence in the status line.
        var conf = (result && typeof result.confidence === 'number') ? result.confidence : 0;
        var msg = 'Auto-map · ' + (result.strategy || 'no winner') +
            ' · ' + nPairs + ' comp pair' + (nPairs === 1 ? '' : 's') +
            ' · ' + nSubFrags + ' sub-frag' + (nSubFrags === 1 ? '' : 's') +
            ' · ' + nMapped + ' atoms mapped · ' +
            nChanges + ' bond change' + (nChanges === 1 ? '' : 's') +
            ' · Confidence = ' + conf.toFixed(5) +
            ' · ' + elapsed + ' ms';
        if (result.timedOut) { msg += ' (timed out)'; }
        if (result.warnings && result.warnings.length > 0) {
            msg += ' [' + result.warnings.join('; ') + ']';
        }
        this.showInfo(msg);
        this._fireCallback('AfterStructureModified');
    };

    /**
     * v1.4.0: toggle per-atom halo rendering (the "Color atoms" UI checkbox).
     * Returns the new state. Leaves componentPairs intact so toggling does
     * not lose the pair data — only suppresses the halos visually.
     */
    MolEditor.prototype.toggleColorAtoms = function() {
        if (!this.renderer) return false;
        this.renderer.colorAtoms = !this.renderer.colorAtoms;
        this.render();
        return this.renderer.colorAtoms;
    };

    /**
     * v1.4.0: clear the mol-mol component-pair overlay. Called when the user
     * presses the "Clear" toolbar button or loads a new molecule. Atom maps
     * themselves are not affected — only the visual pair boxes.
     */
    MolEditor.prototype.clearComponentPairs = function() {
        if (this.renderer) {
            this.renderer.componentPairs = [];
            this.render();
        }
    };

    /**
     * v1.4.0: toggle visibility of the mol-mol component-pair overlay
     * without dropping the underlying pair list. Returns the new state.
     */
    MolEditor.prototype.toggleComponentPairs = function() {
        if (!this.renderer) return false;
        this.renderer.showComponentPairs = !this.renderer.showComponentPairs;
        this.render();
        return this.renderer.showComponentPairs;
    };

    /**
     * Compute Tanimoto similarity between the current molecule and another.
     * Runs entirely in-browser using the bundled SMSDBatch path-fingerprint
     * implementation; no server is contacted.
     * @param {string} smiles2 - SMILES of the second molecule
     * @returns {Promise<Object>} {tanimoto, smiles1, smiles2}
     */
    MolEditor.prototype.similarity = function(smiles2) {
        var smiles1 = this.smiles();
        if (!smiles1) return Promise.reject(new Error('No molecule in editor'));
        if (!smiles2)  return Promise.reject(new Error('Second SMILES required'));

        var SG = (typeof SMSDGraph !== 'undefined') ? (SMSDGraph.SMSDGraph || SMSDGraph) : null;
        var Batch = (typeof SMSDBatch !== 'undefined') ? SMSDBatch : null;
        if (!SG || !Batch || !Batch.generateFingerprint || !Batch.tanimotoFingerprint) {
            return Promise.reject(new Error('Fingerprint similarity not available'));
        }
        try {
            var mol2 = new Molecule();
            SmilesParser.parse(smiles2, mol2);
            var g1 = new SG(this.molecule);
            var g2 = new SG(mol2);
            var fp1 = Batch.generateFingerprint(g1);
            var fp2 = Batch.generateFingerprint(g2);
            var tani = Batch.tanimotoFingerprint(fp1, fp2);
            return Promise.resolve({ tanimoto: tani, smiles1: smiles1, smiles2: smiles2 });
        } catch (e) {
            return Promise.reject(new Error('Similarity failed: ' + (e.message || e)));
        }
    };

    /**
     * Show MCS prompt inline in the status bar.
     * Prompts for a second SMILES and highlights the common substructure
     * using the bundled in-browser SMSDMCS module.
     */
    MolEditor.prototype._promptMCS = function() {
        var self = this;
        var bar = this._statusBar;

        bar.innerHTML = '<div style="display:flex;align-items:center;gap:0.4rem;flex:1">' +
            '<select id="bime-mcs-mode" style="padding:1px 4px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;font-size:10px;background:var(--color-bg,white);color:var(--color-text,#333);cursor:pointer">' +
            '<option value="default">MCS (all)</option>' +
            '<option value="ring">Ring only</option>' +
            '<option value="chain">Chain only</option>' +
            '<option value="connected">Connected</option>' +
            '<option value="disconnected">Disconnected</option>' +
            '<option value="similarity">Similarity</option>' +
            '</select>' +
            '<input id="bime-mcs-input" type="text" placeholder="Second SMILES (e.g. c1ccccc1O) ..." ' +
            'style="flex:1;padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;font-size:11px;font-family:var(--font-mono,monospace);background:var(--color-bg,white);color:var(--color-text,#333)">' +
            '<button id="bime-mcs-btn" style="padding:2px 8px;border:1px solid ' + BRAND_COLOR + ';border-radius:4px;background:' + BRAND_COLOR + ';color:white;font-size:11px;font-weight:600;cursor:pointer;white-space:nowrap">Find MCS</button>' +
            '<span id="bime-mcs-result" style="font-size:11px;font-weight:600"></span>' +
            '<button id="bime-mcs-export" style="display:none;padding:2px 8px;border:1px solid ' + BRAND_COLOR + ';border-radius:4px;background:none;color:' + BRAND_COLOR + ';font-size:10px;font-weight:600;cursor:pointer;white-space:nowrap" title="Export MCS with highlights">Export MCS \u25be</button>' +
            '<button id="bime-mcs-clear" style="padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:10px;font-weight:600" title="Clear highlights">Clear</button>' +
            '<button id="bime-mcs-close" style="padding:2px 6px;border:none;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:13px" title="Close">&times;</button>' +
            '</div>';

        var modeSelect = bar.querySelector('#bime-mcs-mode');
        var input = bar.querySelector('#bime-mcs-input');
        var btn = bar.querySelector('#bime-mcs-btn');
        var result = bar.querySelector('#bime-mcs-result');
        var exportBtn = bar.querySelector('#bime-mcs-export');
        var clearBtn = bar.querySelector('#bime-mcs-clear');
        var closeBtn = bar.querySelector('#bime-mcs-close');

        // Update button label and placeholder when mode changes
        modeSelect.addEventListener('change', function() {
            if (modeSelect.value === 'similarity') {
                btn.textContent = 'Search';
                input.placeholder = 'Leave empty to search library, or enter SMILES...';
            } else {
                btn.textContent = 'Find MCS';
                input.placeholder = 'Second SMILES (e.g. c1ccccc1O) ...';
            }
        });

        function doMCS() {
            var mode = modeSelect.value;
            var smiles1 = self.smiles();
            if (!smiles1) {
                result.style.color = '#dc2626';
                result.textContent = 'Draw a molecule first';
                exportBtn.style.display = 'none';
                return;
            }

            // Similarity search mode
            if (mode === 'similarity') {
                result.style.color = '#ca8a04';
                result.textContent = 'Searching...';
                exportBtn.style.display = 'none';
                try {
                    self._doSimilaritySearch(smiles1, input.value.trim(), result, exportBtn);
                } catch(e) {
                    result.style.color = '#dc2626';
                    result.textContent = 'Error: ' + (e.message || 'unknown');
                }
                return;
            }

            // MCS modes
            var smiles2 = input.value.trim();
            if (!smiles2) {
                result.style.color = '#dc2626';
                result.textContent = 'Enter a SMILES string';
                exportBtn.style.display = 'none';
                return;
            }

            // Build options based on mode
            var chemOpts = { matchAtomType: true, matchBondType: true };
            var mcsOpts = {};
            if (mode === 'ring') {
                chemOpts.ringMatchesRingOnly = true;
            } else if (mode === 'chain') {
                chemOpts.ringMatchesRingOnly = true;
                mcsOpts._chainOnly = true;
            } else if (mode === 'connected') {
                mcsOpts.connectedOnly = true;
            } else if (mode === 'disconnected') {
                mcsOpts.connectedOnly = false;
                mcsOpts.disconnectedMCS = true;
            }

            result.style.color = '#ca8a04';
            result.textContent = 'Computing MCS...';
            exportBtn.style.display = 'none';

            self.findMCS(smiles2, chemOpts, mcsOpts)
                .then(function(data) {
                    var atomCount = data._matchedAtoms || (data.atoms1 ? data.atoms1.length : 0);
                    var tani = data.tanimoto !== undefined ? data.tanimoto : '';
                    var taniStr = tani !== '' ? ', Tanimoto: ' + (typeof tani === 'number' ? tani.toFixed(2) : tani) : '';
                    result.style.color = '#059669';
                    result.textContent = 'MCS: ' + atomCount + ' atoms' + taniStr;
                    if (data.mcs) {
                        result.title = 'MCS SMILES: ' + data.mcs;
                    }
                    exportBtn.style.display = 'inline-flex';
                })
                .catch(function(err) {
                    result.style.color = '#dc2626';
                    result.textContent = 'MCS error: ' + (err.message || 'unknown');
                });
        }

        // Export MCS dropdown menu
        function showExportMenu() {
            // Remove existing menu
            var old = bar.querySelector('#bime-mcs-export-menu');
            if (old) { old.remove(); return; }

            var menu = document.createElement('div');
            menu.id = 'bime-mcs-export-menu';
            menu.style.cssText = 'position:absolute;bottom:100%;right:auto;z-index:200;background:white;border:1px solid var(--color-border,#e2e8f0);border-radius:6px;box-shadow:0 4px 16px rgba(0,0,0,0.12);padding:4px 0;min-width:140px;';
            var items = [
                { label: 'SVG', format: 'svg' },
                { label: 'PNG', format: 'png' },
                { label: 'Copy MCS SMILES', format: 'smiles' }
            ];
            items.forEach(function(item) {
                var opt = document.createElement('button');
                opt.style.cssText = 'display:block;width:100%;padding:6px 14px;border:none;background:none;text-align:left;cursor:pointer;font-size:11px;font-weight:500;color:#334155;';
                opt.textContent = item.label;
                opt.addEventListener('mouseenter', function() { opt.style.background = '#f0fdfa'; });
                opt.addEventListener('mouseleave', function() { opt.style.background = 'none'; });
                opt.addEventListener('click', function() {
                    menu.remove();
                    self.exportMCS(item.format);
                });
                menu.appendChild(opt);
            });
            exportBtn.style.position = 'relative';
            exportBtn.appendChild(menu);

            // Close menu on outside click
            function closeMenu(e) {
                if (!menu.contains(e.target) && e.target !== exportBtn) {
                    menu.remove();
                    document.removeEventListener('click', closeMenu);
                }
            }
            setTimeout(function() { document.addEventListener('click', closeMenu); }, 0);
        }

        function doClear() {
            self.clearHighlights();
            self.render();
            result.textContent = '';
            exportBtn.style.display = 'none';
            self._lastMCSResult = null;
        }

        function doClose() {
            self.clearHighlights();
            self.render();
            self._lastMCSResult = null;
            self._resetStatusBar();
        }

        btn.addEventListener('click', doMCS);
        exportBtn.addEventListener('click', showExportMenu);
        clearBtn.addEventListener('click', doClear);
        closeBtn.addEventListener('click', doClose);
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') doMCS();
            if (e.key === 'Escape') doClose();
        });
        input.focus();
    };

    /**
     * Back-compat stub. The Java SMSD server was removed in v1.0.3; MCS now
     * runs entirely in the browser via the SMSDMCS module. Always resolves
     * with true if the in-browser module is loaded.
     * @returns {Promise<boolean>}
     */
    MolEditor.prototype.checkSmsdServer = function() {
        var ok = (typeof SMSDMCS !== 'undefined' && typeof SMSDGraph !== 'undefined');
        return Promise.resolve(!!ok);
    };

    // =========================================================================
    // MCS Highlighting + Export
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /**
     * Clear all bgColor highlights from atoms and bonds.
     * Removes MCS, SMARTS, and any custom highlighting.
     */
    MolEditor.prototype.clearHighlights = function() {
        for (var i = 0; i < this.molecule.atoms.length; i++) {
            this.molecule.atoms[i].highlighted = false;
            this.molecule.atoms[i].bgColor = null;
            this.molecule.atoms[i].mcsMapNumber = 0;
        }
        for (var i = 0; i < this.molecule.bonds.length; i++) {
            this.molecule.bonds[i].highlighted = false;
            this.molecule.bonds[i].bgColor = null;
        }
    };

    /**
     * Export the current molecule with MCS highlights.
     * @param {string} format - 'svg', 'png', or 'smiles'
     * @returns {string|undefined} MCS SMILES when format='smiles', otherwise triggers download
     */
    MolEditor.prototype.exportMCS = function(format) {
        var self = this;
        var data = this._lastMCSResult;

        if (format === 'smiles') {
            var mcsSmiles = data ? (data.mcs || '') : '';
            // Copy to clipboard if available
            if (mcsSmiles && navigator.clipboard) {
                navigator.clipboard.writeText(mcsSmiles).then(function() {
                    self.showInfo('MCS SMILES copied: ' + mcsSmiles);
                }).catch(function() {
                    self.showInfo('MCS SMILES: ' + mcsSmiles);
                });
            } else {
                this.showInfo(mcsSmiles ? 'MCS SMILES: ' + mcsSmiles : 'No MCS result');
            }
            return mcsSmiles;
        }

        if (this.molecule.isEmpty()) {
            this.showInfo('No molecule to export');
            return;
        }

        // Build filename: "mol1_vs_mol2_mcs.ext"
        var name = this._mcsExportFilename(format);

        if (format === 'svg') {
            if (typeof ImageExport !== 'undefined') {
                var svg = ImageExport.toSVG(this.molecule, this._exportOptions());
                ImageExport.downloadSVG(svg, name);
                this.showInfo('Downloaded ' + name);
            } else {
                this.downloadSVG(name);
            }
        } else if (format === 'png') {
            if (typeof ImageExport !== 'undefined') {
                var opts = this._exportOptions();
                opts.scale = 2;
                this.showInfo('Exporting MCS PNG...');
                ImageExport.toPNG(this.molecule, opts).then(function(blob) {
                    ImageExport.downloadPNG(blob, name);
                    self.showInfo('Downloaded ' + name + ' (2x)');
                }).catch(function(err) {
                    self.showInfo('PNG export failed: ' + err.message);
                });
            } else {
                this.downloadPNG(name, 2);
            }
        }
    };

    /**
     * Generate a filename for MCS export: "mol1_vs_mol2_mcs.ext"
     * Uses molecule names when available, falls back to truncated SMILES.
     * @param {string} ext - file extension ('svg' or 'png')
     * @returns {string}
     */
    MolEditor.prototype._mcsExportFilename = function(ext) {
        var name1 = this.moleculeName || '';
        var name2 = '';

        // Try to get a human name for the query molecule
        if (this._lastMCSQuery) {
            // Check common molecules DB
            if (typeof COMMON_MOLECULES !== 'undefined' && COMMON_MOLECULES) {
                for (var i = 0; i < COMMON_MOLECULES.length; i++) {
                    if (COMMON_MOLECULES[i].smiles === this._lastMCSQuery) {
                        name2 = COMMON_MOLECULES[i].name;
                        break;
                    }
                }
            }
            if (!name2) {
                var q = this._lastMCSQuery;
                name2 = q.length <= 20 ? q.replace(/[^a-zA-Z0-9]/g, '_') : 'mol2';
            }
        }

        if (!name1) {
            var smi = this.smiles();
            name1 = (smi && smi.length <= 20) ? smi.replace(/[^a-zA-Z0-9]/g, '_') : 'mol1';
        }

        // Sanitise and truncate
        name1 = name1.replace(/[^a-zA-Z0-9_\-]/g, '_').substring(0, 24);
        name2 = name2.replace(/[^a-zA-Z0-9_\-]/g, '_').substring(0, 24);

        return name1 + '_vs_' + name2 + '_mcs.' + ext;
    };

    // =========================================================================
    // IUPAC Name Lookup via NCI Chemical Identifier Resolver
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /** Stored molecule name from last lookup */
    MolEditor.prototype.moleculeName = null;

    // -------------------------------------------------------------------------
    // IUPAC / common-name lookup
    //
    // BIME v1.0.0 forwarded structures to https://cactus.nci.nih.gov, exfiltrating
    // every drawn SMILES to a third-party server. v1.0.1 removes that channel
    // entirely for privacy and aligns BIME with a "zero external network calls"
    // posture by default.
    //
    // Public methods are kept so existing callers do not crash; they now resolve
    // with empty / no-op values and a clear message instead of issuing a network
    // request.
    // -------------------------------------------------------------------------

    var _NAME_LOOKUP_DISABLED_MSG =
        'IUPAC / common-name lookup is not available in this build.';

    /** @returns {Promise<string>} always rejects in the free tier. */
    MolEditor.prototype.getIUPACName = function() {
        return Promise.reject(new Error(_NAME_LOOKUP_DISABLED_MSG));
    };

    /** @returns {Promise<string[]>} always resolves with an empty list. */
    MolEditor.prototype.getCommonNames = function() {
        return Promise.resolve([]);
    };

    /** Toolbar handler \u2014 informs the user the feature is not available. */
    MolEditor.prototype._lookupIUPACName = function() {
        this._showNameStatus(_NAME_LOOKUP_DISABLED_MSG, '#6b7280');
    };

    /**
     * Display a molecule name in the status bar with a close button.
     */
    MolEditor.prototype._showNameStatus = function(msg, color) {
        var self = this;
        var bar = this._statusBar;
        bar.innerHTML = '';
        var span = document.createElement('span');
        span.style.cssText = 'font-size:11px;font-weight:600;color:' + (color || '#059669') + ';overflow:hidden;text-overflow:ellipsis;white-space:nowrap;flex:1';
        span.title = msg;
        span.textContent = msg;
        bar.appendChild(span);

        var copyBtn = document.createElement('button');
        copyBtn.style.cssText = 'padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;background:none;color:' + BRAND_COLOR + ';cursor:pointer;font-size:10px;font-weight:600;margin-left:0.5rem;white-space:nowrap';
        copyBtn.textContent = 'Copy';
        copyBtn.title = 'Copy name to clipboard';
        copyBtn.addEventListener('click', function() {
            if (navigator.clipboard) {
                navigator.clipboard.writeText(msg).then(function() {
                    copyBtn.textContent = 'Copied!';
                    setTimeout(function() { copyBtn.textContent = 'Copy'; }, 1500);
                });
            }
        });
        bar.appendChild(copyBtn);

        var closeBtn = document.createElement('button');
        closeBtn.style.cssText = 'border:none;background:none;color:var(--color-text-muted);cursor:pointer;font-size:13px;margin-left:0.3rem';
        closeBtn.title = 'Close';
        closeBtn.textContent = '\u00d7';
        closeBtn.addEventListener('click', function() {
            self._resetStatusBar();
        });
        bar.appendChild(closeBtn);
    };

    // =========================================================================
    // SVG Export
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /**
     * Get the SVG string of the current molecular drawing.
     * @returns {string} SVG markup
     */
    MolEditor.prototype.exportSVG = function() {
        return this.getMolecularAreaGraphicsString();
    };

    /**
     * Trigger a browser download of the molecule as an SVG file.
     * @param {string} [filename='molecule.svg'] download filename
     */
    MolEditor.prototype.downloadSVG = function(filename) {
        filename = filename || 'molecule.svg';
        var svgStr = this.exportSVG();
        if (!svgStr) {
            this.showInfo('No molecule to export');
            return;
        }
        var blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.style.display = 'none';
        document.body.appendChild(a);
        a.click();
        setTimeout(function() {
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }, 100);
        this.showInfo('Downloaded ' + filename);
    };

    // =========================================================================
    // PNG Export
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /**
     * Trigger a browser download of the molecule as a PNG file.
     * @param {string} [filename='molecule.png'] download filename
     * @param {number} [scale=2] resolution multiplier (1, 2, or 4)
     */
    MolEditor.prototype.downloadPNG = function(filename, scale) {
        filename = filename || 'molecule.png';
        scale = scale || 2;
        if (scale < 1) scale = 1;
        if (scale > 4) scale = 4;

        var svgStr = this.exportSVG();
        if (!svgStr) {
            this.showInfo('No molecule to export');
            return;
        }

        var self = this;

        // Parse SVG to get dimensions
        var parser = new DOMParser();
        var svgDoc = parser.parseFromString(svgStr, 'image/svg+xml');
        var svgEl = svgDoc.documentElement;
        var w = parseFloat(svgEl.getAttribute('width')) || this.renderer.width || 500;
        var h = parseFloat(svgEl.getAttribute('height')) || this.renderer.height || 400;

        var canvas = document.createElement('canvas');
        canvas.width = w * scale;
        canvas.height = h * scale;
        var ctx = canvas.getContext('2d');

        // White background
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        // Render SVG to canvas via Image
        var img = new Image();
        var svgBlob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
        var svgUrl = URL.createObjectURL(svgBlob);

        img.onload = function() {
            ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
            URL.revokeObjectURL(svgUrl);

            // Export as PNG
            canvas.toBlob(function(blob) {
                if (!blob) {
                    self.showInfo('PNG export failed');
                    return;
                }
                var pngUrl = URL.createObjectURL(blob);
                var a = document.createElement('a');
                a.href = pngUrl;
                a.download = filename;
                a.style.display = 'none';
                document.body.appendChild(a);
                a.click();
                setTimeout(function() {
                    document.body.removeChild(a);
                    URL.revokeObjectURL(pngUrl);
                }, 100);
                self.showInfo('Downloaded ' + filename + ' (' + scale + 'x)');
            }, 'image/png');
        };

        img.onerror = function() {
            URL.revokeObjectURL(svgUrl);
            self.showInfo('PNG export failed — SVG rendering error');
        };

        img.src = svgUrl;
    };

    // =========================================================================
    // ImageExport integration
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /**
     * Get default export options based on current editor state.
     */
    MolEditor.prototype._exportOptions = function() {
        return {
            width: this.renderer.width || 500,
            height: this.renderer.height || 400,
            title: this.moleculeName || '',
            showHydrogens: this.renderer.showHydrogens,
            watermark: true
        };
    };

    /**
     * Export SVG using ImageExport for standalone, publication-quality output.
     */
    MolEditor.prototype._exportSVGWithImageExport = function() {
        if (typeof ImageExport === 'undefined') {
            this.downloadSVG('molecule.svg');
            return;
        }
        if (this.molecule.isEmpty()) {
            this.showInfo('No molecule to export');
            return;
        }
        var svg = ImageExport.toSVG(this.molecule, this._exportOptions());
        var name = this._exportFilename('svg');
        ImageExport.downloadSVG(svg, name);
        this.showInfo('Downloaded ' + name);
    };

    /**
     * Export PNG using ImageExport for high-resolution output.
     */
    MolEditor.prototype._exportPNGWithImageExport = function() {
        var self = this;
        if (typeof ImageExport === 'undefined') {
            self.downloadPNG('molecule.png', 2);
            return;
        }
        if (self.molecule.isEmpty()) {
            self.showInfo('No molecule to export');
            return;
        }
        var opts = self._exportOptions();
        opts.scale = 2;
        var name = self._exportFilename('png');
        self.showInfo('Exporting PNG...');
        ImageExport.toPNG(self.molecule, opts).then(function(blob) {
            ImageExport.downloadPNG(blob, name);
            self.showInfo('Downloaded ' + name + ' (2x)');
        }).catch(function(err) {
            self.showInfo('PNG export failed: ' + err.message);
        });
    };

    /**
     * Export print-ready SVG using ImageExport.
     */
    MolEditor.prototype._exportPrintSVG = function() {
        if (typeof ImageExport === 'undefined') {
            this.showInfo('ImageExport not loaded');
            return;
        }
        if (this.molecule.isEmpty()) {
            this.showInfo('No molecule to export');
            return;
        }
        var svg = ImageExport.toPrintSVG(this.molecule, {
            width: this.renderer.width || 500,
            height: this.renderer.height || 400,
            title: this.moleculeName || ''
        });
        var name = this._exportFilename('print.svg');
        ImageExport.downloadSVG(svg, name);
        this.showInfo('Downloaded ' + name + ' (print quality)');
    };

    /**
     * Copy molecule as PNG to clipboard using ImageExport.
     */
    MolEditor.prototype._copyToClipboard = function() {
        var self = this;
        if (typeof ImageExport === 'undefined') {
            self.showInfo('ImageExport not loaded');
            return;
        }
        if (self.molecule.isEmpty()) {
            self.showInfo('No molecule to copy');
            return;
        }
        self.showInfo('Copying...');
        ImageExport.copyToClipboard(self.molecule, self._exportOptions(), self._statusBar)
            .then(function() {
                self.showInfo('Copied to clipboard!');
            })
            .catch(function(err) {
                self.showInfo('Copy failed: ' + err.message);
            });
    };

    /**
     * Generate a filename based on molecule name or SMILES.
     */
    MolEditor.prototype._exportFilename = function(ext) {
        var base = this.moleculeName || '';
        if (!base) {
            var smi = this.smiles();
            if (smi && smi.length <= 30) {
                base = smi.replace(/[^a-zA-Z0-9]/g, '_');
            } else {
                base = 'molecule';
            }
        }
        base = base.replace(/[^a-zA-Z0-9_\-]/g, '_').substring(0, 48);
        return base + '.' + ext;
    };

    // =========================================================================
    // Molecule Name Editing, Auto-naming & API
    // Copyright (c) 2026 BioInception PVT LTD. Apache License 2.0.
    // =========================================================================

    /**
     * Get the current molecule name.
     * @returns {string}
     */
    MolEditor.prototype.getMolName = function() {
        return this.molecule.name || '';
    };

    /**
     * Set the molecule name and re-render.
     * @param {string} name
     */
    MolEditor.prototype.setMolName = function(name) {
        this.molecule.name = name || '';
        this._updateNameStatus();
        this.render();
    };

    /**
     * Try to auto-name the molecule: first from common-molecules DB, then NCI Cactus.
     * @returns {Promise<string>} resolves with the name found (or empty string)
     */
    MolEditor.prototype.autoName = function() {
        var self = this;
        var smi = this.smiles();
        if (!smi) return Promise.resolve('');

        // 1. Try local common-molecules database
        var localName = this._lookupCommonName(smi);
        if (localName) {
            this.molecule.name = localName;
            this._updateNameStatus();
            this.render();
            return Promise.resolve(localName);
        }

        // 2. External name lookup is disabled in BIME v1.0.1 for privacy.
        return this._cactusNameLookup(smi).then(function(name) {
            if (name) {
                self.molecule.name = name;
                self._updateNameStatus();
                self.render();
            }
            return name || '';
        });
    };

    // =========================================================================
    // v1.5.0: unified library search.
    //
    // Single dispatch over the 1181-molecule COMMON_MOLECULES library for all
    // four search modes:
    //
    //   exact         — canonical SMILES equality (case-sensitive after
    //                   round-trip through SmilesWriter.write()).
    //   substructure  — query is a SMARTS subgraph of each library entry
    //                   (VF2++ via SmartsMatch.matchSmarts).
    //   mcs           — Maximum Common Substructure size between query and
    //                   each library entry (Bron-Kerbosch via SMSDMCS).
    //                   Score = mcsSize / max(|query|, |entry|).
    //   similarity    — Tanimoto over 1024-bit path fingerprints (existing
    //                   _doSimilaritySearch behaviour, lifted into the API).
    //
    // The library is parsed and indexed on the first call; subsequent calls
    // reuse the cache. Each library entry caches its Molecule, SMSDGraph,
    // canonical SMILES, and fingerprint.
    //
    // Returns a Promise that resolves to an array of hit objects sorted by
    // score descending, capped at options.topN (default 10):
    //
    //   [{ name, smiles, score, mode, mcsSize?, matchCount? }, ...]
    //
    // Promise-based so a future async backend (Web Worker, remote service)
    // can swap in without breaking callers.
    // =========================================================================

    var LIBRARY_CACHE = null; // populated on first searchLibrary call

    function _buildLibraryCache() {
        if (LIBRARY_CACHE) return LIBRARY_CACHE;
        LIBRARY_CACHE = [];
        if (typeof COMMON_MOLECULES === 'undefined' || !COMMON_MOLECULES) {
            return LIBRARY_CACHE;
        }
        if (typeof SmilesParser === 'undefined' || typeof SmilesWriter === 'undefined') {
            return LIBRARY_CACHE;
        }
        var SG = (typeof SMSDGraph !== 'undefined') ? (SMSDGraph.SMSDGraph || SMSDGraph) : null;
        var Batch = (typeof SMSDBatch !== 'undefined') ? SMSDBatch : null;
        for (var i = 0; i < COMMON_MOLECULES.length; i++) {
            var entry = COMMON_MOLECULES[i];
            try {
                var m = new Molecule();
                SmilesParser.parse(entry.smiles, m);
                var canon = '';
                try { canon = SmilesWriter.write(m); } catch (cw) { canon = entry.smiles; }
                var graph = null, fp = null;
                if (SG) {
                    try { graph = new SG(m); } catch (ge) { graph = null; }
                }
                if (graph && Batch && Batch.generateFingerprint) {
                    try { fp = Batch.generateFingerprint(graph); } catch (fe) { fp = null; }
                }
                LIBRARY_CACHE.push({
                    name: entry.name,
                    smiles: entry.smiles,
                    category: entry.category || 'basic',
                    molecule: m,
                    graph: graph,
                    canonical: canon,
                    fingerprint: fp,
                    heavyAtoms: m.atoms.length
                });
            } catch (e) { /* skip invalid entry */ }
        }
        return LIBRARY_CACHE;
    }

    /**
     * Reset the library cache. Useful for tests that mutate COMMON_MOLECULES
     * between runs.
     */
    MolEditor.resetLibraryCache = function() { LIBRARY_CACHE = null; };

    var USER_LIBRARY_LIMIT = 100;

    /**
     * Build a one-off cache for a user-supplied list of [{name, smiles}] entries.
     * Returns the same shape as _buildLibraryCache() so the search loops below
     * don't need a separate code path. Entries with parse errors are silently
     * dropped. The list is capped at USER_LIBRARY_LIMIT (100) to keep MCS
     * search responsive even on the slowest browsers.
     */
    function _buildUserLibraryCache(targets) {
        var out = [];
        if (!Array.isArray(targets) || targets.length === 0) { return out; }
        if (typeof SmilesParser === 'undefined' || typeof SmilesWriter === 'undefined') {
            return out;
        }
        var SG = (typeof SMSDGraph !== 'undefined') ? (SMSDGraph.SMSDGraph || SMSDGraph) : null;
        var Batch = (typeof SMSDBatch !== 'undefined') ? SMSDBatch : null;
        var n = Math.min(targets.length, USER_LIBRARY_LIMIT);
        for (var i = 0; i < n; i++) {
            var entry = targets[i];
            if (!entry || !entry.smiles) { continue; }
            try {
                var m = new Molecule();
                SmilesParser.parse(entry.smiles, m);
                var canon = '';
                try { canon = SmilesWriter.write(m); } catch (cw) { canon = entry.smiles; }
                var graph = null, fp = null;
                if (SG) { try { graph = new SG(m); } catch (ge) { graph = null; } }
                if (graph && Batch && Batch.generateFingerprint) {
                    try { fp = Batch.generateFingerprint(graph); } catch (fe) { fp = null; }
                }
                out.push({
                    name: entry.name || ('entry-' + (i + 1)),
                    smiles: entry.smiles,
                    category: entry.category || 'user',
                    molecule: m,
                    graph: graph,
                    canonical: canon,
                    fingerprint: fp,
                    heavyAtoms: m.atoms.length
                });
            } catch (e) { /* skip invalid entry */ }
        }
        return out;
    }

    /**
     * Maximum number of molecules a user-supplied library may contain.
     * Excess entries beyond this index are silently dropped at parse time.
     * Exposed so the workbench UI can warn users at upload time.
     */
    MolEditor.USER_LIBRARY_LIMIT = USER_LIBRARY_LIMIT;

    /**
     * Search the COMMON_MOLECULES library (or a user-supplied target set)
     * against a query SMILES (or the current canvas if `query` is omitted).
     *
     * @param {string|null} query   SMILES string; null/undefined = use canvas
     * @param {string} mode         'exact' | 'substructure' | 'mcs' | 'similarity'
     * @param {Object} [options]
     *        options.topN          max hits to return (default 10)
     *        options.threshold     similarity / mcs minimum score (default 0.3)
     *        options.timeoutMs     hard cap (default 10000) — only honoured
     *                              by MCS mode (the slow one)
     *        options.targets       (v1.5.0+) user-supplied list of
     *                              [{name, smiles}, ...] to search instead of
     *                              the built-in COMMON_MOLECULES library.
     *                              Capped at MolEditor.USER_LIBRARY_LIMIT (100)
     *                              entries; excess silently dropped to keep
     *                              MCS search responsive.
     * @returns {Promise<Array>}    [{name, smiles, score, mode, ...}, ...]
     */
    MolEditor.prototype.searchLibrary = function(query, mode, options) {
        var self = this;
        options = options || {};
        var topN = options.topN || 10;
        var threshold = (options.threshold !== undefined) ? options.threshold : 0.3;
        var timeoutMs = options.timeoutMs || 10000;
        var modeKey = (mode || 'similarity').toLowerCase();
        if (modeKey !== 'exact' && modeKey !== 'substructure' &&
            modeKey !== 'mcs' && modeKey !== 'similarity') {
            return Promise.reject(new Error('searchLibrary: unknown mode ' + modeKey));
        }

        // Resolve the query SMILES.
        var qSmi = (typeof query === 'string' && query.length > 0) ? query : this.smiles();
        if (!qSmi) { return Promise.reject(new Error('searchLibrary: no query SMILES')); }

        // Pick the target set. options.targets > built-in COMMON_MOLECULES.
        var library;
        if (Array.isArray(options.targets)) {
            library = _buildUserLibraryCache(options.targets);
            if (!library || library.length === 0) {
                return Promise.reject(new Error('searchLibrary: user-supplied targets empty or unparsable'));
            }
        } else {
            library = _buildLibraryCache();
            if (!library || library.length === 0) {
                return Promise.reject(new Error('searchLibrary: COMMON_MOLECULES library not loaded'));
            }
        }

        // Parse the query into a Molecule + graph + canonical + fingerprint.
        var queryMol = new Molecule();
        try { SmilesParser.parse(qSmi, queryMol); }
        catch (pe) { return Promise.reject(new Error('searchLibrary: bad query SMILES: ' + (pe.message || pe))); }
        var SG = (typeof SMSDGraph !== 'undefined') ? (SMSDGraph.SMSDGraph || SMSDGraph) : null;
        var qGraph = SG ? (function() { try { return new SG(queryMol); } catch (e) { return null; } })() : null;
        var qCanon = '';
        try { qCanon = SmilesWriter.write(queryMol); } catch (cw) { qCanon = qSmi; }
        var Batch = (typeof SMSDBatch !== 'undefined') ? SMSDBatch : null;
        var qFp = (qGraph && Batch && Batch.generateFingerprint)
            ? (function() { try { return Batch.generateFingerprint(qGraph); } catch (e) { return null; } })()
            : null;

        var hits = [];
        var t0 = Date.now();

        if (modeKey === 'exact') {
            for (var i = 0; i < library.length; i++) {
                var le = library[i];
                if (le.canonical && le.canonical === qCanon) {
                    hits.push({ name: le.name, smiles: le.smiles, score: 1.0, mode: 'exact' });
                }
            }
        } else if (modeKey === 'substructure') {
            // VF2 subgraph isomorphism: is the query (treated as SMARTS) a
            // substructure of each library entry? Treat the canonical SMILES
            // as the SMARTS pattern (works for non-aromatic queries; aromatic
            // queries gain the standard SMARTS aromatic-perception semantics).
            if (typeof SmartsMatch === 'undefined' || !SmartsMatch.matchSmarts) {
                return Promise.reject(new Error('searchLibrary: SmartsMatch module not loaded'));
            }
            var smartsQuery = qCanon || qSmi;
            for (var j = 0; j < library.length; j++) {
                var lj = library[j];
                try {
                    var matches = SmartsMatch.matchSmarts(lj.molecule, smartsQuery);
                    if (matches && matches.length > 0) {
                        hits.push({
                            name: lj.name,
                            smiles: lj.smiles,
                            score: matches.length / Math.max(lj.heavyAtoms, 1),
                            mode: 'substructure',
                            matchCount: matches.length
                        });
                    }
                } catch (e) { /* skip target on parse error */ }
            }
            // Substructure: sort by (matchCount desc, smiles.length asc, name asc).
            // The trailing name comparison gives a strict total order so two
            // hits that tie on matchCount AND smiles.length stay in a
            // reproducible order independent of source iteration order.
            hits.sort(function(a, b) {
                if (b.matchCount !== a.matchCount) { return b.matchCount - a.matchCount; }
                if (a.smiles.length !== b.smiles.length) { return a.smiles.length - b.smiles.length; }
                return (a.name || '').localeCompare(b.name || '');
            });
        } else if (modeKey === 'mcs') {
            if (typeof SMSDMCS === 'undefined' || !SMSDMCS.findMCS || !qGraph) {
                return Promise.reject(new Error('searchLibrary: SMSDMCS module not loaded'));
            }
            for (var k = 0; k < library.length; k++) {
                if ((Date.now() - t0) > timeoutMs) { break; }
                var lk = library[k];
                if (!lk.graph) { continue; }
                try {
                    var mcs = SMSDMCS.findMCS(qGraph, lk.graph, undefined, { timeoutMs: 1000 });
                    var size = (mcs && mcs.size) ? mcs.size : 0;
                    var denom = Math.max(qGraph.n || 1, lk.heavyAtoms || 1);
                    var score = size / denom;
                    if (score >= threshold) {
                        hits.push({
                            name: lk.name,
                            smiles: lk.smiles,
                            score: score,
                            mode: 'mcs',
                            mcsSize: size
                        });
                    }
                } catch (e) { /* skip on MCS error */ }
            }
            hits.sort(function(a, b) { return b.score - a.score; });
        } else { // similarity
            if (!qFp) {
                return Promise.reject(new Error('searchLibrary: fingerprint not available'));
            }
            for (var l = 0; l < library.length; l++) {
                var ll = library[l];
                if (!ll.fingerprint) { continue; }
                try {
                    var tani = Batch.tanimotoFingerprint(qFp, ll.fingerprint);
                    if (tani >= threshold) {
                        hits.push({
                            name: ll.name,
                            smiles: ll.smiles,
                            score: tani,
                            mode: 'similarity'
                        });
                    }
                } catch (e) { /* skip */ }
            }
            hits.sort(function(a, b) { return b.score - a.score; });
        }

        // Cap at topN; resolve.
        return Promise.resolve(hits.slice(0, topN));
    };

    /**
     * v1.5.0: highlight the part of the currently-loaded molecule that
     * matched a search hit. Called by the workbench after readGenericMolecularInput
     * loads a hit's SMILES so the user sees which atoms made it match.
     *
     * @param {string} mode   'substructure' | 'mcs' | 'exact' | 'similarity'
     * @param {string} querySmiles  the original search query
     * @returns {number}      count of atoms highlighted (0 = no highlight applied)
     *
     * Behaviour by mode:
     *   substructure  — runs SmartsMatch.highlightMatches with the query as the
     *                    SMARTS pattern; sets atom.highlighted on every atom in
     *                    the first match. Returns the matched-atom count.
     *   mcs           — runs SMSDMCS.findMCS between the query and the loaded
     *                    molecule, then highlights MCS atoms via bgColor (uses
     *                    the existing _applyMCSHighlights pipeline).
     *   exact         — sets atom.highlighted on every atom (whole-molecule
     *                    match — every atom contributed).
     *   similarity    — no-op (similarity is a whole-molecule property; per-
     *                    atom contributions to a Tanimoto score aren't
     *                    well-defined for a path fingerprint and would
     *                    mislead more than help).
     */
    MolEditor.prototype.highlightSearchMatch = function(mode, querySmiles) {
        if (!this.molecule || !this.molecule.atoms) { return 0; }
        var modeKey = (mode || '').toLowerCase();
        // v1.5.2 finesse: timestamp the new highlight pass so the Renderer
        // can fire a one-shot reveal animation only on freshly matched
        // atoms/bonds (not on every re-render). The `searchMatchKind` flag
        // lets the renderer paint a richer ring/halo without changing
        // existing bgColor semantics.
        var nowTs = (typeof Date !== 'undefined' && Date.now) ? Date.now() : 0;
        // Reset prior highlights first.
        for (var i = 0; i < this.molecule.atoms.length; i++) {
            var a = this.molecule.atoms[i];
            a.highlighted = false;
            a.bgColor = null;
            a.searchMatchKind = null;
            a.searchMatchTs = 0;
        }
        for (var b = 0; b < this.molecule.bonds.length; b++) {
            this.molecule.bonds[b].highlighted = false;
            this.molecule.bonds[b].bgColor = null;
            this.molecule.bonds[b].searchMatchKind = null;
            this.molecule.bonds[b].searchMatchTs = 0;
        }
        var n = 0;
        if (modeKey === 'similarity') {
            // No per-atom highlight semantics for path-fingerprint Tanimoto;
            // ranking by score on the search panel is the right signal.
            this.render();
            return 0;
        }
        if (modeKey === 'exact') {
            for (var k = 0; k < this.molecule.atoms.length; k++) {
                var aEx = this.molecule.atoms[k];
                aEx.highlighted = true;
                aEx.searchMatchKind = 'exact';
                aEx.searchMatchTs = nowTs;
                n++;
            }
            // Also flag every bond so the bond ring/glow lights up.
            for (var bEx = 0; bEx < this.molecule.bonds.length; bEx++) {
                this.molecule.bonds[bEx].searchMatchKind = 'exact';
                this.molecule.bonds[bEx].searchMatchTs = nowTs;
                this.molecule.bonds[bEx].highlighted = true;
            }
            this.render();
            return n;
        }
        if (modeKey === 'substructure') {
            if (typeof SmartsMatch === 'undefined' || !SmartsMatch.matchSmarts) {
                this.render();
                return 0;
            }
            try {
                var matches = SmartsMatch.matchSmarts(this.molecule, querySmiles);
                if (matches && matches.length > 0) {
                    var first = matches[0];
                    for (var qId in first) {
                        if (!first.hasOwnProperty(qId)) { continue; }
                        var atomId = first[qId];
                        var atom = this.molecule.getAtom(atomId);
                        if (atom) {
                            atom.highlighted = true;
                            atom.bgColor = '#99f6e4';
                            atom.searchMatchKind = 'substructure';
                            atom.searchMatchTs = nowTs;
                            n++;
                        }
                    }
                    // Highlight bonds between matched atoms too.
                    var matchedSet = {};
                    for (var qId2 in first) {
                        if (first.hasOwnProperty(qId2)) { matchedSet[first[qId2]] = true; }
                    }
                    for (var bi = 0; bi < this.molecule.bonds.length; bi++) {
                        var bond = this.molecule.bonds[bi];
                        if (matchedSet[bond.atom1] && matchedSet[bond.atom2]) {
                            bond.bgColor = '#99f6e4';
                            bond.highlighted = true;
                            bond.searchMatchKind = 'substructure';
                            bond.searchMatchTs = nowTs;
                        }
                    }
                }
            } catch (e) { /* leave highlights cleared on error */ }
            this.render();
            return n;
        }
        if (modeKey === 'mcs') {
            if (typeof SMSDMCS === 'undefined' || !SMSDMCS.findMCS ||
                typeof SMSDGraph === 'undefined') {
                this.render();
                return 0;
            }
            try {
                var queryMol = new Molecule();
                SmilesParser.parse(querySmiles, queryMol);
                var SG = SMSDGraph.SMSDGraph || SMSDGraph;
                var gQ = new SG(queryMol);
                var gT = new SG(this.molecule);
                var mcs = SMSDMCS.findMCS(gQ, gT, undefined, { timeoutMs: 2000 });
                var rawMapping = (mcs && mcs.mapping) ? mcs.mapping : {};
                var atomIdMapping = (typeof SMSDMCS.translateToAtomIds === 'function')
                    ? SMSDMCS.translateToAtomIds(rawMapping, gQ, gT)
                    : rawMapping;
                var matchedSetMcs = {};
                var mapCounter = 1;
                for (var qk in atomIdMapping) {
                    if (!atomIdMapping.hasOwnProperty(qk)) { continue; }
                    var tId = +atomIdMapping[qk];
                    var tAtom = this.molecule.getAtom(tId);
                    if (tAtom) {
                        tAtom.bgColor = '#99f6e4';
                        tAtom.mcsMapNumber = mapCounter++;
                        tAtom.highlighted = true;
                        tAtom.searchMatchKind = 'mcs';
                        tAtom.searchMatchTs = nowTs;
                        matchedSetMcs[tId] = true;
                        n++;
                    }
                }
                for (var bj = 0; bj < this.molecule.bonds.length; bj++) {
                    var bnd = this.molecule.bonds[bj];
                    if (matchedSetMcs[bnd.atom1] && matchedSetMcs[bnd.atom2]) {
                        bnd.bgColor = '#99f6e4';
                        bnd.highlighted = true;
                        bnd.searchMatchKind = 'mcs';
                        bnd.searchMatchTs = nowTs;
                    }
                }
            } catch (e) { /* leave highlights cleared on error */ }
            this.render();
            return n;
        }
        this.render();
        return 0;
    };

    /**
     * Look up a SMILES in the COMMON_MOLECULES database (if loaded).
     * Compares canonical SMILES for an exact match.
     * @param {string} inputSmiles
     * @returns {string|null} molecule name or null
     */
    MolEditor.prototype._lookupCommonName = function(inputSmiles) {
        if (typeof COMMON_MOLECULES === 'undefined' || !COMMON_MOLECULES) return null;
        if (!inputSmiles) return null;

        // Build a cached canonical→name index on first call
        if (!MolEditor._cmIndex) {
            MolEditor._cmIndex = {};         // canonical SMILES → name
            MolEditor._cmFpIndex = {};       // "formula|bonds" fingerprint → [name, ...]
            if (typeof SmilesWriter !== 'undefined' && typeof SmilesParser !== 'undefined') {
                for (var k = 0; k < COMMON_MOLECULES.length; k++) {
                    var ent = COMMON_MOLECULES[k];
                    try {
                        var m = new Molecule();
                        SmilesParser.parse(ent.smiles, m);
                        var canon = SmilesWriter.write(m);
                        MolEditor._cmIndex[canon] = ent.name;
                        var fp = _molFingerprint(m);
                        if (!MolEditor._cmFpIndex[fp]) MolEditor._cmFpIndex[fp] = [];
                        MolEditor._cmFpIndex[fp].push({ name: ent.name, smiles: ent.smiles });
                    } catch(e) { /* skip bad entries */ }
                }
            }
        }

        // Canonicalise the input
        var canonInput = inputSmiles.trim();
        var inputMol = null;
        if (typeof SmilesWriter !== 'undefined' && typeof SmilesParser !== 'undefined') {
            try {
                inputMol = new Molecule();
                SmilesParser.parse(canonInput, inputMol);
                canonInput = SmilesWriter.write(inputMol);
            } catch (e) { inputMol = null; }
        }

        // Pass 1: exact canonical match
        if (MolEditor._cmIndex[canonInput]) return MolEditor._cmIndex[canonInput];

        // Pass 2: fingerprint match (formula + bond count) — handles Kekulé vs aromatic mismatch
        if (inputMol) {
            var fp2 = _molFingerprint(inputMol);
            var candidates = MolEditor._cmFpIndex[fp2];
            if (candidates && candidates.length === 1) return candidates[0].name;
            // Multiple candidates: disambiguate by atom count
            if (candidates && candidates.length > 1) {
                for (var c = 0; c < candidates.length; c++) {
                    try {
                        var cm = new Molecule();
                        SmilesParser.parse(candidates[c].smiles, cm);
                        if (cm.atoms.length === inputMol.atoms.length) return candidates[c].name;
                    } catch(e) { /* skip */ }
                }
            }
        }

        return null;
    };

    /** Compute a lightweight fingerprint string for a molecule: "formula|heavyBonds" */
    function _molFingerprint(mol) {
        var counts = {};
        for (var i = 0; i < mol.atoms.length; i++) {
            var sym = (mol.atoms[i] && (mol.atoms[i].symbol || mol.atoms[i].atom)) || 'C';
            counts[sym] = (counts[sym] || 0) + 1;
        }
        var keys = Object.keys(counts).sort(function(a, b) {
            var o = {'C': 0, 'H': 1};
            return ((o[a] !== undefined ? o[a] : 2) - (o[b] !== undefined ? o[b] : 2)) || a.localeCompare(b);
        });
        var formula = keys.map(function(k){ return k + (counts[k] > 1 ? counts[k] : ''); }).join('');
        return formula + '|' + mol.bonds.length;
    }

    /**
     * No-op stub. Was previously a fire-and-forget call to NCI Cactus that
     * exfiltrated drawn SMILES on every load. Kept as a public-API-compatible
     * stub; equivalent functionality may be provided by a future component.
     */
    MolEditor.prototype._asyncCactusName = function() { /* disabled in v1.0.1 */ };

    /**
     * No-op stub (was: NCI Cactus IUPAC lookup). Always resolves with ''.
     * Public callers (e.g. third-party integrations of older BIME) keep
     * working without errors.
     * @returns {Promise<string>}
     */
    MolEditor.prototype._cactusNameLookup = function(/* smi */) {
        return Promise.resolve('');
    };

    /**
     * Toggle molecule name visibility (click) or look up name if none set.
     * Called from the "Name" toolbar button action.
     */
    MolEditor.prototype._toggleOrEditMolName = function() {
        if (this.molecule.name) {
            // Toggle visibility
            this.renderer.showMolName = !this.renderer.showMolName;
            this.render();
            this.showInfo(this.renderer.showMolName ? 'Name shown: ' + this.molecule.name : 'Name hidden');
        } else {
            // No name set yet — do IUPAC lookup (legacy behaviour)
            this._lookupIUPACName();
        }
    };

    /**
     * Update the status bar to show the molecule name (if set).
     */
    MolEditor.prototype._updateNameStatus = function() {
        if (this.molecule.name) {
            this.showInfo(this.molecule.name);
        } else {
            this.showInfo('');
        }
    };

    /**
     * Open an inline text input in the status bar for editing the molecule name.
     */
    MolEditor.prototype._editMolNameInline = function() {
        var self = this;
        var bar = this._statusBar;

        bar.innerHTML = '<div style="display:flex;align-items:center;gap:0.4rem;flex:1">' +
            '<span style="font-size:9px;font-weight:700;letter-spacing:0.04em;text-transform:uppercase;color:var(--color-text-muted,#94a3b8);white-space:nowrap">Name</span>' +
            '<input id="bime-name-input" type="text" placeholder="Enter molecule name..." ' +
            'value="' + escapeHTML(this.molecule.name || '') + '" ' +
            'style="flex:1;padding:2px 6px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;font-size:12px;font-weight:600;background:var(--color-bg,white);color:var(--color-text,#333)">' +
            '<button id="bime-name-ok" style="padding:2px 8px;border:1px solid ' + BRAND_COLOR + ';border-radius:4px;background:' + BRAND_COLOR + ';color:white;font-size:11px;font-weight:600;cursor:pointer">Set</button>' +
            '<button id="bime-name-auto" style="padding:2px 8px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;background:none;color:' + BRAND_COLOR + ';font-size:11px;font-weight:600;cursor:pointer" title="Auto-lookup name">Auto</button>' +
            '<button id="bime-name-close" style="padding:2px 6px;border:none;background:none;color:var(--color-text-muted,#64748b);cursor:pointer;font-size:13px" title="Close">&times;</button>' +
            '</div>';

        var input = bar.querySelector('#bime-name-input');
        var okBtn = bar.querySelector('#bime-name-ok');
        var autoBtn = bar.querySelector('#bime-name-auto');
        var closeBtn = bar.querySelector('#bime-name-close');

        function doSet() {
            self.molecule.name = input.value.trim();
            self.renderer.showMolName = true;
            self.render();
            doClose();
        }

        function doClose() {
            self._resetStatusBar();
            self._updateNameStatus();
        }

        function doAuto() {
            input.value = 'Looking up...';
            self.autoName().then(function(name) {
                input.value = name || '(not found)';
            });
        }

        okBtn.addEventListener('click', doSet);
        autoBtn.addEventListener('click', doAuto);
        closeBtn.addEventListener('click', doClose);
        input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') doSet();
            if (e.key === 'Escape') doClose();
        });
        input.focus();
        input.select();
    };

    // Stub methods (log warning on first call)
    var _stubs = ['setAtomBackgroundColors', 'resetAtomColors', 'setBondBackgroundColors', 'resetBondColors',
        'setAtomToHighLight', 'setBondToHighLight', 'setTemplate', 'setSubstituent', 'addClickHandler',
        'setMenuXShortcuts', 'setMolecularAreaScale', 'setAtomMolecularAreaFontSize',
        'setMolecularAreaLineWidth', 'setAntialias', 'setCopyToClipboardFormat',
        'setMarkerMenuBackGroundColorPalette', 'setAction', 'setAtomAdditionalData',
        'getAtomAdditionalData', 'setBondAdditionalData', 'getBondAdditionalData',
        'nonisomericSmiles', 'getCreationIndex', 'getParentContainer',
        'setUserInterfaceBackgroundColor', 'setVisible', 'isVisible'];

    _stubs.forEach(function(name) {
        if (!MolEditor.prototype[name]) {
            MolEditor.prototype[name] = function() {
                // Silent stub
                return null;
            };
        }
    });

    // =========================================================================
    // v1.4.3: Customize panel — let the user show/hide groups + buttons,
    // reorder them by drag-and-drop (or up/down arrows as fallback), and
    // pick which atom-bar elements appear. Preferences persist to
    // localStorage via ToolbarPrefs.save() and re-apply on every editor
    // build.
    // =========================================================================

    /**
     * Open the Customize panel as a modal overlay. The panel is rebuilt
     * from scratch on each open so it always reflects the latest prefs.
     */
    MolEditor.prototype._openCustomizePanel = function() {
        if (this._customizeOverlay) { return; }
        var self = this;
        var hasPrefs = (typeof ToolbarPrefs !== 'undefined');
        if (!hasPrefs) {
            this.showInfo('Customize: ToolbarPrefs module not loaded');
            return;
        }
        var current = ToolbarPrefs.load() || ToolbarPrefs.snapshot(TOOLBAR_GROUPS, ATOM_BAR);
        // Working copy that gets mutated as the user toggles things.
        var draft = JSON.parse(JSON.stringify(current));

        // Build overlay backdrop.
        var overlay = document.createElement('div');
        overlay.className = 'bime-customize-overlay';
        overlay.setAttribute('role', 'dialog');
        overlay.setAttribute('aria-modal', 'true');
        overlay.setAttribute('aria-label', 'Customize toolbar');
        overlay.style.cssText = 'position:fixed;inset:0;background:rgba(0,0,0,0.45);z-index:9999;display:flex;align-items:center;justify-content:center;padding:1rem;';

        var card = document.createElement('div');
        card.className = 'bime-customize-card';
        card.style.cssText = 'background:var(--color-bg,white);border:1px solid var(--color-border,#e2e8f0);border-radius:8px;box-shadow:0 12px 36px rgba(0,0,0,0.20);max-width:560px;width:100%;max-height:85vh;display:flex;flex-direction:column;overflow:hidden;';

        // Header
        var head = document.createElement('div');
        head.style.cssText = 'padding:12px 16px;border-bottom:1px solid var(--color-border,#e2e8f0);display:flex;align-items:center;justify-content:space-between;background:var(--color-surface,#f8fafc);';
        var title = document.createElement('h3');
        title.style.cssText = 'margin:0;font-size:14px;font-weight:600;color:var(--color-text,#334155);';
        title.textContent = 'Customize toolbar';
        var hint = document.createElement('span');
        hint.style.cssText = 'font-size:11px;color:var(--color-text-muted,#94a3b8);';
        hint.textContent = 'drag rows to reorder · click checkboxes to show/hide';
        head.appendChild(title);
        head.appendChild(hint);
        card.appendChild(head);

        // Scrollable body
        var body = document.createElement('div');
        body.style.cssText = 'padding:12px 16px;overflow-y:auto;flex:1;font-size:12px;color:var(--color-text,#334155);';
        card.appendChild(body);

        self._renderCustomizeBody(body, draft);

        // Footer with Save / Reset / Cancel buttons.
        var foot = document.createElement('div');
        foot.style.cssText = 'padding:10px 16px;border-top:1px solid var(--color-border,#e2e8f0);display:flex;gap:8px;justify-content:flex-end;background:var(--color-surface,#f8fafc);';
        function mkBtn(label, primary) {
            var b = document.createElement('button');
            b.type = 'button';
            b.textContent = label;
            b.style.cssText = 'padding:6px 14px;border-radius:5px;font-size:12px;font-weight:600;cursor:pointer;border:1px solid ' +
                (primary ? BRAND_COLOR : 'var(--color-border,#e2e8f0)') + ';background:' +
                (primary ? BRAND_COLOR : 'var(--color-bg,white)') + ';color:' +
                (primary ? 'white' : 'var(--color-text,#334155)') + ';';
            return b;
        }
        var resetBtn = mkBtn('Reset to defaults', false);
        var cancelBtn = mkBtn('Cancel', false);
        var saveBtn = mkBtn('Save', true);
        foot.appendChild(resetBtn);
        foot.appendChild(cancelBtn);
        foot.appendChild(saveBtn);
        card.appendChild(foot);

        overlay.appendChild(card);
        document.body.appendChild(overlay);
        this._customizeOverlay = overlay;

        // Wire buttons.
        var close = function() {
            if (self._customizeOverlay) {
                self._customizeOverlay.parentNode.removeChild(self._customizeOverlay);
                self._customizeOverlay = null;
            }
            if (self._customizeKeyHandler) {
                document.removeEventListener('keydown', self._customizeKeyHandler);
                self._customizeKeyHandler = null;
            }
        };
        cancelBtn.addEventListener('click', close);
        overlay.addEventListener('click', function(e) {
            if (e.target === overlay) { close(); }
        });
        self._customizeKeyHandler = function(e) {
            if (e.key === 'Escape') { e.preventDefault(); close(); }
        };
        document.addEventListener('keydown', self._customizeKeyHandler);

        resetBtn.addEventListener('click', function() {
            ToolbarPrefs.reset();
            draft = ToolbarPrefs.snapshot(TOOLBAR_GROUPS, ATOM_BAR);
            self._renderCustomizeBody(body, draft);
            self.showInfo('Toolbar customization reset to defaults');
        });
        saveBtn.addEventListener('click', function() {
            var ok = ToolbarPrefs.save(draft);
            if (!ok) {
                self.showInfo('Customize: could not save preferences (storage unavailable)');
                return;
            }
            close();
            self._rebuildToolbar();
            self.showInfo('Toolbar customization saved');
        });
    };

    /**
     * Render the editable rows for groups + items + atom-bar into a scroll
     * body. Every interaction mutates `draft` in place; the caller commits
     * by calling ToolbarPrefs.save(draft).
     */
    MolEditor.prototype._renderCustomizeBody = function(body, draft) {
        var self = this;
        body.innerHTML = '';

        // -------------- Groups + items section --------------
        var hGroups = document.createElement('div');
        hGroups.style.cssText = 'font-size:11px;font-weight:700;color:var(--color-text-muted,#94a3b8);text-transform:uppercase;letter-spacing:0.04em;margin-bottom:6px;';
        hGroups.textContent = 'Groups & buttons';
        body.appendChild(hGroups);

        var groupsList = document.createElement('div');
        groupsList.style.cssText = 'display:flex;flex-direction:column;gap:4px;margin-bottom:14px;';
        body.appendChild(groupsList);

        function rebuildGroups() {
            groupsList.innerHTML = '';
            draft.groupOrder.forEach(function(gid, gIdx) {
                var canon = null;
                for (var i = 0; i < TOOLBAR_GROUPS.length; i++) {
                    if (TOOLBAR_GROUPS[i].id === gid) { canon = TOOLBAR_GROUPS[i]; break; }
                }
                if (!canon) { return; }
                var pg = draft.groups[gid] || { hidden: false, items: [] };
                var row = self._mkCustomizeGroupRow(canon, pg, draft, gIdx, rebuildGroups);
                groupsList.appendChild(row);
            });
        }
        rebuildGroups();

        // -------------- Atom bar section --------------
        var hAtoms = document.createElement('div');
        hAtoms.style.cssText = 'font-size:11px;font-weight:700;color:var(--color-text-muted,#94a3b8);text-transform:uppercase;letter-spacing:0.04em;margin-bottom:6px;';
        hAtoms.textContent = 'Atom bar';
        body.appendChild(hAtoms);

        var atomList = document.createElement('div');
        atomList.style.cssText = 'display:flex;flex-wrap:wrap;gap:6px;';
        body.appendChild(atomList);

        function rebuildAtoms() {
            atomList.innerHTML = '';
            ATOM_BAR.forEach(function(sym) {
                var on = draft.atomBar.indexOf(sym) >= 0;
                var pill = document.createElement('button');
                pill.type = 'button';
                pill.textContent = sym;
                pill.title = on ? 'Click to hide ' + sym : 'Click to show ' + sym;
                pill.setAttribute('aria-pressed', on ? 'true' : 'false');
                pill.style.cssText = 'padding:4px 10px;font-size:12px;font-weight:700;border-radius:14px;cursor:pointer;border:1px solid ' +
                    (on ? BRAND_COLOR : 'var(--color-border,#e2e8f0)') + ';background:' +
                    (on ? BRAND_COLOR : 'var(--color-bg,white)') + ';color:' +
                    (on ? 'white' : 'var(--color-text,#334155)') + ';';
                pill.addEventListener('click', function() {
                    var idx = draft.atomBar.indexOf(sym);
                    if (idx >= 0) { draft.atomBar.splice(idx, 1); }
                    else {
                        // Insert at the position matching ATOM_BAR canonical order.
                        var canonIdx = ATOM_BAR.indexOf(sym);
                        var insertAt = draft.atomBar.length;
                        for (var k = 0; k < draft.atomBar.length; k++) {
                            if (ATOM_BAR.indexOf(draft.atomBar[k]) > canonIdx) {
                                insertAt = k; break;
                            }
                        }
                        draft.atomBar.splice(insertAt, 0, sym);
                    }
                    rebuildAtoms();
                });
                atomList.appendChild(pill);
            });
        }
        rebuildAtoms();
    };

    /**
     * Build a single row in the Customize panel for one toolbar group.
     * Drag-and-drop reorders the row; up/down arrows are an a11y fallback.
     * A summary line shows the buttons in the group with checkboxes for
     * show/hide.
     */
    MolEditor.prototype._mkCustomizeGroupRow = function(canon, pg, draft, gIdx, rebuild) {
        var self = this;
        var row = document.createElement('div');
        row.draggable = true;
        row.style.cssText = 'border:1px solid var(--color-border,#e2e8f0);border-radius:6px;padding:8px 10px;background:var(--color-bg,white);cursor:grab;';

        // Header row: drag handle + group name + visibility checkbox + reorder arrows.
        var head = document.createElement('div');
        head.style.cssText = 'display:flex;align-items:center;gap:8px;';
        var handle = document.createElement('span');
        handle.textContent = '☰'; // hamburger / drag handle glyph
        handle.style.cssText = 'color:var(--color-text-muted,#94a3b8);cursor:grab;user-select:none;font-size:12px;';
        handle.title = 'Drag to reorder';
        head.appendChild(handle);

        var name = document.createElement('strong');
        name.textContent = canon.label || canon.id;
        name.style.cssText = 'font-size:13px;font-weight:600;flex:1;';
        head.appendChild(name);

        // Visibility checkbox for the WHOLE group
        var vis = document.createElement('label');
        vis.style.cssText = 'display:inline-flex;align-items:center;gap:4px;font-size:11px;color:var(--color-text-muted,#94a3b8);cursor:pointer;';
        var visBox = document.createElement('input');
        visBox.type = 'checkbox';
        visBox.checked = !pg.hidden;
        visBox.addEventListener('change', function() {
            pg.hidden = !visBox.checked;
        });
        vis.appendChild(visBox);
        var visLbl = document.createElement('span');
        visLbl.textContent = 'visible';
        vis.appendChild(visLbl);
        head.appendChild(vis);

        function arrowBtn(text, fn, title) {
            var b = document.createElement('button');
            b.type = 'button';
            b.textContent = text;
            b.title = title;
            b.setAttribute('aria-label', title);
            b.style.cssText = 'width:22px;height:22px;border:1px solid var(--color-border,#e2e8f0);border-radius:4px;background:var(--color-bg,white);cursor:pointer;font-size:11px;color:var(--color-text,#334155);';
            b.addEventListener('click', fn);
            return b;
        }
        head.appendChild(arrowBtn('↑', function() {
            if (gIdx > 0) {
                var s = draft.groupOrder.splice(gIdx, 1)[0];
                draft.groupOrder.splice(gIdx - 1, 0, s);
                rebuild();
            }
        }, 'Move group up'));
        head.appendChild(arrowBtn('↓', function() {
            if (gIdx < draft.groupOrder.length - 1) {
                var s = draft.groupOrder.splice(gIdx, 1)[0];
                draft.groupOrder.splice(gIdx + 1, 0, s);
                rebuild();
            }
        }, 'Move group down'));
        row.appendChild(head);

        // Body: per-button checkboxes + reorder arrows.
        var bodyDiv = document.createElement('div');
        bodyDiv.style.cssText = 'margin-top:6px;display:flex;flex-direction:column;gap:3px;padding-left:18px;';
        function rebuildItems() {
            bodyDiv.innerHTML = '';
            if (!Array.isArray(pg.hiddenItems)) { pg.hiddenItems = []; }
            (canon.items || []).forEach(function(item) {
                var on = (pg.hiddenItems.indexOf(item.id) < 0);
                var line = document.createElement('div');
                line.style.cssText = 'display:flex;align-items:center;gap:6px;font-size:11px;';
                var lbl = document.createElement('label');
                lbl.style.cssText = 'display:inline-flex;align-items:center;gap:4px;flex:1;cursor:pointer;';
                var cb = document.createElement('input');
                cb.type = 'checkbox';
                cb.checked = on;
                cb.addEventListener('change', function() {
                    if (cb.checked) {
                        // Visible: ensure in items + remove from hiddenItems.
                        var hk = pg.hiddenItems.indexOf(item.id);
                        if (hk >= 0) { pg.hiddenItems.splice(hk, 1); }
                        if (pg.items.indexOf(item.id) < 0) { pg.items.push(item.id); }
                    } else {
                        // Hidden: add to hiddenItems + remove from items.
                        if (pg.hiddenItems.indexOf(item.id) < 0) { pg.hiddenItems.push(item.id); }
                        var k = pg.items.indexOf(item.id);
                        if (k >= 0) { pg.items.splice(k, 1); }
                    }
                    rebuildItems();
                });
                lbl.appendChild(cb);
                var span = document.createElement('span');
                span.textContent = item.label || item.id;
                lbl.appendChild(span);
                line.appendChild(lbl);
                if (on) {
                    line.appendChild(arrowBtn('↑', function() {
                        var k = pg.items.indexOf(item.id);
                        if (k > 0) {
                            var s = pg.items.splice(k, 1)[0];
                            pg.items.splice(k - 1, 0, s);
                            rebuildItems();
                        }
                    }, 'Move button up'));
                    line.appendChild(arrowBtn('↓', function() {
                        var k = pg.items.indexOf(item.id);
                        if (k >= 0 && k < pg.items.length - 1) {
                            var s = pg.items.splice(k, 1)[0];
                            pg.items.splice(k + 1, 0, s);
                            rebuildItems();
                        }
                    }, 'Move button down'));
                }
                bodyDiv.appendChild(line);
            });
        }
        rebuildItems();
        row.appendChild(bodyDiv);

        // Drag-and-drop reordering of group rows.
        row.addEventListener('dragstart', function(e) {
            row.style.opacity = '0.5';
            try { e.dataTransfer.setData('text/plain', String(gIdx)); } catch (err) {}
            try { e.dataTransfer.effectAllowed = 'move'; } catch (err) {}
        });
        row.addEventListener('dragend', function() { row.style.opacity = '1'; });
        row.addEventListener('dragover', function(e) { e.preventDefault(); try { e.dataTransfer.dropEffect = 'move'; } catch (err) {} });
        row.addEventListener('drop', function(e) {
            e.preventDefault();
            row.style.opacity = '1';
            var srcIdx;
            try { srcIdx = parseInt(e.dataTransfer.getData('text/plain'), 10); } catch (err) { srcIdx = -1; }
            if (isNaN(srcIdx) || srcIdx < 0 || srcIdx === gIdx) { return; }
            var s = draft.groupOrder.splice(srcIdx, 1)[0];
            draft.groupOrder.splice(gIdx, 0, s);
            rebuild();
        });

        return row;
    };

    /**
     * Tear down the existing toolbar DOM and rebuild from current
     * TOOLBAR_GROUPS + ATOM_BAR + saved prefs. Called after the user saves
     * a new layout in the Customize panel.
     */
    MolEditor.prototype._rebuildToolbar = function() {
        if (!this._toolbar || !this._atomBar) { return; }
        // Empty existing children. v1.5.2: also clear the left rail if present.
        while (this._toolbar.firstChild) { this._toolbar.removeChild(this._toolbar.firstChild); }
        if (this._leftRail) {
            while (this._leftRail.firstChild) { this._leftRail.removeChild(this._leftRail.firstChild); }
        }
        while (this._atomBar.firstChild) { this._atomBar.removeChild(this._atomBar.firstChild); }
        this._buildToolbar();
        this._buildAtomBar();
        if (this._currentToolName) {
            // Re-highlight the active tool button if present in the new layout.
            var btn = this._toolbarButtons[this._currentToolName];
            if (btn) { btn.classList.add('active'); btn.setAttribute('aria-pressed', 'true'); }
        }
    };

    // =========================================================================
    // Global compat shim
    // =========================================================================
    if (!global.JSApplet) global.JSApplet = {};
    global.JSApplet.JSME = MolEditor;
    global.MolEditor = MolEditor;

})(typeof window !== 'undefined' ? window : this);
