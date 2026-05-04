# BIME — User & Developer Guide

> **B**io**I**nception **M**olecular **E**ditor — a zero-dependency,
> Apache-2.0 in-browser chemistry editor for SMILES, SMARTS, MOL, RXN,
> reactions, atom-atom mapping, MCS, similarity, and image export.

This guide walks through every component of the BIME user interface and
its public JavaScript API. Skip to whichever section you need:

- [Quick start](#quick-start)
- [The workbench page](#the-workbench-page)
- [Editor toolbar](#editor-toolbar)
- [Atom bar](#atom-bar)
- [Drawing tools](#drawing-tools)
- [Reaction toolkit](#reaction-toolkit)
- [Customize toolbar (v1.4.3+)](#customize-toolbar-v143)
- [SMILES, SMARTS, MOL I/O](#smiles-smarts-mol-io)
- [Auto-mapping & sub-fragment colours (v1.4.0–v1.4.2)](#auto-mapping--sub-fragment-colours)
- [MCS, similarity & substructure search](#mcs-similarity--substructure-search)
- [Image export](#image-export)
- [Keyboard shortcuts](#keyboard-shortcuts)
- [Public API](#public-api)
- [Embedding BIME in your page](#embedding-bime-in-your-page)
- [ToolbarPrefs API (v1.4.3+)](#toolbarprefs-api-v143)
- [RDT API](#rdt-api)
- [Renderer API](#renderer-api)
- [Events & callbacks](#events--callbacks)
- [Accessibility](#accessibility)
- [Versioning & compat](#versioning--compat)

---

## Quick start

BIME is a single bundle, no build step, no npm. Three lines:

```html
<div id="bime-editor" style="width:100%;height:520px"></div>
<script src="dist/bime.min.js"
        integrity="sha384-..."          <!-- copy from dist/SRI.txt -->
        crossorigin="anonymous"></script>
<script>
  var editor = new MolEditor('bime-editor');
  editor.smiles('CC(=O)O.OCC>>CC(=O)OCC.O');   // ester reaction
</script>
```

Works in every browser that ships SVG + ES5 (Edge/Chrome/Firefox/Safari/iOS/Android).
No internet connection required after page load.

---

## The workbench page

`workbench.html` is the live demo and reference layout. It exposes every
feature BIME ships and is a useful starting template for your own pages.

```
┌─────────────────────────────────────────────────────────────┐
│  Toolbar  ←  rebuilt from TOOLBAR_GROUPS + user prefs       │
│  Atoms    ←  rebuilt from ATOM_BAR + user prefs             │
├─────────────────────────────────────────┬───────────────────┤
│                                         │                   │
│         Drawing canvas (SVG)            │  Search type      │
│                                         │  ◯ Exact match    │
│                                         │  ● Substructure   │
│                                         │  ◯ Similarity     │
│                                         │                   │
├─────────────────────────────────────────┴───────────────────┤
│  Status bar: SMILES, atom/bond counts, AAM info             │
├─────────────────────────────────────────────────────────────┤
│  SMILES input:  [_______________________________]  [Load]   │
│  Editor Output: read-only canonical SMILES of canvas        │
│  ┌────────────────────────────────────────────┐             │
│  │  Drag and drop .mol or .sdf files here     │             │
│  │  Or click to select a file                 │             │
│  └────────────────────────────────────────────┘             │
│  [           Paste from Clipboard                ]          │
├─────────────────────────────────────────────────────────────┤
│  Sample molecules · Molecule browser (1181 cats)            │
│  Import structure · Validate                                │
└─────────────────────────────────────────────────────────────┘
```

Open `workbench.html` directly in a browser. No server required for
local development; for production hosting any static-file CDN works.

---

## Editor toolbar

The toolbar is built dynamically from a constant called `TOOLBAR_GROUPS`
at the top of `editor/MolEditor.js`. Each group is a card holding related
buttons. Hover the **(i)** icon at the right edge of any labelled group
for a one-line description.

| Group        | Default contents | What it does |
|--------------|------------------|--------------|
| **Draw**     | Bond, Chain | Place new bonds (Bond) or extend chains (Chain). Bond is the default tool. |
| **Bonds**    | Single, Double, Triple, Wedge | Set the bond order or stereochemistry for new bonds. Click an existing bond to cycle order. |
| **Rings**    | 3, 4, 5, 6, 7 | Drop a ring template at the cursor. The Rings group has a flyout — click the button to expand. |
| **Edit**     | Select, Delete, Move, Undo, Redo, Clear | Selection lasso, eraser, drag-to-move, undo/redo (Ctrl/Cmd+Z, Ctrl/Cmd+Y), clear canvas. |
| **View**     | Layout, Zoom +, Zoom –, Fit, R/S, E/Z, H, Ar | Auto-layout, pan/zoom, toggle CIP stereo labels, show implicit Hs, switch aromatic-ring rendering between dashed-circle and Kekulé. |
| **Reaction** | Arrow, Map, Map #, Auto-map, Colors, Pairs | Draw reaction arrow, manually set atom-atom maps, toggle their display, run the RDT auto-mapper, toggle the per-atom RDT-style colour halos, toggle the molecule-pair overlay. |
| _(actions)_  | Validate, SMARTS, MCS, Sim | Floating action chips. Validate the structure against valence rules; build a SMARTS from selection; run MCS/similarity prompts. |
| **Export**   | Name, SVG, PNG, Print, Copy | Toggle/edit molecule name, save as SVG, save as PNG (up to 4× resolution), print, copy to clipboard. |
| **Custom**   | Customize | Open the Customize panel (v1.4.3+) to show/hide buttons, reorder them, pick atom-bar elements, save to localStorage. |

Buttons styled as **type='action'** are click-once actions; **type='tool'**
toggles the active drawing tool; **type='bond'** sets bond order;
**type='ring'** sets the active ring template.

---

## Atom bar

A horizontal bar of element pills below the toolbar. Click an element to
make it the active atom for the next bond/atom drawn. Default elements:

```
C N O S F Cl Br I P H
```

Click and HOLD any pill to open the **periodic-table popup** (full 18-column
layout with CPK colours). The popup is keyboard-navigable.

Customizable in v1.4.3+ via the Customize panel — pick exactly which
elements appear and in what order.

---

## Drawing tools

| Tool | How to use |
|------|-----------|
| **Bond** | Click two atoms to bond them, or click and drag from one atom to draw a new atom + bond. Click-drag on empty canvas to place a 2-atom chain. |
| **Chain** | Click-drag to extrude a chain of carbons; release to commit. |
| **Select** | Click-drag a lasso to select atoms. Click an atom to toggle. Selected atoms can be moved (Move tool), deleted (Delete key), or wrapped in SMARTS via the SMARTS toolbar action. |
| **Delete** | Click any atom or bond to remove it. Cascades hydrogens automatically. |
| **Move** | Click-drag any atom to reposition. The full molecule recomputes its bounds and re-renders. |

Templates and rings: pick a ring size from the **Rings** group. The
template snaps to the nearest atom under the cursor when clicked, or
drops freestanding at the click location.

---

## Reaction toolkit

To draw a reaction:

1. Draw the reactant molecules on the canvas (separate fragments using
   the dot operator in SMILES, or just leave space between drawings).
2. Pick **Arrow** from the Reaction group. Click and drag to place the
   reaction arrow between reactants and products.
3. Draw the product molecules to the right of the arrow.
4. Optional: pick **Map** to manually set atom-atom map numbers (click
   one reactant atom, then its product partner — they share a number).
5. Or run **Auto-map** (RDT) to compute mapping in one click.

After Auto-map, BIME paints **per-atom RDT-style halos** (v1.4.0+) so the
eye traces chemistry across the arrow:

- **v1.4.0** — whole-molecule pairs share a halo colour (one colour per
  reactant ↔ product component pair).
- **v1.4.1+** — per-MCS-sub-fragment colour groups: each rigid scaffold
  piece that survived a bond-change-free traversal of the reaction gets
  its own colour. A bisubstrate enzyme like dihydrostreptomycin + ATP
  shows three colours (streptomycin scaffold blue, adenosine green,
  γ-phosphate orange).
- **v1.4.2+** — leftover-atom rescue ensures otherwise-orphan atoms
  (e.g. lone Cl in benzene + Cl → chlorobenzene) also get mapped and
  coloured.

Toggle the colour overlay with **Colors**; toggle map-number display
with **Map #**; toggle the molecule-pair box overlay with **Pairs**.

The status bar after Auto-map shows:

```
Auto-map · STRATEGY · N comp pairs · M sub-frags · K atoms mapped
         · P bond change(s) · Confidence = N.NNNNN · ELAPSED ms
```

`comp pairs` = how many reactants pair with how many products.
`sub-frags` = how many rigid scaffold pieces survived without bond
              change — the colour groups visible on screen.

---

## Customize toolbar (v1.4.3+)

Click the **Customize** button in the Custom group (or call
`editor._openCustomizePanel()` from JS) to open the customization modal.

The panel lets you:

- **Show or hide whole groups** via the "visible" checkbox.
- **Show or hide individual buttons** via per-button checkboxes.
- **Reorder groups** by drag-and-drop (drag the **☰** handle on each row)
  or with the **↑** / **↓** arrow buttons.
- **Reorder buttons** within a group via the **↑** / **↓** arrows beside
  each visible button.
- **Pick which atom-bar elements appear** by clicking element pills at
  the bottom (toggle on/off).
- **Reset to defaults** — wipes saved prefs and restores the canonical
  layout.
- **Save** — writes prefs to `localStorage` under
  `bime-toolbar-prefs-v1`. Prefs reload on every editor build.

Forward-compat: when a future BIME version adds a new button, it appears
**by default** in the user's saved layout (appended at the end of its
group). Users explicitly hide buttons via the Customize panel; hidden
buttons are tracked separately so newly-added buttons aren't accidentally
hidden.

Press **Escape** or click the dimmed backdrop to dismiss without saving.

---

## SMILES, SMARTS, MOL I/O

```js
// Read SMILES into the editor
editor.smiles('CC(=O)O.OCC>>CC(=O)OCC.O');

// Get canonical SMILES from the canvas
var smi = editor.smiles();

// Read MOL/V2000/V3000/RXN/MRV (any auto-detected text format)
editor.readGenericMolecularInput(textBlob);

// Write a SMARTS from the current selection
editor.smarts();    // -> '[C][N]'
```

The workbench's **Drag-and-drop zone** routes file contents through
`editor.readGenericMolecularInput()`. **Paste from Clipboard** does the
same with `navigator.clipboard.readText()`.

Validation: click **Validate** in the actions row to run the structural
sanity check. Errors appear inline next to the button.

---

## Auto-mapping & sub-fragment colours

Public API:

```js
// Run RDT atom-atom mapping (synchronous; pure-JS)
var result = RDT.mapReaction(reaction, options);

// result fields:
//   mapping            : { reactantAtomId: productAtomId, ... }
//   bondChanges        : [{ type: 'formed' | 'cleaved' | 'orderChange'
//                                   | 'stereoChange' | 'hydrogenChange',
//                           reactantBond, productBond, atoms,
//                           beforeOrder, afterOrder }, ...]
//   strategy           : 'MIN' | 'MAX' | 'MIXTURE' | 'RING'
//   componentPairs     : [{ reactantCompIdx, productCompIdx, mcsSize,
//                           paletteIndex, rescued? }, ...]
//   confidence         : 0..1 heuristic quality score
//   sides              : { reactants: Molecule[], products: Molecule[] }
//   pairs, score, scoreAfterRescue (v1.4.2+), warnings, timedOut

// Per-MCS-sub-fragment halo groups (v1.4.1+)
var subFrags = RDT.deriveSubFragments(result, { minSize: 1 });
// -> [{ reactantAtomIds, productAtomIds, paletteIndex, size }, ...]

// options:
//   strategies            : ['MIN', 'MAX', 'MIXTURE', 'RING']
//   useBipartitePostPass  : true (default since v1.3.0)
//   useLeftoverRescue     : true (default since v1.4.2)
//   includeHydrogens      : true
//   includeStereo         : true
//   timeoutMs             : 10000
```

Pass `useLeftoverRescue: false` to restore v1.4.1 mapping counts (some
downstream pipelines key into specific atom counts and need byte-equivalent
output).

---

## MCS, similarity & substructure search

### v1.5.0+: unified library search

```js
// Single dispatch over the 1181-molecule built-in library OR a user-
// supplied set, with four search modes:
editor.searchLibrary(query, mode, options).then(function(hits) {
  // hits: [{ name, smiles, score, mode, mcsSize?, matchCount? }, ...]
});

// mode = 'exact' | 'substructure' | 'mcs' | 'similarity'
// options:
//   topN         (default 10)        cap on returned hits
//   threshold    (default 0.3)       similarity / mcs minimum score
//   timeoutMs    (default 10000)     hard cap (MCS only)
//   targets      (default null)      user-supplied [{name, smiles}, ...]
//                                    capped at MolEditor.USER_LIBRARY_LIMIT (100).
//                                    null = use built-in COMMON_MOLECULES.

// Examples:
editor.searchLibrary('c1ccccc1', 'substructure');           // 1181 lib
editor.searchLibrary('c1ccccc1', 'mcs', { topN: 5 });
editor.searchLibrary('CCO', 'similarity', { threshold: 0.5 });
editor.searchLibrary('CCO', 'exact', {                       // user lib
  targets: [
    { name: 'A', smiles: 'CCO' },
    { name: 'B', smiles: 'c1ccccc1' }
  ]
});
```

| Mode | Algorithm | What you get |
|---|---|---|
| `exact` | Canonical SMILES equality (post-`SmilesWriter.write`) | hits with `score = 1.0` |
| `substructure` | VF2++ subgraph isomorphism via `SmartsMatch.matchSmarts` | hits ranked by `matchCount` (number of subgraph matches) |
| `mcs` | Bron-Kerbosch + modular product (`SMSDMCS.findMCS`) | hits ranked by `mcsSize / max(|query|, |target|)`; carry `mcsSize` |
| `similarity` | Tanimoto over 1024-bit path fingerprint (`SMSDBatch`) | hits ranked by Tanimoto coefficient |

### Highlight matched atoms after loading a hit (v1.5.0+)

```js
// Substructure / MCS / Exact: paint the matched atoms on the canvas.
// Similarity: no-op (whole-molecule property, no per-atom contributions).
var matchedCount = editor.highlightSearchMatch(mode, querySmiles);
// Returns the number of atoms highlighted.
```

The workbench panel automatically calls `highlightSearchMatch` when the
user clicks a search hit, so the matched part stays visible after the
hit's full SMILES is loaded.

### User-supplied target library (max 100)

The workbench's "My set (max 100)" toggle lets the user drop a
`.smi` / `.smiles` / `.sdf` / `.txt` file holding their own molecules:

- `.smi` / `.smiles` / `.txt`: one molecule per line, format
  `SMILES <whitespace> name` (rest of line = name; missing → auto-numbered).
- `.sdf`: `$$$$`-separated MOL blocks; the parser extracts the SMILES
  from the first line that lexically looks like one and uses the first
  block line as the name.

Anything beyond the 100-entry cap is silently dropped to keep MCS
search responsive. The cap is exposed as
`MolEditor.USER_LIBRARY_LIMIT`.

### Legacy single-target prompts (still available)

```js
// Maximum Common Substructure between canvas and a single SMILES
editor.findMCS('c1ccccc1OC').then(function(data) {
  console.log(data.mcs, data.atoms1, data.atoms2, data.tanimoto);
});

// Tanimoto similarity (path-fingerprint, in-browser)
editor.similarity('c1ccccc1Cl').then(function(data) {
  console.log(data.tanimoto);
});

// Status-bar prompts (the "Sim" and "MCS" toolbar action buttons)
editor._promptMCS();
editor._promptSimSearch();
```

The workbench's **Search type** radio panel (Exact / Substructure / MCS /
Similarity) chooses the algorithm; the **Search in** source toggle
chooses the target set (built-in 1181 vs user-supplied ≤ 100).

---

## Image export

```js
// Save the canvas as SVG (Apache-2.0, all stylesheets inlined)
var svgString = editor.exportSVG();

// Save as PNG at the given device pixel ratio (1×, 2×, 4×)
editor.exportPNG(2).then(function(pngBlob) { ... });

// Print
editor.exportPrint();

// Copy SVG to clipboard
editor.copyToClipboard();
```

The exported SVG carries the BIME watermark and an
`<metadata>` element with the canonical SMILES + RDT mapping if present.

---

## Keyboard shortcuts

| Shortcut | Action |
|---|---|
| **Ctrl/Cmd + Z** | Undo |
| **Ctrl/Cmd + Y** or **Ctrl/Cmd + Shift + Z** | Redo |
| **Delete** / **Backspace** | Delete selected atoms (Select tool only) |
| **Escape** | Close any open popover / modal (incl. Customize panel) |
| **Tab** / **Shift + Tab** | Navigate toolbar buttons (focus rings always visible) |
| **Enter** / **Space** | Activate the focused button |
| **Ctrl/Cmd + V** | Paste SMILES/MOL text into the editor |

Touch devices: pinch-to-zoom on the canvas. Tap-and-hold an atom button
to open the periodic-table popup.

---

## Public API

The single global is `MolEditor` (also aliased as `JSApplet.JSME` for
legacy consumers). Other namespaces:

| Global | Purpose |
|---|---|
| `MolEditor` | Editor class |
| `Molecule` | Pure-JS molecule model |
| `SmilesParser` / `SmilesWriter` | Canonical SMILES I/O |
| `SmartsParser` / `SmartsMatch` | SMARTS parsing + matching |
| `RDT` | Reaction Decoder Tool — atom-atom mapping, bond-change |
| `SMSDGraph` / `SMSDMCS` / `SMSDVF2` / `SMSDRings` / `SMSDBatch` / `SMSDLayout` | Substructure / MCS / fingerprint algorithms |
| `Renderer` | SVG rendering layer |
| `Layout` | 2-D coordinate generation |
| `Templates` | Ring + scaffold templates (steroid, morphinan, pyranose, etc.) |
| `History` | Undo / redo stack |
| `ImageExport` | SVG / PNG / clipboard helpers |
| `CipStereo` | R/S / E/Z assignment |
| `ToolbarPrefs` | (v1.4.3+) localStorage-backed toolbar customization |

All are read-only after the bundle loads. None require initialization.

---

## Embedding BIME in your page

```html
<!DOCTYPE html>
<html>
<head><meta charset="utf-8"></head>
<body>
  <div id="my-editor" style="width:900px;height:520px"></div>

  <!-- Verified bundle, integrity-checked -->
  <script src="https://asad.github.io/bime-dist/dist/bime.min.js"
          integrity="sha384-..."
          crossorigin="anonymous"></script>

  <script>
    var ed = new MolEditor('my-editor');
    ed.options('depict');                    // read-only mode
    ed.smiles('OC1=CC=CC=C1');               // phenol
    ed.setCallBack('AfterStructureModified', function(e) {
      console.log('canvas changed; new SMILES:', ed.smiles());
    });
  </script>
</body>
</html>
```

Constructor:

```js
new MolEditor(idOrElement[, optionsString])
```

Where `optionsString` is a comma- or space-separated list of flags:

| Flag | Effect |
|---|---|
| `depict` | Read-only — no drawing, just rendering |
| `nohydrogens` | Hide implicit Hs |
| `nodepict` | (default) — full editing enabled |
| `minsubfragsize=N` | (v1.4.2+) Minimum sub-fragment size for halo colouring |

```js
ed.setSize('1200px', '600px');               // resize
ed.reset();                                   // wipe canvas
ed.destroy();                                 // tear down (v1.4.1+ removes ALL listeners)
```

---

## ToolbarPrefs API (v1.4.3+)

Pure persistence — no DOM. Useful if you want to ship a custom layout
with your page, or add your own UI to control the editor's appearance.

```js
// Load from localStorage (returns null if missing or corrupted)
var prefs = ToolbarPrefs.load();

// Snapshot the canonical default layout as a starting point
var snap = ToolbarPrefs.snapshot(MolEditor.TOOLBAR_GROUPS, MolEditor.ATOM_BAR);

// Mutate, then save
snap.groups.export.hidden = true;     // hide the whole Export group
snap.atomBar = ['C','N','O'];         // show only C, N, O
ToolbarPrefs.save(snap);              // returns true on success

// Apply prefs to canonical defaults (returns filtered/reordered copies)
var groups = ToolbarPrefs.applyToGroups(MolEditor.TOOLBAR_GROUPS, prefs);
var atoms  = ToolbarPrefs.applyToAtomBar(MolEditor.ATOM_BAR, prefs);

// Reset to defaults
ToolbarPrefs.reset();

// Validate an arbitrary object against the schema
ToolbarPrefs.validate(someJsonBlob);  // -> bool
```

Storage shape (v1 schema, key `bime-toolbar-prefs-v1`):

```json
{
  "version": 1,
  "groupOrder": ["draw", "bonds", "rings", "edit", "view",
                 "reaction", "actions", "export", "custom"],
  "groups": {
    "draw":  { "hidden": false, "items": ["bond", "chain"], "hiddenItems": [] },
    "bonds": { "hidden": false, "items": ["single", "double", "triple", "stereo"], "hiddenItems": [] }
  },
  "atomBar": ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "H"]
}
```

---

## RDT API

```js
// Run atom-atom mapping
var result = RDT.mapReaction(reaction, options);

// Inspect strategies
RDT.STRATEGIES;  // ['MIN', 'MAX', 'MIXTURE', 'RING']

// Lower-level annotation
var bcs = RDT.annotateBondChanges(reactantSide, productSide, atomMap, opts);

// Per-MCS-sub-fragment halo groups (v1.4.1+)
var subFrags = RDT.deriveSubFragments(result, { minSize: 1 });

// Per-component-pair halo groups (v1.4.0)
var compPairs = RDT.deriveComponentPairs(result, palette);

// Confidence heuristic
var conf = RDT.deriveConfidence(result);
```

See the file-level docstring at the top of `editor/RDT.js` for the full
algorithm description (strategies, fitness, bipartite post-pass,
leftover-atom rescue).

---

## Renderer API

The SVG renderer is exposed as `editor.renderer`. Key knobs:

```js
editor.renderer.showHydrogens   = true;
editor.renderer.showMapNumbers  = true;     // :N atom-atom maps
editor.renderer.showCipRS       = false;    // R/S labels
editor.renderer.showCipEZ       = false;    // E/Z labels
editor.renderer.aromaticStyle   = 'circle'; // or 'kekule'
editor.renderer.colorAtoms      = true;     // RDT halo overlay
editor.renderer.componentPairs  = subFrags; // halo groups
editor.renderer.scale           = 1.0;      // zoom
```

After mutating any flag, call `editor.render()` to re-paint.

---

## Events & callbacks

```js
editor.setCallBack('AfterStructureModified', function(e) {
  // fires after every canvas mutation (draw, delete, undo, paste, AAM)
});

editor.setCallBack('MapHighlight', function(e) {
  // fires when the user clicks a mapped atom (sticky highlight pin)
  // e.atoms = [atomId, atomId, ...] of all atoms sharing that mapNumber
});
```

The full event dispatch path is in `MolEditor._fireCallback`.

---

## Accessibility

BIME targets **WCAG 2.1 Level AA**. Concrete guarantees:

- Every interactive button has a visible focus ring.
- All icon-only buttons carry an `aria-label` matching their tooltip.
- The drawing canvas exposes `<title>` + `<desc>` summarising the molecule.
- Toolbar buttons declare `aria-pressed` for their toggle state.
- All popovers / modals use `role="dialog"` or `role="tooltip"` with
  `aria-modal` and Escape-dismiss.
- The drag-drop zone is `role="button"` and keyboard-activatable
  (Enter / Space).
- Min target size 24 × 24 CSS px throughout.
- `prefers-reduced-motion: reduce` honoured on transitions.

---

## Versioning & compat

BIME follows **Semantic Versioning** strictly. Public API:

- `MolEditor`, its prototype methods, and the constructor signature.
- The `Molecule`, `SmilesParser`, `SmilesWriter`, `RDT`, `Renderer`,
  `ToolbarPrefs` globals and their documented members.
- The bundle's `dist/bime.js` and `dist/bime.min.js` Apache-2.0 banner
  and SRI hash format.

Internal helpers (those starting with an underscore, e.g.
`MolEditor.prototype._buildToolbar`) are subject to change without
notice. They're useful for debugging but should not be relied on by
consumers.

Backwards-compatibility flags:

| Flag | Restores |
|---|---|
| `useBipartitePostPass: false` | v1.2.x greedy component pairing |
| `useLeftoverRescue: false`    | v1.4.1 mapping counts (no rescue) |

Bundle SHA-256 is in `dist/MANIFEST.sha256`; SRI sha384 in `dist/SRI.txt`.
GPG-verifiable tags via `git tag -v v<version>`.

---

## See also

- **`README.md`** — top-level overview, Quick Start, badge row.
- **`CHANGELOG.md`** — full release history with chemistry/algorithm
  rationales for every change.
- **`workbench.html`** — runnable reference layout.
- **`docs.html`** — companion HTML docs page (this guide is the Markdown
  source).
- **`CITATION.cff`** — structured citation for academic use.

If you find a bug or want a new feature, please file an issue at
[github.com/asad/bime-dist/issues](https://github.com/asad/bime-dist/issues).

— BioInception PVT LTD, Cambridge, UK · Apache License 2.0
