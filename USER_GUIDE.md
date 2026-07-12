# BIME User Guide

**BIME** — BioInception Molecular Editor.
A zero-install, zero-dependency browser molecular editor with a headless
command-line companion. Free and Apache-2.0.

This guide takes three audiences from cold start to confident use:

- **End users** — practising chemists, drawing structures and reactions
  for ELN, papers, slides.
- **Students** — biochemistry / organic / pharma courses; using BIME to
  understand structures, reactions, mechanisms, and stereochemistry.
- **Professionals** — cheminformatics / computational-biology engineers
  embedding BIME or scripting it for batch jobs.

Throughout the guide, every command is copy-pasteable. Every CLI example
has been validated against the v3.0.3 sweep.

---

## Table of contents

1. [Getting started](#1-getting-started)
2. [The workbench at a glance](#2-the-workbench-at-a-glance)
3. [Drawing molecules](#3-drawing-molecules)
4. [Building reactions](#4-building-reactions)
5. [Atom–atom mapping (AAM)](#5-atomatom-mapping-aam)
6. [Mechanism arrows](#6-mechanism-arrows)
7. [Highlights — RC / BC / mapped atoms / custom palette](#7-highlights--rc--bc--mapped-atoms--custom-palette)
8. [Exporting](#8-exporting)
9. [The `bime` command-line tool](#9-the-bime-command-line-tool)
10. [For students — a guided tour](#10-for-students--a-guided-tour)
11. [For professionals — scripting recipes](#11-for-professionals--scripting-recipes)
12. [Troubleshooting](#12-troubleshooting)
13. [Glossary](#13-glossary)

---

## 1. Getting started

### Download the public release

```bash
git clone https://github.com/asad/bime-dist.git
cd bime-dist
```

If you do not use `git`, download the public archive instead:

```bash
curl -L https://github.com/asad/bime-dist/archive/refs/heads/main.zip -o bime-dist.zip
unzip bime-dist.zip
cd bime-dist-main
```

### Run BIME in the browser

```bash
python3 -m http.server 8000
# open http://localhost:8000/workbench.html
```

There is no build step. The repo is the deployable artefact. Every page
loads the same `dist/bime.min.js` with an SHA-384 SRI hash baked into
`workbench.html`.

### Build / compile from source

When working from a source checkout, rebuild the browser bundle with:

```bash
node tools/build.js
```

The build writes `dist/bime.js`, `dist/bime.min.js`,
`dist/MANIFEST.sha256`, and `dist/SRI.txt`, then runs the regression suite
against both bundles. The standard build does not require `npm install`.

To verify a downloaded bundle:

```bash
cd dist
shasum -a 256 -c MANIFEST.sha256
```

### Run BIME from the command line

```bash
./bin/bime version
# bime 3.0.3
# RDT.version         3.0.3
# AtomTrace.version   3.0.3
# ToolbarPrefs.version 3.0.3
```

All four lines report the current BIME release and move in lockstep on
every version — `RDT`, `AtomTrace`, and `ToolbarPrefs` are the engine
modules surfaced for provenance. The CLI uses pure Node built-ins (no
`npm install` step). Tested on Node 18+.

> If `./bin/bime` is not executable, run `chmod +x bin/bime`.

### Requirements

| Surface          | Requirement                       |
|------------------|-----------------------------------|
| Browser editor   | Modern browser (Chrome 90+ / Firefox 90+ / Safari 15+) |
| CLI              | Node 18+                          |
| Headless test    | Node 18+                          |

---

## 2. The workbench at a glance

`workbench.html` is a single molecule and reaction editor. It has one drawing
surface — a full toolbar and atom bar over a drawing area — with an **Editor
Output** panel and an **Inspector** docked below it. Both molecules and
reactions live on the *one* editor surface.

- **Drawing area** — draw one molecule, or a full reaction (reactants, a
  reaction arrow, and products), directly on the canvas with the complete
  toolbar (**Draw / Bonds / Rings / Edit / View**) and the atom bar.
- **Editor Output** — read the current structure as **SMILES / SMARTS / MOL /
  SDF / RXN**, copy straight from the textarea, or download via the Export
  buttons.
- **Inspector** — the **Clean & View**, **Reaction Scheme**, and **Templates**
  panels act on the current drawing.
- **Live readout** — a compact strip shows the current **SMILES** (click it to
  copy) and a **warning count** (e.g. *0 warnings* or *⚠ 2 warnings*), updating
  as you draw, and announces updates to screen readers.

```
┌─────────────────────────────────────────────────────────────────┐
│   Toolbar: Draw / Bonds / Rings / Edit / View   ·   atom bar     │
│   ┌──────────────────────────────────────────────────────────┐  │
│   │                                                          │  │
│   │        [ drawing area — a molecule or a reaction ]       │  │
│   │                                                          │  │
│   └──────────────────────────────────────────────────────────┘  │
│   Editor Output (SMILES / SMARTS / MOL / SDF / RXN)  ·  Inspector │
└──────────────────────────────────────────────────────────────────┘
```

---

## 3. Drawing molecules

### From the toolbar

| Tool group | Action |
|---|---|
| **Draw → Bond / Chain** | Click two atoms to draw a bond; chain-tool draws a zig-zag. |
| **Bonds → Single / Double / Triple / Stereo** | Cycle a bond's order or stereo via single-clicks. |
| **Rings → Benzene / Cyclohexane / …** | Click to drop a ring at the cursor. |
| **Edit → Select / Delete / Charge / Isotope** | Modify selected atoms. |
| **View → R/S / E/Z / Atom Labels** | Toggle stereochemistry annotations. |

### From a SMILES paste

Paste any valid SMILES into the **Output → SMILES** box. The editor
parses, lays out, and renders. Round-trip:

```
SMILES in:    Cc1ccccc1
2D rendered:  ⌬ with a methyl
SMILES out:   c1(C)ccccc1   ← canonical reorder
```

### Clean 2D / Best Fit

- **Clean 2D** runs the layout engine (`Layout.layout`) on whatever
  the current coordinates are. Useful after editing.
- **Best Fit** zooms + centres the molecule in the canvas viewport.

### Keyboard shortcuts

The editor follows the hotkey conventions you already know from mainstream
chemistry editors, so drawing stays on the keyboard once a structure is up.
The full list is always one keystroke away — press **?** (or click the
**Shortcuts** button in the Edit toolbar) to open an accessible reference
overlay; **Esc** or a click outside closes it. Every shortcut also has a
matching toolbar button, so the keys are a faster path, never the only one.

| Keys | Action |
|------|--------|
| `1` `2` `3` | Single / double / triple bond tool |
| `C` `N` `O` `S` `P` `F` `I` `H` `B` | Set the element of the atom under the cursor (or the selected atoms); with nothing targeted, arm the atom tool with that element |
| `Shift`+`C` / `Shift`+`B` | Chlorine / bromine |
| `Esc` | Return to the Select tool |
| `Del` / `Backspace` | Delete the selection, or the atom/bond under the cursor |
| `Ctrl`+`Z` | Undo (use `Cmd` on macOS) |
| `Ctrl`+`Shift`+`Z` / `Ctrl`+`Y` | Redo |
| `Ctrl`+`A` | Select all |
| `?` | Open the keyboard-shortcuts overlay |

Shortcuts never fire while you are typing in a text field (a molecule name,
a SMILES box, a label) — so they stay out of the way exactly when you need
the keyboard for something else.

BIME's bonds are single, double, or triple; aromaticity is a property of the
ring (toggle its display with the **Ar** button), so there is no separate
aromatic-bond or any-bond key.

### Build with the keyboard — no mouse needed *(v2.4.11)*

The draw canvas is fully operable from the keyboard, for people who do not use a
mouse and for screen-reader users. **Tab** to the drawing canvas (it is in the
focus order and announces what is on it), then press **Enter** to start *keyboard
build mode*:

- **Arrow keys** grow a bonded atom from the current atom in that direction — up, down, left, or right. Arrow back toward an atom already there and BIME bonds to it instead, closing a ring.
- **C, N, O, …** set the element — the current atom is retyped and the next grown atom uses it (**Shift+C** / **Shift+B** for chlorine / bromine).
- **1, 2, 3** set the bond order (single, double, triple) for the next growth.
- **Delete** removes the current atom; a neighbouring atom becomes current.
- **Escape** leaves keyboard build mode.

The current atom is shown with the selection highlight, every step is announced
aloud (the new atom, the bond, and the running atom and bond count), and every
change is undoable with **Ctrl+Z**. Tidy the geometry afterwards with **Clean
2D**. Build mode only takes over the keyboard while it is active and the canvas
has focus, so the rest of the workbench is unaffected.

---

## 4. Building reactions

A reaction in BIME is one molecule with a *reaction arrow* between
groups of atoms. Two ways to create one:

### Paste a reaction SMILES

```
CCO.CC(=O)O>>CCOC(=O)C.O
```

Reactants are on the left of `>>`, products on the right. BIME parses
the arrow and renders R + arrow + P.

### Draw it interactively

1. Draw the reactant(s).
2. Click the toolbar `Reaction → Arrow` tool.
3. Click the canvas where the arrow should sit.
4. Draw the product(s) to the right of the arrow.

The Inspector → **Reaction Scheme** panel exposes:

- **Arrow type** — forward (→), reverse (←), reversible (⇌), resonance (↔). Each
  is drawn with a distinct *shape* so the four read apart without colour: a single
  head at the end / start, two offset half-harpoons for reversible, and full heads
  at both ends on a **dashed** shaft for resonance.
- **Conditions** — temperature, solvent, catalyst (free-text).
- **Yield** — % yield as a number.
- **Step note** — short caption.

These are reflected in the SVG/PDF export.

---

## 5. Atom–atom mapping (AAM)

### What AAM does

AAM assigns each atom in the reactant side an integer "map number"
that matches one atom on the product side. Same map number = same
atom across the reaction. BIME's atom-mapping engine handles this in
pure JavaScript.

### Running AAM in the editor

1. Build (or paste) a reaction.
2. Click `Reaction → Auto-map`.
3. BIME computes the mapping and:
   - Writes `:n` labels next to mapped atoms.
   - Populates `result.bondChanges` (cleaved / formed / order /
     stereo changes).
   - Surfaces a **Quality** grade (A/B/C/F) in the Inspector.

### Running AAM from the CLI

```bash
$ bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O'
Status:       mapped
Quality:      F
Score:        2
Mapped atoms: 6
Bond changes: 5
Reaction:     substitution
Confidence:   0.659
Decisiveness: 0.000
Timed out:    false

$ bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O' --format json
{
  "status": "mapped",
  "score": 2,
  "quality": "F",
  "mapping": { "0": 0, "1": 1, … },
  "bondChanges": [ … ],
  "diagnostic": { … }
}
```

### Quality grades

| Grade | Meaning |
|---|---|
| **A** | Confidence ≥ 0.95, decisiveness ≥ 0.5, aromatic rings + stereo centres preserved |
| **B** | Confidence ≥ 0.85 |
| **C** | Confidence ≥ 0.65 |
| **F** | Below C threshold — review manually |

### Chemistry-aware tie-breaking

Default ON since v2.0.57. Penalises mappings that break aromatic rings
or drop stereo centres. Golden corpus (138 reactions) still passes
138/138 under chemistry-aware mode. Disable with:

```bash
bime aam '<reaction>' --no-chemistry-aware
```

---

## 6. Mechanism arrows

Curly arrows show electron movement in a mechanism step. Two flavours:

- **Pair (double-barb)** — a lone pair or bond electrons (2 e⁻).
- **Radical (single-barb)** — a single electron (1 e⁻).

To add one:

1. Pick `Curly: pair` or `Curly: radical` in the Inspector.
2. Click the source atom (or bond).
3. Click the target atom (or bond).

Arrows persist with the molecule and round-trip through JSON
serialisation. They render in the SVG / PDF export.

---

## 7. Highlights — RC / BC / mapped atoms / custom palette

After running AAM (Auto-map), three overlays are available to surface
the diagnostic data:

| Overlay | What it draws | Toggle (JS API) |
|---|---|---|
| **Reaction Centre** | Dashed amber halo on every atom in `bondChanges` | `editor.toggleReactionCenter()` |
| **Bond Changes** | Per-bond stroke colour by kind (🟥 cleaved · 🟩 formed · 🟧 order · 🟪 stereo) | `editor.toggleBondChanges()` |
| **Mapped atoms** | Soft teal halo behind every atom with `:n` map number | `editor.toggleMappedAtoms()` |

### Custom palette overrides

```js
editor.setColorOverride('C', '#222222');   // every carbon halo dark
editor.setColorOverride('1', '#ff00ff');   // map :1 atoms magenta
editor.clearColorOverrides();              // back to default
```

Keys may be element symbols (`C`, `N`, `O`, …) or map-number strings
(`'1'`, `'2'`, …).

---

## 8. Exporting

### Every export carries a BIME provenance stamp (v2.0.65+)

Default ON across the workbench AND the CLI. The stamp is:

- **Images (SVG / PNG / PDF)** — a `<metadata id="bime-stamp">` block
  carrying the BIME version, ISO-8601 timestamp, SHA-384 fingerprint,
  copyright, Apache-2.0 licence, and repo URL. Plus a bottom-right
  `<g class="bime-watermark">` with the BIME wordmark + 6-vertex
  glyph at 0.55 opacity (visible, but never competes with the
  structure).
- **Text (MOL / SDF / RXN / SMILES)** — a `# BIME-STAMP-BEGIN` …
  `# BIME-STAMP-END` comment block carrying the same record. MOL/SDF/
  RXN get a leading stamp; SMILES gets a trailing stamp so the
  first non-comment line is still parseable by any SMILES consumer.

#### Opt out

If a downstream tool refuses comment-prefixed input, opt out:

```bash
bime export 'Cc1ccccc1' --format mol --no-stamp
```

#### "Encrypted for public domain" — what that actually means

The stamp is **tamper-evident signing**, not asymmetric encryption.
The two are very different and the BIME design deliberately chose
the former:

- Apache-2.0 artefacts must remain readable by anyone. Encrypting
  the bytes would break the licence.
- What you want is attribution + integrity: any consumer can
  recompute the SHA-384 over the un-stamped content and verify
  the artefact came from BIME unmodified.
- Asymmetric signatures (Ed25519 / RSA) would add cryptographic
  identity-proof but need a private key held only by BioInception.
  They live in a future release alongside the existing stamp, not
  in place of it.

The TL;DR: the stamp **witnesses**, it does not **lock**.
Confidentiality is not promised; **attribution + integrity** are.

### From the browser

The Output panel has SMILES, SMARTS, MOL, SDF, RXN tabs — copy
straight from the textarea, or download via the Export buttons:

- **SVG** — vector, full quality, retains stereo wedges + atom labels.
- **PNG** — rasterised SVG, configurable DPI.
- **PDF** — vector PDF via in-browser PDF generation.
- **MOL / SDF / RXN** — MDL formats.
- **SMILES / SMARTS** — text formats.

### Publication-quality figures *(v2.4.11)*

For journal figures and slides, the **Export Presets** panel (in the **Import / Export** drawer) carries three one-click presets:

- **Publication SVG** — a publication-styled vector figure on a white background, following journal-standard drawing conventions: Helvetica / Arial atom labels, Kekulé double bonds (explicit, no aromatic circles), print-safe element colours, heavier bonds for legibility, and no watermark.
- **Transparent SVG** — the same publication style with a transparent background, so the figure drops cleanly onto any page, poster, or dark slide (the full-canvas background is removed and the atom-label knock-out boxes go clear, so nothing but the structure paints).
- **Transparent PNG** — the publication style rasterised at 4× / ≥1024 px on a transparent background, for tools that need a bitmap.

The same presets are on the workbench command palette — search "publication" or "transparent". All three are opt-in: the plain SVG and PNG download buttons are unchanged.

### From the CLI

```bash
# Export benzene as SVG (default 320×240, transparent background)
bime export 'c1ccccc1' --format svg --out benzene.svg

# Custom canvas
bime export 'Cc1ccccc1' --format svg --width 400 --height 300

# MOL file
bime export 'Cc1ccccc1' --format mol > toluene.mol

# SDF (single record terminated by $$$$)
bime export 'Cc1ccccc1' --format sdf > toluene.sdf

# RXN (requires reaction SMILES)
bime export 'CCO.CC(=O)O>>CCOC(=O)C.O' --format rxn > esterification.rxn

# Canonical SMILES round-trip
bime export 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O' --format smiles
```

PNG / PDF from the CLI need a Canvas polyfill (Node has no DOM); the
CLI politely refuses with a clear message. SVG covers most use cases.

### Publication-quality reaction-mapping figures *(v2.4.13)*

`bime aam <reaction> --format svg` renders a complete, **publication-quality
mapped-reaction figure** — the static, batchable equivalent of the workbench's
on-screen mapping view. It runs the atom-atom mapper, lays the reaction out
(reactants → arrow → products, with `+` separators), and draws four overlays:

- **Colour mol** — soft per-atom halos tint each mapped sub-fragment (the CPK
  element colours stay intact underneath).
- **Atom-atom map numbers** — a reactant atom and its product partner carry the
  same number.
- **Bond changes** — changed bonds are recoloured: **red** cleaved, **green**
  formed, **amber** order change. (A hydrogen-count change has no bond to
  colour, so it marks a reaction-centre atom instead — shown with
  `--reaction-center`.)
- **Trace** — the sub-fragment tint is the same on both sides of the arrow, so
  each atom's fate reads left-to-right.
- **Compound labels** *(v2.4.14)* — each structure is captioned with its
  compound name when recognised: every component is matched against BIME's
  built-in metabolite repertoire by structure (robust to how the SMILES is
  written, and an epimer is never mislabelled as its partner — matching is by
  relative configuration, so the natural biological enantiomer resolves), so a glucose +
  ATP reaction draws "D-Glucose", "ATP", etc. under the structures. Unknown
  components are left uncaptioned. `--ids` captions with the KEGG id instead of
  the name; `--no-labels` turns captions off.

The figure is auto-sized to the reaction (no dead whitespace), uses the
print-quality Helvetica/Arial preset on a white background, and carries the
usual provenance stamp. (No confidence caption is drawn.)

```bash
# A mapped esterification, straight to a figure
bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O' --format svg --out esterification.svg

# Transparent background for a slide, with reaction-centre rings
bime aam 'CCO>>CC=O' --format svg --transparent --reaction-center --out oxidation.svg

# The on-screen look instead of the print preset; pin a size; drop the stamp
bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O' --format svg --screen --width 1000 --height 300 --no-stamp
```

Flags: `--screen` (UI look instead of print), `--transparent`, `--reaction-center`,
`--no-map-numbers`, `--ids` (KEGG ids instead of names), `--no-labels`,
`--width N` / `--height N` (default: auto-fit), `--out file`, `--no-stamp`.

---

## 9. The `bime` command-line tool

### Synopsis

```
bime <command> [args] [flags]
```

### Commands

| Command                | What it does                                                        |
|------------------------|---------------------------------------------------------------------|
| `version`              | BIME + module versions (one per line)                               |
| `parse <smi>`          | SMILES diagnostic (atoms, bonds, elements, charges, stereo)         |
| `clean <smi>`          | `Layout.layout` + canonical SMILES                                  |
| `aam <rxn>`            | `RDT.mapReaction` on a reaction SMILES `A>>B`                       |
| `export <smi>`         | Write SVG / MOL / SDF / RXN / SMILES                                |
| `help [cmd]`           | Top-level help, or per-command flag reference                       |

### Universal flags

- `--in <file>` — read input from a file. `--in -` reads from stdin.
- `--out <file>` — write output to a file (default: stdout).
- `--format <text|json>` — output format for `parse` and `aam`.

### Examples per command

```bash
# version
bime version

# parse
bime parse 'c1ccccc1'
bime parse 'c1ccccc1' --format json
echo 'CCO' | bime parse

# clean
bime clean 'Cc1ccccc1'

# aam
bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O'
bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O' --format json
bime aam 'CCO>>CC=O' --no-chemistry-aware --timeout-ms 1000

# export
bime export 'Cc1ccccc1' --format svg --width 400 --height 300 --out toluene.svg
bime export 'Cc1ccccc1' --format mol > toluene.mol
bime export 'Cc1ccccc1' --format sdf > toluene.sdf
bime export 'CCO>>CC=O' --format rxn > oxidation.rxn

# help
bime help
bime help aam
```

### Exit codes

| Code | Meaning |
|---|---|
| `0` | Success |
| `1` | User error (unknown command, malformed input, missing required flag) or engine failure (parse / AAM error) |

Useful in shell scripts:

```bash
if bime aam "$REACTION" --format json > out.json; then
    echo "AAM OK"
else
    echo "AAM failed: see stderr"
fi
```

---

## 10. For students — a guided tour

**Audience:** undergrad biochemistry / organic / pharma students.

### Day 1 — Draw morphine, get a SMILES

1. Open the workbench. The molecule editor is ready to draw.
2. Paste into the **Output → SMILES** box and click **Load**:
   ```
   CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5
   ```
   (or draw it atom by atom with the toolbar). Try the toolbar **Best Fit**
   to centre it.
3. Read its SMILES/MOL back in **Editor Output** below the drawing area.

### Day 2 — Run AAM on a cross-aldol condensation

```
CC(=O)C.CC(=O)CC(=O)C>>CC(=O)C(=CC(=O)C)C.O
```

(Classical Knoevenagel pairs an active-methylene compound with a
carbonyl; this is a closer fit to a Claisen–Schmidt cross-aldol —
the AAM exercise is identical either way.)

1. Paste into the SMILES box.
2. Click **Reaction → Auto-map**.
3. Observe the map numbers `:n` appear next to atoms.
4. Look at the Inspector → Quality grade. What was it (A/B/C/F)?
5. Look at `result.bondChanges` in the Output → JSON tab.
6. Toggle **Show Reaction Center** to see which atoms moved.

### Day 3 — Mechanism arrows on the SN2 inversion

1. Draw `[C@@H](Br)(C)Cl` in the molecule editor.
2. Click `Curly: pair`.
3. Click the lone pair on the leaving Br, then click the central C.
4. Add a second curly: the C–Cl bond → Br.
5. Save the SVG via Export → SVG.

---

## 11. For professionals — scripting recipes

**Audience:** cheminformatics engineers, computational biologists.

### Recipe 1: batch-AAM a CSV of reactions

```bash
#!/usr/bin/env bash
# in:  reactions.csv  one reaction SMILES per line
# out: aam.jsonl      one JSON record per line
while IFS= read -r rxn; do
    bime aam "$rxn" --format json
done < reactions.csv > aam.jsonl
```

### Recipe 2: embed BIME as a structure widget in your app

```html
<!DOCTYPE html>
<html>
<head>
    <!-- Pin your deployed version's SRI hash from dist/SRI.txt. -->
    <script src="bime.min.js"
            integrity="sha384-…"
            crossorigin="anonymous"></script>
</head>
<body>
    <div id="bime-editor" style="width: 800px; height: 500px"></div>
    <script>
        var editor = new MolEditor({ container: document.getElementById('bime-editor') });
        editor.on('changed', function() {
            console.log('SMILES:', editor.smiles());
        });
    </script>
</body>
</html>
```

The bundle exposes `MolEditor`, `Molecule`, `SmilesParser`, `SmilesWriter`,
`RDT`, `AtomTrace`, `MetaboliteLibrary`, `Layout`, `ImageExport`,
`MolfileWriter`. See [tests/test_v2_0_*](tests/) for a battery of
short usage examples.

### Recipe 3: golden-corpus regression test

```bash
node /PATH/TO/bime/tools/run-tests.js | tail -5
# 1842 passed, 0 failed, total 9XYZ ms
```

Drop a new SMILES into `tests/data/golden_reactions.json` to extend
the corpus. The `test_v2_0_x_golden_reactions.js` runner picks it up
automatically.

---

## 12. Troubleshooting

### "Parse failed" — my SMILES isn't recognised

- Check stereo: BIME accepts `[C@H]`, `[C@@H]`, `/`, `\`. Brackets are
  required for non-organic-subset atoms.
- Lowercase aromatic (`c`, `n`, `o`, `s`) atoms must be inside a ring.
- Reaction SMILES uses `>>`, not `>`.

### AAM status = `unbalanced`

The reaction's heavy-atom counts don't match across `>>`. BIME refuses
to map unbalanced reactions by default. Either:
- Fix the reaction stoichiometry, or
- Add the explicit cofactor on the missing side
  (`CCO.[OH-]>>CC=O.[H2O]`), or
- Disable the check with `--require-balance false` *(CLI flag coming
  in v2.0.65)*.

### Workbench: how do I draw a brand-new structure?

On first load the editor is empty and ready to draw. Build a structure with
the **Draw → Bond / Chain** tools and the atom bar, drop a ring from the
**Rings** group, or paste a SMILES into the **Output → SMILES** box and click
**Load** to have BIME parse, lay out, and render it. Read the result back as
**SMILES / SMARTS / MOL / SDF / RXN** in Editor Output — with a live SMILES and
warning count in the readout strip *(v2.3.8)* — and use **Best Fit** to centre
it or **Clean 2D** to tidy the layout.

### R/S descriptors on symmetric rings *(v2.4.12)*

> **R/S and symmetric rings (v2.4.12).** BIME never guesses a stereodescriptor. A stereocentre whose CIP priorities cannot be resolved unambiguously is left unlabelled rather than shown a possibly-wrong R/S. This is most visible on the symmetric ring cyclitols (the inositol family): each ring carbon's two ring paths are constitutionally identical, so the descriptor cannot be settled without a stereo-aware ring analysis. Acyclic centres — including pseudo-asymmetric ones (lowercase r/s) and meso sugar alcohols (which are labelled with the correct mix of R and S) — resolve normally. The symmetric-ring abstention is a deliberate, documented never-wrong choice: an honest blank, not a wrong letter, pending a stereo-aware ring resolver.

### Accessibility: colour is never the only signal

Everything the editor conveys with colour also carries a **non-colour** cue, so
the workbench reads in greyscale and to colour-blind users:

- **Reaction-arrow type** — the four types differ by *shape* (see §4), not colour.

The **drawing canvas** is reachable and described for screen-reader and
keyboard users: it is in the focus order (**Tab** to it) and announces a live
description of what is drawn — "Molecule drawing canvas, empty canvas" when
nothing is there, otherwise the molecule name (if set) plus an atom/bond tally
(and "reaction scheme" when a reaction arrow is present), updating as you draw.
A focus ring shows when keyboard focus lands on it. (Drawing atoms *with the
keyboard* is planned for a later release; today the canvas is reachable and
self-describing, and editing is via pointer.)

### CLI: `bime export 'X' --format png` is refused

PNG and PDF need a Canvas (browser DOM). Use SVG from the CLI; PNG/PDF
from the browser Export buttons; or pipe SVG through an external
converter:
```bash
bime export 'Cc1ccccc1' --format svg | rsvg-convert -o toluene.png
```

### "Module not found" running the CLI

The CLI must be run from inside the cloned repo (it uses relative
paths to find `tests/shim.js` and `editor/*.js`). Don't `cp bin/bime`
elsewhere. To install globally:
```bash
ln -s "$(pwd)/bin/bime" /usr/local/bin/bime
```

---

## 13. Glossary

| Term | Meaning |
|---|---|
| **AAM** | Atom–atom mapping. Each reactant atom is matched to one product atom (or remains unmapped). |
| **Map number** | Integer label `:1`, `:2`, … attached to an atom indicating its mapping partner. |
| **Bond change** | `{kind, atom1, atom2}` record listing every cleaved / formed / order-change / stereo-change bond in a reaction. |
| **Reaction centre (RC)** | The union of atoms participating in any bond change — the chemically active subset of the structure. |
| **Component pair** | A maximal connected scaffold piece preserved across a reaction — used for sub-fragment colouring. |
| **Curly arrow** | IUPAC mechanism arrow showing electron movement. Single-barb for a single electron, double-barb for a pair. |
| **SRI** | Subresource Integrity. Every BIME bundle is SHA-384-pinned in `workbench.html` and `dist/SRI.txt`. |

---

## Where next?

- **Public release + issues:** <https://github.com/asad/bime-dist>
- **Citation:** see `CITATION.cff` in the repo root — Rahman, S. A., *BIME: BioInception Molecular Editor* (2026).
- **License:** Apache-2.0 — see `LICENSE.txt`.
- **Bug reports:** GitHub Issues. Include a minimal reproducer (SMILES + CLI command or browser steps).

This guide tracks the v2.0.64 feature set. Subsequent releases extend
without removing — older browsers and scripts keep working.
