# BIME for Students

Welcome! This guide is for high-school and undergraduate students
learning chemistry. It assumes no prior programming experience.

By the end of this guide you will know how to:

- Draw molecules in BIME.
- Read and write SMILES — the chemical text language.
- Search for functional groups using SMARTS.
- Compare two molecules and find their shared scaffold (MCS).
- See how atoms move through a reaction (atom-atom mapping).
- Export your drawings for assignments and lab reports.

---

## Getting started in 60 seconds

1. Open the **[BIME workbench](https://asad.github.io/bime-dist/workbench.html)**
   in any browser. **No login. No install.**
2. You'll see a blank drawing canvas in the middle and tool icons on
   the left.
3. Click an atom-bar element (C, N, O, etc.), then click the canvas to
   place atoms. Drag between two atoms to draw a bond.
4. Watch the SMILES update live in the panel below the canvas.

That's it. You're drawing chemistry.

---

## Drawing tools — the cheat sheet

| Icon | Tool | What it does |
|---|---|---|
| **C** | Atom | Click in the canvas to place a carbon. Switch to N, O, etc. via the atom bar. |
| **─** | Single bond | Drag between two atoms to draw a single bond. |
| **=** | Double bond | Drag between two atoms; or click a single bond to upgrade it. |
| **≡** | Triple bond | Drag between two atoms. |
| **○** | Ring | Click a position to drop a benzene/cyclohexane/etc. ring. Switch ring size in the toolbar. |
| **/\\/** | Chain | Click and drag to lay down a zigzag carbon chain. |
| **▲** | Wedge | Drag to add a stereo wedge (out of page). |
| **▼** | Dash | Drag to add a stereo dash (into page). |
| **⨯** | Delete | Click to remove an atom or bond. |
| **↶ ↷** | Undo / Redo | `Ctrl+Z` / `Ctrl+Y` (Cmd-Z / Cmd-Y on Mac). |

Tip: after a few minutes you'll find the **keyboard shortcuts** faster
than the mouse. Hover over any tool to see its shortcut.

---

## SMILES quick reference

SMILES is a one-line text representation of a molecule. It is the
"DNA sequence" of chemistry — compact, copy-pasteable, unambiguous.

### Read these three SMILES strings

```
CCO              # ethanol
c1ccccc1         # benzene
CC(=O)Oc1ccccc1C(=O)O   # aspirin
```

### The basics in five rules

1. **Letters = atoms.** `C` is carbon, `N` is nitrogen, `O` is oxygen,
   `H` is hydrogen, `Cl` is chlorine, `Br` is bromine.
2. **Adjacency = bonds.** `CC` = two carbons single-bonded.
   Like reading sequential atoms.
3. **`=` and `#` = double and triple bonds.** `C=C` is ethene,
   `C#C` is ethyne.
4. **Ring closure with a digit.** `C1CCCCC1` = cyclohexane (the `1`s
   pair up to close the ring).
5. **Lowercase = aromatic.** `c1ccccc1` = benzene with proper aromatic
   delocalisation.

### Branching with parentheses

```
CC(C)C          # isobutane (the middle C has a methyl branch)
CC(=O)C         # acetone (the central C has a =O branch)
CC(=O)O         # acetic acid
CC(=O)OC        # methyl acetate
```

### Special atoms

```
[NH4+]          # ammonium ion (charged)
[O-]            # oxide
[Na+]           # sodium cation
[U]             # uranium (anything not C/N/O/F/S/P/Cl/Br/I/B uses brackets)
[13C]           # carbon-13 isotope
```

### Stereochemistry

```
C[C@@H](N)C(=O)O    # (S)-alanine — anticlockwise around the chiral C
C[C@H](N)C(=O)O     # (R)-alanine — clockwise
C/C=C/C             # (E)-2-butene — trans
C/C=C\C             # (Z)-2-butene — cis
```

The `@` is anticlockwise, `@@` is clockwise (looking from the first
neighbour). The `/` and `\` are bond directions for E/Z double bonds.

### Reactions

```
CC(=O)O.OCC>>CC(=O)OCC.O      # esterification
                              # reactants on the left of >>, products on right
                              # . separates components
```

---

## Common molecules — try these in the workbench

| Molecule | SMILES |
|---|---|
| Water | `O` |
| Ammonia | `N` |
| Methane | `C` |
| Carbon dioxide | `O=C=O` |
| Methanol | `CO` |
| Ethanol | `CCO` |
| Acetic acid | `CC(=O)O` |
| Acetone | `CC(=O)C` |
| Benzene | `c1ccccc1` |
| Phenol | `c1ccc(O)cc1` |
| Aniline | `c1ccc(N)cc1` |
| Toluene | `Cc1ccccc1` |
| Glucose (open) | `OCC(O)C(O)C(O)C(O)C=O` |
| Caffeine | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` |
| Aspirin | `CC(=O)Oc1ccccc1C(=O)O` |
| Paracetamol | `CC(=O)Nc1ccc(O)cc1` |
| Ibuprofen | `CC(C)Cc1ccc(C(C)C(=O)O)cc1` |

Type any of these in the SMILES box at the bottom of the workbench
and press **Load**. The structure appears instantly.

---

## SMARTS — searching for functional groups

SMARTS is like SMILES, but for **patterns**. Use it to ask "where in
this molecule is there a carboxylic acid?"

| SMARTS | Matches |
|---|---|
| `[OH]` | hydroxyl groups |
| `[NX3;H2]` | primary amines |
| `[CX3]=O` | carbonyl groups |
| `[CX3](=O)O` | carboxylic acids |
| `[NX3][CX3](=O)` | amides |
| `[CX3](=O)O[CX4]` | esters |
| `c1ccccc1` | benzene rings |
| `[F,Cl,Br,I]` | halides |
| `[#7]` | any nitrogen |
| `[r6]` | any atom in a 6-ring |

### Try it

1. Load aspirin (`CC(=O)Oc1ccccc1C(=O)O`) in the workbench.
2. Switch to the **Examples → SUB** tab.
3. Type `[CX3](=O)O` in the SMARTS box.
4. Press **Search**. Watch the carboxylic acid light up.
5. Now try `[CX3](=O)O[CX4]` — the ester also matches.

---

## MCS — comparing two molecules

MCS = **Maximum Common Substructure** = the biggest piece two molecules
share.

### Try it

1. Open **Examples → MCS** in the workbench.
2. The page loads aspirin and salicylic acid side by side.
3. The shared scaffold (the salicylate part) is highlighted in colour.
4. The bits that differ — the acetyl group on aspirin — are not
   highlighted.

This is exactly how medicinal chemists ask "what's the same and what's
different between these two drugs?"

Try it on:

- Caffeine vs theobromine — same purine core, different methylation.
- Ibuprofen vs naproxen — same arylpropionate scaffold, different ring.
- Cocaine vs procaine — both local anaesthetics, very different
  scaffolds (BIME's MCS confirms it).

---

## Atom-atom mapping — watching atoms through a reaction

When a reaction happens, every atom from the reactants ends up
**somewhere** in the products. Atom-atom mapping (AAM) traces those
trajectories.

### Try it

1. Open **Examples → AAM**.
2. The page loads `CCO>>CC=O` (ethanol → acetaldehyde).
3. The mapped atoms are connected by coloured halos.
4. Below the structure, a **Bond Changes** panel shows:
   - Which bonds **formed** (none in this case).
   - Which bonds **cleaved** (the C–H and O–H bonds removed during
     oxidation).
   - Which bond **changed order** (the C–O single → double).
   - Which atoms **changed hydrogen count** (the α-C lost 1 H, the
     O lost 1 H).

Try it on esterification, Diels–Alder, SN2 — every reaction tells a
story.

---

## Saving and sharing your work

### Save a molecule

- **File → Save SMILES** — gives a one-line text file.
- **File → Save MOL** — gives a .mol file (V2000) compatible with
  ChemDraw, MarvinSketch, RDKit, OpenBabel.
- **File → Save MOL V3000** — for very large molecules.

### Export an image

- **File → Export SVG** — vector format, scales to any size.
- **File → Export PNG** — pick 1×, 2×, or 4× resolution.
- **File → Copy to Clipboard** — paste straight into Word, Pages, or
  LibreOffice.
- **File → Export print SVG** — CMYK-safe colours, thicker bonds for
  paper.

### Share with a classmate

The SMILES from the SMILES panel **is** your work. Copy it, paste it
into a chat message or email, and your classmate can load the same
structure by pasting it into their workbench.

---

## Common questions

**Q. My SMILES has a red warning. What does it mean?**

A. The atom has too many or too few bonds for its standard valence.
   Click the atom — BIME shows what's wrong and how to fix it. You
   can override with explicit `[N+]`, `[O-]`, or radical notation if
   needed.

**Q. The 2D layout looks tangled.**

A. Click **Tools → Clean Layout** (or press `L`). BIME re-runs the
   layout engine and produces a cleaner depiction.

**Q. I drew a molecule but the SMILES looks weird / not what I expected.**

A. BIME outputs the **canonical** SMILES — the one and only correct
   form for a given molecule. Your hand-written SMILES might describe
   the same molecule, but in a different order. They are equivalent.

**Q. Where do I find more molecules?**

A. The workbench has a built-in browser with **1181 molecules**
   organised by category (drugs, natural products, amino acids,
   sugars, vitamins, etc.). Click **Browse Molecules** to explore.

**Q. Does BIME upload my work anywhere?**

A. **No.** BIME runs entirely in your browser. Open DevTools (F12)
   → Network tab and you'll see zero outgoing requests. Your
   chemistry stays on your computer.

**Q. Can I use BIME on my phone / tablet?**

A. Yes. The interface is fully responsive and supports touch,
   pinch-to-zoom, and stylus input.

**Q. Can I use BIME for my exam answers / lab report?**

A. Absolutely. Export an SVG or PNG and drop it into your document.
   See the **Print SVG** option for paper-friendly output.

**Q. Do I need to credit BIME if I use it?**

A. For internal class use, no. For published work (theses, papers),
   citing BIME is appreciated:

   > Rahman, S. A. *BIME: BioInception Molecular Editor.* BioInception
   > PVT LTD, Cambridge, UK (2026). https://github.com/asad/bime-dist

   See [CITATION.cff](../CITATION.cff) and the README for BibTeX.

---

## Where to go next

- **[USAGE.md](USAGE.md)** — every BIME feature, in detail.
- **[The Examples gallery](https://asad.github.io/bime-dist/examples.html)**
  — interactive, click-through demos of every feature.
- **[The Documentation page](https://asad.github.io/bime-dist/docs.html)**
  — keyboard shortcuts, tool reference, and the JS API.
- **[GitHub Issues](https://github.com/asad/bime-dist/issues)** — ask a
  question, report a bug, or suggest a feature.

Welcome to BIME. Now go draw something chemical.
