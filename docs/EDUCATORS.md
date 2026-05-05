# BIME for Educators

A guide for chemistry teachers, professors, and instructors.

BIME is built for the classroom: zero-install, runs on any device with a
browser, and supports the full spectrum from secondary-school SMILES
notation to graduate-level reaction-mapping pedagogy.

This guide shows how to integrate BIME into your teaching.

---

## Table of contents

1. [Why BIME for the classroom](#1-why-bime-for-the-classroom)
2. [Quick wins — five things you can do today](#2-quick-wins--five-things-you-can-do-today)
3. [Lesson plans by topic](#3-lesson-plans-by-topic)
4. [Pre-built example reactions](#4-pre-built-example-reactions)
5. [Embedding BIME in your LMS](#5-embedding-bime-in-your-lms)
6. [Printable worksheets](#6-printable-worksheets)
7. [Assessment ideas](#7-assessment-ideas)
8. [Accessibility](#8-accessibility)
9. [Common student questions](#9-common-student-questions)
10. [Citing BIME in teaching materials](#10-citing-bime-in-teaching-materials)

---

## 1. Why BIME for the classroom

| Need | BIME's answer |
|---|---|
| **Zero IT setup** | Open `workbench.html` in any browser. No login, no install, no plugin, no flash, no Java. |
| **Works on student laptops** | Pure JavaScript + SVG. Chromebooks, Macs, Windows, Linux, iPads, Android tablets all work. |
| **Works offline** | Drop the folder onto a USB stick or shared drive — runs without internet. |
| **No per-seat licensing** | Apache 2.0. Use it in any class size, public or private. |
| **No data collection** | Zero analytics, zero telemetry, zero tracking. Student SMILES never leaves the browser. |
| **Print-friendly** | SVG export at any DPI. CMYK-safe colours and thicker bonds in print mode. |
| **Honest chemistry** | 1181 PubChem-canonical molecules; aromaticity audit reports 0 issues across the entire database. |

---

## 2. Quick wins — five things you can do today

### 1. Show students how SMILES encodes structure

Open the workbench, click **Examples**, and load **Caffeine**. The page
shows the SMILES `Cn1c(=O)c2c(ncn2C)n(C)c1=O` and the rendered structure
side-by-side. Students can:

- Edit the SMILES live and see the structure update.
- Click an atom and see its valence, charge, hybridisation.
- Toggle aromatic / Kekulé display.

### 2. Demonstrate substructure search

In **Examples → SUB**, students can search for a SMARTS pattern
(`[NX3;H2]` for primary amines, `[CX3](=O)O` for carboxylic acids,
`c1ccccc1` for benzene rings) across a target molecule.

### 3. Draw a reaction and watch atoms map

In **Examples → AAM**, the workbench shows mapped atoms across a
reaction. Esterification, Diels–Alder, SN2, hydration — all worked
out interactively. Excellent for teaching mechanism-tracking.

### 4. Compare two molecules' shared scaffold (MCS)

**Examples → MCS** highlights the maximum common substructure between
caffeine and theobromine, between aspirin and salicylic acid, between
ibuprofen and naproxen. A vivid demonstration of "same scaffold,
different substituent → different drug".

### 5. Export a clean SVG / PNG for handouts

**File → Export SVG** or **Export PNG (4×)**. The output is
publication-quality, prints sharp on A4 / Letter, and is CMYK-safe in
print mode.

---

## 3. Lesson plans by topic

### High school / first-year undergraduate

#### Lesson: SMILES — a chemical text language

- **Objective**: students can read and write simple SMILES strings.
- **Time**: 50 minutes.
- **Activity**:
  1. Show methane `C`, ethane `CC`, propane `CCC` — note the implicit
     hydrogens.
  2. Show benzene `c1ccccc1` — explain ring closure and aromatic
     lowercase.
  3. Show ethanol `CCO` and acetic acid `CC(=O)O` — branches and
     double bonds.
  4. Have students write SMILES for: water, ammonia, methanol,
     formaldehyde, acetone, glucose (open chain), benzaldehyde,
     phenol, aniline.
  5. Each student loads their SMILES in BIME and verifies the
     structure.
- **Assessment**: write down the SMILES for ten student-chosen molecules
  found in their kitchen / bathroom (caffeine, citric acid, vitamin C,
  ibuprofen, aspirin, urea, paracetamol, baking soda's bicarbonate,
  salt, sugar). Verify in BIME.

#### Lesson: Functional groups and SMARTS

- **Objective**: students can recognise a functional group and write
  the SMARTS for it.
- **Time**: 50 minutes.
- **Activity**:
  1. Walk through SMARTS for: hydroxyl `[OH]`, primary amine
     `[NX3;H2]`, carbonyl `[CX3]=O`, carboxylic acid `[CX3](=O)O`,
     amide `[NX3][CX3](=O)`, ester `[CX3](=O)O[CX4]`, halide `[F,Cl,Br,I]`.
  2. Load aspirin in the workbench. Run each SMARTS via **Examples
     → SUB** and observe which atoms light up.
  3. Have students predict, then verify, the SMARTS counts for
     paracetamol, ibuprofen, naproxen, caffeine, glucose.
- **Assessment**: given an unknown SMILES, identify all functional
  groups present and write SMARTS for each.

### Second-year undergraduate / first-year medicinal chemistry

#### Lesson: Stereochemistry and CIP

- **Objective**: students understand R/S and E/Z assignment by the
  Cahn–Ingold–Prelog rules.
- **Time**: 80 minutes.
- **Activity**:
  1. Load (S)-alanine and (R)-alanine: `C[C@@H](N)C(=O)O` and
     `C[C@H](N)C(=O)O`. Click each chiral centre — BIME shows the
     CIP assignment.
  2. Load (E)-2-butene and (Z)-2-butene: `C/C=C/C` and `C/C=C\C`.
     Same exercise.
  3. Load thalidomide: `O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1` and
     its R-enantiomer with `[C@@H]`. Discuss the historical impact.
  4. Have students draw their own chiral molecule, predict the CIP
     assignment, and verify in BIME.
- **Assessment**: given five SMILES with `@` and `@@` markers, predict
  R/S for each chiral centre. Verify.

#### Lesson: Reaction mapping (atom-atom mapping)

- **Objective**: students can trace which atoms map across a reaction.
- **Time**: 80 minutes.
- **Activity**:
  1. Load esterification `CC(=O)O.OCC>>CC(=O)OCC.O` in **Examples →
     AAM**. The halo highlight shows preserved scaffold atoms.
  2. Discuss the mechanism: which oxygen leaves as water? (Answer: the
     acid's hydroxyl oxygen; verified by the classic Roberts–Urey 1939
     ¹⁸O isotopic study.) BIME's MIN strategy with bipartite post-pass
     correctly reproduces this.
  3. Load Diels–Alder `C=CC=C.C=C>>C1CC=CCC1` — students predict the
     six-member ring closure, BIME visualises it.
  4. Load SN2 `CCBr.[OH-]>>CCO.[Br-]` — discuss inversion at the α-C.
- **Assessment**: predict the atom mapping for hydration, aldol, Wittig
  before running them in BIME.

### Graduate / research-level

#### Lesson: Maximum Common Substructure (MCS) for scaffold mining

- **Objective**: students understand how MCS is used to identify shared
  pharmacophores in drug families.
- **Time**: lab-style 2-hour session.
- **Activity**:
  1. Use **Examples → MCS** to compare β-blocker pairs (propranolol,
     atenolol, metoprolol).
  2. Compare statin pairs (atorvastatin, simvastatin, rosuvastatin).
  3. Discuss how the MCS reveals the pharmacophore vs. the
     SAR-modified periphery.
  4. Run a custom comparison: have each student pick two FDA-approved
     drugs from the same therapeutic class, predict the MCS, then
     verify.
- **Assessment**: write a one-page SAR analysis of a chosen drug
  family using BIME's MCS as the structural backbone.

---

## 4. Pre-built example reactions

The **Examples** page (`examples.html`) ships with ready-to-go
demonstrations of every BIME feature. Each example is interactive —
students can edit the input and see results update live.

| Example | Concept demonstrated | Suggested level |
|---|---|---|
| Benzene → chlorobenzene | Aromatic substitution | High school |
| CCO → CC=O | Oxidation, hydrogen accounting | High school / undergrad |
| Esterification | Reaction mapping, mechanism | Undergrad |
| Diels–Alder | Pericyclic, ring formation | Undergrad |
| SN2 inversion | Stereochemistry inversion | Undergrad |
| Aldol condensation | Multi-step mapping | Senior undergrad |
| Aspirin vs salicylic acid | MCS, scaffold + substituent | Undergrad |
| Caffeine vs theobromine | Purine pharmacophore | Undergrad |
| Ibuprofen vs naproxen | Arylpropionate scaffold | Senior undergrad |
| Custom SMARTS | Functional-group recognition | All levels |

---

## 5. Embedding BIME in your LMS

You can embed the workbench inside Moodle, Canvas, Blackboard, or any
LMS that supports HTML or `<iframe>`.

### Option A: Iframe (simplest)

```html
<iframe src="https://asad.github.io/bime-dist/workbench.html"
        width="100%"
        height="700"
        title="BIME molecular editor"
        loading="lazy"
        referrerpolicy="no-referrer">
</iframe>
```

Drop this into a Moodle "HTML block" or Canvas "Embed" page. Students
get the full workbench inside your course.

### Option B: Self-host on the LMS server

Many LMSes have a "static files" area. Upload the BIME folder there and
link to `workbench.html`. See [HOSTING.md](HOSTING.md) for the full
procedure.

### Option C: Direct link

If embedding isn't possible, simply link to the live demo or your
institution's hosted copy:

```markdown
[Open the BIME workbench](https://asad.github.io/bime-dist/workbench.html)
```

For per-question integration (LTI, scoring, gradebook), see
[EMBED.md](EMBED.md) — it covers the JavaScript API for capturing
student-drawn molecules from a parent page.

---

## 6. Printable worksheets

Open any molecule in BIME, then **File → Export SVG (print)**. The
print stylesheet:

- Uses CMYK-safe colours (no fluorescent screen pinks).
- Increases bond width for visibility on paper.
- Removes background gradients.
- Fits naturally onto A4 / Letter at 100% scale.

For a multi-molecule worksheet:

1. Load each molecule.
2. Export each as a numbered SVG.
3. Combine in your favourite document editor (LibreOffice, Word, LaTeX
   with TikZ-include, or Adobe Illustrator).

A pre-made worksheet template (`worksheet.svg`) is available in
`screenshots/` if you want a starting layout.

---

## 7. Assessment ideas

### Short answer

- Given a SMILES, draw the structure on paper.
- Given a structure, write the SMILES.
- Given a SMARTS, identify all matching atoms in a target.
- Given a reaction, predict the atom-atom mapping.

### Practical

- Students draw a molecule in BIME and submit the SMILES + an SVG
  export.
- Students compare two related drugs using MCS and write up the
  pharmacophore.
- Students design a SMARTS query that selectively matches a specific
  scaffold.

### Group work

- Each group is given a drug class (β-blockers, NSAIDs, statins,
  benzodiazepines). They use BIME to map the shared scaffold and the
  varying substituents, then present their findings.

---

## 8. Accessibility

BIME meets WCAG 2.1 AA accessibility standards:

- **Skip-to-content link** for keyboard users.
- **ARIA labels** on every interactive element.
- **Visible focus rings** on every focusable control (no `outline: 0`).
- **Keyboard shortcuts** for every drawing tool — drag-and-drop is
  supported but never required.
- **`prefers-reduced-motion`** is respected — animations disable when
  the OS asks.
- **`forced-colors: active`** — works in Windows High Contrast mode.
- **Print stylesheet** for handouts and exam papers.
- **Touch-friendly** — all controls are ≥ 24 × 24 CSS px (WCAG 2.5.5).
  Tested on iPad and Android tablets.

If you encounter an accessibility barrier, please open a bug report —
accessibility issues are treated as high-priority.

---

## 9. Common student questions

**"Why is benzene written `c1ccccc1` (lowercase) and not `C1=CC=CC=C1`?"**
The lowercase is SMILES's compact notation for aromatic atoms. Both
forms encode the same molecule; canonical SMILES prefers lowercase
because it captures the resonance-delocalised reality.

**"Why does my SMILES show a red warning on this atom?"**
Real-time valence checking. The atom has more bonds than its standard
valence allows. Either correct the structure, or annotate with explicit
charges (`[N+]`, `[O-]`) or radicals.

**"Where do I save my work?"**
File → Save (gives a `.mol` or `.smi` file). BIME does not store
molecules on any server — your work lives on your computer only.

**"Can I share my drawing with classmates?"**
Yes — copy the SMILES from the workbench's SMILES panel and paste it
into a chat message, email, or shared document. Anyone with the SMILES
can load the same structure in BIME.

**"Does BIME work offline?"**
Yes. Once the page is loaded, BIME makes zero network calls. Try
disabling Wi-Fi after loading.

---

## 10. Citing BIME in teaching materials

If you mention BIME in a syllabus, lecture notes, or textbook, please
cite:

> Rahman, S. A. *BIME: BioInception Molecular Editor.* BioInception
> PVT LTD, Cambridge, UK (2026). https://github.com/asad/bime

For a formal academic citation, use the metadata in
[CITATION.cff](../CITATION.cff). For BibTeX, see the bottom of
[README.md](../README.md).

The underlying algorithms are also worth citing in research-track
courses — see the **Citations** section of the README for RDT, EC-BLAST,
SMSD, and the Munkres–Kuhn assignment paper.

---

## More resources

- [USAGE.md](USAGE.md) — comprehensive user guide for advanced features.
- [STUDENTS.md](STUDENTS.md) — share this with your students.
- [EMBED.md](EMBED.md) — embedding into your LMS or course site.
- [HOSTING.md](HOSTING.md) — running BIME on your institution's
  infrastructure.

Questions or classroom case studies you'd like to share? Open an issue
at https://github.com/asad/bime-dist/issues or email
**asad.rahman@bioinceptionlabs.com** with "Educator" in the subject.

We love hearing how BIME is being used in classrooms around the world.
