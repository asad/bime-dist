# BIME for Researchers

A guide for cheminformatics, medicinal-chemistry, and chemical-biology
researchers using BIME in production workflows.

BIME is a 2D molecule editor and cheminformatics workbench with a
**fully exposed JavaScript API**, deterministic algorithms, and
publication-grade chemistry. It runs in the browser with zero
backend — making it ideal for embedded research dashboards,
reproducible analysis, and pipeline-friendly batch operations.

---

## Table of contents

1. [What BIME is good for in research](#1-what-bime-is-good-for-in-research)
2. [Programmatic API — quick tour](#2-programmatic-api--quick-tour)
3. [Atom-atom mapping (AAM)](#3-atom-atom-mapping-aam)
4. [Maximum Common Substructure (MCS)](#4-maximum-common-substructure-mcs)
5. [Substructure search (VF2 + SMARTS)](#5-substructure-search-vf2--smarts)
6. [Canonical SMILES + CIP stereochemistry](#6-canonical-smiles--cip-stereochemistry)
7. [Bond-change annotation](#7-bond-change-annotation)
8. [Building reproducible pipelines](#8-building-reproducible-pipelines)
9. [Determinism and reproducibility](#9-determinism-and-reproducibility)
10. [Citing BIME and its underlying algorithms](#10-citing-bime-and-its-underlying-algorithms)
11. [Comparison with other tools](#11-comparison-with-other-tools)

---

## 1. What BIME is good for in research

| Workflow | BIME's strength |
|---|---|
| **Reaction analysis at scale** | Pure-JS RDT port (Rahman 2016). Process thousands of reactions in a Node script with deterministic output. |
| **Scaffold mining via MCS** | Multi-strategy MCS pipeline (label-frequency / NLF / McSplit / Bron-Kerbosch / McGregor cascade). Tanimoto, Simpson, RASCAL similarity exposed. |
| **SMARTS-driven filtering** | VF2++ matcher with full SMARTS spec (charges, isotopes, ring membership, recursive SMARTS). |
| **Reproducible figures** | SVG/PNG export at any DPI, CMYK-safe print mode. Embeddable in LaTeX, Quarto, or web supplements. |
| **Bond-change classification** | EC-BLAST-style (Rahman 2014) `formed` / `cleaved` / `orderChange` / `stereoChange` / `hydrogenChange` events for every reaction. |
| **In-browser pipelines** | Jupyter notebooks (via JupyterLite or `%%javascript`), Observable, and any web app can call BIME's API directly. No Python install, no RDKit build. |

---

## 2. Programmatic API — quick tour

BIME's runtime is a set of `<script>` tags, each exposing a global
namespace. From any page or notebook with BIME loaded:

```js
// 1. Parse a SMILES into a Molecule
var m = SmilesParser.parse('CC(=O)Oc1ccccc1C(=O)O');     // aspirin

// 2. Generate canonical SMILES
var canon = SmilesWriter.write(m);
console.log(canon);   // → 'CC(=O)Oc1ccccc1C(=O)O' or equivalent canonical

// 3. Assign CIP stereochemistry
CipStereo.assign(m);
m.atoms.filter(function(a){ return a.cipLabel; })
       .map(function(a){ return [a.symbol, a.cipLabel]; });
       // → [ ... ]

// 4. SMARTS substructure match
var query = SmartsParser.parse('[CX3](=O)O');             // carboxylic acid
var matches = SmartsMatch.match(m, query);
console.log(matches.length);                              // → 1

// 5. MCS between two molecules
var a = SmilesParser.parse('CC(=O)Oc1ccccc1C(=O)O');     // aspirin
var b = SmilesParser.parse('Oc1ccccc1C(=O)O');           // salicylic acid
var mcs = SMSDMCS.mcs(a, b);
console.log(mcs.size);                                    // → 9 (atoms in shared scaffold)

// 6. Atom-atom mapping for a reaction
var rxn = SmilesParser.parse('CCO>>CC=O');                // ethanol oxidation
var aam = RDT.mapReaction(rxn);
console.log(aam.bondChanges);
// → [{type:'orderChange',...}, {type:'hydrogenChange',deltaH:-1,...}, ...]
```

The full JS API is in [USAGE.md](USAGE.md). Every public function has
JSDoc comments in the source files.

---

## 3. Atom-atom mapping (AAM)

BIME ships a faithful pure-JS port of the Reaction Decoder Tool
(Rahman *et al.*, *Bioinformatics* 32(13):2065–66, 2016).

### Default behaviour (recommended)

```js
var rxn = SmilesParser.parse('CC(=O)O.OCC>>CC(=O)OCC.O');
var result = RDT.mapReaction(rxn);

console.log(result.strategy);        // → 'MIN' / 'MAX' / 'MIXTURE' / 'RING'
console.log(result.score);           // → fitness, lower is better
console.log(result.mapping);         // → { reactantAtomId: productAtomId, ... }
console.log(result.bondChanges);     // → array of events, see §7
console.log(result.bipartiteApplied);// → true if Munkres-Kuhn post-pass refined
```

### Disabling components

For benchmark reproducibility against the v1.2.x default-off baseline:

```js
RDT.mapReaction(rxn, { useBipartitePostPass: false });
```

To suppress hydrogen-change events (heavy-atom-only output):

```js
RDT.mapReaction(rxn, { includeHydrogens: false });
```

To suppress stereo-change events:

```js
RDT.mapReaction(rxn, { includeStereo: false });
```

### Strategy-specific calls

```js
RDT.runMinPairwise(rxn);       // bond-change minimisation
RDT.runMaxPairwise(rxn);       // maximal MCS coverage
RDT.runMixturePairwise(rxn);   // hybrid
RDT.runRingPairwise(rxn);      // ring-preservation focus
```

Each strategy returns the same shape; BIME's default
`RDT.mapReaction()` runs all four and returns the best by fitness.

---

## 4. Maximum Common Substructure (MCS)

BIME's MCS pipeline cascades through 7 increasingly expensive levels
(L0–L5 + final), each gated by an upper-bound check that prunes
hopeless candidates early.

```js
var a = SmilesParser.parse('Cc1ccc(C(C)C(=O)O)cc1');     // ibuprofen-shaped
var b = SmilesParser.parse('COc1ccc2cc(C(C)C(=O)O)ccc2c1'); // naproxen
var mcs = SMSDMCS.mcs(a, b);

console.log(mcs.size);              // → atoms in MCS
console.log(mcs.atomMap);           // → reactant → query atom map
console.log(mcs.bondMap);           // → reactant → query bond map
console.log(mcs.tanimoto);          // → structural Tanimoto: |MCS| / (|A|+|B|-|MCS|)
console.log(mcs.simpson);           // → |MCS| / min(|A|, |B|)
```

### Tuning MCS

```js
SMSDMCS.mcs(a, b, {
  minSize: 6,                  // advisory lower bound
  ringsOnly: false,            // restrict to ring atoms
  tautomerAware: true,         // try tautomer pairs (slower, more permissive)
  timeoutMs: 10000             // budget per call
});
```

For batch MCS over a fingerprint-clustered library, see
`SMSDBatch.batch(...)` in [USAGE.md](USAGE.md).

---

## 5. Substructure search (VF2 + SMARTS)

```js
var target = SmilesParser.parse('CC(=O)Oc1ccccc1C(=O)O');
var query  = SmartsParser.parse('[CX3](=O)O');

// All matches:
var matches = SmartsMatch.match(target, query);
console.log(matches.length);                 // → 1
console.log(matches[0]);                     // → array of atom IDs

// Boolean check (cheaper):
var hit = SmartsMatch.hasMatch(target, query);
```

VF2++ handles the full SMARTS spec including:

- Atom-property and bond-property constraints (`[NX3;H2]`, `[#6;r6]`).
- Logical operators (`[N,O]`, `[C&D2]`, `[!#6]`).
- Recursive SMARTS (`[$(c1ccccc1)]`).
- Bond expressions (`/`, `\`, `~`, `@`).
- Charge / isotope / ring membership.

Performance on the BIME suite: 16 SMARTS-match tests pass in < 5 ms
end-to-end against medium drug-like molecules.

---

## 6. Canonical SMILES + CIP stereochemistry

```js
var m = SmilesParser.parse('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O');  // glucose
var canon = SmilesWriter.write(m);

// CIP tetrahedral and E/Z:
CipStereo.assign(m);
m.atoms.forEach(function(a){
    if (a.cipLabel) console.log(a.symbol, a.id, '→', a.cipLabel);  // R / S
});
m.bonds.forEach(function(b){
    if (b.cipLabel) console.log('bond', b.id, '→', b.cipLabel);    // E / Z
});
```

The canonical writer uses Morgan (extended connectivity) ranks +
Hückel 4n+2 aromaticity perception. Idempotent across multiple
re-parses; deterministic across runs.

Known limitation: stereo round-trip drift on polycyclic
multi-stereocentre molecules can flip CIP letters across re-parses
while preserving the atom skeleton. The drift is documented and
asserted-against in `tests/test_round_trip.js`. Avoid relying on
exact `[C@H]` / `[C@@H]` round-tripping for fused-ring systems with
stereocentres on the ring junctions.

---

## 7. Bond-change annotation

```js
var rxn = SmilesParser.parse('CCO>>CC=O');           // oxidation
var result = RDT.mapReaction(rxn);
result.bondChanges.forEach(function(ev) {
    console.log(ev.type, ev);
});
```

Event types:

| Type | When emitted | Key fields |
|---|---|---|
| `formed` | A bond exists in the product but not in the reactant. | `productBond`, `atom1`, `atom2`, `mapNumber1`, `mapNumber2` |
| `cleaved` | A bond exists in the reactant but not in the product. | `reactantBond`, `atom1`, `atom2`, `mapNumber1`, `mapNumber2` |
| `orderChange` | The bond order changes (e.g. 1 → 2). | `before`, `after`, `atom1`, `atom2` |
| `stereoChange` | A mapped atom's CIP label flips (R↔S, E↔Z). | `before`, `after`, `atom`, `productAtom` |
| `hydrogenChange` | A mapped atom's implicit-H count differs. | `deltaH`, `beforeH`, `afterH`, `atom`, `productAtom` |

This is the same event taxonomy used in EC-BLAST (Rahman 2014) and is
suitable for downstream classification, pathway-similarity scoring, or
educational visualisation.

---

## 8. Building reproducible pipelines

BIME's editor modules also work under Node.js (with the included test
shim) for headless pipeline use.

### Minimal Node example

```js
// pipeline.js
var shim = require('./tests/shim.js');
shim.loadAll();

var fs = require('fs');
var inputSmiles = fs.readFileSync('input.smi', 'utf8').split('\n');

inputSmiles.forEach(function(smi) {
    if (!smi.trim()) return;
    var m = SmilesParser.parse(smi);
    var canon = SmilesWriter.write(m);
    var hasAcid = SmartsMatch.hasMatch(m, SmartsParser.parse('[CX3](=O)O'));
    process.stdout.write([smi, canon, hasAcid].join('\t') + '\n');
});
```

```bash
node pipeline.js < input.smi > output.tsv
```

### Reactions pipeline

```js
var lines = fs.readFileSync('reactions.smi', 'utf8').split('\n');
lines.forEach(function(rxnSmi) {
    if (!rxnSmi.trim()) return;
    var rxn = SmilesParser.parse(rxnSmi);
    var aam = RDT.mapReaction(rxn);
    var summary = {
        smiles: rxnSmi,
        strategy: aam.strategy,
        score: aam.score,
        nFormed: aam.bondChanges.filter(function(e){return e.type==='formed';}).length,
        nCleaved: aam.bondChanges.filter(function(e){return e.type==='cleaved';}).length,
        nOrderChange: aam.bondChanges.filter(function(e){return e.type==='orderChange';}).length,
        nStereoChange: aam.bondChanges.filter(function(e){return e.type==='stereoChange';}).length,
        nHydrogenChange: aam.bondChanges.filter(function(e){return e.type==='hydrogenChange';}).length
    };
    console.log(JSON.stringify(summary));
});
```

### Browser-side batch

For client-side dashboards, just iterate over an array of SMILES and
call the same APIs. The 410-test regression suite runs in ~430 ms,
giving you a useful upper bound on per-call cost.

---

## 9. Determinism and reproducibility

BIME is fully deterministic for any given input + version + options:

- **No `Math.random`** in algorithm code.
- **No `Date.now`** in the pipelines.
- **Tie-breaks are explicit** — when two strategies post equal fitness,
  BIME breaks ties via the documented `TIEBREAK_ORDER` constant.
- **Munkres-Kuhn returns a unique optimum** under the documented
  column-scan tie-break.
- **The `useBipartitePostPass` strict-improvement gate** (`bipTotal >
  greedyTotal`) means coverage is monotone non-decreasing relative to
  the v1.2.x default-off baseline.

For a reproducible publication, pin the BIME version explicitly:

```bash
# Either pin the tag in your repo's lockfile / submodule:
git submodule add -b v1.4.4 https://github.com/asad/bime-dist.git vendor/bime

# Or vendor the dist/ folder directly into your supplementary materials:
cp -r bime-dist/dist supplementary/bime-1.4.4-dist
```

Then verify integrity in your build:

```bash
cd vendor/bime/dist && shasum -a 256 -c MANIFEST.sha256
```

---

## 10. Citing BIME and its underlying algorithms

When BIME is used in published research, please cite:

```bibtex
@software{bime2026,
  author       = {Rahman, Syed Asad},
  title        = {{BIME}: {BioInception} {Molecular} {Editor}},
  year         = {2026},
  publisher    = {BioInception PVT LTD},
  address      = {Cambridge, UK},
  url          = {https://github.com/asad/bime-dist}
}
```

The underlying algorithms have their own citations — please include the
relevant ones for the methods you use:

**Reaction Decoder Tool (atom-atom mapping):**
Rahman SA, Torrance G, Baldacci L, Cuesta SM, Fenninger F, Gopal N,
Choudhary S, May JW, Holliday GL, Steinbeck C, Thornton JM. *Reaction
Decoder Tool (RDT): extracting features from chemical reactions.*
**Bioinformatics** 32(13):2065-66 (2016).

**EC-BLAST (bond-change classification):**
Rahman SA, Cuesta SM, Furnham N, Holliday GL, Thornton JM. *EC-BLAST:
A Tool to Automatically Search and Compare Enzyme Reactions.* **Nature
Methods** 11:171-174 (2014).

**SMSD (MCS toolkit):**
Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM. *Small
Molecule Subgraph Detector (SMSD) toolkit.* **J. Cheminform.** 1:12
(2009).

**Munkres-Kuhn assignment (component pairing):**
Kuhn HW. *The Hungarian method for the assignment problem.* **Naval
Research Logistics Quarterly** 2:83-97 (1955); Munkres J. *Algorithms
for the assignment and transportation problems.* **SIAM Journal**
5:32-38 (1957).

For full BibTeX entries, see [CITATION.cff](../CITATION.cff) and the
README.

---

## 11. Comparison with other tools

| Need | BIME | RDKit | OpenBabel | ChemDraw | Marvin |
|---|---|---|---|---|---|
| **Pure browser, no install** | ✅ | ⚠️ via Pyodide / wasm builds | ⚠️ via wasm | ❌ | ❌ |
| **Apache 2.0 / FOSS** | ✅ | ✅ BSD | ✅ GPL-2 | ❌ | ❌ |
| **Atom-atom mapping** | ✅ RDT port | ⚠️ via separate RDT bridge | ❌ | ❌ | ✅ |
| **MCS** | ✅ multi-strategy | ✅ FMCS | ⚠️ basic | ❌ | ✅ |
| **SMARTS** | ✅ full spec | ✅ | ✅ | ⚠️ | ✅ |
| **CIP stereochemistry** | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Reaction bond-change events** | ✅ EC-BLAST style | ❌ | ❌ | ❌ | ⚠️ |
| **Embeddable in LMS / web** | ✅ | ❌ (server only) | ❌ | ❌ | ⚠️ commercial |
| **Per-seat cost** | $0 | $0 | $0 | $$$ | $$$ |

BIME does **not** replace RDKit for heavy server-side cheminformatics —
RDKit's force fields, conformer generation, fingerprint variety, and
descriptor catalogue go far beyond BIME's scope. BIME is the right
choice when you need **interactive 2D editing**, **embedded teaching
tools**, **reproducible browser-side analysis**, or **AAM and
bond-change classification with no server**.

The two tools complement each other beautifully — many published
pipelines use RDKit for backend processing and BIME for the
researcher-facing interactive layer.

---

## More resources

- [USAGE.md](USAGE.md) — full programmatic API.
- [HOSTING.md](HOSTING.md) — for deploying BIME on lab infrastructure.
- [EMBED.md](EMBED.md) — for embedding in dashboards or notebooks.
- [GitHub Issues](https://github.com/asad/bime-dist/issues) — bug reports
  and feature requests.

We particularly welcome research-track feature requests: if BIME is 80 %
of the way to your workflow, tell us about the missing 20 %.
