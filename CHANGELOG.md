# Changelog

All notable changes to BIME are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.2] — 2026-05-02

Distribution-only release. Adds a self-contained, zero-dependency build
pipeline that produces tamper-detectable, smaller deployment bundles
without altering the canonical `editor/*.js` source tree. No public-API
changes, no test-suite changes — the existing 410 regression tests pass
identically against the source files and against the new bundle.

### Added — build pipeline

- **`tools/build.js`** — pure-Node build script (no npm dependency).
  Concatenates `editor/*.js` in the canonical browser load order, writes
  `dist/bime.js` (unminified, source-readable) and `dist/bime.min.js`
  (light-minified — whitespace + comment strip only, no symbol mangling).
  Computes SHA-256 manifest and SRI sha384 strings, then re-runs the
  full regression suite against the bundle to verify behaviour
  equivalence. Single command: `node tools/build.js`.
- **`tools/sign-release.sh`** — opt-in GPG release-tag signing helper.
  Prints install hints and exits cleanly if `gpg` is missing, so this
  never blocks a build that simply does not have GPG configured.
- **`tests/run-against-bundle.js`** — runs the existing 18-file
  regression suite against the bundle (`dist/bime.js` or
  `dist/bime.min.js` with `--min`). Used by `tools/build.js` as a
  hard-gate before declaring a successful build.

### Added — distribution artefacts

- **`dist/bime.js`** — concatenated, unminified bundle (~1.0 MB). Source
  remains fully readable; per-file traceability dividers
  (`/*--- editor/<filename>.js ---*/`) keep the origin of each region
  obvious.
- **`dist/bime.min.js`** — light-minified bundle (~480 KB, ~48% of
  source sum). Whitespace collapsed, comments stripped (the Apache 2.0
  bang-comment banner is preserved). **No symbol mangling, no
  control-flow flattening, no string encryption.** The minified file
  produces byte-identical behaviour to the unminified bundle on the
  full 410-test regression suite.
- **`dist/MANIFEST.sha256`** — `shasum -a 256 -c` compatible manifest.
- **`dist/SRI.txt`** — modern SRI sha384 strings ready to paste into
  `<script integrity="...">` attributes.

### Added — documentation

- README "Build and verify" section covering rebuild, manifest
  verification, SRI usage and tag-signature verification.
- `workbench.html` and `examples.html` carry a commented-out alternative
  `<script src="dist/bime.min.js" integrity="...">` block beneath the
  existing `editor/*.js` tags. The default load path is unchanged —
  the bundle is opt-in.

### No obfuscation

This release explicitly **does not** ship an obfuscated bundle. Anyone
can read `editor/*.js` on GitHub anyway, and BIME is Apache 2.0 +
deliberately educational. Tamper-detection is provided by SRI hashes,
the SHA-256 manifest, and (opt-in) signed release tags — not by
hiding the source.

### Changed

- `editor/SMSDVersion.js`, `versions.json`, `CITATION.cff` bumped to
  `1.1.2`.

## [1.1.1] — 2026-05-02

Test-and-example expansion. Five core capabilities each gain a
dedicated regression-test file plus an interactive section in
`examples.html` so students and teachers can discover and verify
them in the browser. Additive only — no API changes, no removals.

### Added — tests (+106 cases)

- **`tests/test_aam.js`** — 25 cases covering atom-atom mapping via
  `RDT.mapReaction`: per-strategy mapping size, strategy
  comparison, determinism over 5 runs, `options.timeoutMs`,
  pre-mapped-atom preservation, empty-side handling, multi-component
  reactions, explicit-hydrogen reactions.
- **`tests/test_mcs_extended.js`** — 15 cases extending the existing
  MCS suite: real-drug pairs (aspirin / salicylic acid, ibuprofen /
  naproxen, caffeine / theobromine), disconnected query, timeout,
  bond-mapping consistency, minSize, Tanimoto from MCS, MCS
  symmetry across 5 pairs, single-atom edge cases.
- **`tests/test_substructure_extended.js`** — 20 cases for
  substructure search: phenyl in toluene / ethylbenzene /
  diphenylmethane, carboxylic-acid SMARTS, amide bond, aromatic-N,
  sulfonamide, recursive carbonyl SMARTS, ring-constraint, negative
  cases, multi-mapping, performance on a 50-atom molecule.
- **`tests/test_bondchange.js`** — 15 cases on bond-change semantics
  via `RDT.annotateBondChanges`: carbonyl reduction, hydration,
  hydrogenation, esterification, hydrolysis, SN2, Diels-Alder,
  no-change reaction, isomerization, ring-opening, ring-closing,
  stereo flip, deterministic event ordering, mirrored-reaction
  symmetry.
- **`tests/test_canon.js`** — 31 cases on canonical SMILES via
  `SmilesWriter.write`: idempotence (10 textbook molecules),
  equivalence under permuted input (6 pairs), determinism across
  5 runs (3 cases), aromatic preference (3 cases), charge / isotope
  / stereo round-trip canonical (3 cases), empty / single-atom edge
  cases.

### Added — interactive examples in `examples.html`

Five new full-width sections (one per capability), each with a
heading, descriptive paragraph, editor canvas, 4-5 example pill
buttons that load pre-canned inputs, and a `<details>` "How does
this work?" block with a 2-3 sentence algorithm summary:

1. **Atom-Atom Mapping (AAM)** — `#aam`. Ethanol oxidation, Fischer
   esterification, SN2, Diels-Alder, aromatic chlorination.
2. **Maximum Common Substructure (MCS)** — `#mcs`. Aspirin /
   salicylic acid, ibuprofen / naproxen, caffeine / theobromine,
   glucose / fructose. Two side-by-side editors; MCS atoms
   highlighted teal.
3. **Substructure Search (SUB) — SMARTS** — `#sub`. Carboxylic acid,
   amide bond, 6-ring, benzene, recursive carbonyl. Match-count
   display + atom highlighting.
4. **Bond Changes (BC)** — `#bc`. Esterification, reduction,
   Diels-Alder, ring-opening. "Auto-map" → "Show Bond Changes"
   pipeline; per-event listing.
5. **Canonical SMILES (Canon)** — `#canon`. Multiple writings of
   the same molecule converge to one canonical form. Includes a
   paired-test button that asserts `canon('CCO') === canon('OCC')`.

A "References" section was added at the bottom of `examples.html`
citing Rahman et al. SMSD 2009 / RDT 2016 / EC-BLAST 2014.

### Tests

- **410 regression tests pass** (was 304 in v1.1.0; +106 new). Plain
  Node, zero deps, ~280 ms total.
- Aromatic audit: 0 issues across 545 molecules.
- All `editor/*.js` parse-clean.

### Known-not-perfect cases pinned (deferred to future)

- Cross-component bond-formed events in `RDT.annotateBondChanges`
  Pass 2 (currently 0 formed events on couplings spanning two
  reactant fragments).
- `stereoChange` events depend on prior CIP perception — `mapReaction`
  does not run CIP in v1.1; flips currently emit 0 events.
- Polycyclic multi-stereocentre canonical drift (`@`/`@@` flip on
  Morphine, DHEA, glucose-with-stereo) — pinned in
  `KNOWN_CANON_DRIFT` allow-list and Canon tests use no-stereo
  variants for equivalence.
- `mcsOpts.minSize` is advisory (does not truncate engine output).
- MIN strategy on aromatic chlorination collapses to a 1-atom
  mapping by design (smallest-match policy of the published RDT
  algorithm).

## [1.1.0] — 2026-05-02

Pure-JavaScript port of the **Reaction Decoder Tool (RDT)** algorithm
(Rahman SA et al., *Bioinformatics* 32(13):2065-66, 2016) running
entirely in the browser. Restores atom-atom mapping (AAM) for drawn
reactions without any backend, completing the v1.0.3 "drop-anywhere
static site" promise.

### Added

- **`editor/RDT.js`** — pure-JS port of the published RDT 2016
  algorithm:
  - 4 mapping strategies: `MIN`, `MAX`, `MIXTURE`, `RING_PERCEPTION`.
  - Bond-change annotator (`formed` / `cleaved` / `orderChange` /
    `stereoChange`).
  - Fitness function: bond-change minimisation with ring-preservation
    bonus and aromatic-atom chemical filter.
  - Selector picks the strategy with lowest fitness; deterministic
    tie-break (fewer events → ring bonus → lex-min strategy name).
  - Honours `options.timeoutMs` with graceful partial-mapping return.
- **"Auto-map" toolbar button** in the Reaction toolbar group of
  `workbench.html` and `test.html`. Enabled when a reaction arrow is
  drawn; calls `RDT.mapReaction(this.molecule)` and re-renders the
  editor with `:N` map-number superscripts.
- **`tests/test_rdt.js`** — 32 regression tests across the four
  strategies, bond-change semantics, fitness behaviour, selector
  determinism, multi-component reactions, pre-mapped-atom
  preservation, unbalanced-reaction edge cases, and in-place mutation
  of the input molecule.

### Public API surface

```js
RDT.mapReaction(reaction, options)        // -> { mapping, bondChanges, score, strategy, timedOut, warnings, strategyResults }
RDT.runMinPairwise(reaction, options)
RDT.runMaxPairwise(reaction, options)
RDT.runMixturePairwise(reaction, options)
RDT.runRingPairwise(reaction, options)
RDT.annotateBondChanges(reactantSide, productSide, atomMap)
RDT.fitness(reactionOrSides, atomMap, bondChanges)
RDT.STRATEGIES   // ['MIN', 'MAX', 'MIXTURE', 'RING']
RDT.version      // '1.1.0'
```

### Tests

- 304 regression tests pass (`node tools/run-tests.js`) — was 272 in
  v1.0.3, +32 RDT cases.
- Aromatic audit: 0 issues across 545 molecules.

### Citations

- Rahman SA, Torrance G, Baldacci L, Martínez Cuesta S, Fenninger F,
  Gopal N, Choudhary S, May JW, Holliday GL, Steinbeck C, Thornton JM.
  *Reaction Decoder Tool (RDT): extracting features from chemical
  reactions.* **Bioinformatics** 32(13):2065-66 (2016).
  doi:10.1093/bioinformatics/btw096
- Rahman SA, Cuesta SM, Furnham N, Holliday GL, Thornton JM. *EC-BLAST:
  a tool to automatically search and compare enzyme reactions.*
  **Nature Methods** 11:171-174 (2014). doi:10.1038/nmeth.2803
- Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM. *Small
  Molecule Subgraph Detector (SMSD) toolkit.* **J. Cheminform.** 1:12
  (2009).

### Known limitations (deferred to v1.2.0)

- Implicit-hydrogen accounting in bond-change events (only heavy-atom
  bonds are walked).
- Greedy-only assignment (matches the published 2016 algorithm).
- EC-class hint derivation from bond-change profile (RDT 2016 has it;
  not yet ported).
- Stereo-change events depend on prior CIP perception.

## [1.0.3] — 2026-05-02

Pure-frontend release. Removed the optional Java backend services
(`rdt-service/`) and their launcher scripts. BIME now runs entirely
in the browser with zero external dependencies and zero external
network calls — the original "drop-anywhere static site" identity
restored.

### Removed

- `rdt-service/RDTServer.java`, `rdt-service/SMSDServer.java`
- `start-bime.sh`, `start-bime.bat`, `start-services.sh`
- `Dockerfile`, `docker-compose.yml`
- Backend-status UI panels in `workbench.html` and `test.html`.
- `MolEditor.prototype._safeRdtUrl`, `_safeSmsdUrl`,
  `_fetchWithFallback`, `_fetchWithTimeout`, `rdtServerUrl`,
  `smsdServerUrl`, `_resolvedRdtUrl`, `_resolvedSmsdUrl`,
  `_backendCheckTimer`, and related plumbing.
- "Auto-map (RDT)" toolbar button.

### Kept (back-compat shims)

- `getIUPACName()`, `getCommonNames()`, `_lookupIUPACName()`,
  `_cactusNameLookup()`, `_asyncCactusName()` remain as no-op
  stubs so v1.0.2 third-party integrations do not crash.

### Privacy / security

- CSP `connect-src` tightened from
  `'self' https://api.bioinceptionlabs.com http://localhost:8766 http://localhost:8767`
  to just `'self'` across all 6 HTML pages. The editor cannot make
  any external HTTP request, even by mistake.

### Roadmap

- v1.1.0 (planned): pure-JS atom-atom mapping (AAM) replacement.
  Educational reactions (esterification, SN2, Diels–Alder, Claisen,
  amide formation, hydration) covered in-browser with no backend.

### Tests

- 272 regression tests still pass (`node tools/run-tests.js`).
- Aromatic audit: 0 issues across 545 molecules.

## [1.0.2] — 2026-05-02

Maintenance release on top of v1.0.1. Bug fixes, robustness, and
accessibility polish. No breaking API changes, no new features. Public
Apache-2.0 source. Pure JavaScript + SVG, zero runtime dependencies.

### Chemistry — algorithm correctness

- **`editor/SMSDMCS.js`**: canonical-hash collisions are now resolved by
  isomorphism check. The previous version mapped two non-isomorphic graphs
  to the same hash bucket in rare cases, producing a missed MCS hit. Added
  per-call `timeoutMs` so MCS is bounded on pathological inputs.
- **`editor/SMSDVF2.js`**: subgraph-containment fast path now goes through
  the `ppx` (predecessor / successor) refinement so trivial-but-large
  queries return immediately. Atom-compatibility cross-check tightened
  to reject false matches across charged / radical / isotope variants.
- **`editor/SMSDBatch.js`**: per-call `timeoutMs` plumbed through the
  batch driver. Long-running pairs no longer stall an entire batch.

### Chemistry — aromaticity & charges

- **`editor/SmilesWriter.js`** + **`editor/SmartsMatch.js`**: aromaticity
  perception is now charge-aware. The cyclopentadienyl anion (Cp⁻) and
  the tropylium cation (C₇H₇⁺) are correctly written aromatic;
  pyridinium (NH⁺) round-trips as `[nH+]`; saturated cyclohexane is no
  longer mis-aromatised on round-trip.
- Saturated-carbon aromaticity rejection: a sp³ ring carbon with two
  hydrogens never participates in an aromatic ring even if its
  neighbours do. Closes a class of false-positive aromatic rings on
  partially saturated polycyclics.

### MDL / RXN / SDF I/O

- **Cross-platform line endings**: `_parseMolV2000`, `_parseMolV3000`,
  `_parseRxnFile`, and the SDF `<NAME>` regex now handle Unix `\n`,
  Windows `\r\n`, and classic-Mac `\r` line endings interchangeably.
  Files saved on one platform open cleanly on the others.
- **MDL property blocks**: `M  CHG`, `M  ISO`, and `M  RAD` are now
  read and written. Mass-difference field is honoured as a fallback
  isotope source. Charge code 4 in the V2000 atom block is correctly
  interpreted as a radical, not a charge.
- **V3000 robustness**: line continuation (`-`) is honoured;
  positional `aamap` (atom-atom-mapping) values are accepted; `$RXN
  V3000` is detected alongside the legacy `$RXN` header.
- **Metadata round-trip**: `mol.name`, `mol.program`, `mol.comment`,
  `> <NAME>`, and `> <COMMENT>` SDF fields preserve through V2000 →
  parse → V3000 → SDF → RXN round-trips. Reaction-header comment
  lines are preserved.
- **Atom-atom mapping (AAM) preservation**: `mapNumber` round-trips
  through MOL V2000 columns 60-62, V3000 `MAP=N`, and RXN. The
  `mapHighlighted` flag now survives undo/redo via `Molecule.toJSON`.
  SMARTS halos preserved via the `bgColor` field.
- **Stereo & CIP preservation in `Molecule.toJSON`**: `chirality`,
  `cipLabel`, `radical`, `aromatic`, `mapHighlighted`, and `bgColor`
  are now serialised so undo / redo no longer drops them after a
  long edit session.

### Editor — data model & memory

- **`editor/Molecule.js`** `clone()` and **`editor/MolEditor.js`**
  `_restoreMolecule` round-trip `reactionPlusSigns`, the arrow `type`,
  and the arrow `conditions` field. Reaction undo / redo now recovers
  the full reaction object, not just the molecule list.
- **`editor/History.js`**: ring-buffer (O(1) push and pop). Snapshots
  are nulled out on overwrite so the GC reclaims large molecules
  promptly. No more long GC pauses on big undo histories.
- **`editor/ImageExport.js`**: measurement-SVG leak fix carried forward;
  null-guards on bond / atom in `editor/Renderer.js` prevent NPE on
  stale renders during fast edits.

### Tests

- **New regression suite at `tools/run-tests.js`**: 272 Node.js tests
  across 12 files (`tests/test_smiles_parser.js`,
  `test_smiles_writer.js`, `test_smarts_match.js`,
  `test_substructure_vf2.js`, `test_mcs.js`, `test_aromaticity.js`,
  `test_hcount.js`, `test_round_trip.js`, `test_history.js`,
  `test_chem_torture.js`, `test_aam_round_trip.js`,
  `test_charge_iso_round_trip.js`). Plain Node, zero dependencies,
  runs in ~80 ms end-to-end. Run with `node tools/run-tests.js`.
- Coverage spans parser / writer round-trip, SMARTS matching, VF2
  subgraph isomorphism, MCS, Hückel aromaticity, hydrogen counting,
  charge / isotope / radical round-trip, atom-atom-mapping
  round-trip, "torture" inputs (zwitterions, polycycles, charged
  aromatics, reaction SMILES) and history (undo / redo with the
  v1.0.2 `reactionPlusSigns` fix).

### Accessibility (WCAG 2.1 AA)

- **Skip-links** on every public page (`index.html`, `workbench.html`,
  `examples.html`, `docs.html`, `screenshots.html`, `test.html`).
- **Landmarks**: `<main>` on every page; nav `aria-label="Primary"`.
- **Toolbar buttons**: `aria-pressed` / `aria-label` for state and
  affordance; SVG icons marked `aria-hidden focusable="false"`.
- **SVG canvas**: `role="img"` plus `<title>` and `<desc>` updated on
  each render so screen readers can identify the current structure.
- **Hamburger menu**: `aria-expanded` mirrored to button state,
  Escape-to-close, focus restored to the toggle on close.
- **Theme toggle**: `aria-pressed` reflects light / dark mode.
- **Touch targets**: WCAG 2.5.5 (44 × 44 CSS px minimum) across
  toolbars and the molecule library.
- **Contrast**: WCAG 1.4.3 AA verified across light and dark themes.

### Design system

- Global `:focus-visible` keyboard ring, distinct from mouse focus.
- Disabled-state styling for buttons that are not yet wired (e.g.
  AAM under construction).
- Print stylesheet: educational handouts print cleanly without nav,
  toolbars, or footers; bonds thicker, fonts larger.
- `prefers-reduced-motion` honoured (animations short-circuited).
- `forced-colors: active` (Windows High Contrast) — explicit
  `currentColor` / `CanvasText` fallbacks added.
- Tighter typography rhythm; `text-wrap: balance` on headings.

### Brand & assets

- `images/svg/smsd-logo.svg`: removed obsolete legacy wordmark; logo
  now ships as a clean glyph. Removed font-family references that
  required system-installed fonts.
- Footer year set to 2026 throughout; consistent BIME / BioInception
  casing across all pages.

### Tooling

- **`tools/audit-aromatic.js`** confirms 0 issues across all 545
  molecules in `common-molecules.js`. Now run as part of pre-commit.
- **`editor/SMSDVersion.js`** simplified — references only the
  upstream SMSD citation (Rahman et al., J. Cheminform. 2009;1:12)
  and the BIME release the file was bundled with.
  `versions.json` simplified to `{ bime: 1.0.2 }`.
- Specific upstream SMSD version numbers stripped from comments and
  UI labels (now just "SMSD"). Cite the published paper, not the
  implementation revision.

### Privacy / security (continuing the v1.0.1 hardening)

- **`js/nav.js`**: `localStorage.getItem` / `setItem` wrapped in
  `try / catch` against `SecurityError` in private-browsing mode and
  sandboxed iframes (theme-toggle no longer crashes init).
- **Java services** (`rdt-service/`): SMILES length cap, body cap,
  30-second process timeout, bounded `ThreadPoolExecutor` with
  `AbortPolicy` (overload returns 503, not OOM); env-driven CORS
  allowlist; `ThreadLocal` CDK parsers (CDK is not thread-safe);
  explicit `StandardCharsets.UTF_8` for child-process I/O.
  Re-confirmed in this release.

## [1.0.1] — 2026-05-01

Maintenance release. Bug fixes only — no breaking API changes, no new
features. Public Apache-2.0 source.

### Privacy & security

- **Removed third-party name-resolver dependency** (`cactus.nci.nih.gov`).
  Earlier versions silently forwarded every drawn SMILES to a third-party
  server for IUPAC / common-name lookup. v1.0.1 removes that channel entirely.
  Public methods (`getIUPACName`, `getCommonNames`,
  `_lookupIUPACName`, `_cactusNameLookup`, `_asyncCactusName`) are kept as
  no-op-compatible stubs so existing integrations do not crash. The built-in
  545-molecule local database remains the authoritative naming source.
- **Removed dead `editor/WasmBridge.js`** (391 LOC of unwired CheerpJ scaffold
  loading a third-party CDN script with no integrity check). Eliminates
  supply-chain risk on a code path that was never reachable.
- **Backend host allowlist** in `_safeRdtUrl()` / `_safeSmsdUrl()`. Previously
  validated only the protocol; now also requires the host to be one of
  `api.bioinceptionlabs.com`, `localhost`, `127.0.0.1`, `0.0.0.0`. Closes a
  potential SSRF / data-exfiltration channel if the URL is ever set from
  untrusted input.
- **Workbench molecule-library XSS hardening**: `buildMolBrowser` migrated
  from string-concatenated `innerHTML` (with insufficient `'` escape only)
  to DOM construction with `textContent` — defends against any future
  `common-molecules.js` entry that contains HTML-significant characters.
- **`document.write` removed from `editor/ImageExport.js`**. PDF print path
  now opens a Blob-URL document, eliminating cross-origin warnings and
  matching modern browser hardening guidance.
- **Content-Security-Policy** meta tag added to every HTML page
  (`workbench.html`, `index.html`, `examples.html`, `docs.html`,
  `screenshots.html`, `test.html`). Restricts `connect-src` to the BIME
  cloud and localhost; forbids `object-src`, `frame-ancestors`, and
  third-party `base-uri`.
- **`test.html` no longer bypasses `_safeRdtUrl()`** — RDT health check
  now goes through the same host allowlist as the rest of the editor.
- **External fetches gain an `AbortController` timeout** via
  `MolEditor.prototype._fetchWithTimeout` (default 8 s).
- **RDT/SMSD Java servers hardened** for non-localhost deployment:
  - SMILES length cap (`BIME_MAX_SMILES_LEN`, default 10 000 chars)
  - Request body cap (`BIME_MAX_BODY_BYTES`, default 64 KB)
  - 30 s process / operation timeout (`BIME_PROC_TIMEOUT_S`)
  - Bounded `ThreadPoolExecutor` (`BIME_POOL_CORE/MAX/QUEUE`,
    `AbortPolicy` returns 503 on overload instead of OOM)
  - CORS allowlist (`BIME_CORS_ORIGINS`, `*` only when unset for local dev)
  - `SMSDServer`: `SmilesParser` / `SmilesGenerator` made `ThreadLocal`
    (CDK parsers are not thread-safe; the prior fixed-pool of 4 raced)
  - `RDTServer`: explicit `StandardCharsets.UTF_8` for child-process I/O
    (avoids Cp1252 mangling on Windows hosts)
  - Process-leak guard in both servers (`destroyForcibly` in `finally`)

### Privacy posture (now accurate)

The previous SECURITY-AUDIT.md §3.1 claim of "zero external network calls"
was inaccurate as of v1.0.0. v1.0.1 restores that posture: the only external
host the public BIME page may contact is `api.bioinceptionlabs.com`, and only
when the user explicitly points the editor at a self-hosted or BioInception
RDT/SMSD endpoint. No molecule data is sent off-host by default.

### Chemistry data restoration

The 552-molecule "canonicalisation" pass in v1.0.0 (commit `984e5a2`) had
mass-lowercased atom symbols without verifying ring aromaticity, producing
chemically invalid encodings (aromatic ribose, aromatic peroxide,
aromatic β-lactams, aromatic cyclohexane, etc.). Plus several stereochemistry
losses had crept in.

This release re-verifies the entire database against PubChem-canonical
SMILES and adds a `tools/audit-aromatic.js` programmatic checker. After
fixes the audit reports **0 issues across 545 molecules**.

Notable named fixes (all now match PubChem canonical SMILES):

- **DHEA** was a byte-identical duplicate of Testosterone; both also had
  wrongly aromatic steroid B/C/D rings. Both replaced with PubChem
  CIDs 5881 / 6013, full stereo restored.
- **Morphine** SMILES contained `[c@]`/`[c@@H]` (chemically impossible —
  aromatic atoms cannot be tetrahedral stereocentres). Replaced with
  PubChem CID 5288826 canonical.
- **Codeine** had a 4-membered ring with aromatic atoms. Replaced with
  PubChem CID 5284371.
- **Atorvastatin** was missing (3R,5R) stereo on the dihydroxyheptanoate
  pharmacophore. Added.
- **Estradiol, Progesterone, Cholesterol** had stereo dropped in the
  canonicalisation pass. PubChem-canonical with full stereo restored.
- **Penicillin G/V, Amoxicillin** had β-lactam + thiazolidine rings written
  as aromatic. Replaced with correct fused bicyclic stereo.
- **Adenosine** had ribose written as aromatic. Replaced.
- **Quinine** had quinuclidine bridge written as aromatic. Replaced.
- **Menthol** was encoded as p-cymen-3-ol (a different molecule).
  Replaced with the correct cyclohexanol with `[C@@H]` stereo.
- **Artemether / Artemisinin** had the 1,2,4-trioxane endoperoxide written
  as aromatic. Replaced with correct trioxane bridge structure.
- **Esomeprazole** previously identical to Omeprazole. Sulfoxide stereo
  added (`[S@]`).
- **Saxagliptin, Budesonide** structurally rewritten to PubChem-canonical.
- **Catechin / Epicatechin** previously shared the same SMILES. Now
  distinguished by C-3 stereochemistry (2R,3S vs 2R,3R).
- **Inositol** was encoded as benzenehexol (aromatic, wrong). Replaced
  with myo-inositol (saturated cyclohexanehexol with full stereo).
- **All 19 chiral L-amino acids** now carry `[C@@H]` (or `[C@H]` for
  L-cysteine due to CIP priority flip).
- **Glucose, Sucrose** carry full Haworth-equivalent stereo annotation.
- **7 duplicate name-cased entries** removed (`Toluene`, `Ethanol`,
  `Methanol`, `Acetone`, `Acetic Acid`, `Urea`, `Glucose` had been listed
  in both the `basic` and `other` categories).

The dataset now ships with a header acknowledging PubChem provenance and a
`tools/audit-aromatic.js` programmatic checker so future regressions can
be caught in CI.

### Editor bug fixes (non-data)

- `editor/History.js`: replaced `Array.shift()` cap-at-100 (O(n) per push)
  with a fixed-capacity ring buffer (O(1) push/pop). Snapshots are nulled
  out on overwrite so the GC can reclaim large molecules.
- `editor/ImageExport.js`: `_measureSVG` measurement node was attached to
  `document.body` and never detached (memory + DOM leak). Rewritten with
  reference-counted attach/detach in `try/finally`. Public
  `ImageExport._cleanup()` added for tests/teardown.
- `editor/ImageExport.js`: `MolfileWriter.toMolfile(mol)` was an undefined
  method (would throw `TypeError`). Changed to `MolfileWriter.write(mol)`.
- `editor/ImageExport.js`: aromatic-circle loop guarded against null
  atoms (stale ring info after a delete previously NPE'd).
- `editor/Molecule.js`: `clone()` and `toJSON()` now serialise
  `reactionPlusSigns` and the arrow's `type` / `conditions` fields, so
  reaction undo/redo round-trips correctly.
- `editor/MolEditor.js`: `_restoreMolecule` now restores
  `reactionPlusSigns` (paired with the Molecule.js change above).
- `editor/Renderer.js`: `_drawAromaticCircle` null-guards atoms.

### SMARTS / SMILES parser & writer

- `editor/SmartsParser.js`:
  - Negative-charge primitive accepted only when followed by digit / `]`
    / EOS — `[C-;R]` silently dropped the `-`. Now accepts at any
    bracket-level terminator.
  - Multi-charge `++` / `+++` (and `--` / `---`) folded into a single
    charge magnitude, mirroring the previous-handled `-` case.
  - Top-level negated bond (`!-`, `!=`, `!:`, `!#`, `!~`) no longer
    silently swallows the follow-on character on invalid input.
- `editor/SmartsWriter.js`: `canWriteUnbracketed` now respects the
  `negate` flag on aromatic / aliphatic constraints (`[!a]` was being
  written unbracketed, dropping the `!`).
- `editor/SmartsMatch.js`: documentary `// TODO(v2.0.0):` comment added at
  `perceiveAromaticity` flagging the divergence from `SmilesWriter.js`
  (PREV+NEXT vs NEXT-only ring-bond inspection). Unification deferred
  to v2.0.0 as it touches algorithmic surface.

### HTML pages

- All `<a target="_blank">` links upgraded to
  `rel="noopener noreferrer"`.
- `js/nav.js`: `localStorage.getItem`/`setItem` wrapped in
  `try/catch`; corrupt or non-permitted values fall back to defaults
  (private browsing / sandboxed iframe no longer crash theme init).

### Legal / licensing

- **`NOTICE` rewritten**: dropped the Apache-incompatible
  "you MUST cite" / "you MUST retain" mandates in favour of a courteous
  citation request. Removed the inline "NOVEL ALGORITHMS" self-disclosure
  block (gratuitous prior-art surface with no licence benefit).
- **`CITATION.cff`** added (replaces the citation block formerly in
  NOTICE).
- **Apache 2.0 header added** to `common-molecules.js` and `js/nav.js`
  (previously missing per §4(c)).
- `common-molecules.js` now carries an explicit provenance comment
  acknowledging PubChem as the source of the SMILES (per
  *Feist v. Rural*; SMILES are factual representations, not
  copyrightable).

### Tooling / build

- `Dockerfile` now installs `curl` so the healthcheck actually runs (the
  base `openjdk:11-jre-slim` image does not ship `curl`).
- `start-bime.sh`, `start-services.sh`: PIDs in trap/cleanup are quoted.
- `tools/audit-aromatic.js`: added (Hückel + ring-size + peroxide
  consistency checker; suppresses the most common false-positive
  patterns: tetrazole / pyrrole-class lone-pair donation, fused-ring
  traversal artefacts, casing-only name duplicates).

### Removed (since the prior public commit)

- `editor/WasmBridge.js` — dead CheerpJ scaffold (never referenced).
- 7 duplicate molecule entries from `common-molecules.js`.
- All four `cactus.nci.nih.gov` external fetch sites.
- `_enableExternalLookup` opt-in flag (no longer needed; lookup removed).
- `MolEditor.prototype.moleculeName` was retained but writes now have no
  network side-effect.

## [1.0.0] — 2026-03-29

Initial public release.
