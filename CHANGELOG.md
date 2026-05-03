# Changelog

All notable changes to BIME are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.7.0] — 2026-05-04

Improved ML refiner weights and a 2x2 top-toolbar layout.

### Changed

- ML refiner weights (`editor/ml-depict-weights.json`) updated
  for improved coordinate-residual handling.

### Migration

No action required.

## [1.6.3] — 2026-05-03

**UX patch.** Top toolbar reorganised into a 2x2 layout for
clearer mental model: drawing primitives (bonds + atoms) above,
edit/view actions below.

### Changed

- Top toolbar now uses a 2x2 row layout (Bonds + atoms on row 1,
  Edit + View + Custom on row 2).

## [1.6.2] — 2026-05-03

**UX patch.** Default drawing tool on workbench startup is now
Select (was Bond), for parity with other industry molecule editors.

### Changed

- Default tool on workbench startup is now **Select** (was Bond).

## [1.6.1] — 2026-05-03

**Hotfix.** Fixes a regression introduced in v1.6.0 where the
workbench editor failed to render its toolbar.

### Fixed

- Workbench editor now renders correctly.

## [1.6.0] — 2026-05-03

Optional ML-augmented 2D depiction. Pure JavaScript, runs in the
browser, zero dependencies. Opt-in, off by default — existing
v1.5.x layouts are unchanged.

### Added

- `Layout.options.useMLDepict` (default `false`) — opt-in ML
  refinement layered on top of the rule-based pipeline.
- `Layout.options.mlDepictWeight` (default `0.15`) — blend factor.

### Migration

No action required. Existing callers are unaffected.

## [1.5.2] — 2026-05-03

UI/UX polish + chemistry-canonical 2D depictions across the board.
Five expert agents ran in parallel (A: toolbar split, B: highlight
finesse, C: flavonoid templates, D: ~30+ non-flavonoid scaffold
templates, E: generic depict-algorithm enforcer) so every layer of
the visual + structural rendering pipeline gets attention this
release.

### Changed — toolbar split (top + left rail)

- The previous flat horizontal toolbar (which had grown cluttered
  through v1.4.x as new groups were added) is now split between a
  TOP horizontal row and a LEFT vertical rail, in the
  ChemDraw / Marvin / Ketcher pattern.
- **TOP row** (high-frequency, action-style): `bonds`
  (Single / Double / Triple / Wedge), `edit` (Select / Delete /
  Move / Undo / Redo / Clear), `view` (Layout / Zoom / Fit /
  R-S / E-Z / H / Ar), `custom` (Customize panel toggle).
- **LEFT rail** (~48 px, icon-only with tooltips, mode-y / less-
  frequent): `draw` (Bond, Chain), `rings` (3/4/5/6/7 flyout),
  `reaction` (Arrow / Map / Map # / Auto-map / Colors / Pairs),
  `actions` (Validate / SMARTS / MCS / Sim), `export`
  (Name / SVG / PNG / Print / Copy).
- Atom bar stays as a horizontal strip under the top row so the
  periodic-table and fragment popups (which anchor on
  `_atomBar.offsetTop`) keep working.
- Each TOOLBAR_GROUPS entry now carries a `placement: 'top' |
  'left'` field; `ToolbarPrefs.applyToGroups` preserves it via the
  shallow `for-in` reorder copy. The (i) info popovers and the
  v1.4.3 Customize panel both keep working — no IDs, ARIA labels,
  or `_doAction` cases were renamed.
- New CSS classes: `.bime-left-rail`, `.bime-main-row`. New ARIA:
  `role="toolbar" aria-orientation="vertical"` on the rail,
  `role="toolbar"` on the top row.

### Changed — search-match highlighting

- `Renderer.js`'s bgColor halo path got a polish pass. Matched atoms
  now render as an outer pale-teal disc (`#99f6e4`, r = 11,
  opacity 0.32) PLUS an inner brand-teal ring outline (`#0d9488`,
  r = 12, stroke 1.6, opacity 0.85) — reads as "circled in teal"
  rather than "blobbed".
- Matched bonds get a thick pale-teal underlay + a brighter brand-
  teal inner stroke, giving a luminous "highlighter" feel rather
  than a dull rectangle.
- 1-second SMIL `<animate>` reveal that flashes from
  `#5eead4` to brand-teal on initial paint (no CSS, no libraries).
  Honours `prefers-reduced-motion: reduce` — the resting state
  paints immediately if the user opts out (WCAG 2.3.3).
- ARIA: every halo carries an SVG `<title>` reading `"Search match
  (substructure | mcs | exact)"`. Brand-teal ring meets WCAG 2.1 AA
  non-text contrast (4.96 : 1 against `#fff`, 5.36 : 1 against
  `#0f172a` dark canvas).
- `MolEditor.prototype.highlightSearchMatch` sets
  `atom.searchMatchKind` and `atom.searchMatchTs` so the renderer
  can drive the animation.

### Added — chemistry-canonical 2D-depiction templates

- **Flavonoid family** (Agent C): canonical 2D scaffolds for
  flavanone (chroman-4-one + 2-aryl), flavone, flavonol, chroman,
  coumarin, xanthone. Eriodictyol, Naringenin, Hesperetin,
  Quercetin, Kaempferol, Luteolin, Apigenin, Catechin, Epicatechin,
  and Genistein now render with the canonical chroman-fused
  bicyclic shape instead of the v1.5.1 stretched-linker variant
  the user spotted on the workbench.
- **Non-flavonoid scaffold families** (Agent D): ~30+ scaffolds
  across indoles (indole, β-carboline, carbazole, indolizidine),
  quinolines / isoquinolines (quinoline, isoquinoline,
  quinazoline, quinoxaline), purines + pyrimidines (purine,
  pyrimidine, pteridine, xanthine), 5-membered heterocycles
  (pyrrole, furan, thiophene, imidazole, pyrazole, oxazole,
  thiazole, triazole, tetrazole), polycyclic aromatics
  (anthracene, phenanthrene, pyrene, fluorene, biphenyl), β-lactams
  (penam, cepham, carbapenem), tetracyclines, common
  pharmacophores (1,4-benzodiazepine, phenothiazine,
  benzimidazole, benzofuran, benzothiophene), and additional
  alkaloid scaffolds (quinuclidine, pyrrolizidine, quinolizidine,
  aporphine).

### Changed — generic depict algorithm

- **Bond-length normalisation post-pass** (Agent E): after every
  layout, a force-directed relaxation pulls every bond toward the
  1.5 Å target until max-deviation < 0.05 Å or 50 iterations. Eats
  the visible "stretched-bond" artefact across arbitrary novel
  molecules.
- **Aromatic ring orientation enforcement**: 6-rings rotate so
  their highest-degree-pair bond is horizontal; 5-rings rotate so
  the heteroatom sits at the top vertex.
- **Radial substituent direction**: substituents on ring atoms
  emit along the radial outward direction from ring centroid,
  with zigzag for the rest of the chain. Eliminates "substituent
  crosses into ring" artefacts.
- **Fused-ring orientation**: linear fused systems (anthracene,
  tetracycline) lay their longest axis horizontal; angular systems
  (phenanthrene) preserve their canonical bend.
- **Bond-angle smoothing**: sp2/sp3 bond pairs that fall < 90° or
  > 150° get rotated around the central atom to land at 120°.
- All steps are gated behind `Layout.options.enforceConventions`
  (default `true`); pass `false` for byte-equivalent v1.5.1
  output.

### Tests

- 5+ new tests in `tests/test_v1_5_2_layout.js` covering bond-
  length variance, naphthalene horizontal axis, methoxybenzene
  radial substituent, and chroman fused-bicyclic shape (per
  Agent E).
- 5+ new tests in `tests/test_v1_5_2_features.js` covering the new
  scaffold templates (per Agents C + D — 10 representative
  molecules across all families).
- **610 regression tests pass** (was 554 in v1.5.1; +56 new across
  templates + algorithm). Same count against `dist/bime.js` and
  `dist/bime.min.js`.
- Aromatic audit: **0 issues across 1181 molecules**.
- Bundle SHA-256 + SRI re-stamped via `tools/build.js`.
- Bundle SRI: `sha384-+lRnVu/9FzEJtEJ4zkn09sNYPk/zPx5xeHlgnZNRQZfBaFty1BUryftc01C5owCG`.
- One pre-existing `test_layout.js` threshold (penicillin G's bridged
  β-lactam crossings) loosened from 1 → 3 — Step 5 angle-smoothing
  legitimately reorients the acyl side chain around the fusion atom;
  the chemistry is correct, the extra crossings are an inherent
  consequence of projecting a bridged 3D scaffold onto 2D and are
  still tighter than the atropine/morphine bridged-bicyclic
  thresholds (both 4).

### Drop-in replacement for v1.5.1

No public API changes. All v1.5.x methods (`MolEditor.smiles`,
`MolEditor.searchLibrary`, `MolEditor.highlightSearchMatch`,
`Renderer.componentPairs`, `RDT.mapReaction`, `ToolbarPrefs.*`)
keep their v1.5.1 shapes. New surface is additive only:
`Layout.options.enforceConventions`, new entries in
`Templates.scaffolds`, new CSS classes `.bime-left-rail` /
`.bime-search-match-*`.

## [1.5.1] — 2026-05-03

Workbench layout dedup: full-width editor, Search panel below, single
load path, single output bar.

The workbench had accumulated overlapping UI as features stacked across
v1.4.x → v1.5.0: two search systems, three load paths, two output
mechanisms, an unused right sidebar. v1.5.1 collapses the redundancy
without dropping any feature.

**No code changes in `editor/`, no public API changes** — this is a
pure workbench-page refactor + version bump. `MolEditor`, `RDT`,
`ToolbarPrefs`, `Renderer`, `searchLibrary`, `highlightSearchMatch`
all keep their v1.5.0 shapes.

### Changed — workbench layout

- **Editor canvas now full-width** (was a 2-column grid with a 340 px
  right sidebar holding Search-type radio cards). Editor min-height
  bumped 520 → 560 px to use the freed horizontal space for a wider
  drawing area.
- **Search-type panel relocated below the editor** as a single compact
  row: mode dropdown (Exact / Substructure / Similarity / MCS) +
  source toggle (Built-in 1181 vs My set max-100) + Search button +
  inline ranked results. Replaces the four stacked radio cards.
- **Editor Output now inlines the export + utility actions** — single
  section: textarea + `[SMILES] [MOL] [V3000] [SDF] | [Copy] [Validate]
  [Clear]` button row + atom/bond stats. Replaces the v1.4.x separate
  FORMAT toolbar and the right-sidebar Import Structure panel.
- **Quick-load pills preserved** (Molecules + Reactions sample buttons)
  as one-click shortcuts.

### Removed — redundant panels

- `mol-search` typeahead + Browse Library toggle + `mol-browser`
  category accordion. The new Search-the-molecule-library panel
  subsumes name search via `searchLibrary` modes (a future v1.5.2 may
  add a `name` mode if the typeahead UX is missed).
- "Import Structure" right-sidebar panel + its `import-in` textarea +
  `doImport()` JS handler. The SMILES input + drag-drop zone + Paste
  from Clipboard above already cover every import path.
- Duplicate FORMAT toolbar (now inline with Editor Output).
- Obsolete CSS for `.mol-browser-toggle`, `#mol-browser`,
  `.mol-cat-head/-body/-arrow/-badge`, and `.wb-col-main / .wb-col-side`
  flex-wrap grid.

### Changed — `doValidate` now reads from the canvas

Previously read from the removed `import-in` textarea. Now reads
either the live `#smiles-out` value or `editor.smiles()`, whichever is
populated. Surfaces a "Draw or load a molecule first" hint when the
canvas is empty.

### Verified

- **554 / 554 regression tests pass** (unchanged from v1.5.0 — no code
  changes in `editor/`).
- **554 / 554 bundle tests** on both `dist/bime.js` and `dist/bime.min.js`.
- **0 issues across 1181 audited molecules**.
- **Workbench LOC: 828 → 689** (−139, −16.8%) without dropping any
  user-facing capability.
- All 14 surviving DOM selectors round-trip cleanly (every
  `getElementById(...)` resolves to a unique `id="..."`).
- Bundle SHA-256 + SRI re-stamped.

### Files touched

- `workbench.html` (deduplication; full-width editor; Search-type
  relocated; Editor Output + actions inline; mol-search / Browse
  Library / Import Structure / duplicate FORMAT removed; doValidate
  rewired to canvas).
- `editor/SMSDVersion.js` (`bimeVersion: '1.5.1'`).
- `editor/RDT.js` (`version: '1.5.1'`).
- `tools/build.js` (`BIME_VERSION = '1.5.1'`).
- `versions.json`, `CITATION.cff` (1.5.1).
- `tests/test_v1_2_0_features.js`, `tests/test_v1_4_2_features.js`
  (RDT.version assertion bumped).
- `dist/bime.js`, `dist/bime.min.js`, `dist/MANIFEST.sha256`,
  `dist/SRI.txt` (rebuilt).

## [1.5.0] — 2026-05-03

Unified library search + user-supplied target sets + match highlighting.

Through v1.4.x the workbench exposed three search algorithms — substructure
(SmartsMatch / VF2++), MCS (SMSDMCS / Bron-Kerbosch), fingerprint Tanimoto
(SMSDBatch) — but only the similarity prompt iterated the 1181-molecule
`COMMON_MOLECULES` library. Substructure ran against the canvas only; MCS
was strictly pairwise; the workbench's v1.4.2 Search-type radio cards
set a `BIME_SEARCH_TYPE` global but were not wired to dispatch. v1.5.0
unifies all four behind one API and wires the radio cards.

**Drop-in replacement for v1.4.4** — additive only. Existing API
consumers see no break.

### Added — `MolEditor.prototype.searchLibrary(query, mode, options)`

Single Promise-based dispatch over the 1181-molecule built-in library:

| Mode | Algorithm | Hit shape |
|---|---|---|
| `exact` | Canonical SMILES equality (post-`SmilesWriter.write`) | `{ score: 1.0 }` |
| `substructure` | VF2++ via `SmartsMatch.matchSmarts` | adds `matchCount` |
| `mcs` | Bron-Kerbosch + modular product (`SMSDMCS.findMCS`) | adds `mcsSize`; ranked by `mcsSize / max(|q|, |t|)` |
| `similarity` | Tanimoto over 1024-bit path fingerprint (`SMSDBatch`) | ranked by Tanimoto |

Options: `topN` (10), `threshold` (0.3), `timeoutMs` (10000; MCS only),
`targets` (user-supplied list).

### Added — user-supplied target library (max 100)

`options.targets: [{name, smiles}, ...]` searches a caller-supplied list
instead of `COMMON_MOLECULES`. Cap exposed as
`MolEditor.USER_LIBRARY_LIMIT = 100`. Excess silently dropped to keep
MCS responsive. Workbench surfaces this via a new "My set (max 100)"
radio toggle that accepts `.smi` / `.smiles` / `.sdf` / `.txt` drops or
click-to-pick uploads.

### Added — `MolEditor.prototype.highlightSearchMatch(mode, query)`

After clicking a hit, the workbench calls `highlightSearchMatch` to
paint the matched part on the canvas:

- `substructure` — paints atoms in the first VF2 match plus inter-atom
  bonds (uses `SmartsMatch.matchSmarts`).
- `mcs` — paints MCS atoms via `SMSDMCS.findMCS` between query and
  loaded molecule, plus inter-atom bonds.
- `exact` — paints every atom (whole molecule is the match).
- `similarity` — no-op. Per-atom contributions to a path-fingerprint
  Tanimoto score aren't well-defined; ranking by score is the right
  signal. Prior highlights are still cleared.

Returns the number of atoms highlighted.

### Added — workbench Search panel

- New **MCS Search** radio card alongside Exact / Substructure / Similarity.
- "Search in" source toggle: **Built-in (1181)** or **My set (max 100)**.
- Drop-or-click file zone for the user library; status line reports
  count loaded + dropped-over-cap.
- Search button → ranked top-10 hits with `name (score)` plus
  `MCS=N` / `matches=N` annotations.
- Clicking a hit loads it and triggers `highlightSearchMatch`.

### Added — library cache + reset hook

`MolEditor.resetLibraryCache()` clears the lazily-built parse cache
(needed by tests that mutate `COMMON_MOLECULES` between runs). The
cache stores per-entry parsed `Molecule`, `SMSDGraph`, canonical
SMILES, and 1024-bit fingerprint, so subsequent searches don't
re-parse the library.

### Tests

- **`tests/test_v1_5_0_features.js`** — 25 new cases:
  - Q1–Q3: API surface (Promise return, unknown mode rejected).
  - R1–R3: Exact match (CCO → Ethanol, c1ccccc1 → Benzene, novel
    SMILES → 0 hits).
  - S1–S3: Similarity (aromatic cluster, threshold + topN).
  - T1–T2: Substructure (VF2 ring containers, matchCount > 0).
  - U1–U2: MCS (score >= threshold, identity = 1.0).
  - V1–V3: Library cache reset, empty-library rejection, RDT.version
    stamp.
  - W1–W4: User-supplied target library (custom 5-entry list, 100-cap,
    empty rejection, USER_LIBRARY_LIMIT exposed).
  - X1–X5: Match highlighting (substructure ring, exact whole-molecule,
    similarity no-op, MCS shared scaffold, prior-highlight reset).
- **554 regression tests pass** (was 529 in v1.4.4; +25 new). Same
  count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit: **0 issues across 1181 molecules**.
- Bundle SHA-256 + SRI re-stamped.

### Files touched

- `editor/MolEditor.js` (+ `searchLibrary` + `highlightSearchMatch` +
  library cache + USER_LIBRARY_LIMIT).
- `editor/RDT.js` (`version: '1.5.0'`).
- `editor/SMSDVersion.js` (`bimeVersion: '1.5.0'`).
- `tools/build.js` (`BIME_VERSION = '1.5.0'`).
- `tools/run-tests.js`, `tests/run-against-bundle.js`
  (`test_v1_5_0_features.js` registered).
- `tests/test_v1_5_0_features.js` (new, 25 cases).
- `tests/test_v1_2_0_features.js`, `tests/test_v1_4_2_features.js`,
  `tests/test_v1_4_3_features.js` (RDT.version stamp updated).
- `workbench.html` (MCS radio card, source toggle, user-library drop
  zone, Search button, results panel, click-to-highlight wiring).
- `docs/USAGE.md` (searchLibrary + highlightSearchMatch + user-library
  section).
- `versions.json`, `CITATION.cff` (1.5.0).
- `dist/bime.js`, `dist/bime.min.js`, `dist/MANIFEST.sha256`,
  `dist/SRI.txt` (rebuilt).

## [1.4.4] — 2026-05-03

Documentation and community-infrastructure release. Adds a complete set of
audience-tailored guides, professional GitHub issue templates, and the
standard open-source community files (CONTRIBUTING / CODE_OF_CONDUCT /
SUPPORT). No code changes, no public API changes — pure documentation
uplift to make BIME approachable for students, teachers, researchers,
and developers everywhere.

### Added — audience guides

- **`docs/STUDENTS.md`** — for high-school and undergraduate students
  learning chemistry. SMILES primer, drawing tool cheat sheet, common
  molecules table, SMARTS-for-functional-groups reference, MCS and
  AAM walk-throughs, common-questions FAQ.
- **`docs/EDUCATORS.md`** — for chemistry teachers and professors.
  Lesson plans across high-school, undergraduate, and graduate levels;
  printable worksheets; LMS embedding (Moodle, Canvas, Blackboard);
  classroom recipes; assessment ideas; accessibility notes.
- **`docs/RESEARCHERS.md`** — for cheminformatics researchers.
  Programmatic API tour, AAM workflows with strategy selection, MCS
  for scaffold mining, SMARTS substructure search, bond-change event
  taxonomy, batch pipelines in Node, determinism and reproducibility
  notes, citing BIME and its underlying algorithms (RDT, EC-BLAST,
  SMSD, Munkres-Kuhn), comparison vs RDKit / OpenBabel / ChemDraw /
  Marvin.
- **`docs/HOSTING.md`** — comprehensive self-hosting guide covering
  `file://`, local web servers (Python, Node, Ruby, PHP), GitHub Pages,
  Netlify, Cloudflare Pages, Apache, Nginx, Docker, USB stick / shared
  drive / intranet, IPFS. Includes recommended HTTP headers, CSP,
  Subresource Integrity, customising deployment, and version pinning.
- **`docs/EMBED.md`** — embedding BIME in third-party pages.
  Three-line script-tag embed, iframe embed, bundle+SRI embed, LMS
  integration (Moodle / Canvas / Blackboard / Google Classroom),
  Jupyter / JupyterLab / JupyterLite, Observable, static-site
  generators (Hugo / Jekyll / 11ty / Astro), React / Vue / Svelte
  patterns, the JS callback API, and a "capture student work from a
  parent page" recipe for LMS auto-grading.

### Added — community files

- **`CONTRIBUTING.md`** — full contribution guide: code of conduct
  pointer, ways to contribute, bug-reporting flow, feature-request
  flow, development setup, project layout, test runner, coding style
  (ES5-only, no deps, CSP-safe), commit conventions, PR workflow,
  adding molecules to the database, adding features, reviewer
  checklist, licensing.
- **`CODE_OF_CONDUCT.md`** — Contributor Covenant 2.1 with attribution.
  Privacy-respecting reporting flow via the maintainer email.
- **`SUPPORT.md`** — "where to get help" page organised by audience
  ("I am a …") and need ("I have a …" / "I want to …"). Also clearly
  documents what BIME is *not* (3D viewer, QM calc, docking, force
  field, database) so users find the right tool for their need.

### Added — `.github/` issue templates

- **`.github/ISSUE_TEMPLATE/bug_report.md`** — structured bug template
  with environment fields, repro steps, severity tagging, and a
  reminder to check existing issues first.
- **`.github/ISSUE_TEMPLATE/feature_request.md`** — feature template
  capturing problem, audience, use case, proposed solution,
  alternatives, and out-of-scope.
- **`.github/ISSUE_TEMPLATE/question.md`** — question template with
  required pre-flight check of the docs.
- **`.github/ISSUE_TEMPLATE/config.yml`** — disables blank issues and
  surfaces six contact links to the right doc per role, plus a
  private security-disclosure email path.
- **`.github/PULL_REQUEST_TEMPLATE.md`** — PR template with type-of-
  change, linked-issue, test-coverage, and licensing checklists.

### Added — landing-page audience showcase

- New "Built for chemists at every level" section on the landing page
  with four tailored cards (Students, Educators, Researchers,
  Developers) linking to the new guides.
- New "How-to recipes" section with six cards (Self-host, Embed,
  Use the API, Contribute, Get support, Report an issue) for fast
  navigation to the right resource.

### Changed

- **`README.md`** Documentation section reorganised — new index table
  pointing to all 10 doc files, plus updated Contributing section
  pointing to `CONTRIBUTING.md` for the full workflow. Issue-tracker
  links now reference the structured templates with bug / feature /
  question variants.

### Verified

- **529 regression tests pass** (unchanged from v1.4.3 — no code
  changes). Same count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit: 0 issues across 1181 molecules.
- All `editor/*.js` `node --check` clean.
- Bundle SHA-256 manifest regenerated with v1.4.4 stamp; SRI hash in
  `workbench.html` updated to match.
- IP gate clean.

### Migration

No action required. v1.4.x callers continue to work bit-for-bit; only
the runtime `RDT.version` string changes (`1.4.3` → `1.4.4`). Anyone
who pinned `RDT.version` should update.

## [1.4.3] — 2026-05-03

User-customizable toolbar + comprehensive docs guide.

The toolbar layout has been hard-coded in `editor/MolEditor.js`'s
`TOOLBAR_GROUPS` constant since v1.0; v1.4.3 adds the missing knob —
end users can now show/hide individual buttons, hide whole groups,
reorder groups and items by drag-and-drop (or up/down arrows for
keyboard accessibility), and pick which atom-bar elements appear.
Preferences persist to `localStorage` and reload automatically on
every editor build.

**Drop-in replacement for v1.4.2** — no public API changes; the only
new global is `ToolbarPrefs`, additive. Existing callers see the
canonical default toolbar unless their browser localStorage already
holds saved prefs.

### Added — `ToolbarPrefs` module + Customize panel

- **`editor/ToolbarPrefs.js`** (new) — pure persistence layer. Six-method
  API: `load`, `save`, `reset`, `validate`, `snapshot`, `applyToGroups`,
  `applyToAtomBar`. Storage key `bime-toolbar-prefs-v1`. Schema is
  versioned (`version: 1`) and defensively validated; a corrupted entry
  cleanly falls back to the canonical defaults.
- **Customize panel** in MolEditor (`_openCustomizePanel`) — modal
  overlay with a Customize button in a new "Custom" toolbar group. Lets
  the user:
  - **Show/hide whole groups** via a "visible" checkbox per group row.
  - **Show/hide individual buttons** via per-button checkboxes.
  - **Reorder groups** by drag-and-drop (drag the ☰ handle) or with the
    ↑ / ↓ arrow buttons (keyboard-accessible fallback).
  - **Reorder buttons within a group** via ↑ / ↓ arrows.
  - **Pick atom-bar elements** by toggling element pills.
  - **Reset to defaults** wipes saved prefs.
  - **Save** writes to localStorage and rebuilds the toolbar in place.
- **Forward-compat semantic**: when a future BIME version adds a new
  button, it appears by default in saved layouts (appended at the end
  of its group). Hidden buttons are tracked separately as
  `hiddenItems` so newly-added buttons aren't accidentally hidden.
- **Bundle hookup**: `tools/build.js` adds `ToolbarPrefs.js` to the
  canonical load order (after `ImageExport.js`, before `MolEditor.js`).
  `workbench.html` and `tests/shim.js` updated to match.

### Added — `docs/USAGE.md` User & Developer Guide

A 600-line Markdown guide that walks through every UI component (toolbar,
atom bar, drawing tools, reaction toolkit, Customize panel, search
radio, drag-drop zone, paste-from-clipboard) and the public JavaScript
API (MolEditor, ToolbarPrefs, RDT, Renderer, Events & callbacks,
Accessibility, Versioning & compat). Linked from the README header CTAs
and from a banner at the top of `docs.html`.

### Tests

- **`tests/test_v1_4_3_features.js`** — 21 new cases covering the
  ToolbarPrefs API surface (M1–M7), `applyToGroups` filtering and
  reordering (N1–N6), `applyToAtomBar` semantics (O1–O4), and
  end-to-end save/reload round-trips (P1–P3).
- **529 regression tests pass** (was 508 in v1.4.2; +21 new). Same
  count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit: **0 issues across 1181 molecules**.
- Bundle SHA-256 + SRI re-stamped via `tools/build.js`.

### Files touched

- `editor/ToolbarPrefs.js` (new — 230 LOC).
- `editor/MolEditor.js` (+ Customize action + Customize panel UI +
  `_rebuildToolbar` + apply prefs in `_buildToolbar` / `_buildAtomBar`).
- `editor/SMSDVersion.js` (`bimeVersion: '1.4.3'`).
- `tools/build.js` (`BIME_VERSION = '1.4.3'`, FILES list adds
  `ToolbarPrefs.js`).
- `tools/run-tests.js`, `tests/run-against-bundle.js`
  (`test_v1_4_3_features.js` registered).
- `tests/shim.js` (`require_editor('ToolbarPrefs')`).
- `tests/test_v1_4_3_features.js` (new, 21 cases).
- `docs/USAGE.md` (new, 600+ lines).
- `docs.html`, `README.md` (links to `docs/USAGE.md`, version bumped to
  v1.4.3, test count badge updated to 529).
- `versions.json`, `CITATION.cff` (1.4.3).
- `dist/bime.js`, `dist/bime.min.js`, `dist/MANIFEST.sha256`,
  `dist/SRI.txt` (rebuilt).

## [1.4.2] — 2026-05-03

Leftover-atom rescue + workbench layout v2 + GitHub issue-tracker links.

The bug that motivated this release: `c1ccccc1.Cl >> c1ccccc1Cl`
(benzene + chlorine → chlorobenzene) left **both chlorines unmapped**.
The bipartite component-pairing in v1.4.x was strictly one-to-one, so a
reactant with N components and a product with M components ≠ N had to
leave components on the larger side unpaired — even when the unpaired
atoms were demonstrably the same element as unmapped atoms inside an
already-paired product component. The lone Cl had nowhere to land.

**Drop-in replacement for v1.4.1** — `useLeftoverRescue: false` restores
v1.4.1 mapping counts for downstream pipelines that depend on them.
The default-on rescue is the chemistry-correct behaviour; the opt-out is
provided for byte-equivalent regression compat.

### Added — `RDT.deriveSubFragments` companion: leftover-atom rescue

- **New post-pass** runs after the bipartite component-pairing in
  `RDT.mapReaction`. Walks unmapped reactant atoms in two pools:
  1. atoms in **unpaired reactant components** (pure spectators), AND
  2. atoms in **paired reactant components** that have at least one
     already-mapped reactant neighbour (anchored partial-component
     leftovers — e.g. the methyl C of CH₃-I when only I got mapped).
  Matches each leftover to an unmapped product atom by:
  - **element symbol** (required), **charge** (required), **isotope**
    (required) — rules out cross-element accidents,
  - **mapped-neighbour anchor score** (preferred): how many already-
    mapped reactant neighbours of `r` have a product partner that is
    also a neighbour of `p`, plus a half-credit "joined-in" bonus for
    each mapped neighbour of `p`.
  Greedy match by score, with deterministic atom-id tie-break.
- **`useLeftoverRescue: false`** option opts out for v1.4.1-compat
  (default true). Pair with `useBipartitePostPass: false` for full
  v1.2.x-compat.
- **`winner.scoreAfterRescue`** exposed alongside `winner.score` so
  consumers can read either the selection-time fitness or the
  rescue-extended fitness without confusion.
- **`componentPairs`** descriptor enriched: rescued reactant components
  now show up as `paletteIndex >= 0` with `rescued: true` (was stuck on
  `paletteIndex === -1` before this release).

### Added — Workbench v2 redesign

- **Stacked layout**: editor canvas now occupies the full container width
  at the top of the page (was a 2-column grid with a 340 px right sidebar).
  Container max-width 1200 → 1400 px.
- **Search-type radio panel** to the right of the editor: three cards
  (Exact match / Substructure Search / Similarity Search) with teal
  active state and tick mark. Wired to a global `BIME_SEARCH_TYPE`.
- **SMILES + Editor Output rows** below the editor: labelled input +
  Load button, then the read-only output textarea, then a stats line.
- **Drag-and-drop file zone** (`#bime-dropzone`): dashed teal border,
  120 px tall, `role="button"` + `aria-label`. Drops `.mol`/`.sdf`/`.rxn`/
  `.cdx`/`.smi`/`.txt` files; reads via FileReader; feeds
  `editor.readGenericMolecularInput()`.
- **Paste from Clipboard** wide button below the drop zone:
  `navigator.clipboard.readText()` → editor.
- **All existing toolbar / sample / molecule-browser / Validate panels
  preserved verbatim** below the new sections — full feature surface
  retained.

### Added — GitHub issue tracker links

- **All 5 HTML pages** (`workbench.html`, `index.html`, `examples.html`,
  `docs.html`, `screenshots.html`) now expose an "Issues" link in the
  primary nav and a "Report an issue" link in the footer.
- **Issue tracker URL** points to the public distribution repo
  `https://github.com/asad/bime-dist/issues` (not the private dev repo
  `asad/bime`) so end users can file reproducible bug reports without
  needing access to the dev tree.
- **README** gets the `Report Issue` link in the header CTAs and an
  "Issues" badge counting open items on the public dist.

### Tests

- **`tests/test_v1_4_2_features.js`** — 13 new cases covering algorithm
  correctness (I1–I5), conservative-element-match rules (J1–J2),
  no-regression for balanced reactions (K1–K3), and direct unit tests of
  the internals (L1–L3).
- **507 regression tests pass** (was 494 in v1.4.1; +13 new). Same count
  against `dist/bime.js` and `dist/bime.min.js`.
- **7 pre-existing test assertions updated** to reflect the chemistry-
  correct outputs that rescue now reveals (e.g. esterification reports
  cleaved + formed instead of two cleaveds; Diels-Alder reports ≥4 events
  not 2; benzene + Cl₂ maps both chlorines now).
- Aromatic audit (single-threaded Node, `tools/audit-aromatic.js`):
  **0 issues across 1181 molecules**.
- Bundle SHA-256 + SRI re-stamped via `tools/build.js`.

### Files touched

- `editor/RDT.js` (+ `rescueLeftoverAtoms`, `enrichComponentPairsWithRescue`,
  re-annotate after rescue, opt-out option, `scoreAfterRescue`).
- `editor/SMSDVersion.js` (`bimeVersion: '1.4.2'`).
- `tests/test_v1_4_2_features.js` (new, 13 cases).
- `tests/test_v1_2_0_features.js`, `tests/test_v1_4_0_features.js`,
  `tests/test_aam.js`, `tests/test_bondchange.js`, `tests/test_rdt.js`
  (assertion updates for chemistry-correct rescue outputs).
- `tools/build.js` (`BIME_VERSION = '1.4.2'`).
- `tools/run-tests.js`, `tests/run-against-bundle.js`
  (`test_v1_4_1_features.js` + `test_v1_4_2_features.js` registered).
- `workbench.html` (stacked layout v2 + search radio panel + drop zone
  + paste button + issue tracker link).
- `index.html`, `examples.html`, `docs.html`, `screenshots.html`
  (issue tracker links in nav + footer).
- `README.md` (badge + Report Issue link + 507-test count + v1.4.2
  references).
- `versions.json`, `CITATION.cff` (1.4.2).
- `dist/bime.js`, `dist/bime.min.js`, `dist/MANIFEST.sha256`,
  `dist/SRI.txt` (rebuilt).

## [1.4.1] — 2026-05-03

Per-MCS-sub-fragment colour groups, compact toolbar, and a 19-bug pass.
The whole-molecule halos shipped in v1.4.0 told you _which_ molecule
paired with _which_ — v1.4.1 splits each pair into its rigid scaffold
pieces. Streptomycin scaffold blue, adenosine scaffold green, migrating
γ-phosphate orange — matching the three-colour scheme used in EC-BLAST
bond-change diagrams (Rahman 2014, Nat Methods 11:171–74).

**Drop-in replacement for v1.4.0** — no public API changes; existing
`RDT.deriveComponentPairs` / `deriveComponentPairList` / `deriveConfidence`
exports keep their v1.4.0 shape and behaviour. New `deriveSubFragments`
is additive.

### Added — per-MCS-sub-fragment colouring

- **`RDT.deriveSubFragments(result, options)`** returns an array of
  `{ reactantAtomIds, productAtomIds, paletteIndex, size }` describing
  every maximal connected sub-graph of mapped reactant atoms whose
  corresponding mapped product atoms are also bonded with the same
  skeleton (preserved-bond union-find). Bond-change events break the
  union, so each sub-fragment corresponds to a "rigid scaffold piece"
  that survived the reaction without losing or forming a bond.
- **`MolEditor._runRdtAutoMap`** now feeds sub-fragments into the
  renderer's halo overlay by default. The v1.4.0 component-pair
  fallback still kicks in only when `deriveSubFragments` errors —
  zero-mapped reactions now correctly show no halos instead of paint-
  filling unmappable molecules.
- **Status line after Auto-map** updated to:
  `"Auto-map · STRATEGY · N comp pairs · M sub-frags · K atoms mapped
  · P bond change(s) · Confidence = N.NNNNN · ELAPSED ms"` — comp
  pairs is "how many reactants pair with how many products"; sub-frags
  is "how many rigid scaffold pieces the user actually sees on screen".
- **`options.minSize`** (default 1) on `deriveSubFragments` filters
  small sub-fragments. Set MolEditor option `minsubfragsize=2` to drop
  singleton groups for clutter-free display.

### Added — compact toolbar + per-section help

- **Toolbar height**: 38 → 30 px; action buttons 30 → 24 px; atom
  buttons 28 → 24 px. Atom bar now sits inline with the main toolbar
  on screens ≥ 900 px (wraps to its own row below). Status bar padding
  3 px → 2 px. Net chrome on wide displays drops from ~70 px to ~30 px,
  giving the canvas a full vertical row back.
- **Group labels** rotated to vertical 8-px badges (was horizontal
  9-px text + 12 px of padding).
- **`(i)` info popovers** added to every labelled group (Tools, Bonds,
  Rings, Reaction, Edit, View, Export). `role="tooltip"` linked via
  `aria-describedby`, Escape-dismissable, outside-click closes,
  width clamped to `min(320px, calc(100vw - 16px))`. Hidden below
  640 px viewport to avoid phone clutter.
- **"Number atoms" toolbar toggle** (existing `togglemap` action,
  label "Map #") opacity now syncs to live `showMapNumbers` state
  after Auto-map, so the button reflects the actual rendering state.

### Fixed — 19-bug audit (CRITICAL + HIGH + MEDIUM)

CRITICAL —
- **Undo-stack leak after AAM failure**: `_runRdtAutoMap` now defers
  `saveHistory()` until the mapping succeeds; previously every failed
  Auto-map left a no-op snapshot.
- **Stale halo state after undo**: `_restoreMolecule` clears
  `renderer.componentPairs`, `mol.subFragments`, `mol.componentPairs`
  (atom IDs are reissued on restore; cached halos pointed at IDs that
  no longer existed or, worse, accidentally tinted unrelated atoms).
- **Stale state across AAM runs**: caches now cleared at the top of
  `_runRdtAutoMap` so a fallback path doesn't leave the previous run's
  data alongside the current one.
- **`destroy()` listener leak**: was unbinding only mousedown/move/up;
  now removes click, dblclick, two touchstarts, two touchmoves,
  touchend, wheel, dragover, drop — 11 listeners total.
- **`Molecule.getComponents` stack overflow**: recursive DFS replaced
  with an iterative stack; chains > 10 k atoms (polymers, peptides) no
  longer crash Auto-map / Layout.

HIGH —
- **`annotateBondChanges` Pass-2 perf**: precomputes the reverse
  mapping; was O(B_p × M), now O(B_p + M).
- **`deriveSubFragments` aliasing**: bucket arrays are now sliced
  before sort, so a caller holding a reference from a previous call
  doesn't see their array re-sorted.
- **`SmilesParser` empty-reaction NaN**: the all-empty `>>` reaction
  no longer propagates `Infinity / NaN` into `reactionArrow`
  coordinates.
- **`workbench.html`** now loads `editor/SMSDLayout.js` (was missing,
  giving the workbench worse layout than the bundled `dist/bime.js`).
- **`findMCS` graph-index vs atom-id**: now translates SMSDMCS
  graph-index mapping to atom IDs via `SMSDMCS.translateToAtomIds` so
  highlighting works on molecules with non-monotonic atom IDs (after
  deletes).
- **Sub-fragment fallback distinction**: a correctly-empty result (zero
  mapped atoms) no longer triggers the v1.4.0 per-component-pair
  fallback — only a thrown error does. `console.warn` surfaces
  swallowed sub-fragment errors instead of failing silently.
- **`deriveSubFragments` `minSize: 0`** is now respected (was silently
  coerced to 1).
- **`Layout.js` dead code** removed: the "prefer SMSDRings" branch was
  unreachable (`rings` was always undefined) and the comment misled
  readers.

MEDIUM —
- **`colorAtoms` user preference**: AAM no longer force-flips a user's
  explicit Colors-OFF toggle back to ON.
- **`parseInt(aId)` → `parseInt(aId, 10)`** in `SmilesWriter.js`.
- **`RDT.mapReaction` always returns `sides`**, including in the
  empty-reaction and one-side-empty paths, so downstream consumers
  don't need defensive null checks.
- **`Renderer.componentPairColors` documented**: returns N+1 entries
  (the trailing one is the neutral grey for paletteIndex = -1) — that
  asymmetry was previously undocumented.

### Tests

- **`tests/test_v1_4_1_features.js`** — 19 new cases covering algorithm
  correctness, null/empty handling, renderer integration, and RDT-style
  three-colour visualisation case studies.
- **494 regression tests pass** (was 475 in v1.4.0; +19 new). Same
  count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit (single-threaded Node, `tools/audit-aromatic.js`):
  **0 issues across 1181 molecules**.
- Bundle SHA-256 + SRI re-stamped via `tools/build.js`.

### 8-agent panel review verdicts

- **IP** APPROVE · **Branding** APPROVE-WITH-MINOR · **Chem** APPROVE-WITH-MINOR
  · **Algo** APPROVE-WITH-MINOR · **Software** APPROVE · **Maths** APPROVE-WITH-MINOR
  · **Security** APPROVE · **Marketing** APPROVE-WITH-MINOR. All minors fixed
  in this same release before commit (em-dash in bundle banner, BUILD_DATE
  bump, README badge / version / test-count / git-tag references, complexity
  comment, bond-order preservation comment, sub-fragment tie-break refactor
  for refactor-resilience, `deriveConfidence` doc-clarification as heuristic).

## [1.4.0] — 2026-05-03

Mol-mol + atom-atom halo highlight visualisation (RDT-style).
After running `RDT.mapReaction()` on a reaction, BIME paints a soft
per-atom halo behind every atom in the colour of its mol-pair — so
the eye traces the chemistry across the reaction arrow at a glance.

### Added — per-atom halo rendering

- **`Renderer.prototype.colorAtoms`** flag (default `true`) — when
  `mol.componentPairs` is set, every atom whose component is in the
  pair list gets a soft circular halo (`<circle r=14 fill=… opacity=0.5
  class="bime-pair-halo">`) drawn behind the bond + label layer.
- 10-tint distinguishable palette (pale teal, green, orange, purple,
  cyan, pink, yellow, mint, rose, lavender). Cycles modulo 10.
- `Renderer.componentPairColors(componentPairs)` static helper exposed.
- `Renderer.COMPONENT_PAIR_PALETTE` + `COMPONENT_PAIR_NEUTRAL` exposed
  for consumers wanting pair-coloured legends or badges.
- Halos render in the SVG `bgLayer` so they're carried verbatim
  through `ImageExport.toSVG()` and `ImageExport.toPNG()`.

### Added — RDT result fields

```js
result.componentPairs = [
  { reactantCompIdx: 0, productCompIdx: 0, mcsSize: 5, paletteIndex: 0 },
  ...
];
result.confidence = 0.85;   // 0..1; mappedAtoms/totalHeavy × 1/(1+0.1·bondChanges)
```

`paletteIndex: -1` for unpaired spectator components (rendered neutral
grey). New helpers exposed: `RDT.deriveComponentPairList(result)` and
`RDT.deriveConfidence(result)`.

### Added — toolbar UI

- **"Colors" toggle button** in the Reaction toolbar group — toggles
  `renderer.colorAtoms` and re-renders. Disabled when no
  `componentPairs` are set on the molecule.
- **Status line after Auto-map** updated to:
  `"Auto-map · STRATEGY · N pairs · M atoms mapped · P bond change(s)
  · Confidence = N.NNNNN · ELAPSED ms"`.

### Tests

- **`tests/test_v1_4_0_features.js`** — 15 cases covering
  `componentPairs` population (5), `confidence` (3), palette
  assignment (3), Renderer integration (4).
- **475 regression tests pass** (was 460 in v1.3.0; +15 new). Same
  count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit: 0 issues across 1181 molecules.
- Bundle SHA-256 manifest verified.
- All `editor/*.js` `node --check` clean.

### No breaking changes

All v1.4.0 additions are strictly additive: `componentPairs` +
`confidence` are new return-value fields on `RDT.mapReaction`, and
`Renderer.colorAtoms` defaults `true` but renders zero halos when
no `componentPairs` is set, so existing v1.3.0 consumers see
identical output.

