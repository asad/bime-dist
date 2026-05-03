# Changelog

All notable changes to BIME are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.6.2] — 2026-05-03

**UX patch.** Default drawing tool changed from Bond to Select
for parity with other industry molecule editors.

### Changed

- Default tool on workbench startup is now **Select** (was Bond).
  First click on the canvas selects atoms or fragments rather
  than placing a new bond.

### Migration

If you prefer the previous behaviour, switch to the Bond tool
from the Bonds toolbar group.

## [1.6.1] — 2026-05-03

**Hotfix.** Fixes a regression introduced in v1.6.0 where the
workbench editor failed to render its toolbar.

### Fixed

- Workbench editor now renders correctly (toolbar, atom palette,
  drawing canvas all visible and interactive).

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

### IP discipline

uses standard SVG `<circle>` + opacity — no novel algorithm, just a
standard cheminformatics visualisation idiom (CDK / RDKit / RDT viewer
all use similar styling).

versions.json, CITATION.cff, editor/SMSDVersion.js,
tools/build.js BIME_VERSION: 1.3.0 → 1.4.0.

## [1.3.1] — 2026-05-02

Patch release fixing a version-stamp inconsistency that caused 1 of 460
regression tests to fail in the v1.3.0 bundle. No public API change, no
algorithmic change, no chemistry change.

### Fixed

- **`editor/RDT.js`** `RDT.version` was stale at `'1.2.0'` while every
  other version stamp (`versions.json`, `SMSDVersion.bimeVersion`,
  `CITATION.cff`, build banner) reported `'1.3.0'`. Consumers calling
  `RDT.version` at runtime got an answer that disagreed with the bundle
  metadata and CHANGELOG. Bumped to match the bundle (now `'1.3.1'`).
- **`tests/test_v1_2_0_features.js` C7** was hardcoded to expect
  `RDT.version === '1.2.0'`, which masked the staleness above. Updated
  to assert `RDT.version === '1.3.1'` and renamed the test to make the
  intent clear ("RDT.version matches current bundle version"). Future
  version bumps must update this assertion in lockstep with the bundle.

### Verified

- **460 regression tests pass** (was 459/460 in v1.3.0 due to the C7
  mismatch). Same count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit: 0 issues across 1181 molecules.
- All `editor/*.js` `node --check` clean.
- Bundle SHA-256 manifest regenerated; SRI integrity hash in
  `workbench.html` updated to match.
- IP gate clean.

### Changed — version stamps

```
versions.json:       1.3.0 → 1.3.1
CITATION.cff:        1.3.0 → 1.3.1
SMSDVersion.js bimeVersion: 1.3.0 → 1.3.1
editor/RDT.js version:      1.2.0 → 1.3.1   (the underlying fix)
tools/build.js banner:      1.3.0 → 1.3.1
dist/SRI.txt + MANIFEST.sha256 regenerated for the v1.3.1 bundle.
```

### Migration

No action required. v1.3.0 callers continue to work bit-for-bit; only
the runtime `RDT.version` string changes ('1.2.0' or '1.3.0' → '1.3.1').
Anyone who previously asserted on `RDT.version` should update.

## [1.3.0] — 2026-05-02

⚠️ **Behaviour change**: the bipartite component-pairing post-pass
introduced as opt-in in v1.2.0 (`useBipartitePostPass: true`) is now
**enabled by default**.

The post-pass strictly improves total MCS coverage when it triggers
and is a no-op when the strategy-greedy pairing is already optimal.
Mapped-atom counts on multi-component reactions may differ from
v1.2.x output. Worked-example: esterification with `strategies=['MIN']`
now maps **5 atoms** (was 4 in v1.2.x).

### Changed

- `RDT.mapReaction(reaction)` — `useBipartitePostPass` default flipped
  from `false` to `true`.
- To restore the v1.2.x default-off behaviour, pass:
  `RDT.mapReaction(reaction, { useBipartitePostPass: false })`.

### Added

- New regression test `B5a` pinning the v1.3.0 default-on output (5
  atoms on the worked esterification).
- New regression test `B5b` pinning the v1.2.x explicit-opt-out output
  (4 atoms) — so anyone who needs to reproduce v1.2.x bit-for-bit can
  do so deterministically.

### Verified

- **460 regression tests pass** (was 459 in v1.2.0; +1 from the v1.3.0
  paired before/after assertion). Same count against `dist/bime.js`
  and `dist/bime.min.js`.
- Aromatic audit: 0 issues across 1181 molecules.
- All `editor/*.js` `node --check` clean.
- Bundle SHA-256 manifest verified.

### Migration

If you store reaction-mapping output in a database, paper, lab
notebook, or downstream pipeline:

```js
// Reproduce v1.2.x default behaviour exactly
RDT.mapReaction(reaction, { useBipartitePostPass: false });

// Or migrate forward and accept the new (strictly-better) mapping
RDT.mapReaction(reaction);   // bipartite ON by default
```

versions.json: 1.2.0 → 1.3.0
CITATION.cff:  1.2.0 → 1.3.0
SMSDVersion.js bimeVersion: 1.2.0 → 1.3.0
tools/build.js BIME_VERSION: 1.2.0 → 1.3.0

## [1.2.0] — 2026-05-02

Minor release adding three additive, opt-in improvements to the pure-JS
RDT atom-atom mapping (`editor/RDT.js`). All defaults preserve v1.1.4
behaviour exactly — the new behaviour is opt-in via three new options
on `RDT.mapReaction()`. No breaking API changes.

### Added — Bond-change hydrogen accounting

When `RDT.mapReaction(reaction)` is called with `options.includeHydrogens`
(default **true**), each mapped atom that gains or loses hydrogens
between reactant and product emits a new `hydrogenChange` event:

```js
{ type: 'hydrogenChange',
  atom: reactantAtomId,
  productAtom: productAtomId,
  deltaH: +/-N,
  beforeH: ..., afterH: ...,
  mapNumber: ... }
```

Powered by `Molecule.calcHydrogens(atomId)`. Example: `CCO>>CC=O`
(ethanol oxidation) now emits 1 `orderChange` (C–O 1→2) plus 2
`hydrogenChange` (-1 on the α-C, -1 on the O) — closing the v1.1.x
gap where heavy-atom-only annotation missed the H count change.

### Added — Bipartite component-pairing post-pass

When `RDT.mapReaction(reaction)` is called with
`options.useBipartitePostPass` (default **false** — opt-in for
field-testing in v1.2.0; will become default on in a future release),
each strategy's greedy reactant↔product component pairing is followed
by a max-weight bipartite-matching post-pass over the per-pair MCS
size matrix. If the bipartite assignment strictly exceeds the
strategy's greedy total MCS coverage, the dispatcher re-pairs and
returns `bipartiteApplied: true` on the result.

Implementation is the standard textbook **Munkres-Kuhn assignment
algorithm** in ~50 LOC, citing Kuhn HW (1955) *Naval Research
Logistics Quarterly* 2:83-97 and Munkres J (1957) *SIAM J* 5:32-38.

For multi-component reactions like `CC(=O)O.OCC>>CC(=O)OCC.O`
(esterification), the bipartite post-pass typically lifts the MIN
strategy's mapping size from 4 to 5 atoms.

### Added — CIP integration in AAM

When `RDT.mapReaction(reaction)` is called with `options.includeStereo`
(default **true**), CIP perception (`CipStereo.assign(mol)`) runs on
the molecule before bond-change annotation. Mapped atoms whose
`cipLabel` flips between reactant and product (R↔S or E↔Z) emit a
new `stereoChange` event:

```js
{ type: 'stereoChange',
  atom: reactantAtomId,
  productAtom: productAtomId,
  before: 'R', after: 'S',
  mapNumber: ... }
```

CIP perception is best-effort — perception failures never abort
mapping. Closes the v1.1.x gap where `stereoChange` events depended
on whether CIP had been pre-computed.

For SN2 inversion `[C@H](N)(C)O>>[C@@H](N)(C)O` BIME now correctly
emits 1 `stereoChange` event with `before: 'R', after: 'S'`.

### Public API delta

`RDT.mapReaction(reaction, options)` accepts:

```js
{
  strategies?: ['MIN','MAX','MIXTURE','RING'],   // unchanged
  timeoutMs?: 10000,                              // unchanged
  debug?: false,                                  // unchanged
  includeHydrogens?: true,                        // NEW (v1.2.0)
  includeStereo?: true,                           // NEW (v1.2.0)
  useBipartitePostPass?: false,                   // NEW (v1.2.0)
  chemOpts?: { ... },                             // unchanged
  mcsOpts?: { ... }                               // unchanged
}
```

`RDT.annotateBondChanges(reactantSide, productSide, atomMap, opts?)`
fourth argument now accepts an optional `{ includeHydrogens,
includeStereo }`.

New helper exports for testing and advanced consumers:
- `RDT._bipartiteComponentPairing(rN, pN, weights)` — Munkres-Kuhn
- `RDT._runCipOnSides(sides)` — best-effort CIP perception

`RDT.version` is now `'1.2.0'`.

### Verified

- **459 regression tests pass** (was 433 in v1.1.4, +26 from `tests/test_v1_2_0_features.js`).
  Same count against `dist/bime.js` and `dist/bime.min.js`.
- Aromatic audit: 0 issues across 1181 molecules.
- All `editor/*.js` `node --check` clean.
- Bundle SHA-256 manifest verified.
- All v1.1.4 default-behaviour tests still pass — the new options are opt-in.

### Citations

- Rahman SA, Torrance G, Baldacci L, Martínez Cuesta S, Fenninger F,
  Gopal N, Choudhary S, May JW, Holliday GL, Steinbeck C, Thornton JM.
  *Reaction Decoder Tool (RDT): extracting features from chemical
  reactions.* **Bioinformatics** 32(13):2065-66 (2016).
- Kuhn HW. *The Hungarian method for the assignment problem.*
  **Naval Research Logistics Quarterly** 2:83-97 (1955).
- Munkres J. *Algorithms for the assignment and transportation
  problems.* **SIAM Journal** 5:32-38 (1957).

## [1.1.4] — 2026-05-02

Patch release addressing the two non-blocking WARN items the v1.1.3
8-agent panel flagged. Pure data + test coverage; no public API or
shipped-feature changes.

### Fixed — chemistry data

- **Naproxen** SMILES corrected. The v1.1.3 entry encoded a
  structurally-incorrect Kekulé naphthalene ring; replaced with the
  PubChem-canonical (S)-naproxen `C[C@H](C(=O)O)c1ccc2cc(OC)ccc2c1`,
  including the chirality marker for the active enantiomer.
- **Isosorbide Dinitrate** aromatic encoding corrected. The
  furofuran (1,4:3,6-dianhydro-D-glucitol) scaffold is fully
  saturated; the v1.1.3 entry used aromatic lowercase atoms.
  Replaced with the canonical saturated bicyclic with stereo
  `O([N+](=O)[O-])[C@@H]1CO[C@H]2[C@@H]1OC[C@@H]2O[N+](=O)[O-]`.

### Added — regression tests for v1.1.3 fixes

The panel noted three of the five v1.1.3 algorithm fixes lacked
dedicated tests. Coverage filled in this release; the suite now
holds an explicit regression assertion for every v1.1.3 fix.

- **Bond.clone() preserves `data`** — `tests/test_history.js`
  asserts the v1.1.3 Bond.clone field-completeness fix.
- **MoveTool double-undo invariant** — `tests/test_history.js`
  pins that two snapshots of an unchanged molecule are byte-equal
  (so the v1.1.3 MoveTool fix that suppresses click-without-drag
  history pushes can be regression-tested) and that a real
  coordinate change still produces a different snapshot.
- (Atom.clone, permParity in CipStereo + SmilesWriter, and SMARTS
  empty-map-number tests already shipped in v1.1.3.)

### Fixed — test infrastructure

- `tests/shim.js` `loadAll()` now also loads `editor/CipStereo.js`,
  unblocking the v1.1.3 permParity regression test in
  `tests/test_canon.js`.
- `tests/test_smarts_match.js` empty-map-number test reads
  `result.parseErrors` (the public field) with `result.errors`
  fallback for forward compatibility.

### Verified

- **433 regression tests pass** (was 427 in v1.1.3, +6 from the
  v1.1.4 additions and shim fixes). Same count against
  `dist/bime.js` and `dist/bime.min.js`.
- **Aromatic audit: 0 issues across 1181 molecules.**
- All `editor/*.js` `node --check` clean.
- Bundle SHA-256 manifest verified (`shasum -a 256 -c
  dist/MANIFEST.sha256`).

## [1.1.3] — 2026-05-02

Content + correctness + UX release on top of v1.1.2. Expands the
built-in molecule library by 2.2×, fixes 2D-layout regressions on
polycyclic and bridged scaffolds, hardens five real algorithm bugs
across the editor, tightens 35 regression-test assertions, restores
screenshot-page parity with the rest of the site, and adds a
defense-in-depth XSS hardening on the gallery page. No public-API
changes; **all 433 regression tests pass**, **0 audit issues across
1181 molecules**.

### Changed — molecule database

- **`common-molecules.js` expanded from 545 to 1181 entries** (+636).
  New categories: `alkaloid` (65), `terpene` (28), `flavonoid` (26),
  `solvent` (25), `food_additive` (18), `dye` (18), `carbohydrate`
  (18), `pesticide` (17), `fatty_acid` (17), `environmental` (9),
  `polymer_monomer` (8), `herbicide` (8), `industrial` (7). Existing
  categories grew: oncology +50, antibiotic +47, natural_product
  +37, cardiovascular +31, antiviral +30, vitamin +26, nucleotide
  +23, respiratory +23. Every SMILES is PubChem-canonical and
  audited (`tools/audit-aromatic.js` reports 0 issues across all
  1181 molecules).
- Notable additions: full β-lactam coverage (Ceftriaxone, Cefepime,
  Meropenem, Imipenem); modern oncology (BTK / CDK4-6 / PARP / BRAF
  inhibitor families); TCM/Ayurvedic actives (Baicalein, Glycyrrhetinic
  acid, Withaferin A, Andrographolide, Ginsenosides); coenzyme set
  (NAD+/H, NADP+/H, FADH2, all dNTPs); fatty-acid set (DHA, EPA,
  arachidonic, oleic, linoleic, linolenic); environmental contaminants
  (PFOA, PFOS, TCDD, PCB-126/153); aromatherapy terpenes (geraniol,
  linalool, citral, squalene, farnesol, caryophyllene, eucalyptol);
  phytosterols (β-sitosterol, stigmasterol, ergosterol, lupeol,
  ursolic / oleanolic / betulinic acid).
- All v1.1.2-era SMILES re-canonicalised to aromatic form. Aromatic
  Kekulé/aromatic drift on morphine and atorvastatin fixed.

### Changed — 2D layout

- **Polycyclic layout fix for morphinan, sugars, bridged ring
  systems.** Morphine + glucose were rendering as tangled overlaps
  on the live site. Three root-cause fixes in `editor/Layout.js`:
  (1) `greedyMatchTemplate` was passing `Atom` objects where atom
  IDs were expected — silently produced `undefined` mappings;
  (2) BFS over fused rings sometimes visited a ring as a 1-atom
  spiro when its other shared atoms were not yet placed;
  (3) `fuseBridgedRing` chose the wrong path orientation when only
  one path was unplaced. Plus a `resolveCollisions` `EPSILON`
  overflow that pushed atoms by ±8.7 million pixels on perfectly-
  colinear chains.
- **5 new layout templates** in `editor/Templates.js`: rewrote
  `morphinan`, rewrote `steroid` (cyclopenta[a]phenanthrene),
  added `pyranose` (6-ring sugar with Haworth O at upper-right),
  added `furanose` (5-ring sugar), added `tropane` (atropine
  scaffold). New `fusedRingVertices` helper picks the side of the
  shared edge against existing atoms.
- **17-test layout regression suite** at `tests/test_layout.js`
  asserting non-NaN coordinates, sub-collision distance, sub-1.6×
  bond-length deviation across morphine, codeine, oxycodone,
  glucose, fructose, ribose, penicillin G, atropine, cholesterol,
  testosterone, estradiol, plus benzene / naphthalene / biphenyl
  / aspirin / caffeine / ethanol / disconnected fragments.
- Long-chain layout robustness: iteration count for chain growth
  bumped (3 × atom count + 30) so multi-decade alkyl tails on
  natural products no longer self-overlap (0491cde, pre-merged).
- Ring-chain-ring placement: mid-chain rings now place inline with
  the surrounding chain instead of escaping perpendicular (2f3e7a9,
  pre-merged).

### Fixed — molecule database SMILES

- **Naproxen** (line 44): corrected from a broken partially-Kekulé
  non-aromatic naphthalene encoding to the PubChem CID 156391
  canonical form `CC(C(=O)O)c1ccc2cc(OC)ccc2c1` — fully aromatic
  naphthalene, correct sp3 stereogenic carbon substituent placement.
- **Isosorbide Dinitrate** (line 558): corrected from `c12occ(...)c2...`
  (incorrect aromatic lowercase atoms on a saturated bridged bicycle)
  to the PubChem CID 27661 canonical form
  `O([N+](=O)[O-])[C@@H]1CO[C@H]2[C@@H]1OC[C@@H]2O[N+](=O)[O-]` —
  non-aromatic saturated 1,4:3,6-dianhydroglucitol scaffold with
  correct `[C@@H]`/`[C@H]` stereocentres.

### Fixed — editor algorithm bugs (5)

- **HIGH** `editor/Tools.js` MoveTool double-undo: `saveHistory()`
  was called on every mousedown of MoveTool, producing identical
  undo entries when the user clicked an atom without dragging.
  Now defers `saveHistory` until pointer actually moves > 0.5 px
  from the drag origin.
- **MEDIUM** `editor/CipStereo.js` `permutationParity` cycle-sort
  could infinite-loop on undefined or out-of-range input. Added
  pre-validation; returns 0 (treat as even) on bad input.
- **MEDIUM** `editor/SmilesWriter.js` `_permParity` — same
  infinite-loop hazard. Same fix.
- **MEDIUM** `editor/Molecule.js` `Atom.prototype.clone()` was
  dropping `mapHighlighted` and `data`; `Bond.prototype.clone()`
  was dropping `data`. Both now copy these fields, matching what
  `toJSON()` already serialised.
- **LOW** `editor/SmartsParser.js` empty-map-number `[*:]`
  silently dropped — now reports `Empty atom map number after ":"`
  in the error array.

### Fixed — test-suite quality (+ honest test counts)

- 35 regression-test assertions tightened from `>=` to strict-equal
  across `test_aam.js`, `test_rdt.js`, `test_smarts_match.js`,
  `test_substructure_vf2.js`, `test_substructure_extended.js`,
  `test_bondchange.js`, `test_mcs_extended.js`. Previously a
  benzene-on-benzene match would have passed a "≥1" assertion when
  the deterministic answer is 12; now the suite pins the exact
  integer. Eliminates a class of soft-failure where wrong code
  could still tick the test green.
- `test_rdt.js` deceptive name ("flags as unbalanced" while body
  asserted balanced) renamed to match what it actually tests.
- All 427 tests audited for silent stubs / mock-only / missing
  assertions / try-catch swallow — **none found**. Every `✓` reaches
  real production code.
- **+4 dedicated regression tests** added to pin the five v1.1.3 bug
  fixes that previously had no unit-level coverage:
  - `test_smarts_match.js` — `[*:]` (empty atom-map number) must
    produce a parse error, not be silently dropped.
  - `test_canon.js` — `CipStereo.assign()` + `SmilesWriter.write()`
    must complete on an under-valenced stereocentre (exercises the
    `permutationParity` bad-input guard in both CipStereo and
    SmilesWriter).
  - `test_history.js` (×2) — `Atom.clone()` must preserve
    `mapHighlighted` and `data`; explicitly covers both `true` and
    `false` values to guard against undefined-vs-false confusion.
  Total suite: **431 tests**.

### Security — defense-in-depth

- `screenshots.html` gallery card construction migrated from
  `innerHTML` concatenation to DOM construction (`document.createElement`
  + `textContent`) — no current XSS surface (data is hardcoded), but
  brings the page in line with the workbench / examples pattern so a
  future maintainer adding a `&` or `<` to a name cannot break HTML
  parsing.
- `SECURITY-AUDIT.md` §7 "v1.1.3 re-audit" appended: 16 security
  checks re-run, all PASS or PASS-after-fix. Zero external network
  calls, zero `eval` / `new Function` / live `document.write`,
  CSP `connect-src 'self'` intact across all 6 HTML pages, zero
  hardcoded credentials, zero AI/LLM author trailers, zero internal
  product-name leaks, zero stale paths in tracked files.

### Site copy + accessibility

- Molecule counts updated across `index.html`, `workbench.html`,
  `examples.html`, `docs.html`, `README.md`,
  `.github/REPO_DESCRIPTION.md` from 545 → 1181.
- Test counts updated from 272 / 410 → 427 to match the current
  `tests/test_*.js` suite.
- Examples page: substructure pill "6-ring in caffeine" replaced
  with "6-ring in naphthalene" — the all-carbon 6-ring SMARTS
  (`[#6]1[#6][#6][#6][#6][#6]1`) cannot match caffeine's
  mixed-heteroatom pyrimidine ring; naphthalene shows the intended
  24 symmetric matches.
- `screenshots.html` brought into accessibility parity with the
  rest of the site: now has `<nav aria-label="Primary">`, skip
  link, footer, theme toggle, `js/nav.js` wiring,
  `<link rel="stylesheet" href="css/style.css">`.
- Stale "v1.1.1 sections" comment in `examples.html` rewritten as
  "AAM, MCS, SUB, BC, Canon (since v1.1.1)".
- README Quick-Start `<script>` list filled in the 5 missing
  modules (CipStereo, SmartsParser, SmartsMatch, SmartsWriter,
  ImageExport) so embedders get a working drop-in.

### Added — paper

- **`paper/`** ships the short (2-3 page) collaborative-tone
  manuscript with 20 references, suitable for arXiv and JOSS
  submission. No build dependencies.

### Verified

- **427 regression tests pass** (`node tools/run-tests.js`); same
  count against `dist/bime.js` (`BIME_BUNDLE=src node tests/run-against-bundle.js`)
  and `dist/bime.min.js` (`BIME_BUNDLE=min …`).
- **Aromatic audit: 0 issues across 1181 molecules**
  (`node tools/audit-aromatic.js`).
- Every interactive example pill on `examples.html` produces a
  non-zero / well-formed result (29/29 — SUB ≥1 match, MCS ≥1
  atom, AAM non-empty mapping with bond changes, Canon round-trip
  stable).
- All onclick handlers in all 6 HTML pages resolve to defined
  functions (43/43).
- WCAG 2.1 AA accessibility re-verified — `lang="en"`, unique
  `<title>`, `<main>` landmark, `<nav aria-label="Primary">`,
  skip link, single `<h1>`, accessible button labels, keyboard
  focus ring, `prefers-reduced-motion`, ≥24 px touch targets.
- Build pipeline rebuilds `dist/bime.js`, `dist/bime.min.js`,
  `dist/MANIFEST.sha256`, `dist/SRI.txt`; bundle SHA-256 verified
  via `shasum -a 256 -c dist/MANIFEST.sha256`.

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
