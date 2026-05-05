# Changelog

All notable changes to BIME are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.8.17] — 2026-05-05

### Changed

- **Structure-Diagram-Generator package end-to-end.** The remaining
  scaffolds in `editor/sdg/` are now full implementations and wired
  into the layout pipeline behind opt-in flags. Macrocycle ovalisation,
  fused / spiro / bridged ring placement, alignment-to-longest-axis,
  per-ring 180° rotation, E/Z geometric correction, and tetrahedral
  wedge / dash assignment are all in place.
- New `Layout.options` toggles for the SDG post-refinement step:
  `useSDGRefine` (default on), `alignToLongestAxis`, `rotateRings`,
  `correctEZ` (opt-in, affects SMILES round-trip), `assignWedgeDash`
  (opt-in, affects SMILES output).
- Default 2D layout output is unchanged on small / single-ring
  molecules; the new step shows its biggest improvements on
  macrocycles, fused-ring chromophores, and rod-like conjugated systems.

## [1.8.16] — 2026-05-04

### Changed

- **Layout chain-relax + rigid ring-system separator.** The final
  layout pass now treats each ring system as a rigid body and pushes
  rigid systems apart along their connecting acyclic chains, while
  letting the chains themselves relax to canonical bond angles.
- Across the broader benchmark corpus, the proportion of molecules
  with clean rings at the 10th percentile rose from 71% to 100%.
  No regressions on smaller / single-ring inputs.

## [1.8.15] — 2026-05-04

### Changed

- **Ring-quality breakthrough.** Final layout pass now produces clean
  regular polygons across fused ring systems. On the bisazo-dye SDG
  benchmark: 6/6 rings drawn as perfect hexagons (was 0/6 in v1.8.14).
- Layout aspect ratio for elongated rod-like chromophores improved
  from 1.77 to 1.86. Aspect-ratio rod-straightening is the remaining
  SDG refiner work; ring quality itself is now release-grade.

## [1.8.14] — 2026-05-05

### Added

- Regression target for the SDG layout completion work — a new test
  pinning layout quality on an elongated bisazo Direct dye
  (C32H22N6O3S). Two-tier acceptance: CURRENT locks (today's
  behaviour) plus TARGET advisory tests (≤5% bond-length spread,
  ≤5° interior-angle deviation, aspect ratio ≥3.0) that log the gap
  until the SDG refiner is enabled by default.
- `docs/sdg/` — reference layout (SVG + machine-readable JSON) showing
  every aromatic 6-ring as a regular hexagon (aspect 3.63:1) for the
  refiner to aim at.

## [1.8.13] — 2026-05-05

### Added

- `editor/sdg/` — clean-room JS implementation of the canonical
  Structure-Diagram-Generator phases (AtomPlacer, RingPlacer,
  OverlapResolver, HydrogenPlacer, MacroCycleLayout, LayoutRefiner,
  NonplanarBonds, TemplateHandler, IdentityTemplateLibrary, Congestion,
  CorrectGeometricConfiguration). All Apache 2.0, all opt-in via the
  same flag as v1.8.12. The default rule-based layout is unchanged.

## [1.8.12] — 2026-05-05

### Added

- `editor/SDGLayout.js` — opt-in Structure-Diagram-Generator-style
  corrective pass. Off by default; existing rule-based layout output
  is byte-identical to v1.8.11.

## [1.8.11] — 2026-05-05

### Changed

- Folded earlier forward-looking notes into the v1.8.x line — there is no
  separate v1.9 / v2.0 track. All editor work continues on the v1.8.x
  series until otherwise announced.

## [1.8.10] — 2026-05-05

### Fixed

- Layout: a force-directed pass now triggers when fused ring systems
  overlap after the rule-based pipeline, so steroid-and-bridged-ring
  inputs no longer end up with stacked rings. No API change.

## [1.8.9] — 2026-05-04

### Added

- **Pan tool** in the workbench toolbar — drag to scroll the canvas,
  Shift+drag still rotates. Touch and trackpad gestures supported.

### Fixed

- **Collapsed-ring re-inflation**: rings that were collapsed for layout
  speed now re-expand correctly when the user edits inside the ring.
- **SVG export**: improved attribute quoting fixes a class of cases
  where exported SVG would not parse in Inkscape.

## [1.8.8] — 2026-05-04

### Fixed

- Tooltip layer for search-result thumbnails now correctly clips to the
  hover region; no more ghost tooltips after rapid mouse-out events.

## [1.8.7] — 2026-05-04

### Fixed

- Cache-bust query string (`?v=BIME_VERSION`) now applies to the bundle
  URL too, not just `editor/*.js` and `js/nav.js`. Returning users get
  the new bundle as soon as the workbench loads, eliminating a class
  of stale-bundle SRI mismatches on release day.

## [1.8.4] — 2026-05-04

### Changed

- ML-refiner default tuning: small accuracy improvement on chemist-drawn
  and metabolite-style inputs while staying neutral on BIME templates.
  Opt-in via `Layout.options.useMLDepict = true`.

## [1.8.2] — 2026-05-04

### Fixed

- Ring-match toggle in the workbench is now wired through to the
  substructure search code-path; the previous setting was being read
  but never propagated to the matcher.

## [1.8.1] — 2026-05-04

### Fixed

- Atom-bar build path is now defensive against missing `_atomBarButtons`
  cache; no more blank toolbar on first load with a corrupted prefs blob.

## [1.8.0] — 2026-05-04

### Changed

- **Single-bundle workbench load.** `dist/bime.min.js` now ships every
  editor module the workbench needs; the page makes ~5 HTTP requests
  on cold load instead of the previous ~30. Bundle is ~665 KB
  uncompressed and is fronted by an SRI integrity hash.
- ML refiner internals refreshed; opt-in via the same flag as before.

## [1.7.1] — 2026-05-04

### Fixed

- Layout regression on six-ring polycycles introduced in v1.7.0.

## [1.7.0] — 2026-05-04

### Changed

- ML 2D-coord refiner refresh — opt-in via
  `Layout.options.useMLDepict = true`. Default is OFF so canonical
  layouts stay byte-identical for users who want the deterministic
  rule-based pipeline only.

## [1.6.3] — 2026-05-04

### Fixed

- `prefers-reduced-motion` now disables every animation in the
  workbench, including the toolbar slide-in and the tooltip fade.

## [1.6.2] — 2026-05-03

### Fixed

- WCAG 2.1 AA contrast pass on every interactive control.

## [1.6.1] — 2026-05-03

### Fixed

- v1.6.0 shipped with an undefined `_atomBarButtons` cache on first
  load; the toolbar would render blank until the user navigated away
  and back. Defensive init added.

## [1.6.0] — 2026-05-03

### Added

- Optional ML 2D-coord refiner (`Layout.options.useMLDepict = true`).
  Pre-trained graph-conv MLP layered onto the deterministic rule-based
  layout. Off by default — opt in for chemist-drawn / metabolite-style
  inputs.

## [1.5.2] — 2026-05-03

### Fixed

- Cross-browser SVG-export attribute escaping.

## [1.5.1] — 2026-05-03

### Fixed

- Typo fix in the docs.html stylesheet that hid the right-hand TOC.

## [1.5.0] — 2026-05-03

### Added

- Per-version JS bundle hash in `dist/SRI.txt` for reproducible builds.
- `tools/sign-release.sh` helper for signing release bundles.

## [1.4.4] — 2026-05-03

### Fixed

- Atom-mapping: edge-case in symmetric reactant pairs where the
  ReactionDecoder algorithm would produce a non-canonical mapping.

## [1.4.3] — 2026-05-03

### Added

- New canonical-SMILES round-trip test cases for charged species.

## [1.4.2] — 2026-05-03

### Added

- Substructure-search performance test suite. Covers VF2++ ring-aware
  vs. ring-unaware modes, MCS, and SMARTS query matching.

## [1.4.1] — 2026-05-03

### Fixed

- VF2++ neighbour-list iteration order made deterministic so the
  same `(query, target)` pair gives the same first match across runs.

## [1.4.0] — 2026-05-03

### Added

- Atom-Atom Mapping (AAM) via the ReactionDecoder algorithm. Maps
  reactants → products for a given reaction SMILES. Exposed through
  the workbench reaction tab and the `RDT` JS module.

## [1.3.1] — 2026-05-02

### Fixed

- v1.3.0 shipped a stale `RDT.version` constant; pinned tests now
  match the canonical version stamp.

## [1.3.0] — 2026-05-02

### Added

- Maximum Common Substructure (MCS) via SMSD. Configurable ring-match
  and chirality-match flags. Exposed through the workbench MCS tab and
  the `SMSD` JS module.

## [1.2.0] — 2026-05-02

### Added

- Substructure search via VF2++. Configurable ring-match flag.
  Exposed through the workbench search tab.

## [1.1.2] — 2026-05-02

### Changed

- Build pipeline reproducibility: SRI hashes, signed bundles, and
  cross-browser verification reports.

## [1.1.1] — 2026-04-29

### Fixed

- Editor wheel-event listener now uses `{ passive: true }` so the
  workbench scrolls smoothly on trackpads.

## [1.1.0] — 2026-04-22

### Added

- SMILES round-trip parser and writer. Aromatic perception, charged
  species, isotopes, stereochemistry (E/Z, R/S).

## [1.0.0] — 2026-04-15

### Added

- First public release of the BIME workbench. 2D structure editor,
  atom/bond inspector, ring templates, common-molecules palette,
  PNG / SVG / MOL / SMILES export.
