# Contributing to BIME

Thank you for your interest in BIME — the BioInception Molecular Editor.
We welcome contributions from chemists, software engineers, educators,
students, and researchers around the world.

This guide explains how to contribute effectively. Please read it before
opening an issue or pull request.

---

## Table of contents

1. [Code of conduct](#code-of-conduct)
2. [Ways to contribute](#ways-to-contribute)
3. [Reporting bugs](#reporting-bugs)
4. [Suggesting enhancements](#suggesting-enhancements)
5. [Asking questions](#asking-questions)
6. [Development setup](#development-setup)
7. [Project layout](#project-layout)
8. [Running the tests](#running-the-tests)
9. [Coding style](#coding-style)
10. [Commit conventions](#commit-conventions)
11. [Pull request workflow](#pull-request-workflow)
12. [Adding a new molecule to the database](#adding-a-new-molecule-to-the-database)
13. [Adding a new feature](#adding-a-new-feature)
14. [Reviewing changes](#reviewing-changes)
15. [Licensing](#licensing)

---

## Code of conduct

By participating in this project you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md). In short: be respectful, be
inclusive, and assume good intent.

---

## Ways to contribute

You do **not** need to write code to contribute. Useful contributions
include:

- **Reporting bugs** with a minimal reproducible example.
- **Suggesting features** with a clear use case and proposed UX.
- **Improving documentation** — typo fixes, clarifications, examples,
  translations of `docs/USAGE.md` and the audience guides.
- **Adding molecules** to `common-molecules.js` (PubChem-canonical SMILES,
  with provenance).
- **Writing tests** for edge cases that are currently uncovered.
- **Reviewing pull requests** from others.
- **Sharing screenshots / classroom case studies** that we can include in
  `docs/EDUCATORS.md`, `docs/STUDENTS.md`, or `docs/RESEARCHERS.md`.

---

## Reporting bugs

**Before opening an issue**, please:

1. Check the [open issues](https://github.com/asad/bime-dist/issues) and
   [closed issues](https://github.com/asad/bime-dist/issues?q=is%3Aissue+is%3Aclosed)
   to see if your bug has already been reported.
2. Try the latest `main` build at the
   [live demo](https://asad.github.io/bime/workbench.html) to confirm the
   bug still reproduces.
3. Open a browser console (F12 / Cmd-Opt-I) and copy any error messages.

**When opening an issue**, use the **Bug report** template and provide:

- BIME version (`RDT.version` in the browser console, or check
  `versions.json` if self-hosted).
- Browser + OS + version.
- Exact steps to reproduce (numbered list).
- Expected behaviour vs actual behaviour.
- A minimal SMILES / SMARTS / reaction string that triggers the bug.
- Console errors (if any).
- Screenshot or short screen-recording, if a UX bug.

The more reproducible your report, the faster we can fix it.

---

## Suggesting enhancements

Use the **Feature request** template. A great enhancement request includes:

- **The problem** you are trying to solve (one paragraph, no jargon).
- **The use case** — who benefits? Students, teachers, researchers,
  industry?
- **Proposed solution** — UI sketch, API shape, or behaviour description.
- **Alternatives considered** — what would happen if BIME did nothing?
- **Out of scope** — what this enhancement is *not* asking for.

Enhancement issues that are well-specified and have community
agreement are more likely to be merged quickly.

---

## Asking questions

Use the **Question** template. Before asking, please check:

- [docs/USAGE.md](docs/USAGE.md) — comprehensive user guide.
- [docs/EDUCATORS.md](docs/EDUCATORS.md) — for teachers/professors.
- [docs/STUDENTS.md](docs/STUDENTS.md) — for students.
- [docs/RESEARCHERS.md](docs/RESEARCHERS.md) — for research workflows.
- [docs/HOSTING.md](docs/HOSTING.md) — for deployment / self-hosting.
- [docs/EMBED.md](docs/EMBED.md) — for embedding BIME in your own page.
- [README.md](README.md) — quick-start and feature overview.

Questions answered in the docs may be closed with a pointer.

---

## Development setup

BIME has **zero runtime dependencies** and **zero build dependencies**
beyond Node.js (for the test runner and bundler). You do not need npm,
webpack, vite, or anything else.

### Prerequisites

- **Node.js** ≥ 14 (any modern LTS works; we use only ES5-compatible
  syntax in source files for browser compatibility).
- **A modern browser** (Chrome, Firefox, Safari, Edge — current and
  previous major version).
- **Git** to clone the repo.

### Clone and run locally

```bash
git clone https://github.com/asad/bime-dist.git
cd bime

# Open directly in a browser — no server needed:
open workbench.html              # macOS
xdg-open workbench.html          # Linux
start workbench.html             # Windows

# Or serve locally if your browser blocks file:// for some features:
python3 -m http.server 8080
# then open http://localhost:8080/workbench.html
```

That's it. No `npm install`, no transpile step, no watch mode.

---

## Project layout

```
bime/
├── README.md                 Public README, quick start
├── CONTRIBUTING.md           This file
├── CODE_OF_CONDUCT.md        Community standards
├── SUPPORT.md                Where to get help
├── LICENSE.txt               Apache 2.0
├── NOTICE                    Required Apache 2.0 attribution notice
├── CITATION.cff              Citation metadata
├── CHANGELOG.md              All notable changes per Keep-a-Changelog
├── versions.json             Single source of truth for bundle version
│
├── index.html                Landing page
├── workbench.html            Main workbench app
├── examples.html             Interactive examples gallery
├── screenshots.html          Visual gallery
├── docs.html                 In-browser API docs
├── test.html                 Browser test harness
│
├── docs/
│   ├── USAGE.md              Comprehensive user guide
│   ├── EDUCATORS.md          For teachers/professors
│   ├── STUDENTS.md           For students learning chemistry
│   ├── RESEARCHERS.md        For research workflows
│   ├── HOSTING.md            How to host BIME yourself
│   └── EMBED.md              How to embed BIME in your page
│
├── editor/                   Source modules (loaded as <script> tags)
│   ├── Molecule.js           Atom/Bond/Molecule data model
│   ├── SmilesParser.js       SMILES recursive-descent parser
│   ├── SmilesWriter.js       Canonical SMILES writer (Morgan)
│   ├── SmartsParser.js       SMARTS pattern parser
│   ├── SmartsMatch.js        VF2 substructure matcher
│   ├── CipStereo.js          R/S and E/Z assignment
│   ├── Layout.js             2D coordinate generation
│   ├── Templates.js          Ring/scaffold templates
│   ├── Renderer.js           SVG renderer
│   ├── MolEditor.js          Top-level editor controller
│   ├── History.js            Undo/redo ring buffer
│   ├── Tools.js              Drawing tools (atom, bond, ring, …)
│   ├── ImageExport.js        SVG/PNG/clipboard export
│   ├── MolfileWriter.js      MOL V2000/V3000 writer
│   ├── RDT.js                Reaction Decoder Tool — AAM port
│   ├── SMSDGraph.js          Aromaticity + ring perception
│   ├── SMSDVF2.js            VF2 graph isomorphism
│   ├── SMSDMCS.js            Maximum Common Substructure
│   ├── SMSDRings.js          SSSR ring set
│   ├── SMSDLayout.js         MCS layout
│   ├── SMSDBatch.js          Batch MCS over fingerprint clusters
│   └── SMSDVersion.js        Version metadata
│
├── tests/
│   ├── shim.js               Shared Node test scaffolding
│   ├── run-against-bundle.js Re-runs suite against dist/ bundle
│   └── test_*.js             Per-module regression tests (~21 files)
│
├── tools/
│   ├── build.js              Pure-Node bundle builder
│   ├── run-tests.js          Test runner
│   ├── audit-aromatic.js     Aromaticity audit over all molecules
│   └── sign-release.sh       Optional GPG signing helper
│
├── dist/                     Built bundles (regenerated by tools/build.js)
│   ├── bime.js               Unminified concat bundle
│   ├── bime.min.js           Light-minified bundle
│   ├── MANIFEST.sha256       SHA-256 manifest
│   └── SRI.txt               SRI sha384 strings for HTML integrity
│
├── images/                   Brand and toolbar SVG icons
├── css/                      Stylesheet
└── js/                       Top-level page scripts (nav, theme)
```

---

## Running the tests

```bash
# Run the full regression suite (~500 tests, < 1 second):
node tools/run-tests.js

# Audit aromaticity across all 1181 built-in molecules:
node tools/audit-aromatic.js

# Build the dist bundle and re-run tests against it:
node tools/build.js
```

Every pull request must pass `node tools/run-tests.js` cleanly. Add new
tests for any new behaviour or bug fix.

### Browser tests

Open `test.html` directly in any browser. The page will load every test
file and show pass/fail counts in the page.

### Bundle parity tests

`node tools/run-tests.js` runs only against source files. To verify the
minified bundle behaves identically to source:

```bash
node tests/run-against-bundle.js
node tests/run-against-bundle.js --min   # against dist/bime.min.js
```

Both must produce the same pass count as the source-file run.

---

## Coding style

- **ES5-compatible syntax** in all `editor/*.js` source files. No `let`,
  no `const`, no arrow functions, no template literals, no `class`,
  no destructuring, no spread, no `async`/`await`. Reason: BIME runs
  unmodified on every browser back to IE 11 and on `file://` URIs.
- **No build step.** Source files are loaded directly via `<script>`
  tags; they must work as-is.
- **No external dependencies.** No npm packages at runtime.
  `tools/*.js` may use Node built-ins only (`fs`, `path`, `crypto`).
- **4-space indent** for `editor/*.js`. 2-space indent for HTML/CSS.
- **JSDoc-style comments** at the top of each module describing its
  purpose, public API, and any algorithmic citations.
- **Apache 2.0 header** at the top of every `editor/*.js` and
  `tests/*.js` file.
- **Strict equality** (`===` / `!==`) — never `==`.
- **No `eval`, `new Function`, `document.write`, `innerHTML` on
  user-controlled data.** All DOM construction goes through
  `createElement` + `textContent` or `createElementNS` for SVG.
- **CSP-safe.** The default workbench CSP forbids inline scripts at
  runtime (only static page-level `<script>` tags); your code must work
  under that policy.
- **Deterministic.** Two runs with identical inputs must produce
  identical outputs. Avoid `Math.random`, `Date.now` in algorithm code.
- **Algorithm citations** belong in source comments. If you implement a
  published algorithm, cite the paper.

---

## Commit conventions

- One logical change per commit. Bug fixes and refactors should be
  separate commits.
- Subject line: imperative mood, ≤ 72 chars. Examples:
  - `fix Naproxen SMILES to PubChem CID 156391 canonical form`
  - `add Bond.clone() preserves data field test`
  - `refactor SMSDRings to use iterative DFS for >10k atom chains`
- Body (optional, ≤ 80 chars per line) explains *why* the change is
  needed. Mention any issue / discussion the commit closes.
- Author: use your real name and a working email address. We do not
  accept anonymous or pseudonymous commits.

---

## Pull request workflow

1. **Fork** the repo on GitHub.
2. **Branch** from `main` with a descriptive name:
   `git checkout -b fix-naproxen-smiles`
3. **Make your changes** with the coding-style rules above.
4. **Run the tests**: `node tools/run-tests.js` must report 0 failures.
5. **Update CHANGELOG.md** under an `[Unreleased]` heading or under the
   next version's section.
6. **Update CITATION.cff and versions.json** only if you are bumping
   the bundle version (typically a maintainer task).
7. **Commit** with a clear message (see above).
8. **Push** to your fork and **open a pull request** against
   `asad/bime-dist:main`.
9. **Fill in the PR template** — describe what changes, why, how
   tested, and any linked issues.
10. **Respond to review comments** promptly. Maintainers may request
    changes; this is normal and not personal.

PRs should be **small and focused**. A 30-line fix with a regression test
is easier to review and merge than a 1000-line refactor.

---

## Adding a new molecule to the database

`common-molecules.js` is a curated list of 1181 PubChem-canonical SMILES.
To add a molecule:

1. Find the **PubChem CID** for the molecule.
2. Take the **canonical SMILES** from PubChem (use the **Canonical SMILES**
   field, or the **Isomeric SMILES** field if stereochemistry matters).
3. Add an entry to the appropriate category array in
   `common-molecules.js` with the form:

   ```js
   {"name":"YourMolecule","smiles":"...","category":"existing_category"}
   ```

4. Run `node tools/audit-aromatic.js` and ensure your addition passes
   (0 issues across all 1181+ molecules).
5. Run `node tools/run-tests.js` to confirm regression suite still
   passes.
6. Commit with a message like:
   `add Tubocurarine (PubChem CID 6000) to alkaloid category`.
7. Open a PR.

If your molecule fits into an entirely new category, add the category
array and update the in-code provenance comment at the top of the file.

---

## Adding a new feature

For non-trivial features:

1. **Open a feature-request issue first** to discuss design before
   writing code. This avoids wasted effort if the maintainers prefer a
   different approach.
2. Once the design is agreed:
   - Implement the feature in the relevant `editor/*.js` module.
   - Add tests in `tests/test_<module>.js` (or a dedicated
     `test_<feature>.js`).
   - Update `docs/USAGE.md` with usage examples.
   - Add a `[Unreleased]` CHANGELOG entry.
3. **Default behaviour must not change** for existing callers. New
   features are additive — opt-in via options where possible.
4. **Document the public API delta** in the CHANGELOG.

---

## Reviewing changes

Pull request reviewers should check:

- **Tests pass**: `node tools/run-tests.js` — 0 failed.
- **Style conformant**: ES5-only syntax, no new dependencies.
- **CHANGELOG updated**: a clear entry under `[Unreleased]`.
- **No IP-sensitive content**: do not merge anything that leaks
  internal codenames, unpublished research, or BioInception confidential
  material.
- **Citations present**: any new algorithm cites its source paper in
  source comments.
- **No regressions**: existing tests continue to pass.

---

## Licensing

By submitting a contribution to BIME, you agree that:

- Your contribution will be licensed under the
  [Apache License, Version 2.0](LICENSE.txt) — the same licence as the
  rest of BIME.
- You have the right to submit it (you wrote it, or you have the
  copyright owner's permission to submit it under Apache 2.0).
- You are not knowingly introducing patented or otherwise restricted
  algorithms without the patent owner's permission.

If you are unsure whether your contribution is acceptable, open an
issue and ask **before** investing significant time.

---

## Contact

- **Issues**: https://github.com/asad/bime-dist/issues
- **Maintainer**: Syed Asad Rahman, BioInception PVT LTD, Cambridge, UK

Thank you for contributing to open-source chemistry. ❤️
