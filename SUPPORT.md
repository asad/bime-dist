# Getting help with BIME

Welcome! BIME is the BioInception Molecular Editor — a free, open-source,
browser-based tool for drawing, analysing, and exporting chemical
structures, reactions, and SMARTS queries.

This page tells you **where to look first** depending on what you need.

---

## I am a …

### 🎓 Student learning chemistry

- Read [docs/STUDENTS.md](docs/STUDENTS.md) — a guide tailored for
  undergraduates and high-school students.
- Open the [live workbench](https://asad.github.io/bime-dist/workbench.html)
  and try the **Examples** tab.
- Browse the [examples gallery](https://asad.github.io/bime-dist/examples.html)
  for interactive walk-throughs of MCS, AAM, substructure search, and
  bond-change analysis.
- For SMILES syntax help, see
  [docs/USAGE.md § SMILES Quick Reference](docs/USAGE.md#smiles-quick-reference).

### 👩‍🏫 Teacher / professor

- Read [docs/EDUCATORS.md](docs/EDUCATORS.md) — classroom recipes,
  printable worksheets, and lecture-friendly examples.
- Embed BIME in your LMS (Moodle, Canvas, Blackboard) using
  [docs/EMBED.md](docs/EMBED.md).
- Host BIME on your university web space using
  [docs/HOSTING.md](docs/HOSTING.md).

### 🔬 Researcher

- Read [docs/RESEARCHERS.md](docs/RESEARCHERS.md) — workflows for
  reaction analysis, MCS-based scaffold mining, atom-atom mapping, and
  pipeline integration.
- The full programmatic API is in [docs/USAGE.md](docs/USAGE.md).
- Cite BIME using the metadata in [CITATION.cff](CITATION.cff).

### 🛠️ Developer / integrator

- Read [docs/EMBED.md](docs/EMBED.md) — how to embed BIME in your own
  page, web app, Jupyter notebook, or Observable notebook.
- Read [docs/HOSTING.md](docs/HOSTING.md) — how to host BIME on your
  own server, GitHub Pages, Netlify, Cloudflare Pages, or even via
  IPFS / `file://`.
- Read [CONTRIBUTING.md](CONTRIBUTING.md) — how to contribute code.

---

## I have a …

### 🐛 Bug to report

→ [Open a Bug report](https://github.com/asad/bime-dist/issues/new?template=bug_report.md)

Please use the template — it asks for the BIME version, browser, exact
steps to reproduce, and a minimal SMILES that triggers the problem. The
more reproducible your report, the faster we can fix it.

### 💡 Feature idea

→ [Open a Feature request](https://github.com/asad/bime-dist/issues/new?template=feature_request.md)

Describe the **problem** first (one paragraph), then the **proposed
solution**, then the **alternatives**. Feature ideas with a clear use
case and proposed UX are most likely to be merged.

### ❓ Question about how to use BIME

→ [Open a Question](https://github.com/asad/bime-dist/issues/new?template=question.md)

Before asking, please check the docs (we link them at the top of every
template). If your question is already answered there, we may close
the issue with a pointer.

### 🔒 Security vulnerability

**Please do not open a public issue for security bugs.**

Email the maintainer privately: **asad.rahman@bioinceptionlabs.com**

We will acknowledge receipt within 5 working days, investigate, and
coordinate a fix. Public disclosure of the vulnerability will be timed
so users have a chance to update first.

---

## I want to …

### Run BIME on my computer / school / lab without internet

BIME runs entirely in the browser with **zero backend** and **zero
external network calls** by default. You can:

- Open `workbench.html` directly via `file://` — no server needed.
- Host it on a school intranet using any static-file server.
- Drop it onto a USB stick or shared drive — it will run from there.

See [docs/HOSTING.md](docs/HOSTING.md) for details.

### Embed BIME in my own webpage / blog / LMS

You only need a few `<script>` tags and a `<div>`. See
[docs/EMBED.md](docs/EMBED.md) for copy-pasteable snippets and a working
HTML example.

### Translate the UI to my language

The current UI is English-only. Translation contributions are very
welcome — open an issue to discuss the approach (i18n keys vs. per-page
translation), and we'll guide you through. The user-facing strings live
in `editor/MolEditor.js`, `index.html`, `workbench.html`,
`docs.html`, and `examples.html`.

### Cite BIME in a paper / thesis

Use [CITATION.cff](CITATION.cff). For BibTeX, see the citation block at
the bottom of [README.md](README.md).

The underlying algorithms are also citable — RDT (Rahman 2016),
EC-BLAST (Rahman 2014), SMSD (Rahman 2009). Full citations are in the
README and CHANGELOG.

### Sponsor / fund BIME development

BIME is volunteer-maintained under Apache 2.0 and will remain so. If
your organisation benefits from BIME and would like to support
maintenance, please email **asad.rahman@bioinceptionlabs.com** with
"BIME Sponsorship" in the subject line. We are particularly open to
sponsorships from universities, research institutes, and educational
non-profits.

---

## What BIME is **not**

- **Not** a 3D molecule viewer (use 3Dmol.js, Mol\*, or NGL Viewer).
- **Not** a quantum chemistry calculator (use ORCA, Gaussian, Psi4).
- **Not** a docking tool (use AutoDock Vina, Smina, GNINA).
- **Not** a force-field minimiser (use OpenMM, RDKit, Open Babel).
- **Not** a database (use ChEMBL, PubChem, ChEBI).

BIME is a 2D molecule editor + cheminformatics workbench focused on
**structure drawing, SMARTS search, MCS, atom-atom mapping, and
bond-change analysis**.

---

## Response times

This is a community-maintained project. Issues are typically triaged
within 1-2 weeks. Pull requests with passing tests and a clear
description are reviewed faster.

If your issue is time-sensitive (e.g. teaching deadline, paper
submission), please mention that in your issue and we will do our best.

---

Thank you for using BIME, and for being part of the open-source
chemistry community.

— Maintainer: Syed Asad Rahman, BioInception PVT LTD, Cambridge, UK
