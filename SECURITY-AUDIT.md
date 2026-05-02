# BIME Security, License, and Privacy Audit

**Last updated:** 2026-05-01 (v1.0.1 release)
**Original date:** 2026-03-29
**Auditor:** Automated security review (8-agent panel)
**Scope:** All files in the BIME repository

> **v1.0.1 update (2026-05-01).** This document was originally written for
> v1.0.0 and has been amended to reflect the changes in the v1.0.1 release.
> The summary of v1.0.1 changes is at the **bottom** of this file under
> **§6 — v1.0.1 update**. The rest of the document remains in its original
> form for historical traceability.

---

## 1. SECURITY FINDINGS

### 1.1 XSS Vulnerabilities (Fixed)

**SEV-01 (HIGH) -- `_showRDTStatus` innerHTML injection in MolEditor.js (line ~998)**
- The `msg` parameter was injected directly into `innerHTML` without escaping.
- The `msg` value originates from RDT server JSON responses (`data.error`, `data.mapped`).
  If a user connects to a malicious RDT server, the server could return HTML/script payloads
  in the `error` or `mapped` fields, achieving XSS in the editor context.
- **Fix applied:** Replaced innerHTML with `document.createElement` and `textContent` assignment.
  The close button now uses `addEventListener` instead of an inline `onclick` attribute.

**SEV-02 (HIGH) -- test.html RDT health response innerHTML injection (line ~138)**
- Server response fields `d.engine` and `d.version` were interpolated into `innerHTML`.
- A malicious RDT server could inject HTML via these fields.
- **Fix applied:** Replaced with DOM construction using `textContent` and `String()` coercion.

**SEV-03 (LOW) -- `updateStats` innerHTML in test.html and workbench.html**
- Used innerHTML to display atom/bond counts. While the values are integers from internal
  methods (not user-controllable), innerHTML is unnecessary here.
- **Fix applied:** Replaced with `textContent` in both files.

### 1.2 URL Injection / Protocol Abuse (Fixed)

**SEV-04 (MEDIUM) -- User-controllable RDT server URL without protocol validation**
- In workbench.html, the RDT URL input field value was used directly in `fetch()` calls.
  A `javascript:` or `data:` URL could potentially be abused (though `fetch()` itself
  rejects non-HTTP protocols, the validation adds defense in depth).
- In MolEditor.js, the `rdtServerUrl` property was used directly without validation.
- **Fix applied:** Added `_safeRdtUrl()` method that validates the URL starts with
  `http://` or `https://`, falling back to `http://localhost:8766`. Applied the same
  validation in workbench.html's `rdtGetUrl()` function.

### 1.3 SMILES Parser DoS Prevention (Fixed)

**SEV-05 (MEDIUM) -- No input length limit on SMILES parser**
- The SMILES parser had no guard against extremely long input strings. A SMILES string
  with thousands of atoms could cause the parser, layout engine, and renderer to consume
  excessive CPU time, effectively freezing the browser tab.
- The parser itself does NOT contain infinite loops -- it is a single-pass linear scanner
  with well-defined termination (the `while (pos < len)` loop always advances `pos`).
  Ring closures and branches use finite stacks. The layout algorithm uses BFS with a
  visited set. No algorithmic infinite-loop risk was found.
- **Fix applied:** Added `MAX_SMILES_LENGTH = 10000` guard at the `parse()` entry point.
  Strings exceeding this limit produce a parse error immediately.

### 1.4 Items Reviewed -- No Issues Found

**No `eval()`, `new Function()`, or `document.write()` calls anywhere in the codebase.**

**No unsafe DOM manipulation with user-supplied data:** All remaining `innerHTML` uses
either inject static HTML templates (toolbar construction, status bar reset) or use
developer-controlled constants (icon SVG strings, BRAND_COLOR). The toolbar buttons
use `item.label` which comes from hardcoded TOOLBAR_GROUPS -- not user input.

**No credential leaks in fetch calls:** The RDT fetch calls do not send cookies,
authentication headers, or credentials. They use plain JSON POST with
`Content-Type: application/json`. No `credentials: 'include'` option is set.

**localStorage usage is safe:** Only `bime-theme` key is stored (value: `'light'` or
`'dark'`). No sensitive data, molecule data, or user information is persisted.

### 1.5 RDT Server (RDTServer.java) -- Advisory Findings

**CORS wildcard (Advisory):** The Java server sets `Access-Control-Allow-Origin: *`.
This is acceptable for a local development server running on localhost. If this server
is ever deployed publicly, the CORS origin should be restricted to specific domains.

**No input sanitisation on SMILES (Advisory):** The server passes user-supplied SMILES
directly to the ReactionDecoder CLI via ProcessBuilder. While this is not a direct
shell injection (ProcessBuilder does not invoke a shell), unusually crafted input could
potentially trigger bugs in ReactionDecoder itself. The `escapeJson()` method properly
escapes output. For production deployment, consider adding SMILES input validation
and a processing timeout.

**Thread pool size (Advisory):** The server uses a fixed thread pool of 4. A burst of
slow requests could exhaust the pool. Consider adding request timeouts for production.

---

## 2. LICENSE COMPLIANCE

### 2.1 Apache 2.0 Header Coverage -- PASS

All 11 JavaScript source files in `editor/` contain the Apache 2.0 header:
```
Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
All rights reserved. Licensed under the Apache License, Version 2.0 -- see LICENSE.txt
```

Files verified:
- editor/Molecule.js
- editor/Layout.js
- editor/Templates.js
- editor/SmilesParser.js
- editor/SmilesWriter.js
- editor/MolfileWriter.js
- editor/Renderer.js
- editor/History.js
- editor/Tools.js
- editor/CipStereo.js
- editor/MolEditor.js
- editor/SmartsParser.js
- editor/SmartsMatch.js
- editor/SmartsWriter.js
- rdt-service/RDTServer.java
- common-molecules.js

The `js/nav.js` file does NOT contain a license header. This should be added for
completeness, though it is a small utility file.

### 2.2 LICENSE.txt -- PASS

Contains the full Apache License 2.0 text with correct copyright:
```
Copyright 2026 BioInception PVT LTD
```

### 2.3 LICENSE-HEADER.txt -- PASS

Contains the short-form Apache 2.0 boilerplate with dual copyright:
```
Copyright (c) 2026 BioInception PVT LTD
Algorithm Copyright (c) 2026 Syed Asad Rahman
```

### 2.4 NOTICE File -- PASS with Advisory

The NOTICE file meets Apache 2.0 Section 4(d) requirements:
- Contains project name and copyright notice
- Contains attribution requirement for academic and downstream use
- Lists novel algorithms
- Lists related tools

**Advisory:** The NOTICE file includes a TRADEMARK NOTICE section and an
"ATTRIBUTION REQUIREMENT" section that goes beyond standard Apache 2.0.
The Apache 2.0 license does not require citation in academic work -- that is
a request, not a license obligation. The current wording ("you MUST cite")
could be confusing. Consider rephrasing to "we request that you cite" or
adding this as a separate CITATION.cff file to avoid conflating license
obligations with citation requests.

### 2.5 Third-Party Code Attribution -- PASS

**No third-party code was found that requires separate attribution.**

Specific checks performed:
- **CDK (Chemistry Development Kit):** Layout.js references "Helson's SDG review"
  and "CDK StructureDiagramGenerator architecture" as algorithmic inspiration in
  comments. The implementation is original JavaScript code, not a port or copy of
  CDK Java source. CDK is LGPL-licensed; no LGPL code was detected.
- **RDKit:** No RDKit code patterns found. The SMILES parser, Morgan algorithm,
  and SSSR implementation are original.
- **OpenBabel:** No OpenBabel code patterns found.
- **JSME (Peter Ertl / Novartis):** MolEditor.js mentions "Drop-in API-compatible
  replacement for JSApplet.JSME" in a comment, indicating API compatibility but
  the implementation is entirely original.
- **PubChem/DrugBank/ChEMBL:** The common-molecules.js file sources SMILES strings
  from these databases and properly attributes them in comments. SMILES strings
  are factual chemical representations and not copyrightable.

### 2.6 README License Section -- PASS

The README correctly states:
- Apache 2.0 License
- Copyright holders (BioInception PVT LTD and Syed Asad Rahman)
- Links to LICENSE.txt
- Citation information

---

## 3. PRIVACY FINDINGS

### 3.1 External Network Requests -- PASS (Zero)

**The site makes NO external network requests.** Verified:

- No CDN links (no Google Fonts, no Bootstrap CDN, no jQuery CDN)
- No analytics scripts (no Google Analytics, no Matomo, no Plausible)
- No tracking pixels or beacons
- No external CSS imports
- No external JavaScript imports
- All resources (CSS, JS, SVG icons) are served from the same origin
- The only external links are navigation links to GitHub (which are standard
  `<a>` tags with `target="_blank" rel="noopener"` -- they do not load resources)

The `index.html` Quick Start section shows example `<script>` tags pointing to
`https://asad.github.io/bime/...` but these are inside a `<pre>` code block
(display only, not executed).

### 3.2 Tracking / Telemetry -- PASS (None)

No tracking, telemetry, fingerprinting, or analytics code found anywhere.

### 3.3 RDT Integration Privacy -- PASS with Advisory

**Molecule data sent to RDT stays on localhost by default.** The default
`rdtServerUrl` is `http://localhost:8766`. SMILES strings are sent via
`POST /api/map` as `{"smiles": "..."}`.

**Advisory:** If a user changes the RDT URL to point to a remote server,
molecule data (reaction SMILES) would be transmitted to that server. This
is by design and under user control. Consider adding a visual indicator
or warning when the RDT URL is not localhost, so users are aware their
molecule data is leaving the local machine.

### 3.4 localStorage -- PASS

Only one key is stored: `bime-theme` with value `'light'` or `'dark'`.
No molecule data, SMILES strings, user identifiers, or session data are
persisted in localStorage or sessionStorage.

### 3.5 Cookies -- PASS (None)

No cookies are set or read anywhere in the codebase.

---

## 4. SUMMARY OF CHANGES MADE

| File | Change | Severity |
|------|--------|----------|
| editor/MolEditor.js | `_showRDTStatus`: replaced innerHTML with textContent DOM construction | HIGH |
| editor/MolEditor.js | Added `_safeRdtUrl()` with http/https protocol validation | MEDIUM |
| editor/MolEditor.js | Replaced `this.rdtServerUrl` with `this._safeRdtUrl()` in all fetch calls | MEDIUM |
| editor/SmilesParser.js | Added `MAX_SMILES_LENGTH = 10000` guard against DoS | MEDIUM |
| test.html | RDT health display: replaced innerHTML with textContent DOM construction | HIGH |
| test.html | `updateStats`: replaced innerHTML with textContent | LOW |
| workbench.html | `updateStats`: replaced innerHTML with textContent | LOW |
| workbench.html | `rdtGetUrl()`: added http/https protocol validation | MEDIUM |

---

## 5. RECOMMENDATIONS

1. **Add license header to js/nav.js** for completeness.
2. **Rephrase NOTICE citation requirement** from "you MUST cite" to "we request"
   to avoid conflating Apache 2.0 obligations with academic citation norms.
   Consider adding a CITATION.cff file instead.
3. **For production RDT deployment:** restrict CORS to specific origins, add
   SMILES input length limits server-side, add request timeouts.
4. **Consider Content Security Policy (CSP):** Adding a CSP meta tag to the HTML
   pages would provide an additional layer of XSS defense. A suitable policy:
   `default-src 'self'; style-src 'self' 'unsafe-inline'; script-src 'self' 'unsafe-inline'`
5. **Visual indicator for remote RDT:** Warn users when RDT URL is not localhost.

---

## 6. v1.0.1 update (2026-05-01)

The following items from §5 are addressed in BIME v1.0.1. The full release
notes are in `CHANGELOG.md`.

| Recommendation | v1.0.1 status |
|---|---|
| §5.1 License header on `js/nav.js` | **Done.** Apache 2.0 header added. Also added to `common-molecules.js`. |
| §5.2 Rephrase NOTICE; add CITATION.cff | **Done.** NOTICE rewritten — "MUST cite/retain" mandates removed; "NOVEL ALGORITHMS" self-disclosure block deleted (gratuitous prior-art surface). New `CITATION.cff` created. |
| §5.3 Production RDT/SMSD hardening | **Done.** Both Java servers gain SMILES length cap, request body cap, 30 s process/operation timeout, bounded `ThreadPoolExecutor` (`AbortPolicy` → 503 on overload), and a CORS allowlist via `BIME_CORS_ORIGINS`. SMSDServer parsers made `ThreadLocal` (CDK is not thread-safe). RDTServer charset switched to explicit `StandardCharsets.UTF_8`. |
| §5.4 Content-Security-Policy | **Done.** CSP meta tag added to all six HTML pages: `default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self' https://api.bioinceptionlabs.com http://localhost:8766 http://localhost:8767; object-src 'none'; base-uri 'self'; frame-ancestors 'none'`. Inline `'unsafe-inline'` retained for now because of inline page bootstraps; tightening is a v1.0.2 candidate. |
| §5.5 Remote-RDT visual indicator | Deferred. The host allowlist in `_safeRdtUrl()` / `_safeSmsdUrl()` mitigates the original concern (only `api.bioinceptionlabs.com` and localhost can be configured). |

### 6.1 Privacy correction

The original §3.1 claim of *"zero external network calls"* was incorrect for
v1.0.0 — three external origins were actually contacted: `cactus.nci.nih.gov`
(NCI Chemical Identifier Resolver, used for IUPAC / common-name lookup on
every drawn SMILES), `api.bioinceptionlabs.com` (configured RDT/SMSD cloud
endpoints), and `https://cjrtnc.leaningtech.com/3.0/cj3loader.js` (CheerpJ
loader, only reachable through the dead `WasmBridge.js`).

**v1.0.1 restores the privacy posture:**

- `cactus.nci.nih.gov` — **removed**. All four call sites stripped from
  `editor/MolEditor.js`. The four public methods (`getIUPACName`,
  `getCommonNames`, `_lookupIUPACName`, `_cactusNameLookup`) are kept as
  no-op-compatible stubs so existing third-party integrations do not crash.
- `cjrtnc.leaningtech.com` — **removed**. `editor/WasmBridge.js` deleted
  outright (391 LOC of dead CheerpJ scaffold, never wired in by any caller,
  no SRI integrity attribute, no JARs in the repo).
- `api.bioinceptionlabs.com` — **retained**, but now subject to a host
  allowlist in `_safeRdtUrl()` / `_safeSmsdUrl()` (protocol AND host both
  validated). Default endpoints unchanged. No molecule data is sent
  off-host unless the user explicitly enables RDT or SMSD cloud mode.

### 6.2 Additional hardening shipped in v1.0.1

- **`workbench.html` molecule-library XSS hardening.** `buildMolBrowser`
  was concatenating `innerHTML` from `common-molecules.js` data with only
  a single-quote escape; an entry containing `<`, `>`, `&`, or `"` would
  break parsing or open stored-XSS. Migrated to DOM construction with
  `textContent` (mirrors the `buildMolBrowser` pattern at the molecule
  list site).
- **`editor/ImageExport.js` `document.write` removed.** PDF print path
  rewritten to use a Blob URL — eliminates cross-origin warnings and
  matches modern browser hardening guidance.
- **`test.html` `_safeRdtUrl()` bypass closed.** RDT health probe now
  routes through `editor._safeRdtUrl()` + `_fetchWithTimeout()` like the
  rest of the editor.
- **`AbortController`-based `_fetchWithTimeout`** added (default 8 s).
  Used for any external call so a slow / unreachable server cannot hang
  the browser tab.
- **Memory leaks:** `ImageExport._measureSVG` was attached to
  `document.body` and never detached; rewritten with reference-counted
  attach/detach. `editor/History.js` moved from `Array.shift()` cap to a
  ring-buffer (O(1) push/pop, no large-snapshot retention).

### 6.3 Items deferred to v1.0.2 / v2.0.0

- **Aromaticity perception duplicated** between `editor/SmilesWriter.js`
  (`PREV+NEXT` ring-bond inspection) and `editor/SmartsMatch.js`
  (`NEXT`-only). Documented with a `// TODO(v2.0.0):` comment at the
  divergent site. Consolidation is deferred to v2.0.0 because it touches
  algorithmic surface.
- **Inline event handlers (`onclick=`) on hardcoded SMILES literals** in
  `workbench.html`, `examples.html`, `index.html`, `test.html`. No XSS
  surface today (no user input flows in), but a `script-src 'self'`
  CSP without `'unsafe-inline'` would forbid them. Wholesale migration
  is a v1.0.2 candidate.

## 7. v1.1.3 re-audit (2026-05-02)

Re-audit performed at HEAD (`5593e57`) on the v1.1.3 working tree by the
Security + Privacy panel agent. Scope: full `editor/*.js` (24 files
including v1.1.0 `RDT.js` and v1.1.3 `Templates.js` updates), `js/nav.js`,
all six HTML pages, `tools/build.js`, `tools/sign-release.sh`,
`tools/audit-aromatic.js`, `tools/run-tests.js`, `dist/bime.js` and
`dist/bime.min.js`.

### 7.1 Re-verified clean (no regressions since v1.0.3)

- **No external network calls in shipped code.** `grep` for
  `fetch(`, `XMLHttpRequest`, `WebSocket`, `EventSource`,
  `navigator.sendBeacon`, `new Image(...http` across `editor/`, `js/`,
  and all six HTML pages returns zero hits. The dist bundle has zero
  hits as well. Editor remains 100% offline. (Tools `audit-aromatic.js`
  and `run-tests.js` also have zero `fetch`/`XHR` — they are pure
  Node.js scripts.)
- **No dynamic code execution.** Zero hits for `eval(`, `new Function(`,
  string-form `setTimeout`/`setInterval`. `document.write` appears once
  as a comment in `editor/ImageExport.js:217` describing why it was
  removed; no live call site.
- **CSP coverage.** All six HTML pages (`index.html`, `workbench.html`,
  `examples.html`, `docs.html`, `screenshots.html`, `test.html`) carry
  the identical locked-down meta CSP: `default-src 'self'; script-src
  'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src
  'self' data: blob:; connect-src 'self'; object-src 'none'; base-uri
  'self'; frame-ancestors 'none'`. No regression from v1.0.3.
- **No external `<script src="https://...">` anywhere.** SRI not
  applicable. No third-party CDN dependency.
- **Storage.** `localStorage`/`sessionStorage`/IndexedDB/cookies sweep:
  the only key is `bime-theme` in `js/nav.js` (dark/light mode), with
  try/catch around both read and write to handle private-browsing.
  No PII, no molecule data, no auth tokens stored anywhere.
- **No credentials.** Sweep for `api_key|secret|password|token|bearer|
  authorization|x-api-key` across all shipped code, docs, manifests, and
  shell scripts: zero real hits. (Only false-positive: the word
  "tokenise" in SMILES tokeniser code and "forbidden token" in
  `tools/build.js`'s self-guard.)
- **No `postMessage` / cross-origin handlers, no `<iframe>`.** Zero
  hits.
- **No unsafe URL navigation.** Zero `window.location.href = ...` or
  `location.replace` / `location.assign` call sites.
- **Build pipeline.** `tools/build.js` reads `editor/*.js`, writes only
  to `dist/`, and only spawns `node --check` for syntax verification.
  No `eval`, no `new Function`. Includes a self-guard that throws if any
  of `chemxpert`, `tinyrank`, `commercial`, `paid` survives into the
  minified bundle (forbidden marketing-token list at line 510).
- **GitHub Actions workflows.** `.github/workflows/` does not exist.
  Only `.github/REPO_DESCRIPTION.md`. `update-deps.yml` deletion (v1.1.2)
  confirmed. No hosted-runner spend triggered by any push, tag, or
  schedule. Aligns with BIME repo cost guardrails.
- **No stale local paths** (`/Users/bioinception`, `/private/tmp`,
  `.claude/`) in any tracked file.
- **No AI/LLM author trailers** (`Co-Authored-By: Claude`,
  `Generated with Claude`, `Generated by Anthropic`, etc.) in any md /
  cff / json file.
- **No internal product names** (`chemxpert`, `tinyrank`, `MULTI_ANCHOR`,
  `aam_kernel`, `BIME Pro`, `paid`, `commercial`, `tenant`,
  `future server`) anywhere outside the `tools/build.js` forbidden-token
  guard.
- **dist bundle re-check.** Zero hits for `chemxpert|tinyrank|paid|
  commercial|api_key|password` in `dist/bime.js` or `dist/bime.min.js`.
  Bundle has only `innerHTML` calls that mirror the source (all reviewed
  in this re-audit and accepted).

### 7.2 v1.1.0 / v1.1.3 surface review (new code)

- **`editor/RDT.js` (v1.1.0, 963 LOC)** — pure-JS atom-atom mapping.
  Contains no DOM access, no `innerHTML`, no `fetch`, no `eval`, no
  `localStorage`. Pure functions over molecule objects. Clean.
- **`editor/Templates.js` (v1.1.3 update, 1284 LOC)** — pure-JS template
  / fragment library. Same posture: no DOM, no network, no eval. Clean.

### 7.3 v1.1.3 fixes applied in this audit

- **`screenshots.html`** — defense-in-depth. `MOLECULES` array data is
  hardcoded inside the file (no external input), but the gallery card
  `label` was being built by `innerHTML` concatenation of `entry.name`
  and `entry.smiles`. Migrated to DOM construction with `textContent`,
  matching the pattern already used in `workbench.html:324-326` for
  the molecule library. Severity: MEDIUM (defense-in-depth — no
  exploit today, but hardens against any future entry containing
  HTML-significant characters such as ampersand, less-than,
  double-quote in a generated SMILES). Public API unchanged. Six
  showcase molecules still render identically.

### 7.4 Items still deferred (re-confirmed)

- Aromaticity perception consolidation between `SmilesWriter.js` and
  `SmartsMatch.js` (originally flagged in §6.3) is still deferred to
  v2.0.0.
- Inline `onclick=` handlers in HTML pages on hardcoded SMILES literals
  (originally flagged in §6.3) — still acceptable under current CSP
  (`script-src 'self' 'unsafe-inline'`); wholesale migration deferred.

### 7.5 Regression test status post-fix

`node tools/run-tests.js` — **427 passed, 0 failed** in 579 ms.
`node tools/build.js` — bime.js 1036171 bytes, bime.min.js 497330 bytes,
bundle test 427 passed, min bundle 427 passed. SRI regenerated in
`dist/SRI.txt`.

### 7.6 Files modified

- `screenshots.html` — replaced one `innerHTML` concat block with DOM
  construction (lines 194-209 in working tree).
- `dist/bime.js`, `dist/bime.min.js`, `dist/MANIFEST.sha256`,
  `dist/SRI.txt` — regenerated by `tools/build.js`.
- `SECURITY-AUDIT.md` — appended this §7.

### 7.7 Ancillary observation (out of audit scope)

`tools/build.js:61` still hard-codes `BIME_VERSION = '1.1.2'` while
`versions.json` declares `bime: 1.1.3`. This is a release-pipeline
drift, not a security issue, but the owner may want to bump
`BIME_VERSION` to `1.1.3` before tagging.

