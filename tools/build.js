#!/usr/bin/env node
/**
 * tools/build.js - BIME v1.1.4 build pipeline (zero npm dep).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Concatenates editor/*.js in the canonical browser load order, writes
 * dist/bime.js (unminified, source-readable) and dist/bime.min.js
 * (light-minified - whitespace + comment strip, no symbol mangling).
 * Computes SHA-256 manifest and SRI sha384 strings, then re-runs the
 * regression suite against the bundle to verify behaviour
 * equivalence.
 *
 * Usage:    node tools/build.js
 * Outputs:  dist/bime.js, dist/bime.min.js, dist/MANIFEST.sha256, dist/SRI.txt
 *
 * No external dependencies. No network. No symbol mangling. No
 * obfuscation. The Apache 2.0 banner is preserved as a /*! bang
 * comment so any downstream minifier respects it.
 */
'use strict';

var fs = require('fs');
var path = require('path');
var crypto = require('crypto');
var child_process = require('child_process');

var ROOT = path.join(__dirname, '..');
var EDITOR_DIR = path.join(ROOT, 'editor');
var DIST_DIR = path.join(ROOT, 'dist');

// Canonical browser load order, taken from workbench.html.
// SMSDLayout is appended after SMSDBatch (used by Layout.js when present).
var FILES = [
    'Molecule.js',
    'Layout.js',
    'Templates.js',
    'SmilesParser.js',
    'SmilesWriter.js',
    'MolfileWriter.js',
    'Renderer.js',
    'History.js',
    'Tools.js',
    'CipStereo.js',
    'SmartsParser.js',
    'SmartsMatch.js',
    'SmartsWriter.js',
    'ImageExport.js',
    'ToolbarPrefs.js',
    'MolEditor.js',
    'SMSDVersion.js',
    'SMSDGraph.js',
    'SMSDVF2.js',
    'SMSDMCS.js',
    'SMSDRings.js',
    'SMSDBatch.js',
    'SMSDLayout.js',
    // v1.6.0: MLDepict + bundled weights (loaded eagerly so
    // Layout.options.useMLDepict can flip on at runtime).
    'MLDepict.js',
    'ml-depict-weights.js',
    'RDT.js'
];

var BIME_VERSION = '1.8.1';
var BUILD_DATE = '2026-05-03';

// ---------------------------------------------------------------------------
// Banner: combined Apache 2.0 header. /*! bang comment survives minifiers.
// ---------------------------------------------------------------------------

function makeBanner() {
    // Use BUILD_DATE (not Date.now) so the bundle is reproducible: the same
    // source + same date should yield the same bytes and the same SRI hash.
    // Override with BIME_BUILD_TIMESTAMP env var if you want to embed a more
    // precise build moment (and accept that SRI will change accordingly).
    var stamp = process.env.BIME_BUILD_TIMESTAMP || (BUILD_DATE + 'T00:00:00Z');
    return [
        '/*! bundled BIME v' + BIME_VERSION + ' */',
        '/*!',
        ' * BIME — BioInception Molecular Editor',
        ' * Version: v' + BIME_VERSION + '   Build: ' + stamp,
        ' *',
        ' * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.',
        ' * Algorithm copyright (c) 2009-2026 Syed Asad Rahman (SMSD-derived modules).',
        ' * Licensed under the Apache License, Version 2.0 - see LICENSE.txt',
        ' *',
        ' * Source remains readable in the editor/*.js files; this bundle is a',
        ' * deployment optimisation only. No symbol mangling, no obfuscation.',
        ' * Citation: Rahman SA. BIME: BioInception Molecular Editor (2026).',
        ' */',
        ''
    ].join('\n');
}

// ---------------------------------------------------------------------------
// Strip per-file leading JSDoc-style header (/** ... */ at file start).
// ---------------------------------------------------------------------------

function stripLeadingHeader(src) {
    // Remove a single /** ... */ at the very start (allow trailing newline).
    var m = src.match(/^\s*\/\*\*[\s\S]*?\*\/\s*/);
    if (m) return src.slice(m[0].length);
    return src;
}

// ---------------------------------------------------------------------------
// Concatenate sources with traceability dividers.
// ---------------------------------------------------------------------------

function concat(files) {
    var out = [makeBanner()];
    for (var i = 0; i < files.length; i++) {
        var f = files[i];
        var fp = path.join(EDITOR_DIR, f);
        var src = fs.readFileSync(fp, 'utf8');
        src = stripLeadingHeader(src);
        out.push('/*--- editor/' + f + ' ---*/');
        out.push(src);
    }
    // Bottom shim: when loaded in Node via require(), expose global symbols
    // as module.exports so a single require('dist/bime.js') populates them.
    out.push('/*--- node-bundle export shim ---*/');
    out.push([
        '(function(){',
        '  if (typeof module !== "undefined" && typeof module.exports === "object") {',
        '    module.exports = (typeof window !== "undefined") ? window : globalThis;',
        '  }',
        '})();'
    ].join('\n'));
    return out.join('\n');
}

// ---------------------------------------------------------------------------
// Light tokenizer-aware minifier.
//
// Goals:
//   1. Strip /* ... */ comments EXCEPT keep /*! ... */ banners.
//   2. Strip // ... line comments to EOL.
//   3. Collapse runs of whitespace (incl. newlines) to a single space.
//
// Rules respected:
//   - String literals (' " `) are scanned and copied verbatim.
//   - Template literals preserve their contents (no interpolation tokenizing
//     beyond brace counting).
//   - Regex literals are detected after expression-context tokens and kept
//     verbatim. Conservative: when ambiguous, treat / as division.
//   - Whitespace adjacent to identifiers/keywords is preserved as a single
//     space; whitespace adjacent to punctuation is dropped.
//   - We never rename identifiers, never rearrange tokens, never alter
//     semantics. The only changes are whitespace and comment removal.
// ---------------------------------------------------------------------------

function isIdent(ch) {
    return (ch >= 'a' && ch <= 'z') ||
           (ch >= 'A' && ch <= 'Z') ||
           (ch >= '0' && ch <= '9') ||
           ch === '_' || ch === '$';
}

function isWs(ch) {
    return ch === ' ' || ch === '\t' || ch === '\n' || ch === '\r' ||
           ch === ' ' || ch === ' ' || ch === ' ';
}

// Tokens after which a `/` starts a regex (vs. division).
var REGEX_PRECEDING = {
    '(': 1, ',': 1, '=': 1, ':': 1, '[': 1, '!': 1, '&': 1, '|': 1,
    '?': 1, '{': 1, '}': 1, ';': 1, '+': 1, '-': 1, '*': 1, '%': 1,
    '<': 1, '>': 1, '^': 1, '~': 1
};
// Keywords after which a `/` starts a regex.
var REGEX_KEYWORDS = {
    'return': 1, 'typeof': 1, 'instanceof': 1, 'in': 1, 'of': 1,
    'new': 1, 'delete': 1, 'void': 1, 'throw': 1, 'case': 1,
    'do': 1, 'else': 1, 'yield': 1, 'await': 1
};

function lastSignificantToken(out) {
    // Walk back over trailing whitespace to find the last non-ws non-comment
    // character or identifier.
    for (var i = out.length - 1; i >= 0; i--) {
        var c = out[i];
        if (!isWs(c)) return c;
    }
    return '';
}

function lastWordBeforeSlash(out) {
    // Find the trailing identifier/keyword if any, else return ''.
    var i = out.length - 1;
    while (i >= 0 && isWs(out[i])) i--;
    var end = i + 1;
    while (i >= 0 && isIdent(out[i])) i--;
    var start = i + 1;
    if (start === end) return '';
    var word = out.slice(start, end).join('');
    // Reject if it's all digits (numeric literal, not keyword).
    if (/^[0-9]/.test(word)) return '';
    return word;
}

function minify(src) {
    var n = src.length;
    var out = [];   // char array we push onto
    var i = 0;

    while (i < n) {
        var ch = src[i];
        var ch2 = src[i + 1];

        // ----- block comment /* ... */ -----
        if (ch === '/' && ch2 === '*') {
            var isBang = src[i + 2] === '!';
            var end = src.indexOf('*/', i + 2);
            if (end === -1) {
                // Unterminated - keep as-is, abort minification of this run.
                out.push(src.slice(i));
                i = n;
                break;
            }
            if (isBang) {
                // Preserve banner. Push verbatim incl. delimiters.
                out.push(src.slice(i, end + 2));
                // Add a newline so successive banners stay on separate lines.
                out.push('\n');
            }
            // else: drop comment entirely (no whitespace replacement; the
            // caller is responsible for separating tokens around the comment
            // - in practice JS comments live between statements / inside
            // declarations where surrounding whitespace already exists).
            i = end + 2;
            continue;
        }

        // ----- line comment // ... -----
        if (ch === '/' && ch2 === '/') {
            // Skip to end of line. Preserve the newline so statement
            // boundaries (ASI) are not broken.
            var nl = src.indexOf('\n', i + 2);
            if (nl === -1) { i = n; break; }
            i = nl;  // Land on \n; the whitespace branch below will handle.
            continue;
        }

        // ----- string literals ' " ` -----
        if (ch === '"' || ch === "'" || ch === '`') {
            var quote = ch;
            var start = i;
            i++;
            while (i < n) {
                var c = src[i];
                if (c === '\\') { i += 2; continue; }
                if (quote === '`' && c === '$' && src[i + 1] === '{') {
                    // Template-literal interpolation. Skip into the JS expr,
                    // tracking braces. Inside it, strings/comments are again
                    // possible - cheapest correct path is to walk char-by-char
                    // honouring nested strings/template-literals/braces.
                    i += 2;
                    var depth = 1;
                    while (i < n && depth > 0) {
                        var cc = src[i];
                        if (cc === '{') { depth++; i++; continue; }
                        if (cc === '}') { depth--; i++; continue; }
                        if (cc === '"' || cc === "'" || cc === '`') {
                            var iq = cc;
                            i++;
                            while (i < n) {
                                var ccc = src[i];
                                if (ccc === '\\') { i += 2; continue; }
                                if (ccc === iq) { i++; break; }
                                i++;
                            }
                            continue;
                        }
                        if (cc === '/' && src[i + 1] === '/') {
                            var nl2 = src.indexOf('\n', i + 2);
                            i = (nl2 === -1) ? n : nl2;
                            continue;
                        }
                        if (cc === '/' && src[i + 1] === '*') {
                            var ce = src.indexOf('*/', i + 2);
                            i = (ce === -1) ? n : ce + 2;
                            continue;
                        }
                        i++;
                    }
                    continue;
                }
                if (c === quote) { i++; break; }
                i++;
            }
            out.push(src.slice(start, i));
            continue;
        }

        // ----- regex literal -----
        if (ch === '/') {
            // Decide regex vs. division based on previous significant token.
            var lastCh = lastSignificantToken(out);
            var isRegex = false;
            if (lastCh === '' || REGEX_PRECEDING[lastCh] === 1) {
                isRegex = true;
            } else if (isIdent(lastCh)) {
                // identifier or keyword - distinguish by lookup.
                var word = lastWordBeforeSlash(out);
                if (REGEX_KEYWORDS[word] === 1) isRegex = true;
            }
            if (isRegex) {
                var rstart = i;
                i++;
                var inClass = false;
                while (i < n) {
                    var rc = src[i];
                    if (rc === '\\') { i += 2; continue; }
                    if (rc === '[') { inClass = true; i++; continue; }
                    if (rc === ']') { inClass = false; i++; continue; }
                    if (rc === '/' && !inClass) { i++; break; }
                    if (rc === '\n') { break; }  // safety: not a regex
                    i++;
                }
                // Trailing flags (gimsuy).
                while (i < n && /[a-z]/.test(src[i])) i++;
                out.push(src.slice(rstart, i));
                continue;
            }
            // else fall through, treat as division operator.
        }

        // ----- whitespace -----
        if (isWs(ch)) {
            // Collapse the run.
            var j = i;
            while (j < n && isWs(src[j])) j++;
            // Decide: do we need to keep one space?
            var prev = out.length ? out[out.length - 1] : '';
            var next = src[j] || '';
            var keepSpace = false;
            if (prev && next) {
                // Identifier-identifier boundary requires a space.
                var prevIdent = isIdent(prev);
                var nextIdent = isIdent(next);
                if (prevIdent && nextIdent) keepSpace = true;
                // `+ +` / `- -` / `+ ++` etc. require a space to avoid `++`.
                if ((prev === '+' || prev === '-') && prev === next) keepSpace = true;
                // Avoid joining `/` and `/` (would start a comment) - rare
                // since we handle comments above, but defensive.
                if (prev === '/' && next === '/') keepSpace = true;
            }
            if (keepSpace) out.push(' ');
            i = j;
            continue;
        }

        // ----- default: emit char -----
        out.push(ch);
        i++;
    }

    return out.join('');
}

// ---------------------------------------------------------------------------
// Hashing & SRI.
// ---------------------------------------------------------------------------

function sha256Hex(buf) {
    return crypto.createHash('sha256').update(buf).digest('hex');
}

function sha384Sri(buf) {
    var h = crypto.createHash('sha384').update(buf).digest('base64');
    return 'sha384-' + h;
}

// ---------------------------------------------------------------------------
// Build.
// ---------------------------------------------------------------------------

function ensureDist() {
    if (!fs.existsSync(DIST_DIR)) fs.mkdirSync(DIST_DIR, { recursive: true });
}

function writeFile(p, content) {
    fs.writeFileSync(p, content);
}

function fileSize(p) {
    return fs.statSync(p).size;
}

// ---------------------------------------------------------------------------
// Bundle test runner.
// ---------------------------------------------------------------------------

function runBundleTests(useMin) {
    var bundleTestRunner = path.join(ROOT, 'tests', 'run-against-bundle.js');
    if (!fs.existsSync(bundleTestRunner)) {
        throw new Error('tests/run-against-bundle.js missing - cannot validate bundle');
    }
    var args = [bundleTestRunner];
    if (useMin) args.push('--min');
    var res = child_process.spawnSync(process.execPath, args, {
        cwd: ROOT,
        encoding: 'utf8',
        stdio: ['ignore', 'pipe', 'pipe']
    });
    var stdout = res.stdout || '';
    var stderr = res.stderr || '';
    if (res.status !== 0) {
        process.stderr.write(stdout);
        process.stderr.write(stderr);
        throw new Error('bundle test run failed (exit ' + res.status + ')');
    }
    // Extract pass count from the FINAL summary line (after the dashes).
    var summaryRe = /(\d+)\s+passed,\s+(\d+)\s+failed,\s+total\s+\d+\s+ms/;
    var m = stdout.match(summaryRe);
    if (!m) throw new Error('could not parse bundle test summary');
    if (parseInt(m[2], 10) !== 0) {
        process.stderr.write(stdout);
        throw new Error('bundle tests had failures');
    }
    return parseInt(m[1], 10);
}

// ---------------------------------------------------------------------------
// Main.
// ---------------------------------------------------------------------------

function main() {
    ensureDist();

    // Verify all sources exist.
    for (var i = 0; i < FILES.length; i++) {
        var fp = path.join(EDITOR_DIR, FILES[i]);
        if (!fs.existsSync(fp)) {
            throw new Error('missing source: ' + fp);
        }
    }

    // Sum source bytes for the size-reduction check.
    var srcBytes = 0;
    for (var k = 0; k < FILES.length; k++) {
        srcBytes += fileSize(path.join(EDITOR_DIR, FILES[k]));
    }

    var bundle = concat(FILES);
    var bundlePath = path.join(DIST_DIR, 'bime.js');
    writeFile(bundlePath, bundle);

    // Quick syntax check via Node's --check.
    var syntax1 = child_process.spawnSync(process.execPath, ['--check', bundlePath], { encoding: 'utf8' });
    if (syntax1.status !== 0) {
        process.stderr.write(syntax1.stderr || '');
        throw new Error('syntax error in dist/bime.js');
    }

    var minified = minify(bundle);
    var minPath = path.join(DIST_DIR, 'bime.min.js');
    writeFile(minPath, minified);

    var syntax2 = child_process.spawnSync(process.execPath, ['--check', minPath], { encoding: 'utf8' });
    if (syntax2.status !== 0) {
        process.stderr.write(syntax2.stderr || '');
        throw new Error('syntax error in dist/bime.min.js after minification');
    }

    // Hashes & SRI.
    var bundleBuf = fs.readFileSync(bundlePath);
    var minBuf = fs.readFileSync(minPath);
    var sha256bundle = sha256Hex(bundleBuf);
    var sha256min = sha256Hex(minBuf);
    var sriBundle = sha384Sri(bundleBuf);
    var sriMin = sha384Sri(minBuf);

    // The two-column lines below are valid `shasum -c` input. The footer
    // line is prefixed with `# ` so shasum treats it as a comment and
    // doesn't warn. Keeps the file usable both as a manifest and as a
    // human-readable note.
    var manifest = [
        '# BIME v' + BIME_VERSION + ' - built ' + BUILD_DATE + ' by tools/build.js',
        '# verify with: shasum -a 256 -c MANIFEST.sha256',
        sha256bundle + '  bime.js',
        sha256min + '  bime.min.js',
        ''
    ].join('\n');
    writeFile(path.join(DIST_DIR, 'MANIFEST.sha256'), manifest);

    var sri = [
        '# BIME v' + BIME_VERSION + ' SRI integrity strings',
        '# Generated by tools/build.js on ' + BUILD_DATE + '. Regenerated on every release.',
        '',
        '<script src="bime.js" integrity="' + sriBundle + '" crossorigin="anonymous"></script>',
        '<script src="bime.min.js" integrity="' + sriMin + '" crossorigin="anonymous"></script>',
        ''
    ].join('\n');
    writeFile(path.join(DIST_DIR, 'SRI.txt'), sri);

    // ---------------------------------------------------------------------
    // Cache-bust local resources in HTML pages so users see a new version
    // immediately on hard-reload instead of waiting for max-age=600 to
    // expire on cached editor/*.js. We bump the `?v=...` suffix on every
    // <script src="editor/...">, <script src="js/...">,
    // <script src="common-molecules.js">, and <link href="css/...">.
    // ---------------------------------------------------------------------
    var HTML_PAGES = ['index.html', 'workbench.html', 'examples.html',
                      'docs.html', 'screenshots.html', 'test.html'];
    HTML_PAGES.forEach(function (page) {
        var p = path.join(ROOT, page);
        if (!fs.existsSync(p)) { return; }
        var html = fs.readFileSync(p, 'utf8');
        var before = html;
        // Replace existing ?v=... or add it if missing.
        html = html.replace(
            /(src="(?:editor|js)\/[^"?]*\.js)(?:\?v=[^"]*)?(")/g,
            '$1?v=' + BIME_VERSION + '$2'
        );
        html = html.replace(
            /(src="common-molecules\.js)(?:\?v=[^"]*)?(")/g,
            '$1?v=' + BIME_VERSION + '$2'
        );
        html = html.replace(
            /(href="css\/[^"?]*\.css)(?:\?v=[^"]*)?(")/g,
            '$1?v=' + BIME_VERSION + '$2'
        );
        if (html !== before) {
            fs.writeFileSync(p, html);
        }
    });

    // Run regression suite against the bundle, then against the minified bundle.
    var bundlePassed = runBundleTests(false);
    var bundleMinPassed = runBundleTests(true);
    if (bundleMinPassed !== bundlePassed) {
        throw new Error('minified bundle test count mismatch: ' + bundleMinPassed + ' vs. ' + bundlePassed);
    }

    // Sanity: banner survived minification.
    var minHead = minified.slice(0, 4096);
    if (minHead.indexOf('/*! bundled BIME v' + BIME_VERSION + ' */') === -1) {
        throw new Error('banner did not survive minification');
    }

    var unminSize = bundleBuf.length;
    var minSize = minBuf.length;
    var minPct = Math.round((minSize / unminSize) * 100);
    var srcPct = Math.round((minSize / srcBytes) * 100);

    console.log('');
    console.log('BIME v' + BIME_VERSION + ' build complete');
    console.log('---------------------------------------');
    console.log('bime.js      ' + unminSize + ' bytes  sha256 ' + sha256bundle.slice(0, 12) + '...');
    console.log('bime.min.js  ' + minSize + ' bytes  sha256 ' + sha256min.slice(0, 12) + '...   (' + minPct + '% of unminified, ' + srcPct + '% of source sum)');
    console.log('SRI:         ' + sriMin);
    console.log('Bundle test: ' + bundlePassed + ' passed');
    console.log('Min bundle:  ' + bundleMinPassed + ' passed');
    console.log('');

    if (srcPct >= 70) {
        // Soft warning - the spec asks for at least 30% reduction.
        console.warn('warning: bundle is only ' + (100 - srcPct) + '% smaller than the source sum (target >= 30%)');
    }
}

try {
    main();
    process.exit(0);
} catch (err) {
    console.error('build failed: ' + (err && err.message ? err.message : err));
    if (err && err.stack) console.error(err.stack);
    process.exit(1);
}
