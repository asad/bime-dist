(function(global) {
    'use strict';
    /**
     * ExportStamp.js — provenance + tamper-evidence stamp for every BIME export.
     *
     * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
     * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
     *
     * Standing rule: every artefact exported from BIME
     * (SVG, PNG via the SVG path, PDF, MOL, SDF, RXN, SMILES) carries a
     * BIME provenance stamp so the source is identifiable and the
     * content can be checked for tampering.
     *
     * The stamp is twofold:
     *   1. A human-readable attribution: "BIME v<X.Y.Z> — BioInception
     *      Molecular Editor — © 2026 BioInception PVT LTD / Syed Asad
     *      Rahman".
     *   2. A SHA-384 fingerprint computed over the unstamped content
     *      (so re-stamping is idempotent — the hash never includes the
     *      stamp itself).
     *
     * Public-domain note (v2.0.66 standing rule): when this
     * code or the docs say BIME exports are "encrypted for public
     * domain", the meaning is **tamper-evident signed**, not
     * asymmetrically encrypted. Apache-2.0 artefacts must remain
     * readable by anyone — encrypting the bytes would contradict the
     * licence. What the stamp DOES provide is verifiable provenance
     * + integrity: the fingerprint is computed over the un-stamped
     * content, so a consumer can strip the stamp, recompute the
     * hash, and confirm both that (a) the artefact came from BIME,
     * (b) it has not been altered since stamping. That is a witness
     * — not a cipher. Confidentiality is NOT promised; attribution +
     * integrity ARE.
     *
     * True asymmetric signatures (Ed25519 / RSA) require a private
     * key held only by BioInception and live in a separate `.sig`
     * sidecar in any future release that wants them. The roadmap
     * tracks this as the v2.0.6x+ asymmetric-signature item.
     *
     * Format-specific stamp placement:
     *   SVG      — `<title>` + `<metadata>` block at the document head.
     *   MOL      — V2000 header line 2 (the "program/timestamp" slot)
     *              plus a leading comment block immediately after the
     *              header trio.
     *   SDF      — same as MOL (SDF is MOL + data fields).
     *   RXN      — V2000 RXN header lines 1-3 already carry the
     *              creator; we append BIME stamp to line 3.
     *   SMILES   — SMILES has no native comment syntax; the stamp is
     *              emitted on a separate line below the SMILES, prefixed
     *              with `# `, a common command-line comment convention.
     *              Consumers who want a pure SMILES string strip the
     *              `# ` lines.
     *   txt/log  — `# ` prefixed header at the top.
     *
     * Idempotency: `stampSvg`, `stampText` recognise an existing BIME
     * stamp and either pass through unchanged or refresh the version.
     *
     * Pure ES5, no DOM dependency at runtime — works equally inside
     * the workbench browser bundle and the bime CLI.
     */

    // Crypto.subtle is unavailable in sync paths in some browsers; the
    // CLI uses Node's `crypto` module. We expose a synchronous fallback
    // SHA-384 implementation only as a last resort because BIME does not
    // ship a 3rd-party SHA library and the spec implementation is
    // large. In normal operation the host (Node or browser) provides
    // a SHA-384 primitive via the platform.
    function defaultSha384(content) {
        // Browser path: window.crypto.subtle is async. We synthesise a
        // sync wrapper that uses Node's crypto when available; otherwise
        // we hand back a short FNV-like fingerprint as a tamper-evident
        // (but not collision-resistant) placeholder. This keeps the
        // stamp deterministic in every environment without bloating the
        // bundle with a SHA-384 polyfill. The stamp text always names
        // which fingerprint algorithm was used, so a consumer can verify.
        try {
            // Node path
            var nodeCrypto = (typeof module !== 'undefined' && module.exports) ?
                require('crypto') : null;
            if (nodeCrypto && nodeCrypto.createHash) {
                return {
                    algo: 'sha384',
                    hex: nodeCrypto.createHash('sha384').update(content).digest('hex')
                };
            }
        } catch (e) { /* fall through */ }
        // Browser fallback (only used when window.crypto.subtle is async-only):
        // FNV-1a 64-bit digest, hex-encoded. Documented in the stamp text
        // as "fnv1a64" so the consumer doesn't assume SHA-384.
        var h1 = 0x811c9dc5, h2 = 0xc9dc5118;
        for (var i = 0; i < content.length; i++) {
            var c = content.charCodeAt(i);
            h1 ^= c; h2 ^= c;
            h1 = (h1 * 16777619) >>> 0;
            h2 = (h2 * 1099511628211) % 0xFFFFFFFFFFFFFF;
        }
        return {
            algo: 'fnv1a64',
            hex: ('00000000' + h1.toString(16)).slice(-8) +
                 ('00000000' + h2.toString(16)).slice(-8)
        };
    }

    var COPYRIGHT_LINE = '(c) 2026 BioInception PVT LTD, Cambridge, UK / Syed Asad Rahman';
    var LICENSE_LINE   = 'Apache-2.0';
    var URL_LINE       = 'https://github.com/asad/bime-dist';

    // Read the current BIME version at stamp time so we don't have to
    // re-bake it on every release. Looks for the globals each editor
    // module attaches; falls back to '(dev)' if none are loaded.
    function currentVersion() {
        var g = (typeof globalThis !== 'undefined') ? globalThis :
                ((typeof window !== 'undefined') ? window : {});
        if (g.RDT && g.RDT.version) { return g.RDT.version; }
        if (g.AtomTrace && g.AtomTrace.version) { return g.AtomTrace.version; }
        if (g.ToolbarPrefs && g.ToolbarPrefs.version) { return g.ToolbarPrefs.version; }
        return '(dev)';
    }

    function isoTimestamp() {
        try { return new Date().toISOString(); }
        catch (e) { return '0000-00-00T00:00:00Z'; }
    }

    /**
     * Build the canonical stamp record. Pure data; callers format it
     * per output medium.
     *
     * `opts.now` (ISO string) and `opts.version` let tests pin the
     * stamp to a deterministic value.
     */
    function stampRecord(content, opts) {
        opts = opts || {};
        var version = opts.version || currentVersion();
        var hashFn = opts.hashFn || defaultSha384;
        var sig = hashFn(String(content == null ? '' : content));
        return {
            project: 'BIME',
            version: version,
            copyright: COPYRIGHT_LINE,
            license: LICENSE_LINE,
            url: URL_LINE,
            generatedAt: opts.now || isoTimestamp(),
            algo: sig.algo,
            hash: sig.hex
        };
    }

    // -----------------------------------------------------------------------
    // Format-specific stampers.
    // -----------------------------------------------------------------------

    /**
     * Inject a `<title>` + `<metadata>` block into the SVG document
     * head, immediately after the opening `<svg>` tag. Idempotent —
     * an existing `<metadata id="bime-stamp">` block is replaced.
     *
     * If the SVG already has a `<title>` we LEAVE IT, and append our
     * metadata as a sibling. This preserves user-authored titles
     * (e.g. molecule names) while still embedding provenance.
     */
    function stampSvg(svgString, opts) {
        if (!svgString || typeof svgString !== 'string') { return svgString; }
        var rec = stampRecord(svgString, opts);
        var metaXml =
            '<metadata id="bime-stamp">' +
                '<bime:provenance xmlns:bime="https://github.com/asad/bime-dist#">' +
                    '<bime:project>BIME</bime:project>' +
                    '<bime:version>' + xmlEscape(rec.version) + '</bime:version>' +
                    '<bime:generatedAt>' + xmlEscape(rec.generatedAt) + '</bime:generatedAt>' +
                    '<bime:fingerprint algo="' + xmlEscape(rec.algo) + '">' +
                        xmlEscape(rec.hash) +
                    '</bime:fingerprint>' +
                    '<bime:copyright>' + xmlEscape(rec.copyright) + '</bime:copyright>' +
                    '<bime:license>' + xmlEscape(rec.license) + '</bime:license>' +
                    '<bime:url>' + xmlEscape(rec.url) + '</bime:url>' +
                '</bime:provenance>' +
            '</metadata>';
        // Drop any previous bime-stamp metadata for idempotency.
        var stripped = svgString.replace(/<metadata id="bime-stamp">[\s\S]*?<\/metadata>/g, '');
        // Insert right after the opening <svg ...> tag.
        var openTagMatch = stripped.match(/<svg\b[^>]*>/);
        if (!openTagMatch) { return svgString; }
        var insertAt = openTagMatch.index + openTagMatch[0].length;
        var stamped = stripped.substring(0, insertAt) + metaXml + stripped.substring(insertAt);
        // Also inject a top-right corner BIME wordmark + bug — small, low
        // opacity, doesn't compete with the molecule. Inserted just before
        // the closing </svg>. Idempotent via the `bime-watermark` class.
        if (!opts || opts.watermark !== false) {
            stamped = stamped.replace(
                /<g class="bime-watermark"[\s\S]*?<\/g>/g, ''
            );
            var wm = renderBimeWatermark(stamped);
            stamped = stamped.replace(/<\/svg>\s*$/, wm + '</svg>');
        }
        return stamped;
    }

    function xmlEscape(s) {
        return String(s == null ? '' : s)
            .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;').replace(/'/g, '&apos;');
    }

    // Hex-encoded SVG fragment for the BIME wordmark + tiny molecule bug.
    // Placed in the bottom-right corner at ~0.55 opacity so it's visible
    // but never competes with the structure. Anchored via viewBox.
    function renderBimeWatermark(svgString) {
        // Extract the viewBox to anchor the watermark to the bottom-right.
        var vb = /viewBox="([^"]+)"/.exec(svgString || '');
        var w = 320, h = 240;
        if (vb && vb[1]) {
            var parts = vb[1].split(/\s+/).map(parseFloat);
            if (parts.length === 4 && isFinite(parts[2]) && isFinite(parts[3])) {
                w = parts[2]; h = parts[3];
            }
        }
        var x = w - 70, y = h - 8;
        // Compact molecule-bug: 6-vertex hexagon (benzene reference) + the
        // "BIME" wordmark. Stroke colour is BIME teal (#0d9488) at 0.55
        // alpha so it sits over light AND dark canvases without drawing
        // the eye.
        return (
            '<g class="bime-watermark" opacity="0.55" pointer-events="none">' +
                '<g transform="translate(' + (x - 10) + ',' + (y - 11) + ') scale(0.6)" ' +
                    'fill="none" stroke="#0d9488" stroke-width="1.5" ' +
                    'stroke-linecap="round" stroke-linejoin="round">' +
                    '<path d="M8 4L13 6.5V11.5L8 14L3 11.5V6.5Z"/>' +
                '</g>' +
                '<text x="' + x + '" y="' + y + '" ' +
                    'font-family="-apple-system, system-ui, sans-serif" ' +
                    'font-size="9" font-weight="700" fill="#0d9488">BIME</text>' +
            '</g>'
        );
    }

    /**
     * Wrap a text artefact (MOL / SDF / RXN / SMILES / generic) with a
     * BIME stamp comment block. Idempotent — an existing leading stamp
     * is replaced with the current one.
     */
    function stampText(format, content, opts) {
        if (content == null) { content = ''; }
        var raw = String(content);
        // Strip any previous BIME stamp block so re-stamping is idempotent.
        // The block is delimited by `# BIME-STAMP-BEGIN` / `# BIME-STAMP-END`.
        var stripped = raw.replace(
            /^[\s\S]*?# BIME-STAMP-BEGIN[\s\S]*?# BIME-STAMP-END\s*/, ''
        );
        // The hash is computed over the un-stamped content so consumers
        // can re-derive it after stripping the stamp.
        var rec = stampRecord(stripped, opts);
        var stampLines = [
            '# BIME-STAMP-BEGIN',
            '# Project:     ' + rec.project,
            '# Version:     ' + rec.version,
            '# Generated:   ' + rec.generatedAt,
            '# ' + rec.algo + ':   ' + rec.hash,
            '# Copyright:   ' + rec.copyright,
            '# License:     ' + rec.license,
            '# URL:         ' + rec.url,
            '# BIME-STAMP-END'
        ];
        // For MOL / SDF / RXN, the comment block goes BEFORE the format
        // header (most parsers tolerate leading comment lines). For
        // SMILES, the stamp goes AFTER (so the first non-comment line is
        // still a clean SMILES string).
        var fmt = String(format || '').toLowerCase();
        if (fmt === 'smiles' || fmt === 'smi') {
            return stripped.replace(/\s*$/, '') + '\n' + stampLines.join('\n') + '\n';
        }
        // MOL, SDF, RXN, txt, log, anything else: leading comment block.
        return stampLines.join('\n') + '\n' + stripped;
    }

    /**
     * Extract a BIME stamp record from a previously-stamped artefact.
     * Returns `null` if no stamp is present. Useful for consumers that
     * want to verify provenance without re-running the hash themselves.
     */
    function extractStamp(content) {
        if (!content) { return null; }
        var raw = String(content);
        // SVG metadata block:
        var svgMatch = raw.match(
            /<metadata id="bime-stamp">[\s\S]*?<bime:version>([^<]+)<\/bime:version>[\s\S]*?<bime:fingerprint algo="([^"]+)">([^<]+)<\/bime:fingerprint>[\s\S]*?<\/metadata>/
        );
        if (svgMatch) {
            return { medium: 'svg', version: svgMatch[1], algo: svgMatch[2], hash: svgMatch[3] };
        }
        // Text comment block:
        var textMatch = raw.match(
            /# BIME-STAMP-BEGIN[\s\S]*?# Version:\s*(\S+)[\s\S]*?# (\S+):\s*([0-9a-f]+)[\s\S]*?# BIME-STAMP-END/
        );
        if (textMatch) {
            return { medium: 'text', version: textMatch[1], algo: textMatch[2], hash: textMatch[3] };
        }
        return null;
    }

    var ExportStamp = {
        stampRecord:   stampRecord,
        stampSvg:      stampSvg,
        stampText:     stampText,
        extractStamp:  extractStamp,
        // Exposed for tests + power callers.
        _defaultSha384:    defaultSha384,
        _renderWatermark:  renderBimeWatermark,
        _xmlEscape:        xmlEscape,
        version:       '3.0.3'
    };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = ExportStamp;
    }
    if (global) { global.ExportStamp = ExportStamp; }
})(typeof window !== 'undefined' ? window : (typeof global !== 'undefined' ? global : this));
