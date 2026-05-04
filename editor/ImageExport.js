/**
 * ImageExport.js -- Publication-quality image export for BIME molecular structures
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 -- see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Provides SVG, PNG, print-ready SVG, batch export, clipboard copy,
 * and download helpers for molecular structures rendered by BIME.
 */
(function(global) {
    'use strict';

    var ELEMENTS = Molecule.ELEMENTS;
    var BOND_LENGTH = Molecule.BOND_LENGTH;

    // ---------------------------------------------------------------------------
    // Default rendering constants (screen)
    // ---------------------------------------------------------------------------
    var DEFAULTS = {
        width: 500,
        height: 400,
        padding: 40,
        background: '#ffffff',
        fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
        fontSize: 13,
        bondWidth: 1.8,
        doubleBondGap: 3,
        tripleBondGap: 4,
        wedgeWidth: 5,
        labelPad: 4,
        showHydrogens: true,
        showAtomNumbers: false,
        showAromaticCircles: true,
        watermark: false,
        title: ''
    };

    // Print overrides
    var PRINT = {
        fontSize: 16,
        bondWidth: 2.5
    };

    // CMYK-safe CPK colours for print
    var PRINT_COLORS = {
        'C': '#000000', 'N': '#0000cc', 'O': '#cc0000', 'S': '#ccaa00',
        'P': '#cc6600', 'F': '#009900', 'Cl': '#009900', 'Br': '#660099',
        'I': '#660099', 'H': '#000000'
    };

    // ---------------------------------------------------------------------------
    // True text measurement via hidden SVG element
    //
    // FIX: previously the measurement <svg> was lazily attached to document.body
    // on first use and never detached, leaking a DOM node for the lifetime of
    // the page (and held a strong ref to <text> children created during builds).
    // We now attach it for the duration of a build via _beginMeasureSession() /
    // _endMeasureSession() so it leaves no residue on document.body when idle.
    // ---------------------------------------------------------------------------
    var _measureSVG = null;
    var _measureDepth = 0; // re-entrant session counter

    function _beginMeasureSession() {
        if (typeof document === 'undefined') return;
        _measureDepth++;
        if (!_measureSVG) {
            _measureSVG = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
            _measureSVG.setAttribute('style', 'position:absolute;visibility:hidden;width:0;height:0;pointer-events:none');
        }
        if (!_measureSVG.parentNode) {
            document.body.appendChild(_measureSVG);
        }
    }

    function _endMeasureSession() {
        if (typeof document === 'undefined') return;
        if (_measureDepth > 0) _measureDepth--;
        if (_measureDepth === 0 && _measureSVG && _measureSVG.parentNode) {
            _measureSVG.parentNode.removeChild(_measureSVG);
        }
    }

    function measureText(text, fSize, fFamily) {
        if (typeof document === 'undefined') return text.length * fSize * 0.6;
        // If no session is open (defensive: ad-hoc measurement) attach lazily
        // and detach immediately so we still don't leak.
        var ownSession = false;
        if (!_measureSVG || !_measureSVG.parentNode) {
            _beginMeasureSession();
            ownSession = true;
        }
        var t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        t.setAttribute('font-size', fSize);
        t.setAttribute('font-family', fFamily || 'sans-serif');
        t.textContent = text;
        _measureSVG.appendChild(t);
        var w = t.getComputedTextLength();
        _measureSVG.removeChild(t);
        if (ownSession) _endMeasureSession();
        return w;
    }

    // ---------------------------------------------------------------------------
    // ImageExport namespace
    // ---------------------------------------------------------------------------
    var ImageExport = {};

    // ---------------------------------------------------------------------------
    // 1. SVG Export (publication quality)
    // ---------------------------------------------------------------------------

    /**
     * Render a Molecule to a self-contained SVG string.
     * @param {Molecule} mol
     * @param {Object} [options]
     * @returns {string} standalone SVG markup
     */
    ImageExport.toSVG = function(mol, options) {
        var opts = _mergeOpts(DEFAULTS, options);
        return _buildSVG(mol, opts);
    };

    // ---------------------------------------------------------------------------
    // 2. PNG Export (high resolution)
    // ---------------------------------------------------------------------------

    /**
     * Render a Molecule to a PNG Blob via offscreen canvas.
     * @param {Molecule} mol
     * @param {Object} [options]  scale: 1|2|4, width, height, background
     * @returns {Promise<Blob>}
     */
    ImageExport.toPNG = function(mol, options) {
        var opts = _mergeOpts(DEFAULTS, options);
        var scale = Math.max(1, Math.min(10, opts.scale || 2));
        var svgStr = _buildSVG(mol, opts);

        var w = opts.width;
        var h = opts.height;
        // Enforce minimum 1024px at 4x for publication quality
        if (scale === 4 && w * scale < 1024) {
            var factor = 1024 / (w * scale);
            w = Math.ceil(w * factor);
            h = Math.ceil(h * factor);
            svgStr = _buildSVG(mol, _mergeOpts(opts, { width: w, height: h }));
        }

        return new Promise(function(resolve, reject) {
            var canvas = document.createElement('canvas');
            canvas.width = w * scale;
            canvas.height = h * scale;
            var ctx = canvas.getContext('2d');
            ctx.imageSmoothingEnabled = true;
            ctx.imageSmoothingQuality = 'high';

            // Background
            if (opts.background && opts.background !== 'transparent') {
                ctx.fillStyle = opts.background;
                ctx.fillRect(0, 0, canvas.width, canvas.height);
            }

            var img = new Image();
            var blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
            var url = URL.createObjectURL(blob);

            img.onload = function() {
                ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
                URL.revokeObjectURL(url);
                canvas.toBlob(function(pngBlob) {
                    if (pngBlob) resolve(pngBlob);
                    else reject(new Error('PNG encoding failed'));
                }, 'image/png');
            };
            img.onerror = function() {
                URL.revokeObjectURL(url);
                reject(new Error('SVG rasterisation failed'));
            };
            img.src = url;
        });
    };

    // ---------------------------------------------------------------------------
    // 3. PDF-ready SVG (print optimised)
    // ---------------------------------------------------------------------------

    /**
     * Return SVG optimised for print/PDF embedding.
     * @param {Molecule} mol
     * @param {Object} [options]
     * @returns {string} SVG string
     */
    ImageExport.toPrintSVG = function(mol, options) {
        var printOpts = _mergeOpts(DEFAULTS, {
            fontSize: PRINT.fontSize,
            bondWidth: PRINT.bondWidth,
            background: '#ffffff',
            watermark: false,
            printMode: true
        });
        printOpts = _mergeOpts(printOpts, options);
        return _buildSVG(mol, printOpts);
    };

    // ---------------------------------------------------------------------------
    // 3b. PDF Export (via browser print dialog)
    // ---------------------------------------------------------------------------

    /**
     * Open a print dialog for "Save as PDF" -- works in all browsers
     * without external dependencies.
     * @param {Molecule} mol
     * @param {Object} [options]
     */
    ImageExport.toPDF = function(mol, options) {
        var svg = ImageExport.toPrintSVG(mol, options);
        // Build a complete HTML document and open via a Blob URL.
        // Avoids document.write (deprecated, triggers cross-origin warnings, and was
        // flagged in the BIME security audit as defense-in-depth regression).
        var html =
            '<!DOCTYPE html><html><head><meta charset="utf-8"><title>BIME Export</title>' +
            '<style>@page{margin:0.5in}html,body{margin:0;padding:0}' +
            'body{display:flex;justify-content:center;align-items:center;min-height:100vh}</style>' +
            '</head><body>' + svg + '</body></html>';
        var blob = new Blob([html], { type: 'text/html;charset=utf-8' });
        var url = URL.createObjectURL(blob);
        var win = window.open(url, '_blank');
        if (!win) { URL.revokeObjectURL(url); return; }
        // Print once the new window has loaded the Blob document.
        var printed = false;
        var doPrint = function() {
            if (printed) return;
            printed = true;
            try { win.print(); } catch (e) { /* ignore */ }
            // Revoke after a delay so the print preview can finish reading the URL.
            setTimeout(function() { URL.revokeObjectURL(url); }, 60000);
        };
        if (win.document && win.document.readyState === 'complete') {
            doPrint();
        } else {
            win.addEventListener('load', doPrint);
            // Fallback in case the load event fires before the listener attaches.
            setTimeout(doPrint, 1500);
        }
    };

    // ---------------------------------------------------------------------------
    // 4. Batch Export
    // ---------------------------------------------------------------------------

    /**
     * Batch-convert an array of {name, smiles} to SVG strings.
     * @param {Array<{name:string, smiles:string}>} molecules
     * @param {Object} [options]
     * @returns {Array<string>}
     */
    ImageExport.batchSVG = function(molecules, options) {
        var results = [];
        for (var i = 0; i < molecules.length; i++) {
            var entry = molecules[i];
            var mol = _parseSMILES(entry.smiles);
            if (!mol) { results.push(''); continue; }
            var opts = _mergeOpts(DEFAULTS, options);
            opts.title = entry.name || '';
            results.push(_buildSVG(mol, opts));
        }
        return results;
    };

    /**
     * Batch-convert an array of {name, smiles} to PNG Blobs.
     * @param {Array<{name:string, smiles:string}>} molecules
     * @param {Object} [options]
     * @returns {Promise<Array<Blob>>}
     */
    ImageExport.batchPNG = function(molecules, options) {
        var promises = [];
        for (var i = 0; i < molecules.length; i++) {
            var entry = molecules[i];
            var mol = _parseSMILES(entry.smiles);
            if (!mol) {
                promises.push(Promise.resolve(null));
                continue;
            }
            var opts = _mergeOpts(DEFAULTS, options);
            opts.title = entry.name || '';
            promises.push(ImageExport.toPNG(mol, opts));
        }
        return Promise.all(promises);
    };

    // ---------------------------------------------------------------------------
    // 5. Download Helpers
    // ---------------------------------------------------------------------------

    /**
     * Trigger browser download of an SVG string.
     * @param {string} svg
     * @param {string} [filename]
     */
    ImageExport.downloadSVG = function(svg, filename) {
        filename = filename || 'molecule.svg';
        _downloadBlob(new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }), filename);
    };

    /**
     * Trigger browser download of a PNG Blob.
     * @param {Blob} blob
     * @param {string} [filename]
     */
    ImageExport.downloadPNG = function(blob, filename) {
        filename = filename || 'molecule.png';
        _downloadBlob(blob, filename);
    };

    /**
     * Download all molecules as images (ZIP if JSZip available, else sequential).
     * @param {Array<{name:string, smiles:string}>} molecules
     * @param {string} [format]  'svg' or 'png'
     * @param {Object} [options]
     * @returns {Promise<void>}
     */
    ImageExport.downloadAll = function(molecules, format, options) {
        format = (format || 'svg').toLowerCase();

        if (format === 'svg') {
            var svgs = ImageExport.batchSVG(molecules, options);
            if (typeof JSZip !== 'undefined') {
                return _zipDownload(molecules, svgs, null, 'svg');
            }
            // Sequential download fallback
            for (var i = 0; i < svgs.length; i++) {
                if (svgs[i]) {
                    var name = _safeName(molecules[i].name || ('molecule_' + (i + 1)));
                    ImageExport.downloadSVG(svgs[i], name + '.svg');
                }
            }
            return Promise.resolve();
        }

        // PNG
        return ImageExport.batchPNG(molecules, options).then(function(blobs) {
            if (typeof JSZip !== 'undefined') {
                return _zipDownload(molecules, null, blobs, 'png');
            }
            for (var i = 0; i < blobs.length; i++) {
                if (blobs[i]) {
                    var name = _safeName(molecules[i].name || ('molecule_' + (i + 1)));
                    ImageExport.downloadPNG(blobs[i], name + '.png');
                }
            }
        });
    };

    // ---------------------------------------------------------------------------
    // 6. Clipboard
    // ---------------------------------------------------------------------------

    /**
     * Copy molecule as PNG to clipboard.
     * @param {Molecule} mol
     * @param {Object} [options]
     * @param {HTMLElement} [statusBar] optional element to show "Copied!" confirmation
     * @returns {Promise<void>}
     */
    ImageExport.copyToClipboard = function(mol, options, statusBar) {
        return ImageExport.toPNG(mol, _mergeOpts({ scale: 2 }, options)).then(function(blob) {
            if (!navigator.clipboard || !navigator.clipboard.write) {
                throw new Error('Clipboard API not available');
            }
            // Include SMILES (with name if present) as plain text alongside the PNG
            var smilesText = '';
            if (typeof SmilesWriter !== 'undefined') {
                var rawSmiles = SmilesWriter.write ? SmilesWriter.write(mol) : '';
                // Strip any name already appended by SmilesWriter.write, then re-add cleanly
                var baseSmiles = rawSmiles.split(/[ \t]/)[0];
                smilesText = baseSmiles + (mol.name ? ' ' + mol.name : '');
            }
            var clipItems = { 'image/png': blob };
            if (smilesText && typeof ClipboardItem !== 'undefined') {
                clipItems['text/plain'] = new Blob([smilesText], { type: 'text/plain' });
            }
            var item = new ClipboardItem(clipItems);
            return navigator.clipboard.write([item]);
        }).then(function() {
            if (statusBar) {
                var info = statusBar.querySelector('#bime-info');
                if (info) {
                    info.textContent = 'Copied!';
                    setTimeout(function() { info.textContent = ''; }, 2000);
                }
            }
        });
    };

    /**
     * Copy molecule as MOL V2000 text to clipboard.
     * @param {Molecule} mol
     * @returns {Promise<void>}
     */
    ImageExport.copyMOL = function(mol) {
        if (typeof MolfileWriter === 'undefined') {
            return Promise.reject(new Error('MolfileWriter not available'));
        }
        // FIX: MolfileWriter exposes write() / writeSDF(); there is no toMolfile().
        // The previous typo caused copyMOL to throw TypeError.
        var molText = MolfileWriter.write(mol);
        return navigator.clipboard.writeText(molText);
    };

    // =========================================================================
    // Internal: SVG builder
    // =========================================================================

    /**
     * Build a complete, self-contained SVG string for a molecule.
     */
    function _buildSVG(mol, opts) {
        if (!mol || mol.isEmpty()) {
            return _emptySVG(opts);
        }
        _beginMeasureSession();
        try {
            return _buildSVGImpl(mol, opts);
        } finally {
            _endMeasureSession();
        }
    }

    function _buildSVGImpl(mol, opts) {

        var fontSize = opts.fontSize || DEFAULTS.fontSize;
        var bondWidth = opts.bondWidth || DEFAULTS.bondWidth;
        var padding = opts.padding !== undefined ? opts.padding : DEFAULTS.padding;
        var doubleBondGap = opts.doubleBondGap || DEFAULTS.doubleBondGap;
        var tripleBondGap = opts.tripleBondGap || DEFAULTS.tripleBondGap;
        var wedgeWidth = opts.wedgeWidth || DEFAULTS.wedgeWidth;
        var labelPad = opts.labelPad || DEFAULTS.labelPad;
        var showH = opts.showHydrogens !== false;
        var showAtomNums = !!opts.showAtomNumbers;
        var showAromCircles = opts.showAromaticCircles !== false;
        var printMode = !!opts.printMode;
        var aromaticStyle = opts.aromaticStyle || 'dashed';

        // Print mode uses journal-standard font stack.
        var fontFamily = printMode
            ? '"Helvetica Neue", Helvetica, Arial, sans-serif'
            : (opts.fontFamily || DEFAULTS.fontFamily);

        // v1.8.6: emit fontFamily into SVG attributes safely. The default
        // stack contains "Segoe UI" with literal double quotes, which —
        // when interpolated into a double-quoted attribute (font-family="..."
        // or style="...") — produces malformed XML that strict parsers
        // reject (e.g. browsers loading via <object>, <img>, or copy-paste
        // into another tool). The bug manifested as "broken" depicts on
        // search-result hover tooltips for any molecule, because the
        // tooltip's <object>-style consumption tripped the strict path.
        // Fix: replace inner " with &quot; for SVG-attribute use; the
        // raw value is still emitted into <style> elements where the
        // CSS parser tolerates it. See also: docs/USAGE.md "SVG export".
        var fontFamilyAttr = String(fontFamily).replace(/"/g, '&quot;');

        // Compute molecule bounds and auto-scale
        var bounds = mol.getBounds();
        var molW = bounds.w || 1;
        var molH = bounds.h || 1;

        var targetW = opts.width || DEFAULTS.width;
        var targetH = opts.height || DEFAULTS.height;
        var innerW = targetW - padding * 2;
        var innerH = targetH - padding * 2;

        var scaleX = innerW / molW;
        var scaleY = innerH / molH;
        var scale;
        if (opts.fixedScale) {
            // Fixed scale: pixels per bond length (e.g., 30 = 30px per bond)
            scale = opts.fixedScale / BOND_LENGTH;
        } else {
            scale = Math.min(scaleX, scaleY, 3);
        }
        scale = Math.max(scale, 0.2);

        var offsetX = (targetW / scale - molW) / 2 - bounds.x;
        var offsetY = (targetH / scale - molH) / 2 - bounds.y;

        function tx(x) { return (x + offsetX) * scale; }
        function ty(y) { return (y + offsetY) * scale; }

        function atomColor(symbol) {
            if (printMode && PRINT_COLORS[symbol]) return PRINT_COLORS[symbol];
            var elem = ELEMENTS[symbol];
            return elem ? elem.color : '#333333';
        }

        function estimateTextWidth(text) {
            return measureText(text, fontSize, fontFamily);
        }

        function labelRadius(atom) {
            if (atom.symbol === 'C' && atom.charge === 0 && atom.isotope === 0 && mol.degree(atom.id) > 0) return 0;
            return estimateTextWidth(atom.symbol) / 2 + labelPad;
        }

        // Collect ring info for double bond offset and aromatic circles
        var ringInfo = (typeof Layout !== 'undefined' && Layout.getRingInfo) ? Layout.getRingInfo(mol) : [];

        function findRingForBond(bond) {
            for (var i = 0; i < ringInfo.length; i++) {
                var ring = ringInfo[i];
                var atoms = ring.atoms;
                var idx1 = atoms.indexOf(bond.atom1);
                var idx2 = atoms.indexOf(bond.atom2);
                if (idx1 >= 0 && idx2 >= 0) {
                    var diff = Math.abs(idx1 - idx2);
                    if (diff === 1 || diff === atoms.length - 1) return ring;
                }
            }
            return null;
        }

        // Hydrogen placement direction
        function hydrogenDirection(atom) {
            var neighbors = mol.getNeighbors(atom.id);
            if (neighbors.length === 0) return 'right';
            var sumX = 0, sumY = 0;
            for (var i = 0; i < neighbors.length; i++) {
                var n = mol.getAtom(neighbors[i]);
                if (!n) continue;
                var dx = n.x - atom.x;
                var dy = n.y - atom.y;
                var len = Math.sqrt(dx * dx + dy * dy);
                if (len > 0) { sumX += dx / len; sumY += dy / len; }
            }
            var angle = Math.atan2(-sumY, -sumX);
            if (angle > -Math.PI / 4 && angle <= Math.PI / 4) return 'right';
            if (angle > Math.PI / 4 && angle <= 3 * Math.PI / 4) return 'down';
            if (angle > -3 * Math.PI / 4 && angle <= -Math.PI / 4) return 'up';
            return 'left';
        }

        var parts = [];

        // --- Aromatic circles ---
        if (showAromCircles) {
            for (var ri = 0; ri < ringInfo.length; ri++) {
                var ring = ringInfo[ri];
                if (!ring.aromatic) continue;
                var rcx = tx(ring.center.x);
                var rcy = ty(ring.center.y);
                // FIX: guard against stale ringInfo referencing a removed atom
                var totalDist = 0;
                var counted = 0;
                for (var ai = 0; ai < ring.atoms.length; ai++) {
                    var ra = mol.getAtom(ring.atoms[ai]);
                    if (!ra) continue;
                    var rdx = tx(ra.x) - rcx;
                    var rdy = ty(ra.y) - rcy;
                    totalDist += Math.sqrt(rdx * rdx + rdy * rdy);
                    counted++;
                }
                if (counted === 0) continue;
                var avgR = totalDist / counted;
                var innerR = avgR * 0.6;
                var circColor = printMode ? '#666666' : '#999999';
                if (aromaticStyle !== 'none') {
                    var dashAttr = aromaticStyle === 'solid' ? '' : ' stroke-dasharray="3,2"';
                    parts.push('<circle cx="' + _r(rcx) + '" cy="' + _r(rcy) + '" r="' + _r(innerR) +
                        '" fill="none" stroke="' + circColor + '" stroke-width="1"' + dashAttr + '/>');
                }
            }
        }

        // --- Bonds ---
        for (var bi = 0; bi < mol.bonds.length; bi++) {
            var bond = mol.bonds[bi];
            var a1 = mol.getAtom(bond.atom1);
            var a2 = mol.getAtom(bond.atom2);
            if (!a1 || !a2) continue;

            var x1 = tx(a1.x), y1 = ty(a1.y);
            var x2 = tx(a2.x), y2 = ty(a2.y);

            var dx = x2 - x1, dy = y2 - y1;
            var len = Math.sqrt(dx * dx + dy * dy);
            if (len < 1) continue;
            var nx = dx / len, ny = dy / len;

            // Trim bonds at atom label boundaries
            var trim1 = labelRadius(a1);
            var trim2 = labelRadius(a2);
            x1 += nx * trim1; y1 += ny * trim1;
            x2 -= nx * trim2; y2 -= ny * trim2;

            var bColor = printMode ? '#000000' : '#333333';

            if (bond.stereo === Molecule.STEREO_WEDGE) {
                // Wedge bond (filled triangle)
                var wdx = x2 - x1, wdy = y2 - y1;
                var wlen = Math.sqrt(wdx * wdx + wdy * wdy);
                if (wlen < 0.5) continue;
                var wnx = -wdy / wlen, wny = wdx / wlen;
                var w = wedgeWidth;
                parts.push('<polygon points="' +
                    _r(x1) + ',' + _r(y1) + ' ' +
                    _r(x2 + wnx * w) + ',' + _r(y2 + wny * w) + ' ' +
                    _r(x2 - wnx * w) + ',' + _r(y2 - wny * w) +
                    '" fill="' + bColor + '" stroke="' + bColor + '" stroke-width="0.5" stroke-linejoin="round"/>');

            } else if (bond.stereo === Molecule.STEREO_DASH) {
                // Dashed wedge bond
                var ddx = x2 - x1, ddy = y2 - y1;
                var dlen = Math.sqrt(ddx * ddx + ddy * ddy);
                if (dlen < 0.5) continue;
                var dnx = -ddy / dlen, dny = ddx / dlen;
                var segments = 7;
                for (var si = 1; si <= segments; si++) {
                    var t = si / segments;
                    var cx = x1 + ddx * t, cy = y1 + ddy * t;
                    var dw = wedgeWidth * t;
                    parts.push('<line x1="' + _r(cx - dnx * dw) + '" y1="' + _r(cy - dny * dw) +
                        '" x2="' + _r(cx + dnx * dw) + '" y2="' + _r(cy + dny * dw) +
                        '" stroke="' + bColor + '" stroke-width="1.5" stroke-linecap="round"/>');
                }

            } else if (bond.type === Molecule.BOND_DOUBLE) {
                var ringForBond = findRingForBond(bond);
                if (ringForBond) {
                    // Ring double bond: main + offset
                    var rcx2 = tx(ringForBond.center.x);
                    var rcy2 = ty(ringForBond.center.y);
                    var mx = (x1 + x2) / 2, my = (y1 + y2) / 2;
                    var perpX = -ny, perpY = nx;
                    var dot = (rcx2 - mx) * perpX + (rcy2 - my) * perpY;
                    var sign = dot >= 0 ? 1 : -1;
                    var offX = perpX * doubleBondGap * 1.5 * sign;
                    var offY = perpY * doubleBondGap * 1.5 * sign;

                    parts.push(_svgLine(x1, y1, x2, y2, bColor, bondWidth));
                    var trimFrac = 0.15;
                    var bdx2 = x2 - x1, bdy2 = y2 - y1;
                    var blen2 = Math.sqrt(bdx2 * bdx2 + bdy2 * bdy2);
                    parts.push(_svgLine(
                        x1 + offX + nx * blen2 * trimFrac, y1 + offY + ny * blen2 * trimFrac,
                        x2 + offX - nx * blen2 * trimFrac, y2 + offY - ny * blen2 * trimFrac,
                        bColor, bondWidth * 0.8));
                } else {
                    // Chain double bond: symmetric
                    var px2 = -ny * doubleBondGap, py2 = nx * doubleBondGap;
                    parts.push(_svgLine(x1 + px2, y1 + py2, x2 + px2, y2 + py2, bColor, bondWidth));
                    parts.push(_svgLine(x1 - px2, y1 - py2, x2 - px2, y2 - py2, bColor, bondWidth));
                }

            } else if (bond.type === Molecule.BOND_TRIPLE) {
                var px3 = -ny * tripleBondGap, py3 = nx * tripleBondGap;
                parts.push(_svgLine(x1, y1, x2, y2, bColor, bondWidth));
                parts.push(_svgLine(x1 + px3, y1 + py3, x2 + px3, y2 + py3, bColor, bondWidth * 0.7));
                parts.push(_svgLine(x1 - px3, y1 - py3, x2 - px3, y2 - py3, bColor, bondWidth * 0.7));

            } else {
                // Single bond
                parts.push(_svgLine(x1, y1, x2, y2, bColor, bondWidth));
            }
        }

        // --- Atoms ---
        for (var ai2 = 0; ai2 < mol.atoms.length; ai2++) {
            var atom = mol.atoms[ai2];
            var ax = tx(atom.x);
            var ay = ty(atom.y);
            var elemColor = atomColor(atom.symbol);

            var showLabel = atom.symbol !== 'C' || atom.charge !== 0 ||
                atom.isotope > 0 || mol.degree(atom.id) === 0;

            if (showLabel) {
                var tw = estimateTextWidth(atom.symbol);

                // White background behind label
                parts.push('<rect x="' + _r(ax - tw / 2 - labelPad) + '" y="' + _r(ay - fontSize / 2 - 2) +
                    '" width="' + _r(tw + labelPad * 2) + '" height="' + _r(fontSize + 4) +
                    '" fill="' + opts.background + '" stroke="none"/>');

                // Atom label
                parts.push('<text x="' + _r(ax) + '" y="' + _r(ay + fontSize * 0.35) +
                    '" fill="' + elemColor + '" font-size="' + fontSize + 'px" font-family="' + fontFamilyAttr +
                    '" font-weight="bold" text-anchor="middle">' + _esc(atom.symbol) + '</text>');

                // Charge
                if (atom.charge !== 0) {
                    var chargeStr = atom.charge > 0
                        ? (atom.charge === 1 ? '+' : atom.charge + '+')
                        : (atom.charge === -1 ? '\u2212' : Math.abs(atom.charge) + '\u2212');
                    parts.push('<text x="' + _r(ax + tw / 2 + 2) + '" y="' + _r(ay - fontSize * 0.2) +
                        '" fill="' + elemColor + '" font-size="9px" font-family="' + fontFamilyAttr +
                        '" text-anchor="start">' + _esc(chargeStr) + '</text>');
                }

                // Implicit hydrogens
                if (showH) {
                    var hCount = mol.calcHydrogens(atom.id);
                    if (hCount > 0) {
                        var hStr = 'H';
                        var hCountStr = hCount > 1 ? '' + hCount : '';
                        var hDir = hydrogenDirection(atom);
                        var hFontSize = fontSize * 0.77;
                        var hTextW = estimateTextWidth(hStr) + (hCountStr ? estimateTextWidth(hCountStr) * 0.6 : 0);
                        var hx, hy, hAnchor;
                        if (hDir === 'left') {
                            hx = ax - tw / 2 - (atom.charge !== 0 ? 10 : 2) - hTextW / 2;
                            hy = ay + fontSize * 0.35;
                            hAnchor = 'end';
                        } else if (hDir === 'up') {
                            hx = ax; hy = ay - fontSize * 0.6; hAnchor = 'middle';
                        } else if (hDir === 'down') {
                            hx = ax; hy = ay + fontSize * 1.2; hAnchor = 'middle';
                        } else {
                            hx = ax + tw / 2 + (atom.charge !== 0 ? 10 : 2) + hTextW / 2;
                            hy = ay + fontSize * 0.35;
                            hAnchor = 'start';
                        }

                        // Background behind H label
                        if (hDir === 'left' || hDir === 'right') {
                            var hBgX = hDir === 'right'
                                ? ax + tw / 2 + (atom.charge !== 0 ? 10 : 2)
                                : ax - tw / 2 - (atom.charge !== 0 ? 10 : 2) - hTextW;
                            parts.push('<rect x="' + _r(hBgX - 1) + '" y="' + _r(ay - fontSize / 2 - 2) +
                                '" width="' + _r(hTextW + 2) + '" height="' + _r(fontSize + 4) +
                                '" fill="' + opts.background + '" stroke="none"/>');
                        }

                        var hLabel = hStr;
                        if (hCountStr) {
                            // Use tspan for subscript count
                            hLabel = 'H<tspan font-size="' + Math.round(hFontSize * 0.75) +
                                'px" dy="3">' + hCountStr + '</tspan>';
                        }
                        parts.push('<text x="' + _r(hx) + '" y="' + _r(hy) +
                            '" fill="#666666" font-size="' + _r(hFontSize) + 'px" font-family="' + fontFamilyAttr +
                            '" text-anchor="' + hAnchor + '">' + hLabel + '</text>');
                    }
                }

                // Isotope
                if (atom.isotope > 0) {
                    parts.push('<text x="' + _r(ax - tw / 2 - 4) + '" y="' + _r(ay - fontSize * 0.2) +
                        '" fill="' + elemColor + '" font-size="9px" font-family="' + fontFamilyAttr +
                        '" text-anchor="end">' + atom.isotope + '</text>');
                }
            }

            // Atom numbers (optional)
            if (showAtomNums) {
                parts.push('<text x="' + _r(ax) + '" y="' + _r(ay + fontSize + 8) +
                    '" fill="#999999" font-size="8px" font-family="' + fontFamilyAttr +
                    '" text-anchor="middle">' + (ai2 + 1) + '</text>');
            }

            // Atom-atom mapping numbers
            if (atom.mapNumber > 0) {
                var mapY = ay + fontSize + 6;
                parts.push('<text x="' + _r(ax) + '" y="' + _r(mapY + 12 * 0.3) +
                    '" fill="#0d9488" font-size="12px" font-family="' + fontFamilyAttr +
                    '" font-weight="bold" text-anchor="middle">' + atom.mapNumber + '</text>');
            }
        }

        // --- Reaction arrow ---
        if (mol.reactionArrow) {
            var arrow = mol.reactionArrow;
            var ax1 = tx(arrow.x1), ay1 = ty(arrow.y1);
            var ax2 = tx(arrow.x2), ay2 = ty(arrow.y2);
            var arrowColor = printMode ? '#000000' : '#333333';
            var adx = ax2 - ax1, ady = ay2 - ay1;
            var alen = Math.sqrt(adx * adx + ady * ady);
            if (alen > 1) {
                var anx = adx / alen, any = ady / alen;
                var apx = -any, apy = anx;
                var hs = 12;
                var isRetro = (arrow.type === 'retro');

                if (isRetro) {
                    var shaftEnd = alen - hs;
                    parts.push(_svgLine(ax1, ay1, ax1 + anx * shaftEnd, ay1 + any * shaftEnd, arrowColor, 2.2));
                    var hw = hs * 0.5;
                    parts.push(_svgLine(ax2 - anx * hs + apx * hw, ay2 - any * hs + apy * hw, ax2, ay2, arrowColor, 2.2));
                    parts.push(_svgLine(ax2 - anx * hs - apx * hw, ay2 - any * hs - apy * hw, ax2, ay2, arrowColor, 2.2));
                } else {
                    var shaftEnd2 = alen - hs * 0.85;
                    parts.push(_svgLine(ax1, ay1, ax1 + anx * shaftEnd2, ay1 + any * shaftEnd2, arrowColor, 2.2));
                    var hw2 = hs * 0.45;
                    var baseX = ax2 - anx * hs, baseY = ay2 - any * hs;
                    parts.push('<polygon points="' +
                        _r(ax2) + ',' + _r(ay2) + ' ' +
                        _r(baseX + apx * hw2) + ',' + _r(baseY + apy * hw2) + ' ' +
                        _r(baseX - apx * hw2) + ',' + _r(baseY - apy * hw2) +
                        '" fill="' + arrowColor + '" stroke="' + arrowColor +
                        '" stroke-width="0.5" stroke-linejoin="round"/>');
                }

                // Conditions text
                if (arrow.conditions) {
                    var cmx = (ax1 + ax2) / 2;
                    var cmy = (ay1 + ay2) / 2 - 12;
                    parts.push('<text x="' + _r(cmx) + '" y="' + _r(cmy) +
                        '" fill="#666666" font-size="10px" font-family="' + fontFamilyAttr +
                        '" font-style="italic" text-anchor="middle">' + _esc(arrow.conditions) + '</text>');
                }
            }
        }

        // --- Plus signs ---
        if (mol.reactionPlusSigns) {
            for (var pi = 0; pi < mol.reactionPlusSigns.length; pi++) {
                var pos = mol.reactionPlusSigns[pi];
                parts.push('<text x="' + _r(tx(pos.x)) + '" y="' + _r(ty(pos.y) + 18 * 0.35) +
                    '" fill="#555555" font-size="18px" font-family="' + fontFamilyAttr +
                    '" font-weight="bold" text-anchor="middle">+</text>');
            }
        }

        // --- Watermark ---
        if (opts.watermark) {
            parts.push('<text x="' + (targetW - 8) + '" y="' + (targetH - 6) +
                '" fill="#cccccc" font-size="8px" font-family="' + fontFamilyAttr +
                '" text-anchor="end">Generated by BIME</text>');
        }

        // --- Molecule name label at the bottom ---
        var molDisplayName = mol.name || opts.title || '';
        if (molDisplayName) {
            var cx = targetW / 2;
            parts.push('<text x="' + _r(cx) + '" y="' + _r(targetH - 10) +
                '" text-anchor="middle" font-family="' + fontFamilyAttr +
                '" font-size="12" fill="#333333">' + _esc(molDisplayName) + '</text>');
        }

        // --- Assemble SVG ---
        var bgRect = '<rect width="' + targetW + '" height="' + targetH + '" fill="' + (opts.background || '#ffffff') + '"/>';

        // <title> from mol.name or opts.title
        var titleStr = mol.name || opts.title || 'Molecule';
        var metadata = '<title>' + _esc(titleStr) + '</title>';
        if (opts.title || mol.name) {
            var dcTitle = mol.name || opts.title;
            metadata += '<metadata><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"' +
                ' xmlns:dc="http://purl.org/dc/elements/1.1/">' +
                '<rdf:Description><dc:title>' + _esc(dcTitle) + '</dc:title>' +
                '<dc:creator>BIME - BioInception Molecular Editor</dc:creator>' +
                '</rdf:Description></rdf:RDF></metadata>';
        }

        var svgOpen = '<svg xmlns="http://www.w3.org/2000/svg" width="' + targetW + '" height="' + targetH +
            '" viewBox="0 0 ' + targetW + ' ' + targetH + '"' +
            ' style="font-family:' + fontFamilyAttr + '">';

        return svgOpen + metadata + bgRect + parts.join('') + '</svg>';
    }

    // =========================================================================
    // Internal helpers
    // =========================================================================

    function _emptySVG(opts) {
        var w = opts.width || DEFAULTS.width;
        var h = opts.height || DEFAULTS.height;
        var bg = opts.background || '#ffffff';
        return '<svg xmlns="http://www.w3.org/2000/svg" width="' + w + '" height="' + h +
            '" viewBox="0 0 ' + w + ' ' + h + '">' +
            '<rect width="' + w + '" height="' + h + '" fill="' + bg + '"/></svg>';
    }

    function _svgLine(x1, y1, x2, y2, color, width) {
        return '<line x1="' + _r(x1) + '" y1="' + _r(y1) + '" x2="' + _r(x2) + '" y2="' + _r(y2) +
            '" stroke="' + color + '" stroke-width="' + _r(width) + '" stroke-linecap="round"/>';
    }

    function _r(n) {
        return Math.round(n * 100) / 100;
    }

    function _esc(str) {
        return String(str).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    }

    function _mergeOpts(defaults, overrides) {
        var result = {};
        for (var k in defaults) {
            if (defaults.hasOwnProperty(k)) result[k] = defaults[k];
        }
        if (overrides) {
            for (var k2 in overrides) {
                if (overrides.hasOwnProperty(k2) && overrides[k2] !== undefined) result[k2] = overrides[k2];
            }
        }
        return result;
    }

    function _parseSMILES(smiles) {
        if (!smiles || typeof SmilesParser === 'undefined') return null;
        try {
            var mol = new Molecule();
            SmilesParser.parse(smiles.trim(), mol);
            if (mol.isEmpty()) return null;
            if (typeof Layout !== 'undefined' && Layout.layout) {
                Layout.layout(mol);
            }
            return mol;
        } catch (e) {
            return null;
        }
    }

    function _downloadBlob(blob, filename) {
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.style.display = 'none';
        document.body.appendChild(a);
        a.click();
        setTimeout(function() {
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }, 100);
    }

    function _safeName(name) {
        return name.replace(/[^a-zA-Z0-9_\-]/g, '_').substring(0, 64);
    }

    function _zipDownload(molecules, svgs, pngBlobs, ext) {
        var zip = new JSZip();
        for (var i = 0; i < molecules.length; i++) {
            var name = _safeName(molecules[i].name || ('molecule_' + (i + 1)));
            if (ext === 'svg' && svgs && svgs[i]) {
                zip.file(name + '.svg', svgs[i]);
            } else if (ext === 'png' && pngBlobs && pngBlobs[i]) {
                zip.file(name + '.png', pngBlobs[i]);
            }
        }
        return zip.generateAsync({ type: 'blob' }).then(function(content) {
            _downloadBlob(content, 'bime-molecules.' + ext + '.zip');
        });
    }

    // =========================================================================
    // Export
    // =========================================================================
    /**
     * Detach the hidden measurement <svg> from document.body if attached.
     * Useful for tests / teardown so the DOM is empty after exports complete.
     */
    ImageExport._cleanup = function() {
        if (_measureSVG && _measureSVG.parentNode) {
            _measureSVG.parentNode.removeChild(_measureSVG);
        }
        _measureDepth = 0;
    };

    global.ImageExport = ImageExport;

})(typeof window !== 'undefined' ? window : this);
