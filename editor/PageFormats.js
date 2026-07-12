/**
 * editor/PageFormats.js — print page-format sizes for the pathway canvas
 * (v2.0.38).
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Closes R2 (A4/Letter auto-scaling canvas). All sizes are in CSS pixels
 * (1 in = 96 px). The dimensions match the standard paper-size
 * definitions, NOT BIME's existing 1200×620 historical viewBox — so any
 * pathway authored for a specific page format renders to exact print
 * dimensions when exported to SVG / PDF.
 *
 *   a4 portrait  =  794 × 1123  (210 × 297 mm)
 *   a4 landscape = 1123 × 794
 *   letter portrait  =  816 × 1056  (8.5 × 11 in)
 *   letter landscape = 1056 × 816
 *
 * Plus the legacy default 'screen' = 1200 × 620 which preserves the
 * pre-v2.0.38 viewBox so existing exports are byte-identical when the
 * user has not changed the format toggle.
 *
 * No external dependencies; pure ES5; no DOM access.
 */
(function (global) {
    'use strict';

    var FORMATS = {
        'screen': {
            label: 'Screen (legacy)',
            portrait:  { w: 1200, h:  620 },
            landscape: { w: 1200, h:  620 }
        },
        'a4': {
            label: 'A4',
            portrait:  { w:  794, h: 1123 },
            landscape: { w: 1123, h:  794 }
        },
        'letter': {
            label: 'US Letter',
            portrait:  { w:  816, h: 1056 },
            landscape: { w: 1056, h:  816 }
        }
    };

    var DEFAULTS = { format: 'screen', orientation: 'landscape' };

    function sizeFor(format, orientation) {
        var f = FORMATS[format] || FORMATS[DEFAULTS.format];
        var o = (orientation === 'portrait' || orientation === 'landscape') ? orientation : DEFAULTS.orientation;
        var s = f[o] || f.landscape;
        return { w: s.w, h: s.h };
    }

    function viewBoxFor(format, orientation) {
        var s = sizeFor(format, orientation);
        return '0 0 ' + s.w + ' ' + s.h;
    }

    function formatList() {
        var out = [];
        for (var key in FORMATS) {
            if (Object.prototype.hasOwnProperty.call(FORMATS, key)) {
                out.push({ key: key, label: FORMATS[key].label });
            }
        }
        return out;
    }

    function isKnown(format, orientation) {
        if (!FORMATS[format]) { return false; }
        return orientation === 'portrait' || orientation === 'landscape';
    }

    var PageFormats = {
        FORMATS: FORMATS,
        DEFAULTS: DEFAULTS,
        sizeFor: sizeFor,
        viewBoxFor: viewBoxFor,
        formatList: formatList,
        isKnown: isKnown,
        version: '2.0.38'
    };

    global.PageFormats = PageFormats;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = PageFormats;
    }

})(typeof window !== 'undefined' ? window : this);
