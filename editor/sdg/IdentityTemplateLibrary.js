/**
 * editor/sdg/IdentityTemplateLibrary.js — CDK IdentityTemplateLibrary port.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Reference (Java original, ~18 KB):
 *   https://github.com/cdk/cdk/blob/main/tool/sdg/src/main/java/org/openscience/cdk/layout/IdentityTemplateLibrary.java
 *
 * Stores pre-computed coordinates for known scaffolds keyed by canonical
 * SMILES. Used by TemplateHandler / StructureDiagramGenerator to skip
 * algorithmic placement when a faster "memo" lookup is available.
 *
 * Status: SCAFFOLD. BIME has its own templates in editor/Templates.js
 * (76 scaffolds in v1.5.2). v1.8.x will:
 *   1. Wrap editor/Templates.js entries as IdentityTemplate objects
 *   2. Optionally augment with the CDK distribution's template library
 *      (which ships ~3000+ scaffolds via the cdk-data jar — license-
 *      compatible since CDK is LGPL-2.1 and BIME-derived weights
 *      remain usable under Apache 2.0 as long as the original CDK
 *      distribution is not redistributed verbatim).
 */
(function (global) {
    'use strict';

    var IdentityTemplateLibrary = {};

    // ---------------------------------------------------------------------
    // store / retrieve / size — match CDK API.
    // Stub: in-memory map. v1.8.x will route through editor/Templates.js.
    // ---------------------------------------------------------------------
    var _byCanon = {};

    IdentityTemplateLibrary.store = function (canonicalSmiles, atomCoords) {
        _byCanon[canonicalSmiles] = atomCoords;
    };

    IdentityTemplateLibrary.retrieve = function (canonicalSmiles) {
        return _byCanon[canonicalSmiles] || null;
    };

    IdentityTemplateLibrary.size = function () {
        return Object.keys(_byCanon).length;
    };

    IdentityTemplateLibrary.has = function (canonicalSmiles) {
        return Object.prototype.hasOwnProperty.call(_byCanon, canonicalSmiles);
    };

    global.SDG = global.SDG || {};
    global.SDG.IdentityTemplateLibrary = IdentityTemplateLibrary;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = IdentityTemplateLibrary;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
