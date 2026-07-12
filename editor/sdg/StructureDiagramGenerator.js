/**
 * editor/sdg/StructureDiagramGenerator.js — SDG orchestrator.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * Entry point for callers that expect the
 * `StructureDiagramGenerator.generateCoordinates(mol)` API. BIME's own
 * code uses `Layout.layout(mol)` directly; Layout.layout is the canonical
 * pipeline and already runs every built-in SDG step:
 *
 *   1–8, 10  ring/chain/substituent placement (Layout.layout)
 *      11    E/Z correction         (SDG.CorrectGeometricConfiguration)
 *      12    wedge/dash assignment  (SDG.NonplanarBonds)
 *      13    refinement             (SDG.LayoutRefiner)
 *      14    final overlap resolve  (SDG.OverlapResolver)
 *
 * This shim adds step 9 (explicit H placement via SDG.HydrogenPlacer),
 * which Layout.layout does not run by default, so API callers get
 * a fully decorated 2D depiction in one call.
 */
(function (global) {
    'use strict';

    var SDG_NS = global.SDG = global.SDG || {};

    var StructureDiagramGenerator = {};

    StructureDiagramGenerator.generateCoordinates = function (mol, options) {
        if (!mol || !mol.atoms || mol.atoms.length === 0) return;
        options = options || {};

        // Steps 1–8, 10–14 — the full Layout pipeline. SDG.LayoutRefiner,
        // SDG.OverlapResolver, SDG.CorrectGeometricConfiguration, and
        // SDG.NonplanarBonds are all called from inside Layout.layout
        // (see editor/Layout.js step 16); calling them again here would
        // be a redundant second pass.
        if (typeof Layout !== 'undefined' && Layout.layout) {
            Layout.layout(mol);
        }

        // Step 9 — explicit H placement (Layout.layout skips this by
        // default; this API includes it for a fully decorated depiction).
        if (SDG_NS.HydrogenPlacer && options.placeHydrogens !== false) {
            SDG_NS.HydrogenPlacer.placeHydrogens2D(mol);
        }
    };

    // Convenience aliases for familiar overload-style usage.
    StructureDiagramGenerator.setMolecule = function (mol) { this._mol = mol; };
    StructureDiagramGenerator.getMolecule = function () { return this._mol; };
    StructureDiagramGenerator.layout = StructureDiagramGenerator.generateCoordinates;

    SDG_NS.StructureDiagramGenerator = StructureDiagramGenerator;
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = StructureDiagramGenerator;
    }
})(typeof globalThis !== 'undefined' ? globalThis :
   typeof window !== 'undefined' ? window : this);
