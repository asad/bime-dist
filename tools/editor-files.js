/**
 * Canonical browser/editor source load order for BIME.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 *
 * This list feeds source-mode pages and dist bundles. Keep ordering changes
 * here so the build, pages, and release tests cannot drift independently.
 */
'use strict';

var FILES = [
    // v2.0.33: shared helpers loaded BEFORE the modules that depend on them.
    'XmlUtil.js',
    'CurlyArrowEndpoint.js',
    'Molecule.js',
    'Layout.js',
    'Templates.js',
    'MetaboliteLibrary.js',
    'SmilesParser.js',
    'SmilesWriter.js',
    'MolfileWriter.js',
    'Renderer.js',
    'History.js',
    'Tools.js',
    'CIPStereo.js',
    'SmartsParser.js',
    'SmartsMatch.js',
    'SmartsWriter.js',
    'ImageExport.js',
    'ExportStamp.js',
    'ToolbarPrefs.js',
    'MolEditor.js',
    'SMSDVersion.js',
    'SMSDGraph.js',
    'SMSDVF2.js',
    'SMSDMCS.js',
    'SMSDRings.js',
    'SMSDBatch.js',
    'SMSDLayout.js',
    'SDGLayout.js',
    'sdg/Congestion.js',
    'sdg/LayoutQuality.js',
    'sdg/HydrogenPlacer.js',
    'sdg/OverlapResolver.js',
    'sdg/AtomPlacer.js',
    'sdg/RingPlacer.js',
    'sdg/MacroCycleLayout.js',
    'sdg/IdentityTemplateLibrary.js',
    'sdg/TemplateHandler.js',
    'sdg/LayoutRefiner.js',
    'sdg/CorrectGeometricConfiguration.js',
    'sdg/NonplanarBonds.js',
    'sdg/StructureDiagramGenerator.js',
    'MLDepict.js',
    'ml-depict-weights.js',
    'RDT.js',
    'AtomTrace.js',
    'PageFormats.js',
    'CanvasView.js',
    'CanvasSurface.js',
    'FocusLens.js',
];

function scriptTags(version) {
    return FILES.map(function (file) {
        var suffix = version ? '?v=' + version : '';
        return '<script src="editor/' + file + suffix + '"></script>';
    }).join('\n');
}

module.exports = {
    FILES: FILES,
    scriptTags: scriptTags
};
