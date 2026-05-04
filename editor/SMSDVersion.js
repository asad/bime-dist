(function() {
    'use strict';
    /**
     * Version metadata for the SMSD-derived modules used inside BIME's
     * client-side editor. The original SMSD toolkit (Rahman et al.,
     * J. Cheminform. 2009;1:12) is the upstream Java reference; the
     * `editor/SMSD*.js` files are the BIME JavaScript port.
     *
     * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK
     * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
     * Licensed under the Apache License, Version 2.0
     */
    window.SMSD_VERSION = {
        // BIME release this version metadata was bundled with.
        bimeVersion: '1.8.7',
        modules: {
            SMSDGraph:  { source: 'mol_graph', status: 'ported' },
            SMSDVF2:    { source: 'vf2pp',     status: 'ported' },
            SMSDMCS:    { source: 'mcs',       status: 'ported' },
            SMSDRings:  { source: 'rings',     status: 'ported' },
            SMSDBatch:  { source: 'batch',     status: 'ported' },
            SMSDLayout: { source: 'layout',    status: 'ported' }
        },
        // Detailed change history is maintained in BIME's top-level
        // CHANGELOG.md rather than duplicated here.
        citation: 'Rahman SA et al. Small Molecule Subgraph Detector (SMSD) toolkit. J Cheminform. 2009;1:12.'
    };

    // Back-compat aliases for v1.0.1 and earlier consumers that referenced
    // the camel-cased symbols (Smsd...). The canonical product-name spelling
    // is the all-caps acronym (SMSD); the camel-cased aliases will be
    // removed in a future major release. Window-level re-exports below
    // delegate to whichever modules have already loaded by the time this
    // file runs (which is last in the standard load order).
    var W = window;
    if (W.SMSDGraph  && !W.SmsdGraph)  W.SmsdGraph  = W.SMSDGraph;
    if (W.SMSDVF2    && !W.SmsdVF2)    W.SmsdVF2    = W.SMSDVF2;
    if (W.SMSDMCS    && !W.SmsdMCS)    W.SmsdMCS    = W.SMSDMCS;
    if (W.SMSDRings  && !W.SmsdRings)  W.SmsdRings  = W.SMSDRings;
    if (W.SMSDBatch  && !W.SmsdBatch)  W.SmsdBatch  = W.SMSDBatch;
    if (W.SMSDLayout && !W.SmsdLayout) W.SmsdLayout = W.SMSDLayout;
})();
