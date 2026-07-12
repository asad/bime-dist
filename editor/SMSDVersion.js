(function() {
    'use strict';
    /**
     * Version metadata for BIME's substructure / MCS editor modules.
     *
     * Copyright (c) 2018-2026 BioInception PVT LTD, Cambridge, UK
     * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
     * Licensed under the Apache License, Version 2.0
     */
    window.SMSD_VERSION = {
        // BIME release this version metadata was bundled with.
        bimeVersion: '3.0.3',
        modules: {
            SMSDGraph:  { source: 'mol_graph', status: 'ready' },
            SMSDVF2:    { source: 'subgraph',  status: 'ready' },
            SMSDMCS:    { source: 'mcs',       status: 'ready' },
            SMSDRings:  { source: 'rings',     status: 'ready' },
            SMSDBatch:  { source: 'batch',     status: 'ready' },
            SMSDLayout: { source: 'layout',    status: 'ready' }
        }
        // Detailed change history is maintained in BIME's top-level
        // CHANGELOG.md rather than duplicated here.
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
