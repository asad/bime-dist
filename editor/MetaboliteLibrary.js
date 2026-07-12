(function(global) {
    'use strict';
    /**
     * MetaboliteLibrary.js — a built-in repertoire of common metabolites and a
     * few canonical pathway templates, so pathways can be drawn with real
     * molecular structures and names instead of bare text labels.
     *
     * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
     * Licensed under the Apache License, Version 2.0 - see LICENSE.txt
     *
     * Design notes
     * ------------
     * - Small intermediates carry `display: 'structure'` and are drawn from
     *   their SMILES; large cofactors (ATP, NAD+, CoA, ...) carry
     *   `display: 'label'` and render as labelled nodes by default (the
     *   standard pathway-map convention), with SMILES kept for on-demand drawing
     *   where a compact one exists.
     * - Every SMILES here has been validated against BIME's own SmilesParser
     *   (see tests/test_v2_0_19_metabolite_library.js). KEGG COMPOUND ids are
     *   included for cross-referencing with KGML imports.
     * - Reaction facts (substrate -> product, enzyme, EC) are uncopyrightable
     *   data; the glycolysis template mirrors KEGG map00010.
     */

    var METABOLITES = [
        // ---- glycolysis intermediates (drawn structures) ----
        { name: 'D-Glucose',            smiles: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',                 kegg: 'C00031', display: 'structure', aliases: ['Glucose', 'Glc', 'alpha-D-Glucose'] },
        { name: 'Glucose 6-phosphate',  smiles: 'O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',           kegg: 'C00668', display: 'structure', aliases: ['G6P', 'Glucose-6-P'] },
        { name: 'Fructose 6-phosphate', smiles: 'OCC1(O)O[C@H](COP(O)(O)=O)[C@@H](O)[C@@H]1O',             kegg: 'C05345', display: 'structure', aliases: ['F6P', 'Fructose-6-P'] },
        { name: 'Fructose 1,6-bisphosphate', smiles: 'O[C@H]1[C@H](O)C(O)(COP(O)(O)=O)O[C@@H]1COP(O)(O)=O', kegg: 'C00354', display: 'structure', aliases: ['F1,6BP', 'FBP', 'F16BP'] },
        { name: 'Dihydroxyacetone phosphate', smiles: 'OCC(=O)COP(O)(O)=O',                                 kegg: 'C00111', display: 'structure', aliases: ['DHAP', 'Glycerone phosphate'] },
        { name: 'Glyceraldehyde 3-phosphate', smiles: 'O=C[C@H](O)COP(O)(O)=O',                             kegg: 'C00118', display: 'structure', aliases: ['GAP', 'G3P', 'Glyceraldehyde-3-P'] },
        { name: '1,3-Bisphosphoglycerate', smiles: 'O[C@H](COP(O)(O)=O)C(=O)OP(O)(O)=O',                     kegg: 'C00236', display: 'structure', aliases: ['1,3-BPG', '1,3BPG', 'BPG'] },
        { name: '3-Phosphoglycerate',   smiles: 'O[C@H](COP(O)(O)=O)C(O)=O',                                 kegg: 'C00197', display: 'structure', aliases: ['3-PG', '3PG'] },
        { name: '2-Phosphoglycerate',   smiles: 'OC[C@@H](OP(O)(O)=O)C(O)=O',                                kegg: 'C00631', display: 'structure', aliases: ['2-PG', '2PG'] },
        { name: 'Phosphoenolpyruvate',  smiles: 'OC(=O)C(=C)OP(O)(O)=O',                                     kegg: 'C00074', display: 'structure', aliases: ['PEP'] },
        { name: 'Pyruvate',             smiles: 'CC(=O)C(O)=O',                                              kegg: 'C00022', display: 'structure', aliases: ['Pyr'] },
        { name: 'Lactate',              smiles: 'C[C@H](O)C(O)=O',                                           kegg: 'C00186', display: 'structure', aliases: ['L-Lactate', '(S)-Lactate'] },
        { name: 'Acetaldehyde',         smiles: 'CC=O',                                                     kegg: 'C00084', display: 'structure', aliases: ['Ethanal'] },
        // ---- small inorganics / side metabolites ----
        { name: 'Phosphate',            smiles: 'OP(O)(O)=O',                                               kegg: 'C00009', display: 'label', aliases: ['Pi', 'Orthophosphate', 'Phosphate (Pi)'] },
        { name: 'Water',                smiles: 'O',                                                        kegg: 'C00001', display: 'label', aliases: ['H2O'] },
        { name: 'Carbon dioxide',       smiles: 'O=C=O',                                                    kegg: 'C00011', display: 'label', aliases: ['CO2'] },
        { name: 'Ammonia',              smiles: 'N',                                                        kegg: 'C00014', display: 'label', aliases: ['NH3'] },
        { name: 'Proton',               smiles: '[H+]',                                                     kegg: 'C00080', display: 'label', aliases: ['H+'] },
        // ---- large cofactors (labelled nodes) ----
        { name: 'ATP',                  smiles: '', kegg: 'C00002', display: 'label', aliases: ['Adenosine triphosphate'] },
        { name: 'ADP',                  smiles: '', kegg: 'C00008', display: 'label', aliases: ['Adenosine diphosphate'] },
        { name: 'AMP',                  smiles: '', kegg: 'C00020', display: 'label', aliases: ['Adenosine monophosphate'] },
        { name: 'NAD+',                 smiles: '', kegg: 'C00003', display: 'label', aliases: ['NAD'] },
        { name: 'NADH',                 smiles: '', kegg: 'C00004', display: 'label', aliases: [] },
        { name: 'NADP+',                smiles: '', kegg: 'C00006', display: 'label', aliases: ['NADP'] },
        { name: 'NADPH',                smiles: '', kegg: 'C00005', display: 'label', aliases: [] },
        { name: 'CoA',                  smiles: '', kegg: 'C00010', display: 'label', aliases: ['Coenzyme A', 'CoA-SH'] },
        { name: 'Acetyl-CoA',           smiles: '', kegg: 'C00024', display: 'label', aliases: ['Acetyl coenzyme A'] },
        { name: 'FAD',                  smiles: '', kegg: 'C00016', display: 'label', aliases: [] },
        { name: 'FADH2',                smiles: '', kegg: 'C01352', display: 'label', aliases: [] },
        { name: 'GTP',                  smiles: '', kegg: 'C00044', display: 'label', aliases: [] },
        { name: 'GDP',                  smiles: '', kegg: 'C00035', display: 'label', aliases: [] },
        // ---- v2.0.23 TCA cycle intermediates (drawn structures) ----
        { name: 'Oxaloacetate',         smiles: 'OC(=O)CC(=O)C(O)=O',                                     kegg: 'C00036', display: 'structure', aliases: ['OAA', 'Oxalacetic acid'] },
        { name: 'Citrate',              smiles: 'OC(=O)CC(O)(C(O)=O)CC(O)=O',                             kegg: 'C00158', display: 'structure', aliases: ['Citric acid'] },
        { name: 'cis-Aconitate',        smiles: 'OC(=O)CC(=CC(O)=O)C(O)=O',                                kegg: 'C00417', display: 'structure', aliases: ['cis-Aconitic acid'] },
        { name: 'Isocitrate',           smiles: 'O=C(O)C[C@H](C(=O)O)[C@@H](O)C(=O)O',                     kegg: 'C00311', display: 'structure', aliases: ['Isocitric acid', '(2R,3S)-Isocitrate'] },
        { name: '2-Oxoglutarate',       smiles: 'OC(=O)CCC(=O)C(O)=O',                                    kegg: 'C00026', display: 'structure', aliases: ['alpha-Ketoglutarate', 'aKG', '2-Ketoglutarate', 'Oxoglutaric acid'] },
        { name: 'Succinyl-CoA',         smiles: 'OC(=O)CCC(=O)SC',                                         kegg: 'C00091', display: 'structure', aliases: ['Succinyl coenzyme A'] },
        { name: 'Succinate',            smiles: 'OC(=O)CCC(O)=O',                                          kegg: 'C00042', display: 'structure', aliases: ['Succinic acid', 'Butanedioic acid'] },
        { name: 'Fumarate',             smiles: 'OC(=O)/C=C/C(O)=O',                                       kegg: 'C00122', display: 'structure', aliases: ['Fumaric acid', '(E)-Butenedioic acid'] },
        { name: 'L-Malate',             smiles: 'OC(=O)[C@@H](O)CC(O)=O',                                  kegg: 'C00149', display: 'structure', aliases: ['Malate', 'Malic acid', '(S)-Malate'] },
        // ---- v2.0.25 Urea cycle intermediates (Krebs-Henseleit) ----
        { name: 'L-Ornithine',          smiles: 'NCCC[C@H](N)C(O)=O',                                       kegg: 'C00077', display: 'structure', aliases: ['Ornithine', 'Orn', '(S)-Ornithine'] },
        { name: 'L-Citrulline',         smiles: 'NC(=O)NCCC[C@H](N)C(O)=O',                                 kegg: 'C00327', display: 'structure', aliases: ['Citrulline', 'Cit'] },
        { name: 'L-Argininosuccinate',  smiles: 'N=C(NCCC[C@H](N)C(O)=O)N[C@@H](CC(O)=O)C(O)=O',            kegg: 'C03406', display: 'structure', aliases: ['Argininosuccinate', 'ASA'] },
        { name: 'L-Arginine',           smiles: 'N=C(N)NCCC[C@H](N)C(O)=O',                                 kegg: 'C00062', display: 'structure', aliases: ['Arginine', 'Arg'] },
        { name: 'Carbamoyl phosphate',  smiles: 'NC(=O)OP(O)(O)=O',                                         kegg: 'C00169', display: 'structure', aliases: ['Carbamyl phosphate', 'CP'] },
        { name: 'Urea',                 smiles: 'NC(N)=O',                                                  kegg: 'C00086', display: 'structure', aliases: ['Carbamide'] },
        // ---- v2.0.26 Pentose Phosphate Pathway intermediates ----
        { name: '6-Phosphoglucono-1,5-lactone', smiles: 'O=C1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H](O)[C@@H]1O', kegg: 'C01236', display: 'structure', aliases: ['6PGL', '6-Phosphogluconolactone'] },
        { name: '6-Phosphogluconate',   smiles: 'OC(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)COP(O)(O)=O',       kegg: 'C00345', display: 'structure', aliases: ['6PG', 'Phosphogluconate'] },
        { name: 'D-Ribulose 5-phosphate', smiles: 'OCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O',                      kegg: 'C00199', display: 'structure', aliases: ['Ru5P', 'Ribulose-5-P'] },
        { name: 'D-Ribose 5-phosphate', smiles: 'O=C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O',                   kegg: 'C00117', display: 'structure', aliases: ['R5P', 'Ribose-5-P'] },
        { name: 'D-Xylulose 5-phosphate', smiles: 'OCC(=O)[C@@H](O)[C@H](O)COP(O)(O)=O',                     kegg: 'C00231', display: 'structure', aliases: ['Xu5P', 'Xylulose-5-P'] },
        // ---- v2.0.27 β-oxidation intermediates (S-methyl thioester CoA abbreviation per v2.0.23 Succinyl-CoA convention) ----
        { name: 'Butyryl-CoA',                     smiles: 'CCCC(=O)SC',         kegg: 'C00136', display: 'structure', aliases: ['Butanoyl-CoA', 'BuCoA'] },
        { name: 'trans-Crotonyl-CoA',              smiles: 'C/C=C/C(=O)SC',      kegg: 'C00877', display: 'structure', aliases: ['Crotonyl-CoA', '(E)-Crotonoyl-CoA', 'Δ²-trans-Enoyl-CoA'] },
        { name: 'L-3-Hydroxybutyryl-CoA',          smiles: 'C[C@H](O)CC(=O)SC',  kegg: 'C01144', display: 'structure', aliases: ['(S)-3-Hydroxybutanoyl-CoA', 'L-β-Hydroxybutyryl-CoA', '3-OH-Butyryl-CoA'] },
        { name: 'Acetoacetyl-CoA',                 smiles: 'CC(=O)CC(=O)SC',     kegg: 'C00332', display: 'structure', aliases: ['3-Ketobutyryl-CoA', '3-Oxobutanoyl-CoA'] },
        // ---- v2.0.76 pathway-template expansion ----
        // Gluconeogenesis + Cori cycle reuse the glycolysis / TCA intermediates
        // already above (zero new metabolites). These six add what the
        // fermentation, glyoxylate, and ketogenesis templates need. All are
        // small, well-known structures; the ketone-body D-β-hydroxybutyrate
        // is the (R)-3-hydroxybutanoate — verified by CIP (priority
        // OH > CH2COOH > CH3 > H) against this file's (S)-L-Lactate anchor,
        // which shares the identical stereocentre neighbour order: the (R)
        // token is therefore '@@' (the opposite of L-Lactate's '@').
        { name: 'Ethanol',              smiles: 'CCO',              kegg: 'C00469', display: 'structure', aliases: ['EtOH', 'Ethyl alcohol'] },
        { name: 'Glyoxylate',           smiles: 'O=CC(O)=O',        kegg: 'C00048', display: 'structure', aliases: ['Glyoxylic acid', 'Oxoacetate'] },
        { name: 'HMG-CoA',              smiles: '',                 kegg: 'C00356', display: 'label',     aliases: ['3-Hydroxy-3-methylglutaryl-CoA', 'Hydroxymethylglutaryl-CoA'] },
        { name: 'Acetoacetate',         smiles: 'CC(=O)CC(O)=O',    kegg: 'C00164', display: 'structure', aliases: ['Acetoacetic acid', '3-Oxobutanoate'] },
        { name: 'D-β-Hydroxybutyrate',  smiles: 'C[C@@H](O)CC(O)=O', kegg: 'C01089', display: 'structure', aliases: ['β-Hydroxybutyrate', '3-Hydroxybutyrate', 'D-3-Hydroxybutyrate', 'BHB', 'R-3-Hydroxybutanoate'] },
        { name: 'Acetone',              smiles: 'CC(C)=O',          kegg: 'C00207', display: 'structure', aliases: ['Propan-2-one', 'Dimethyl ketone'] }
    ];

    // Case-insensitive lookup by canonical name or any alias.
    var _index = null;
    function buildIndex() {
        _index = {};
        for (var i = 0; i < METABOLITES.length; i++) {
            var m = METABOLITES[i];
            _index[m.name.toLowerCase()] = m;
            for (var a = 0; a < m.aliases.length; a++) {
                var key = m.aliases[a].toLowerCase();
                if (!_index[key]) { _index[key] = m; }
            }
        }
    }
    function find(name) {
        if (!name) { return null; }
        if (!_index) { buildIndex(); }
        return _index[String(name).toLowerCase()] || null;
    }

    // v2.4.14: structure lookup by canonical SMILES — powers name/KEGG-id
    // auto-labelling of figures. Lazily canonicalises every library entry's
    // SMILES (via the shared SmilesParser + SmilesWriter) into an index, then
    // matches a query structure by its canonical form. Returns the metabolite
    // entry ({ name, kegg, ... }) or null. Requires SmilesParser + SmilesWriter;
    // returns null when they (or a parseable structure) are unavailable.
    // A raw canonical-SMILES string is NOT a reliable structure key: the
    // constitutional-only canonical ranking emits tetrahedral @/@@ in an
    // atom-order-dependent way (the SmilesWriter KNOWN_ORBIT_DRIFT note), so the
    // SAME molecule obtained two ways can serialise with a flipped @/@@. We
    // therefore key on (constitutional canonical SMILES, stereo-stripped, which
    // IS order-invariant) PLUS an order-invariant, stereo-aware CIP R/S multiset.
    // That is robust to the drift yet still distinguishes epimers/diastereomers.
    // NOTE: the SORTED R/S multiset is invariant under reflection, so it does
    // NOT distinguish enantiomers (L-glucose keys the same as D-glucose). This
    // is harmless here: the curated library has no enantiomeric pairs and the
    // metabolic reactions this names are single-enantiomer, and the
    // unambiguous-single-match guard below never guesses on a tie.
    var _structIndex = null;
    function _stripStereo(s) { return s ? s.replace(/@+/g, '').replace(/[\/\\]/g, '') : s; }
    function _structureKey(smi) {
        if (!smi || !global.SmilesParser || !global.SmilesWriter) { return null; }
        try {
            var m = global.SmilesParser.parse(smi);
            if (!m || !m.atoms || !m.atoms.length) { return null; }
            if (m.parseErrors && m.parseErrors.length) { return null; }
            var constitutional = _stripStereo(global.SmilesWriter.write(m));
            var stereo = '';
            if (global.CIPStereo && global.CIPStereo.assignRS) {
                global.CIPStereo.assignRS(m);
                var labels = [];
                for (var i = 0; i < m.atoms.length; i++) {
                    if (m.atoms[i].cipLabel) { labels.push(m.atoms[i].cipLabel); }
                }
                labels.sort();
                stereo = labels.join('');
            }
            // v2.4.16: also fold in the canonical E/Z multiset (from the parsed
            // directional bonds, not geometry) so cis/trans isomers — maleate vs
            // fumarate — get distinct keys and are never conflated.
            var ez = '';
            if (global.CIPStereo && global.CIPStereo.ezDescriptors) {
                ez = global.CIPStereo.ezDescriptors(m).join('');
            }
            return constitutional + '|' + stereo + '|' + ez;
        } catch (e) { return null; }
    }
    function buildStructIndex() {
        _structIndex = {};
        for (var i = 0; i < METABOLITES.length; i++) {
            var m = METABOLITES[i];
            if (!m.smiles) { continue; }
            var k = _structureKey(m.smiles);
            if (!k) { continue; }
            if (!_structIndex[k]) { _structIndex[k] = []; }
            _structIndex[k].push(m);
        }
    }
    function findBySmiles(smiles) {
        var k = _structureKey(smiles);
        if (!k) { return null; }
        if (!_structIndex) { buildStructIndex(); }
        var hits = _structIndex[k];
        // Only an UNAMBIGUOUS match (exactly one library compound for this
        // constitution + stereo signature) is returned — never mislabel.
        return (hits && hits.length === 1) ? hits[0] : null;
    }

    /**
     * Canonical pathway templates (reaction facts mirroring KEGG maps). Each
     * step is one reaction: from -> to, an enzyme + EC, a reaction-arrow type
     * ('forward' = irreversible, 'reversible' = equilibrium), and optional
     * cofactors consumed/produced (shown alongside the enzyme label).
     */
    var PATHWAYS = {
        glycolysis: {
            id: 'glycolysis',
            name: 'Glycolysis (Embden–Meyerhof–Parnas)',
            kegg: 'map00010',
            shape: 'ladder',
            steps: [
                { from: 'D-Glucose', to: 'Glucose 6-phosphate', enzyme: 'Hexokinase', ec: '2.7.1.1', arrowType: 'forward', cofactors: 'ATP → ADP' },
                { from: 'Glucose 6-phosphate', to: 'Fructose 6-phosphate', enzyme: 'Phosphoglucose isomerase', ec: '5.3.1.9', arrowType: 'reversible', cofactors: '' },
                { from: 'Fructose 6-phosphate', to: 'Fructose 1,6-bisphosphate', enzyme: 'Phosphofructokinase-1', ec: '2.7.1.11', arrowType: 'forward', cofactors: 'ATP → ADP' },
                { from: 'Fructose 1,6-bisphosphate', to: 'Glyceraldehyde 3-phosphate', enzyme: 'Aldolase', ec: '4.1.2.13', arrowType: 'reversible', cofactors: '' },
                { from: 'Fructose 1,6-bisphosphate', to: 'Dihydroxyacetone phosphate', enzyme: 'Aldolase', ec: '4.1.2.13', arrowType: 'reversible', cofactors: '' },
                { from: 'Dihydroxyacetone phosphate', to: 'Glyceraldehyde 3-phosphate', enzyme: 'Triosephosphate isomerase', ec: '5.3.1.1', arrowType: 'reversible', cofactors: '' },
                { from: 'Glyceraldehyde 3-phosphate', to: '1,3-Bisphosphoglycerate', enzyme: 'GAPDH', ec: '1.2.1.12', arrowType: 'reversible', cofactors: 'NAD+ + Pi → NADH' },
                { from: '1,3-Bisphosphoglycerate', to: '3-Phosphoglycerate', enzyme: 'Phosphoglycerate kinase', ec: '2.7.2.3', arrowType: 'reversible', cofactors: 'ADP → ATP' },
                { from: '3-Phosphoglycerate', to: '2-Phosphoglycerate', enzyme: 'Phosphoglycerate mutase', ec: '5.4.2.11', arrowType: 'reversible', cofactors: '' },
                { from: '2-Phosphoglycerate', to: 'Phosphoenolpyruvate', enzyme: 'Enolase', ec: '4.2.1.11', arrowType: 'reversible', cofactors: '− H2O' },
                { from: 'Phosphoenolpyruvate', to: 'Pyruvate', enzyme: 'Pyruvate kinase', ec: '2.7.1.40', arrowType: 'forward', cofactors: 'ADP → ATP' }
            ]
        },
        // v2.0.23: TCA / citric-acid / Krebs cycle. Canonical closed 9-reaction
        // cycle starting at oxaloacetate (clockwise from the 12 o'clock position
        // under the cycle-layout primitive). Acetyl-CoA enters as the cofactor
        // partner of the first step (citrate synthase). Per turn: 3 NADH, 1
        // FADH2, 1 GTP, 2 CO2. Source: KEGG map00020 + Lehninger ch. 16.
        tca: {
            id: 'tca',
            name: 'TCA cycle (Citric acid / Krebs)',
            kegg: 'map00020',
            shape: 'cycle',
            steps: [
                { from: 'Oxaloacetate',         to: 'Citrate',         enzyme: 'Citrate synthase',                    ec: '2.3.3.1',  arrowType: 'forward',    cofactors: '+ Acetyl-CoA + H2O → CoA-SH' },
                { from: 'Citrate',              to: 'cis-Aconitate',   enzyme: 'Aconitase',                            ec: '4.2.1.3',  arrowType: 'reversible', cofactors: '− H2O' },
                { from: 'cis-Aconitate',        to: 'Isocitrate',      enzyme: 'Aconitase',                            ec: '4.2.1.3',  arrowType: 'reversible', cofactors: '+ H2O' },
                { from: 'Isocitrate',           to: '2-Oxoglutarate',  enzyme: 'Isocitrate dehydrogenase',             ec: '1.1.1.41', arrowType: 'forward',    cofactors: 'NAD+ → NADH + CO2' },
                { from: '2-Oxoglutarate',       to: 'Succinyl-CoA',    enzyme: '2-Oxoglutarate dehydrogenase complex', ec: '1.2.4.2',  arrowType: 'forward',    cofactors: 'NAD+ + CoA-SH → NADH + CO2' },
                { from: 'Succinyl-CoA',         to: 'Succinate',       enzyme: 'Succinyl-CoA synthetase',              ec: '6.2.1.4',  arrowType: 'reversible', cofactors: 'GDP + Pi → GTP + CoA-SH' },
                { from: 'Succinate',            to: 'Fumarate',        enzyme: 'Succinate dehydrogenase',              ec: '1.3.5.1',  arrowType: 'reversible', cofactors: 'FAD → FADH2' },
                { from: 'Fumarate',             to: 'L-Malate',        enzyme: 'Fumarase',                             ec: '4.2.1.2',  arrowType: 'reversible', cofactors: '+ H2O' },
                { from: 'L-Malate',             to: 'Oxaloacetate',    enzyme: 'Malate dehydrogenase',                 ec: '1.1.1.37', arrowType: 'reversible', cofactors: 'NAD+ → NADH' }
            ]
        },
        // v2.0.25: Urea cycle (Krebs–Henseleit). Four-step closed cycle on the
        // ornithine backbone; physiologically irreversible in the forward
        // direction. Carbamoyl phosphate enters as a cofactor at step 0 (built
        // upstream from NH3 + CO2 + 2 ATP by carbamoyl-P synthetase I,
        // EC 6.3.4.16, in the mitochondrion); urea exits at the final step;
        // aspartate enters / fumarate exits at the lyase step. KEGG map00220.
        urea: {
            id: 'urea',
            name: 'Urea cycle (Krebs–Henseleit)',
            kegg: 'map00220',
            shape: 'cycle',
            steps: [
                { from: 'L-Ornithine',         to: 'L-Citrulline',         enzyme: 'Ornithine transcarbamylase',  ec: '2.1.3.3',  arrowType: 'forward', cofactors: '+ Carbamoyl-P → Pi' },
                { from: 'L-Citrulline',        to: 'L-Argininosuccinate',  enzyme: 'Argininosuccinate synthetase', ec: '6.3.4.5',  arrowType: 'forward', cofactors: '+ Aspartate + ATP → AMP + PPi' },
                { from: 'L-Argininosuccinate', to: 'L-Arginine',           enzyme: 'Argininosuccinate lyase',      ec: '4.3.2.1',  arrowType: 'forward', cofactors: '→ Fumarate' },
                { from: 'L-Arginine',          to: 'L-Ornithine',          enzyme: 'Arginase',                    ec: '3.5.3.1',  arrowType: 'forward', cofactors: '+ H2O → Urea' }
            ]
        },
        // v2.0.26: Pentose Phosphate Pathway. Oxidative phase (G6P → 6PG-lactone
        // → 6PG → Ribulose-5-P, 2 NADPH + CO2 produced) is the linear trunk;
        // ribulose-5-P then branches into the non-oxidative phase via the
        // isomerase (→ Ribose-5-P) and the epimerase (→ Xylulose-5-P). KEGG
        // map00030. The downstream sugar-shuffle (transketolase / transaldolase
        // shuttling C2/C3 units to/from glycolytic F6P + GAP) is intentionally
        // not included in this minimal template.
        ppp: {
            id: 'ppp',
            name: 'Pentose Phosphate Pathway',
            kegg: 'map00030',
            shape: 'branched',
            trunk: ['Glucose 6-phosphate', '6-Phosphoglucono-1,5-lactone', '6-Phosphogluconate', 'D-Ribulose 5-phosphate'],
            steps: [
                { from: 'Glucose 6-phosphate',          to: '6-Phosphoglucono-1,5-lactone', enzyme: 'Glucose-6-phosphate dehydrogenase', ec: '1.1.1.49', arrowType: 'forward',    cofactors: 'NADP+ → NADPH' },
                { from: '6-Phosphoglucono-1,5-lactone', to: '6-Phosphogluconate',           enzyme: 'Gluconolactonase',                   ec: '3.1.1.31', arrowType: 'forward',    cofactors: '+ H2O' },
                { from: '6-Phosphogluconate',           to: 'D-Ribulose 5-phosphate',       enzyme: '6-Phosphogluconate dehydrogenase',   ec: '1.1.1.44', arrowType: 'forward',    cofactors: 'NADP+ → NADPH + CO2' },
                { from: 'D-Ribulose 5-phosphate',       to: 'D-Ribose 5-phosphate',         enzyme: 'Ribose-5-phosphate isomerase',       ec: '5.3.1.6',  arrowType: 'reversible', cofactors: '' },
                { from: 'D-Ribulose 5-phosphate',       to: 'D-Xylulose 5-phosphate',       enzyme: 'Ribulose-5-phosphate 3-epimerase',  ec: '5.1.3.1',  arrowType: 'reversible', cofactors: '' }
            ]
        },
        pyruvate_fates: {
            id: 'pyruvate_fates',
            name: 'Pyruvate fates',
            kegg: 'map00620',
            shape: 'fanout',
            hub: 'Pyruvate',
            steps: [
                { from: 'Pyruvate', to: 'Lactate',      enzyme: 'Lactate dehydrogenase',          ec: '1.1.1.27', arrowType: 'reversible', cofactors: 'NADH → NAD+' },
                { from: 'Pyruvate', to: 'Acetaldehyde', enzyme: 'Pyruvate decarboxylase',         ec: '4.1.1.1',  arrowType: 'forward',    cofactors: '→ CO2' },
                { from: 'Pyruvate', to: 'Acetyl-CoA',   enzyme: 'Pyruvate dehydrogenase system',  ec: '1.2.1.104', arrowType: 'forward',   cofactors: 'CoA + NAD+ → CO2 + NADH' },
                { from: 'Pyruvate', to: 'Oxaloacetate', enzyme: 'Pyruvate carboxylase',           ec: '6.4.1.1',   arrowType: 'forward',   cofactors: 'ATP + HCO3- → ADP + Pi' }
            ]
        },
        // v2.0.27: β-oxidation spiral (Lynen helix). One full turn shown on the
        // shortest meaningful substrate — butyryl-CoA (C4) → 2× acetyl-CoA. The
        // four enzyme steps are the canonical FAD-dependent dehydrogenation,
        // hydration, NAD-dependent dehydrogenation, and thiolytic cleavage; the
        // closing arrow from acetoacetyl-CoA back to butyryl-CoA carries the
        // thiolysis cofactor text "+ CoA-SH → 2× Acetyl-CoA". For longer chains
        // (e.g. palmitoyl C16), the same square geometry repeats with the
        // closure feeding a shortened acyl-CoA — the "iteration" semantics
        // encoded by the `loop-iterative` shape primitive. Source: KEGG
        // map00071 + Lehninger ch. 17.
        beta_oxidation: {
            id: 'beta_oxidation',
            name: 'β-Oxidation (fatty-acid spiral, one turn)',
            kegg: 'map00071',
            shape: 'loop-iterative',
            steps: [
                { from: 'Butyryl-CoA',              to: 'trans-Crotonyl-CoA',         enzyme: 'Acyl-CoA dehydrogenase',           ec: '1.3.8.7',  arrowType: 'forward',    cofactors: 'FAD → FADH2' },
                { from: 'trans-Crotonyl-CoA',       to: 'L-3-Hydroxybutyryl-CoA',     enzyme: 'Enoyl-CoA hydratase',              ec: '4.2.1.17', arrowType: 'reversible', cofactors: '+ H2O' },
                { from: 'L-3-Hydroxybutyryl-CoA',   to: 'Acetoacetyl-CoA',            enzyme: '3-Hydroxyacyl-CoA dehydrogenase',  ec: '1.1.1.35', arrowType: 'reversible', cofactors: 'NAD+ → NADH' },
                { from: 'Acetoacetyl-CoA',          to: 'Butyryl-CoA',                enzyme: 'β-Ketothiolase',                    ec: '2.3.1.16', arrowType: 'forward',    cofactors: '+ CoA-SH → 2× Acetyl-CoA' }
            ]
        },
        // ---- v2.0.76 pathway-template expansion (5 new pathways) ----
        // Gluconeogenesis: the anabolic glucose-synthesis ladder. It reverses
        // glycolysis except at the three physiologically-irreversible bypasses
        // (pyruvate carboxylase + PEPCK; FBPase-1; G6Pase). Reuses every
        // glycolysis intermediate — zero new metabolites. KEGG map00010.
        gluconeogenesis: {
            id: 'gluconeogenesis',
            name: 'Gluconeogenesis',
            kegg: 'map00010',
            shape: 'ladder',
            steps: [
                { from: 'Pyruvate',                  to: 'Oxaloacetate',           enzyme: 'Pyruvate carboxylase',          ec: '6.4.1.1',  arrowType: 'forward',    cofactors: 'ATP + HCO3- → ADP + Pi' },
                { from: 'Oxaloacetate',              to: 'Phosphoenolpyruvate',    enzyme: 'PEP carboxykinase',             ec: '4.1.1.32', arrowType: 'forward',    cofactors: 'GTP → GDP + CO2' },
                { from: 'Phosphoenolpyruvate',       to: '2-Phosphoglycerate',     enzyme: 'Enolase',                       ec: '4.2.1.11', arrowType: 'reversible', cofactors: '+ H2O' },
                { from: '2-Phosphoglycerate',        to: '3-Phosphoglycerate',     enzyme: 'Phosphoglycerate mutase',       ec: '5.4.2.11', arrowType: 'reversible', cofactors: '' },
                { from: '3-Phosphoglycerate',        to: '1,3-Bisphosphoglycerate',enzyme: 'Phosphoglycerate kinase',       ec: '2.7.2.3',  arrowType: 'reversible', cofactors: 'ATP → ADP' },
                { from: '1,3-Bisphosphoglycerate',   to: 'Glyceraldehyde 3-phosphate', enzyme: 'GAPDH',                     ec: '1.2.1.12', arrowType: 'reversible', cofactors: 'NADH → NAD+ + Pi' },
                { from: 'Glyceraldehyde 3-phosphate',to: 'Fructose 1,6-bisphosphate', enzyme: 'Aldolase',                  ec: '4.1.2.13', arrowType: 'reversible', cofactors: '+ DHAP' },
                { from: 'Fructose 1,6-bisphosphate', to: 'Fructose 6-phosphate',   enzyme: 'Fructose-1,6-bisphosphatase',   ec: '3.1.3.11', arrowType: 'forward',    cofactors: '+ H2O → Pi' },
                { from: 'Fructose 6-phosphate',      to: 'Glucose 6-phosphate',    enzyme: 'Phosphoglucose isomerase',      ec: '5.3.1.9',  arrowType: 'reversible', cofactors: '' },
                { from: 'Glucose 6-phosphate',       to: 'D-Glucose',              enzyme: 'Glucose-6-phosphatase',         ec: '3.1.3.9',  arrowType: 'forward',    cofactors: '+ H2O → Pi' }
            ]
        },
        // Cori cycle: the inter-organ lactate↔glucose loop. Muscle runs
        // glycolysis to lactate; liver runs gluconeogenesis back to glucose.
        // Shown compactly as a 3-node cycle with the multi-step legs as summary
        // arrows. Reuses glucose / pyruvate / lactate — zero new metabolites.
        cori: {
            id: 'cori',
            name: 'Cori cycle (lactate–glucose)',
            kegg: 'map00010',
            shape: 'cycle',
            steps: [
                { from: 'D-Glucose', to: 'Pyruvate',  enzyme: 'Glycolysis (muscle)',       ec: '', arrowType: 'forward',    cofactors: 'ADP + Pi → ATP' },
                { from: 'Pyruvate',  to: 'Lactate',   enzyme: 'Lactate dehydrogenase',     ec: '1.1.1.27', arrowType: 'reversible', cofactors: 'NADH → NAD+' },
                { from: 'Lactate',   to: 'D-Glucose', enzyme: 'Gluconeogenesis (liver)',   ec: '', arrowType: 'forward',    cofactors: 'ATP/GTP → ADP/GDP' }
            ]
        },
        // Ethanol fermentation: pyruvate → acetaldehyde → ethanol, the anaerobic
        // NAD+-regenerating fate of pyruvate in yeast. KEGG map00010.
        fermentation: {
            id: 'fermentation',
            name: 'Ethanol fermentation',
            kegg: 'map00010',
            shape: 'ladder',
            steps: [
                { from: 'Pyruvate',     to: 'Acetaldehyde', enzyme: 'Pyruvate decarboxylase', ec: '4.1.1.1', arrowType: 'forward',    cofactors: '→ CO2' },
                { from: 'Acetaldehyde', to: 'Ethanol',      enzyme: 'Alcohol dehydrogenase',  ec: '1.1.1.1', arrowType: 'reversible', cofactors: 'NADH → NAD+' }
            ]
        },
        // Glyoxylate cycle: the plant/microbe anaplerotic bypass that skips the
        // two decarboxylations of the TCA cycle, letting net carbon be built
        // from acetyl-CoA. Isocitrate lyase + malate synthase are the signature
        // enzymes; the rest of the ring shares TCA enzymes. KEGG map00630.
        glyoxylate: {
            id: 'glyoxylate',
            name: 'Glyoxylate cycle',
            kegg: 'map00630',
            shape: 'cycle',
            steps: [
                { from: 'Isocitrate',   to: 'Glyoxylate',    enzyme: 'Isocitrate lyase',     ec: '4.1.3.1',  arrowType: 'forward',    cofactors: '→ Succinate' },
                { from: 'Glyoxylate',   to: 'L-Malate',      enzyme: 'Malate synthase',      ec: '2.3.3.9',  arrowType: 'forward',    cofactors: '+ Acetyl-CoA + H2O → CoA-SH' },
                { from: 'L-Malate',     to: 'Oxaloacetate',  enzyme: 'Malate dehydrogenase', ec: '1.1.1.37', arrowType: 'reversible', cofactors: 'NAD+ → NADH' },
                { from: 'Oxaloacetate', to: 'Citrate',       enzyme: 'Citrate synthase',     ec: '2.3.3.1',  arrowType: 'forward',    cofactors: '+ Acetyl-CoA + H2O → CoA-SH' },
                { from: 'Citrate',      to: 'cis-Aconitate', enzyme: 'Aconitase',            ec: '4.2.1.3',  arrowType: 'reversible', cofactors: '− H2O' },
                { from: 'cis-Aconitate',to: 'Isocitrate',    enzyme: 'Aconitase',            ec: '4.2.1.3',  arrowType: 'reversible', cofactors: '+ H2O' }
            ]
        },
        // Ketogenesis: hepatic synthesis of ketone bodies from acetyl-CoA via
        // the HMG-CoA branch. Acetoacetate is the hub: reduced to
        // D-β-hydroxybutyrate (the transport form) or spontaneously
        // decarboxylated to acetone. KEGG map00072.
        ketogenesis: {
            id: 'ketogenesis',
            name: 'Ketogenesis (ketone bodies)',
            kegg: 'map00072',
            shape: 'ladder',
            steps: [
                { from: 'Acetyl-CoA',        to: 'Acetoacetyl-CoA',       enzyme: 'Acetyl-CoA acetyltransferase (thiolase)', ec: '2.3.1.9',  arrowType: 'reversible', cofactors: '+ Acetyl-CoA → CoA-SH' },
                { from: 'Acetoacetyl-CoA',   to: 'HMG-CoA',               enzyme: 'HMG-CoA synthase',                        ec: '2.3.3.10', arrowType: 'forward',    cofactors: '+ Acetyl-CoA + H2O → CoA-SH' },
                { from: 'HMG-CoA',           to: 'Acetoacetate',          enzyme: 'HMG-CoA lyase',                           ec: '4.1.3.4',  arrowType: 'forward',    cofactors: '→ Acetyl-CoA' },
                { from: 'Acetoacetate',      to: 'D-β-Hydroxybutyrate',   enzyme: 'β-Hydroxybutyrate dehydrogenase',         ec: '1.1.1.30', arrowType: 'reversible', cofactors: 'NADH → NAD+' },
                { from: 'Acetoacetate',      to: 'Acetone',               enzyme: 'Acetoacetate decarboxylase (spontaneous)',ec: '4.1.1.4',  arrowType: 'forward',    cofactors: '→ CO2' }
            ]
        }
    };

    global.MetaboliteLibrary = {
        metabolites: METABOLITES,
        pathways: PATHWAYS,
        find: find,
        findBySmiles: findBySmiles,
        version: '2.1.2'
    };
})(typeof window !== 'undefined' ? window : this);
