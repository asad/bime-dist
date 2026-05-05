/**
 * test-cases.js -- SMILES parser / writer QA test suite for BIME
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 *
 * Contains real-world drug SMILES, reaction SMILES, and known edge cases
 * drawn from public chemistry references and parser-spec documentation.
 *
 * Each entry records:
 *   name     - human-readable label
 *   smiles   - the SMILES string to parse
 *   atoms    - expected heavy-atom count (H excluded)
 *   bonds    - expected bond count (explicit bonds only)
 *   expected - "valid" or "error" (whether the parser should accept it)
 *   notes    - optional description of parser-relevant features
 */
var TEST_CASES = {

    // ==================================================================
    // 1.  Real drug molecules  (20+ entries)
    // ==================================================================
    molecules: [
        {
            name:  "Aspirin",
            smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
            atoms: 13, bonds: 13,
            expected: "valid",
            notes: "Kekule phenyl ring, ester + carboxylic acid"
        },
        {
            name:  "Ibuprofen",
            smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            atoms: 15, bonds: 15,
            expected: "valid",
            notes: "Aromatic ring (lowercase), branching, carboxylic acid"
        },
        {
            name:  "Paracetamol (Acetaminophen)",
            smiles: "CC(=O)Nc1ccc(O)cc1",
            atoms: 11, bonds: 11,
            expected: "valid",
            notes: "Aromatic ring, amide, phenol"
        },
        {
            name:  "Caffeine",
            smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
            atoms: 14, bonds: 16,
            expected: "valid",
            notes: "Fused aromatic bicyclic (purine), N-methyls, two carbonyls"
        },
        {
            name:  "Metformin",
            smiles: "CN(C)C(=N)NC(N)=N",
            atoms: 7, bonds: 7,
            expected: "valid",
            notes: "Guanidine groups, all organic-subset atoms"
        },
        {
            name:  "Omeprazole",
            smiles: "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",
            atoms: 23, bonds: 25,
            expected: "valid",
            notes: "Benzimidazole ([nH]), sulfinyl S(=O), pyridine, methoxy groups"
        },
        {
            name:  "Warfarin",
            smiles: "CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O",
            atoms: 21, bonds: 23,
            expected: "valid",
            notes: "Coumarin scaffold, two aromatic rings, ketone, enol"
        },
        {
            name:  "Diazepam",
            smiles: "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc31",
            atoms: 19, bonds: 21,
            expected: "valid",
            notes: "Seven-membered ring, Cl substituent, aromatic ring fusion"
        },
        {
            name:  "Morphine",
            smiles: "CN1CCC23C4OC5=C(O)C=CC(=C25)C(O)C1C3C=C4",
            atoms: 21, bonds: 24,
            expected: "valid",
            notes: "Complex bridged/fused polycyclic, multiple ring closures"
        },
        {
            name:  "Penicillin G",
            smiles: "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
            atoms: 23, bonds: 25,
            expected: "valid",
            notes: "Beta-lactam (4-ring), thiazolidine, amide, carboxylic acid"
        },
        {
            name:  "Amoxicillin",
            smiles: "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
            atoms: 25, bonds: 27,
            expected: "valid",
            notes: "Like penicillin G with amino-hydroxyphenyl group"
        },
        {
            name:  "Ciprofloxacin",
            smiles: "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
            atoms: 24, bonds: 27,
            expected: "valid",
            notes: "Fluoroquinolone, cyclopropyl, piperazine, F substituent"
        },
        {
            name:  "Atorvastatin",
            smiles: "CC(C)c1n(CC(O)CC(O)CC(=O)O)c(c2ccc(F)cc2)c(c1c3ccccc3)C(=O)Nc4ccccc4",
            atoms: 42, bonds: 45,
            expected: "valid",
            notes: "Statin, pyrrole core, multiple aromatic rings, amide, dihydroxy acid chain"
        },
        {
            name:  "Losartan",
            smiles: "CCCCc1nc(Cl)c(n1Cc2ccc(cc2)c3ccccc3c4nn[nH]n4)CO",
            atoms: 31, bonds: 34,
            expected: "valid",
            notes: "Imidazole, biphenyl, tetrazole ([nH]), Cl"
        },
        {
            name:  "Amlodipine",
            smiles: "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c2ccccc2Cl)C(=O)OC",
            atoms: 27, bonds: 28,
            expected: "valid",
            notes: "Dihydropyridine, ester groups, Cl, ether chain"
        },
        {
            name:  "Sertraline",
            smiles: "CNC1CCC(c2ccc(Cl)c(Cl)c2)c3ccccc31",
            atoms: 21, bonds: 23,
            expected: "valid",
            notes: "Naphthalene fused ring, dichlorophenyl, secondary amine"
        },
        {
            name:  "Fluoxetine (Prozac)",
            smiles: "CNCCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2",
            atoms: 22, bonds: 22,
            expected: "valid",
            notes: "Trifluoromethyl CF3, ether linkage, two aromatic rings"
        },
        {
            name:  "Clopidogrel",
            smiles: "COC(=O)C(c1ccccc1Cl)N2CCc3sccc3C2",
            atoms: 20, bonds: 22,
            expected: "valid",
            notes: "Thiophene (s in aromatic subset), seven-membered ring, Cl"
        },
        {
            name:  "Sildenafil (Viagra)",
            smiles: "CCCC1=NN(C)C2=C1NC(=NC2=O)c3cc(ccc3OCC)S(=O)(=O)N4CCN(CC4)C",
            atoms: 33, bonds: 36,
            expected: "valid",
            notes: "Pyrazolopyrimidine, sulfonamide, piperazine, ethoxy"
        },
        {
            name:  "Lisinopril",
            smiles: "NCCCC(NC(CCc1ccccc1)C(=O)O)C(=O)N2CCCC2C(=O)O",
            atoms: 27, bonds: 28,
            expected: "valid",
            notes: "ACE inhibitor, proline ring, two carboxylic acids, amine"
        },
        {
            name:  "Metoprolol",
            smiles: "CC(C)NCC(O)COc1ccc(CCOC)cc1",
            atoms: 18, bonds: 18,
            expected: "valid",
            notes: "Beta-blocker, aromatic ring, ether, secondary amine/alcohol"
        },
        {
            name:  "Doxycycline",
            smiles: "CC1C2C(O)C3C(=C(O)c4cccc(O)c4C3=O)C(=O)C2(O)C(N(C)C)C(=O)C1(C)O",
            atoms: 30, bonds: 33,
            expected: "valid",
            notes: "Tetracycline tetracyclic scaffold, multiple stereocentres, enol"
        },
        {
            name:  "Naproxen",
            smiles: "COc1ccc2cc(CC(C)C(=O)O)ccc2c1",
            atoms: 17, bonds: 18,
            expected: "valid",
            notes: "Naphthalene, methoxy, propionic acid"
        },
        {
            name:  "Glucose (open-chain)",
            smiles: "OCC(O)C(O)C(O)C(O)C=O",
            atoms: 12, bonds: 11,
            expected: "valid",
            notes: "Linear sugar, no rings, many hydroxyl groups"
        },
        {
            name:  "Cholesterol",
            smiles: "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4CC(O)CCC34C)C",
            atoms: 27, bonds: 30,
            expected: "valid",
            notes: "Steroid tetracyclic core, many ring closures (1-4)"
        }
    ],

    // ==================================================================
    // 2.  Reaction SMILES  (10+ entries with atom mapping)
    // ==================================================================
    reactions: [
        {
            name:    "Fischer esterification",
            smiles:  "[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]>>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5].[OH2:4]",
            mapped:  true,
            notes:   "Acetic acid + ethanol -> ethyl acetate + water"
        },
        {
            name:    "SN2 displacement (iodide)",
            smiles:  "[I-:1].[Na+].[CH2:2]([CH:3]=[CH2:4])[Br:5]>>[Na+].[Br-:5].[CH2:2]([CH:3]=[CH2:4])[I:1]",
            mapped:  true,
            notes:   "Allyl bromide + NaI -> allyl iodide + NaBr"
        },
        {
            name:    "Suzuki coupling",
            smiles:  "c1ccc([B:1](O)O)cc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1",
            mapped:  false,
            notes:   "Phenylboronic acid + bromobenzene -> biphenyl (no mapping)"
        },
        {
            name:    "Aldol condensation",
            smiles:  "[CH3:1][CH:2]=[O:3].[CH3:4][C:5](=[O:6])[CH3:7]>>[CH3:1][CH:2]([OH:3])[CH2:7][C:5](=[O:6])[CH3:4]",
            mapped:  true,
            notes:   "Acetaldehyde + acetone -> beta-hydroxy ketone"
        },
        {
            name:    "Amide bond formation",
            smiles:  "[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>>[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]",
            mapped:  true,
            notes:   "Acetic acid + methylamine -> N-methylacetamide + water"
        },
        {
            name:    "Diels-Alder cycloaddition",
            smiles:  "C=CC=C.C=C>>C1CC=CCC1",
            mapped:  false,
            notes:   "1,3-butadiene + ethylene -> cyclohexene"
        },
        {
            name:    "Grignard addition",
            smiles:  "[CH3:1][Mg]Br.[CH3:2][CH:3]=[O:4]>>[CH3:1][CH:3]([OH:4])[CH3:2]",
            mapped:  true,
            notes:   "MeMgBr + acetaldehyde -> 2-propanol; bracket [Mg] needed"
        },
        {
            name:    "Williamson ether synthesis",
            smiles:  "[CH3:1][O-:2].[Na+].[CH3:3][Br:4]>>[CH3:1][O:2][CH3:3].[Na+].[Br-:4]",
            mapped:  true,
            notes:   "Sodium methoxide + methyl bromide -> dimethyl ether"
        },
        {
            name:    "Saponification",
            smiles:  "[CH3:1][C:2](=[O:3])[O:4][CH2:5][CH3:6].[OH-:7]>>[CH3:1][C:2](=[O:3])[O-:7].[HO:4][CH2:5][CH3:6]",
            mapped:  true,
            notes:   "Ethyl acetate + hydroxide -> acetate + ethanol"
        },
        {
            name:    "Nucleophilic aromatic substitution",
            smiles:  "c1cc([N+](=O)[O-])ccc1[F:1].[NH2:2]c1ccccc1>>c1cc([N+](=O)[O-])ccc1[NH:2]c1ccccc1.[F-:1]",
            mapped:  true,
            notes:   "Fluoronitrobenzene + aniline -> diarylamine + F-; charged N"
        },
        {
            name:    "Beckmann rearrangement",
            smiles:  "[CH3:1][C:2](=[N:3][OH:4])[CH3:5]>>[CH3:1][NH:3][C:2](=[O:4])[CH3:5]",
            mapped:  true,
            notes:   "Acetone oxime -> N-methylacetamide"
        },
        {
            name:    "Wittig reaction",
            smiles:  "[CH2:1]=[P](c1ccccc1)(c2ccccc2)c3ccccc3.[CH3:2][CH:3]=[O:4]>>[CH2:1]=[CH:3][CH3:2].O=[P](c1ccccc1)(c2ccccc2)c3ccccc3",
            mapped:  true,
            notes:   "Methylenetriphenylphosphorane + acetaldehyde -> propene + Ph3PO; bracket [P] not needed (organic subset)"
        }
    ],

    // ==================================================================
    // 3.  Edge cases  (tricky SMILES that parsers commonly fail on)
    // ==================================================================
    edge_cases: [
        // --- Charged atoms ---
        {
            name:     "Ammonium cation",
            smiles:   "[NH4+]",
            expected: "valid",
            notes:    "Bracket atom: N with 4 explicit H and +1 charge"
        },
        {
            name:     "Hydroxide anion",
            smiles:   "[OH-]",
            expected: "valid",
            notes:    "Bracket atom: O with 1 H and -1 charge"
        },
        {
            name:     "Sodium cation",
            smiles:   "[Na+]",
            expected: "valid",
            notes:    "Metal ion, no organic-subset entry"
        },
        {
            name:     "Iron(III) cation",
            smiles:   "[Fe+3]",
            expected: "valid",
            notes:    "Multi-digit charge: +3 on transition metal"
        },
        {
            name:     "Nitro group (charge-separated)",
            smiles:   "C[N+]([O-])=O",
            expected: "valid",
            notes:    "Charge-separated nitro; both + and - charges in one molecule"
        },
        {
            name:     "Zwitterion (glycine)",
            smiles:   "[NH3+]CC([O-])=O",
            expected: "valid",
            notes:    "Amino acid zwitterion form; opposite charges"
        },

        // --- Isotopes ---
        {
            name:     "Deuterium (D)",
            smiles:   "[2H]",
            expected: "valid",
            notes:    "Isotope 2 on hydrogen"
        },
        {
            name:     "Carbon-13",
            smiles:   "[13C]",
            expected: "valid",
            notes:    "Isotope prefix inside brackets"
        },
        {
            name:     "Carbon-14 in aromatic ring",
            smiles:   "[14cH]1ccccc1",
            expected: "valid",
            notes:    "Isotope + aromatic bracket atom + Kekule ring"
        },
        {
            name:     "Deuterochloroform",
            smiles:   "[2H]C(Cl)(Cl)Cl",
            expected: "valid",
            notes:    "Isotopic H bonded to carbon with three Cl"
        },
        {
            name:     "Tritium water",
            smiles:   "[3H]O[3H]",
            expected: "valid",
            notes:    "Both hydrogens are tritium"
        },

        // --- Aromatic heterocycles ---
        {
            name:     "Pyrrole",
            smiles:   "c1cc[nH]c1",
            expected: "valid",
            notes:    "[nH] required for pyrrole nitrogen (lone pair contributes 2 e-)"
        },
        {
            name:     "Imidazole",
            smiles:   "c1cnc[nH]1",
            expected: "valid",
            notes:    "Two nitrogens: one pyridine-type (n), one pyrrole-type ([nH])"
        },
        {
            name:     "Furan",
            smiles:   "c1ccoc1",
            expected: "valid",
            notes:    "Aromatic oxygen; 'o' is in AROMATIC_ORGANIC"
        },
        {
            name:     "Thiophene",
            smiles:   "c1ccsc1",
            expected: "valid",
            notes:    "Aromatic sulfur; 's' is in AROMATIC_ORGANIC"
        },
        {
            name:     "Pyridine",
            smiles:   "c1ccncc1",
            expected: "valid",
            notes:    "Aromatic nitrogen without H (pyridine-type)"
        },
        {
            name:     "Indole",
            smiles:   "c1ccc2[nH]ccc2c1",
            expected: "valid",
            notes:    "Fused bicyclic aromatic heterocycle"
        },
        {
            name:     "Selenophene",
            smiles:   "[se]1cccc1",
            expected: "valid",
            notes:    "Aromatic selenium; [se] bracket required (not in organic subset)"
        },
        {
            name:     "Purine (adenine scaffold)",
            smiles:   "c1nc2c([nH]1)ncnc2N",
            expected: "valid",
            notes:    "Fused imidazole + pyrimidine with exocyclic amine"
        },

        // --- Ring edge cases ---
        {
            name:     "Biphenyl (explicit single bond)",
            smiles:   "c1ccccc1-c2ccccc2",
            expected: "valid",
            notes:    "Explicit '-' between aromatic rings prevents aromatic bond"
        },
        {
            name:     "Naphthalene",
            smiles:   "c1ccc2ccccc2c1",
            expected: "valid",
            notes:    "Fused bicyclic with ring closures 1 and 2"
        },
        {
            name:     "Cubane",
            smiles:   "C12C3C4C1C5C3C4C25",
            expected: "valid",
            notes:    "Multiple ring closures on single atoms (4 digits in a row)"
        },
        {
            name:     "Adamantane",
            smiles:   "C1C2CC3CC1CC(C2)C3",
            expected: "valid",
            notes:    "Bridged tricyclic cage structure; ring closures 1-3"
        },
        {
            name:     "Spiro compound (spiro[4.5]decane)",
            smiles:   "C1CCC2(CC1)CCCCC2",
            expected: "valid",
            notes:    "Spiro center with two ring closures on one atom"
        },
        {
            name:     "Large ring %nn notation",
            smiles:   "C%10CCCCCCCCCCCCC%10",
            expected: "valid",
            notes:    "14-membered ring using percent notation for ring number 10"
        },

        // --- Bond types ---
        {
            name:     "Acetylene (triple bond)",
            smiles:   "C#C",
            expected: "valid",
            notes:    "Triple bond symbol #"
        },
        {
            name:     "Allene",
            smiles:   "C=C=C",
            expected: "valid",
            notes:    "Cumulated double bonds; carbon with two double bonds"
        },
        {
            name:     "E/Z double bond (trans-2-butene)",
            smiles:   "C/C=C/C",
            expected: "valid",
            notes:    "Directional bonds / for E/Z stereochemistry"
        },
        {
            name:     "Cis-2-butene",
            smiles:   "C/C=C\\C",
            expected: "valid",
            notes:    "Mixed / and \\ directional bonds"
        },

        // --- Chirality ---
        {
            name:     "L-Alanine (@ chirality)",
            smiles:   "N[C@@H](C)C(=O)O",
            expected: "valid",
            notes:    "Tetrahedral @@ chirality on bracket atom with explicit H"
        },
        {
            name:     "D-Alanine (@ chirality)",
            smiles:   "N[C@H](C)C(=O)O",
            expected: "valid",
            notes:    "Tetrahedral @ chirality (opposite hand)"
        },

        // --- Disconnected fragments ---
        {
            name:     "Sodium chloride (disconnected)",
            smiles:   "[Na+].[Cl-]",
            expected: "valid",
            notes:    "Dot-separated ionic pair"
        },
        {
            name:     "Ethanol in water",
            smiles:   "CCO.O",
            expected: "valid",
            notes:    "Two disconnected fragments"
        },

        // --- Empty / degenerate ---
        {
            name:     "Single atom (methane)",
            smiles:   "C",
            expected: "valid",
            notes:    "Simplest SMILES; implicit 4 H"
        },
        {
            name:     "Hydrogen molecule",
            smiles:   "[H][H]",
            expected: "valid",
            notes:    "Two bracket hydrogens bonded together"
        },

        // --- Known parser stress tests ---
        {
            name:     "Wildcard atom",
            smiles:   "[*]",
            expected: "valid",
            notes:    "Generic R-group / wildcard"
        },
        {
            name:     "Atom with map number only",
            smiles:   "[C:1]",
            expected: "valid",
            notes:    "Map number without reaction context"
        },
        {
            name:     "Multiple charges (deprecated syntax)",
            smiles:   "[Fe+++]",
            expected: "valid",
            notes:    "Three consecutive + signs = +3 charge (deprecated but valid OpenSMILES)"
        },
        {
            name:     "Quadruple bond (GaAs)",
            smiles:   "[Ga+]$[As-]",
            expected: "error",
            notes:    "BIME does not support '$' quadruple bond -- Unexpected character"
        },
        {
            name:     "Bare wildcard (no brackets)",
            smiles:   "C*C",
            expected: "error",
            notes:    "Wildcard outside brackets is not supported"
        },
        {
            name:     "Unclosed bracket",
            smiles:   "[NH4+",
            expected: "error",
            notes:    "Missing ']' -- parser must report unclosed bracket"
        },
        {
            name:     "Unclosed ring",
            smiles:   "C1CCC",
            expected: "error",
            notes:    "Ring 1 opened but never closed"
        },
        {
            name:     "Unbalanced parentheses",
            smiles:   "C(CC",
            expected: "error",
            notes:    "Branch opened but not closed"
        },
        {
            name:     "Empty SMILES",
            smiles:   "",
            expected: "error",
            notes:    "Empty string should be rejected by validateSmiles"
        }
    ]
};

// ==================================================================
// Self-test runner (works in browser or Node.js)
// ==================================================================
(function() {
    'use strict';

    // Guard: only run when SmilesParser is available
    if (typeof SmilesParser === 'undefined' || typeof Molecule === 'undefined') {
        if (typeof console !== 'undefined') {
            console.log('[test-cases] SmilesParser or Molecule not loaded -- skipping self-test.');
        }
        return;
    }

    var pass = 0, fail = 0, total = 0;

    function runSection(label, cases) {
        console.log('\n=== ' + label + ' ===');
        for (var i = 0; i < cases.length; i++) {
            var tc = cases[i];
            total++;
            var result = SmilesParser.validateSmiles(tc.smiles);
            var expectValid = (tc.expected !== 'error');

            // For molecules and reactions, also check atom/bond counts
            var ok = (result.valid === expectValid);
            var detail = '';

            if (ok && expectValid && tc.atoms !== undefined) {
                var mol = SmilesParser.parse(tc.smiles);
                if (mol.atoms.length !== tc.atoms) {
                    ok = false;
                    detail = ' [atom count: got ' + mol.atoms.length + ', expected ' + tc.atoms + ']';
                }
                if (ok && tc.bonds !== undefined && mol.bonds.length !== tc.bonds) {
                    ok = false;
                    detail = ' [bond count: got ' + mol.bonds.length + ', expected ' + tc.bonds + ']';
                }
            }

            if (ok) {
                pass++;
            } else {
                fail++;
                var status = result.valid ? 'parsed OK' : 'errors: ' + result.errors.join('; ');
                console.log('  FAIL: ' + tc.name + ' (' + tc.smiles + ') -- '
                    + status + detail
                    + (result.warnings.length ? ' [warnings: ' + result.warnings.join('; ') + ']' : ''));
            }
        }
    }

    runSection('Molecules', TEST_CASES.molecules);
    runSection('Reactions', TEST_CASES.reactions);
    runSection('Edge Cases', TEST_CASES.edge_cases);

    console.log('\n--- Results: ' + pass + '/' + total + ' passed, ' + fail + ' failed ---');
})();
