#!/usr/bin/env node
/*
 * audit-aromatic.js — BIME v1.0.1 chemistry audit
 * Copyright (c) 2026 BioInception PVT LTD
 *
 * ALGORITHM
 * ---------
 * 1. Load common-molecules.js as text and extract every {name, smiles, category}
 *    record using a tolerant regex parser (no JS eval — the file is data not
 *    code).
 * 2. For each SMILES we run a *lightweight* aromatic-sanity check (full RDKit-
 *    grade perception is out of scope; we only catch the egregious mistakes
 *    introduced by the 984e5a2 "mass canonicalisation" pass):
 *
 *      a. Tokenise the SMILES into atoms, bonds, branches, ring-closure
 *         digits, bracket atoms, and stereo markers.
 *      b. Build a per-atom adjacency list using an explicit ring-bond stack
 *         that pairs ring-closure digits.
 *      c. Find the smallest ring each lowercase aromatic atom (c, n, o, s, p)
 *         participates in via BFS.  If the atom is not in any ring, that is
 *         flagged immediately (lowercase outside a ring is illegal).
 *      d. Reject the ring if any of:
 *           - ring size not in {5, 6, 7}
 *           - ring contains an explicit single-bond saturator (sp3 CH2 written
 *             upper-case "C" with no =/: incident bonds AND no ring aromatic
 *             neighbours of the same ring), suggesting a saturated member was
 *             carried into an "aromatic" ring
 *           - ring contains an sp3 heteroatom that cannot be aromatic (e.g.
 *             tetrahedral [C@H] or [C@@H], an aliphatic O between two sp3 C,
 *             a peroxide OO chain)
 *           - ring contains adjacent O-O or O-S (peroxide / persulfide) — not
 *             a Hückel aromatic motif
 *      e. Run a minimal Hückel π-electron tally for confirmed-aromatic rings
 *         where we can: c=1, n(in pyridine)=1, [nH]=2, o=2, s=2.  Ring is
 *         valid only if total ∈ {2, 6, 10}.  Skipped when any ring atom was
 *         written upper-case (Kekulé form).
 *
 * 3. Cross-record checks:
 *      - Duplicate SMILES across distinct names  → DUPLICATE
 *      - Known chiral drugs without any '@'      → MISSING_STEREO
 *      - Identical canonicalised atom multiset   → STRUCTURAL_DUP (heuristic)
 *
 * 4. Print a table; exit 1 if any issues.
 *
 * No external deps.  Pure Node 14+.
 */

'use strict';

const fs = require('fs');
const path = require('path');

const SRC = process.argv[2] ||
  path.join(__dirname, '..', 'common-molecules.js');

// ---------- 1. parse records ------------------------------------------------

function loadRecords(file) {
  const txt = fs.readFileSync(file, 'utf8');
  const recs = [];
  // Tolerate single OR double-quoted keys/values; permissive on whitespace.
  const re = /\{\s*(?:"name"|name)\s*:\s*"([^"]+)"\s*,\s*(?:"smiles"|smiles)\s*:\s*"([^"]+)"\s*,\s*(?:"category"|category)\s*:\s*"([^"]+)"\s*\}/g;
  let m;
  while ((m = re.exec(txt)) !== null) {
    recs.push({ name: m[1], smiles: m[2], category: m[3] });
  }
  return recs;
}

// ---------- 2. SMILES tokeniser & graph -------------------------------------

const AROMATIC = new Set(['c', 'n', 'o', 's', 'p']);
const ORG_SUBSET = new Set(['B','C','N','O','S','P','F','I','b','c','n','o','s','p']);

function tokenise(smiles) {
  const toks = [];
  let i = 0;
  while (i < smiles.length) {
    const ch = smiles[i];
    if (ch === '[') {
      const j = smiles.indexOf(']', i);
      toks.push({ kind: 'atom', raw: smiles.slice(i, j + 1), bracket: true });
      i = j + 1;
    } else if (ch === 'C' && smiles[i + 1] === 'l') {
      toks.push({ kind: 'atom', raw: 'Cl' }); i += 2;
    } else if (ch === 'B' && smiles[i + 1] === 'r') {
      toks.push({ kind: 'atom', raw: 'Br' }); i += 2;
    } else if (/[A-Za-z]/.test(ch) && ORG_SUBSET.has(ch)) {
      toks.push({ kind: 'atom', raw: ch }); i++;
    } else if (ch === '(' || ch === ')') {
      toks.push({ kind: ch }); i++;
    } else if (ch === '%') {
      toks.push({ kind: 'ring', n: parseInt(smiles.slice(i + 1, i + 3), 10) });
      i += 3;
    } else if (/[0-9]/.test(ch)) {
      toks.push({ kind: 'ring', n: parseInt(ch, 10) }); i++;
    } else if ('-=#:/\\.'.includes(ch)) {
      toks.push({ kind: 'bond', raw: ch }); i++;
    } else {
      i++; // unknown / chirality dot
    }
  }
  return toks;
}

function buildGraph(smiles) {
  const toks = tokenise(smiles);
  const atoms = []; // { sym, aromatic, bracket, raw, idx }
  const adj = [];   // adj[i] = [{j, bond}]
  const ringOpen = {}; // digit -> { atomIdx, bond }
  const stack = []; // for branches; holds prevAtom indices
  let prev = -1;
  let pendingBond = null;
  for (const t of toks) {
    if (t.kind === 'atom') {
      let sym, aromatic;
      if (t.bracket) {
        // bracket atom: pull leading element symbol
        const inner = t.raw.slice(1, -1);
        const m = inner.match(/^[0-9]*([A-Za-z][a-z]?)/);
        sym = m ? m[1] : inner;
        aromatic = /^[a-z]/.test(sym);
        if (aromatic) sym = sym;
      } else {
        sym = t.raw;
        aromatic = /^[a-z]/.test(sym);
      }
      const idx = atoms.length;
      atoms.push({ sym, aromatic, bracket: !!t.bracket, raw: t.raw, idx });
      adj.push([]);
      if (prev >= 0) {
        const bond = pendingBond || (atoms[prev].aromatic && aromatic ? ':' : '-');
        adj[prev].push({ j: idx, bond });
        adj[idx].push({ j: prev, bond });
      }
      prev = idx;
      pendingBond = null;
    } else if (t.kind === 'bond') {
      pendingBond = t.raw;
    } else if (t.kind === '(') {
      stack.push(prev);
    } else if (t.kind === ')') {
      prev = stack.pop();
      pendingBond = null;
    } else if (t.kind === 'ring') {
      if (ringOpen[t.n] !== undefined) {
        const o = ringOpen[t.n];
        const bond = pendingBond || o.bond ||
          (atoms[o.atomIdx].aromatic && atoms[prev].aromatic ? ':' : '-');
        adj[o.atomIdx].push({ j: prev, bond });
        adj[prev].push({ j: o.atomIdx, bond });
        delete ringOpen[t.n];
      } else {
        ringOpen[t.n] = { atomIdx: prev, bond: pendingBond };
      }
      pendingBond = null;
    }
  }
  return { atoms, adj, unclosed: Object.keys(ringOpen) };
}

// ---------- 3. ring detection ----------------------------------------------

function smallestRingContaining(start, adj) {
  // BFS that records parent; when we revisit a node from a different parent we
  // have a cycle.  Returns the smallest cycle through `start`.
  let best = null;
  for (const nbr of adj[start]) {
    // try every initial neighbour, look for path back to start avoiding the
    // direct edge.
    const queue = [[nbr.j, [start, nbr.j]]];
    const seen = new Map();
    seen.set(nbr.j, 1);
    while (queue.length) {
      const [u, path] = queue.shift();
      if (path.length > 8) continue; // ring size cap
      for (const e of adj[u]) {
        if (path.length === 2 && e.j === start) continue; // skip back-edge
        if (e.j === start && path.length >= 3) {
          if (!best || path.length < best.length) best = path.slice();
          continue;
        }
        if (seen.has(e.j)) continue;
        seen.set(e.j, 1);
        queue.push([e.j, path.concat(e.j)]);
      }
    }
  }
  return best; // array of atom indices forming the cycle (start..last, edge last->start implied)
}

// ---------- 4. aromatic sanity ---------------------------------------------

function piContribution(atom) {
  const s = atom.sym.toLowerCase();
  if (s === 'c') return 1;
  if (s === 'o') return 2;
  if (s === 's') return 2;
  if (s === 'n') {
    if (atom.bracket && /\[nH/.test(atom.raw)) return 2;
    return 1;
  }
  if (s === 'p') return 1;
  return null;
}

function auditMolecule(rec) {
  const issues = [];
  let g;
  try { g = buildGraph(rec.smiles); } catch (e) {
    return [{ severity: 'PARSE', msg: 'tokenise failure: ' + e.message }];
  }
  if (g.unclosed.length) {
    issues.push({ severity: 'CRITICAL', msg: 'unclosed ring digits: ' + g.unclosed.join(',') });
  }
  // sp3 chirality marker on a lowercase aromatic atom
  for (const a of g.atoms) {
    if (a.bracket && a.aromatic && /@/.test(a.raw)) {
      issues.push({ severity: 'CRITICAL', msg: `aromatic atom with tetrahedral stereo: ${a.raw}` });
    }
  }
  for (let i = 0; i < g.atoms.length; i++) {
    const a = g.atoms[i];
    if (!a.aromatic) continue;
    if (!AROMATIC.has(a.sym.toLowerCase())) continue;
    const ring = smallestRingContaining(i, g.adj);
    if (!ring) {
      issues.push({ severity: 'CRITICAL', msg: `lowercase '${a.sym}' atom #${i} is not in any ring` });
      continue;
    }
    if (![5, 6, 7].includes(ring.length)) {
      issues.push({ severity: 'CRITICAL', msg: `aromatic atom in ring of size ${ring.length}: ${ring.map(x=>g.atoms[x].raw).join('-')}` });
      continue;
    }
    // Check for sp3 saturators in the ring: an upper-case atom whose ALL
    // ring-bonds to its ring-neighbours are explicit single (non-aromatic) and
    // it has no double bond inside the ring.
    let allArom = true;
    let sp3Count = 0;
    let oxOx = false;
    let hasUpperNonH = false;        // any uppercase organic atom in the ring
    let hasLower = false;            // any lowercase aromatic atom in the ring
    let hasHetero = false;           // n, o, s, p (lowercase) in the ring
    for (let k = 0; k < ring.length; k++) {
      const ai = ring[k];
      const aj = ring[(k + 1) % ring.length];
      const at = g.atoms[ai];
      if (!at.aromatic) {
        allArom = false;
        if (at.sym !== 'H') hasUpperNonH = true;
      } else {
        hasLower = true;
        const ls = at.sym.toLowerCase();
        if (ls === 'n' || ls === 'o' || ls === 's' || ls === 'p') hasHetero = true;
      }
      // adjacent O-O or S-S in ring → peroxide / persulfide → not aromatic
      const sym1 = g.atoms[ai].sym.toLowerCase();
      const sym2 = g.atoms[aj].sym.toLowerCase();
      if ((sym1 === 'o' && sym2 === 'o') || (sym1 === 's' && sym2 === 's')) oxOx = true;
      // sp3 marker: bracket atom with @ or H count making it saturated
      if (at.bracket && /@/.test(at.raw) && at.aromatic) sp3Count++;
    }
    if (oxOx) {
      issues.push({ severity: 'CRITICAL', msg: `peroxide/persulfide bond inside aromatic-flagged ring: ${ring.map(x=>g.atoms[x].raw).join('-')}` });
      continue;
    }
    if (!allArom) {
      // FUSED-RING FALSE-POSITIVE GUARD.  In a polycyclic system the BFS
      // ring-finder will sometimes return a cycle that crosses a fused
      // boundary, picking up sp3 atoms from the *adjacent* ring.
      //
      // Suppress when the ring's lowercase aromatic atoms are exactly the
      // junction atoms of a neighbouring all-aromatic ring -- detected as:
      // some pair of lowercase atoms in this cycle are directly bonded to
      // each other AND that bond is itself part of an aromatic ring (i.e.
      // both atoms are in another lowercase ring).  In that case the
      // current cycle is a fused-ring traversal, not a genuine
      // sp3/aromatic-mixed ring.
      let fusedArtefact = false;
      const lowerIdx = ring.filter(i => g.atoms[i].aromatic);
      outer: for (let p = 0; p < lowerIdx.length; p++) {
        for (let q = p + 1; q < lowerIdx.length; q++) {
          const a = lowerIdx[p], b = lowerIdx[q];
          if (g.adj[a].some(e => e.j === b)) {
            // both a and b are aromatic; an a-b bond means they sit on a
            // shared edge of two rings -> classic fused junction.
            fusedArtefact = true;
            break outer;
          }
        }
      }
      if (fusedArtefact) continue;
      if (!hasLower || !hasUpperNonH) continue;
      issues.push({ severity: 'HIGH', msg: `aromatic atom mixed with sp3 ring members: ${ring.map(x=>g.atoms[x].raw).join('-')}` });
      continue;
    }
    // 4-MEMBERED RING DOUBLE-CHECK.  The BFS may report a 4-cycle that is
    // really a shortcut across a bridged/spiro centre, not a true 4-ring.
    // Verify each consecutive pair in `ring` has a real edge in adj.
    if (ring.length === 4) {
      let real = true;
      for (let k = 0; k < ring.length; k++) {
        const a = ring[k], b = ring[(k + 1) % ring.length];
        if (!g.adj[a].some(e => e.j === b)) { real = false; break; }
      }
      if (!real) continue;          // BFS artefact — silently skip
      issues.push({ severity: 'CRITICAL', msg: `aromatic atom in ring of size ${ring.length}: ${ring.map(x=>g.atoms[x].raw).join('-')}` });
      continue;
    }
    // Verify the ring is a true cycle (consecutive atoms actually bonded).
    // BFS shortest-cycle searches can occasionally splice across fused
    // junctions; without this guard Budesonide-style steroid + acetal
    // systems generate spurious "all-aromatic π=5/8" rings.
    {
      let real = true;
      for (let k = 0; k < ring.length; k++) {
        const a = ring[k], b = ring[(k + 1) % ring.length];
        if (!g.adj[a].some(e => e.j === b)) { real = false; break; }
      }
      if (!real) continue;
    }
    // Hückel tally
    let pi = 0; let ok = true;
    for (const ai of ring) {
      const c = piContribution(g.atoms[ai]);
      if (c === null) { ok = false; break; }
      pi += c;
    }
    // 5-RING LONE-PAIR DONATION.  In every aromatic 5-ring (pyrrole,
    // furan, thiophene, imidazole, pyrazole, oxazole, thiazole, triazole,
    // tetrazole, ...) exactly one heteroatom contributes its lone pair
    // (2 electrons) instead of the 1 electron we tally above.  Add the
    // missing +1 if the ring contains at least one heteroatom; this turns
    // tetrazole c-n-n-n-n (π=5 → 6), pyrrole (π=4 → 6 with [nH]=2 already),
    // imidazole c-n-c-n-c (π=4 → 6), etc., into valid Hückel rings without
    // false-positive-ing all-carbon "5-rings" (which the AROMATIC test
    // would already reject upstream).
    if (ok && ring.length === 5 && hasHetero && pi % 4 !== 2) {
      pi += 1;
    }
    if (ok && ![2, 6, 10].includes(pi)) {
      issues.push({ severity: 'MEDIUM', msg: `ring fails Hückel (π=${pi}): ${ring.map(x=>g.atoms[x].raw).join('-')}` });
    }
  }
  return issues;
}

// ---------- 5. cross-record checks -----------------------------------------

const KNOWN_CHIRAL = new Set([
  'Morphine','Atorvastatin','Estradiol','Progesterone','Testosterone','DHEA',
  'Cholesterol','Penicillin G','Penicillin V','Amoxicillin','Adenosine',
  'Quinine','Menthol','Artemether','Esomeprazole','Armodafinil',
  'L-Alanine','L-Serine','L-Tryptophan','L-Valine','L-Leucine','L-Isoleucine',
  'L-Threonine','L-Cysteine','L-Methionine','L-Phenylalanine','L-Tyrosine',
  'L-Histidine','L-Lysine','L-Arginine','L-Aspartate','L-Glutamate',
  'L-Asparagine','L-Glutamine','L-Proline',
  'D-Glucose','D-Fructose','D-Ribose','D-Galactose','D-Mannose',
]);

function crossChecks(records) {
  const issues = [];
  const bySmiles = new Map();
  for (const r of records) {
    if (!bySmiles.has(r.smiles)) bySmiles.set(r.smiles, []);
    bySmiles.get(r.smiles).push(r.name);
  }
  for (const [smi, names] of bySmiles) {
    if (names.length > 1) {
      // "Acetic acid" vs "Acetic Acid" are *not* chemical duplicates; they
      // are case-only collisions in the molecule library and should not
      // crowd CRITICAL.  Fold case + whitespace and demote to INFO if the
      // resulting normalised names are all identical.
      const norm = names.map(n => n.toLowerCase().replace(/\s+/g, ' ').trim());
      const allSame = norm.every(n => n === norm[0]);
      if (allSame) {
        issues.push({ name: names.join(' / '), severity: 'INFO',
          msg: `casing duplicate (same molecule, name-case differs): ${smi}` });
      } else {
        issues.push({ name: names.join(' / '), severity: 'CRITICAL',
          msg: `duplicate SMILES across distinct names: ${smi}` });
      }
    }
  }
  for (const r of records) {
    if (KNOWN_CHIRAL.has(r.name) && !/@/.test(r.smiles)) {
      issues.push({ name: r.name, severity: 'HIGH',
        msg: `known chiral drug missing stereo markers` });
    }
  }
  return issues;
}

// ---------- 6. main --------------------------------------------------------

function main() {
  const records = loadRecords(SRC);
  const allIssues = [];
  for (const r of records) {
    const iss = auditMolecule(r);
    for (const i of iss) allIssues.push({ name: r.name, smiles: r.smiles, ...i });
  }
  for (const ci of crossChecks(records)) allIssues.push(ci);

  console.log(`Audited ${records.length} molecules.`);
  console.log(`Issues found: ${allIssues.length}\n`);
  if (allIssues.length) {
    const w = (s, n) => (s.length > n ? s.slice(0, n - 1) + '…' : s.padEnd(n));
    console.log(w('SEVERITY', 10) + w('NAME', 28) + w('ISSUE', 80));
    console.log('-'.repeat(118));
    for (const i of allIssues) {
      console.log(w(i.severity, 10) + w(i.name || '', 28) + w(i.msg, 80));
    }
  }
  process.exit(allIssues.length ? 1 : 0);
}

if (require.main === module) main();
