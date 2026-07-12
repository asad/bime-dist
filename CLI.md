# BIME Command-Line Tool

`bime` is a headless command-line wrapper around the same engine that
powers the browser workbench — SMILES parsing, canonical layout,
atom-atom mapping, and image/structure export, all
offline with zero network calls.

- **No build step.** Pure Node, zero runtime dependencies.
- **Node 18+** for the source/npm path; the standalone binaries need no
  Node at all.
- **Apache-2.0**, free.

---

## Table of contents

1. [Install](#install)
2. [Quick start](#quick-start)
3. [Commands](#commands)
   - [`version`](#version)
   - [`parse`](#parse)
   - [`clean`](#clean)
   - [`aam`](#aam)
   - [`export`](#export)
   - [`help`](#help)
4. [Global flags & I/O](#global-flags--io)
5. [How-to recipes](#how-to-recipes)
6. [Scripting & pipelines](#scripting--pipelines)
7. [Exit codes](#exit-codes)

---

## Install

### Option A — standalone binary (no Node required)

Download the binary for your platform from the
[latest release](https://github.com/asad/bime-dist/releases/latest):

| Platform | Asset |
|---|---|
| macOS (Apple Silicon) | `bime-<ver>-macos-arm64` |
| macOS (Intel) | `bime-<ver>-macos-x64` |
| Linux x86-64 | `bime-<ver>-linux-x64` |
| Linux ARM64 | `bime-<ver>-linux-arm64` |
| Windows x64 | `bime-<ver>-win-x64.exe` |

```bash
curl -L -o bime https://github.com/asad/bime-dist/releases/latest/download/bime-<ver>-macos-arm64
chmod +x bime
./bime version
```

Verify the download against the published checksums:

```bash
shasum -a 256 -c bime-<ver>-checksums.txt
```

### Option B — npm (every platform, needs Node 18+)

```bash
npm install -g bime-dist
bime version
```

### Option C — from a source checkout

```bash
git clone https://github.com/asad/bime-dist.git
cd bime-dist
./bin/bime version          # chmod +x bin/bime if needed
```

---

## Quick start

```bash
bime version                                   # print engine versions
bime parse 'CC(=O)Oc1ccccc1C(=O)O'             # diagnose a SMILES (aspirin)
bime clean 'C1=CC=CC=C1' --format json         # canonicalise + layout
bime aam 'CCO>>CC=O' --format json             # atom-atom map a reaction
bime export 'c1ccccc1' --format svg --out benzene.svg
bime help aam                                  # per-command help
```

---

## Commands

### `version`

Print the bundle and module versions (all move in lockstep per release).

```bash
$ bime version
bime 3.0.3
RDT.version         3.0.3
AtomTrace.version   3.0.3
ToolbarPrefs.version 3.0.3
```

### `parse`

SMILES diagnostic — atom/bond counts, charges, stereocentres, element
histogram, reaction detection.

```bash
$ bime parse 'CC(=O)Oc1ccccc1C(=O)O'
SMILES:        CC(=O)Oc1ccccc1C(=O)O
Atoms:         13
Bonds:         13
Charges:       0
Stereo atoms:  0
Reaction:      no
Elements:
  C: 9
  O: 4
```

Add `--format json` for machine-readable output.

### `clean`

Run the 2D layout engine and print the canonical SMILES.

```bash
bime clean 'C1=CC=CC=C1'              # → c1ccccc1
bime clean 'Cc1ccccc1' --format json
```

### `aam`

Atom-atom mapping for a reaction SMILES (`reactants>>products`). Returns
the atom correspondence, bond changes, a confidence score, and a quality
grade (A/B/C/F).

```bash
$ bime aam 'CCO>>CC=O' --format json
{
  "status": "mapped",
  "score": 0.5,
  "confidence": 0.75,
  "quality": "B",
  "mappedCount": 3,
  "mapping": { "7": 10, "8": 11, "9": 12 },
  "bondChanges": [ { "type": "orderChange", ... } ]
}
```

Flags:

| Flag | Effect |
|---|---|
| `--format text\|json\|svg` | output format (default `text`) |
| `--no-chemistry-aware` | disable the chemistry-aware fitness penalty (atom-radii / electronegativity / stereo consistency) |
| `--timeout-ms N` | per-call solver timeout (default 2000) |
| `--keep-mapping` | use the input `:n` atom-map numbers as-is — render/report a curated or externally-supplied mapping **without re-solving** |

#### Publication-quality mapped-reaction figure (`--format svg`)

`bime aam <reaction> --format svg` renders the mapped reaction as a
publication-ready figure: reactants, arrow, and products laid out as a
scheme, with atom-atom map numbers, bond-change highlighting (red cleaved,
green formed, amber order change), sub-fragment trace colouring, and a
per-component compound-name caption. Auto-sized, white background by default.

```bash
bime aam 'CCO.CC(=O)O>>CCOC(=O)C.O' --format svg --out esterification.svg
bime aam 'CCO>>CC=O' --format svg --reaction-center --out oxidation.svg
bime aam 'CCO>>CC=O' --format svg --transparent --no-labels --out bare.svg
bime aam 'CCO>>CC=O' --format svg --ids        # caption with KEGG ids

# Already have a mapped reaction? Render it WITHOUT re-mapping —
# the input :n atom-map numbers are used verbatim:
bime aam '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]' --format svg --keep-mapping --out fig.svg
```

| Flag | Effect |
|---|---|
| `--screen` | on-screen look instead of the print preset |
| `--transparent` | transparent background (drop onto any page) |
| `--reaction-center` | add dashed reaction-centre rings |
| `--no-map-numbers` | omit the atom-atom map numbers |
| `--ids` | caption components with the KEGG id instead of the name |
| `--no-labels` | omit the per-component compound-name captions |
| `--width N` `--height N` | pin the canvas size (default: auto-fit) |
| `--out file` | write to a file (default stdout) |
| `--no-stamp` | omit the BIME provenance stamp |

### `export`

Export a structure to SVG, MOL (V2000), SDF, RXN, or canonical SMILES.

```bash
bime export 'c1ccccc1' --format svg --width 800 --height 600 --out benzene.svg
bime export 'CCO' --format mol --out ethanol.mol
bime export 'c1ccccc1C' --format smiles --no-stamp     # clean SMILES, no provenance header
```

Flags:

| Flag | Effect |
|---|---|
| `--format svg\|mol\|sdf\|rxn\|smiles` | output format |
| `--width N` `--height N` | SVG canvas size (SVG only) |
| `--out file` | write to a file (default stdout) |
| `--no-stamp` | omit the BIME provenance stamp (version + SHA-384 fingerprint + Apache-2.0 line embedded by default) |

### `help`

```bash
bime help            # top-level command list
bime help aam        # per-command flags + usage
```

---

## Global flags & I/O

| Flag | Meaning |
|---|---|
| `--in <file\|->` | read input from a file, or `-` for stdin |
| `--out <file>` | write output to a file (default stdout) |
| `--format <text\|json>` | output format where applicable |

```bash
echo 'CCO>>CC=O' | bime aam --in - --format json
bime parse --in molecule.smi
```

---

## How-to recipes

**Batch-canonicalise a list of SMILES**

```bash
while read smi; do bime clean "$smi"; done < input.smi > canonical.smi
```

**Render a folder of SMILES to SVG**

```bash
i=0
while read smi; do
  bime export "$smi" --format svg --out "mol_$((i++)).svg"
done < library.smi
```

**Grade a batch of reactions by AAM quality**

```bash
while read rxn; do
  q=$(bime aam "$rxn" --format json | grep -o '"quality": "[A-F]"')
  echo "$q  $rxn"
done < reactions.rxn
```

---

## Scripting & pipelines

- All commands are deterministic — same input → same output, byte for
  byte. Safe to cache and diff.
- `--format json` output is stable and suitable for piping into `jq`:

  ```bash
  bime aam 'CCO>>CC=O' --format json | jq '.quality, .mappedCount'
  ```
- Everything runs offline. No telemetry, no network calls.
- The provenance stamp on exports is opt-out (`--no-stamp`) when you need
  byte-clean output for downstream tooling.

---

## Exit codes

| Code | Meaning |
|---|---|
| `0` | success |
| `1` | bad input (invalid SMILES, malformed flags) |
| `2` | usage error (unknown command/flag) |

---

For the browser editor and the full API, see
[`USER_GUIDE.md`](USER_GUIDE.md) and the live
[docs](https://asad.github.io/bime-dist/docs.html).
