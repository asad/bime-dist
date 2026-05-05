# SDG Layout Benchmarks

Reference artefacts for the Structure-Diagram-Generator (SDG) layout
completion work that started in v1.8.12 / v1.8.13.

## Benchmark molecule: bisazo Direct dye

A mono-sulfonated bisazo dye, C32H22N6O3S, MW ≈ 570.6 — representative
of the elongated chromophore class that BIME's default rule-based
layout currently mis-renders.

```
SMILES: c5(cc6ccc(\N=N\c1cc2ccc(cc2cc1N)c4ccc(\N=N\c3ccccc3)cc4)cc6cc5N)S(=O)(=O)O

Topology: [phenyl]-N=N-[phenylene]-[aminonaphthyl]-N=N-[aminonaphthyl-SO3H]
```

## Files

- `bisazo-dye-reference.svg` — hand-placed reference layout. Every
  ring is a regular hexagon (bond-length CV 0 %, max interior-angle
  deviation 0°). Aspect ratio 3.63 : 1.
- `bisazo-dye-reference-coords.json` — atom-id → (x, y) in BL units
  (BL = 1.0). Use as a template / target layout for the SDG refiner.

## Acceptance gate

Pinned in `tests/test_v1_8_x_chromophore_layout.js`. Two tiers:

| Metric                       | Current (locked) | SDG target |
|------------------------------|------------------|------------|
| Ring bond-length spread      | ≤ 20 %           | ≤ 5 %      |
| Ring interior-angle dev      | ≤ 90°            | ≤ 5°       |
| Layout aspect ratio          | ≥ 1.5            | ≥ 3.0      |
| All 6-rings regular hexagons | n/a              | yes        |

The TARGET tests are advisory — they log the gap rather than fail the
suite — until the SDG refiner is on by default. Switch them from
advisory to hard once the refiner lands.
