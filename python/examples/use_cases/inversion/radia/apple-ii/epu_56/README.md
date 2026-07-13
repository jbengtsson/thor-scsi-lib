# Daresbury HU56 refined prototype

Run from this directory:

```bash
./generate_epu_56_phase_0.sh
```

Optional output directory:

```bash
./generate_epu_56_phase_0.sh runs/my_phase_0
```

The script uses the source tree directly through the project-local `python/` source directory; installing the
`apple-ii` package is not required.

## Implemented published parameters

- 56 mm magnetic period
- 17 periods
- 21 mm nominal gap
- four blocks per period
- 40 x 40 mm transverse block dimensions
- 1.25 T remanence
- 0.5 mm horizontal gap between adjacent left/right arrays
- 6.95 mm end-block length

The end blocks use full remanence in this prototype. That is an explicit modeling
assumption, not a reported measured end-block magnetization.

## Generated files

- `field.csv`: combined on-axis field
- `analysis.json`: harmonics, field ellipse, quality gate, and field integrals
- `diagnostics.json`: exact generated block geometry and magnetization plus
  independent field analysis for UL, UR, LL, and LR and a superposition check

## Not yet represented

- the published S1/S2/S3 end-spacing topology
- 5 x 5 mm mounting notches
- measured individual block strengths and magnetization-direction errors
- block sorting and local shimming

This output is a refined diagnostic prototype, not a validated engineering
reproduction of the Daresbury device.
