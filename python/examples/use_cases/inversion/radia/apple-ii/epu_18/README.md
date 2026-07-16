# AQUA APPLE-X 18 mm RADIA model — Cartesian kick-map grid

The second-order map is computed with `rad.FldFocKickPer` on one regular
Cartesian grid. The default bounds are derived from the model's rotated-square
magnetic opening:

```text
|x| + |y| <= gap_mm + inner_chamfer_mm
```

At the nominal 1.5 mm gap this gives axis vertices at approximately
`x,y = +/-3.889087 mm`. The rectangular CSV includes every grid point and flags
whether each point is inside the physical rotated-square aperture.

Nominal run:

```bash
./epu_18_gen.sh runs/production_cartesian
```

Explicit Figure-7-like bounds can be requested with:

```bash
./epu_18_gen.sh runs/cartesian_pm4 \
  --field-x-min-mm -4 --field-x-max-mm 4 --field-x-points 61 \
  --field-y-min-mm -4 --field-y-max-mm 4 --field-y-points 61
```

The central symplectic polynomial fit uses all Cartesian points satisfying
`r <= --kick-fit-r-max-mm`; it does not use a radial output subset.

## RADIA compatibility

`FldFocKickPer` is called with formatted output mode `"fix"`. The `"tab"` mode
caused native heap corruption in the tested macOS Python binding. Numerical
kick matrices are separate return elements.
