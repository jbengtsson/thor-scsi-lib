# APPLE-II RADIA toolkit v0.2.4

A parameterized APPLE-II reference model for RADIA. Version 0.2.4 preserves the
v0.2.3 modeling and validation workflow and corrects a macOS-only test-harness
PATH defect discovered during target-machine validation. It includes a persistent
four-row device, bounded phase optimization, numerical quality gates, field
integrals, higher-harmonic fitting, and reference-data comparison.

This remains a **generic reference model**, not a certified replica of a specific
ESRF or MAX IV insertion device. Device-faithful conclusions require exact
magnet geometry and material data plus comparison with published or measured
field maps at the relevant gaps and phases.

## What changed from v0.1

- Tests are real pytest functions with positive and failure-case coverage.
- All public numeric inputs are validated for finiteness and applicable bounds.
- Generic end blocks are now contiguous by default. `end_clearance_mm` makes any
  intended longitudinal gap explicit.
- Optimization rejects points below a minimum total fundamental amplitude or
  above a maximum combined relative fit residual.
- `scan_phase()` is a reusable function and `optimize_phase()` performs bounded
  coarse-to-fine refinement.
- `Apple2Device` retains handles to all four rows; `set_phase()` translates the
  existing row objects with RADIA transformations instead of deleting and
  reconstructing the model.
- Harmonic orders are configurable; first and second numerical field integrals
  are reported.
- `apple2-reference-compare` compares calculated metrics with a user-supplied
  CSV reference set.

## Install

```bash
python -c "import radia; print(radia.__file__)"
cd apple2_radia
python -m pip install -e .
```

For tests:

```bash
python -m pip install -e '.[test]'
python -m pytest
```

## Analyze one operating point

```bash
apple2-radia --period 40 --periods 10 --gap 12 --phase-mm 10 \
  --harmonics 1 3 5 --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10
```

Outputs `apple2_field.csv` and `apple2_analysis.json`.

## Scan and optimize phase

```bash
apple2-phase-scan --period 40 --gap 12 \
  --phase-min -20 --phase-max 20 --phase-steps 41 \
  --refinement-levels 3 --refinement-steps 11
```

The device is built once. Each evaluation calls `Apple2Device.set_phase()`, which
applies only the required incremental row translations. The best point must pass
the amplitude and residual gate; normalized circularity alone is insufficient.

## End-block convention

Let `dz` be the regular block length and `edz` the end-block length. The default
end centers are positioned at half the sum of the adjacent lengths from the
nearest regular center, making the surfaces contiguous. A nonzero
`end_clearance_mm` adds a documented surface-to-surface gap at each end.

The magnetization magnitude and length fractions remain generic placeholders.
Replace them with the exact termination sequence before engineering comparison.

## Harmonics and field integrals

The toolkit fits requested harmonic orders independently for Bx and By over the
central analysis interval. Numerical field integrals use the sampled field:

- first integral: `∫ B dz`, in T·mm;
- second integral: `∫(∫B dz) dz`, in T·mm², with the cumulative first integral
  set to zero at the first sample.

## Reference-data comparison

Prepare calculated and reference CSV files with:

```text
gap_mm,phase_mm,Bx1_T,By1_T,phase_deg,circularity
```

Then run:

```bash
apple2-reference-compare calculated.csv reference.csv
```

The tool requires a unique calculated point at every reference `(gap, phase)`
and reports RMSE for Bx, By, wrapped relative phase, and circularity. The included
`examples/reference_template.csv` is an **illustrative schema only**, not validated
published data.

## Interpretation limits

`|s3|`, `s1`, and `s2` are magnetic-field ellipse metrics derived from fitted
fundamental components. They are not photon-radiation Stokes parameters. A high
`|s3|` is accepted only when the field-amplitude and fit-residual criteria pass.

## Remaining external validation

The package must still be run with the user's installed macOS RADIA extension and
compared with authoritative field data. This release contains no claimed ESRF or
MAX IV agreement and no radiation calculation.

## Version 0.2.4 test-harness correction

The failure-path regression test now preserves the host command-search path and
overrides only `PYTHON_BIN`. The previous test replaced `PATH` with the directory
containing `shasum`; on macOS that excluded `/bin/bash`, so the script failed at
its `#!/usr/bin/env bash` launcher with exit status 127 before the intended
missing-pytest path could be exercised. No runtime validation or modeling
algorithm changed.

## Target-Mac validation and evidence capture

On the Intel Mac where `import radia` succeeds, run:

```bash
cd apple2_radia
./scripts/validate_macos_radia.sh
```

The script compiles the exact source tree, runs the tests, executes one real
RADIA operating point, performs the bounded phase optimization, and—after the
log stream is closed—writes exact-byte SHA-256 evidence for every output. `pytest` is a declared
test-only dependency and must be installed in the selected interpreter:

```bash
python3 -m pip install -e '.[test]'
```

When it is missing, the script now stops before any RADIA calculation and prints
that exact corrective command instead of emitting a Python traceback. An explicit
`AUTO_INSTALL_TEST_DEPS=1` opt-in lets the script install the declared test extra.
The test step disables automatic loading of third-party pytest plugins, clears
`PYTEST_ADDOPTS`, prints each test as it runs, and enforces a configurable timeout.
Override the default 120-second timeout with, for example,
`PYTEST_TIMEOUT_SECONDS=300`.
To include reference comparison:

```bash
REFERENCE_CSV=/absolute/path/reference.csv ./scripts/validate_macos_radia.sh
```


The final evidence chain avoids a self-referential log hash:

- `validation.log.sha256` verifies the closed final log;
- `BYTES.txt` records byte lengths for the scientific outputs, `validation.log`,
  and `validation.log.sha256`;
- `SHA256SUMS.txt` uses the standard `digest  filename` format and verifies those
  artifacts plus `BYTES.txt`;
- `SHA256SUMS.txt.sha256` verifies that checksum manifest.

From the generated output directory, verify the chain with:

```bash
shasum -a 256 -c validation.log.sha256
shasum -a 256 -c SHA256SUMS.txt.sha256
shasum -a 256 -c SHA256SUMS.txt
```

Successful script completion establishes only the checks recorded in its log;
device fidelity still depends on the authority and applicability of the supplied
reference data.
