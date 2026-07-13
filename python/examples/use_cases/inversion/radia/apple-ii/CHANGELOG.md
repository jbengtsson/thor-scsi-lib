# Changelog

## 0.3.0

- Added an explicit horizontal inter-array gap parameter and CLI option.
- Added an explicit end-block length parameter that overrides the fractional default.
- Added exact per-block geometry and magnetization records.
- Added optional per-row field and superposition diagnostics through `--diagnostics-json`.
- Added a refined Daresbury HU56 phase-zero prototype script using the published
  56 mm period, 17 periods, 21 mm gap, 40 x 40 mm blocks, 0.5 mm adjacent-array
  gap, and 6.95 mm end-block length.
- Kept S1/S2/S3 end spacing, mounting notches, and measured block errors explicitly
  outside the prototype scope.


## 0.2.5

- Corrected `SHA256SUMS.txt` to use the standard `digest  filename` format accepted by `shasum -c`.
- Moved byte lengths into a separate, checksummed `BYTES.txt` ledger.
- Added a regression assertion that runs `shasum -a 256 -c SHA256SUMS.txt` directly.
- Left all APPLE-II geometry, field, analysis, optimization, and reference algorithms unchanged.

## 0.2.4

- Corrected the macOS failure-path regression test so it preserves the host `PATH`.
- The test now overrides only `PYTHON_BIN` and `AUTO_INSTALL_TEST_DEPS`.
- Eliminated the false exit-127 failure caused by excluding `/bin/bash` from `PATH`.
- Left the validation script and all APPLE-II modeling algorithms unchanged.

## 0.2.3

- Closed `validation.log` before calculating its digest.
- Moved the final log digest to `validation.log.sha256`, eliminating self-reference.
- Added a post-log `SHA256SUMS.txt` covering outputs, the final log, and its sidecar.
- Added `SHA256SUMS.txt.sha256` so the final manifest can be verified independently.
- Preserved validation exit status while producing evidence for successful or failed runs.
- Synchronized `apple_ii.__version__` with the 0.2.3 project metadata.
- Left all APPLE-II geometry, analysis, optimization, and reference algorithms unchanged.

## 0.2.2

- Added a fail-fast dependency preflight for the declared pytest test extra.
- Replaced the missing-pytest traceback with an exact installation command.
- Added opt-in `AUTO_INSTALL_TEST_DEPS=1` installation; the default remains non-mutating.
- Ensured no RADIA calculation starts when the test runner is unavailable.

## 0.2.1

- Hardened the macOS validation test step against host-specific pytest plugin stalls.
- Disabled automatic third-party pytest plugin loading during validation.
- Cleared inherited `PYTEST_ADDOPTS`, enabled per-test progress output, and added a configurable 120-second timeout.
- Added explicit pytest/NumPy preflight version reporting and actionable timeout diagnostics.

## 0.2.0

- Replaced module-scope assertions with a collected pytest suite.
- Added strict finite/range validation and explicit CLI argument constraints.
- Corrected generic end-block placement and added explicit end clearance.
- Added amplitude/residual quality gating and clipped normalized metrics.
- Added reusable scan and bounded coarse-to-fine phase optimization.
- Added persistent `Apple2Device` and `RowHandle` objects with in-place RADIA row translation.
- Added higher-harmonic fitting, numerical first/second field integrals, and CSV reference comparison.
- Expanded all scan geometry and analysis controls.

## 0.1.0

- Initial reusable four-row reference geometry, fundamental fit, ellipse metrics, and uniform phase scan.
