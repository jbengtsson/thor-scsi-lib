# Validation scope

This source release is designed for two validation layers:

1. **Environment-independent:** Python compilation, package build/install, unit tests,
   synthetic harmonic recovery, geometry bookkeeping with a RADIA test double,
   optimizer behavior, field-integral numerics, and reference comparison.
2. **Environment-dependent:** import and execution against the user's real macOS RADIA
   extension, actual field calculations, and comparison with authoritative device data.

The first layer can be run in a standard Python environment. The second layer is not
established until executed on the target Mac with the installed RADIA binary.


## macOS pytest isolation

The validation script executes the repository tests in a subprocess with
`PYTEST_DISABLE_PLUGIN_AUTOLOAD=1` and an empty `PYTEST_ADDOPTS`. The project tests
require no third-party pytest plugins. Per-test output is enabled, and the test
subprocess is terminated with exit status 124 if it exceeds
`PYTEST_TIMEOUT_SECONDS` (default: 120 seconds).

## Test dependency preflight

`pytest>=8` is declared in the `test` optional dependency group. The macOS
validation script checks that `pytest` is importable by the selected
`PYTHON_BIN` before importing it for version reporting. If it is absent, the
script exits with status 2, prints `python -m pip install -e '.[test]'`, and does
not begin a RADIA calculation. Environment mutation is opt-in only through
`AUTO_INSTALL_TEST_DEPS=1`.



## Post-log evidence finalization

Version 0.2.3 does not place the final hash of `validation.log` inside the log
itself. The validation body and its terminal status are first written through
`tee`; the pipeline then closes. Only after closure does the parent shell write:

1. `validation.log.sha256`, containing the exact-byte digest of the stable log;
2. `BYTES.txt`, containing byte lengths for the output artifacts available when
   evidence finalization begins;
3. `SHA256SUMS.txt`, using the standard `digest  filename` syntax accepted by
   `shasum -c`, and covering the outputs plus `BYTES.txt`;
4. `SHA256SUMS.txt.sha256`, containing the exact-byte digest of the checksum manifest.

The script returns the original validation status after producing this evidence,
so evidence finalization does not turn a failed validation into a successful one.

## Failure-path test environment

Version 0.2.4 corrects the regression-test environment used to exercise missing
`pytest`. The test inherits the normal host `PATH`, ensuring `/usr/bin/env bash`
can resolve the shell on macOS, and changes only `PYTHON_BIN` plus the explicit
non-mutating `AUTO_INSTALL_TEST_DEPS=0` setting. This tests the intended exit-2
preflight path without introducing an unrelated launcher failure.


## Version 0.2.5 manifest-format correction

Version 0.2.4 mixed a byte-length column into `SHA256SUMS.txt`. Although the
digests were correct, macOS `shasum -c` interpreted the size field as part of
the filename. Version 0.2.5 separates sizes into `BYTES.txt` and restores the
standard two-field checksum grammar.

## Version 0.3.0 Daresbury prototype validation

Environment-independent validation additionally covers:

- explicit horizontal inter-array gap placement;
- explicit end-block length override;
- exact geometry and magnetization records;
- per-row field analysis and combined-field superposition using a RADIA test double;
- source-only CLI argument exposure for the new parameters.

The refined HU56 prototype must still be executed against the user's real RADIA
extension. Agreement with the published Daresbury field map is not established by
the source tests or synthetic diagnostic smoke test.
