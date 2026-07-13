#!/usr/bin/env bash
# Run on the target Intel Mac with the Python interpreter that imports RADIA.
set -uo pipefail

export LC_ALL=C

ROOT=$(cd "$(dirname "$0")/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}
PYTEST_TIMEOUT_SECONDS=${PYTEST_TIMEOUT_SECONDS:-120}
AUTO_INSTALL_TEST_DEPS=${AUTO_INSTALL_TEST_DEPS:-0}
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
OUT=${1:-"$ROOT/validation_macos_$STAMP"}
REFERENCE_CSV=${REFERENCE_CSV:-}
LOG="$OUT/validation.log"
LOG_SIDECAR="$OUT/validation.log.sha256"
MANIFEST="$OUT/SHA256SUMS.txt"
MANIFEST_SIDECAR="$OUT/SHA256SUMS.txt.sha256"
mkdir -p "$OUT"
rm -f "$LOG" "$LOG_SIDECAR" "$MANIFEST" "$MANIFEST_SIDECAR"

sha256_digest() {
  shasum -a 256 "$1" | awk '{print $1}'
}

write_sha256_sidecar() {
  local target=$1
  local sidecar=$2
  printf '%s  %s\n' "$(sha256_digest "$target")" "$(basename "$target")" > "$sidecar"
}

write_output_manifest() {
  local path name size digest
  : > "$MANIFEST"
  for path in "$OUT"/*; do
    [[ -f "$path" ]] || continue
    name=$(basename "$path")
    case "$name" in
      SHA256SUMS.txt|SHA256SUMS.txt.sha256) continue ;;
    esac
    size=$(wc -c < "$path" | tr -d '[:space:]')
    digest=$(sha256_digest "$path")
    printf '%s  %12d  %s\n' "$digest" "$size" "$name" >> "$MANIFEST"
  done
  write_sha256_sidecar "$MANIFEST" "$MANIFEST_SIDECAR"
}

validation_body() {
  printf 'APPLE-II RADIA macOS validation\n'
  printf 'UTC stamp: %s\n' "$STAMP"
  printf 'source root: %s\n' "$ROOT"
  printf 'output directory: %s\n' "$OUT"
  printf 'python command: %s\n' "$PYTHON_BIN"
  printf 'pytest timeout: %s seconds\n' "$PYTEST_TIMEOUT_SECONDS"
  printf 'auto-install test dependencies: %s\n' "$AUTO_INSTALL_TEST_DEPS"

  "$PYTHON_BIN" --version

  # Pytest is a test-only dependency. Detect it before the general Python preflight so
  # a missing package produces an actionable error instead of a traceback.
  if ! "$PYTHON_BIN" -c 'import pytest' >/dev/null 2>&1; then
    printf 'pytest: NOT INSTALLED for %s\n' "$PYTHON_BIN"
    if [[ "$AUTO_INSTALL_TEST_DEPS" == "1" ]]; then
      printf 'Installing the project test extra into the selected Python environment...\n'
      "$PYTHON_BIN" -m pip install -e "${ROOT}[test]"
    else
      cat >&2 <<EOF2
ERROR: pytest is required for validation step [2/4], but it is not installed
for the selected interpreter: $PYTHON_BIN

Install the declared test dependencies with:

  cd "$ROOT"
  "$PYTHON_BIN" -m pip install -e '.[test]'

Then rerun this script. To allow this script to perform that installation
explicitly, rerun with:

  AUTO_INSTALL_TEST_DEPS=1 "$0"${1:+ "$1"}

No RADIA calculation was attempted.
EOF2
      return 2
    fi
  fi

  # Verify that the optional installation actually supplied pytest.
  if ! "$PYTHON_BIN" -c 'import pytest' >/dev/null 2>&1; then
    printf 'ERROR: pytest remains unavailable after dependency installation.\n' >&2
    return 2
  fi

  "$PYTHON_BIN" - <<'PY'
import hashlib
import pathlib
import platform

import numpy
import pytest
import radia

path = pathlib.Path(radia.__file__).resolve()
print(f"platform: {platform.platform()}")
print(f"machine: {platform.machine()}")
print(f"NumPy: {numpy.__version__}")
print(f"pytest: {pytest.__version__}")
print(f"radia path: {path}")
print(f"radia SHA-256: {hashlib.sha256(path.read_bytes()).hexdigest()}")
PY

  printf '\n[1/4] Compile exact source tree\n'
  "$PYTHON_BIN" -m compileall -q "$ROOT/apple2" "$ROOT/tests"

  printf '\n[2/4] Run isolated source tests with progress\n'
  "$PYTHON_BIN" - "$ROOT" "$PYTEST_TIMEOUT_SECONDS" <<'PY'
from __future__ import annotations

import os
from pathlib import Path
import shlex
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
try:
    timeout = float(sys.argv[2])
except ValueError as exc:
    raise SystemExit("PYTEST_TIMEOUT_SECONDS must be a positive number") from exc
if not timeout > 0:
    raise SystemExit("PYTEST_TIMEOUT_SECONDS must be a positive number")

env = os.environ.copy()
env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
env["PYTEST_ADDOPTS"] = ""
env["PYTHONPATH"] = str(root) + (
    os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
)
command = [
    sys.executable,
    "-m",
    "pytest",
    "-vv",
    "-s",
    "--maxfail=1",
    "--strict-config",
    "--strict-markers",
]
print("pytest plugin autoload: disabled", flush=True)
print("pytest command:", " ".join(shlex.quote(part) for part in command), flush=True)
try:
    result = subprocess.run(command, cwd=root, env=env, timeout=timeout, check=False)
except subprocess.TimeoutExpired:
    print(
        f"ERROR: source tests exceeded {timeout:g} seconds. "
        "The process was terminated. Re-run the printed command manually, or set "
        "PYTEST_TIMEOUT_SECONDS to a larger value if the machine is unusually slow.",
        file=sys.stderr,
        flush=True,
    )
    raise SystemExit(124)
raise SystemExit(result.returncode)
PY

  printf '\n[3/4] Run real RADIA operating point\n'
  PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" -m apple2.cli \
    --period 40 --periods 10 --gap 12 --phase-mm 10 \
    --samples 1201 --central-periods 6 --harmonics 1 3 5 \
    --csv "$OUT/apple2_field.csv" --json "$OUT/apple2_analysis.json"

  printf '\n[4/4] Run real RADIA bounded phase optimization\n'
  PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" -m apple2.scan \
    --period 40 --periods 10 --gap 12 \
    --phase-min -20 --phase-max 20 --phase-steps 41 \
    --refinement-levels 3 --refinement-steps 11 \
    --samples 801 --central-periods 6 --harmonics 1 3 5 \
    --csv "$OUT/apple2_phase_scan.csv" --json "$OUT/apple2_phase_optimization.json"

  if [[ -n "$REFERENCE_CSV" ]]; then
    printf '\n[optional] Compare against supplied reference CSV\n'
    PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" -m apple2.reference \
      "$OUT/apple2_phase_scan.csv" "$REFERENCE_CSV" \
      --json "$OUT/apple2_reference_comparison.json"
  fi
}

run_logged_validation() {
  local status
  set +e
  (set -e; validation_body)
  status=$?
  if [[ $status -eq 0 ]]; then
    printf '\nValidation completed without command failure.\n'
  else
    printf '\nValidation stopped with exit status %d.\n' "$status"
  fi
  return "$status"
}

# The validation body is the only producer of validation.log. The tee pipeline is
# allowed to close before any digest is calculated, so the final log hash cannot
# be invalidated by later log output.
set +e
run_logged_validation 2>&1 | tee "$LOG"
validation_status=${PIPESTATUS[0]}
set -e

# Exact-byte evidence is deliberately written after validation.log is closed.
# The log hash is held in a separate sidecar, avoiding self-reference. The output
# manifest includes validation.log and its sidecar, then receives its own sidecar.
write_sha256_sidecar "$LOG" "$LOG_SIDECAR"
write_output_manifest

printf '\nFinal exact-byte evidence (written after validation.log closed)\n'
cat "$MANIFEST"
printf '%s\n' "$(cat "$MANIFEST_SIDECAR")"
printf 'Verify from the output directory with:\n'
printf '  shasum -a 256 -c validation.log.sha256\n'
printf '  shasum -a 256 -c SHA256SUMS.txt.sha256\n'

exit "$validation_status"
