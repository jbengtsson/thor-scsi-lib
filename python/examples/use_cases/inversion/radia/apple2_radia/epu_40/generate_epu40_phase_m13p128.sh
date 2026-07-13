#!/usr/bin/env bash
#
# Generate the EPU-40 dataset for phase phase_m13p128.
#
# This is a source-only workflow:
#   PYTHONPATH=<repo root> python3 -m apple2.cli
#
# No pip installation of apple2-radia is required.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ROOT="${APPLE2_ROOT:-$SCRIPT_DIR}"
OUTPUT_SUBDIR="epu_40/runs/phase_m13p128"

PYTHON_BIN="${PYTHON_BIN:-python3}"
FORCE=0
DRY_RUN=0

usage() {
    cat <<'EOF'
Usage:
  generate_epu40_phase_m13p128.sh [options]

Options:
  -o, --output-subdir DIR  Output directory, relative to the source root unless
                           an absolute path is supplied.
                           Default: epu_40/runs/phase_m13p128
      --root DIR           APPLE-II toolkit source root.
                           Default: directory containing this script, or APPLE2_ROOT.
      --python COMMAND     Python interpreter. Default: PYTHON_BIN or python3
  -f, --force              Replace existing generated files in the target.
      --dry-run            Print the resolved command without running RADIA.
  -h, --help               Show this help.

Generated files:
  field.csv
  analysis.json
  generation_command.txt
  BYTES.txt
  SHA256SUMS.txt
  SHA256SUMS.txt.sha256
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output-subdir)
            [[ $# -ge 2 ]] || die "$1 requires a directory"
            OUTPUT_SUBDIR="$2"
            shift 2
            ;;
        --root)
            [[ $# -ge 2 ]] || die "$1 requires a directory"
            ROOT="$2"
            shift 2
            ;;
        --python)
            [[ $# -ge 2 ]] || die "$1 requires a command"
            PYTHON_BIN="$2"
            shift 2
            ;;
        -f|--force)
            FORCE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown option: $1"
            ;;
    esac
done

ROOT="$(cd "$ROOT" 2>/dev/null && pwd -P)" ||
    die "source root does not exist: $ROOT"

[[ -f "$ROOT/pyproject.toml" ]] || die "missing $ROOT/pyproject.toml"
[[ -d "$ROOT/apple2" ]] || die "missing $ROOT/apple2"
[[ -f "$ROOT/apple2/cli.py" ]] || die "missing $ROOT/apple2/cli.py"

command -v "$PYTHON_BIN" >/dev/null 2>&1 ||
    die "Python interpreter not found: $PYTHON_BIN"

if [[ "$OUTPUT_SUBDIR" = /* ]]; then
    OUTPUT_DIR="$OUTPUT_SUBDIR"
else
    OUTPUT_DIR="$ROOT/$OUTPUT_SUBDIR"
fi

FIELD_CSV="$OUTPUT_DIR/field.csv"
ANALYSIS_JSON="$OUTPUT_DIR/analysis.json"
COMMAND_FILE="$OUTPUT_DIR/generation_command.txt"
BYTES_FILE="$OUTPUT_DIR/BYTES.txt"
MANIFEST_FILE="$OUTPUT_DIR/SHA256SUMS.txt"
MANIFEST_HASH_FILE="$OUTPUT_DIR/SHA256SUMS.txt.sha256"

if [[ "$FORCE" -ne 1 ]]; then
    for candidate in \
        "$FIELD_CSV" \
        "$ANALYSIS_JSON" \
        "$COMMAND_FILE" \
        "$BYTES_FILE" \
        "$MANIFEST_FILE" \
        "$MANIFEST_HASH_FILE"
    do
        [[ ! -e "$candidate" ]] ||
            die "target already exists: $candidate (use --force to replace)"
    done
fi

cmd=(
    "$PYTHON_BIN" -m apple2.cli
    --period 40
    --periods 10
    --gap 12
    --phase-mm -13.128
    --samples 1201
    --central-periods 6
    --harmonics 1 3 5
    --min-total-amplitude 1e-4
    --max-relative-residual 0.10
    --csv "$FIELD_CSV"
    --json "$ANALYSIS_JSON"
)

echo "EPU-40 dataset generation"
echo "phase:             -13.128 mm"
echo "source root:       $ROOT"
echo "output directory:  $OUTPUT_DIR"
echo "python command:    $PYTHON_BIN"
echo

printf 'Resolved command:\n  '
printf '%q ' env "PYTHONPATH=$ROOT${PYTHONPATH:+:$PYTHONPATH}" "${cmd[@]}"
printf '\n'

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo
    echo "Dry run only; no files were generated."
    exit 0
fi

mkdir -p "$OUTPUT_DIR"

if [[ "$FORCE" -eq 1 ]]; then
    rm -f \
        "$FIELD_CSV" \
        "$ANALYSIS_JSON" \
        "$COMMAND_FILE" \
        "$BYTES_FILE" \
        "$MANIFEST_FILE" \
        "$MANIFEST_HASH_FILE"
fi

{
    printf 'working_directory=%q\n' "$ROOT"
    printf 'PYTHONPATH=%q\n' "$ROOT${PYTHONPATH:+:$PYTHONPATH}"
    printf 'command='
    printf '%q ' "${cmd[@]}"
    printf '\n'
} > "$COMMAND_FILE"

(
    cd "$ROOT"
    PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" "${cmd[@]}"
)

[[ -s "$FIELD_CSV" ]] || die "field CSV was not created or is empty"
[[ -s "$ANALYSIS_JSON" ]] || die "analysis JSON was not created or is empty"

"$PYTHON_BIN" - "$ANALYSIS_JSON" "-13.128" <<'PY'
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

path = Path(sys.argv[1])
expected_phase = float(sys.argv[2])

data = json.loads(path.read_text(encoding="utf-8"))
parameters = data["parameters"]

expected = {
    "period_mm": 40.0,
    "n_periods": 10,
    "gap_mm": 12.0,
    "phase_mm": expected_phase,
}

errors = []
for name, wanted in expected.items():
    actual = float(parameters[name])
    if not math.isclose(actual, float(wanted), rel_tol=0.0, abs_tol=1e-9):
        errors.append(f"{name}: emitted={actual!r}, expected={wanted!r}")

if errors:
    raise SystemExit("Generated metadata mismatch:\n  " + "\n  ".join(errors))
PY

(
    cd "$OUTPUT_DIR"

    wc -c \
        field.csv \
        analysis.json \
        generation_command.txt \
        > BYTES.txt

    shasum -a 256 \
        field.csv \
        analysis.json \
        generation_command.txt \
        BYTES.txt \
        > SHA256SUMS.txt

    shasum -a 256 SHA256SUMS.txt > SHA256SUMS.txt.sha256

    shasum -a 256 -c SHA256SUMS.txt.sha256
    shasum -a 256 -c SHA256SUMS.txt
)

echo
echo "Dataset generated successfully:"
echo "  $OUTPUT_DIR"
