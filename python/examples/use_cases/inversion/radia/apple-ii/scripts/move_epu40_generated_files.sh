#!/usr/bin/env bash
#
# Move generated outputs for the 40 mm APPLE-II/EPU model into:
#   <apple2_radia source root>/epu_40
#
# Safety:
#   - Default mode is a preview; nothing is moved.
#   - Use --apply to perform the moves.
#   - Source files under apple2/, tests/, scripts/, and documentation are untouched.
#   - A generated directory/file is moved only when its JSON metadata reports
#     parameters.period_mm == 40.0.
#
# Typical use:
#   cd /path/to/apple2_radia
#   ./move_epu40_generated_files.sh
#   ./move_epu40_generated_files.sh --apply
#
set -euo pipefail

MODE="preview"
ROOT="${APPLE2_ROOT:-$PWD}"
TARGET=""
PYTHON_BIN="${PYTHON_BIN:-python3}"

usage() {
    cat <<'EOF'
Usage:
  move_epu40_generated_files.sh [--apply] [--root DIR] [--target DIR]

Options:
  --apply         Perform the moves. Without this option, only preview actions.
  --root DIR      APPLE-II toolkit source root. Default: current directory.
  --target DIR    Destination. Default: <root>/epu_40
  -h, --help      Show this help.

Environment:
  APPLE2_ROOT     Alternative default source root.
  PYTHON_BIN      Python command used to inspect JSON metadata. Default: python3.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)
            MODE="apply"
            shift
            ;;
        --root)
            [[ $# -ge 2 ]] || { echo "ERROR: --root requires a directory." >&2; exit 2; }
            ROOT="$2"
            shift 2
            ;;
        --target)
            [[ $# -ge 2 ]] || { echo "ERROR: --target requires a directory." >&2; exit 2; }
            TARGET="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

ROOT=$(cd "$ROOT" 2>/dev/null && pwd -P) || {
    echo "ERROR: source root does not exist: $ROOT" >&2
    exit 2
}

if [[ ! -f "$ROOT/pyproject.toml" || ! -d "$ROOT/apple2" ]]; then
    cat >&2 <<EOF
ERROR: this does not look like the apple2_radia source root:

  $ROOT

Expected:
  $ROOT/pyproject.toml
  $ROOT/apple2/
EOF
    exit 2
fi

if [[ -z "$TARGET" ]]; then
    TARGET="$ROOT/epu_40"
elif [[ "$TARGET" != /* ]]; then
    TARGET="$ROOT/$TARGET"
fi

# Resolve the parent even when TARGET does not yet exist.
TARGET_PARENT=$(cd "$(dirname "$TARGET")" 2>/dev/null && pwd -P) || {
    echo "ERROR: target parent does not exist: $(dirname "$TARGET")" >&2
    exit 2
}
TARGET="$TARGET_PARENT/$(basename "$TARGET")"

if [[ "$TARGET" == "$ROOT" ]]; then
    echo "ERROR: target must not equal the source root." >&2
    exit 2
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "ERROR: Python command not found: $PYTHON_BIN" >&2
    exit 2
fi

period_is_40_mm() {
    local json_file=$1

    [[ -f "$json_file" ]] || return 1

    "$PYTHON_BIN" - "$json_file" <<'PY'
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

path = Path(sys.argv[1])
try:
    data = json.loads(path.read_text(encoding="utf-8"))
    value = float(data["parameters"]["period_mm"])
except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
    raise SystemExit(1)

raise SystemExit(0 if math.isclose(value, 40.0, rel_tol=0.0, abs_tol=1e-9) else 1)
PY
}

moved_count=0
skipped_count=0

move_item() {
    local source=$1
    local destination_dir=$2
    local destination="$destination_dir/$(basename "$source")"

    [[ -e "$source" ]] || return 0

    if [[ -e "$destination" ]]; then
        echo "SKIP: destination already exists:"
        echo "  $destination"
        skipped_count=$((skipped_count + 1))
        return 0
    fi

    echo "MOVE:"
    echo "  $source"
    echo "  -> $destination"

    if [[ "$MODE" == "apply" ]]; then
        mkdir -p "$destination_dir"
        mv "$source" "$destination"
    fi

    moved_count=$((moved_count + 1))
}

echo "APPLE-II generated-output organizer"
echo "mode:        $MODE"
echo "source root: $ROOT"
echo "target:      $TARGET"
echo "period gate: parameters.period_mm == 40.0"
echo

# 1. Move matching interactive run directories, preserving runs/<name>.
if [[ -d "$ROOT/runs" ]]; then
    for run_dir in "$ROOT"/runs/*; do
        [[ -d "$run_dir" ]] || continue
        [[ "$run_dir" == "$TARGET"* ]] && continue

        metadata="$run_dir/analysis.json"
        if period_is_40_mm "$metadata"; then
            move_item "$run_dir" "$TARGET/runs"
        else
            echo "SKIP: no verified 40 mm metadata: $run_dir"
            skipped_count=$((skipped_count + 1))
        fi
    done
fi

# 2. Move matching macOS validation directories as complete evidence bundles.
for validation_dir in "$ROOT"/validation_macos_*; do
    [[ -d "$validation_dir" ]] || continue

    metadata="$validation_dir/apple2_analysis.json"
    if period_is_40_mm "$metadata"; then
        move_item "$validation_dir" "$TARGET/validation"
    else
        echo "SKIP: no verified 40 mm metadata: $validation_dir"
        skipped_count=$((skipped_count + 1))
    fi
done

# 3. Move root-level single-operating-point outputs as a pair.
if period_is_40_mm "$ROOT/apple2_analysis.json"; then
    move_item "$ROOT/apple2_analysis.json" "$TARGET/root_outputs"
    move_item "$ROOT/apple2_field.csv" "$TARGET/root_outputs"
fi

# 4. Move root-level phase-optimization outputs as a pair.
if period_is_40_mm "$ROOT/apple2_phase_optimization.json"; then
    move_item "$ROOT/apple2_phase_optimization.json" "$TARGET/root_outputs"
    move_item "$ROOT/apple2_phase_scan.csv" "$TARGET/root_outputs"
fi

# Remove the now-empty runs directory only after a real move.
if [[ "$MODE" == "apply" && -d "$ROOT/runs" ]]; then
    rmdir "$ROOT/runs" 2>/dev/null || true
fi

echo
if [[ "$MODE" == "preview" ]]; then
    echo "Preview complete: $moved_count item(s) would be moved; $skipped_count skipped."
    echo "Run again with --apply to perform the moves."
else
    echo "Move complete: $moved_count item(s) moved; $skipped_count skipped."
    echo "Destination: $TARGET"
fi
