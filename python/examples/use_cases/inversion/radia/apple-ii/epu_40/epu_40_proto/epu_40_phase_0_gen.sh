#!/usr/bin/env bash
set -euo pipefail

# Uniformly scaled EPU-40 prototype derived from the attached EPU-56 scripts.
#
# Geometric scale factor: 40 / 56 = 5 / 7.
# The following EPU-56 lengths are scaled by 5/7:
#   period:                   56.00 -> 40.000000 mm
#   minimum magnetic gap:     21.00 -> 15.000000 mm
#   block width and height:   40.00 -> 28.571429 mm
#   inter-array x clearance:   0.50 ->  0.357143 mm
#   end-block length:          6.95 ->  4.964286 mm
#   circular/helical phase:   18.16 -> 12.971429 mm
#
# Dimensionless and non-geometric inputs are retained from the EPU-56
# prototype: 17 periods, Br=1.25 T, 430 MeV, sampling and fit settings.
# Set PERIODS=24 when invoking this script to obtain approximately the same
# active magnetic length as the 17-period EPU-56 (960 mm versus 952 mm).
#
# Usage:
#   ./epu_40_scaled_gen.sh [helical|phase0] [output-directory]
#
# Examples:
#   ./epu_40_scaled_gen.sh helical
#   ./epu_40_scaled_gen.sh phase0 ./runs/epu40_phase0
#   PERIODS=24 ./epu_40_scaled_gen.sh helical ./runs/epu40_24p_helical
#
# APPLE_II_ROOT may be set explicitly. Otherwise this script looks for the
# repository either in its own directory or one directory above it.
# RADIA must already be importable by the selected Python interpreter.
#
# This preserves the prototype's known simplifications. It does not add the
# Daresbury S1/S2/S3 end spacings, mounting notches, measured block errors,
# block sorting, or a physically optimized end-field correction.

MODE="${1:-helical}"
USER_OUT_DIR="${2:-}"

case "$MODE" in
  helical)
    PHASE_MM="12.9714285714"
    DEFAULT_RUN_NAME="epu40_hel_phase_12p9714"
    ;;
  phase0|phase-0|zero)
    PHASE_MM="0"
    DEFAULT_RUN_NAME="epu40_phase_0"
    ;;
  *)
    printf 'Error: mode must be "helical" or "phase0"; received: %s\n' "$MODE" >&2
    exit 2
    ;;
esac

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "${APPLE_II_ROOT:-}" ]]; then
  ROOT="$(cd -- "$APPLE_II_ROOT" && pwd)"
elif [[ -d "$SCRIPT_DIR/python/apple_ii" ]]; then
  ROOT="$SCRIPT_DIR"
elif [[ -d "$SCRIPT_DIR/../python/apple_ii" ]]; then
  ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
else
  cat >&2 <<'MSG'
Error: could not locate python/apple_ii.
Place this script in the repository root or one directory below it, or set:
  APPLE_II_ROOT=/absolute/path/to/repository
MSG
  exit 2
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
PERIODS="${PERIODS:-17}"
CENTRAL_PERIODS="${CENTRAL_PERIODS:-12}"
SAMPLES="${SAMPLES:-2401}"
OUT_DIR="${USER_OUT_DIR:-$SCRIPT_DIR/runs/$DEFAULT_RUN_NAME}"

mkdir -p "$OUT_DIR"

printf 'EPU-40 uniformly scaled prototype\n'
printf '  mode:             %s\n' "$MODE"
printf '  periods:          %s\n' "$PERIODS"
printf '  phase:            %s mm\n' "$PHASE_MM"
printf '  output directory: %s\n' "$OUT_DIR"

PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" -m apple_ii.cli \
  --period 40 \
  --periods "$PERIODS" \
  --gap 15 \
  --width 28.5714285714 \
  --height 28.5714285714 \
  --remanence 1.25 \
  --inter-array-gap-x-mm 0.3571428571 \
  --motion elliptical \
  --phase-mm "$PHASE_MM" \
  --end-block-fraction 0 \
  --end-block-length-mm 4.9642857143 \
  --end-magnetization-fraction 1.0 \
  --end-clearance-mm 0 \
  --samples "$SAMPLES" \
  --central-periods "$CENTRAL_PERIODS" \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --energy-ev 430e6 \
  --field-plot "$OUT_DIR/field.png" \
  --csv "$OUT_DIR/field.csv" \
  --json "$OUT_DIR/analysis.json" \
  --diagnostics-json "$OUT_DIR/diagnostics.json" \
  --kickmap-prefix "$OUT_DIR/kickmap"
