#!/usr/bin/env bash
set -euo pipefail

# MAX IV HIPPIE EPU53 periodic-core prototype at minimum gap and zero phase.
# This script resolves the apple_ii source package relative to its own location.
# The RADIA extension must already be importable, for example through the user's
# existing PYTHONPATH entry for Radia/env/radia_python.
#
# Published nominal body parameters represented here:
#   - 53 mm magnetic period
#   - 70 full periods
#   - 11 mm minimum gap
#   - 30 x 30 x 13.25 mm magnet-block envelope
#   - NdFeB remanence:
#       vertically magnetized blocks:     1.28 T
#       longitudinally magnetized blocks: 1.25 T
#
# Current-toolkit approximation:
#   apple_ii.cli presently accepts one common remanence. Because a standard
#   four-block Halbach period contains equal numbers of vertical and
#   longitudinal blocks, this prototype uses their arithmetic mean, 1.265 T.
#   Override it with EPU53_COMMON_REMANENCE_T when testing sensitivity.
#
# Still not implemented / not established from the public design information:
#   - orientation-dependent remanence in a single model
#   - actual termination-block geometry
#   - exact horizontal spacing between adjacent left/right arrays
#   - mounting cut-outs, keepers, sorting, shimming, and measured block errors
#   - authoritative measured field-map comparison
#
# The internal quality gate checks harmonic-fit quality only. It does not
# establish agreement with the as-built MAX IV device.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/runs/epu_53_phase_0_prototype}"

COMMON_REMANENCE_T="${EPU53_COMMON_REMANENCE_T:-1.265}"
INTER_ARRAY_GAP_X_MM="${EPU53_INTER_ARRAY_GAP_X_MM:-0.0}"

mkdir -p "$OUT_DIR"

PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" python3 -m apple_ii.cli \
  --period 53 \
  --periods 70 \
  --gap 11 \
  --width 30 \
  --height 30 \
  --remanence "$COMMON_REMANENCE_T" \
  --inter-array-gap-x-mm "$INTER_ARRAY_GAP_X_MM" \
  --motion elliptical \
  --phase-mm 0 \
  --end-block-fraction 0 \
  --end-clearance-mm 0 \
  --samples 8401 \
  --central-periods 12 \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --csv "$OUT_DIR/field.csv" \
  --json "$OUT_DIR/analysis.json" \
  --diagnostics-json "$OUT_DIR/diagnostics.json"
