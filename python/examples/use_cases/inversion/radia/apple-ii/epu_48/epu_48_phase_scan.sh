#!/usr/bin/env bash
set -euo pipefail

# BESSY II EMIL UE48 APPLE-II phase scan.
#
# Scans one half magnetic period:
#   phase = 0 ... 24 mm for a 48 mm period
#
# Nominal periodic-core parameters:
#   - 48 mm period
#   - 29 full periods
#   - 15 mm nominal minimum gap
#   - 40 x 40 mm magnet cross-section
#   - 0.8 mm slit between adjacent magnetic rows
#
# Current-toolkit approximation:
#   The published A- and B-block remanences are 1.28 T and 1.33 T.
#   apple_ii.scan currently accepts one common value, so the default is their
#   equal-count mean, 1.305 T.
#
# Environment overrides:
#   UE48_GAP_MM
#   UE48_COMMON_REMANENCE_T
#   UE48_INTER_ARRAY_GAP_X_MM
#
# The scan's numerical quality criteria assess the fitted field only.
# They do not establish agreement with the as-built UE48 or its end fields.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
RUN_DIR="${1:-$SCRIPT_DIR/runs/ue_48_phase_scan_prototype}"

GAP_MM="${UE48_GAP_MM:-15}"
COMMON_REMANENCE_T="${UE48_COMMON_REMANENCE_T:-1.305}"
INTER_ARRAY_GAP_X_MM="${UE48_INTER_ARRAY_GAP_X_MM:-0.8}"

mkdir -p "$RUN_DIR"

PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" python3 -m apple_ii.scan \
  --period 48 \
  --periods 29 \
  --gap "$GAP_MM" \
  --width 40 \
  --height 40 \
  --remanence "$COMMON_REMANENCE_T" \
  --inter-array-gap-x-mm "$INTER_ARRAY_GAP_X_MM" \
  --phase-min 0 \
  --phase-max 24 \
  --phase-steps 25 \
  --refinement-levels 2 \
  --refinement-steps 11 \
  --samples 4001 \
  --central-periods 12 \
  --csv "$RUN_DIR/phase_scan.csv" \
  --json "$RUN_DIR/phase_scan.json"
