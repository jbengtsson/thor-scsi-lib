#!/usr/bin/env bash
set -euo pipefail

# BESSY II EMIL UE48 APPLE-II periodic-core prototype.
#
# Default operating point:
#   horizontal-linear mode
#   phase = 0 mm
#   gap   = 15 mm
#
# Published nominal body parameters represented here:
#   - 48 mm magnetic period
#   - 29 full periods
#   - 15 mm nominal minimum magnetic gap
#   - 40 x 40 mm magnet transverse dimensions
#   - 0.8 mm gap/slit between adjacent magnetic rows
#   - Nd2Fe14B permanent magnets
#   - longitudinally magnetized A blocks: Br = 1.28 T
#   - vertically magnetized B blocks:     Br = 1.33 T
#
# Published magnetic benchmark at the nominal minimum gap:
#   horizontal-linear effective field: 0.791 T
#
# Current-toolkit approximation:
#   apple_ii.cli accepts one common remanence. A standard four-block Halbach
#   period contains equal numbers of A and B blocks, so this prototype uses
#   their arithmetic mean:
#
#       (1.28 T + 1.33 T) / 2 = 1.305 T
#
#   Override it with UE48_COMMON_REMANENCE_T when testing sensitivity.
#
# Still not implemented / not established numerically in the public paper:
#   - orientation-dependent remanence in one model
#   - the published symmetric BESSY end-pole configuration
#   - glued-pair machining details
#   - magic-finger corrections, sorting, shimming, and measured block errors
#   - an authoritative machine-readable field map
#
# The internal quality gate checks harmonic-fit quality only. It does not
# establish agreement with the as-built BESSY II UE48.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/runs/ue_48_phase_0_prototype}"

GAP_MM="${UE48_GAP_MM:-15}"
COMMON_REMANENCE_T="${UE48_COMMON_REMANENCE_T:-1.305}"
INTER_ARRAY_GAP_X_MM="${UE48_INTER_ARRAY_GAP_X_MM:-0.8}"

mkdir -p "$OUT_DIR"

PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" python3 -m apple_ii.cli \
  --period 48 \
  --periods 29 \
  --gap "$GAP_MM" \
  --width 40 \
  --height 40 \
  --remanence "$COMMON_REMANENCE_T" \
  --inter-array-gap-x-mm "$INTER_ARRAY_GAP_X_MM" \
  --motion elliptical \
  --phase-mm 0 \
  --end-block-fraction 0 \
  --end-clearance-mm 0 \
  --samples 4001 \
  --central-periods 12 \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --csv "$OUT_DIR/field.csv" \
  --json "$OUT_DIR/analysis.json" \
  --diagnostics-json "$OUT_DIR/diagnostics.json"
