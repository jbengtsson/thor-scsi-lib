#!/usr/bin/env bash
set -euo pipefail

# MAX III I3 / MAX IV ID-lab EPU69.1 periodic-core prototype.
#
# Default operating point:
#   phase = 0 mm
#   gap   = 18 mm
#
# The 18 mm gap is selected because the 2021 MAX IV pulsed-wire thesis
# publishes directly comparable Hall-probe results at that gap. The historic
# I3 source had a 69.1 mm period, 57 poles, a 1.962 m total length, and a
# 16 mm minimum gap.
#
# Published quantities represented here:
#   - 69.1 mm design period
#   - 57 poles, represented as 28 complete magnetic periods
#   - 18 mm comparison gap
#   - Nd-Fe-B APPLE-II / EPU architecture
#
# Prototype assumptions not established by the public sources:
#   - 30 x 30 mm transverse block dimensions
#   - 1.25 T common remanence
#   - zero horizontal inter-array spacing
#   - four equal 17.275 mm blocks per period
#   - no explicit end termination, sorting, shimming, or measured block errors
#
# Override the uncertain assumptions without editing this file:
#   EPU69_BLOCK_WIDTH_MM=...
#   EPU69_BLOCK_HEIGHT_MM=...
#   EPU69_COMMON_REMANENCE_T=...
#   EPU69_INTER_ARRAY_GAP_X_MM=...
#   EPU69_GAP_MM=...
#
# Public Hall-probe benchmark for the same large-period test EPU:
#   measured period:       69.109 mm
#   effective field Beff:  0.7513 T
#   RMS phase error:       2.986 degrees
#   first integral:       -222.2 G.cm
#   second integral:       -39.440 kG.cm^2
#
# The thesis benchmark table does not explicitly repeat the gap. The device's
# closest reported measurement setting is 18 mm, whose plotted pole field is
# consistent with roughly 0.75 T; use this as a comparison hypothesis rather
# than as an independently proven table-to-gap identity.
#
# The apple_ii quality gate checks harmonic-fit quality only. It does not
# establish agreement with the measured MAX IV field data.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/runs/epu_69_phase_0_prototype}"

GAP_MM="${EPU69_GAP_MM:-18}"
BLOCK_WIDTH_MM="${EPU69_BLOCK_WIDTH_MM:-30}"
BLOCK_HEIGHT_MM="${EPU69_BLOCK_HEIGHT_MM:-30}"
COMMON_REMANENCE_T="${EPU69_COMMON_REMANENCE_T:-1.25}"
INTER_ARRAY_GAP_X_MM="${EPU69_INTER_ARRAY_GAP_X_MM:-0.0}"

mkdir -p "$OUT_DIR"

PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" python3 -m apple_ii.cli \
  --period 69.1 \
  --periods 28 \
  --gap "$GAP_MM" \
  --width "$BLOCK_WIDTH_MM" \
  --height "$BLOCK_HEIGHT_MM" \
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
