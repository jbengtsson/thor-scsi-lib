#!/usr/bin/env bash
set -e

# Refined Daresbury HU56 prototype at minimum gap and zero phase.
# Run from:
#   apple2_radia/epu_56_daresbury
#
# Implemented refinements:
#   - 0.5 mm horizontal gap between adjacent left/right arrays
#   - explicit 6.95 mm end-block length
#   - full-remanence end-block prototype assumption
#   - geometry and per-row field diagnostics
#
# Still not implemented:
#   - published S1/S2/S3 end spacings
#   - 5 x 5 mm mounting notches
#   - measured/sorted individual magnet-block data
#
# Optional first argument changes the output directory.

OUT_DIR="${1:-runs/phase_0_refined_prototype}"

mkdir -p "$OUT_DIR"

python3 -m apple_ii.cli \
  --period 56 \
  --periods 17 \
  --gap 21 \
  --width 40 \
  --height 40 \
  --remanence 1.25 \
  --inter-array-gap-x-mm 0.5 \
  --motion elliptical \
  --phase-mm 0 \
  --end-block-fraction 0 \
  --end-block-length-mm 6.95 \
  --end-magnetization-fraction 1.0 \
  --end-clearance-mm 0 \
  --samples 2401 \
  --central-periods 12 \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --csv "$OUT_DIR/field.csv" \
  --json "$OUT_DIR/analysis.json" \
  --diagnostics-json "$OUT_DIR/diagnostics.json"
