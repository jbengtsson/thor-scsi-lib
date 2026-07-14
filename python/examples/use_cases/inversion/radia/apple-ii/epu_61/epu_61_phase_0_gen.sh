#!/usr/bin/env bash
set -euo pipefail

# Prototype MAX IV / MAX-lab EPU61 APPLE-II model at minimum mechanical gap
# and zero longitudinal phase. This script mirrors epu_56_phase_0_gen.sh and
# resolves the apple_ii source package relative to its own location.
#
# The RADIA extension must already be importable, for example through the
# user's existing PYTHONPATH entry for Radia/env/radia_python.
#
# Published nominal inputs represented here:
#   - 61 mm magnetic period
#   - 41.5 periods / 83 full-size poles
#   - 14 mm minimum mechanical gap
#   - 30 x 30 mm block cross-section
#   - 15.20 mm individual block thickness
#   - NdFeB remanence of at least 1.25 T
#
# Toolkit mapping for the published 41.5-period length:
#   - apple_ii currently accepts an integer number of full periods
#   - 41 full periods create 164 quarter-period body blocks per row
#   - two explicit 15.20 mm, full-remanence continuation blocks produce
#     166 modeled blocks per row, corresponding nominally to 41.5 periods
#
# Prototype assumptions and limitations:
#   - horizontal inter-array spacing is not established from the public source;
#     this script uses the toolkit's zero-spacing default explicitly
#   - body-block length remains period/4 = 15.25 mm in the current toolkit;
#     repeated glue layers and the exact 30.40 mm glued-pair construction are
#     therefore not represented independently
#   - the two generic continuation blocks do not reproduce the real EPU61 end
#     sections, whose complete block sequence and dimensions are unpublished
#   - 5 x 5 mm corner cut-outs, measured block errors, magnet sorting,
#     alignment offsets, shimming, and measured field maps are not implemented
#
# This is a diagnostic prototype, not an as-built engineering model.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/runs/phase_0_prototype}"

mkdir -p "$OUT_DIR"

PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" \
python3 -m apple_ii.cli \
  --period 61 \
  --periods 41 \
  --gap 14 \
  --width 30 \
  --height 30 \
  --remanence 1.25 \
  --inter-array-gap-x-mm 0.0 \
  --motion elliptical \
  --phase-mm 0 \
  --end-block-fraction 0 \
  --end-block-length-mm 15.20 \
  --end-magnetization-fraction 1.0 \
  --end-clearance-mm 0 \
  --samples 4801 \
  --central-periods 32 \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --csv "$OUT_DIR/field.csv" \
  --json "$OUT_DIR/analysis.json" \
  --diagnostics-json "$OUT_DIR/diagnostics.json"
