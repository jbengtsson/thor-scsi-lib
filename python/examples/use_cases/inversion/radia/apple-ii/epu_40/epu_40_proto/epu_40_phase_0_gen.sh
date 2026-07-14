#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/runs/phase_0}"

mkdir -p "$OUT_DIR"

python3 -m apple_ii.cli \
  --period 40 \
  --periods 17 \
  --gap 15 \
  --width 28.5714285714 \
  --height 28.5714285714 \
  --remanence 1.25 \
  --inter-array-gap-x-mm 0.3571428571 \
  --motion elliptical \
  --phase-mm 0 \
  --end-block-fraction 0 \
  --end-block-length-mm 4.9642857143 \
  --end-magnetization-fraction 1.0 \
  --end-clearance-mm 0 \
  --samples 2401 \
  --central-periods 12 \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --energy-ev 430e6 \
  --field-plot "$OUT_DIR/field.png" \
  --csv "$OUT_DIR/field.csv" \
  --json "$OUT_DIR/analysis.json" \
  --diagnostics-json "$OUT_DIR/diagnostics.json" \
  --kickmap-prefix "$OUT_DIR/kickmap"
