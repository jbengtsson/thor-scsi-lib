#!/usr/bin/env bash
set -e

mkdir -p runs/phase_0

python3 -m apple_ii.cli \
  --period 40 \
  --periods 10 \
  --gap 12 \
  --width 30 \
  --height 12 \
  --remanence 1.20 \
  --motion elliptical \
  --phase-mm 0 \
  --end-block-fraction 0.5 \
  --end-magnetization-fraction 0.5 \
  --end-clearance-mm 0 \
  --samples 1201 \
  --central-periods 6 \
  --harmonics 1 3 5 \
  --min-total-amplitude 1e-4 \
  --max-relative-residual 0.10 \
  --csv runs/phase_0/field.csv \
  --json runs/phase_0/analysis.json
