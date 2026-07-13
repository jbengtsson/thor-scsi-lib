#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$SCRIPT_DIR/runs/phase_0_refined_prototype"

mkdir -p "$RUN_DIR"

export PYTHONPATH="$SCRIPT_DIR/../python${PYTHONPATH:+:$PYTHONPATH}"

python3 -m apple_ii.scan \
  --period 56 \
  --periods 17 \
  --gap 21 \
  --width 40 \
  --height 40 \
  --remanence 1.25 \
  --inter-array-gap-x-mm 0.5 \
  --end-block-length-mm 6.95 \
  --end-magnetization-fraction 1.0 \
  --phase-min 0 \
  --phase-max 28 \
  --phase-steps 29 \
  --refinement-levels 2 \
  --refinement-steps 11 \
  --samples 2401 \
  --central-periods 12 \
  --csv "$RUN_DIR/phase_scan.csv" \
  --json "$RUN_DIR/phase_scan.json"
