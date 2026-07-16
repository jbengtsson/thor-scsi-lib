#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR="${1:-runs/fig3}"

python3 epu_18_fig3.py \
    "${OUTPUT_DIR}" \
    --model-script epu_18_gen.py \
    --mode circular+ \
    --gap-min-mm 1.5 \
    --gap-max-mm 5.0 \
    --gap-points 15 \
    --period-mm 18 \
    --periods 110 \
    --energy-gev 1.0 \
    --fit-periods 20 \
    --samples-per-period 64 \
    --target-wavelength-nm 5.8 4.0 3.5 \
    --target-colors red limegreen red
