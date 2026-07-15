#!/usr/bin/env bash
set -euo pipefail

# AQUA-style APPLE-X, nominal helical/circular mode (FEL2022 WEP38 reference).

OUTPUT_DIR="${1:-runs/hel}"

exec python3 epu_18_gen.py \
    circular+ \
    "${OUTPUT_DIR}" \
    --period-mm 18.0 \
    --periods 110 \
    --gap-mm 1.5 \
    --br-t 1.35 \
    --block-xy-mm 18.0 \
    --central-hole-mm 5.5 \
    --chamfer-mm 1.5 \
    --samples 4001 \
    --fit-periods 20 \
    --energy-gev 1.0
