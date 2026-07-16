#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${1:-runs/fig3}"
if [[ $# -gt 0 ]]; then
    shift
fi

python3 "${SCRIPT_DIR}/epu_18_fig3.py" "${OUTPUT_DIR}" \
    --model-script "${SCRIPT_DIR}/epu_18_gen.py" \
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
    --target-colors red limegreen red \
    "$@"
