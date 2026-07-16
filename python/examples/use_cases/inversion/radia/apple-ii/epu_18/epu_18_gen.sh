#!/usr/bin/env bash
set -euo pipefail

# Nominal WEP38 AQUA APPLE-X device in positive-helicity circular mode.
# RADIA FldFocKickPer computes both the periodic focusing-potential matrix and
# the second-order kick matrices on one regular Cartesian mesh. The Python
# script writes a Figure-7-style focusing-potential plot masked to the valid
# diamond aperture, a full-mesh diagnostic potential plot, and a separate central
# fit-domain kick plot. By default the x/y bounds are the axis vertices of the
# rotated-square magnetic aperture derived from the model geometry.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${1:-runs/hel}"
if (( $# > 0 )); then
    shift
fi

exec python3 "${SCRIPT_DIR}/epu_18_gen.py" \
    circular+ \
    "${OUTPUT_DIR}" \
    --period-mm 18.0 \
    --periods 110 \
    --gap-mm 1.5 \
    --br-t 1.35 \
    --block-xy-mm 18.0 \
    --central-aperture-across-flats-mm 5.5 \
    --outer-chamfer-mm 5.0 \
    --samples 4001 \
    --fit-periods 20 \
    --energy-gev 1.0 \
    --susceptibility \
    --kickmap \
    --field-x-points 101 \
    --field-y-points 101 \
    --kickmap-max-harmonic 5 \
    --kickmap-points-per-period 32 \
    --kickmap-derivative-step-mm 0.0 \
    --kick-fit-r-max-mm 0.8 \
    "$@"
