
# Superperiod Optimizer

This script optimizes a higher-order achromat superperiod.

## Features
- Live interactive $\chi^2$ plot during optimization  
- Final dashboard figure (`dashboard.png`) with $\chi^2$ and constraints  
- Logs progress to console

## Requirements
- Python 3.7+
- Packages: `numpy`, `scipy`, `matplotlib`
- `thor_scsi` installed and a valid `.lat` lattice file

## Usage

python3 optimize_superperiod.py my_lattice_file

where `my_lattice_file.lat` is your lattice file (without `.lat` extension).

You can specify output dashboard filename with:

python3 optimize_superperiod.py my_lattice_file --dashboard results.png

## Output
- Interactive live plot of $\chi^2$
- Saved dashboard figure: `dashboard.png`

---

(c) 2025
