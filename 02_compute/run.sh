#!/usr/bin/env bash
set -euo pipefail
python compute_solubility.py \
  --smiles_column  "${smiles_column:-smiles}" \
  --solvents       "${solvents:-top50}" \
  --temperature_c  "${temperature_c:-25.0}"
