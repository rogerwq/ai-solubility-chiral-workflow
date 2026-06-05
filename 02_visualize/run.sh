#!/usr/bin/env bash
set -euo pipefail
python visualize.py \
  --top_n_solvents "${top_n_solvents:-20}"
