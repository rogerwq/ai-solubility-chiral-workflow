#!/usr/bin/env bash
set -euo pipefail
python resolve_validate.py \
  --smiles_column "${smiles_column:-smiles}" \
  --name_column   "${name_column:-name}"
