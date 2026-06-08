# ai-solubility-chiral-workflow

A [Silva](https://github.com/chiral-data/silva) workflow for predicting organic solubility in multiple solvents, powered by **[fastsolv](https://github.com/JacksonBurns/fastsolv)** (MIT, 2025).

Designed for experimentalists in chemical and pharmaceutical industries who need to shortlist the best solvents for purification without exhaustively screening a long solvent list. For a given set of molecules, the pipeline predicts log(S) solubility (mol/L) across up to 136 solvents at any temperature from −30 to 130 °C, ranks them, and delivers an interactive HTML report with a heatmap and per-molecule bar charts.

**Interactive demo:** https://huggingface.co/spaces/AliGhiami/Solubility

---

## Quick test

A single-molecule input is provided in `test/aspirin_test.csv` to verify the pipeline is working correctly before running your full compound set.

Use it as the workflow input in Silva. The expected output for Aspirin at 25 °C is:

| Best solvent | logS (approx.) | Category |
|---|---|---|
| water | ~ −2.0 | 10–1000 mmol/L |
| ethanol | ~ −1.5 | 10–1000 mmol/L |
| methanol | ~ −1.5 | 10–1000 mmol/L |

> Values are approximate — small run-to-run variation (~±0.01) is normal due to PyTorch non-determinism.

---

## Workflow structure

```
00_resolve_validate  →  01_compute  →  02_visualize
```

| Node | Input | Output | Description |
|------|-------|--------|-------------|
| `00_resolve_validate` | `molecules.csv` | `standardized_molecules.csv`, `validation_report.json` | Resolves missing SMILES via PubChem, canonicalizes with RDKit, deduplicates |
| `01_compute` | `standardized_molecules.csv` | `solubility_predictions.csv` | Runs fastsolv ensemble (4 checkpoints) in long format: one row per molecule–solvent pair |
| `02_visualize` | `solubility_predictions.csv` | `solubility_report.html` | Interactive heatmap + per-molecule ranked bar charts |

---

## Docker image

Silva auto-builds the image from `apps/solubility_pipeline_2026_06_05/Dockerfile`.

The Dockerfile installs `fastsolv`, `rdkit`, `pandas`, `numpy`, and `requests`, then pre-downloads the 4 ensemble checkpoints from Zenodo (~900 MB total) so predictions run fully offline.

To build manually:

```bash
docker build -t solubility_pipeline:2026_06_05 apps/solubility_pipeline_2026_06_05/
```

---

## Parameters (`global_params.json`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `smiles_column` | `"smiles"` | Column name for SMILES in the input CSV |
| `name_column` | `"name"` | Column name for molecule names |
| `solvents` | `"top50"` | Solvent set: `"top50"` (50 common pharma/lab solvents) or `"all"` (136 BigSolDB solvents) |
| `temperature_c` | `25.0` | Prediction temperature in °C (valid range: −30 to 130) |
| `top_n_solvents` | `20` | Number of solvents shown in the heatmap overview (node 02 only) |

---

## Input format

A CSV with at least one of `name` or `smiles` columns:

```csv
name,smiles
Aspirin,CC(=O)Oc1ccccc1C(=O)O
Ibuprofen,CC(C)Cc1ccc(C(C)C(=O)O)cc1
```

SMILES can be omitted — the pipeline resolves them from PubChem using the name.

---

## Output: `solubility_predictions.csv`

Long format — one row per molecule–solvent pair:

| Column | Description |
|--------|-------------|
| `name` | Molecule name |
| `smiles` | Canonical SMILES |
| `solvent` | Solvent name |
| `solvent_smiles` | Solvent SMILES |
| `temperature_c` | Prediction temperature (°C) |
| `logS` | Predicted log(S) in mol/L |
| `logS_stdev` | Ensemble standard deviation (uncertainty) |

---

## Model

**fastsolv** uses a graph neural network trained on BigSolDB v1 (54,273 solubility measurements, 839 solutes, 136 solvents). An ensemble of 4 checkpoints provides uncertainty estimates via the `logS_stdev` column.

Reference: Krasnov et al. 2022, doi:10.5281/zenodo.6984601  
Checkpoints: https://zenodo.org/records/13943074
