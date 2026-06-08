#!/usr/bin/env python3
"""
Node 00 — Resolve & Validate

Accepts a CSV with a 'name' column, a 'smiles' column, or both.
- If SMILES is missing for a row, resolves it from PubChem using the name.
- Validates and canonicalizes all SMILES strings using RDKit.
- Removes duplicates by canonical SMILES (keeps first occurrence).
- Outputs standardized_molecules.csv and validation_report.json.
"""

import argparse
import json
import os
import sys

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not available — falling back to basic string validation.", flush=True)


def _resolve_smiles_pubchem(name: str) -> str | None:
    """Resolve a molecule name to SMILES via PubChem REST API."""
    try:
        import requests
        from urllib.parse import quote
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
               f"{quote(name)}/property/CanonicalSMILES/JSON")
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            return r.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except Exception:
        pass
    return None


def standardize_smiles(smiles: str) -> str | None:
    """Canonicalize and standardize a SMILES string using RDKit."""
    if not RDKIT_AVAILABLE:
        return smiles.strip() if smiles else None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        mol = rdMolStandardize.Cleanup(mol)
    except Exception:
        pass
    return Chem.MolToSmiles(mol)


def main():
    parser = argparse.ArgumentParser(description="Resolve and validate molecule inputs.")
    parser.add_argument("--smiles_column", default="smiles")
    parser.add_argument("--name_column", default="name")
    args = parser.parse_args()

    # --- Load input ---
    input_path = "./inputs/molecules.csv"
    if not os.path.exists(input_path):
        input_path = "../input_files/sample_molecules.csv"
    if not os.path.exists(input_path):
        print("ERROR: Input file not found at ./inputs/molecules.csv", flush=True)
        sys.exit(1)

    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows from {input_path}.", flush=True)

    smiles_col = args.smiles_column if args.smiles_column in df.columns else None
    name_col = args.name_column if args.name_column in df.columns else None

    if smiles_col is None and name_col is None:
        print(
            f"ERROR: Neither '{args.smiles_column}' nor '{args.name_column}' column found. "
            f"Available columns: {list(df.columns)}",
            flush=True,
        )
        sys.exit(1)

    # --- Resolve & validate ---
    results = []
    invalid_count = 0
    resolved_count = 0

    for _, row in df.iterrows():
        smiles = str(row[smiles_col]).strip() if smiles_col and pd.notna(row.get(smiles_col)) else ""
        name = str(row[name_col]).strip() if name_col and pd.notna(row.get(name_col)) else ""

        # Attempt PubChem resolution when SMILES is absent
        if not smiles and name:
            resolved = _resolve_smiles_pubchem(name)
            if resolved:
                smiles = resolved
                resolved_count += 1
                print(f"  Resolved '{name}' via PubChem.", flush=True)
            else:
                print(f"  Could not resolve '{name}' — skipping.", flush=True)

        if not smiles:
            invalid_count += 1
            continue

        canonical = standardize_smiles(smiles)
        if canonical is None:
            invalid_count += 1
            print(f"  Invalid SMILES skipped: {smiles[:60]}", flush=True)
            continue

        results.append({"name": name, "smiles": canonical})

    # --- Deduplicate by canonical SMILES ---
    out_df = pd.DataFrame(results)
    if out_df.empty:
        print("ERROR: No valid molecules after validation.", flush=True)
        sys.exit(1)

    before_dedup = len(out_df)
    out_df = out_df.drop_duplicates(subset=["smiles"], keep="first").reset_index(drop=True)
    duplicates_removed = before_dedup - len(out_df)

    # --- Save outputs ---
    os.makedirs("./outputs", exist_ok=True)
    out_df.to_csv("./outputs/standardized_molecules.csv", index=False)

    report = {
        "input_count": len(df),
        "invalid_removed": invalid_count,
        "duplicates_removed": duplicates_removed,
        "pubchem_resolved": resolved_count,
        "output_count": len(out_df),
    }
    with open("./outputs/validation_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print(
        f"Validation complete: {len(out_df)} valid molecules "
        f"({invalid_count} invalid, {duplicates_removed} duplicates removed, "
        f"{resolved_count} resolved via PubChem).",
        flush=True,
    )


if __name__ == "__main__":
    main()
