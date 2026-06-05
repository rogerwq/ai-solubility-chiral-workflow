#!/usr/bin/env python3
"""
Node 01 — Compute Solubility

Runs fastsolv predictions on the validated molecule set.
For each molecule, predicts log(S) solubility (mol/L) in the selected
solvent set at the specified temperature.

Output: solubility_predictions.csv — long format with columns:
  name, smiles, solvent, solvent_smiles, logS, logS_stdev
"""

import argparse
import os
import sys

import pandas as pd

# ── Top-50 solvent set (most common pharma/lab solvents from BigSolDB) ──────
TOP50_SOLVENTS = {
    'ethanol':              'CCO',
    'methanol':             'CO',
    'isopropanol':          'CC(C)O',
    'ethyl acetate':        'CCOC(C)=O',
    'acetone':              'CC(C)=O',
    'n-propanol':           'CCCO',
    'water':                'O',
    'n-butanol':            'CCCCO',
    'acetonitrile':         'CC#N',
    'toluene':              'Cc1ccccc1',
    'Dimethylformamide':    'CN(C)C=O',
    '1,4-dioxane':          'C1COCCO1',
    'isobutanol':           'CC(C)CO',
    '2-butanone':           'CCC(C)=O',
    'Tetrahydrofuran':      'C1CCOC1',
    'cyclohexane':          'C1CCCCC1',
    'n-pentanol':           'CCCCCO',
    'methyl acetate':       'COC(C)=O',
    'sec-butanol':          'CCC(C)O',
    'N-Methyl-2-pyrrolidone': 'CN1CCCC1=O',
    'chloroform':           'ClC(Cl)Cl',
    'n-hexane':             'CCCCCC',
    'ethylene glycol':      'OCCO',
    'Dimethyl sulfoxide':   'CS(C)=O',
    'n-octanol':            'CCCCCCCCO',
    'dichloromethane':      'ClCCl',
    'acetic acid':          'CC(=O)O',
    'n-butyl acetate':      'CCCCOC(C)=O',
    'cyclohexanone':        'O=C1CCCCC1',
    'benzene':              'c1ccccc1',
    '2-ethoxyethanol':      'CCOCCO',
    'n-propyl acetate':     'CCCOC(C)=O',
    'isopentanol':          'CC(C)CCO',
    'n-heptane':            'CCCCCCC',
    'n-hexanol':            'CCCCCCO',
    '1,2-dichloroethane':   'ClCCCl',
    'tetrachloromethane':   'ClC(Cl)(Cl)Cl',
    'propylene glycol':     'CC(O)CO',
    'Dimethylacetamide':    'CC(=O)N(C)C',
    'ethyl formate':        'CCOC=O',
    'isopropyl acetate':    'CC(=O)OC(C)C',
    'n-heptanol':           'CCCCCCCO',
    '2-methoxyethanol':     'COCCO',
    'Methyl tert-butyl ether': 'COC(C)(C)C',
    'isobutyl acetate':     'CC(=O)OCC(C)C',
    'p-xylene':             'Cc1ccc(C)cc1',
    'ethylbenzene':         'CCc1ccccc1',
    'diethyl ether':        'CCOCC',
    'o-xylene':             'Cc1ccccc1C',
    'chlorobenzene':        'Clc1ccccc1',
}

# ── All 136 BigSolDB solvents ────────────────────────────────────────────────
ALL_SOLVENTS = {
    **TOP50_SOLVENTS,
    'methanol':                     'CO',
    'n-propanol':                   'CCCO',
    'n-butanol':                    'CCCCO',
    'n-octanol':                    'CCCCCCCCO',
    '1,2-dichlorobenzene':          'Clc1ccccc1Cl',
    '1,2,4-trichlorobenzene':       'Clc1ccc(Cl)c(Cl)c1',
    'Dimethylformamide':            'CN(C)C=O',
    'N-Methyl-2-pyrrolidone':       'CN1CCCC1=O',
    'n-propyl acetate':             'CCCOC(C)=O',
    'diethylene glycol':            'OCCOCCO',
    'diethyl ether':                'CCOCC',
    'tetrachloromethane':           'ClC(Cl)(Cl)Cl',
    'benzene':                      'c1ccccc1',
    'n-butyl acetate':              'CCCCOC(C)=O',
    'benzylalcohol':                'OCc1ccccc1',
    'tert-butanol':                 'CC(C)(C)O',
    'n-heptane':                    'CCCCCCC',
    'n-hexane':                     'CCCCCC',
    'isopentanol':                  'CC(C)CCO',
    'ethyl formate':                'CCOC=O',
    'chlorobenzene':                'Clc1ccccc1',
    'n-heptanol':                   'CCCCCCCO',
    'isopropyl acetate':            'CC(=O)OC(C)C',
    'methyl propionate':            'CCC(=O)OC',
    'Methyl tert-butyl ether':      'COC(C)(C)C',
    'Dimethylacetamide':            'CC(=O)N(C)C',
    'acetic acid':                  'CC(=O)O',
    'formic acid':                  'O=CO',
    'acetylacetone':                'CC(=O)CC(C)=O',
    'propionic acid':               'CCC(=O)O',
    'acrylic acid':                 'C=CC(=O)O',
    'formamide':                    'NC=O',
    'Methyl isobutyl ketone':       'CC(=O)CC(C)C',
    '3-pentanone':                  'CCC(=O)CC',
    'dipropyl ether':               'CCCOCCC',
    'methylcyclohexane':            'CC1CCCCC1',
    'triethyl phosphate':           'CCOP(=O)(OCC)OCC',
    'propylene carbonate':          'CC1COC(=O)O1',
    'cyclopentanone':               'O=C1CCCC1',
    'pyridine':                     'c1ccncc1',
    '2-pentanol':                   'CCCC(C)O',
    'tert-amyl alcohol':            'CCC(C)(C)O',
    '2-methoxyethanol':             'COCCO',
    '2-ethoxyethanol':              'CCOCCO',
    '2-propoxyethanol':             'CCCOCCO',
    'n-hexanol':                    'CCCCCCO',
    '2-butoxyethanol':              'CCCCOCCO',
    'n-pentyl acetate':             'CCCCCOC(C)=O',
    '2-pentanone':                  'CCCC(C)=O',
    '1-methoxy-2-propanol':         'COCC(C)O',
    'ethylbenzene':                 'CCc1ccccc1',
    'p-cymene':                     'Cc1ccc(C(C)C)cc1',
    'p-xylene':                     'Cc1ccc(C)cc1',
    'isooctane':                    'CC(C)CC(C)(C)C',
    'tert-butylamine':              'CC(C)(C)N',
    'n-pentane':                    'CCCCC',
    'isooctanol':                   'CC(C)CCCCCO',
    'morpholine':                   'C1COCCN1',
    'anisole':                      'COc1ccccc1',
    'tributyl phosphate':           'CCCCOP(=O)(OCCCC)OCCCC',
    'diacetone alcohol':            'CC(=O)CC(C)(C)O',
    '2-ethyl-n-hexanol':            'CCCCC(CC)CO',
    '2-methyltetrahydrofuran':      'CC1CCCO1',
    'm-xylene':                     'Cc1cccc(C)c1',
    'benzyl alcohol':               'OCc1ccccc1',
    'isobutyl acetate':             'CC(=O)OCC(C)C',
    'transcutol':                   'CCOCCOCCO',
    '1-propoxy-2-propanol':         'CCCOCC(C)O',
    '1-methoxy-2-propyl acetate':   'COCC(C)OC(C)=O',
    '1,1-dichloroethane':           'CC(Cl)Cl',
    'n-octane':                     'CCCCCCCC',
    'diisopropyl ether':            'CC(C)OC(C)C',
    '2-isopropoxyethanol':          'CC(C)OCCO',
    'n-dodecanol':                  'CCCCCCCCCCCCO',
    'n-decanol':                    'CCCCCCCCCCO',
    'aniline':                      'Nc1ccccc1',
    'acetophenone':                 'CC(=O)c1ccccc1',
    '1,2-diethoxyethane':           'CCOCCOCC',
    'trichloroethylene':            'ClC=C(Cl)Cl',
    '1-bromopropane':               'CCCBr',
    'n-butyric acid':               'CCCC(=O)O',
    '4-methylpyridine':             'Cc1ccncc1',
    'sec-pentanol':                 'CCCC(C)O',
    'sec-butyl acetate':            'CCC(C)OC(C)=O',
    'tert-butyl acetate':           'CC(=O)OC(C)(C)C',
    'methyl formate':               'COC=O',
    'decalin':                      'C1CCC2CCCCC2C1',
    'sulfolane':                    'O=S1(=O)CCCC1',
    '2-(2-butoxyethoxy)ethanol':    'CCCCOCCOCCO',
    '2-hexanone':                   'CCCCC(C)=O',
    'Methyl ethyl ketone':          'CCC(C)=O',
    'nonan-1-ol':                   'CCCCCCCCCO',
    'n-dodecane':                   'CCCCCCCCCCCC',
    'n-hexadecane':                 'CCCCCCCCCCCCCCCC',
    '2-aminoethanol':               'NCCO',
    '1,3-propanediol':              'OCCCO',
    '2-octanol':                    'CCCCCCC(C)O',
    'cumene':                       'CC(C)c1ccccc1',
}


TEMP_MIN_K = -30.0 + 273.15
TEMP_MAX_K = 130.0 + 273.15


def main():
    parser = argparse.ArgumentParser(description="Compute fastsolv solubility predictions.")
    parser.add_argument("--smiles_column", default="smiles")
    parser.add_argument("--solvents", default="top50", choices=["top50", "all"])
    parser.add_argument("--temperature_c", type=float, default=25.0)
    args = parser.parse_args()

    # --- Validate temperature ---
    temp_c = args.temperature_c
    if not (-30.0 <= temp_c <= 130.0):
        print(f"WARNING: temperature_c={temp_c} is outside fastsolv validated range "
              f"(-30 to 130 °C). Clamping.", flush=True)
        temp_c = max(-30.0, min(130.0, temp_c))
    temp_k = temp_c + 273.15

    # --- Load input ---
    input_path = "./inputs/standardized_molecules.csv"
    if not os.path.exists(input_path):
        input_path = "../input_files/sample_molecules.csv"
    if not os.path.exists(input_path):
        print("ERROR: standardized_molecules.csv not found in inputs/", flush=True)
        import sys; sys.exit(1)

    df = pd.read_csv(input_path)
    if args.smiles_column not in df.columns:
        print(f"ERROR: Column '{args.smiles_column}' not found.", flush=True)
        import sys; sys.exit(1)

    smiles_list = df[args.smiles_column].tolist()
    names_list  = df["name"].tolist() if "name" in df.columns else [s[:20] for s in smiles_list]
    print(f"Running fastsolv on {len(smiles_list)} molecules "
          f"({args.solvents} solvent set, {temp_c:.1f} °C) …", flush=True)

    # --- Select solvent set ---
    pool = TOP50_SOLVENTS if args.solvents == "top50" else ALL_SOLVENTS
    solvent_names  = list(pool.keys())
    solvent_smiles = list(pool.values())
    n_solvents = len(solvent_names)
    print(f"Solvent set: {n_solvents} solvents.", flush=True)

    # --- Build fastsolv input DataFrame (one row per molecule–solvent pair) ---
    from fastsolv import fastsolv

    rows = []
    for mol_smiles, mol_name in zip(smiles_list, names_list):
        df_in = pd.DataFrame({
            "solute_smiles":  [mol_smiles]  * n_solvents,
            "solvent_smiles": solvent_smiles,
            "temperature":    [temp_k]      * n_solvents,
        })
        try:
            result = fastsolv(df_in).reset_index(drop=True)
            preds = result["predicted_logS"].tolist()
            stds  = result["predicted_logS_stdev"].tolist()
        except Exception as e:
            print(f"  WARNING: fastsolv failed for {mol_name} ({mol_smiles[:40]}): {e}", flush=True)
            preds = [None] * n_solvents
            stds  = [None] * n_solvents

        for sv_name, sv_smiles, logS, logS_std in zip(solvent_names, solvent_smiles, preds, stds):
            rows.append({
                "name":          mol_name,
                "smiles":        mol_smiles,
                "solvent":       sv_name,
                "solvent_smiles": sv_smiles,
                "temperature_c": temp_c,
                "logS":          round(float(logS), 4) if logS is not None else None,
                "logS_stdev":    round(float(logS_std), 4) if logS_std is not None else None,
            })

        print(f"  Done: {mol_name}", flush=True)

    # --- Save output ---
    out_df = pd.DataFrame(rows)
    os.makedirs("./outputs", exist_ok=True)
    out_df.to_csv("./outputs/solubility_predictions.csv", index=False)
    print(f"Saved solubility_predictions.csv: "
          f"{len(smiles_list)} molecules × {n_solvents} solvents = {len(out_df)} rows.", flush=True)


if __name__ == "__main__":
    main()
