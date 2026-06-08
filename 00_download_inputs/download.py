#!/usr/bin/env python3
"""
Node 00 — Download Sample Input

Downloads the bundled sample molecules CSV from the workflow repository so the
pipeline can run end-to-end without a manually supplied input file.

To screen your own compounds: skip this node and place molecules.csv directly
in the workflow's input_files/ directory before running node 01.
"""

import os
import urllib.request

URL = (
    "https://raw.githubusercontent.com/AliGhiami/ai-solubility-chiral-workflow"
    "/main/input_files/sample_molecules.csv"
)
DST = "outputs/molecules.csv"

os.makedirs("outputs", exist_ok=True)
print(f"Downloading sample molecules from:\n  {URL}", flush=True)
urllib.request.urlretrieve(URL, DST)
print(f"Saved {DST}", flush=True)
