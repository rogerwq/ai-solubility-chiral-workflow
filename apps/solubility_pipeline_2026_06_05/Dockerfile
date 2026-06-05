FROM continuumio/miniconda3

# fastsolv requires Python >=3.9; pin to 3.12 for consistency
RUN conda install -y python=3.12 setuptools && conda clean -afy

ENV USER=silva \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

# Install Python dependencies
RUN pip install --no-cache-dir fastsolv rdkit pandas numpy requests

# Pre-download the 4 fastsolv ensemble checkpoints from Zenodo (~224 MB each)
# Baking them into the image makes predictions fully offline / reproducible.
RUN python - <<'EOF'
import importlib.util, os, requests
from pathlib import Path

spec = importlib.util.find_spec("fastsolv")
pkg_dir = Path(spec.submodule_search_locations[0])
ckpt_dir = pkg_dir / "checkpoints"
ckpt_dir.mkdir(exist_ok=True)

base_url = "https://zenodo.org/records/13943074/files"
for i in range(1, 5):
    dest = ckpt_dir / f"fastsolv_1701_{i}.ckpt"
    if dest.exists() and dest.stat().st_size >= 100_000_000:
        print(f"Checkpoint {i}/4 already present.", flush=True)
        continue
    url = f"{base_url}/fastsolv_1701_{i}.ckpt"
    print(f"Downloading checkpoint {i}/4 …", flush=True)
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                f.write(chunk)
    print(f"Checkpoint {i}/4 done ({dest.stat().st_size // 1_000_000} MB).", flush=True)
EOF

WORKDIR /workspace
