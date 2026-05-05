"""
Extract spatial coordinates for the two comparison samples (OIL + CORT)
into a small CSV used by the interactive hemisphere selector.

Output: data/processed/hemisphere_coords.csv
Columns: bead_id, sample, treatment, x, y

Usage:
    uv run python scripts/extract/extract_hemisphere_coords.py
    uv run python scripts/extract/extract_hemisphere_coords.py --cort B03 --oil B14
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"
OUTPUT_PATH = ROOT / "data" / "processed" / "hemisphere_coords.csv"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cort", default=None)
    parser.add_argument("--oil",  default=None)
    args = parser.parse_args()

    cfg = load_config()
    h5ad_path = ROOT / cfg["h5ad_path"]

    samples_cfg = {s["label"]: s["id"] for s in cfg["samples"]}
    cort_id = args.cort or samples_cfg.get("CORT", "B03")
    oil_id  = args.oil  or samples_cfg.get("OIL",  "B14")

    keep_samples = {cort_id, oil_id}

    print(f"Reading {h5ad_path.name} …")
    adata = ad.read_h5ad(h5ad_path, backed="r")

    mask = adata.obs["sample"].isin(keep_samples).values
    xy   = adata.obsm["X_spatial"][mask]
    obs  = adata.obs.loc[mask, ["sample", "treatment"]].copy()

    df = pd.DataFrame({
        "bead_id":   obs.index,
        "sample":    obs["sample"].values,
        "treatment": obs["treatment"].values,
        "x":         xy[:, 0].astype("float32"),
        "y":         xy[:, 1].astype("float32"),
    })

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)
    print(f"Saved {len(df):,} beads → {OUTPUT_PATH}")
    print(df.groupby(["sample", "treatment"]).size().to_string())


if __name__ == "__main__":
    main()
