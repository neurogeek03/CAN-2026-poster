"""
Extract spatial coordinates + cell class/subclass for the two comparison
samples (OIL + CORT) into a CSV used by the brain subclass figure.

Output: data/processed/brain_subclass_coords.csv
Columns: bead_id, x, y, sample, treatment, cell_class, cell_subclass

Usage:
    uv run python scripts/extract/extract_brain_subclass.py
"""

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

ROOT        = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"
OUTPUT_PATH = ROOT / "data" / "processed" / "brain_subclass_coords.csv"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def derive_class(cell_type: str) -> tuple[str, str]:
    body   = re.sub(r"^\d+_", "", str(cell_type))
    suffix = body.split("_")[-1]
    if suffix in ("Glut", "Gaba", "IMN"):
        cls = suffix
    elif suffix == "NN":
        cls = "Non-Neuron"
    else:
        cls = "Other"
    return cls, body


def main():
    cfg = load_config()
    h5ad_path    = ROOT / cfg["h5ad_path"]
    keep_samples = {s["id"] for s in cfg["samples"]}

    print(f"Reading {h5ad_path.name} …")
    adata = ad.read_h5ad(h5ad_path, backed="r")

    mask = adata.obs["sample"].isin(keep_samples).values
    obs  = adata.obs.loc[mask, ["sample", "treatment", "RCTD_first_type_rat"]].copy()
    xy   = adata.obsm["X_spatial"][mask]

    classes, subclasses = zip(
        *obs["RCTD_first_type_rat"].astype(str).apply(derive_class)
    )

    df = pd.DataFrame({
        "bead_id":      obs.index,
        "x":            xy[:, 0].astype("float32"),
        "y":            xy[:, 1].astype("float32"),
        "sample":       obs["sample"].values,
        "treatment":    obs["treatment"].values,
        "cell_class":   list(classes),
        "cell_subclass": list(subclasses),
    })

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)
    print(f"Saved {len(df):,} beads → {OUTPUT_PATH}")
    print(df.groupby(["sample", "treatment", "cell_class"]).size().to_string())


if __name__ == "__main__":
    main()
