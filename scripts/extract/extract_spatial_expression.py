"""
Extract spatial expression data for a given gene from the slide-seq h5ad.

Outputs a CSV to data/processed/ with per-bead spatial coordinates,
expression value, quality-filter flag, and cell class / subclass labels.

Usage:
    uv run python scripts/extract/extract_spatial_expression.py <gene>
    uv run python scripts/extract/extract_spatial_expression.py Sgk1
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def derive_class(cell_type: str) -> tuple[str, str]:
    """Return (class_label, subclass_label) from RCTD_first_type_rat string.

    e.g. '017_CA3_Glut' → ('Glut', 'CA3_Glut')
         '319_Astro_TE_NN' → ('Non-Neuron', 'Astro_TE_NN')
    """
    # Strip leading number prefix (e.g. '017_')
    body = re.sub(r"^\d+_", "", cell_type)
    suffix = body.split("_")[-1]
    if suffix in ("Glut", "Gaba", "IMN"):
        cls = suffix
    elif suffix == "NN":
        cls = "Non-Neuron"
    else:
        cls = "Other"
    return cls, body


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gene", help="Gene symbol to extract (e.g. Sgk1)")
    args = parser.parse_args()

    cfg = load_config()
    cfg["gene"] = args.gene

    h5ad_path = ROOT / cfg["h5ad_path"]
    gene_symbol = cfg["gene"]
    samples = [s["id"] for s in cfg["samples"]]
    neuron_thresh = cfg["neuron_score_threshold"]
    nn_thresh = cfg["non_neuron_score_threshold"]

    out_path = ROOT / cfg["processed_csv"].format(gene=gene_symbol)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading {h5ad_path.name} …")
    adata = sc.read_h5ad(h5ad_path)

    # ── Look up Ensembl ID for the requested gene symbol ─────────────────────
    gene_row = adata.var[adata.var["gene_symbol"] == gene_symbol]
    if gene_row.empty:
        sys.exit(f"Gene '{gene_symbol}' not found in var['gene_symbol'].")
    ensembl_id = gene_row.index[0]
    print(f"  {gene_symbol} → {ensembl_id}")

    # ── Subset to target samples ──────────────────────────────────────────────
    adata = adata[adata.obs["sample"].isin(samples)].copy()
    print(f"  Beads in target samples: {adata.n_obs:,}")

    # ── Derive class / subclass from RCTD_first_type_rat ─────────────────────
    classes, subclasses = zip(
        *adata.obs["RCTD_first_type_rat"].astype(str).apply(derive_class)
    )
    adata.obs["cell_class"] = list(classes)
    adata.obs["cell_subclass"] = list(subclasses)

    # ── Quality filter flags ──────────────────────────────────────────────────
    is_neuron = adata.obs["cell_class"].isin(["Glut", "Gaba", "IMN"])
    is_nn = ~is_neuron

    neuron_pass = is_neuron & (adata.obs["RCTD_singlet_score_rat"] > neuron_thresh)
    nn_pass = (
        is_nn
        & (adata.obs["RCTD_spot_class_rat"] == "singlet")
        & (adata.obs["RCTD_singlet_score_rat"] > nn_thresh)
    )
    adata.obs["quality_pass"] = (neuron_pass | nn_pass).values

    print(
        f"  Quality-filtered beads: {adata.obs['quality_pass'].sum():,} "
        f"/ {adata.n_obs:,}"
    )

    # ── Extract gene expression ───────────────────────────────────────────────
    gene_idx = list(adata.var_names).index(ensembl_id)
    expr_mat = adata.X[:, gene_idx]
    if hasattr(expr_mat, "toarray"):
        expr = expr_mat.toarray().flatten()
    else:
        expr = np.asarray(expr_mat).flatten()

    # ── Spatial coordinates ───────────────────────────────────────────────────
    xy = adata.obsm["X_spatial"]

    # ── Assemble CSV ──────────────────────────────────────────────────────────
    df = pd.DataFrame(
        {
            "bead_id": adata.obs_names,
            "x": xy[:, 0],
            "y": xy[:, 1],
            "sample": adata.obs["sample"].values,
            "treatment": adata.obs["treatment"].values,
            "cell_type": adata.obs["RCTD_first_type_rat"].values,
            "cell_class": adata.obs["cell_class"].values,
            "cell_subclass": adata.obs["cell_subclass"].values,
            "is_neuron": is_neuron.values,
            "RCTD_spot_class_rat": adata.obs["RCTD_spot_class_rat"].values,
            "RCTD_singlet_score_rat": adata.obs["RCTD_singlet_score_rat"].values,
            "quality_pass": adata.obs["quality_pass"].values,
            f"{gene_symbol}_expr": expr,
        }
    )

    df.to_csv(out_path, index=False)
    n_expressing = (df[f"{gene_symbol}_expr"] > 0).sum()
    print(f"  Saved {len(df):,} beads → {out_path}")
    print(f"  Beads expressing {gene_symbol}: {n_expressing:,} ({100*n_expressing/len(df):.1f}%)")


if __name__ == "__main__":
    main()
