"""
Extract UMAP coordinates and metadata from the Slide-tags h5ad file.

Output: data/processed/umap_coords.csv
Columns: cell_id, subclass_name, treatment, UMAP1, UMAP2
"""

import anndata as ad
import pandas as pd
from pathlib import Path

RAW = Path("data/raw/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad")
OUT = Path("data/processed/umap_coords.csv")

print(f"Loading {RAW} ...")
adata = ad.read_h5ad(RAW, backed="r")

# Verify expected columns exist
required_obs = ["subclass_name", "treatment"]
missing = [c for c in required_obs if c not in adata.obs.columns]
if missing:
    print(f"WARNING: missing obs columns: {missing}")
    print(f"Available obs columns:\n{list(adata.obs.columns)}")

# Verify UMAP embedding exists
if "X_umap" not in adata.obsm:
    print(f"WARNING: X_umap not found in obsm.")
    print(f"Available obsm keys: {list(adata.obsm.keys())}")
    raise KeyError("X_umap not found — check obsm keys above and update script.")

# Build output dataframe
df = pd.DataFrame(
    adata.obsm["X_umap"],
    index=adata.obs_names,
    columns=["UMAP1", "UMAP2"],
)
df.index.name = "cell_id"

for col in required_obs:
    if col in adata.obs.columns:
        df[col] = adata.obs[col].values

OUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT)

print(f"Saved {len(df):,} cells to {OUT}")
print(df.head())
print(f"\nsubclass_name value counts:\n{df['subclass_name'].value_counts()}")
print(f"\ntreatment value counts:\n{df['treatment'].value_counts()}")
