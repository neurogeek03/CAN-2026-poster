"""
Figure 01 — UMAP colored by subclass_name (Slide-tags dataset)

Input:  data/processed/umap_coords.csv
        data/cluster_annotation_term.csv
Output: figures/fig01_umap.svg
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
UMAP_CSV    = Path("data/processed/umap_coords.csv")
COLOR_CSV   = Path("data/cluster_annotation_term.csv")
OUT_SVG     = Path("figures/fig01_umap.svg")

# ── Config ─────────────────────────────────────────────────────────────────────
COLOR_COL   = "subclass_name"
POINT_SIZE  = 0.8       # adjust for density
POINT_ALPHA = 0.7
FIG_W, FIG_H = 10, 8   # inches — tune to poster panel dimensions
LEGEND_NCOL = 2
LEGEND_FONTSIZE = 6

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(UMAP_CSV)

# ── Build color palette ────────────────────────────────────────────────────────
color_df = pd.read_csv(COLOR_CSV, usecols=["name", "cluster_annotation_term_set_name", "color_hex_triplet"])
subclass_colors = (
    color_df[color_df["cluster_annotation_term_set_name"] == "subclass"]
    .set_index("name")["color_hex_triplet"]
    .to_dict()
)

# ── Determine plot order: sort by numeric prefix so legend is ordered ──────────
def _num_prefix(label):
    try:
        return int(str(label).split(" ")[0])
    except ValueError:
        return 9999

present = df[COLOR_COL].dropna().unique().tolist()
present_sorted = sorted(present, key=_num_prefix)
df[COLOR_COL] = pd.Categorical(df[COLOR_COL], categories=present_sorted, ordered=True)
df = df.sort_values(COLOR_COL)

# ── Plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

for subclass in present_sorted:
    mask = df[COLOR_COL] == subclass
    color = subclass_colors.get(subclass, "#808080")
    ax.scatter(
        df.loc[mask, "UMAP1"],
        df.loc[mask, "UMAP2"],
        c=color,
        s=POINT_SIZE,
        alpha=POINT_ALPHA,
        linewidths=0,
        rasterized=True,   # rasterize dots inside the SVG for file size
    )

ax.set_xlabel("UMAP1", fontsize=10)
ax.set_ylabel("UMAP2", fontsize=10)
ax.set_title(f"Slide-tags — {COLOR_COL} (n = {len(df):,})", fontsize=11)
ax.axis("equal")
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# ── Legend ─────────────────────────────────────────────────────────────────────
legend_elements = [
    Patch(facecolor=subclass_colors.get(s, "#808080"), label=s)
    for s in present_sorted
]
ax.legend(
    handles=legend_elements,
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    ncol=LEGEND_NCOL,
    fontsize=LEGEND_FONTSIZE,
    frameon=False,
    markerscale=1.5,
    handlelength=1,
)

# ── Save ───────────────────────────────────────────────────────────────────────
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_SVG, format="svg", bbox_inches="tight")
print(f"Saved {OUT_SVG}")
plt.close()
