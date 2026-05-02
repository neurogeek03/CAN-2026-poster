"""
Figure 01 v2 — UMAP colored by subclass_name, top 30 cell types labeled at centroids.

All cell types retain their original colors (no greying).
Top 30 by cell count get abbreviated labels placed at cluster centroids.
Label style: white text on a background matching the cluster color.
Overlapping labels are repelled using adjustText.

Input:  data/processed/umap_coords.csv
        data/cluster_annotation_term.csv
Output: figures/fig01_umap_v2_centroid_labels.svg

Dependencies: matplotlib, pandas, adjustText
  conda install -c conda-forge adjusttext
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from adjustText import adjust_text
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
UMAP_CSV  = Path("data/processed/umap_coords.csv")
COLOR_CSV = Path("data/cluster_annotation_term.csv")
OUT_SVG   = Path("figures/fig01_umap_v2_centroid_labels.svg")

# ── Config ─────────────────────────────────────────────────────────────────────
COLOR_COL      = "subclass_name"
N_LABELED      = 30         # top N cell types (by cell count) to receive labels
FORCE_LABEL    = [          # always label these regardless of rank
    "025 CA2-FC-IG Glut",
]
POINT_SIZE     = 0.8
POINT_ALPHA    = 0.7
FIG_W, FIG_H   = 12, 9
LABEL_FONTSIZE = 9
LEGEND_NCOL    = 2
LEGEND_FONTSIZE = 5.5

# ── Abbreviation rule ──────────────────────────────────────────────────────────
# Keep numeric prefix + next 2 words: "006 L4/5 IT CTX Glut" → "006 L4/5 IT"
def abbreviate(name: str, n_words: int = 3) -> str:
    return " ".join(name.split()[:n_words])

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(UMAP_CSV)

# ── Build color palette (subclass rows only) ───────────────────────────────────
color_df = pd.read_csv(
    COLOR_CSV,
    usecols=["name", "cluster_annotation_term_set_name", "color_hex_triplet"],
)
subclass_colors = (
    color_df[color_df["cluster_annotation_term_set_name"] == "subclass"]
    .set_index("name")["color_hex_triplet"]
    .to_dict()
)

# ── Sort categories by numeric prefix ─────────────────────────────────────────
def _num_prefix(label):
    try:
        return int(str(label).split(" ")[0])
    except ValueError:
        return 9999

present = df[COLOR_COL].dropna().unique().tolist()
present_sorted = sorted(present, key=_num_prefix)
df[COLOR_COL] = pd.Categorical(df[COLOR_COL], categories=present_sorted, ordered=True)
df = df.sort_values(COLOR_COL)

# ── Determine top N cell types to label ───────────────────────────────────────
top_n = (
    df[COLOR_COL]
    .value_counts()
    .head(N_LABELED)
    .index.tolist()
)
top_n_set = set(top_n) | set(FORCE_LABEL)

# ── Compute centroids for labeled cell types ───────────────────────────────────
centroids = (
    df[df[COLOR_COL].isin(top_n_set)]
    .groupby(COLOR_COL, observed=True)[["UMAP1", "UMAP2"]]
    .median()
)

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
        rasterized=True,
    )

ax.set_xlabel("UMAP1", fontsize=10)
ax.set_ylabel("UMAP2", fontsize=10)
ax.set_title(
    f"Slide-tags — {COLOR_COL} (n = {len(df):,}), top {N_LABELED} labeled",
    fontsize=11,
)
ax.axis("equal")
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# ── Place centroid labels ──────────────────────────────────────────────────────
texts = []
for subclass, row in centroids.iterrows():
    color = subclass_colors.get(subclass, "#808080")
    label = abbreviate(subclass)
    t = ax.text(
        row["UMAP1"],
        row["UMAP2"],
        label,
        fontsize=LABEL_FONTSIZE,
        color="white",
        fontweight="bold",
        ha="center",
        va="center",
        path_effects=[
            pe.withStroke(linewidth=0.6, foreground="black"),
        ],
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor=color,
            edgecolor="none",
            alpha=1.0,
        ),
    )
    texts.append(t)

# Repel overlapping labels; arrowprops draws thin leader lines when labels move
adjust_text(
    texts,
    ax=ax,
    expand=(1.3, 1.5),
    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
)

ax.legend_.remove() if ax.get_legend() else None

# ── Save ───────────────────────────────────────────────────────────────────────
OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_SVG, format="svg", bbox_inches="tight")
print(f"Saved {OUT_SVG}")
plt.close()
