"""
Figure: two brain sections (BC28 = OIL, BC3 = CORT) colored by cell subclass,
stacked in a single column (2 rows × 1 col), no legend.

Colors loaded from data/cluster_annotation_term.csv (Allen Brain Atlas hex colors).

Input:  data/Spatial/slide_tags/coords_BC28_BC3_score0.5.csv
Output: pdfs/fig_brain_subclass_column.pdf

Usage:
    uv run python scripts/figures/fig_brain_subclass_column.py
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.rcParams["font.family"] = "Liberation Sans"

ROOT      = Path(__file__).resolve().parents[2]
CSV       = ROOT / "data" / "Spatial" / "slide_tags" / "coords_BC28_BC3_score0.5.csv"
COLOR_CSV = ROOT / "data" / "cluster_annotation_term.csv"
OUT       = ROOT / "pdfs" / "fig_brain_subclass_column.pdf"

FIG_W          = 7
FIG_H          = 14
HSPACE         = 0.04
POINT_SIZE     = 1.5
POINT_ALPHA    = 1.0
FALLBACK_COLOR = "#aaaaaa"
TITLE_FONTSIZE = 30
TITLE_FONTWEIGHT = "bold"

SAMPLE_TITLES = {"BC28": "OIL", "BC3": "CORT"}


def _normalize(name: str) -> str:
    name = re.sub(r"^\d+ ", "", name)
    return name.replace("/", "_").replace("-", "_").replace(" ", "_")


def load_color_map() -> dict[str, str]:
    df = pd.read_csv(COLOR_CSV, usecols=["name", "cluster_annotation_term_set_name", "color_hex_triplet"])
    sub = df[df["cluster_annotation_term_set_name"] == "subclass"]
    return {_normalize(row["name"]): row["color_hex_triplet"] for _, row in sub.iterrows()}


def plot_section(ax, df: pd.DataFrame, color_map: dict, title: str):
    df = df.dropna(subset=["subclass_name"])
    for sub in sorted(df["subclass_name"].unique()):
        mask = df["subclass_name"] == sub
        color = color_map.get(_normalize(sub), FALLBACK_COLOR)
        ax.scatter(
            df.loc[mask, "x_um"],
            df.loc[mask, "y_um"],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            c=[color],
            linewidths=0,
            rasterized=True,
        )

    ax.set_aspect("equal")
    ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight=TITLE_FONTWEIGHT, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_facecolor("white")


def main():
    df = pd.read_csv(CSV)
    color_map = load_color_map()

    fig, axes = plt.subplots(2, 1, figsize=(FIG_W, FIG_H), facecolor="white")
    fig.subplots_adjust(hspace=HSPACE, left=0.02, right=0.98, top=0.97, bottom=0.02)

    for ax, sample_id in zip(axes, ["BC28", "BC3"]):
        sub_df = df[df["sample"] == sample_id].reset_index(drop=True)
        plot_section(ax, sub_df, color_map, SAMPLE_TITLES[sample_id])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, format="pdf", bbox_inches="tight", facecolor="white", dpi=600)
    print(f"Saved → {OUT}")


if __name__ == "__main__":
    main()
