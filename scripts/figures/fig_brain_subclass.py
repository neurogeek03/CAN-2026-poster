"""
Figure: two full brains (OIL and CORT) colored by cell subclass, no legend.

Colors are loaded from data/cluster_annotation_term.csv (Allen Brain Atlas
canonical hex colors per subclass). Alignment is read from conf/brain_params.json
(set interactively via scripts/align_brains.py).

Output: figures/fig_brain_subclass.svg

Usage:
    uv run python scripts/figures/fig_brain_subclass.py
"""

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = "Liberation Sans"

ROOT         = Path(__file__).resolve().parents[2]
CSV          = ROOT / "data" / "processed" / "brain_subclass_coords.csv"
BRAIN_PARAMS = ROOT / "conf" / "brain_params.json"
COLOR_CSV    = ROOT / "data" / "cluster_annotation_term.csv"
OUT          = ROOT / "pdfs" / "fig_brain_subclass.pdf"

# ── Config ────────────────────────────────────────────────────────────────────
FIG_W        = 14       # figure width  (inches)
FIG_H        = 7        # figure height (inches)
WSPACE       = 0.04     # horizontal gap between panels
POINT_SIZE   = 1.5      # scatter marker size (pt²)
POINT_ALPHA  = 1.0      # marker opacity
FALLBACK_COLOR = "#aaaaaa"  # color for subclasses missing from atlas
TITLE_FONTSIZE = 30     # panel title font size (pt)
TITLE_FONTWEIGHT = "bold"
# ─────────────────────────────────────────────────────────────────────────────


def _normalize(name: str) -> str:
    """Convert atlas subclass names to the format used in brain_subclass_coords."""
    name = re.sub(r"^\d+ ", "", name)          # strip leading "001 "
    return name.replace("/", "_").replace("-", "_").replace(" ", "_")


def load_color_map() -> dict[str, str]:
    df = pd.read_csv(COLOR_CSV, usecols=["name", "cluster_annotation_term_set_name", "color_hex_triplet"])
    sub = df[df["cluster_annotation_term_set_name"] == "subclass"]
    return {_normalize(row["name"]): row["color_hex_triplet"] for _, row in sub.iterrows()}


def rotate(x, y, angle_deg):
    cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return cx + c*(x-cx) - s*(y-cy), cy + s*(x-cx) + c*(y-cy)


def apply_alignment(df_oil: pd.DataFrame, df_cort: pd.DataFrame,
                    params: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    gap  = params["alignment"]["x_gap_um"]
    yoff = params["alignment"]["y_offset_um"]

    def transform(df, rot):
        x = df["x"].to_numpy(float)
        y = df["y"].to_numpy(float)
        x, y = rotate(x, y, rot)
        x -= x.mean(); y -= y.mean()
        return x, y

    xo, yo = transform(df_oil,  params["oil"]["rotation_deg"])
    xc, yc = transform(df_cort, params["cort"]["rotation_deg"])

    # normalize both brains to the same bounding box so they appear the same size
    def _max_span(x, y):
        return max(x.max() - x.min(), y.max() - y.min())

    ref_span = _max_span(xo, yo)
    sc = ref_span / _max_span(xc, yc)
    xc, yc = xc * sc, yc * sc

    half_w_o = (xo.max() - xo.min()) / 2
    half_w_c = (xc.max() - xc.min()) / 2

    xo = xo - half_w_o - gap / 2
    xc = xc + half_w_c + gap / 2
    yc = yc + yoff

    df_oil  = df_oil.copy();  df_oil["x"]  = xo;  df_oil["y"]  = yo
    df_cort = df_cort.copy(); df_cort["x"] = xc;  df_cort["y"] = yc
    return df_oil, df_cort


def plot_brain(ax, df: pd.DataFrame, color_map: dict, title: str):
    for sub in sorted(df["cell_subclass"].unique()):
        mask = df["cell_subclass"] == sub
        ax.scatter(
            df.loc[mask, "x"],
            df.loc[mask, "y"],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            c=[color_map.get(sub, FALLBACK_COLOR)],
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
    with open(BRAIN_PARAMS) as f:
        params = json.load(f)

    df        = pd.read_csv(CSV)
    color_map = load_color_map()

    oil_id  = params["oil"]["sample"]
    cort_id = params["cort"]["sample"]

    df_oil  = df[df["sample"] == oil_id].reset_index(drop=True)
    df_cort = df[df["sample"] == cort_id].reset_index(drop=True)

    df_oil, df_cort = apply_alignment(df_oil, df_cort, params)

    fig, axes = plt.subplots(1, 2, figsize=(FIG_W, FIG_H), facecolor="white")
    fig.subplots_adjust(wspace=WSPACE, left=0.02, right=0.98, top=0.93, bottom=0.02)

    plot_brain(axes[0], df_oil,  color_map, "OIL")
    plot_brain(axes[1], df_cort, color_map, "CORT")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, format="pdf", bbox_inches="tight", facecolor="white", dpi=600)
    print(f"Saved → {OUT}")


if __name__ == "__main__":
    main()
