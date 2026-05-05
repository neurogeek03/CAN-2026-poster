"""
Test plot: apply hemisphere_params.json (line split + rotation + scale + alignment)
and show both hemispheres side-by-side in a single panel.

This mirrors the exact transform used in select_hemisphere.py's AlignSession so the
output should match what was saved interactively.

Usage:
    uv run python scripts/figures/fig_hemisphere_test_v2.py
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT   = Path(__file__).resolve().parents[2]
CSV    = ROOT / "data" / "processed" / "hemisphere_coords.csv"
PARAMS = ROOT / "conf" / "hemisphere_params.json"
OUT    = ROOT / "figures" / "fig_hemisphere_test_v2.svg"


def line_split_mask(x, y, pt1, pt2, side):
    x1, y1 = pt1
    x2, y2 = pt2
    cross = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    return (cross * side) > 0


def rotate(x, y, angle_deg):
    cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return cx + c*(x-cx) - s*(y-cy), cy + s*(x-cx) + c*(y-cy)


def flip_x(x):
    return 2 * x.mean() - x


def place(x, y, rot, do_flip, scale, x_anchor, y_shift):
    """Mirror of AlignSession._place() in select_hemisphere.py."""
    x, y = rotate(x, y, rot)
    if do_flip:
        x = flip_x(x)
    x = x * scale
    y = y * scale
    if x_anchor >= 0:   # CORT on right: left edge → x_anchor
        x = x - x.min() + x_anchor
    else:               # OIL on left: right edge → x_anchor
        x = x - x.max() + x_anchor
    y = y - y.mean() + y_shift
    return x, y


def load_hemisphere(df, params, key):
    p   = params[key]
    ls  = params["line_split"][key]
    sub = df[df["sample"] == p["sample"]]
    x, y = sub["x"].to_numpy(float), sub["y"].to_numpy(float)
    mask = line_split_mask(x, y, ls["pt1"], ls["pt2"], ls["side"])
    return x[mask], y[mask]


def main():
    with open(PARAMS) as f:
        params = json.load(f)

    df = pd.read_csv(CSV)

    x_oil,  y_oil  = load_hemisphere(df, params, "oil")
    x_cort, y_cort = load_hemisphere(df, params, "cort")

    p_oil  = params["oil"]
    p_cort = params["cort"]
    aln    = params["alignment"]
    gap    = aln["x_gap_um"]
    yoff   = aln.get("y_offset_um", 0.0)

    x_oil_p, y_oil_p = place(
        x_oil, y_oil,
        rot=p_oil["rotation_deg"],
        do_flip=p_oil["flip_x"],
        scale=p_oil["scale"],
        x_anchor=-gap / 2,
        y_shift=0,
    )
    x_cort_p, y_cort_p = place(
        x_cort, y_cort,
        rot=p_cort["rotation_deg"],
        do_flip=p_cort["flip_x"],
        scale=p_cort["scale"],
        x_anchor=+gap / 2,
        y_shift=yoff,
    )

    fig, ax = plt.subplots(figsize=(12, 7), facecolor="white")
    ax.scatter(x_oil_p,  y_oil_p,  s=0.5, c="#4393C3", linewidths=0,
               rasterized=True, label="OIL")
    ax.scatter(x_cort_p, y_cort_p, s=0.5, c="#D6604D", linewidths=0,
               rasterized=True, label="CORT")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.legend(loc="upper right", markerscale=8, frameon=False, fontsize=9)
    ax.set_facecolor("white")

    plt.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"Saved → {OUT}")


if __name__ == "__main__":
    main()
