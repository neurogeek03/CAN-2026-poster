"""
Test plot: apply hemisphere_params.json (line split + rotation) and show
one hemisphere per sample side by side.

Usage:
    uv run python scripts/figures/fig_hemisphere_test_v1.py
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT   = Path(__file__).resolve().parents[2]
CSV    = ROOT / "data" / "processed" / "hemisphere_coords.csv"
PARAMS = ROOT / "conf" / "hemisphere_params.json"
OUT    = ROOT / "figures" / "fig_hemisphere_test_v1.svg"


def line_split_mask(x, y, pt1, pt2, side):
    x1, y1 = pt1
    x2, y2 = pt2
    cross = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    return (cross * side) > 0


def rotate(x, y, angle_deg):
    cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    xr = cx + c * (x - cx) - s * (y - cy)
    yr = cy + s * (x - cx) + c * (y - cy)
    return xr, yr


def main():
    df = pd.read_csv(CSV)
    with open(PARAMS) as f:
        params = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor="white")

    for ax, (label, key) in zip(axes, [("OIL", "oil"), ("CORT", "cort")]):
        p = params[key]
        ls = params["line_split"][key]
        sid = p["sample"]

        sub = df[df["sample"] == sid]
        x, y = sub["x"].to_numpy(float), sub["y"].to_numpy(float)

        mask = line_split_mask(x, y, ls["pt1"], ls["pt2"], ls["side"])
        x, y = x[mask], y[mask]
        x, y = rotate(x, y, p["rotation_deg"])

        ax.scatter(x, y, s=0.5, c="#AAAAAA", linewidths=0, rasterized=True)
        ax.set_title(f"{label}  ({sid})", fontsize=12, fontweight="bold")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    plt.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"Saved → {OUT}")


if __name__ == "__main__":
    main()
