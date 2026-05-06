"""
Interactive rectangular crop tool for the aligned hemisphere view.

Renders both hemispheres using the alignment already saved in
conf/hemisphere_params.json, then lets you define a rectangular
crop region with four sliders. Press Save to write the crop bounds
back into the JSON under the "crop" key.

Usage:
    uv run python scripts/crop_hemispheres.py
"""

import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

_backend = matplotlib.get_backend()
if _backend.lower() in ("agg", "svg", "pdf", "ps", "cairo"):
    sys.exit(
        f"Non-interactive backend ({_backend}). Run locally, not on a headless server."
    )

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, Slider

ROOT   = Path(__file__).resolve().parents[1]
CSV    = ROOT / "data" / "processed" / "hemisphere_coords.csv"
PARAMS = ROOT / "conf" / "hemisphere_params.json"

PREVIEW_N = 6000


def line_split_mask(x, y, pt1, pt2, side):
    x1, y1 = pt1
    x2, y2 = pt2
    cross = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    return (cross * side) > 0


def rotate(x, y, angle_deg):
    cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return cx + c * (x - cx) - s * (y - cy), cy + s * (x - cx) + c * (y - cy)


def load_hemisphere(df, params, key):
    p  = params[key]
    ls = params["line_split"][key]
    sub = df[df["sample"] == p["sample"]]
    x, y = sub["x"].to_numpy(float), sub["y"].to_numpy(float)
    mask = line_split_mask(x, y, ls["pt1"], ls["pt2"], ls["side"])
    x, y = x[mask], y[mask]
    x, y = rotate(x, y, p["rotation_deg"])
    x -= x.mean()
    y -= y.mean()
    return x, y


def subsample(x, y, n):
    if len(x) > n:
        idx = np.random.choice(len(x), n, replace=False)
        return x[idx], y[idx]
    return x, y


def apply_alignment(x_oil, y_oil, x_cort, y_cort, aln):
    so   = aln.get("scale_oil",   1.0)
    sc   = aln.get("scale_cort",  1.0)
    gap  = aln.get("x_gap_um",    670.0)
    yoff = aln.get("y_offset_um", 0.0)

    oil_half_w  = (x_oil.max()  - x_oil.min())  * so / 2
    cort_half_w = (x_cort.max() - x_cort.min()) * sc / 2

    x_o = x_oil  * so - oil_half_w  - gap / 2
    y_o = y_oil  * so
    x_c = x_cort * sc + cort_half_w + gap / 2
    y_c = y_cort * sc + yoff
    return x_o, y_o, x_c, y_c


def main():
    with open(PARAMS) as f:
        params = json.load(f)

    df = pd.read_csv(CSV)

    x_oil,  y_oil  = load_hemisphere(df, params, "oil")
    x_cort, y_cort = load_hemisphere(df, params, "cort")

    aln = params.get("alignment", {})
    xo_full, yo_full, xc_full, yc_full = apply_alignment(
        x_oil, y_oil, x_cort, y_cort, aln
    )

    # subsample for display
    xo, yo = subsample(xo_full, yo_full, PREVIEW_N)
    xc, yc = subsample(xc_full, yc_full, PREVIEW_N)

    all_x = np.concatenate([xo_full, xc_full])
    all_y = np.concatenate([yo_full, yc_full])
    x_data_min, x_data_max = all_x.min(), all_x.max()
    y_data_min, y_data_max = all_y.min(), all_y.max()
    x_pad = (x_data_max - x_data_min) * 0.05
    y_pad = (y_data_max - y_data_min) * 0.05

    # Load any previously saved crop, or default to full extent
    saved_crop = params.get("crop", {})
    init_left   = saved_crop.get("x_min", x_data_min)
    init_right  = saved_crop.get("x_max", x_data_max)
    init_bottom = saved_crop.get("y_min", y_data_min)
    init_top    = saved_crop.get("y_max", y_data_max)

    # ── figure layout ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 9))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.46)
    fig.suptitle("Crop hemispheres — adjust sliders, then Save", fontsize=11)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_facecolor("white")

    ax.scatter(xo, yo, s=0.8, c="#4393C3", linewidths=0, rasterized=True, label="OIL")
    ax.scatter(xc, yc, s=0.8, c="#D6604D", linewidths=0, rasterized=True, label="CORT")
    ax.legend(loc="upper right", markerscale=8, frameon=False, fontsize=9)

    # Crop rectangle
    rect = Rectangle(
        (init_left, init_bottom),
        init_right - init_left,
        init_top   - init_bottom,
        linewidth=1.5,
        edgecolor="#222222",
        facecolor="none",
        linestyle="--",
        zorder=10,
    )
    ax.add_patch(rect)

    # Grey shading outside crop (4 rectangles — top, bottom, left, right)
    SHADE = dict(facecolor="#999999", alpha=0.25, linewidth=0, zorder=5)
    r_top    = ax.add_patch(Rectangle((0, 0, ), 1, 1, **SHADE))
    r_bottom = ax.add_patch(Rectangle((0, 0, ), 1, 1, **SHADE))
    r_left   = ax.add_patch(Rectangle((0, 0, ), 1, 1, **SHADE))
    r_right  = ax.add_patch(Rectangle((0, 0, ), 1, 1, **SHADE))

    margin = 1.10   # how far the shading extends beyond the data

    def _shade_bounds():
        xl, xr = x_data_min - x_pad * margin, x_data_max + x_pad * margin
        yb, yt = y_data_min - y_pad * margin, y_data_max + y_pad * margin
        L = sl_left.val
        R = sl_right.val
        B = sl_bottom.val
        T = sl_top.val
        r_left.set_bounds(xl, yb, L - xl, yt - yb)
        r_right.set_bounds(R, yb, xr - R, yt - yb)
        r_bottom.set_bounds(L, yb, R - L, B - yb)
        r_top.set_bounds(L, T, R - L, yt - T)

    # fixed view — show all data with padding
    ax.set_xlim(x_data_min - x_pad, x_data_max + x_pad)
    ax.set_ylim(y_data_min - y_pad, y_data_max + y_pad)

    # ── sliders ────────────────────────────────────────────────────────────────
    x_range = x_data_max - x_data_min
    y_range = y_data_max - y_data_min
    x_step  = round(x_range / 400, 1)
    y_step  = round(y_range / 400, 1)

    def make_slider(left, bottom, label, vmin, vmax, vinit, step):
        a = fig.add_axes([left, bottom, 0.84, 0.03])
        return Slider(a, label, vmin, vmax, valinit=vinit, valstep=step)

    sl_left   = make_slider(0.08, 0.38, "Left   (x_min)",
                            x_data_min - x_pad, x_data_max + x_pad, init_left,   x_step)
    sl_right  = make_slider(0.08, 0.31, "Right  (x_max)",
                            x_data_min - x_pad, x_data_max + x_pad, init_right,  x_step)
    sl_bottom = make_slider(0.08, 0.24, "Bottom (y_min)",
                            y_data_min - y_pad, y_data_max + y_pad, init_bottom, y_step)
    sl_top    = make_slider(0.08, 0.17, "Top    (y_max)",
                            y_data_min - y_pad, y_data_max + y_pad, init_top,    y_step)

    def update(_val=None):
        L = sl_left.val
        R = sl_right.val
        B = sl_bottom.val
        T = sl_top.val
        # clamp so rectangle is never inverted
        R = max(R, L + x_step)
        T = max(T, B + y_step)
        rect.set_xy((L, B))
        rect.set_width(R - L)
        rect.set_height(T - B)
        _shade_bounds()
        fig.canvas.draw_idle()

    sl_left.on_changed(update)
    sl_right.on_changed(update)
    sl_bottom.on_changed(update)
    sl_top.on_changed(update)
    update()

    # ── reset button ───────────────────────────────────────────────────────────
    ax_reset = fig.add_axes([0.20, 0.08, 0.18, 0.04])
    btn_reset = Button(ax_reset, "Reset crop")

    def reset(_event):
        sl_left.set_val(x_data_min)
        sl_right.set_val(x_data_max)
        sl_bottom.set_val(y_data_min)
        sl_top.set_val(y_data_max)

    btn_reset.on_clicked(reset)

    # ── save button ────────────────────────────────────────────────────────────
    ax_save = fig.add_axes([0.62, 0.08, 0.18, 0.04])
    btn_save = Button(ax_save, "Save crop")

    def save(_event):
        params["crop"] = {
            "x_min": float(sl_left.val),
            "x_max": float(sl_right.val),
            "y_min": float(sl_bottom.val),
            "y_max": float(sl_top.val),
        }
        with open(PARAMS, "w") as f:
            json.dump(params, f, indent=2)
        print(f"Saved crop → {PARAMS}")
        print(f"  x: [{sl_left.val:.1f}, {sl_right.val:.1f}]")
        print(f"  y: [{sl_bottom.val:.1f}, {sl_top.val:.1f}]")
        plt.close(fig)

    btn_save.on_clicked(save)
    plt.show(block=True)


if __name__ == "__main__":
    main()
