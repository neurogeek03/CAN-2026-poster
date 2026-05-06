"""
Interactive rectangular crop tool for the aligned hemisphere view.

Renders both hemispheres using the SAME coordinate system as
fig_spatial_expression_v3.py (line-split → rotate → straighten cut edge →
gap anchor), then lets you define a rectangular crop with four sliders.
Press Save to write the crop bounds into conf/hemisphere_params.json under
the "crop" key.

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


def straighten_cut_edge(x, y, hkey):
    """
    Shear-correct so the inner cut edge is perfectly vertical.
    Returns (x_corrected, y_corrected, cut_edge_x).
    Mirrors _straighten_cut_edge in fig_spatial_expression_v3.py.
    """
    edge_pct = 95 if hkey == "oil" else 5
    n_bins = 30
    edges = np.percentile(y, np.linspace(0, 100, n_bins + 1))
    x_edge, y_mid = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (y >= lo) & (y <= hi)
        if m.sum() >= 5:
            x_edge.append(np.percentile(x[m], edge_pct))
            y_mid.append((lo + hi) / 2)

    if len(x_edge) < 3:
        return x, y, (x.max() if hkey == "oil" else x.min())

    slope = np.polyfit(y_mid, x_edge, 1)[0]
    y_ref = float(np.median(y))
    x_corr = x - slope * (y - y_ref)
    cut_edge_x = float(np.mean(x_edge) + slope * (y_ref - np.mean(y_mid)))
    return x_corr, y, cut_edge_x


def load_and_place_hemisphere(df, params, hkey):
    """
    Apply the full pipeline (line-split → rotate → straighten → gap anchor)
    matching fig_spatial_expression_v3.py exactly.
    Returns (x, y) in the figure coordinate space.
    """
    p  = params[hkey]
    ls = params["line_split"][hkey]
    aln = params.get("alignment", {})
    gap  = aln.get("x_gap_um",    670.0)
    yoff = aln.get("y_offset_um", 0.0)

    sub = df[df["sample"] == p["sample"]]
    x = sub["x"].to_numpy(float)
    y = sub["y"].to_numpy(float)

    mask = line_split_mask(x, y, ls["pt1"], ls["pt2"], ls["side"])
    x, y = x[mask], y[mask]

    x, y = rotate(x, y, p["rotation_deg"])

    scale = p.get("scale", 1.0)
    x = x * scale
    y = y * scale

    x, y, cut_edge_x = straighten_cut_edge(x, y, hkey)

    x_anchor = gap / 2 if hkey == "cort" else -gap / 2
    y_shift  = yoff    if hkey == "cort" else 0.0
    x += -cut_edge_x + x_anchor
    y += -y.mean() + y_shift

    return x, y


def subsample(x, y, n):
    if len(x) > n:
        idx = np.random.choice(len(x), n, replace=False)
        return x[idx], y[idx]
    return x, y


def main():
    with open(PARAMS) as f:
        params = json.load(f)

    df = pd.read_csv(CSV)

    x_oil,  y_oil  = load_and_place_hemisphere(df, params, "oil")
    x_cort, y_cort = load_and_place_hemisphere(df, params, "cort")

    all_x = np.concatenate([x_oil, x_cort])
    all_y = np.concatenate([y_oil, y_cort])
    x_data_min, x_data_max = all_x.min(), all_x.max()
    y_data_min, y_data_max = all_y.min(), all_y.max()
    x_pad = (x_data_max - x_data_min) * 0.05
    y_pad = (y_data_max - y_data_min) * 0.05

    # Subsample for display
    xo, yo = subsample(x_oil,  y_oil,  PREVIEW_N)
    xc, yc = subsample(x_cort, y_cort, PREVIEW_N)

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

    SHADE = dict(facecolor="#999999", alpha=0.25, linewidth=0, zorder=5)
    r_top    = ax.add_patch(Rectangle((0, 0), 1, 1, **SHADE))
    r_bottom = ax.add_patch(Rectangle((0, 0), 1, 1, **SHADE))
    r_left   = ax.add_patch(Rectangle((0, 0), 1, 1, **SHADE))
    r_right  = ax.add_patch(Rectangle((0, 0), 1, 1, **SHADE))

    margin = 1.10

    def _shade_bounds():
        xl = x_data_min - x_pad * margin
        xr = x_data_max + x_pad * margin
        yb = y_data_min - y_pad * margin
        yt = y_data_max + y_pad * margin
        L, R, B, T = sl_left.val, sl_right.val, sl_bottom.val, sl_top.val
        r_left.set_bounds(xl, yb, L - xl, yt - yb)
        r_right.set_bounds(R, yb, xr - R, yt - yb)
        r_bottom.set_bounds(L, yb, R - L, B - yb)
        r_top.set_bounds(L, T, R - L, yt - T)

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
        L, R = sl_left.val, max(sl_right.val, sl_left.val + x_step)
        B, T = sl_bottom.val, max(sl_top.val, sl_bottom.val + y_step)
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
