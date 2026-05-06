"""
Interactive alignment tool — scale and position the two hemispheres
relative to each other. Updates conf/hemisphere_params.json on Save.

Usage:
    uv run python scripts/align_hemispheres.py

Requirements: active display (run locally).
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
from matplotlib.widgets import Button, Slider

ROOT   = Path(__file__).resolve().parents[1]
CSV    = ROOT / "data" / "processed" / "hemisphere_coords.csv"
PARAMS = ROOT / "conf" / "hemisphere_params.json"

PREVIEW_N = 5000


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


def load_hemisphere(df, params, key):
    p  = params[key]
    ls = params["line_split"][key]
    sub = df[df["sample"] == p["sample"]]
    x, y = sub["x"].to_numpy(float), sub["y"].to_numpy(float)
    mask = line_split_mask(x, y, ls["pt1"], ls["pt2"], ls["side"])
    x, y = x[mask], y[mask]
    x, y = rotate(x, y, p["rotation_deg"])
    # Centre at origin
    x -= x.mean()
    y -= y.mean()
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

    x_oil,  y_oil  = load_hemisphere(df, params, "oil")
    x_cort, y_cort = load_hemisphere(df, params, "cort")

    # Subsample for smooth interaction
    xo, yo = subsample(x_oil,  y_oil,  PREVIEW_N)
    xc, yc = subsample(x_cort, y_cort, PREVIEW_N)

    # Current saved alignment (defaults if missing)
    aln = params.get("alignment", {})
    init_scale_oil  = aln.get("scale_oil",  1.0)
    init_scale_cort = aln.get("scale_cort", 1.0)
    init_gap        = aln.get("x_gap_um",   670.0)
    init_y_offset   = aln.get("y_offset_um", 0.0)

    # ── figure layout ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 9))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.42)
    fig.suptitle("Align hemispheres — adjust sliders, then Save", fontsize=11)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_facecolor("white")

    sc_oil  = ax.scatter([], [], s=0.8, c="#4393C3", linewidths=0, rasterized=True, label="OIL")
    sc_cort = ax.scatter([], [], s=0.8, c="#D6604D", linewidths=0, rasterized=True, label="CORT")
    ax.legend(loc="upper right", markerscale=8, frameon=False, fontsize=9)

    # ── sliders ───────────────────────────────────────────────────────────────
    span = max(x_oil.max()-x_oil.min(), x_cort.max()-x_cort.min())

    def make_slider(left, bottom, label, vmin, vmax, vinit, step=None):
        a = fig.add_axes([left, bottom, 0.35, 0.03])
        kw = dict(valmin=vmin, valmax=vmax, valinit=vinit)
        if step:
            kw["valstep"] = step
        return Slider(a, label, **kw)

    sl_scale_oil  = make_slider(0.08, 0.29, "OIL scale",   0.5, 2.0, init_scale_oil,  0.01)
    sl_scale_cort = make_slider(0.57, 0.29, "CORT scale",  0.5, 2.0, init_scale_cort, 0.01)

    ax_gap = fig.add_axes([0.08, 0.22, 0.84, 0.03])
    sl_gap = Slider(ax_gap, "Gap (µm)", 0, span * 1.5, valinit=init_gap, valstep=10)

    ax_yoff = fig.add_axes([0.08, 0.15, 0.84, 0.03])
    sl_yoff = Slider(ax_yoff, "Y offset (µm)", -span*0.5, span*0.5, valinit=init_y_offset, valstep=10)

    ax_zoom = fig.add_axes([0.08, 0.08, 0.84, 0.03])
    sl_zoom = Slider(ax_zoom, "Zoom", 0.5, 8.0, valinit=1.0, valstep=0.1)

    # ── update ────────────────────────────────────────────────────────────────
    def update(_val=None):
        so = sl_scale_oil.val
        sc = sl_scale_cort.val
        gap = sl_gap.val
        yoff = sl_yoff.val
        zoom = sl_zoom.val

        # OIL on left, CORT on right — use full-resolution spans for positioning
        # so a heavily-masked CORT still lands at the right x_anchor
        oil_half_w  = (x_oil.max()  - x_oil.min())  * so / 2
        cort_half_w = (x_cort.max() - x_cort.min()) * sc / 2

        x_oil_plot  = xo * so - oil_half_w - gap/2
        y_oil_plot  = yo * so

        x_cort_plot = xc * sc + cort_half_w + gap/2
        y_cort_plot = yc * sc + yoff

        sc_oil.set_offsets(np.column_stack([x_oil_plot, y_oil_plot]))
        sc_cort.set_offsets(np.column_stack([x_cort_plot, y_cort_plot]))

        # View: center on combined extent, apply zoom
        all_x = np.concatenate([x_oil_plot, x_cort_plot])
        all_y = np.concatenate([y_oil_plot, y_cort_plot])
        cx = (all_x.min() + all_x.max()) / 2
        cy = (all_y.min() + all_y.max()) / 2
        half_w = (all_x.max() - all_x.min()) * 0.55 / zoom
        half_h = (all_y.max() - all_y.min()) * 0.55 / zoom
        ax.set_xlim(cx - half_w, cx + half_w)
        ax.set_ylim(cy - half_h, cy + half_h)
        fig.canvas.draw_idle()

    sl_scale_oil.on_changed(update)
    sl_scale_cort.on_changed(update)
    sl_gap.on_changed(update)
    sl_yoff.on_changed(update)
    sl_zoom.on_changed(update)
    update()

    # ── save button ───────────────────────────────────────────────────────────
    ax_save = fig.add_axes([0.40, 0.02, 0.20, 0.04])
    btn = Button(ax_save, "Save params")

    def save(_event):
        params["alignment"] = {
            "scale_oil":    float(sl_scale_oil.val),
            "scale_cort":   float(sl_scale_cort.val),
            "x_gap_um":     float(sl_gap.val),
            "y_offset_um":  float(sl_yoff.val),
        }
        with open(PARAMS, "w") as f:
            json.dump(params, f, indent=2)
        print(f"Saved → {PARAMS}")
        print(f"  scale_oil={sl_scale_oil.val:.2f}  scale_cort={sl_scale_cort.val:.2f}")
        print(f"  gap={sl_gap.val:.0f} µm   y_offset={sl_yoff.val:.0f} µm")
        plt.close(fig)

    btn.on_clicked(save)
    plt.show(block=True)


if __name__ == "__main__":
    main()
