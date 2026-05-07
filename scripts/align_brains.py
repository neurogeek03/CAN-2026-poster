"""
Interactive alignment tool — rotate and position two brain sections
relative to each other. Saves alignment params to a JSON file in conf/.

Usage:
    uv run python scripts/align_brains.py                # Slide-seq brains
    uv run python scripts/align_brains.py --slidetags    # Slide-tags BC28/BC3

Requirements: active display (run locally).
"""

import argparse
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

ROOT = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser()
parser.add_argument("--slidetags", action="store_true",
                    help="Align Slide-tags sections (BC28=OIL, BC3=CORT)")
args = parser.parse_args()

if args.slidetags:
    CSV        = ROOT / "data" / "Spatial" / "slide_tags" / "coords_BC28_BC3_score0.5.csv"
    PARAMS     = ROOT / "conf" / "slidetags_coords.json"
    X_COL      = "x_um"
    Y_COL      = "y_um"
    OIL_SAMPLE  = "BC28"
    CORT_SAMPLE = "BC3"
else:
    CSV        = ROOT / "data" / "processed" / "hemisphere_coords.csv"
    PARAMS     = ROOT / "conf" / "brain_params.json"
    X_COL      = "x"
    Y_COL      = "y"
    OIL_SAMPLE  = "B14"
    CORT_SAMPLE = "B03"

PREVIEW_N = 5000

print(f"Mode     : {'slide-tags (BC28/BC3)' if args.slidetags else 'slide-seq (B14/B03)'}")
print(f"CSV      : {CSV}")
print(f"Params   : {PARAMS}")


def rotate(x, y, angle_deg):
    cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return cx + c*(x-cx) - s*(y-cy), cy + s*(x-cx) + c*(y-cy)


def subsample(x, y, n):
    if len(x) > n:
        idx = np.random.choice(len(x), n, replace=False)
        return x[idx], y[idx]
    return x, y


def main():
    # Load or initialise params
    if PARAMS.exists():
        with open(PARAMS) as f:
            params = json.load(f)
    else:
        params = {}

    df = pd.read_csv(CSV)

    sub_oil  = df[df["sample"] == OIL_SAMPLE]
    sub_cort = df[df["sample"] == CORT_SAMPLE]

    x_oil_raw  = sub_oil[X_COL].to_numpy(float)
    y_oil_raw  = sub_oil[Y_COL].to_numpy(float)
    x_cort_raw = sub_cort[X_COL].to_numpy(float)
    y_cort_raw = sub_cort[Y_COL].to_numpy(float)

    # Centre each brain at origin before subsampling (rotation applied in update)
    x_oil_raw  -= x_oil_raw.mean()
    y_oil_raw  -= y_oil_raw.mean()
    x_cort_raw -= x_cort_raw.mean()
    y_cort_raw -= y_cort_raw.mean()

    xo, yo = subsample(x_oil_raw,  y_oil_raw,  PREVIEW_N)
    xc, yc = subsample(x_cort_raw, y_cort_raw, PREVIEW_N)

    # Defaults from saved params
    aln = params.get("alignment", {})
    init_rot_oil  = params.get("oil",  {}).get("rotation_deg", 0.0)
    init_rot_cort = params.get("cort", {}).get("rotation_deg", 0.0)
    init_gap      = aln.get("x_gap_um",    670.0)
    init_yoff     = aln.get("y_offset_um",   0.0)

    span = max(x_oil_raw.max() - x_oil_raw.min(),
               x_cort_raw.max() - x_cort_raw.min())

    # ── figure layout ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.35)
    mode_label = "Slide-tags BC28/BC3" if args.slidetags else "Slide-seq B14/B03"
    fig.suptitle(f"Align brains [{mode_label}] — rotate with sliders, then Save", fontsize=11)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_facecolor("white")

    sc_oil  = ax.scatter([], [], s=0.8, c="#4393C3", linewidths=0, rasterized=True, label="OIL")
    sc_cort = ax.scatter([], [], s=0.8, c="#D6604D", linewidths=0, rasterized=True, label="CORT")
    ax.legend(loc="upper right", markerscale=8, frameon=False, fontsize=9)

    # ── sliders ────────────────────────────────────────────────────────────────
    def make_slider(left, bottom, width, label, vmin, vmax, vinit, step=None):
        a = fig.add_axes([left, bottom, width, 0.03])
        kw = dict(valmin=vmin, valmax=vmax, valinit=vinit)
        if step:
            kw["valstep"] = step
        return Slider(a, label, **kw)

    sl_rot_oil  = make_slider(0.08, 0.28, 0.35, "OIL rotation (°)",  -180, 180, init_rot_oil,  0.5)
    sl_rot_cort = make_slider(0.57, 0.28, 0.35, "CORT rotation (°)", -180, 180, init_rot_cort, 0.5)

    sl_gap  = make_slider(0.08, 0.18, 0.84, "Gap (µm)",       0,         span * 1.5, init_gap,  10)
    sl_yoff = make_slider(0.08, 0.10, 0.84, "Y offset (µm)", -span * 0.5, span * 0.5, init_yoff, 10)

    # ── update ─────────────────────────────────────────────────────────────────
    def update(_val=None):
        rot_oil  = sl_rot_oil.val
        rot_cort = sl_rot_cort.val
        gap      = sl_gap.val
        yoff     = sl_yoff.val

        rx_o, ry_o = rotate(xo, yo, rot_oil)
        rx_c, ry_c = rotate(xc, yc, rot_cort)

        # Re-centre after rotation
        rx_o -= rx_o.mean(); ry_o -= ry_o.mean()
        rx_c -= rx_c.mean(); ry_c -= ry_c.mean()

        half_w_o = (rx_o.max() - rx_o.min()) / 2
        half_w_c = (rx_c.max() - rx_c.min()) / 2

        x_oil_plot  = rx_o - half_w_o - gap / 2
        y_oil_plot  = ry_o

        x_cort_plot = rx_c + half_w_c + gap / 2
        y_cort_plot = ry_c + yoff

        sc_oil.set_offsets(np.column_stack([x_oil_plot, y_oil_plot]))
        sc_cort.set_offsets(np.column_stack([x_cort_plot, y_cort_plot]))

        all_x = np.concatenate([x_oil_plot, x_cort_plot])
        all_y = np.concatenate([y_oil_plot, y_cort_plot])
        pad_x = (all_x.max() - all_x.min()) * 0.05
        pad_y = (all_y.max() - all_y.min()) * 0.05
        ax.set_xlim(all_x.min() - pad_x, all_x.max() + pad_x)
        ax.set_ylim(all_y.min() - pad_y, all_y.max() + pad_y)
        fig.canvas.draw_idle()

    sl_rot_oil.on_changed(update)
    sl_rot_cort.on_changed(update)
    sl_gap.on_changed(update)
    sl_yoff.on_changed(update)
    update()

    # ── save button ────────────────────────────────────────────────────────────
    ax_save = fig.add_axes([0.40, 0.02, 0.20, 0.05])
    btn = Button(ax_save, "Save params")

    def save(_event):
        out = {
            "oil":  {"sample": OIL_SAMPLE,  "rotation_deg": float(sl_rot_oil.val)},
            "cort": {"sample": CORT_SAMPLE, "rotation_deg": float(sl_rot_cort.val)},
            "alignment": {
                "x_gap_um":    float(sl_gap.val),
                "y_offset_um": float(sl_yoff.val),
            },
        }
        PARAMS.parent.mkdir(parents=True, exist_ok=True)
        with open(PARAMS, "w") as f:
            json.dump(out, f, indent=2)
        print(f"Saved → {PARAMS}")
        print(f"  rot_oil={sl_rot_oil.val:.1f}°  rot_cort={sl_rot_cort.val:.1f}°")
        print(f"  gap={sl_gap.val:.0f} µm   y_offset={sl_yoff.val:.0f} µm")
        plt.close(fig)

    btn.on_clicked(save)
    plt.show(block=True)


if __name__ == "__main__":
    main()
