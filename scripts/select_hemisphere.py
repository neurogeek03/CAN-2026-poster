"""
Interactive hemisphere selector — run once to define crop regions and alignment
transforms, saved to conf/hemisphere_params.json.

All downstream figure scripts load that JSON to apply the same crop + rotation.

Phase 1: For each sample (OIL then CORT), draw a lasso around the hemisphere
         you want to keep. Click "Confirm" to proceed.
Phase 2: Both cropped hemispheres shown side by side. Use sliders to set
         rotation and flip per sample, then click "Save".

Usage:
    uv run python scripts/select_hemisphere.py
    uv run python scripts/select_hemisphere.py --cort B03 --oil B14

Requirements: active display (run locally or via X forwarding / VNC).
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

# Fail fast if there's no display before loading heavy dependencies.
_backend = matplotlib.get_backend()
if _backend.lower() in ("agg", "svg", "pdf", "ps", "cairo"):
    import sys as _sys
    _sys.exit(
        f"Matplotlib selected a non-interactive backend ({_backend}).\n"
        "This script requires a display — run it locally, not on a headless server."
    )

import numpy as np
import pandas as pd
import yaml
from matplotlib.path import Path as MplPath
from matplotlib.widgets import Button, CheckButtons, Slider, LassoSelector

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"
COORDS_CSV  = ROOT / "data" / "processed" / "hemisphere_coords.csv"
OUTPUT_PATH = ROOT / "conf" / "hemisphere_params.json"

# Use a subsample size for the live alignment preview (sliders).
PREVIEW_N = 6000


# ── helpers ──────────────────────────────────────────────────────────────────

def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_coords(sample_id: str):
    """Return (x, y) arrays for a single sample from the pre-extracted CSV."""
    if not COORDS_CSV.exists():
        sys.exit(
            f"Coords CSV not found: {COORDS_CSV}\n"
            "Run on the HPC first:\n"
            "  uv run python scripts/extract/extract_hemisphere_coords.py"
        )
    df = pd.read_csv(COORDS_CSV)
    sub = df[df["sample"] == sample_id]
    if sub.empty:
        sys.exit(f"Sample '{sample_id}' not found in {COORDS_CSV}.")
    return sub["x"].to_numpy(float), sub["y"].to_numpy(float)


def apply_lasso(x, y, poly_verts):
    """Return boolean mask of points inside the lasso polygon."""
    path = MplPath(poly_verts)
    pts = np.column_stack([x, y])
    return path.contains_points(pts)


def rotate(x, y, angle_deg, cx=None, cy=None):
    """Rotate (x, y) around (cx, cy). Defaults to centroid."""
    if cx is None:
        cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    xr = cx + c * (x - cx) - s * (y - cy)
    yr = cy + s * (x - cx) + c * (y - cy)
    return xr, yr


def transform(x, y, rotation_deg, flip_x):
    """Apply flip then rotation (around centroid of current coords)."""
    if flip_x:
        cx = x.mean()
        x = 2 * cx - x
    x, y = rotate(x, y, rotation_deg)
    return x, y


# ── Phase 1: lasso selection ─────────────────────────────────────────────────

class LassoSession:
    """One interactive lasso window for a single sample."""

    def __init__(self, x, y, title):
        self.x = x
        self.y = y
        self.verts = None
        self._done = False

        # Subsample for display speed (scatter with >100 K pts can lag)
        if len(x) > 40_000:
            idx = np.random.choice(len(x), 40_000, replace=False)
            self._xd, self._yd = x[idx], y[idx]
        else:
            self._xd, self._yd = x, y

        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.fig.patch.set_facecolor("white")
        self.fig.suptitle(
            f"{title}\n"
            "Draw a lasso around the hemisphere to keep, then click Confirm.",
            fontsize=10,
        )
        self._sc_all = self.ax.scatter(
            self._xd, self._yd, s=0.8, c="#BBBBBB", linewidths=0, rasterized=True
        )
        self._sc_sel = self.ax.scatter(
            [], [], s=1.2, c="steelblue", linewidths=0, rasterized=True, zorder=3
        )
        self.ax.set_aspect("equal")
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)

        ax_btn = self.fig.add_axes([0.35, 0.01, 0.3, 0.05])
        self._btn = Button(ax_btn, "Confirm selection")
        self._btn.on_clicked(self._confirm)

        self._lasso = LassoSelector(self.ax, self._on_select, useblit=True)

    def _on_select(self, verts):
        self.verts = verts
        mask = apply_lasso(self._xd, self._yd, verts)
        self._sc_sel.set_offsets(np.column_stack([self._xd[mask], self._yd[mask]]))
        self.fig.canvas.draw_idle()

    def _confirm(self, _event):
        if self.verts is None:
            print("  Draw a lasso first.")
            return
        self._done = True
        plt.close(self.fig)

    def run(self):
        plt.show(block=True)
        if not self._done:
            sys.exit("Lasso window closed without confirming — exiting.")
        return self.verts


# ── Phase 2: alignment ────────────────────────────────────────────────────────

class AlignSession:
    """Side-by-side view with rotation sliders and flip toggles."""

    def __init__(self, x_oil, y_oil, x_cort, y_cort):
        # Subsample each for live preview
        def _sub(x, y, n):
            if len(x) > n:
                idx = np.random.choice(len(x), n, replace=False)
                return x[idx], y[idx]
            return x, y

        self.xo, self.yo = _sub(x_oil, y_oil, PREVIEW_N)
        self.xc, self.yc = _sub(x_cort, y_cort, PREVIEW_N)

        self.rot_oil = 0.0
        self.rot_cort = 0.0
        self.flip_oil = False
        self.flip_cort = False
        self.gap = 500.0
        self._saved = False

        self.fig = plt.figure(figsize=(13, 7))
        self.fig.patch.set_facecolor("white")
        self.fig.suptitle(
            "Adjust rotation / flip / gap, then click Save.",
            fontsize=10,
        )

        # Two axes for the hemispheres
        self.ax_oil = self.fig.add_axes([0.03, 0.22, 0.44, 0.72])
        self.ax_crt = self.fig.add_axes([0.53, 0.22, 0.44, 0.72])
        for ax, lbl in [(self.ax_oil, "OIL"), (self.ax_crt, "CORT")]:
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(lbl, fontsize=10, fontweight="bold")
            for sp in ax.spines.values():
                sp.set_visible(False)

        self._sc_oil = self.ax_oil.scatter([], [], s=1.0, c="#BBBBBB", linewidths=0)
        self._sc_crt = self.ax_crt.scatter([], [], s=1.0, c="#BBBBBB", linewidths=0)

        # ── sliders ───────────────────────────────────────────────────────────
        sl_kw = dict(valmin=-180, valmax=180, valinit=0, valstep=0.5)

        ax_ro = self.fig.add_axes([0.07, 0.14, 0.35, 0.03])
        self.sl_rot_oil = Slider(ax_ro, "OIL rot°", **sl_kw)
        self.sl_rot_oil.on_changed(self._update)

        ax_rc = self.fig.add_axes([0.57, 0.14, 0.35, 0.03])
        self.sl_rot_crt = Slider(ax_rc, "CORT rot°", **sl_kw)
        self.sl_rot_crt.on_changed(self._update)

        ax_gap = self.fig.add_axes([0.25, 0.09, 0.5, 0.03])
        x_span = max(self.xo.max() - self.xo.min(), self.xc.max() - self.xc.min())
        self.sl_gap = Slider(
            ax_gap, "Gap (µm)", valmin=0, valmax=x_span,
            valinit=x_span * 0.05, valstep=10,
        )
        self.sl_gap.on_changed(self._update)

        # ── flip toggles ──────────────────────────────────────────────────────
        ax_fo = self.fig.add_axes([0.07, 0.02, 0.15, 0.06])
        self._chk_oil = CheckButtons(ax_fo, ["Flip OIL (mirror X)"], [False])
        self._chk_oil.on_clicked(self._toggle_flip_oil)

        ax_fc = self.fig.add_axes([0.57, 0.02, 0.15, 0.06])
        self._chk_crt = CheckButtons(ax_fc, ["Flip CORT (mirror X)"], [False])
        self._chk_crt.on_clicked(self._toggle_flip_cort)

        # ── save button ───────────────────────────────────────────────────────
        ax_save = self.fig.add_axes([0.40, 0.02, 0.18, 0.06])
        self._btn = Button(ax_save, "Save params")
        self._btn.on_clicked(self._save)

        self._update(None)

    def _toggle_flip_oil(self, _label):
        self.flip_oil = not self.flip_oil
        self._update(None)

    def _toggle_flip_cort(self, _label):
        self.flip_cort = not self.flip_cort
        self._update(None)

    def _transformed(self):
        xo, yo = transform(self.xo.copy(), self.yo.copy(),
                           self.sl_rot_oil.val, self.flip_oil)
        xc, yc = transform(self.xc.copy(), self.yc.copy(),
                           self.sl_rot_crt.val, self.flip_cort)
        return xo, yo, xc, yc

    def _update(self, _val):
        xo, yo, xc, yc = self._transformed()

        self._sc_oil.set_offsets(np.column_stack([xo, yo]))
        self._sc_crt.set_offsets(np.column_stack([xc, yc]))

        for ax, x, y in [(self.ax_oil, xo, yo), (self.ax_crt, xc, yc)]:
            pad = (x.max() - x.min()) * 0.05
            ax.set_xlim(x.min() - pad, x.max() + pad)
            ax.set_ylim(y.min() - pad, y.max() + pad)

        self.fig.canvas.draw_idle()

    def _save(self, _event):
        self._saved = True
        self.final_rot_oil = float(self.sl_rot_oil.val)
        self.final_rot_crt = float(self.sl_rot_crt.val)
        self.final_gap = float(self.sl_gap.val)
        plt.close(self.fig)

    def run(self):
        plt.show(block=True)
        if not self._saved:
            sys.exit("Alignment window closed without saving — exiting.")
        return {
            "rotation_deg": self.final_rot_oil,
            "flip_x": self.flip_oil,
        }, {
            "rotation_deg": self.final_rot_crt,
            "flip_x": self.flip_cort,
        }, self.final_gap


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cort", default=None, help="CORT sample ID (default: from config)")
    parser.add_argument("--oil",  default=None, help="OIL sample ID (default: from config)")
    args = parser.parse_args()

    cfg = load_config()

    # Resolve sample IDs from args or config
    samples = {s["label"]: s["id"] for s in cfg["samples"]}
    cort_id = args.cort or samples.get("CORT", "B03")
    oil_id  = args.oil  or samples.get("OIL",  "B14")

    print(f"Loading OIL  sample {oil_id} …")
    x_oil, y_oil = load_coords(oil_id)
    print(f"  {len(x_oil):,} beads")

    print(f"Loading CORT sample {cort_id} …")
    x_cort, y_cort = load_coords(cort_id)
    print(f"  {len(x_cort):,} beads")

    # ── Phase 1: lasso ────────────────────────────────────────────────────────
    print("\nPhase 1 — OIL hemisphere selection")
    oil_verts = LassoSession(x_oil, y_oil, f"OIL  ({oil_id})").run()
    oil_mask  = apply_lasso(x_oil, y_oil, oil_verts)
    print(f"  Selected {oil_mask.sum():,} / {len(oil_mask):,} beads")

    print("\nPhase 1 — CORT hemisphere selection")
    cort_verts = LassoSession(x_cort, y_cort, f"CORT  ({cort_id})").run()
    cort_mask  = apply_lasso(x_cort, y_cort, cort_verts)
    print(f"  Selected {cort_mask.sum():,} / {len(cort_mask):,} beads")

    # ── Phase 2: alignment ────────────────────────────────────────────────────
    print("\nPhase 2 — alignment")
    align_oil, align_cort, gap = AlignSession(
        x_oil[oil_mask], y_oil[oil_mask],
        x_cort[cort_mask], y_cort[cort_mask],
    ).run()

    # ── Save ──────────────────────────────────────────────────────────────────
    params = {
        "oil": {
            "sample": oil_id,
            "lasso_poly": [[float(v[0]), float(v[1])] for v in oil_verts],
            **align_oil,
        },
        "cort": {
            "sample": cort_id,
            "lasso_poly": [[float(v[0]), float(v[1])] for v in cort_verts],
            **align_cort,
        },
        "alignment": {
            "x_gap_um": gap,
        },
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(params, f, indent=2)

    print(f"\nSaved → {OUTPUT_PATH}")
    print("  OIL :", align_oil)
    print("  CORT:", align_cort)
    print(f"  gap : {gap:.0f} µm")


if __name__ == "__main__":
    main()
