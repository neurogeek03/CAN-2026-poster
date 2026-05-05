"""
Interactive hemisphere selector + alignment — run once locally.
Saves all parameters to conf/hemisphere_params.json.

Phase 1: For each sample (OIL then CORT), select one hemisphere:
  - Default (lasso): draw a freehand polygon around the region to keep.
  - --line-split: click two points to draw a straight dividing line;
    click "Flip sides" to switch, then "Confirm".

Phase 2: Both hemispheres shown together in one view. Adjust rotation,
  scale, gap (negative = overlap), and vertical offset. Click "Save".

Usage:
    uv run python scripts/select_hemisphere.py
    uv run python scripts/select_hemisphere.py --line-split
    uv run python scripts/select_hemisphere.py --cort B03 --oil B14

Requirements: active display — run locally, not on HPC.
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
        f"Non-interactive backend ({_backend}).\n"
        "Run this script locally, not on a headless server."
    )

import numpy as np
import pandas as pd
import yaml
from matplotlib.path import Path as MplPath
from matplotlib.widgets import Button, CheckButtons, Slider, LassoSelector

ROOT        = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"
COORDS_CSV  = ROOT / "data" / "processed" / "hemisphere_coords.csv"
OUTPUT_PATH = ROOT / "conf" / "hemisphere_params.json"

PREVIEW_N = 6000   # points used in the live alignment preview


# ── helpers ──────────────────────────────────────────────────────────────────

def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_coords(sample_id: str):
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
    path = MplPath(poly_verts)
    return path.contains_points(np.column_stack([x, y]))


def apply_line_split(x, y, pt1, pt2, side):
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


def subsample(x, y, n=PREVIEW_N):
    if len(x) > n:
        idx = np.random.choice(len(x), n, replace=False)
        return x[idx], y[idx]
    return x.copy(), y.copy()


# ── Phase 1a: line-split selection ───────────────────────────────────────────

class LineSplitSession:
    """Click two points to draw a straight dividing line; keep one side."""

    def __init__(self, x, y, title):
        self.x, self.y = x, y
        self._pts  = []
        self._side = 1
        self._line = None
        self._done = False

        if len(x) > 40_000:
            idx = np.random.choice(len(x), 40_000, replace=False)
            self._xd, self._yd = x[idx], y[idx]
        else:
            self._xd, self._yd = x, y

        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.fig.patch.set_facecolor("white")
        self.fig.suptitle(
            f"{title}\nClick two points to define the dividing line, then Confirm.",
            fontsize=10,
        )
        self._sc_all = self.ax.scatter(
            self._xd, self._yd, s=0.8, c="#BBBBBB", linewidths=0, rasterized=True
        )
        self._sc_sel = self.ax.scatter(
            [], [], s=1.2, c="steelblue", linewidths=0, rasterized=True, zorder=3
        )
        self._pt_markers, = self.ax.plot([], [], "rx", ms=10, mew=2, zorder=5)
        self.ax.set_aspect("equal")
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)

        ax_flip = self.fig.add_axes([0.15, 0.01, 0.25, 0.05])
        self._btn_flip = Button(ax_flip, "Flip sides")
        self._btn_flip.on_clicked(self._flip)

        ax_btn = self.fig.add_axes([0.55, 0.01, 0.25, 0.05])
        self._btn_confirm = Button(ax_btn, "Confirm selection")
        self._btn_confirm.on_clicked(self._confirm)

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

    def _side_mask(self):
        if len(self._pts) < 2:
            return np.zeros(len(self._xd), dtype=bool)
        (x1, y1), (x2, y2) = self._pts
        cross = (x2-x1)*(self._yd-y1) - (y2-y1)*(self._xd-x1)
        return (cross * self._side) > 0

    def _draw_line(self):
        if self._line is not None:
            self._line.remove()
            self._line = None
        if len(self._pts) < 2:
            return
        (x1, y1), (x2, y2) = self._pts
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        dx, dy = x2-x1, y2-y1
        if abs(dx) < 1e-9 and abs(dy) < 1e-9:
            return
        ts = []
        if abs(dx) > 1e-9:
            ts += [(xlim[0]-x1)/dx, (xlim[1]-x1)/dx]
        if abs(dy) > 1e-9:
            ts += [(ylim[0]-y1)/dy, (ylim[1]-y1)/dy]
        t0, t1 = min(ts), max(ts)
        self._line, = self.ax.plot(
            [x1+t0*dx, x1+t1*dx], [y1+t0*dy, y1+t1*dy],
            "r-", lw=1.5, zorder=4
        )

    def _refresh(self):
        mask = self._side_mask()
        self._sc_sel.set_offsets(
            np.column_stack([self._xd[mask], self._yd[mask]]) if mask.any()
            else np.empty((0, 2))
        )
        self._pt_markers.set_data([p[0] for p in self._pts], [p[1] for p in self._pts])
        self._draw_line()
        self.fig.canvas.draw_idle()

    def _on_click(self, event):
        if event.inaxes is not self.ax or event.button != 1:
            return
        if len(self._pts) >= 2:
            self._pts = []
        self._pts.append((event.xdata, event.ydata))
        self._refresh()

    def _flip(self, _e):
        self._side *= -1
        self._refresh()

    def _confirm(self, _e):
        if len(self._pts) < 2:
            print("  Click two points first.")
            return
        self._done = True
        plt.close(self.fig)

    def run(self):
        plt.show(block=True)
        if not self._done:
            sys.exit("Window closed without confirming.")
        return self._pts[0], self._pts[1], self._side


# ── Phase 1b: lasso selection ─────────────────────────────────────────────────

class LassoSession:
    def __init__(self, x, y, title):
        self.x, self.y = x, y
        self.verts = None
        self._done = False

        if len(x) > 40_000:
            idx = np.random.choice(len(x), 40_000, replace=False)
            self._xd, self._yd = x[idx], y[idx]
        else:
            self._xd, self._yd = x, y

        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.fig.patch.set_facecolor("white")
        self.fig.suptitle(
            f"{title}\nDraw a lasso around the hemisphere to keep, then Confirm.",
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
        self._btn_confirm = Button(ax_btn, "Confirm selection")
        self._btn_confirm.on_clicked(self._confirm)

        self._lasso = LassoSelector(self.ax, self._on_select, useblit=True)

    def _on_select(self, verts):
        self.verts = verts
        mask = apply_lasso(self._xd, self._yd, verts)
        self._sc_sel.set_offsets(
            np.column_stack([self._xd[mask], self._yd[mask]])
        )
        self.fig.canvas.draw_idle()

    def _confirm(self, _e):
        if self.verts is None:
            print("  Draw a lasso first.")
            return
        self._done = True
        plt.close(self.fig)

    def run(self):
        plt.show(block=True)
        if not self._done:
            sys.exit("Window closed without confirming.")
        return self.verts


# ── Phase 2: combined alignment view ─────────────────────────────────────────

class AlignSession:
    """Both hemispheres in one view. Sliders: rotation, scale, gap, y-offset."""

    def __init__(self, x_oil, y_oil, x_cort, y_cort):
        self.xo_full, self.yo_full = x_oil,  y_oil
        self.xc_full, self.yc_full = x_cort, y_cort
        self.xo, self.yo = subsample(x_oil,  y_oil)
        self.xc, self.yc = subsample(x_cort, y_cort)
        self._saved = False

        span = max(
            x_oil.max()  - x_oil.min(),
            x_cort.max() - x_cort.min(),
        )

        self.fig = plt.figure(figsize=(13, 8))
        self.fig.patch.set_facecolor("white")
        self.fig.suptitle(
            "Adjust alignment, then click Save", fontsize=11
        )
        self.fig.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.42)

        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect("equal")
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for sp in self.ax.spines.values():
            sp.set_visible(False)
        self.ax.set_facecolor("white")

        self._sc_o = self.ax.scatter([], [], s=1.0, c="#4393C3",
                                     linewidths=0, rasterized=True, label="OIL")
        self._sc_c = self.ax.scatter([], [], s=1.0, c="#D6604D",
                                     linewidths=0, rasterized=True, label="CORT")
        self.ax.legend(loc="upper right", markerscale=8, frameon=False, fontsize=9)

        # ── sliders ───────────────────────────────────────────────────────────
        def sl(left, bottom, label, vmin, vmax, vinit, step=None):
            a = self.fig.add_axes([left, bottom, 0.38, 0.03])
            kw = dict(valmin=vmin, valmax=vmax, valinit=vinit)
            if step is not None:
                kw["valstep"] = step
            s = Slider(a, label, **kw)
            s.on_changed(self._update)
            return s

        self.sl_rot_o  = sl(0.06, 0.34, "OIL rot°",    -180, 180,       0,    0.5)
        self.sl_rot_c  = sl(0.56, 0.34, "CORT rot°",   -180, 180,       0,    0.5)
        self.sl_scl_o  = sl(0.06, 0.28, "OIL scale",    0.3,   2.0,     1.0,  0.01)
        self.sl_scl_c  = sl(0.56, 0.28, "CORT scale",   0.3,   2.0,     1.0,  0.01)

        ax_gap = self.fig.add_axes([0.06, 0.22, 0.88, 0.03])
        self.sl_gap = Slider(ax_gap, "Gap (µm)",
                             valmin=-span * 0.8, valmax=span * 0.5,
                             valinit=0, valstep=10)
        self.sl_gap.on_changed(self._update)

        ax_yoff = self.fig.add_axes([0.06, 0.15, 0.88, 0.03])
        self.sl_yoff = Slider(ax_yoff, "Y offset (µm)",
                              valmin=-span * 0.5, valmax=span * 0.5,
                              valinit=0, valstep=10)
        self.sl_yoff.on_changed(self._update)

        # ── flip checkboxes ───────────────────────────────────────────────────
        self._flip_o = False
        self._flip_c = False

        ax_fo = self.fig.add_axes([0.06, 0.07, 0.22, 0.06])
        self._chk_o = CheckButtons(ax_fo, ["Flip OIL (mirror X)"], [False])
        self._chk_o.on_clicked(lambda _: self._toggle("o"))

        ax_fc = self.fig.add_axes([0.56, 0.07, 0.22, 0.06])
        self._chk_c = CheckButtons(ax_fc, ["Flip CORT (mirror X)"], [False])
        self._chk_c.on_clicked(lambda _: self._toggle("c"))

        ax_save = self.fig.add_axes([0.40, 0.01, 0.20, 0.06])
        self._btn_save = Button(ax_save, "Save params")
        self._btn_save.on_clicked(self._save)

        self._update()

    def _toggle(self, which):
        if which == "o":
            self._flip_o = not self._flip_o
        else:
            self._flip_c = not self._flip_c
        self._update()

    def _place(self, x, y, rot, do_flip, scale, x_anchor, y_shift):
        """Rotate → flip → scale → translate so x_anchor is the edge toward centre."""
        x, y = rotate(x, y, rot)
        if do_flip:
            x = flip_x(x)
        x = x * scale
        y = y * scale
        # Shift so the inner edge sits at x_anchor
        if x_anchor >= 0:   # CORT on the right: left edge → x_anchor
            x = x - x.min() + x_anchor
        else:               # OIL on the left: right edge → x_anchor
            x = x - x.max() + x_anchor
        y = y - y.mean() + y_shift
        return x, y

    def _update(self, _v=None):
        gap   = self.sl_gap.val
        yoff  = self.sl_yoff.val

        xo, yo = self._place(
            self.xo, self.yo,
            self.sl_rot_o.val, self._flip_o, self.sl_scl_o.val,
            x_anchor=-gap/2, y_shift=0,
        )
        xc, yc = self._place(
            self.xc, self.yc,
            self.sl_rot_c.val, self._flip_c, self.sl_scl_c.val,
            x_anchor=+gap/2, y_shift=yoff,
        )

        self._sc_o.set_offsets(np.column_stack([xo, yo]))
        self._sc_c.set_offsets(np.column_stack([xc, yc]))

        all_x = np.concatenate([xo, xc])
        all_y = np.concatenate([yo, yc])
        px = (all_x.max() - all_x.min()) * 0.04
        py = (all_y.max() - all_y.min()) * 0.04
        self.ax.set_xlim(all_x.min()-px, all_x.max()+px)
        self.ax.set_ylim(all_y.min()-py, all_y.max()+py)
        self.fig.canvas.draw_idle()

    def _save(self, _e):
        self._saved = True
        self.result = {
            "oil":  {
                "rotation_deg": float(self.sl_rot_o.val),
                "flip_x":       self._flip_o,
                "scale":        float(self.sl_scl_o.val),
            },
            "cort": {
                "rotation_deg": float(self.sl_rot_c.val),
                "flip_x":       self._flip_c,
                "scale":        float(self.sl_scl_c.val),
            },
            "alignment": {
                "x_gap_um":    float(self.sl_gap.val),
                "y_offset_um": float(self.sl_yoff.val),
            },
        }
        plt.close(self.fig)

    def run(self):
        plt.show(block=True)
        if not self._saved:
            sys.exit("Alignment window closed without saving.")
        return self.result


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cort",       default=None, help="CORT sample ID")
    parser.add_argument("--oil",        default=None, help="OIL sample ID")
    parser.add_argument("--line-split", action="store_true",
                        help="Use straight-line splitting instead of lasso")
    args = parser.parse_args()

    cfg     = load_config()
    samples = {s["label"]: s["id"] for s in cfg["samples"]}
    cort_id = args.cort or samples.get("CORT", "B03")
    oil_id  = args.oil  or samples.get("OIL",  "B14")

    print(f"Loading OIL  {oil_id} …")
    x_oil,  y_oil  = load_coords(oil_id)
    print(f"  {len(x_oil):,} beads")

    print(f"Loading CORT {cort_id} …")
    x_cort, y_cort = load_coords(cort_id)
    print(f"  {len(x_cort):,} beads")

    # ── Phase 1 ───────────────────────────────────────────────────────────────
    split_meta = {}

    if args.line_split:
        print("\nPhase 1 — OIL (line split)")
        pt1o, pt2o, so = LineSplitSession(x_oil,  y_oil,  f"OIL  ({oil_id})").run()
        oil_mask = apply_line_split(x_oil, y_oil, pt1o, pt2o, so)
        print(f"  {oil_mask.sum():,} / {len(oil_mask):,} beads kept")

        print("\nPhase 1 — CORT (line split)")
        pt1c, pt2c, sc = LineSplitSession(x_cort, y_cort, f"CORT  ({cort_id})").run()
        cort_mask = apply_line_split(x_cort, y_cort, pt1c, pt2c, sc)
        print(f"  {cort_mask.sum():,} / {len(cort_mask):,} beads kept")

        split_meta = {
            "oil":  {"pt1": list(pt1o), "pt2": list(pt2o), "side": so},
            "cort": {"pt1": list(pt1c), "pt2": list(pt2c), "side": sc},
        }
        oil_verts  = [list(pt1o), list(pt2o)]
        cort_verts = [list(pt1c), list(pt2c)]
    else:
        print("\nPhase 1 — OIL (lasso)")
        oil_verts  = LassoSession(x_oil,  y_oil,  f"OIL  ({oil_id})").run()
        oil_mask   = apply_lasso(x_oil, y_oil, oil_verts)
        print(f"  {oil_mask.sum():,} / {len(oil_mask):,} beads kept")

        print("\nPhase 1 — CORT (lasso)")
        cort_verts = LassoSession(x_cort, y_cort, f"CORT  ({cort_id})").run()
        cort_mask  = apply_lasso(x_cort, y_cort, cort_verts)
        print(f"  {cort_mask.sum():,} / {len(cort_mask):,} beads kept")

    # ── Phase 2 ───────────────────────────────────────────────────────────────
    print("\nPhase 2 — alignment")
    result = AlignSession(
        x_oil[oil_mask],   y_oil[oil_mask],
        x_cort[cort_mask], y_cort[cort_mask],
    ).run()

    # ── Save ──────────────────────────────────────────────────────────────────
    params = {
        "oil": {
            "sample":     oil_id,
            "lasso_poly": [[float(v[0]), float(v[1])] for v in oil_verts],
            **result["oil"],
        },
        "cort": {
            "sample":     cort_id,
            "lasso_poly": [[float(v[0]), float(v[1])] for v in cort_verts],
            **result["cort"],
        },
        "alignment": result["alignment"],
    }
    if split_meta:
        params["line_split"] = split_meta

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(params, f, indent=2)

    print(f"\nSaved → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
