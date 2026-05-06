"""
Figure: spatial expression map for a gene across two brain samples,
with one hemisphere per sample selected and aligned via hemisphere_params.json.

Layout: side-by-side panels (OIL left, CORT right). Each panel shows only
the selected hemisphere, cropped by line-split, rotated, flipped, and scaled
to match the interactive alignment set in select_hemisphere.py.

Layers (bottom to top):
  1. Gray         — all quality-filtered beads in the hemisphere
  2. fg_zero_color — selected cell type, expression = 0
  3. Purple cmap  — selected cell type, expression > 0 (white → #732B8B)

colormap in spatial_expression.yaml accepts either a matplotlib colormap name
or a hex color string (e.g. "#732B8B"), which builds a white→hex linear colormap.

Usage:
    uv run python scripts/figure_code/fig_spatial_expression_v3.py <cell_type> <gene>
    uv run python scripts/figure_code/fig_spatial_expression_v3.py Oligo_NN Sgk1
    uv run python scripts/figure_code/fig_spatial_expression_v3.py hippo Sgk1
"""

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

matplotlib.rcParams["font.family"] = "Liberation Sans"

ROOT        = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"
HEM_PARAMS  = ROOT / "conf" / "hemisphere_params.json"

CELL_TYPE_ALIASES = {
    "hippo": [
        "DG_Glut",
        "CA1_ProS_Glut",
        "CA2_FC_IG_Glut",
        "CA3_Glut",
        "SUB_ProS_Glut",
        "NP_SUB_Glut",
        "CT_SUB_Glut",
        "HPF_CR_Glut",
    ],
    "ca": [
        "CA1_ProS_Glut",
        "CA2_FC_IG_Glut",
        "CA3_Glut",
        "DG_Glut",
    ],
    "all": None,
}


# ── hemisphere transforms (mirrors select_hemisphere.py) ─────────────────────

def _line_split_mask(x, y, pt1, pt2, side):
    x1, y1 = pt1
    x2, y2 = pt2
    cross = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    return (cross * side) > 0


def _rotate(x, y, angle_deg):
    cx, cy = x.mean(), y.mean()
    rad = np.deg2rad(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return cx + c*(x-cx) - s*(y-cy), cy + s*(x-cx) + c*(y-cy)


def _flip_x(x):
    return 2 * x.mean() - x


def apply_hemisphere_transform(df: pd.DataFrame, key: str, params: dict) -> pd.DataFrame:
    """
    Crop df to one hemisphere and apply rotation/flip/scale.
    Returns a new df with x, y replaced by transformed coordinates.
    """
    p  = params[key]
    ls = params["line_split"][key]

    x = df["x"].to_numpy(float)
    y = df["y"].to_numpy(float)

    mask = _line_split_mask(x, y, ls["pt1"], ls["pt2"], ls["side"])
    x, y = x[mask], y[mask]

    x, y = _rotate(x, y, p["rotation_deg"])
    if p.get("flip_x", False):
        x = _flip_x(x)

    scale = p.get("scale", 1.0)
    x = x * scale
    y = y * scale

    out = df[mask].copy()
    out["x"] = x
    out["y"] = y
    return out


def _straighten_cut_edge(df_bg: pd.DataFrame, df_fg: pd.DataFrame, hkey: str):
    """
    Shear-correct a section so its inner cut edge is perfectly vertical.
    Fits a line to the inner-edge beads of df_bg (rightmost for 'oil',
    leftmost for 'cort'), then removes the slope from both df_bg and df_fg.
    Returns (df_bg, df_fg, cut_edge_x) where cut_edge_x is the median inner
    edge x after correction — use this to anchor the gap, not a global percentile.
    """
    x = df_bg["x"].to_numpy(float)
    y = df_bg["y"].to_numpy(float)

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
        return df_bg, df_fg, (x.max() if hkey == "oil" else x.min())

    slope = np.polyfit(y_mid, x_edge, 1)[0]   # dx/dy of the inner edge
    y_ref = float(np.median(y))

    df_bg = df_bg.copy()
    df_fg = df_fg.copy()
    for df in (df_bg, df_fg):
        df["x"] = df["x"].to_numpy(float) - slope * (df["y"].to_numpy(float) - y_ref)

    # After shear, inner-edge beads all land at: mean(x_edge) + slope*(y_ref - mean(y_mid))
    cut_edge_x = float(np.mean(x_edge) + slope * (y_ref - np.mean(y_mid)))
    return df_bg, df_fg, cut_edge_x


# ── plotting helpers (unchanged from v1) ─────────────────────────────────────

def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def make_colormap(base_cmap: str):
    if ":" in base_cmap:
        start_hex, end_hex = base_cmap.split(":", 1)
        start_rgb = np.array(mcolors.to_rgb(start_hex))
        end_rgb   = np.array(mcolors.to_rgb(end_hex))
        t = np.linspace(0.0, 1.0, 256)[:, None]
        colors = start_rgb + t * (end_rgb - start_rgb)
        return mcolors.ListedColormap(colors)
    if base_cmap.startswith("#"):
        end_rgb = mcolors.to_rgb(base_cmap)
        colors = np.array(
            [(1.0, 1.0, 1.0)] +
            [tuple(1.0 - (1.0 - c) * t for c in end_rgb)
             for t in np.linspace(0.05, 1.0, 255)]
        )
        return mcolors.ListedColormap(colors)
    cmap = plt.get_cmap(base_cmap)
    colors = cmap(np.linspace(0.05, 1.0, 256))
    return mcolors.ListedColormap(colors)


def resolve_fg(df_filtered: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    subclass   = cfg.get("cell_subclass")
    cell_class = cfg.get("cell_class")
    if subclass:
        if isinstance(subclass, list):
            return df_filtered[df_filtered["cell_subclass"].isin(subclass)].copy()
        return df_filtered[df_filtered["cell_subclass"] == subclass].copy()
    if cell_class:
        return df_filtered[df_filtered["cell_class"] == cell_class].copy()
    return df_filtered.copy()


def plot_sample(ax, df_bg, df_fg, gene: str, cfg: dict, vmax: float, cmap):
    ax.set_facecolor("white")
    ax.set_aspect("equal")
    expr_col = f"{gene}_expr"

    ax.scatter(
        df_bg["x"], df_bg["y"],
        s=cfg["dot_size_bg"], c=cfg["bg_color"],
        linewidths=0, rasterized=True, zorder=1,
    )

    if df_fg is not None and len(df_fg) > 0:
        zero       = df_fg[df_fg[expr_col] == 0]
        expressing = df_fg[df_fg[expr_col] > 0]

        if len(zero) > 0:
            ax.scatter(
                zero["x"], zero["y"],
                s=cfg["dot_size_fg"], c=cfg["fg_zero_color"],
                linewidths=0, rasterized=True, zorder=2,
            )

        if len(expressing) > 0:
            sc_obj = ax.scatter(
                expressing["x"], expressing["y"],
                s=cfg["dot_size_fg"], c=expressing[expr_col],
                cmap=cmap, vmin=0, vmax=vmax,
                linewidths=0, rasterized=True, zorder=3,
            )
            return sc_obj

    return ax.scatter([], [], c=[], cmap=cmap, vmin=0, vmax=vmax)


def add_scale_bar(ax, scale_um: float, df: pd.DataFrame):
    x_range  = df["x"].max() - df["x"].min()
    y_range  = df["y"].max() - df["y"].min()
    margin_x = x_range * 0.03
    margin_y = y_range * 0.03
    x0 = df["x"].max() - margin_x - scale_um
    x1 = df["x"].max() - margin_x
    y0 = df["y"].min() + margin_y

    ax.plot([x0, x1], [y0, y0], color="black", lw=2, solid_capstyle="butt")
    label = f"{int(scale_um/1000)} mm" if scale_um >= 1000 else f"{int(scale_um)} µm"
    ax.text(
        (x0 + x1) / 2, y0 + y_range * 0.015,
        label, ha="center", va="bottom", fontsize=16, color="black",
    )


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cell_type", help="Cell subclass to highlight (e.g. Oligo_NN, hippo, all)")
    parser.add_argument("gene",      help="Gene symbol (e.g. Sgk1)")
    parser.add_argument("--dot-size", type=float, default=None,
                        help="Override dot size for bg and fg")
    args = parser.parse_args()

    cfg = load_config()
    cfg["cell_subclass"] = None
    cfg["cell_class"]    = None
    if args.dot_size is not None:
        cfg["dot_size_bg"] = args.dot_size
        cfg["dot_size_fg"] = args.dot_size
    cfg["gene"]           = args.gene
    cfg["_cell_type_arg"] = args.cell_type

    gene    = cfg["gene"]
    samples = cfg["samples"]

    with open(HEM_PARAMS) as f:
        hparams = json.load(f)

    # sample id → hemisphere key ("oil" or "cort")
    hem_key_map = {
        hparams["oil"]["sample"]:  "oil",
        hparams["cort"]["sample"]: "cort",
    }

    csv_path = ROOT / cfg["processed_csv"].format(gene=gene)
    if not csv_path.exists():
        print(f"CSV not found for {gene} — running extraction…")
        subprocess.run(
            ["uv", "run", "python",
             str(ROOT / "scripts" / "extract" / "extract_spatial_expression.py"),
             gene],
            check=True, cwd=ROOT,
        )
    else:
        print(f"CSV found for {gene}, skipping extraction.")

    df_all      = pd.read_csv(csv_path)
    df_filtered = df_all[df_all["quality_pass"]].copy()

    cell_type_arg = cfg["_cell_type_arg"]
    if cell_type_arg in CELL_TYPE_ALIASES:
        cfg["cell_subclass"] = CELL_TYPE_ALIASES[cell_type_arg]
        cfg["_alias"] = cell_type_arg
    elif cell_type_arg in df_filtered["cell_subclass"].values:
        cfg["cell_subclass"] = cell_type_arg
    elif cell_type_arg in df_filtered["cell_class"].values:
        cfg["cell_class"] = cell_type_arg
    else:
        available_sub = sorted(df_filtered["cell_subclass"].dropna().unique())
        available_cls = sorted(df_filtered["cell_class"].dropna().unique())
        raise ValueError(
            f"'{cell_type_arg}' not found in cell_subclass or cell_class.\n"
            f"  cell_class values: {available_cls}\n"
            f"  cell_subclass values (first 10): {available_sub[:10]} ..."
        )

    df_fg_all = resolve_fg(df_filtered, cfg)

    expr_col = f"{gene}_expr"
    vmax = df_fg_all[expr_col].quantile(0.99)
    vmax = max(vmax, 0.01)
    cmap = make_colormap(cfg["colormap"])

    fig, ax = plt.subplots(
        1, 1,
        figsize=(cfg["fig_width"], cfg["fig_height"]),
        facecolor="white",
    )

    sc_obj       = None
    df_bg_placed = []   # list of (label, df_bg)

    for sample_info in samples:
        sid   = sample_info["id"]
        label = sample_info["label"]
        hkey  = hem_key_map.get(sid)

        df_bg_raw = df_filtered[df_filtered["sample"] == sid]
        df_fg_raw = df_fg_all[df_fg_all["sample"] == sid]

        if hkey is not None and "line_split" in hparams:
            df_bg = apply_hemisphere_transform(df_bg_raw, hkey, hparams)
            # Derive fg from the already-transformed bg so rotation/flip centers are identical
            df_fg = df_bg[df_bg.index.isin(df_fg_raw.index)].copy()
            df_bg, df_fg, cut_edge_x = _straighten_cut_edge(df_bg, df_fg, hkey)

            if "alignment" in hparams:
                aln      = hparams["alignment"]
                gap      = aln["x_gap_um"]
                yoff     = aln.get("y_offset_um", 0.0)
                x_anchor = gap / 2 if hkey == "cort" else -gap / 2
                y_shift  = yoff    if hkey == "cort" else 0.0

                y_bg = df_bg["y"].to_numpy(float)
                # Anchor using the bin-based inner edge x from _straighten_cut_edge
                x_off = -cut_edge_x + x_anchor
                y_off = -y_bg.mean() + y_shift

                df_bg = df_bg.copy()
                df_bg["x"] += x_off
                df_bg["y"] += y_off
                df_fg = df_fg.copy()
                df_fg["x"] += x_off
                df_fg["y"] += y_off
        else:
            df_bg = df_bg_raw.copy()
            df_fg = df_fg_raw.copy()

        result = plot_sample(ax, df_bg, df_fg, gene, cfg, vmax, cmap)
        if result is not None:
            sc_obj = result

        df_bg_placed.append((label, df_bg))

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    df_bg_combined = pd.concat([df for _, df in df_bg_placed])
    add_scale_bar(ax, cfg["scale_bar_um"], df_bg_combined)

    # OIL / CORT labels — symmetric distance from gap centre (x=0)
    x_range   = df_bg_combined["x"].max() - df_bg_combined["x"].min()
    margin    = x_range * 0.02
    hem_left  = next((lbl, df) for lbl, df in df_bg_placed if df["x"].mean() < 0)
    hem_right = next((lbl, df) for lbl, df in df_bg_placed if df["x"].mean() >= 0)
    x_sym   = max(abs(hem_left[1]["x"].min()), abs(hem_right[1]["x"].max())) + margin
    x_label = x_sym / 4 + 600
    y_label = max(hem_left[1]["y"].max(), hem_right[1]["y"].max())
    ax.text(-x_label, y_label, hem_left[0],
            ha="right", va="top", fontsize=24, fontweight="bold")
    ax.text(+x_label, y_label, hem_right[0],
            ha="left",  va="top", fontsize=24, fontweight="bold")

    # Force gap centre (x=0) to be the figure centre
    x_half = max(abs(df_bg_combined["x"].min()), abs(df_bg_combined["x"].max()), x_sym)
    y_all  = df_bg_combined["y"]
    y_mid  = y_all.mean()
    y_half = max(y_all.max() - y_mid, y_mid - y_all.min())
    ax.set_xlim(-x_half, x_half)
    ax.set_ylim(y_mid - y_half, y_mid + y_half)

    subclass   = cfg.get("cell_subclass")
    cell_class = cfg.get("cell_class")
    alias      = cfg.get("_alias")
    fg_label   = alias or (subclass if isinstance(subclass, str) else None) or cell_class

    if sc_obj is not None:
        cbar_ax = fig.add_axes([0.35, 0.01, 0.3, 0.06])
        cb = fig.colorbar(sc_obj, cax=cbar_ax, orientation="horizontal")
        cb.set_label(f"{gene} expression", fontsize=20)
        cb.ax.tick_params(labelsize=20)
        if fg_label:
            cbar_ax.set_title(fg_label, fontsize=20, pad=8, fontweight="bold")

    plt.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.07, wspace=0)

    cell_type_slug = (
        alias or (subclass if isinstance(subclass, str) else None) or cell_class or "all"
    ).replace(" ", "_")
    out_path = ROOT / "figures" / "genes_half_half" / f"fig_spatial_{gene}_{cell_type_slug}_v3.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    print(f"Saved → {out_path}")


if __name__ == "__main__":
    main()
