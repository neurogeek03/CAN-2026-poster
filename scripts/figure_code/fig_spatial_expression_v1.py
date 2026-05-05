"""
Figure: spatial expression map for a gene across two brain samples.

Layout: side-by-side panels (left = first sample in config, right = second).

Layers (bottom to top):
  1. Gray       — all quality-filtered beads (brain shape)
  2. fg_zero_color — selected cell type, expression = 0
  3. Reds cmap  — selected cell type, expression > 0

Usage:
    uv run python scripts/figures/fig_spatial_expression_v1.py <cell_type> <gene>
    uv run python scripts/figures/fig_spatial_expression_v1.py Oligo_NN Sgk1
    uv run python scripts/figures/fig_spatial_expression_v1.py hippo Sgk1
"""

import argparse
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "conf" / "spatial_expression.yaml"

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


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def make_colormap(base_cmap: str):
    cmap = plt.get_cmap(base_cmap)
    colors = cmap(np.linspace(0.05, 1.0, 256))
    return mcolors.ListedColormap(colors)


def resolve_fg(df_filtered: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Return the foreground subset according to cell_subclass > cell_class > all."""
    subclass = cfg.get("cell_subclass")
    cell_class = cfg.get("cell_class")
    if subclass:
        if isinstance(subclass, list):
            return df_filtered[df_filtered["cell_subclass"].isin(subclass)].copy()
        return df_filtered[df_filtered["cell_subclass"] == subclass].copy()
    if cell_class:
        return df_filtered[df_filtered["cell_class"] == cell_class].copy()
    return df_filtered.copy()


def plot_sample(ax, df_bg, df_fg, gene: str, cfg: dict, vmax: float, cmap):
    """Draw one panel. Returns a scatter mappable for the colorbar."""
    ax.set_facecolor("white")
    ax.set_aspect("equal")
    expr_col = f"{gene}_expr"

    # Layer 1 — all quality-filtered beads (gray)
    ax.scatter(
        df_bg["x"], df_bg["y"],
        s=cfg["dot_size_bg"],
        c=cfg["bg_color"],
        linewidths=0,
        rasterized=True,
        zorder=1,
    )

    if df_fg is not None and len(df_fg) > 0:
        zero = df_fg[df_fg[expr_col] == 0]
        expressing = df_fg[df_fg[expr_col] > 0]

        # Layer 2 — selected type, not expressing
        if len(zero) > 0:
            ax.scatter(
                zero["x"], zero["y"],
                s=cfg["dot_size_fg"],
                c=cfg["fg_zero_color"],
                linewidths=0,
                rasterized=True,
                zorder=2,
            )

        # Layer 3 — selected type, expressing (colored by level)
        if len(expressing) > 0:
            sc_obj = ax.scatter(
                expressing["x"], expressing["y"],
                s=cfg["dot_size_fg"],
                c=expressing[expr_col],
                cmap=cmap,
                vmin=0,
                vmax=vmax,
                linewidths=0,
                rasterized=True,
                zorder=3,
            )
            return sc_obj

    return ax.scatter([], [], c=[], cmap=cmap, vmin=0, vmax=vmax)


def add_scale_bar(ax, scale_um: float, df: pd.DataFrame):
    x_range = df["x"].max() - df["x"].min()
    y_range = df["y"].max() - df["y"].min()
    margin_x = x_range * 0.03
    margin_y = y_range * 0.03
    x0 = df["x"].max() - margin_x - scale_um
    x1 = df["x"].max() - margin_x
    y0 = df["y"].min() + margin_y

    ax.plot([x0, x1], [y0, y0], color="black", lw=2, solid_capstyle="butt")
    label = f"{int(scale_um/1000)} mm" if scale_um >= 1000 else f"{int(scale_um)} µm"
    ax.text(
        (x0 + x1) / 2, y0 + y_range * 0.015,
        label, ha="center", va="bottom",
        fontsize=8, color="black",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cell_type", help="Cell subclass to highlight (e.g. Oligo_NN, hippo, all)")
    parser.add_argument("gene", help="Gene symbol (e.g. Sgk1)")
    parser.add_argument("--dot-size", type=float, default=None, help="Override dot size for both bg and fg (default: from config)")
    args = parser.parse_args()

    cfg = load_config()
    cfg["cell_subclass"] = None
    cfg["cell_class"] = None
    if args.dot_size is not None:
        cfg["dot_size_bg"] = args.dot_size
        cfg["dot_size_fg"] = args.dot_size
    cfg["gene"] = args.gene
    cfg["_cell_type_arg"] = args.cell_type  # resolved after CSV is loaded

    gene = cfg["gene"]
    samples = cfg["samples"]

    csv_path = ROOT / cfg["processed_csv"].format(gene=gene)
    if not csv_path.exists():
        print(f"CSV not found for {gene} — running extraction…")
        subprocess.run(
            ["uv", "run", "python",
             str(ROOT / "scripts" / "extract" / "extract_spatial_expression.py"),
             gene],
            check=True,
            cwd=ROOT,
        )
    else:
        print(f"CSV found for {gene}, skipping extraction.")

    df_all = pd.read_csv(csv_path)
    df_filtered = df_all[df_all["quality_pass"]].copy()

    cell_type_arg = cfg["_cell_type_arg"]
    if cell_type_arg in CELL_TYPE_ALIASES:
        cfg["cell_subclass"] = CELL_TYPE_ALIASES[cell_type_arg]  # None = show all
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

    n_panels = len(samples)
    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(cfg["fig_width"], cfg["fig_height"]),
        facecolor="white",
    )
    if n_panels == 1:
        axes = [axes]

    sc_obj = None
    for ax, sample_info in zip(axes, samples):
        sid = sample_info["id"]
        label = sample_info["label"]

        df_bg = df_filtered[df_filtered["sample"] == sid]
        df_fg = df_fg_all[df_fg_all["sample"] == sid]

        result = plot_sample(ax, df_bg, df_fg, gene, cfg, vmax, cmap)
        if result is not None:
            sc_obj = result

        ax.set_title(label, fontsize=11, fontweight="bold", pad=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        add_scale_bar(ax, cfg["scale_bar_um"], df_bg if len(df_bg) > 0 else df_all)

    # ── Horizontal colorbar ───────────────────────────────────────────────────
    if sc_obj is not None:
        cbar_ax = fig.add_axes([0.2, 0.06, 0.6, 0.025])
        cb = fig.colorbar(sc_obj, cax=cbar_ax, orientation="horizontal")
        cb.set_label(f"{gene} expression (SCT)", fontsize=9)
        cb.ax.tick_params(labelsize=8)

    # ── Legend for zero-expression cell type color ────────────────────────────
    subclass = cfg.get("cell_subclass")
    cell_class = cfg.get("cell_class")
    alias = cfg.get("_alias")
    fg_label = alias or (subclass if isinstance(subclass, str) else None) or cell_class
    if fg_label:
        patch = mpatches.Patch(color=cfg["fg_zero_color"], label=f"{fg_label} (not expressed)")
        fig.legend(
            handles=[patch],
            loc="lower right",
            bbox_to_anchor=(0.98, 0.02),
            frameon=False,
            fontsize=8,
        )

    # ── Title ─────────────────────────────────────────────────────────────────
    fg_label = fg_label or "all cells"
    fig.suptitle(
        f"{gene}  ·  {fg_label}   |   Slide-seq",
        fontsize=12, y=0.97, fontweight="bold",
    )

    plt.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.12, wspace=0.04)

    cell_type_slug = (alias or (subclass if isinstance(subclass, str) else None) or cell_class or "all").replace(" ", "_")
    out_path = ROOT / cfg["output_svg"].format(gene=gene, cell_type=cell_type_slug)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=cfg["dpi"], bbox_inches="tight", facecolor="white")
    print(f"Saved → {out_path}")


if __name__ == "__main__":
    main()
