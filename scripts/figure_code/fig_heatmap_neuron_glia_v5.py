"""
Figure v5 — Split layout:
  Panel A (themes 1–4): all 19 cell types  → figures/fig_heatmap_themes1to4_v5.svg
  Panel B (themes 5–8): neurons only       → figures/fig_heatmap_themes5to8_neurons_v5.svg

Input:  data/processed/heatmap_neuron_glia.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

IN_CSV = Path("data/processed/heatmap_neuron_glia.csv")

# ── Shared config ──────────────────────────────────────────────────────────────
FDR_THRESH = 0.05
VMAX       = 4.0
CELL_W     = 0.42
CELL_H     = 0.38
CMAP       = "RdBu_r"
NOT_TESTED = "#d0d0d0"

THEME_BAND_W   = 0.42
THEME_BAND_GAP = 0.08
LEFT_MARGIN    = 0.9
TOP_MARGIN     = 0.9
BOT_MARGIN     = 0.9
RIGHT_MARGIN   = THEME_BAND_GAP + THEME_BAND_W + 0.15

ZONE_COLORS = {
    "Hippocampal\nneurons":    "#1a6b72",
    "Cortical\nglutamatergic": "#5a3472",
    "Interneurons":            "#7b3f1e",
    "Non-neuronal":            "#3d3d6b",
}

THEME_COLORS = {
    "1. GR desensitization":   "#4e79a7",
    "2. Neurotrophic collapse": "#59a14f",
    "3. Circadian disruption":  "#9467bd",
    "4. Neuroinflammation":     "#e15759",
    "5. E/I balance":           "#f28e2b",
    "7. cAMP/CREB plasticity":  "#439894",
    "8. Wnt/synaptic adhesion": "#b07a35",
}

# ── Panel definitions ──────────────────────────────────────────────────────────
PANELS = [
    {
        "out":    Path("figures/fig_heatmap_themes1to4_v5.svg"),
        "themes": [
            "1. GR desensitization",
            "2. Neurotrophic collapse",
            "3. Circadian disruption",
            "4. Neuroinflammation",
        ],
        "col_zones": {
            "Hippocampal\nneurons":       ["CA1-ProS", "CA3", "DG"],
            "Cortical\nglutamatergic":    ["L2/3 IT", "L4/5 IT", "L5 ET", "L6 IT", "L6 CT"],
            "Interneurons":               ["Pvalb", "Sst"],
            "Non-neuronal":               ["Astro-NT", "Astro-TE", "OPC", "Oligo",
                                           "Microglia", "Endo", "Peri", "VLMC", "CHOR"],
        },
    },
    {
        "out":    Path("figures/fig_heatmap_themes5to8_neurons_v5.svg"),
        "themes": [
            "5. E/I balance",
            "7. cAMP/CREB plasticity",
            "8. Wnt/synaptic adhesion",
        ],
        "col_zones": {
            "Hippocampal\nneurons":       ["CA1-ProS", "CA3", "DG"],
            "Cortical\nglutamatergic":    ["L2/3 IT", "L4/5 IT", "L5 ET", "L6 IT", "L6 CT"],
            "Interneurons":               ["Pvalb", "Sst"],
        },
    },
]

# ── Drawing function ───────────────────────────────────────────────────────────
def draw_panel(df_all, theme_order, col_zones, out_path):
    col_order = [ct for zone_cols in col_zones.values() for ct in zone_cols]

    gene_order = []
    theme_map  = {}
    for theme in theme_order:
        genes = list(df_all[df_all["theme"] == theme].drop_duplicates("gene")["gene"])
        gene_order.extend(genes)
        for g in genes:
            theme_map[g] = theme

    df = df_all[df_all["gene"].isin(gene_order)]

    logfc_wide = (
        df.pivot(index="gene", columns="cell_type", values="logFC")
        .reindex(index=gene_order, columns=col_order)
    )
    fdr_wide = (
        df.pivot(index="gene", columns="cell_type", values="FDR")
        .reindex(index=gene_order, columns=col_order)
    )

    n_rows = len(gene_order)
    n_cols = len(col_order)

    hmap_w = n_cols * CELL_W
    hmap_h = n_rows * CELL_H
    fig_w  = LEFT_MARGIN + hmap_w + RIGHT_MARGIN
    fig_h  = TOP_MARGIN  + hmap_h + BOT_MARGIN

    fig = plt.figure(figsize=(fig_w, fig_h))

    ax_l = LEFT_MARGIN / fig_w
    ax_b = BOT_MARGIN  / fig_h
    ax_w = hmap_w / fig_w
    ax_h = hmap_h / fig_h

    ax = fig.add_axes([ax_l, ax_b, ax_w, ax_h])

    # heatmap
    mat      = logfc_wide.values.astype(float)
    tested   = ~np.isnan(mat)
    mat_disp = np.where(tested, mat, 0.0)

    norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)
    im   = ax.imshow(mat_disp, aspect="auto", cmap=CMAP, norm=norm,
                     interpolation="nearest")

    for r in range(n_rows):
        for c in range(n_cols):
            if not tested[r, c]:
                ax.add_patch(mpatches.Rectangle(
                    (c - 0.5, r - 0.5), 1, 1, color=NOT_TESTED, zorder=2))

    sig_r, sig_c = np.where((fdr_wide.values < FDR_THRESH) & tested)
    ax.scatter(sig_c, sig_r, s=18, c="black", zorder=3, linewidths=0)

    # theme row separators
    theme_boundaries = []
    prev_theme = None
    for i, gene in enumerate(gene_order):
        t = theme_map[gene]
        if t != prev_theme:
            theme_boundaries.append(i)
            prev_theme = t

    for boundary in theme_boundaries[1:]:
        ax.axhline(boundary - 0.5, color="white", linewidth=2.5, zorder=4)

    # column zone separators
    zone_starts = {}
    col_idx = 0
    for zone, cols in col_zones.items():
        zone_starts[zone] = col_idx
        col_idx += len(cols)

    for zone in list(col_zones.keys())[1:]:
        ax.axvline(zone_starts[zone] - 0.5, color="white", linewidth=3, zorder=4)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_order, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(gene_order, fontstyle="italic", fontsize=8)
    ax.tick_params(length=0)

    # zone header bands
    band_h_frac = THEME_BAND_W / fig_h
    for zone, cols in col_zones.items():
        z_start = zone_starts[zone]
        z_ncols = len(cols)
        x0_frac = ax_l + (z_start / n_cols) * ax_w
        w_frac  = (z_ncols / n_cols) * ax_w
        band_ax = fig.add_axes([x0_frac, ax_b + ax_h + 0.004, w_frac, band_h_frac])
        band_ax.set_xlim(0, 1); band_ax.set_ylim(0, 1)
        band_ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, color=ZONE_COLORS[zone]))
        band_ax.text(0.5, 0.5, zone, ha="center", va="center",
                     fontsize=7, color="white", fontweight="bold",
                     transform=band_ax.transAxes)
        band_ax.axis("off")

    # theme annotation band (right)
    theme_ax_x = ax_l + ax_w + THEME_BAND_GAP / fig_w
    theme_ax_w = THEME_BAND_W / fig_w
    theme_ax   = fig.add_axes([theme_ax_x, ax_b, theme_ax_w, ax_h])
    theme_ax.set_xlim(0, 1)
    theme_ax.set_ylim(-0.5, n_rows - 0.5)
    theme_ax.invert_yaxis()
    theme_ax.axis("off")

    for ti, theme in enumerate(theme_order):
        start = theme_boundaries[ti]
        end   = theme_boundaries[ti + 1] if ti + 1 < len(theme_boundaries) else n_rows
        mid   = (start + end - 1) / 2
        theme_ax.add_patch(mpatches.Rectangle(
            (0, start - 0.5), 1, end - start,
            color=THEME_COLORS[theme], zorder=2))
        theme_ax.text(0.5, mid, theme, ha="center", va="center",
                      fontsize=6.5, color="white", fontweight="bold", rotation=90)

    # colorbar
    cbar_y  = 0.12 / fig_h
    cbar_ax = fig.add_axes([ax_l, cbar_y, ax_w * 0.28, 0.025])
    cb = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cb.set_label(f"log₂ fold change (CORT vs OIL), capped ±{VMAX}", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    fig.text(
        ax_l + ax_w * 0.31, 0.38 / fig_h,
        "• FDR < 0.05    grey = not tested",
        fontsize=6.5, va="center", color="#333333",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    print(f"Saved {out_path}")
    plt.close()


# ── Run ────────────────────────────────────────────────────────────────────────
df_all = pd.read_csv(IN_CSV)
for panel in PANELS:
    draw_panel(df_all, panel["themes"], panel["col_zones"], panel["out"])
