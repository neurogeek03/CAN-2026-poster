"""
Bipartite network graph of neuron ↔ glia CCC, split into 4 theme panels.
2×2 grid. Each panel sits inside a soft-rounded #732B8B container.
Edge labels = LR pair, UMAP label style (white bold on colored bg).
Nodes drawn with scatter (display-space circles, always circular).
"""

import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages

PANEL_COLOR = "#732B8B"

# ── Data ──────────────────────────────────────────────────────────────────────

df = pd.read_csv("data/LIANA+/lr_res.csv", index_col=0)
sig = df[df.interaction_padj < 0.05].copy()

NEURONS = {"016 CA1-ProS Glut", "017 CA3 Glut", "025 CA2-FC-IG Glut", "037 DG Glut",
           "052 Pvalb Gaba"}
GLIA    = {"318 Astro-NT NN", "319 Astro-TE NN", "326 OPC NN", "327 Oligo NN",
           "334 Microglia NN", "325 CHOR NN", "330 VLMC NN", "331 Peri NN", "333 Endo NN"}

sig = sig[
    (sig.source.isin(NEURONS) & sig.target.isin(GLIA)) |
    (sig.source.isin(GLIA)    & sig.target.isin(NEURONS))
].copy()

THEMES = {
    "Axon guidance": [
        "Sema6a^Plxna2", "Sema6d^Kdr_Plxna1", "Sema3d^Nrp2_Plxna2", "Sema3c^Nrp2_Plxna2",
        "Nell2^Robo3", "Ncam1^Robo1", "Ncam1^Robo3",
        "Slit2^Robo1", "Slit2^Robo2", "Slit3^Robo1", "Efna5^Ephb1", "Slit2^Dcc",
    ],
    "Cell adhesion / CAMs": [
        "Cntn4^Ptprg", "Cntn1^Ptprz1", "Ncam1^Ptprz1", "L1cam^Ptprz1", "L1cam^Cntn1",
        "Nfasc^Cntn1_Cntnap1", "Alcam^Chl1", "Fstl5^Chl1",
        "Ntng1^Lrrc4c", "Ntng2^Lrrc4c",
    ],
    "Neurotrophic / growth factors": [
        "Bdnf^Ntrk2", "Bdnf^Sort1",
        "Ptn^Ncl", "Ptn^Ptprz1", "Ptn^Itgav_Itgb3", "Ptn^Alk", "Ptn^Cdh10",
        "Mdk^Sorl1", "Mdk^Ncl", "Mdk^Alk", "Mdk^Ptprz1", "Fstl1^Dip2a",
    ],
    "ECM–Integrin remodeling": [
        "Fn1^Itgav_Itgb3", "Fn1^Itgav_Itgb6", "Fn1^Itga9", "Fn1^Itgav_Itgb8", "Fn1^Nt5e",
        "Tnr^Itgav_Itgb3", "Tnr^Itgav_Itgb6",
        "Col4a3^Itgav_Itgb8", "Col6a1^Itgav_Itgb8", "Col18a1^Kdr",
        "Lama1^Nt5e", "Lama4^Itgav_Itgb8", "Lamb1^Itgav_Itgb8",
        "Vwf^Itgav_Itgb3", "Vwf^Itga9", "Gpc3^Igf1r", "Gpc3^Unc5c", "Gpc3^Unc5d",
    ],
}

NEURON_ORDER = ["CA1-ProS", "CA2-FC-IG", "DG", "CA3", "Pvalb"]
GLIA_ORDER   = ["Endo", "Oligo", "Astro-TE", "OPC", "Microglia", "CHOR", "Astro-NT", "VLMC"]

# ── Helpers ───────────────────────────────────────────────────────────────────

def short_name(ct):
    core = " ".join(ct.split(" ")[1:])
    return core.replace(" Glut", "").replace(" Gaba", "").replace(" NN", "")

sig["src_short"] = sig.source.map(short_name)
sig["tgt_short"] = sig.target.map(short_name)

ann = pd.read_csv("data/cluster_annotation_term.csv")
subclass_colors = (
    ann[ann.cluster_annotation_term_set_name == "subclass"]
    .set_index("name")["color_hex_triplet"]
)
all_cts = sorted(set(sig.source) | set(sig.target))
node_colors = {short_name(ct): subclass_colors.get(ct, "#888888") for ct in all_cts}

UP_COLOR     = "#d62728"
DN_COLOR     = "#1f77b4"
LWIDTH_SCALE = 5.6
NODE_SIZE    = 1760

def bipartite_pos(left_nodes, right_nodes):
    pos = {}
    for i, n in enumerate(left_nodes):
        pos[n] = np.array([0.0, 1 - i / max(len(left_nodes)  - 1, 1)])
    for i, n in enumerate(right_nodes):
        pos[n] = np.array([1.0, 1 - i / max(len(right_nodes) - 1, 1)])
    return pos

def draw_panel(ax, theme_name, theme_lr_pairs):
    subset = sig[sig.interaction.isin(theme_lr_pairs)]
    if subset.empty:
        ax.text(0.5, 0.5, "No interactions", ha="center", va="center",
                transform=ax.transAxes, color="#cccccc", fontsize=9)
        ax.set_title(theme_name, fontsize=18, fontweight="bold", pad=8,
                     color="white",
                     bbox=dict(boxstyle="square,pad=0.45", facecolor=PANEL_COLOR,
                               edgecolor="none"))
        ax.axis("off")
        return

    present      = set(subset.src_short) | set(subset.tgt_short)
    neuron_nodes = [n for n in NEURON_ORDER if n in present]
    glia_nodes   = [n for n in GLIA_ORDER   if n in present]
    pos          = bipartite_pos(neuron_nodes, glia_nodes)

    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-0.08, 1.12)
    ax.axis("off")

    STRAIGHT_ARROWS = {"Tnr^Itgav_Itgb3"}

    # ── edges ──────────────────────────────────────────────────────────────────
    pair_counts = defaultdict(int)
    for _, row in subset.iterrows():
        pair_counts[(row.src_short, row.tgt_short)] += 1
    pair_drawn = defaultdict(int)

    def get_rad(key, count):
        idx = pair_drawn[key]
        if count == 1:
            return 0.10
        return float(np.linspace(-0.25, 0.25, count)[idx])

    edge_labels = []

    for _, row in subset.sort_values("interaction_stat").iterrows():
        s, t  = row.src_short, row.tgt_short
        color = UP_COLOR if row.interaction_stat > 0 else DN_COLOR
        lw    = max(1.0, abs(row.interaction_stat) * LWIDTH_SCALE)
        key   = (s, t)
        rad   = 0.0 if row.interaction in STRAIGHT_ARROWS else get_rad(key, pair_counts[key])
        pair_drawn[key] += 1

        x0, y0 = pos[s]
        x1, y1 = pos[t]

        ax.add_patch(FancyArrowPatch(
            posA=(x0, y0), posB=(x1, y1),
            arrowstyle="-|>",
            connectionstyle=f"arc3,rad={rad:.3f}",
            color=color, linewidth=lw,
            mutation_scale=32, shrinkA=36, shrinkB=36,
            zorder=2, alpha=1.0,
        ))

        my = (y0 + y1) / 2 - rad * (x1 - x0) * 0.5
        label = row.interaction.replace("^", " : ")
        edge_labels.append((0.5, my, label, color))

    # ── edge labels ────────────────────────────────────────────────────────────
    MANUAL_LABEL_POS = {
        "Gpc3 : Unc5d": (0.5, 0.0),
    }

    texts = []
    anchors = []
    for mx, my, label, color in edge_labels:
        x_init, y_init = MANUAL_LABEL_POS.get(label, (mx, my))
        t = ax.text(
            x_init, y_init, label,
            fontsize=14, color="white", fontweight="bold",
            ha="center", va="center", zorder=5,
            path_effects=[pe.withStroke(linewidth=0.5, foreground="black")],
            bbox=dict(boxstyle="round,pad=0.22", facecolor=color,
                      edgecolor="none", alpha=0.92),
        )
        if label not in MANUAL_LABEL_POS:
            texts.append(t)
            anchors.append(my)

    # ── nodes ─────────────────────────────────────────────────────────────────
    scatter_artists = []
    for node, (x, y) in pos.items():
        color = node_colors.get(node, "#888888")
        sc = ax.scatter(x, y, s=NODE_SIZE, c=color, zorder=6,
                        linewidths=1.5, edgecolors="white")
        scatter_artists.append(sc)

    # nudge overlapping labels upward in anchor-y order, max ±0.12 from anchor
    if len(texts) > 1:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        MIN_GAP = 0.015
        MAX_SHIFT = 0.12
        order = np.argsort(anchors)
        for _ in range(60):
            moved = False
            for i in range(len(order) - 1):
                a, b = order[i], order[i + 1]
                ya = texts[a].get_position()[1]
                yb = texts[b].get_position()[1]
                overlap = MIN_GAP - (yb - ya)
                if overlap > 0:
                    shift = overlap / 2
                    new_ya = ya - shift
                    new_yb = yb + shift
                    if abs(new_ya - anchors[a]) <= MAX_SHIFT:
                        texts[a].set_position((0.5, new_ya))
                    if abs(new_yb - anchors[b]) <= MAX_SHIFT:
                        texts[b].set_position((0.5, new_yb))
                    moved = True
            if not moved:
                break

    # ── node labels ───────────────────────────────────────────────────────────
    for node, (x, y) in pos.items():
        ha     = "right" if x == 0 else "left"
        offset = -0.14 if x == 0 else 0.14
        ax.text(x + offset, y, node, ha=ha, va="center",
                fontsize=18, fontweight="bold", color="#222222", zorder=7)

    # title: white text on a #732B8B rectangle
    ax.set_title(theme_name, fontsize=18, fontweight="bold", pad=8,
                 color="white",
                 bbox=dict(boxstyle="square,pad=0.45", facecolor=PANEL_COLOR,
                           edgecolor="none"))

# ── Figure ────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(20, 12))
fig.patch.set_facecolor("white")

for ax in axes.flat:
    ax.set_facecolor("white")

theme_items = list(THEMES.items())
for ax, (theme_name, lr_pairs) in zip(axes.flat, theme_items):
    draw_panel(ax, theme_name, lr_pairs)

# shared legend
dir_handles = [
    mpatches.Patch(color=UP_COLOR, label="Upregulated in PPD"),
    mpatches.Patch(color=DN_COLOR, label="Downregulated in PPD"),
]
lw_handles = [
    Line2D([0], [0], color="#aaaaaa", linewidth=0.5 * LWIDTH_SCALE, label="|score| = 0.5"),
    Line2D([0], [0], color="#aaaaaa", linewidth=1.0 * LWIDTH_SCALE, label="|score| = 1.0"),
    Line2D([0], [0], color="#aaaaaa", linewidth=1.5 * LWIDTH_SCALE, label="|score| = 1.5"),
]
fig.legend(handles=dir_handles + lw_handles,
           loc="lower center", ncol=5, frameon=False,
           fontsize=11, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(
    "Neuron ↔ glia cell-cell communication  (LIANA+, padj < 0.05)",
    fontsize=15, fontweight="bold", y=1.01,
)
plt.tight_layout(w_pad=3, h_pad=3)


out_pdf = "figures/pdfs/fig_network_ccc_v8.pdf"
out_svg = "figures/fig_network_ccc_v8.svg"
with PdfPages(out_pdf) as pdf:
    pdf.savefig(fig, bbox_inches="tight")
fig.savefig(out_svg, bbox_inches="tight")

print(f"Saved: {out_pdf}")
print(f"Saved: {out_svg}")
