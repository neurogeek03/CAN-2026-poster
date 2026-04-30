"""
compose_poster.py — Assemble all figure SVGs into a single poster preview.

Usage:
    python scripts/compose_poster.py             # writes figures/poster_preview.svg + .png
    python scripts/compose_poster.py --svg-only  # skip PNG export (no cairosvg needed)

Dependencies:
    uv pip install svgutils cairosvg
"""

import argparse
from pathlib import Path
import svgutils.transform as sg

# ---------------------------------------------------------------------------
# Poster canvas (inches → points at 96 dpi for SVG coordinate space)
# Standard academic poster: 48" × 36"
# ---------------------------------------------------------------------------
PX_PER_INCH = 96
POSTER_W = 48 * PX_PER_INCH   # 4608 px
POSTER_H = 36 * PX_PER_INCH   # 3456 px

FIGURES_DIR = Path("figures")
OUT_SVG = FIGURES_DIR / "poster_preview.svg"
OUT_PNG = FIGURES_DIR / "poster_preview.png"

# ---------------------------------------------------------------------------
# Layout config — edit this to reposition panels
#
# Each entry: (figure_file, x, y, scale)
#   x, y   — top-left corner in poster pixels
#   scale  — scale factor applied to the SVG (1.0 = native size)
#
# Grid guide (rough thirds × halves):
#   Column 1: x ≈ 50          Column 2: x ≈ 1580       Column 3: x ≈ 3110
#   Row top:  y ≈ 200          Row bottom: y ≈ 1850
# ---------------------------------------------------------------------------
COL = [50, 1580, 3110]
ROW = [200, 1850]

LAYOUT = [
    # (filename,                    x,       y,       scale)
    ("fig01_umap.svg",          COL[0], ROW[0], 1.0),   # Panel A — UMAP
    ("fig02_markers.svg",       COL[1], ROW[0], 1.0),   # Panel B — Marker plot
    ("fig03_cumulative_degs.svg", COL[2], ROW[0], 1.0), # Panel C — Cumulative DEGs
    ("fig04_volcano.svg",       COL[0], ROW[1], 1.0),   # Panel D — Volcano
    ("fig05_spatial.svg",       COL[1], ROW[1], 1.0),   # Panel E — Spatial expression
    ("fig06_ccc.svg",           COL[2], ROW[1], 1.0),   # Panel F — Cell-cell communication
]

PANEL_LABELS = ["A", "B", "C", "D", "E", "F"]
LABEL_OFFSET_X = 0    # px relative to panel x
LABEL_OFFSET_Y = -30  # px above panel top edge
LABEL_SIZE = 48       # font size in px


def load_svg(path: Path) -> sg.SVGFigure | None:
    if not path.exists():
        print(f"  [skip] {path.name} not found")
        return None
    return sg.fromfile(str(path))


def compose(svg_only: bool = False) -> None:
    FIGURES_DIR.mkdir(exist_ok=True)

    canvas = sg.SVGFigure(width=f"{POSTER_W}px", height=f"{POSTER_H}px")
    all_elements = []

    for (fname, x, y, scale), label in zip(LAYOUT, PANEL_LABELS):
        fig_path = FIGURES_DIR / fname
        fig = load_svg(fig_path)
        if fig is None:
            continue

        root = fig.getroot()
        root.moveto(x, y, scale)
        all_elements.append(root)

        # Panel label (A, B, C …)
        txt = sg.TextElement(
            x + LABEL_OFFSET_X,
            y + LABEL_OFFSET_Y,
            label,
            size=LABEL_SIZE,
            weight="bold",
            font="Helvetica",
        )
        all_elements.append(txt)

    if not all_elements:
        print("No figures found — nothing to compose. Run figure scripts first.")
        return

    canvas.set_elements(all_elements)
    canvas.save(str(OUT_SVG))
    print(f"Saved: {OUT_SVG}")

    if not svg_only:
        try:
            import cairosvg
            cairosvg.svg2png(
                url=str(OUT_SVG),
                write_to=str(OUT_PNG),
                output_width=POSTER_W // 2,   # half-res PNG for quick review
                output_height=POSTER_H // 2,
            )
            print(f"Saved: {OUT_PNG}")
        except ImportError:
            print("cairosvg not installed — skipping PNG export. Run with --svg-only to suppress this.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--svg-only", action="store_true", help="Skip PNG export")
    args = parser.parse_args()
    compose(svg_only=args.svg_only)
