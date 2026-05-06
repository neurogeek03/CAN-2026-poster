#!/usr/bin/env python3
"""
Generate individual spatial-expression figures and assemble them into a grid panel.

Calls fig_spatial_expression_v3.py for each (cell_type, gene) entry,
then stitches the resulting PDFs into an NxM panel with fitz (PyMuPDF).

Usage:
    uv run python scripts/make_panel_colored.py
"""

import subprocess
import sys
from pathlib import Path

import fitz  # PyMuPDF

ROOT = Path(__file__).resolve().parents[1]
FIG_SCRIPT = ROOT / "scripts" / "figure_code" / "fig_spatial_expression_v3.py"

# ── Panel definition ──────────────────────────────────────────────────────────
# (cell_type, gene) — order is left-to-right, top-to-bottom
# Row 1: Sgk1 | Ddit4 | Fkbp5
# Row 2: Nfkbia | Pdk4 | Ptn

CELL_TYPE = "Oligo_NN"

GENES = [
    "Sgk1",
    "Ddit4",
    "Fkbp5",
    "Nfkbia",
    "Pdk4",
    "Ptn",
]

NCOLS              = 3
NROWS              = 2
CELL_TYPE_LABEL    = False  # True to show cell type above colorbar, False to hide
FIG_DIR  = ROOT / "figures" / "genes_half_half" / CELL_TYPE
OUTFILE  = FIG_DIR / f"panel_{CELL_TYPE}_{NCOLS}x{NROWS}.pdf"


def figure_path(cell_type: str, gene: str) -> Path:
    return ROOT / "figures" / "genes_half_half" / f"fig_spatial_{gene}_{cell_type}_v3.pdf"


# ── Step 1: generate figures ──────────────────────────────────────────────────

def generate_figures():
    for gene in GENES:
        out = figure_path(CELL_TYPE, gene)
        print(f"Generating {gene} ({CELL_TYPE})…")
        cmd = ["uv", "run", "python", str(FIG_SCRIPT), CELL_TYPE, gene]
        if not CELL_TYPE_LABEL:
            cmd.append("--no-cell-type-label")
        subprocess.run(cmd, check=True, cwd=ROOT)
        if not out.exists():
            print(f"  ERROR: expected output not found: {out}", file=sys.stderr)
            sys.exit(1)
        print(f"  Saved → {out}")


# ── Step 2: assemble panel ────────────────────────────────────────────────────

def assemble_panel():
    pdfs = [figure_path(CELL_TYPE, gene) for gene in GENES]

    cell_w, cell_h = None, None
    for p in pdfs:
        doc = fitz.open(p)
        r = doc[0].rect
        doc.close()
        if cell_w is None:
            cell_w, cell_h = r.width, r.height

    panel_w = cell_w * NCOLS
    panel_h = cell_h * NROWS

    out_doc = fitz.open()
    page = out_doc.new_page(width=panel_w, height=panel_h)

    for idx, pdf_path in enumerate(pdfs):
        col = idx % NCOLS
        row = idx // NCOLS
        x0 = col * cell_w
        y0 = row * cell_h
        cell_rect = fitz.Rect(x0, y0, x0 + cell_w, y0 + cell_h)

        src = fitz.open(pdf_path)
        page.show_pdf_page(cell_rect, src, 0)
        src.close()

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out_doc.save(str(OUTFILE), garbage=4, deflate=True)
    out_doc.close()
    print(f"\nPanel saved → {OUTFILE}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    generate_figures()
    assemble_panel()
