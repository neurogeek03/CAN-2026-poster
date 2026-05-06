#!/usr/bin/env python3
"""
Generate individual spatial-expression figures and assemble them into a grid panel.

Layout: 2 rows × 4 columns
  Row 1 (Hippocampus): Nr3c1 | Nr3c2 | Nr4a3 | Ifngr1
  Row 2 (DG_Glut):     Sgk1  | Bdnf  | Pcdh7 | Pcdh15

Calls fig_spatial_expression_v3.py for each (cell_type, gene) entry,
then stitches the resulting PDFs into a panel with fitz (PyMuPDF).

Usage:
    uv run python scripts/make_panel_hippo.py
"""

import subprocess
import sys
from pathlib import Path

import fitz  # PyMuPDF

ROOT       = Path(__file__).resolve().parents[1]
FIG_SCRIPT = ROOT / "scripts" / "figure_code" / "fig_spatial_expression_v3.py"
OUTFILE    = ROOT / "figures" / "genes_half_half" / "panel_hippo_4x2.pdf"

# ── Panel definition ──────────────────────────────────────────────────────────
# Each entry: (cell_type, [gene, ...]) — defines one row, left to right.

ROWS = [
    ("Hippocampus", ["Nr3c1", "Nr3c2", "Nr4a3", "Ifngr1"]),
    ("DG_Glut",     ["Sgk1",  "Bdnf",  "Pcdh7", "Pcdh15"]),
]

NCOLS           = 4
NROWS           = len(ROWS)
CELL_TYPE_LABEL = False


def figure_path(cell_type: str, gene: str) -> Path:
    return (
        ROOT / "figures" / "genes_half_half" / cell_type
        / f"fig_spatial_{gene}_{cell_type}_v3.pdf"
    )


# ── Step 1: generate figures ──────────────────────────────────────────────────

def generate_figures():
    for cell_type, genes in ROWS:
        for gene in genes:
            out = figure_path(cell_type, gene)
            print(f"Generating {gene} ({cell_type})…")
            cmd = ["uv", "run", "python", str(FIG_SCRIPT), cell_type, gene]
            if not CELL_TYPE_LABEL:
                cmd.append("--no-cell-type-label")
            subprocess.run(cmd, check=True, cwd=ROOT)
            if not out.exists():
                print(f"  ERROR: expected output not found: {out}", file=sys.stderr)
                sys.exit(1)
            print(f"  Saved → {out}")


# ── Step 2: assemble panel ────────────────────────────────────────────────────

def assemble_panel():
    pdfs = [figure_path(ct, gene) for ct, genes in ROWS for gene in genes]

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

    OUTFILE.parent.mkdir(parents=True, exist_ok=True)
    out_doc.save(str(OUTFILE), garbage=4, deflate=True)
    out_doc.close()
    print(f"\nPanel saved → {OUTFILE}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    generate_figures()
    assemble_panel()
