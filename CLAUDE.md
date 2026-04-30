# CAN 2026 Poster — Project Instructions

## What this project is
Reproducible figure generation for a conference poster on postpartum depression (PPD),
using two spatial transcriptomics datasets: Slide-tags and Slide-seq.

## Environment
- Python venv managed with `uv`, located at `.venv/`
- Activate: `source .venv/bin/activate`
- Reproduce: `uv pip install -r venv.log`
- After adding any new package: `uv pip freeze > venv.log`

## Rules
- All figure scripts live in `scripts/figures/` — one script per figure.
- All extraction scripts live in `scripts/extract/` — one script per dataset.
- Figure scripts must only read from `data/processed/` (CSVs). Never load raw `.h5ad` files in figure scripts.
- All figures output to `figures/` as SVG.
- Python only — no R.
- Shared style settings (fonts, colors, dimensions) go in `scripts/figures/style.py`.

## Key files
- `docs/plan.md` — figure inventory, data requirements, status tracking
- `scripts/run_all.py` — wrapper that runs the full pipeline

## Data
- `data/raw/` — large source files (h5ad, etc.), gitignored
- `data/processed/` — extracted CSVs, gitignored
- Both directories are on-disk only; not tracked in git.
