#!/usr/bin/env bash
# Usage: ./make_panel.sh <dir> <cols>x<rows> [outfile]
#   dir    — directory containing PDFs to panel
#   NxM    — grid layout (e.g. 3x2)
#   outfile — optional output filename (default: panel_<dirname>_<NxM>.pdf, saved in <dir>)
#
# Example:
#   ./make_panel.sh figures/genes_half_half/Oligo_NN 3x2

set -euo pipefail

DIR="${1:?Usage: $0 <dir> <NxM> [outfile]}"
NUP="${2:?Usage: $0 <dir> <NxM> [outfile]}"
DIRNAME=$(basename "$DIR")
OUTFILE="${3:-${DIR}/panel_${DIRNAME}_${NUP}.pdf}"

PDFS=("$DIR"/*.pdf)

if [[ ${#PDFS[@]} -eq 0 ]]; then
    echo "Error: no PDFs found in $DIR" >&2
    exit 1
fi

echo "Arranging ${#PDFS[@]} PDFs in a ${NUP} grid -> $OUTFILE"

pdfjam --nup "$NUP" --landscape --outfile "$OUTFILE" "${PDFS[@]}"

echo "Done: $OUTFILE"
