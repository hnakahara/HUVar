#!/usr/bin/env bash
# Main and supplementary figures (Arial). Run from this folder: bash make_figs.sh <CGP supplementary table (xlsx)> [workdir]
# Writes PDF + PNG of every figure to <workdir>/out and the A4 "Fig. N" pages to <workdir>/final.
set -euo pipefail
S3=$1; W=${2:-$PWD/figure_work}; HERE=$(cd "$(dirname "$0")" && pwd)
rm -rf "$W"; mkdir -p "$W"; cp -r "$HERE"/. "$W"/; cd "$W"; mkdir -p out final
for s in make_criteria_barplot.py make_fig3_combined.py make_code_ratio_barplot.py; do python3 "$s"; done
python3 make_fig_cgp_barplot.py --s3 "$S3" --results-dir cgp --out out/cgp_barplot_hg38
python3 make_fig_cgp_sankey.py --s3 "$S3" --results-dir cgp --f1cdx cgp/f1cdx_report_classes.tsv --out out/cgp_sankey3_hg38
for b in hg38 hg19; do (cd confusion_$b && python3 ../make_fig_gene_concordance.py --summary confusion_gene_summary.tsv --config revel_manual --min-n 10 --out ../out/fig_gene_concordance_revel_manual_$b); done
python3 compose.py
for f in final/*.pdf; do pdftoppm -png -r 300 -singlefile "$f" "${f%.pdf}" 2>/dev/null || true; done
ls out final
