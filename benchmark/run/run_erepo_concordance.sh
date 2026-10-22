#!/usr/bin/env bash
# Agreement with the eRepo classifications (Reviewer 3, Comment 15), all eRepo records (none excluded).
#   bash run_erepo_concordance.sh <inputs dir> <out dir> <HUVar dir> [builds...]
# Example: bash run_erepo_concordance.sh ../inputs ../results/concordance /path/to/HUVar hg38 hg19
# <inputs dir> must contain erepo-original.csv, erepo_<build>.vcf and erepo_<build>_huhvar_<cfg>_results.tsv.
# Output in <out>/<build>/: concordance_*, diag_*, codes_reclass/ (for reference: eRepo records whose
#       classification is not reproduced from their own codes; not excluded from the analysis), confusion/, discordant/,
#       categories_counts.tsv, fig_gene_concordance_revel_manual_<build>.{png,pdf}
set -euo pipefail
IN=$1; OUT=$2; HUVAR=$3; shift 3
BUILDS=${*:-hg38 hg19}
HERE=$(cd "$(dirname "$0")/../analysis" && pwd)   # analysis scripts
CF="default default_manual revel revel_manual spliceai spliceai_manual"
for B in $BUILDS; do
  case $B in hg38) A=GRCh38;; hg19) A=GRCh37;; *) echo "unknown build $B"; exit 1;; esac
  O=$OUT/$B; mkdir -p "$O"
  ER="$HUVAR/resources/$A/erepo_classifications_$B.tsv"
  echo "== $B"
  python "$HERE/erepo_concordance.py" "$IN" "$O" "$B" > "$O/concordance.log" 2>&1
  python "$HERE/erepo_concordance_diag.py" "$IN" "$O" revel_manual "$B" > "$O/diag_revel_manual.txt" 2>&1
  python "$HERE/erepo_codes_reclass.py" --build "$B" --dir "$O" --erepo "$ER" \
      --huvar-src "$HUVAR/src" --results-dir "$IN" --out "$O/codes_reclass" > "$O/codes_reclass.log" 2>&1
  python "$HERE/erepo_confusion.py" --build "$B" --dir "$O" --configs $CF \
      --diag "$O/diag_revel_manual_variants.tsv" --out "$O/confusion" > "$O/confusion.log" 2>&1
  python "$HERE/erepo_discordant_table.py" --build "$B" --dir "$O" --results-dir "$IN" \
      --erepo "$ER" --out "$O/discordant" > "$O/discordant.log" 2>&1
  python "$HERE/erepo_category_counts.py" --build "$B" --dir "$O" \
      --out "$O/categories_counts.tsv" > "$O/categories_counts.log" 2>&1
  (cd "$O/confusion" && python "$HERE/make_fig_gene_concordance.py" --summary confusion_gene_summary.tsv \
      --config revel_manual --min-n 10 --out "fig_gene_concordance_revel_manual_$B") > /dev/null 2>&1 || true
  cat "$O/categories_counts.log"
done
