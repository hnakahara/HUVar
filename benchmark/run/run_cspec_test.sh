#!/usr/bin/env bash
# Reviewer 2, Comment 6: effect of the disease-specific ClinGen CSpecs for genes
# with several released specifications (RYR1: malignant hyperthermia [MH VCEP],
# congenital myopathy AD/AR [Congenital Myopathies VCEP]; ACTA1: congenital
# myopathy AD/AR; VWF: VWD AD [types 1/2] / AR [type 3]).
# The eRepo variants of these genes are classified (GRCh38) without --cspec (the
# conservative merged row = HUVar default) and under each specification
# (--cspec GENE=ID). Three runs cover every specification:
#   run A: RYR1=malignant_hyperthermia  ACTA1=congenital_myopathy_ad  VWF=vwd_ad
#   run B: RYR1=congenital_myopathy_ad  ACTA1=congenital_myopathy_ar  VWF=vwd_ar
#   run C: RYR1=congenital_myopathy_ar
# Configurations: revel_manual (REVEL+SpliceAI+BayesDel/CADD + eRepo supplement),
# default_manual (ESM1b+OpenSpliceAI + supplement), revel (automated only).
# Then:  python compare_cspec.py --dir $OUT
#   (from a HUVar checkout installed with pip install -e .)
#   nohup bash run_cspec_test.sh > cspec.log 2>&1 &
set -euo pipefail
DATA=${DATA:?set DATA to the HUVar data directory (python scripts/setup_data.py --data-dir ...)}
SPLICEAI_DIR=${SPLICEAI_DIR:?set SPLICEAI_DIR to the directory with the precomputed SpliceAI scores}
REPO=$(cd "$(dirname "$0")/../.." && pwd)   # HUVar repository root
IN=${IN:-$REPO/data}   # erepo_hg38.vcf / erepo_hg19.vcf
OUT=${OUT:-$PWD/erepo_cspec_results}
WORKERS=${WORKERS:-16}
FORCE=${FORCE:-0}
mkdir -p "$OUT"
VCF="$OUT/erepo_hg38_multispec_genes.vcf"
{ grep '^#' "$IN/erepo_hg38.vcf"; grep -v '^#' "$IN/erepo_hg38.vcf" | awk -F'\t' '$3 ~ /^(RYR1|ACTA1|VWF):/'; } > "$VCF"
echo "variants: $(grep -vc '^#' "$VCF")  commit: $(git -C "${HUVAR_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}" rev-parse --short HEAD 2>/dev/null)"

declare -A RUNS=(
  [merged]=""
  [A]="--cspec RYR1=malignant_hyperthermia --cspec ACTA1=congenital_myopathy_ad --cspec VWF=vwd_ad"
  [B]="--cspec RYR1=congenital_myopathy_ad --cspec ACTA1=congenital_myopathy_ar --cspec VWF=vwd_ar"
  [C]="--cspec RYR1=congenital_myopathy_ar"
)
for cfg in revel_manual default_manual revel; do
  case $cfg in
    default_manual) tools="--insilico-tool esm1b --splice-tool openspliceai --openspliceai-flanking-size 80 --supplement $DATA/shared/erepo_manual_criteria_hg38.tsv" ;;
    revel_manual)   tools="--insilico-tool revel --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-bayesdel --with-cadd --supplement $DATA/shared/erepo_manual_criteria_hg38.tsv" ;;
    revel)          tools="--insilico-tool revel --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-bayesdel --with-cadd --exclude-self-expert-panel" ;;
  esac
  for run in merged A B C; do
    name="cspec_${cfg}_${run}.tsv"
    if [ "$FORCE" != 1 ] && [ -s "$OUT/$name" ]; then echo "skip (exists): $name"; continue; fi
    echo "[$(date '+%F %T')] $name"
    # shellcheck disable=SC2086
    acmg-classify classify "$VCF" -o "$OUT/$name" --data-dir "$DATA" --workers "$WORKERS" \
      --assembly GRCh38 --no-progress $tools ${RUNS[$run]}
  done
done
ls -la "$OUT"
