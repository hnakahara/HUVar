#!/usr/bin/env bash
# Reviewer 1: PS1/PM5 with ClinVar comparators of >=2 stars instead of >=1.
# Re-runs the eRepo configurations (GRCh38) with ACMG_PS1_MIN_STARS=2
# ACMG_PM5_MIN_STARS=2; compare with the standard runs from run_benchmark.sh.
#   MODE=manual (default): with the eRepo curated-evidence supplement
#                          (erepo_hg38_huhvar_{default,revel}_manual_stars2_results.tsv)
#   MODE=auto:             automated only (--exclude-self-expert-panel)
# Existing result files are skipped (FORCE=1 to re-run).
#   (from a HUVar checkout installed with pip install -e .)
#   nohup bash run_stars2_test.sh > stars2.log 2>&1 &
set -euo pipefail
DATA=${DATA:?set DATA to the HUVar data directory (python scripts/setup_data.py --data-dir ...)}
SPLICEAI_DIR=${SPLICEAI_DIR:?set SPLICEAI_DIR to the directory with the precomputed SpliceAI scores}
REPO=$(cd "$(dirname "$0")/../.." && pwd)   # HUVar repository root
IN=${IN:-$REPO/inputs}   # erepo_hg38_may_2026.vcf
OUT=${OUT:-$PWD/erepo_stars2_$(date +%Y%m%d)}
WORKERS=${WORKERS:-36}
MODE=${MODE:-manual}
FORCE=${FORCE:-0}
mkdir -p "$OUT"
export ACMG_PS1_MIN_STARS=2 ACMG_PM5_MIN_STARS=2
echo "[$(date '+%F %T')] host=$(hostname) commit=$(git -C "${HUVAR_DIR:-$PWD}" rev-parse --short HEAD 2>/dev/null || echo unknown) out=$OUT"
for cfg in default revel; do
  case $cfg in
    default) tools="--insilico-tool esm1b --splice-tool openspliceai --openspliceai-flanking-size 80" ;;
    revel)   tools="--insilico-tool revel --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-bayesdel --with-cadd --with-tp53-functional" ;;
  esac
  if [ "$MODE" = manual ]; then
    extra="--supplement $DATA/shared/erepo_manual_criteria_hg38.tsv"
    name="erepo_hg38_huhvar_${cfg}_manual_stars2_results.tsv"
  else
    extra="--exclude-self-expert-panel"
    name="erepo_hg38_huhvar_${cfg}_stars2_results.tsv"
  fi
  if [ "$FORCE" != 1 ] && [ -s "$OUT/$name" ]; then echo "skip (exists): $name"; continue; fi
  echo "[$(date '+%F %T')] start: $name"
  # shellcheck disable=SC2086
  acmg-classify classify "$IN/erepo_hg38_may_2026.vcf" -o "$OUT/$name" \
    --data-dir "$DATA" --workers "$WORKERS" --assembly GRCh38 --no-progress \
    $extra $tools
  echo "[$(date '+%F %T')] done:  $name"
done
ls -la "$OUT"
