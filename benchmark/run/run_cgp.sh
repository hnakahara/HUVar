#!/usr/bin/env bash
# Re-run the CGP (FoundationOne CDx) cohort with the revised HUVar for Reviewer 3, comment 19.
# Same commands as the original analysis (eRepo supplement included), written to a new folder.
#   (from a HUVar checkout installed with pip install -e .)
#   nohup bash run_cgp_revision.sh GRCh37 > cgp_GRCh37.log 2>&1 &
#   nohup bash run_cgp_revision.sh GRCh38 > cgp_GRCh38.log 2>&1 &
# Existing result files are skipped (FORCE=1 to re-run).
set -euo pipefail
DATA=${DATA:?set DATA to the HUVar data directory (python scripts/setup_data.py --data-dir ...)}
SPLICEAI_DIR=${SPLICEAI_DIR:?set SPLICEAI_DIR to the directory with the precomputed SpliceAI scores}
IN=${IN:?set IN to the folder with f1_acmg_<build>.vcf (CGP cohort VCF; not distributed)}
OUT=${OUT:-$PWD/cgp_results}
WORKERS=${WORKERS:-24}
FORCE=${FORCE:-0}
mkdir -p "$OUT"
ASSEMBLIES=("$@"); [ $# -eq 0 ] && ASSEMBLIES=(GRCh37 GRCh38)
echo "[$(date '+%F %T')] host=$(hostname) commit=$(git -C "${HUVAR_DIR:-$PWD}" rev-parse --short HEAD 2>/dev/null || echo unknown) out=$OUT"

run() {  # run <hg> <asm> <cfg>
  local hg=$1 asm=$2 cfg=$3 tools name tbase
  case $cfg in
    default)       tools="--insilico-tool esm1b --splice-tool openspliceai --openspliceai-flanking-size 80" ;;
    alphamissense) tools="--insilico-tool alphamissense --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-tp53-functional" ;;
    revel)         tools="--insilico-tool revel --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-bayesdel --with-cadd --with-tp53-functional" ;;
  esac
  name="f1_acmg_${hg}_${cfg}_results.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$OUT/$name" ]; then echo "skip (exists): $name"; return 0; fi
  tbase="$OUT/.tmp_${name%.tsv}"; rm -f "$tbase".*
  echo "[$(date '+%F %T')] start: $name"
  # shellcheck disable=SC2086
  acmg-classify classify "$IN/f1_acmg_${hg}.vcf" -o "$tbase.tsv" \
    --data-dir "$DATA" --workers "$WORKERS" --assembly "$asm" --no-progress \
    --supplement "$DATA/shared/erepo_manual_criteria_${hg}.tsv" $tools
  for f in "$tbase".*; do
    [ "$f" = "$tbase.tsv" ] && continue
    [ -e "$f" ] && mv -f "$f" "$OUT/${name%.tsv}${f#"$tbase"}"
  done
  mv -f "$tbase.tsv" "$OUT/$name"
  echo "[$(date '+%F %T')] done:  $name"
}
for asm in "${ASSEMBLIES[@]}"; do
  case $asm in GRCh37) hg=hg19 ;; GRCh38) hg=hg38 ;; *) echo "unknown assembly: $asm" >&2; exit 1 ;; esac
  for cfg in default alphamissense revel; do run "$hg" "$asm" "$cfg"; done
done
ls -la "$OUT"
