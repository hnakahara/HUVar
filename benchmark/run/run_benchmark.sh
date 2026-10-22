#!/usr/bin/env bash
# Re-run the eRepo benchmark with the revised HUVar (Genome Medicine revision).
#
# Same commands as the original benchmark, plus --exclude-self-expert-panel on the
# automated-only runs. Run from the HUVar checkout:
#   (from a HUVar checkout installed with pip install -e .)
#   nohup bash run_benchmark.sh GRCh37 > bench_GRCh37.log 2>&1 &
#   nohup bash run_benchmark.sh GRCh38 > bench_GRCh38.log 2>&1 &
# With no argument both assemblies run one after the other.
#
# Runs whose result file already exists are skipped, so an interrupted run can
# simply be restarted. Each run writes to a temporary name and is renamed only
# when it finishes, so a half-written file is never mistaken for a result.
# FORCE=1 re-runs everything.
#
# Configurations (names unchanged so the analysis scripts pick them up):
#   default    ESM1b + OpenSpliceAI (80 nt)
#   spliceai   AlphaMissense + SpliceAI (no BayesDel/CADD)
#   revel      REVEL + SpliceAI + BayesDel/CADD auxiliary rules (VCEP-oriented)
#   *_manual   same + eRepo curated-evidence supplement (ClinVar expert-panel criteria kept)
#
# --with-tp53-functional (TP53 PS3/BS3 from the NCI TP53 Database assay classes,
# non-commercial terms) is used only in the non-commercial configurations
# (spliceai = AlphaMissense, revel), like BayesDel/CADD; the default configuration
# stays commercial-use ready.
# Automated-only runs use --exclude-self-expert-panel: criteria the ClinGen expert
# panel applied to the variant itself (PS3/PS4/PP1/BS2 in ClinVar >=3-star records)
# are not imported, so the evaluation does not read the reference curation.
set -euo pipefail

DATA=${DATA:?set DATA to the HUVar data directory (python scripts/setup_data.py --data-dir ...)}
SPLICEAI_DIR=${SPLICEAI_DIR:?set SPLICEAI_DIR to the directory with the precomputed SpliceAI scores}
REPO=$(cd "$(dirname "$0")/../.." && pwd)   # HUVar repository root
IN=${IN:-$REPO/inputs}   # erepo_hg38_may_2026.vcf / erepo_hg19_may_2026.vcf
OUT=${OUT:-$PWD/erepo_results}
WORKERS=${WORKERS:-36}
FORCE=${FORCE:-0}
mkdir -p "$OUT"
ASSEMBLIES=("$@")
[ $# -eq 0 ] && ASSEMBLIES=(GRCh37 GRCh38)

# Record the environment so that runs on different servers can be checked to be identical.
echo "[$(date '+%F %T')] host=$(hostname) assemblies=${ASSEMBLIES[*]} out=$OUT"
echo "  HUVar commit: $(git -C "${HUVAR_DIR:-$PWD}" rev-parse --short HEAD 2>/dev/null || echo unknown) (${HUVAR_DIR:-$PWD})"
for f in "$DATA"/GRCh3?/clinvar/clinvar_ps1_pm5_GRCh3?.sqlite "$DATA"/shared/clingen_gene_disease_validity.csv; do
  [ -e "$f" ] && echo "  $(ls -l --time-style=+%F_%T "$f" | awk '{print $5, $6, $7}')"
done
echo "  python: $(command -v python)  acmg-classify: $(command -v acmg-classify)"
echo "  data: $DATA"

run() {  # run <hg> <asm> <cfg> <manual:0|1>
  local hg=$1 asm=$2 cfg=$3 manual=$4 tools extra name
  case $cfg in
    default)  tools="--insilico-tool esm1b --splice-tool openspliceai --openspliceai-flanking-size 80" ;;
    spliceai) tools="--insilico-tool alphamissense --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-tp53-functional" ;;
    revel)    tools="--insilico-tool revel --splice-tool spliceai --spliceai-dir $SPLICEAI_DIR --with-bayesdel --with-cadd --with-tp53-functional" ;;
  esac
  if [ "$manual" = 1 ]; then
    extra="--supplement $DATA/shared/erepo_manual_criteria_${hg}.tsv"
    name="erepo_${hg}_huhvar_${cfg}_manual_results.tsv"
  else
    extra="--exclude-self-expert-panel"
    name="erepo_${hg}_huhvar_${cfg}_results.tsv"
  fi
  if [ "$FORCE" != 1 ] && [ -s "$OUT/$name" ]; then
    echo "[$(date '+%F %T')] skip (exists): $name"
    return 0
  fi
  local tmp="$OUT/.tmp_${name}"
  rm -f "$tmp" "$tmp.provenance.json"
  echo "[$(date '+%F %T')] start: $name"
  # shellcheck disable=SC2086
  acmg-classify classify "$IN/erepo_${hg}_may_2026.vcf" -o "$tmp" \
    --data-dir "$DATA" --workers "$WORKERS" --assembly "$asm" --no-progress \
    $tools $extra
  [ -e "$tmp.provenance.json" ] && mv -f "$tmp.provenance.json" "$OUT/$name.provenance.json"
  mv -f "$tmp" "$OUT/$name"
  for side in "$tmp"*; do  # other sidecars (e.g. skipped records)
    [ -e "$side" ] && mv -f "$side" "$OUT/${name}${side#"$tmp"}"
  done
  echo "[$(date '+%F %T')] done:  $name"
}

for asm in "${ASSEMBLIES[@]}"; do
  case $asm in
    GRCh37) hg=hg19 ;;
    GRCh38) hg=hg38 ;;
    *) echo "unknown assembly: $asm (use GRCh37 or GRCh38)" >&2; exit 1 ;;
  esac
  for cfg in default spliceai revel; do
    run "$hg" "$asm" "$cfg" 0
    run "$hg" "$asm" "$cfg" 1
  done
done

echo "[$(date '+%F %T')] done -> $OUT"
ls -la "$OUT"
