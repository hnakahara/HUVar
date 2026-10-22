#!/usr/bin/env bash
# Hardware requirements (Reviewer 2, Comment 2; README "Measured resource usage").
# Classifies the 11,545 eRepo variants (GRCh38) with the VCEP-oriented configuration
# (REVEL + SpliceAI + BayesDel/CADD, TP53 functional data, eRepo supplement) with
# 1-16 worker processes, without a memory limit and with the memory limited to
# 16 GB and 8 GB (cgroup MemoryMax via systemd-run). Wall time, CPU usage, maximum
# resident set size and major page faults are recorded with GNU time (/usr/bin/time -v).
#
#   DATA=/path/to/data SPLICEAI_DIR=/path/to/spliceai bash measure_time_memory.sh > measure.log 2>&1
set -uo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)   # HUVar repository root
DATA=${DATA:?set DATA to the HUVar data directory (python scripts/setup_data.py --data-dir ...)}
SPLICEAI_DIR=${SPLICEAI_DIR:?set SPLICEAI_DIR to the directory with the precomputed SpliceAI scores}
IN=${IN:-$REPO/data}
OUT=${OUT:-$PWD/resource_measurement}
WORKER_LIST=${WORKER_LIST:-"1 2 4 8 16"}
LIMITS=${LIMITS:-"none 16G 8G"}
mkdir -p "$OUT"
command -v /usr/bin/time >/dev/null || { echo "GNU time (/usr/bin/time) is required" >&2; exit 1; }

for lim in $LIMITS; do
  for w in $WORKER_LIST; do
    pre=()
    [ "$lim" != none ] && pre=(systemd-run --user --scope --quiet -p "MemoryMax=$lim" -p MemorySwapMax=0)
    echo "== workers=$w memory=$lim"
    "${pre[@]}" /usr/bin/time -v acmg-classify classify "$IN/erepo_hg38.vcf" \
        -o "$OUT/erepo_hg38_w${w}_${lim}.tsv" --data-dir "$DATA" --workers "$w" --assembly GRCh38 \
        --insilico-tool revel --splice-tool spliceai --spliceai-dir "$SPLICEAI_DIR" \
        --supplement "$DATA/shared/erepo_manual_criteria_hg38.tsv" \
        --with-bayesdel --with-cadd --with-tp53-functional --no-progress \
      2> "$OUT/time_w${w}_${lim}.txt"
    grep -E "Elapsed|Percent of CPU|Maximum resident|Major" "$OUT/time_w${w}_${lim}.txt"
  done
done
