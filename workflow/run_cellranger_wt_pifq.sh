#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(pwd)}"
REF="${PROJECT_DIR}/references/Arabidopsis_thaliana_Araport11_CellRanger_reference"
FASTQS="${PROJECT_DIR}/data/raw/CRA010863"
OUTDIR="${PROJECT_DIR}/cellranger"

if ! command -v cellranger >/dev/null 2>&1; then
  echo "ERROR: cellranger is not available in PATH." >&2
  echo "Install Cell Ranger from 10x Genomics or load it on the server, then rerun." >&2
  exit 1
fi

for required in "$REF/reference.json" "$FASTQS/ScWT_S1_L001_R1_001.fastq.gz" "$FASTQS/Scpifq_S1_L001_R1_001.fastq.gz"; do
  if [ ! -e "$required" ]; then
    echo "ERROR: missing required input: $required" >&2
    exit 1
  fi
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"

run_if_needed() {
  local id="$1"
  local sample="$2"

  if [ -s "$id/outs/metrics_summary.csv" ]; then
    echo "[$(date '+%F %T')] $id already finished; skipping."
    return 0
  fi

  echo "[$(date '+%F %T')] Running Cell Ranger for $id"
  cellranger count \
    --localcores=80 \
    --id="$id" \
    --fastqs="$FASTQS" \
    --sample="$sample" \
    --transcriptome="$REF" \
    --no-bam
}

run_if_needed Scpifq Scpifq
run_if_needed ScWT ScWT
