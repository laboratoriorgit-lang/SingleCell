#!/usr/bin/env bash
# OPTIONAL: runs the workflow with the values of params.yaml.
# The step scripts in workflow/ (and the report) contain plain values and run on
# their own. This script first generates, with make_param_steps.py, copies of them
# in results/params_steps/ whose values are read from params.yaml, and runs those.
#
#   bash workflow/run_all.sh              # Steps 1, 2 and 3
#   bash workflow/run_all.sh "2 3"        # only some steps
#   bash workflow/run_all.sh "0 1 2 3"    # also Step 0 (Cell Ranger)
#
# Run it inside the container, with the repository mounted at /workspace.
# Step 0 is NOT part of the default: it needs a signed, account-tied 10x
# download URL and Cell Ranger on PATH (see step0_cellranger.sh), and the
# Cell Ranger matrices are already in the repository.
# step1_supplementary.R is not run here: it continues the R session of Step 1.
# Each step writes its log to results/logs/.
set -euo pipefail

STEPS="${1:-1 2 3}"
cd "$(dirname "$0")/.."
[ "$(pwd)" = "/workspace" ] || echo "WARNING: repository is at $(pwd), the scripts expect /workspace" >&2
export PARAMS_FILE="${PARAMS_FILE:-$(pwd)/params.yaml}"
[ -f "$PARAMS_FILE" ] || { echo "params file not found: $PARAMS_FILE" >&2; exit 1; }
mkdir -p results/logs
python3 workflow/make_param_steps.py results/params_steps

run() {   # run <step> <command...>
  local step="$1"; shift
  echo "=== Step ${step}: $* ==="
  "$@" 2>&1 | tee "results/logs/step${step}.log"
  echo "=== Step ${step} done ==="
}

for s in $STEPS; do
  case "$s" in
    0) run 0 bash    workflow/step0_cellranger.sh ;;
    1) run 1 Rscript results/params_steps/step1_singlecell.R ;;
    2) run 2 Rscript results/params_steps/step2_degs.R ;;
    3) run 3 python3 results/params_steps/step3_pseudotime.py ;;
    *) echo "Unknown step: $s (use 0, 1, 2 or 3)" >&2; exit 1 ;;
  esac
done
echo "All requested steps finished."
