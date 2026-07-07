#!/bin/bash
LOG=/workspace/pipeline_v2.log
echo "[V2] Starting $(date)" | tee $LOG
Rscript /workspace/cap1_v2.R >> $LOG 2>&1 && echo "[V2] Cap1 DONE $(date)" | tee -a $LOG || { echo "[V2] Cap1 ERROR" | tee -a $LOG; exit 1; }
Rscript /workspace/cap2_v2.R >> $LOG 2>&1 && echo "[V2] Cap2 DONE $(date)" | tee -a $LOG || { echo "[V2] Cap2 ERROR" | tee -a $LOG; exit 1; }
python3 /workspace/cap3_v2.py >> $LOG 2>&1 && echo "[V2] Cap3 DONE $(date)" | tee -a $LOG || { echo "[V2] Cap3 ERROR" | tee -a $LOG; exit 1; }
echo "[V2] PIPELINE COMPLETE" | tee -a $LOG
