#!/bin/bash
LOG=/workspace/pipeline_v1.log
echo "[V1] Starting $(date)" | tee $LOG
Rscript /workspace/cap1_v1.R >> $LOG 2>&1 && echo "[V1] Cap1 DONE $(date)" | tee -a $LOG || { echo "[V1] Cap1 ERROR" | tee -a $LOG; exit 1; }
Rscript /workspace/cap2_v1.R >> $LOG 2>&1 && echo "[V1] Cap2 DONE $(date)" | tee -a $LOG || { echo "[V1] Cap2 ERROR" | tee -a $LOG; exit 1; }
python3 /workspace/cap3_v1.py >> $LOG 2>&1 && echo "[V1] Cap3 DONE $(date)" | tee -a $LOG || { echo "[V1] Cap3 ERROR" | tee -a $LOG; exit 1; }
echo "[V1] PIPELINE COMPLETE" | tee -a $LOG
