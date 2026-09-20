#!/bin/bash

set -u
set -o pipefail

# ============================================================
# PATHS
# ============================================================

HARD="/Volumes/Toshiba-Ext/raw-data"
OUT="$HOME/Desktop/start-analysis/FA-connectomes"
SUBJECT_FILE="$OUT/subjects.txt"
LOG="$OUT/FA_connectome_batch.log"
FAIL="$OUT/FA_connectome_failures.txt"

mkdir -p "$OUT"

# Start fresh logs
: > "$LOG"
: > "$FAIL"

# ============================================================
# CHECK SUBJECT FILE
# ============================================================

if [[ ! -f "$SUBJECT_FILE" ]]; then
    echo "ERROR: Missing subject list:"
    echo "$SUBJECT_FILE"
    exit 1
fi

mapfile_cmd_available=false

# Count non-empty lines
N_SUBJECTS=$(grep -v '^[[:space:]]*$' "$SUBJECT_FILE" | wc -l | tr -d ' ')

echo "============================================================" | tee -a "$LOG"
echo "FA CONNECTOME BATCH"                                         | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "Subjects: $N_SUBJECTS"                                      | tee -a "$LOG"
echo "Input:    $HARD"                                            | tee -a "$LOG"
echo "Output:   $OUT"                                             | tee -a "$LOG"
echo "Started:  $(date)"                                          | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo

# ============================================================
# PROCESS SUBJECTS
# ============================================================

CURRENT=0

while IFS= read -r SID || [[ -n "$SID" ]]; do

    # Remove whitespace / CR characters
    SID=$(echo "$SID" | tr -d '\r' | xargs)

    # Skip blank lines
    [[ -z "$SID" ]] && continue

    CURRENT=$((CURRENT + 1))

    echo
    echo "============================================================" | tee -a "$LOG"
    echo "[$CURRENT/$N_SUBJECTS] $SID"                                | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"

    TRACKS="$HARD/$SID/longi_FT/template/tracks_10M.tck"
    NODES="$HARD/$SID/longi_FT/template/nodes.mif"

    BL_FA="$HARD/$SID/BL/dmri/preprocessing/dti-metrics/fa.mif"
    FU_FA="$HARD/$SID/FU/dmri/preprocessing/dti-metrics/fa.mif"

    SUBJECT_OUT="$OUT/$SID"
    WEIGHTS="$SUBJECT_OUT/weights"
    CONNECTOMES="$SUBJECT_OUT/connectomes"

    mkdir -p "$WEIGHTS" "$CONNECTOMES"

    # --------------------------------------------------------
    # INPUT VALIDATION
    # --------------------------------------------------------

    MISSING=0

    for FILE in "$TRACKS" "$NODES" "$BL_FA" "$FU_FA"; do

        if [[ ! -f "$FILE" ]]; then
            echo "ERROR: Missing input:" | tee -a "$LOG"
            echo "  $FILE"             | tee -a "$LOG"
            MISSING=1
        fi

    done

    if [[ "$MISSING" -eq 1 ]]; then
        echo "$SID : missing input" >> "$FAIL"
        echo "SKIPPING $SID" | tee -a "$LOG"
        continue
    fi

    # ========================================================
    # BASELINE
    # ========================================================

    echo "[BL] tcksample..." | tee -a "$LOG"

    if ! tcksample \
        "$TRACKS" \
        "$BL_FA" \
        "$WEIGHTS/BL_FA.csv" \
        -stat_tck mean \
        -force \
        >> "$LOG" 2>&1
    then
        echo "$SID : BL tcksample failed" >> "$FAIL"
        echo "ERROR: BL tcksample failed." | tee -a "$LOG"
        continue
    fi

    echo "[BL] tck2connectome..." | tee -a "$LOG"

    if ! tck2connectome \
        "$TRACKS" \
        "$NODES" \
        "$CONNECTOMES/BL_FA.csv" \
        -scale_file "$WEIGHTS/BL_FA.csv" \
        -stat_edge mean \
        -force \
        >> "$LOG" 2>&1
    then
        echo "$SID : BL tck2connectome failed" >> "$FAIL"
        echo "ERROR: BL tck2connectome failed." | tee -a "$LOG"
        continue
    fi

    # ========================================================
    # FOLLOW-UP
    # ========================================================

    echo "[FU] tcksample..." | tee -a "$LOG"

    if ! tcksample \
        "$TRACKS" \
        "$FU_FA" \
        "$WEIGHTS/FU_FA.csv" \
        -stat_tck mean \
        -force \
        >> "$LOG" 2>&1
    then
        echo "$SID : FU tcksample failed" >> "$FAIL"
        echo "ERROR: FU tcksample failed." | tee -a "$LOG"
        continue
    fi

    echo "[FU] tck2connectome..." | tee -a "$LOG"

    if ! tck2connectome \
        "$TRACKS" \
        "$NODES" \
        "$CONNECTOMES/FU_FA.csv" \
        -scale_file "$WEIGHTS/FU_FA.csv" \
        -stat_edge mean \
        -force \
        >> "$LOG" 2>&1
    then
        echo "$SID : FU tck2connectome failed" >> "$FAIL"
        echo "ERROR: FU tck2connectome failed." | tee -a "$LOG"
        continue
    fi

    # ========================================================
    # BASIC OUTPUT CHECK
    # ========================================================

    if [[ ! -s "$CONNECTOMES/BL_FA.csv" || ! -s "$CONNECTOMES/FU_FA.csv" ]]; then
        echo "$SID : empty connectome output" >> "$FAIL"
        echo "ERROR: Empty connectome output." | tee -a "$LOG"
        continue
    fi

    BL_ROWS=$(wc -l < "$CONNECTOMES/BL_FA.csv" | tr -d ' ')
    FU_ROWS=$(wc -l < "$CONNECTOMES/FU_FA.csv" | tr -d ' ')

    echo "BL matrix rows: $BL_ROWS" | tee -a "$LOG"
    echo "FU matrix rows: $FU_ROWS" | tee -a "$LOG"

    if [[ "$BL_ROWS" -ne 164 || "$FU_ROWS" -ne 164 ]]; then
        echo "$SID : unexpected matrix dimensions" >> "$FAIL"
        echo "WARNING: Expected 164 rows." | tee -a "$LOG"
        continue
    fi

    echo "✓ $SID COMPLETE" | tee -a "$LOG"

done < "$SUBJECT_FILE"


# ============================================================
# FINAL REPORT
# ============================================================

echo
echo "============================================================" | tee -a "$LOG"
echo "BATCH FINISHED"                                               | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "Finished: $(date)"                                          | tee -a "$LOG"

if [[ -s "$FAIL" ]]; then

    echo
    echo "Subjects requiring attention:" | tee -a "$LOG"
    cat "$FAIL" | tee -a "$LOG"

else

    echo "No failures recorded." | tee -a "$LOG"

fi

echo
echo "Outputs:"
echo "$OUT"
echo
echo "Log:"
echo "$LOG"
echo "============================================================"
