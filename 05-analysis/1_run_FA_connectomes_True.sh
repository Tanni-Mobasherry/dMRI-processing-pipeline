#!/bin/bash

set -u
set -o pipefail

# ============================================================
# FA-WEIGHTED CONNECTOMES
# ============================================================
#
# CHANGE vs the original version of this script
# ---------------------------------------------
# The original passed the NATIVE-space FA maps
#     <SID>/{BL,FU}/dmri/preprocessing/dti-metrics/fa.mif
# directly to tcksample, while the streamlines live in the longitudinal
# FOD-template space. Those are different spaces:
#     fa.mif        104 x 104 x  72 @ 2.00 mm, rotated ~7 deg
#     fod_template  141 x 158 x 136 @ 1.25 mm
# tcksample does not error on a space mismatch - it samples wherever the
# streamline coordinates land - so the run completed cleanly while
# returning anatomically wrong values. The symptom: the BL and FU
# connectomes, built from IDENTICAL streamlines, correlated at only
# r = 0.20 across subjects (15/36 negative), where FBC on the same
# streamlines gives r = 0.999.
#
# Fix: warp each session's FA into template space first, using the rigid
# session->template transforms that population_template already produced
# (template/transforms/{BL,FU}.txt). Verified on YTH001: BL-FU connectome
# correlation goes from -0.329 to +0.988.
#
# A hard guard below aborts the subject if the warped FA grid does not
# match the template, so this cannot silently recur.
# ============================================================

HARD="/Volumes/Toshiba-Ext/raw-data"
START="$HOME/Desktop/start-analysis"
CODE="$START/FA-connectomes"          # scripts, logs, QC live here
OUT="$START/FA"                        # results, one folder per participant
SUBJECT_FILE="$CODE/subjects.txt"
LOG="$CODE/FA_connectome_batch.log"
FAIL="$CODE/FA_connectome_failures.txt"

mkdir -p "$OUT"
: > "$LOG"
: > "$FAIL"

if [[ ! -f "$SUBJECT_FILE" ]]; then
    echo "ERROR: Missing subject list: $SUBJECT_FILE"
    exit 1
fi

N_SUBJECTS=$(grep -v '^[[:space:]]*$' "$SUBJECT_FILE" | wc -l | tr -d ' ')

echo "============================================================" | tee -a "$LOG"
echo "FA CONNECTOME BATCH (template-space corrected)"              | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "Subjects: $N_SUBJECTS"                                       | tee -a "$LOG"
echo "Input:    $HARD"                                             | tee -a "$LOG"
echo "Output:   $OUT"                                              | tee -a "$LOG"
echo "Started:  $(date)"                                           | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

CURRENT=0

while IFS= read -r SID || [[ -n "$SID" ]]; do

    SID=$(echo "$SID" | tr -d '\r' | xargs)
    [[ -z "$SID" ]] && continue
    CURRENT=$((CURRENT + 1))

    echo                                                            | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"
    echo "[$CURRENT/$N_SUBJECTS] $SID"                              | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"

    TDIR="$HARD/$SID/longi_FT/template"
    TRACKS="$TDIR/tracks_10M.tck"
    NODES="$TDIR/nodes.mif"
    FOD_TEMPLATE="$TDIR/fod_template.mif"

    SUBJECT_OUT="$OUT/$SID"
    WEIGHTS="$SUBJECT_OUT/weights"
    CONNECTOMES="$SUBJECT_OUT"          # connectomes sit at the top level,
                                        # mirroring start-analysis/FBC/<ID>/
    FA_TEMPL="$SUBJECT_OUT/fa_template_space"
    mkdir -p "$WEIGHTS" "$CONNECTOMES" "$FA_TEMPL"

    # --------------------------------------------------------
    # INPUT VALIDATION
    # --------------------------------------------------------
    MISSING=0
    for FILE in "$TRACKS" "$NODES" "$FOD_TEMPLATE" \
                "$TDIR/transforms/BL.txt" "$TDIR/transforms/FU.txt" \
                "$HARD/$SID/BL/dmri/preprocessing/dti-metrics/fa.mif" \
                "$HARD/$SID/FU/dmri/preprocessing/dti-metrics/fa.mif"; do
        if [[ ! -f "$FILE" ]]; then
            echo "ERROR: Missing input: $FILE" | tee -a "$LOG"
            MISSING=1
        fi
    done
    if [[ "$MISSING" -eq 1 ]]; then
        echo "$SID : missing input" >> "$FAIL"
        echo "SKIPPING $SID" | tee -a "$LOG"
        continue
    fi

    REF_SIZE=$(mrinfo "$FOD_TEMPLATE" -size | cut -d' ' -f1-3 | tr -s ' ')

    SUBJECT_FAILED=0

    for TP in BL FU; do

        FA_NATIVE="$HARD/$SID/$TP/dmri/preprocessing/dti-metrics/fa.mif"
        FA_IN_TEMPLATE="$FA_TEMPL/${TP}_fa_template_space.mif"

        # ----------------------------------------------------
        # 1. WARP FA INTO TEMPLATE SPACE  (the corrected step)
        # ----------------------------------------------------
        echo "[$TP] mrtransform FA -> template space..." | tee -a "$LOG"

        if ! mrtransform "$FA_NATIVE" \
            -linear "$TDIR/transforms/$TP.txt" \
            -template "$FOD_TEMPLATE" \
            "$FA_IN_TEMPLATE" \
            -force \
            >> "$LOG" 2>&1
        then
            echo "$SID : $TP mrtransform failed" >> "$FAIL"
            echo "ERROR: $TP mrtransform failed." | tee -a "$LOG"
            SUBJECT_FAILED=1; break
        fi

        # ----------------------------------------------------
        # 2. GUARD: warped FA must sit on the template grid
        # ----------------------------------------------------
        GOT_SIZE=$(mrinfo "$FA_IN_TEMPLATE" -size | cut -d' ' -f1-3 | tr -s ' ')
        if [[ "$GOT_SIZE" != "$REF_SIZE" ]]; then
            echo "$SID : $TP FA grid mismatch ($GOT_SIZE vs $REF_SIZE)" >> "$FAIL"
            echo "ERROR: $TP warped FA grid [$GOT_SIZE] != template [$REF_SIZE]" | tee -a "$LOG"
            SUBJECT_FAILED=1; break
        fi
        echo "  grid OK: $GOT_SIZE" | tee -a "$LOG"

        # ----------------------------------------------------
        # 3. SAMPLE ALONG STREAMLINES  (supervisor's recipe)
        # ----------------------------------------------------
        echo "[$TP] tcksample..." | tee -a "$LOG"
        if ! tcksample "$TRACKS" "$FA_IN_TEMPLATE" "$WEIGHTS/${TP}_FA.csv" \
            -stat_tck mean -force >> "$LOG" 2>&1
        then
            echo "$SID : $TP tcksample failed" >> "$FAIL"
            echo "ERROR: $TP tcksample failed." | tee -a "$LOG"
            SUBJECT_FAILED=1; break
        fi

        echo "[$TP] tck2connectome..." | tee -a "$LOG"
        if ! tck2connectome "$TRACKS" "$NODES" "$CONNECTOMES/${TP}_FA.csv" \
            -scale_file "$WEIGHTS/${TP}_FA.csv" -stat_edge mean -force \
            >> "$LOG" 2>&1
        then
            echo "$SID : $TP tck2connectome failed" >> "$FAIL"
            echo "ERROR: $TP tck2connectome failed." | tee -a "$LOG"
            SUBJECT_FAILED=1; break
        fi

    done

    [[ "$SUBJECT_FAILED" -eq 1 ]] && continue

    # --------------------------------------------------------
    # 4. OUTPUT CHECKS
    # --------------------------------------------------------
    if [[ ! -s "$CONNECTOMES/BL_FA.csv" || ! -s "$CONNECTOMES/FU_FA.csv" ]]; then
        echo "$SID : empty connectome output" >> "$FAIL"
        echo "ERROR: Empty connectome output." | tee -a "$LOG"
        continue
    fi

    BL_ROWS=$(wc -l < "$CONNECTOMES/BL_FA.csv" | tr -d ' ')
    FU_ROWS=$(wc -l < "$CONNECTOMES/FU_FA.csv" | tr -d ' ')
    echo "BL matrix rows: $BL_ROWS | FU matrix rows: $FU_ROWS" | tee -a "$LOG"

    if [[ "$BL_ROWS" -ne 164 || "$FU_ROWS" -ne 164 ]]; then
        echo "$SID : unexpected matrix dimensions" >> "$FAIL"
        echo "WARNING: Expected 164 rows." | tee -a "$LOG"
        continue
    fi

    # --------------------------------------------------------
    # 5. GUARD: BL and FU come from identical streamlines, so the
    #    two connectomes must be strongly correlated. This is the
    #    check that would have caught the original space bug.
    # --------------------------------------------------------
    R=$(python3 -c "
import numpy as np
iu = np.triu_indices(164, 1)
b = np.loadtxt('$CONNECTOMES/BL_FA.csv', delimiter=',')[iu]
f = np.loadtxt('$CONNECTOMES/FU_FA.csv', delimiter=',')[iu]
m = (b > 0) & (f > 0)
print(f'{np.corrcoef(b[m], f[m])[0,1]:.4f}')
" 2>/dev/null)

    echo "BL-FU connectome correlation: $R" | tee -a "$LOG"

    if awk "BEGIN{exit !($R < 0.9)}"; then
        echo "$SID : LOW BL-FU correlation ($R) - check registration" >> "$FAIL"
        echo "WARNING: BL-FU correlation $R is below 0.9." | tee -a "$LOG"
        continue
    fi

    echo "OK $SID COMPLETE (r = $R)" | tee -a "$LOG"

done < "$SUBJECT_FILE"

echo                                                                | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "BATCH FINISHED"                                               | tee -a "$LOG"
echo "Finished: $(date)"                                            | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

if [[ -s "$FAIL" ]]; then
    echo; echo "Subjects requiring attention:" | tee -a "$LOG"
    cat "$FAIL" | tee -a "$LOG"
else
    echo "No failures recorded." | tee -a "$LOG"
fi

echo
echo "Next: re-run the labelling step, then"
echo "  python3 $CODE/qc/check_FA_connectome_validity.py"
echo "which should now report FA r > 0.9 and VERDICT: PASS."
