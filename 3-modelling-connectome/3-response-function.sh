
# ================================================================
# Response-Function Estimation
#
# Purpose:
# Estimate tissue-specific response functions for white matter,
# grey matter, and cerebrospinal fluid using the unsupervised
# Dhollander algorithm for subsequent MSMT-CSD.
#
# Subject:
#   YTH001_BL
#
# Input:
#
#   dwi_eddy_BA.mif
#   Bias-field-corrected diffusion image generated in:
#   1-diffusion-preprocessing/09_bias_field_correction.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/"
#
# Outputs:
#
#   wm.txt
#   gm.txt
#   csf.txt
#   rf_voxels.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/response-function/"
# ================================================================

# ------------------------------------------------------------
# Step 1: Estimate WM, GM and CSF response functions
# ------------------------------------------------------------

dwi2response dhollander \
    dwi_eddy_BA.mif \
    wm.txt \
    gm.txt \
    csf.txt \
    -voxels rf_voxels.mif \
    -force
