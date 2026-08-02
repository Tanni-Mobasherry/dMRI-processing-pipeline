#!/bin/bash

# ============================================================
# Step 9: DWI Bias-Field Correction
# Test subject: YTH001 / BL
#
# Inputs:
#   dwi_eddy.nii.gz
#   dwi_eddy.eddy_rotated_bvecs
#   dwi_cat.bval
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/eddy/"
#
# Outputs:
#   dwi_eddy.mif
#   dwi_eddy_BA.mif
#   biasfield.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/"
# ============================================================


# ------------------------------------------------------------
# Step 1: Convert the eddy-corrected DWI from NIfTI to MRtrix
#          format and import the gradient information
# ------------------------------------------------------------

mrconvert \
dwi_eddy.nii.gz \
dwi_eddy.mif \
-fslgrad \
dwi_eddy.eddy_rotated_bvecs \
dwi_cat.bval \
-force


# ------------------------------------------------------------
# Step 2: Perform ANTs N4 bias-field correction on the
#          eddy-corrected DWI
# ------------------------------------------------------------

dwibiascorrect ants \
dwi_eddy.mif \
dwi_eddy_BA.mif \
-bias biasfield.mif \
-force
