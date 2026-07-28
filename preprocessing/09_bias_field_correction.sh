#!/bin/bash

# DWI Bias-Field Correction
# Example dataset: YTH001/BL
# Paths below are relative to the YTH001 directory.
# Run this script from inside the YTH001 directory.

DATA="BL"

# 1. EDDY Output -> MRtrix Conversion
mrconvert \
  "$DATA/dmri/preproc/eddy/dwi_eddy.nii.gz" \
  "$DATA/dmri/preproc/bias-field-correction/dwi_eddy.mif" \
  -fslgrad \
  "$DATA/dmri/preproc/eddy/dwi_eddy.eddy_rotated_bvecs" \
  "$DATA/dmri/preproc/eddy/dwi_cat.bval" \
  -force

# 2. DWI Bias-Field Correction — ANTs N4
dwibiascorrect ants \
  "$DATA/dmri/preproc/bias-field-correction/dwi_eddy.mif" \
  "$DATA/dmri/preproc/bias-field-correction/dwi_eddy_BA.mif" \
  -bias "$DATA/dmri/preproc/bias-field-correction/biasfield.mif" \
  -force
