#!/usr/bin/env bash
# ================================================================
# YTH001_BL: Parcellation → Diffusion Space (NiftyReg)
#
# Purpose:
# Register the longitudinal FreeSurfer T1 image to diffusion
# space and transform the cortical parcellation into diffusion
# space.
#
# Subject:
#   YTH001_BL
#
# Inputs:
#   dwi_eddy_BA.mif/ "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preproc/bias-field-correction/dwi_eddy_BA.mif"
#   Brain.mgz (longitudinal FreeSurfer) "/Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/Brain.mgz"
#   parcels_a2009s.mif "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation"
#
# Final output:
#   parcels_a2009s_dwi.mif
# ================================================================

# ------------------------------------------------------------
# Step 1: Create mean bias-corrected b0 image
# ------------------------------------------------------------

dwiextract \
dwi_eddy_BA.mif \
- \
-bzero \
-quiet | \
mrmath \
- \
mean \
mean_b0_BA.mif \
-axis 3 \
-force \
-quiet
# ------------------------------------------------------------
# Step 2: Convert mean b0 to NIfTI
# ------------------------------------------------------------

mrconvert \
    mean_b0_BA.mif \
    mean_b0_BA.nii.gz\
    -force \
    -quiet

# ------------------------------------------------------------
# Step 3: Convert longitudinal FreeSurfer T1 (brain.mgz) to NIfTI
# ------------------------------------------------------------
mrconvert \
brain.mgz \
T1_brain.nii.gz \
-force \
-quiet
# ------------------------------------------------------------
# Step 4: Affine registration
# ------------------------------------------------------------
reg_aladin \
-ref mean_b0_BA.nii.gz \
-flo T1_brain.nii.gz \
-aff aff.txt \
-res t1_aff.nii.gz \
> reg_aladin_console.log 2>&1

# ------------------------------------------------------------
# Step 5: Nonlinear registration
# ------------------------------------------------------------

reg_f3d \
-ref mean_b0_BA.nii.gz \
-flo T1_brain.nii.gz \
-aff aff.txt \
-cpp cpp.nii.gz \
-res t1_warp.nii.gz \
> reg_f3d_console.log 2>&1

# ------------------------------------------------------------
# Step 6: Convert parcellation to NIfTI
# ------------------------------------------------------------
mrconvert \
parcels_a2009s.mif \
parcels_a2009s.nii.gz \
-force \
-quiet
# ------------------------------------------------------------
# Step 7: Resample parcellation into diffusion space
#
# Nearest-neighbour interpolation (-inter 0) is used because
# the parcellation contains discrete integer labels.
# ------------------------------------------------------------
reg_resample \
-ref mean_b0_BA.nii.gz \
-flo parcels_a2009s.nii.gz \
-trans cpp.nii.gz \
-inter 0 \
-res parcels_a2009s_dwi.nii.gz \
> reg_resample_console.log 2>&1

# ------------------------------------------------------------
# Step 8: Convert back to MRtrix format
# ------------------------------------------------------------
mrconvert \
parcels_a2009s_dwi.nii.gz \
parcels_a2009s_dwi.mif \
-force \
-quiet
