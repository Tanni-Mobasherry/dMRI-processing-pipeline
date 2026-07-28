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
#   dwi_eddy_BA.mif
#   T1.mgz (longitudinal FreeSurfer)
#   parcels_a2009s.mif
#
# Final output:
#   parcels_a2009s_dwi.mif
# ================================================================

ROOT="/Volumes/Toshiba-Ext/raw-data/YTH001/BL"

WORKDIR="$ROOT/dmri/anatomy-parcellation"

DWI_INPUT="$ROOT/dmri/preproc/bias-field-correction/dwi_eddy_BA.mif"

T1_MGZ="/Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/T1.mgz"

PARCELLATION_MIF="$WORKDIR/parcels_a2009s.mif"

# ------------------------------------------------------------
# Step 1: Create mean bias-corrected b0 image
# ------------------------------------------------------------

dwiextract \
    "$DWI_INPUT" \
    - \
    -bzero \
    -quiet | \
mrmath \
    - \
    mean \
    "$WORKDIR/mean_b0_BA.mif" \
    -axis 3 \
    -force \
    -quiet

# ------------------------------------------------------------
# Step 2: Convert mean b0 to NIfTI
# ------------------------------------------------------------

mrconvert \
    "$WORKDIR/mean_b0_BA.mif" \
    "$WORKDIR/mean_b0_BA.nii.gz" \
    -force \
    -quiet

# ------------------------------------------------------------
# Step 3: Convert longitudinal FreeSurfer T1 to NIfTI
# ------------------------------------------------------------

mrconvert \
    "$T1_MGZ" \
    "$WORKDIR/T1_FS_longi.nii.gz" \
    -force \
    -quiet

# ------------------------------------------------------------
# Step 4: Affine registration
#
# Reference image:
#   mean_b0_BA.nii.gz
#
# Floating image:
#   T1_FS_longi.nii.gz
# ------------------------------------------------------------

reg_aladin \
    -ref "$WORKDIR/mean_b0_BA.nii.gz" \
    -flo "$WORKDIR/T1_FS_longi.nii.gz" \
    -aff "$WORKDIR/aff.txt" \
    -res "$WORKDIR/t1_aff.nii.gz" \
    > "$WORKDIR/reg_aladin_console.log" 2>&1

# ------------------------------------------------------------
# Step 5: Nonlinear registration
# ------------------------------------------------------------

reg_f3d \
    -ref "$WORKDIR/mean_b0_BA.nii.gz" \
    -flo "$WORKDIR/T1_FS_longi.nii.gz" \
    -aff "$WORKDIR/aff.txt" \
    -cpp "$WORKDIR/cpp.nii.gz" \
    -res "$WORKDIR/t1_warp.nii.gz" \
    > "$WORKDIR/reg_f3d_console.log" 2>&1

# ------------------------------------------------------------
# Step 6: Convert parcellation to NIfTI
# ------------------------------------------------------------

mrconvert \
    "$PARCELLATION_MIF" \
    "$WORKDIR/parcels_a2009s.nii.gz" \
    -force \
    -quiet

# ------------------------------------------------------------
# Step 7: Resample parcellation into diffusion space
#
# Nearest-neighbour interpolation (-inter 0) is used because
# the parcellation contains discrete integer labels.
# ------------------------------------------------------------

reg_resample \
    -ref "$WORKDIR/mean_b0_BA.nii.gz" \
    -flo "$WORKDIR/parcels_a2009s.nii.gz" \
    -trans "$WORKDIR/cpp.nii.gz" \
    -inter 0 \
    -res "$WORKDIR/parcels_a2009s_dwi.nii.gz" \
    > "$WORKDIR/reg_resample_console.log" 2>&1

# ------------------------------------------------------------
# Step 8: Convert back to MRtrix format
# ------------------------------------------------------------

mrconvert \
    "$WORKDIR/parcels_a2009s_dwi.nii.gz" \
    "$WORKDIR/parcels_a2009s_dwi.mif" \
    -force \
    -quiet
