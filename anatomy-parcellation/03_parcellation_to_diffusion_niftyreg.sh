
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
#   brain.mgz (longitudinal FreeSurfer) "/Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/brain.mgz"
#   parcels_a2009s.mif "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation"
#
# Final output:
#   parcels_a2009s_dwi.mif
# ================================================================

# ------------------------------------------------------------
# Step 1: Create mean bias-corrected b0 image
# ------------------------------------------------------------
dwiextract dwi_eddy_BA.mif - -bzero -quiet | \
mrmath - mean mean_b0_BA.mif -axis 3 -force -quiet

# ------------------------------------------------------------
#Step 2: Create registration masks
# ------------------------------------------------------------
dwi2mask dwi_eddy_BA.mif dwi_mask.mif -force

# ------------------------------------------------------------
# Step 3: Convert inputs to NIfTI
# ------------------------------------------------------------

# Mean b0 (brain)
mrconvert \
    mean_b0_BA.mif \
    mean_b0_BA.nii.gz \
    -force \
    -quiet

# Longitudinal FreeSurfer brain
mrconvert \
    brain.mgz \
    T1_brain.nii.gz \
    -force \
    -quiet

# Parcellation
mrconvert \
    parcels_a2009s.mif \
    parcels_a2009s.nii.gz \
    -force \
    -quiet
#refrence mask
mrconvert \
    dwi_mask.mif \
    dwi_mask.nii.gz \
    -datatype uint8 \
    -force \
    -quiet

# ------------------------------------------------------------
# Step 4: Affine T1-to-b0 registration + (reference mask applied)
# ------------------------------------------------------------
reg_aladin \
   -ref mean_b0_BA.nii.gz \
   -flo T1_brain.nii.gz \
   -rmask dwi_mask.nii.gz \
   -aff T1_to_b0_aff.txt \
   -res T1_in_b0_aff.nii.gz \
   > reg_aladin_console.log 2>&1
# ------------------------------------------------------------
# Step 5: Nonlinear registration + (reference mask applied)
# ------------------------------------------------------------
reg_f3d \
    -ref mean_b0_BA.nii.gz \
    -flo T1_brain.nii.gz \
    -rmask dwi_mask.nii.gz \
    -aff T1_to_b0_aff.txt \
    -cpp T1_to_b0_cpp.nii.gz \
    -res T1_in_b0_nonlinear_t1wrap.nii.gz \
    > reg_f3d_console.log 2>&1

# ---------------------------------------------------------------------------------
# Step 6: Warp the parcellation to diffusion space using the nonlinear transform
# 
# Nearest-neighbour interpolation (-inter 0) is used because
# the parcellation contains discrete integer labels.
# ----------------------------------------------------------------------------------
reg_resample \
    -ref mean_b0_BA.nii.gz \
    -flo parcels_a2009s.nii.gz \
    -trans T1_to_b0_cpp.nii.gz \
    -inter 0 \
    -res parcels_a2009s_dwi.nii.gz \
    > reg_resample_console.log 2>&1

# ------------------------------------------------------------
# Step 7: Convert output to MRtrix format
# ------------------------------------------------------------
mrconvert \
parcels_a2009s_dwi.nii.gz \
parcels_a2009s_dwi.mif \
-force \
-quiet
