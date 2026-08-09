# ================================================================
# DWI Upsampling 
#
# Purpose:
# Upsample the bias-field-corrected DWI from 2 mm to
# 1.25 mm isotropic spatial resolution before FOD estimation,
# and generate a corresponding brain mask from the upsampled DWI.
#
# Subject:
# YTH001_BL
#
# Input:
#
# dwi_eddy_BA.mif
# Bias-field-corrected diffusion image generated in:
# 1-diffusion-preprocessing/09_bias_field_correction.sh
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/"
#
# Outputs:
#
# dwi_eddy_BA_upsampled.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/"
# ------------------------------------------------------------
# Step 1: Upsample the bias-field-corrected DWI
# ------------------------------------------------------------
#
# Resample the spatial voxel size:
#
# 2 × 2 × 2 mm
# →
# 1.25 × 1.25 × 1.25 mm
#
# The diffusion dimension and diffusion-gradient information
# are retained.
# ------------------------------------------------------------

mrgrid \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/dwi_eddy_BA.mif \
regrid \
-vox 1.25 \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif

# ------------------------------------------------------------
# Step 2: Check upsampled DWI geometry
# ------------------------------------------------------------
#
# Expected input:
# Dimensions: 104 × 104 × 72 × 198
# Voxel size: 2 × 2 × 2 × 3.5 mm
#
# Expected output:
# Dimensions: 166 × 166 × 115 × 198
# Voxel size: 1.25 × 1.25 × 1.25 × 3.5 mm
# ------------------------------------------------------------

mrinfo \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-size \
-spacing


# ------------------------------------------------------------
# Step 3: Check diffusion shells and gradient count
# ------------------------------------------------------------
#
# Expected shells:
# 0, 1500 and 3000 s/mm²
#
# Expected number of gradient rows:
# 198
# ------------------------------------------------------------

mrinfo \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-shell_bvalues

mrinfo \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-dwgrad | wc -l

# ------------------------------------------------------------
# Step 4: Check DWI mean intensity before and after upsampling
# ------------------------------------------------------------
# Upsampling changes the voxel grid through spatial interpolation.
# With the same physical field of view and comparable brain mask, the
# mean or median signal in each diffusion volume should remain broadly
# similar, but exact equality is not expected.
#
# Large changes may indicate an FOV/mask/resampling problem or a globally
# abnormal volume. Combine this check with finite-voxel counts, slice-wise
# coverage, robust intensity statistics, and visual QC.
# The commands below display the first 10 volumes as a quick
# inspection. Full content QC across all 198 volumes is
# performed separately using:
#
# upsampling_mean_intensity_QC.sh
# ------------------------------------------------------------

echo "=== BEFORE UPSAMPLING: FIRST 10 VOLUMES ==="

mrstats \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/dwi_eddy_BA.mif \
-output min \
-output max \
-output mean | sed -n '1,10p'

echo "=== AFTER UPSAMPLING: FIRST 10 VOLUMES ==="

mrstats \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-output min \
-output max \
-output mean | sed -n '1,10p'
# ------------------------------------------------------------



