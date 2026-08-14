# ================================================================

# Upsampled DWI Brain-Mask Generation and QC
# Purpose:
# Generate a brain mask from the 1.25 mm upsampled,
# bias-field-corrected DWI and perform geometry, integrity, volume-intensity, and visual quality-control checks.

# This step identifies mask holes or excluded brain regions but does not modify the mask. Any required mask correction is performed in a separate step.

# Subject:

# YTH001_BL

# Input:

# dwi_eddy_BA_upsampled.mif

# Upsampled bias-field-corrected diffusion image generated in:

# 03-modelling-connectome/05_dwi_upsampling.md

# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/"

# Primary output:

# dwi_mask_upsampled.mif

# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/"

# QC outputs:

# YTH001_BL_masked_volume_statistics.csv

# YTH001_BL_mean_b0.mif

# YTH001_BL_mean_b1500.mif

# YTH001_BL_mean_b3000.mif

#YTH001_BL_upsampled_mask_QC/"

# ================================================================


# ------------------------------------------------------------
# Step 1: Generate the brain mask from the upsampled DWI
# ------------------------------------------------------------

# The default dwi2mask procedure is used to generate a binary brain mask directly from the 1.25 mm upsampled DWI.
# ------------------------------------------------------------

dwi2mask \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-force


# ------------------------------------------------------------
# Step 2: Check the brain-mask geometry
# ------------------------------------------------------------
# Expected dimensions:
# 166 × 166 × 115
# Expected voxel size:
# 1.25 × 1.25 × 1.25 mm
# The mask geometry should correspond to the first three spatial
# dimensions of the upsampled DWI.
# ------------------------------------------------------------

mrinfo \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-size \
-spacing


# ------------------------------------------------------------
# Step 3: Check that the mask is binary and non-empty
# ------------------------------------------------------------
# Expected intensity range:
# Minimum = 0
# Maximum = 1
# The mean represents the proportion of the complete image grid
# included within the mask. It is recorded as a screening metric and does not replace visual inspection.

# ------------------------------------------------------------

mrstats \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-output min \
-output max \
-output mean


# ------------------------------------------------------------
# Step 4: Check for NaN or infinite values
# ------------------------------------------------------------
# The finite operator identifies whether each DWI voxel contains a finite numerical value.
# The resulting binary values are counted across all 198 volumes.
# Expected result:
# Non-finite voxels: 0

# ------------------------------------------------------------

mrcalc \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-finite \
0 \
-eq \
- |
mrstats \
- \
-ignorezero \
-output count |
awk '
{
    total += $1
}
END {
    print "Non-finite voxels:", total+0
}'

# ------------------------------------------------------------
# Step 5: Calculate masked statistics for all 198 DWI volumes
# ------------------------------------------------------------
# Calculate the minimum, maximum, mean, and median intensity
# within the upsampled brain mask for every diffusion volume.
# The first column records the zero-based volume number.
# Output:
# YTH001_BL_masked_volume_statistics.csv
# -----------------------------------------------------------

mrstats \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-mask \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-output min \
-output max \
-output mean \
-output median |
awk '
BEGIN {
    print "volume,min,max,mean,median"
}
{
    print NR-1 "," $1 "," $2 "," $3 "," $4
}' > \
YTH001_BL_upsampled_mask_QC/YTH001_BL_masked_volume_statistics.csv

# ------------------------------------------------------------

# Step 6: Screen for all-zero or near-zero DWI volumes

# ------------------------------------------------------------
# Flag any diffusion volume with a masked mean intensity below 1.
# This is a gross-corruption screening threshold. Any flagged
# volume requires further quantitative and visual investigation.
# Expected result:
# No zero or near-zero volumes detected.
# ------------------------------------------------------------

awk -F, '
NR > 1 && ($4 == 0 || $4 < 1) {
    print "Suspicious volume:", $1, "| masked mean:", $4
    flagged++
}
END {
    if (flagged == 0) {
        print "No zero or near-zero volumes detected."
    }
}' \
/YTH001_BL_upsampled_mask_QC/YTH001_BL_masked_volume_statistics.csv


# ------------------------------------------------------------
# Step 7: Generate the mean b=0 image
# ------------------------------------------------------------
# Extract all b=0 volumes and calculate their mean for visual
# assessment of anatomy, signal coverage, and mask boundaries.

# ------------------------------------------------------------
dwiextract \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-bzero \
- |
mrmath \
- \
mean \
/Users/tanni/Desktop/YTH001_BL_upsampled_mask_QC/YTH001_BL_mean_b0.mif \
-axis 3 \
-force

# ------------------------------------------------------------
# Step 8: Generate the mean b=1500 image
# ------------------------------------------------------------
# Extract all b=1500 s/mm² volumes and calculate their mean for
# visual assessment of diffusion-image integrity.
# ------------------------------------------------------------

dwiextract \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-shell 1500 \
- |
mrmath \
- \
mean \
/Users/tanni/Desktop/YTH001_BL_upsampled_mask_QC/YTH001_BL_mean_b1500.mif \
-axis 3 \
-force


# ------------------------------------------------------------
# Step 9: Generate the mean b=3000 image
# ------------------------------------------------------------
# Extract all b=3000 s/mm² volumes and calculate their mean for
# visual assessment of high-b-value diffusion-image integrity.
# ------------------------------------------------------------

dwiextract \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-shell 3000 \
- |
mrmath \
- \
mean \
/Users/tanni/Desktop/YTH001_BL_upsampled_mask_QC/YTH001_BL_mean_b3000.mif \
-axis 3 \
-force


# ------------------------------------------------------------
# Step 9: Visually inspect the mask over the mean b=0 image
# ------------------------------------------------------------
# Display the binary mask over the mean b=0 image with partial transparency so that the underlying anatomy remains visible.
# Inspect all axial, coronal, and sagittal slices.
# Confirm inclusion of:
# Cerebral cortex
# Intended white matter
# Cerebellum
# Brainstem
# Superior and inferior brain regions

# Check for:
# Internal holes
# Excluded brain tissue
# Cropping
# Boundary errors
# Substantial non-brain inclusion

# Any mask requiring correction should be documented and handled
# in the separate mask-correction step.

# ------------------------------------------------------------

mrview \
/Users/tanni/Desktop/YTH001_BL_upsampled_mask_QC/YTH001_BL_mean_b0.mif \
-overlay.load \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-overlay.opacity 0.3


# ------------------------------------------------------------

# Step 11: Visually inspect the mean b=1500 image

# ------------------------------------------------------------

# Inspect for missing or corrupted slices, unexpected signal loss,
# cropping, severe blurring, ringing, or other abnormalities.
# ------------------------------------------------------------

mrview \
/Users/tanni/Desktop/YTH001_BL_upsampled_mask_QC/YTH001_BL_mean_b1500.mif


# ------------------------------------------------------------
# Step 12: Visually inspect the mean b=3000 image
# ------------------------------------------------------------
# Inspect for missing or corrupted slices, unexpected signal loss,
# cropping, severe blurring, ringing, or other abnormalities.

# ------------------------------------------------------------

mrview \
/Users/tanni/Desktop/YTH001_BL_upsampled_mask_QC/YTH001_BL_mean_b3000.mif


# ------------------------------------------------------------
# QC interpretation
# ------------------------------------------------------------
# The automated QC checks pass when:
# The mask has the expected geometry.
# The mask is binary and non-empty.
# No NaN or infinite voxels are detected.
# No all-zero or near-zero DWI volumes are detected.
# The masked volume statistics show plausible signal intensities.
# Final acceptance still requires visual inspection of all three
# anatomical planes and the diffusion-shell mean images.
# Hole correction is not included in this step.

# ------------------------------------------------------------
