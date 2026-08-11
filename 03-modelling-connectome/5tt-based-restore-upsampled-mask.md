# 5TT-Based Restoration of an Upsampled DWI Mask

```bash
#!/bin/bash

# ============================================================
# Supplementary QC: 5TT-Based Restoration of an Upsampled
# DWI Mask
# Test subject: YTH001 / BL
#
# Restoration rule:
# Restore a candidate voxel when:
#
# WM > 0.25 OR CSF < 0.5
#
# Inputs:
#
# mask.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/"
#
# dwi_eddy_BA_upsampled.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/"
#
# dwi_mask_upsampled.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/"
#
# 5tt_coreg.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
#
# Intermediate outputs:
#
# old_mask_on_upsampled_grid.mif
# old_mask_on_upsampled_grid_eroded.mif
# 5tt_on_upsampled_grid.mif
# 5tt_cortical_GM.mif
# 5tt_subcortical_GM.mif
# 5tt_all_GM.mif
# 5tt_WM.mif
# 5tt_CSF.mif
# 5tt_restoration_rule.mif
# candidate_GM_ge_0.5.mif
# candidate_GM_ge_0.75.mif
# voxels_actually_added.mif
# restoration_mismatch.mif
#
# Main outputs:
#
# candidate_excluded_voxels.mif
# voxels_supported_by_5tt.mif
# voxels_not_supported_by_5tt.mif
# dwi_mask_upsampled_5tt_corrected.mif
#
# Output directory:
# "/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/"
# ============================================================

# ------------------------------------------------------------
# Step 1: Create the output directory
# ------------------------------------------------------------

mkdir -p \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration

# ------------------------------------------------------------
# Step 2: Resample the original 2 mm binary mask onto the
# upsampled 1.25 mm DWI grid
#
# Nearest-neighbour interpolation is used to retain binary values.
# ------------------------------------------------------------

mrgrid \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/mask.mif \
regrid \
-template \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-interp nearest \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/old_mask_on_upsampled_grid.mif \
-force

# ------------------------------------------------------------
# Step 3: Erode the resampled original mask by one voxel
#
# This reduces superficial differences caused by changing from
# the 2 mm grid to the 1.25 mm grid.
# ------------------------------------------------------------

maskfilter \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/old_mask_on_upsampled_grid.mif \
erode \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/old_mask_on_upsampled_grid_eroded.mif \
-npass 1 \
-force

# ------------------------------------------------------------
# Step 4: Identify candidate excluded voxels
#
# candidate = max(eroded_original_mask - new_mask, 0)
#
# These voxels require tissue-based assessment and are not yet
# confirmed as genuine brain-tissue exclusion.
# ------------------------------------------------------------

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/old_mask_on_upsampled_grid_eroded.mif \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-sub \
0 \
-max \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-datatype bit \
-force

# Count all candidate excluded voxels.

mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-ignorezero \
-output count

# ------------------------------------------------------------
# Step 5: Resample the co-registered 5TT image onto the exact
# upsampled DWI grid
#
# Linear interpolation is used because the 5TT volumes contain
# continuous tissue fractions.
# ------------------------------------------------------------

mrgrid \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/5tt_coreg.mif \
regrid \
-template \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-interp linear \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_on_upsampled_grid.mif \
-force

# ------------------------------------------------------------
# Step 6: Extract the 5TT tissue compartments
#
# Standard MRtrix 5TT volume order:
# 0 = cortical GM
# 1 = subcortical GM
# 2 = WM
# 3 = CSF
# 4 = pathological tissue
# ------------------------------------------------------------

mrconvert \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_on_upsampled_grid.mif \
-coord 3 0 \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_cortical_GM.mif \
-force

mrconvert \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_on_upsampled_grid.mif \
-coord 3 1 \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_subcortical_GM.mif \
-force

mrconvert \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_on_upsampled_grid.mif \
-coord 3 2 \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_WM.mif \
-force

mrconvert \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_on_upsampled_grid.mif \
-coord 3 3 \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_CSF.mif \
-force

# Combine cortical and subcortical GM for supplementary QC.

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_cortical_GM.mif \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_subcortical_GM.mif \
-add \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_all_GM.mif \
-force

# ------------------------------------------------------------
# Step 7: Confirm that all tissue maps have the same grid as
# the upsampled DWI and candidate-exclusion map
# ------------------------------------------------------------

mrinfo \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-size \
-spacing

mrinfo \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-size \
-spacing

mrinfo \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_WM.mif \
-size \
-spacing

mrinfo \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_CSF.mif \
-size \
-spacing

# ------------------------------------------------------------
# Step 8: Create the 5TT restoration rule
#
# Restore when:
# WM > 0.25 OR CSF < 0.5
# ------------------------------------------------------------

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_WM.mif \
0.25 \
-gt \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_CSF.mif \
0.5 \
-lt \
-or \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_restoration_rule.mif \
-datatype bit \
-force

# ------------------------------------------------------------
# Step 9: Divide all candidate voxels into two complementary
# outputs
#
# Supported output:
# candidate AND restoration_rule
#
# Not-supported output:
# candidate MINUS supported
# ------------------------------------------------------------

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_restoration_rule.mif \
-mult \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_supported_by_5tt.mif \
-datatype bit \
-force

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_supported_by_5tt.mif \
-sub \
0 \
-max \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_not_supported_by_5tt.mif \
-datatype bit \
-force

# ------------------------------------------------------------
# Step 10: Create supplementary GM-overlap maps
#
# These maps are for interpretation only and are not separate
# restoration rules.
# ------------------------------------------------------------

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_all_GM.mif \
0.5 \
-ge \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-mult \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_GM_ge_0.5.mif \
-datatype bit \
-force

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_all_GM.mif \
0.75 \
-ge \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-mult \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_GM_ge_0.75.mif \
-datatype bit \
-force

# ------------------------------------------------------------
# Step 11: Add only the 5TT-supported voxels to the new mask
#
# corrected_mask = new_mask OR supported_voxels
#
# The original upsampled mask is not overwritten.
# ------------------------------------------------------------

mrcalc \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_supported_by_5tt.mif \
-or \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/dwi_mask_upsampled_5tt_corrected.mif \
-datatype bit \
-force

# ------------------------------------------------------------
# Step 12: Count the original, candidate, supported,
# not-supported, and corrected-mask voxels
# ------------------------------------------------------------

echo "Original upsampled mask:"
mrstats \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-ignorezero \
-output count \
-quiet

echo "Candidate excluded voxels:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-ignorezero \
-output count \
-quiet

echo "Voxels supported by 5TT:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_supported_by_5tt.mif \
-ignorezero \
-output count \
-quiet

echo "Voxels not supported by 5TT:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_not_supported_by_5tt.mif \
-ignorezero \
-output count \
-quiet

echo "Corrected mask:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/dwi_mask_upsampled_5tt_corrected.mif \
-ignorezero \
-output count \
-quiet

echo "Candidate voxels with GM >= 0.5:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_GM_ge_0.5.mif \
-ignorezero \
-output count \
-quiet

echo "Candidate voxels with GM >= 0.75:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_GM_ge_0.75.mif \
-ignorezero \
-output count \
-quiet

# ------------------------------------------------------------
# Step 13: Calculate the mean and maximum tissue fractions
# inside all candidate excluded voxels
# ------------------------------------------------------------

echo "GM fraction within candidate voxels: mean and maximum"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_all_GM.mif \
-mask \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-output mean \
-output max \
-quiet

echo "WM fraction within candidate voxels: mean and maximum"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_WM.mif \
-mask \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-output mean \
-output max \
-quiet

echo "CSF fraction within candidate voxels: mean and maximum"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/5tt_CSF.mif \
-mask \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/candidate_excluded_voxels.mif \
-output mean \
-output max \
-quiet

# ------------------------------------------------------------
# Step 14: Confirm which voxels were actually added to the mask
# ------------------------------------------------------------

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/dwi_mask_upsampled_5tt_corrected.mif \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-sub \
0 \
-max \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_actually_added.mif \
-datatype bit \
-force

echo "Voxels actually added to the corrected mask:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_actually_added.mif \
-ignorezero \
-output count \
-quiet

# ------------------------------------------------------------
# Step 15: Confirm that the actually added voxels exactly match
# the 5TT-supported restoration map
#
# A correct result must contain:
# 0 mismatched voxels
# maximum mismatch = 0
# ------------------------------------------------------------

mrcalc \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_actually_added.mif \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_supported_by_5tt.mif \
-sub \
-abs \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/restoration_mismatch.mif \
-force

echo "Number of mismatched voxels:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/restoration_mismatch.mif \
-ignorezero \
-output count \
-quiet

echo "Maximum mismatch:"
mrstats \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/restoration_mismatch.mif \
-output max \
-quiet

# ------------------------------------------------------------
# Step 16: Perform visual QC of both complementary outputs
#
# In mrview:
# turn interpolation off for the binary overlays;
# inspect axial, coronal, and sagittal views;
# scroll through adjacent slices; and
# pay particular attention to the corpus callosum and nearby CSF.
# ------------------------------------------------------------

mrview \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-overlay.load \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_supported_by_5tt.mif \
-overlay.opacity 0.8 \
-overlay.threshold_min 0.5 \
-overlay.threshold_max 1 \
-overlay.interpolation 0 \
-overlay.load \
/Users/tanni/Desktop/loss-check-upsampled/YTH001/BL/5tt-based-restoration/voxels_not_supported_by_5tt.mif \
-overlay.opacity 0.8 \
-overlay.threshold_min 0.5 \
-overlay.threshold_max 1 \
-overlay.interpolation 0
```

## Output interpretation

| Output | Meaning |
| --- | --- |
| `candidate_excluded_voxels.mif` | Voxels inside the eroded original mask but absent from the new upsampled mask |
| `voxels_supported_by_5tt.mif` | Candidate voxels satisfying `WM > 0.25 OR CSF < 0.5`; these are restored |
| `voxels_not_supported_by_5tt.mif` | Candidate voxels satisfying neither condition; these remain excluded and are saved for QC |
| `dwi_mask_upsampled_5tt_corrected.mif` | Original upsampled mask plus only the 5TT-supported candidate voxels |
| `voxels_actually_added.mif` | Exact voxels added to the corrected mask |
| `restoration_mismatch.mif` | Difference between intended and actual restoration; its count and maximum must both be zero |

## Important notes

- The alignment of `5tt_coreg.mif` with the upsampled DWI must be visually confirmed before using this procedure.
- The restoration rule is applied only inside `candidate_excluded_voxels.mif`, not across the whole image.
- The original `dwi_mask_upsampled.mif` is not overwritten.
- No automatic `fillh`, `fillh26`, or dilation-and-erosion closing is applied to the final mask.
- The corrected mask should remain a separate QC output until it has passed visual inspection.
