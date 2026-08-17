
```bash

# 5TT-Based Restoration of an Upsampled DWI Mask
# subject: YTH001 / BL
#
# Restoration rule:
# Restore a candidate voxel when: WM > 0.25 OR CSF < 0.5
#
# Inputs:
#
# mask.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/"
#
# dwi_eddy_BA_upsampled.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/updated-method/upsampling/"
#
# dwi_mask_upsampled.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/updated-method/mask/"
#
# 5tt_coreg.mif
# "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
#

# INTERMEDIATE OUTPUTS
# old_mask_on_upsampled_grid.mif
# all_mask_loss_without_erosion.mif
# old_mask_on_upsampled_grid_eroded.mif
# 5tt_on_upsampled_grid.mif
# 5tt_WM.mif
# 5tt_CSF.mif
# 5tt_total_fraction.mif
# valid_5tt_coverage.mif
# 5tt_restoration_rule_raw.mif
# 5tt_restoration_rule_valid.mif
# candidate_reconstructed_from_two_groups.mif
# candidate_partition_mismatch.mif

# FINAL QC OUTPUTS
# candidate_excluded_voxels.mif
# voxels_to_restore_5tt_supported.mif
# candidate_not_selected_for_restoration.mif
# candidate_without_valid_5tt.mif
# candidate_valid_but_criterion_not_met.mif

# FINAL CORRECTED MASK
# dwi_mask_upsampled_5tt_corrected.mif

# /Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/5tt-based-restoration/

# ============================================================
# ------------------------------------------------------------
# Step 1: Resample the original 2 mm binary mask onto the exact
# upsampled 1.25 mm DWI grid
# Nearest-neighbour interpolation retains binary mask values.
# ------------------------------------------------------------

mrgrid \
mask.mif \
regrid \
-template \
dwi_eddy_BA_upsampled.mif \
-interp nearest \
old_mask_on_upsampled_grid.mif \
-force

# ------------------------------------------------------------
# Step 2: Calculate all old-mask voxels absent from the new mask
# before erosion
#
# all loss = max(resampled old mask - new mask, 0)
#
# This output includes superficial grid-boundary differences and
# is retained for QC only.
# ------------------------------------------------------------

mrcalc \
old_mask_on_upsampled_grid.mif \
dwi_mask_upsampled.mif \
-sub \
0 \
-max \
all_mask_loss_without_erosion.mif \
-datatype bit \
-force

echo "All differences before erosion:"
mrstats \
all_mask_loss_without_erosion.mif \
-ignorezero \
-output count \
-quiet


# ------------------------------------------------------------
# Step 3: Erode the resampled original mask by one voxel
#
# This reduces superficial differences caused by changing from
# the native 2 mm grid to the upsampled 1.25 mm grid.
# ------------------------------------------------------------

maskfilter \
old_mask_on_upsampled_grid.mif \
erode \
old_mask_on_upsampled_grid_eroded.mif \
-npass 1 \
-force


# ------------------------------------------------------------
# Step 4: Identify candidate excluded voxels after erosion
#
# candidate = max(eroded resampled old mask - new mask, 0)
#
# These are candidates for tissue-based assessment; they are not
# automatically interpreted as genuine brain-tissue exclusion.
# ------------------------------------------------------------

mrcalc \
old_mask_on_upsampled_grid_eroded.mif \
dwi_mask_upsampled.mif \
-sub \
0 \
-max \
candidate_excluded_voxels.mif \
-datatype bit \
-force

echo "Candidate excluded voxels after one-voxel erosion:"
mrstats \
candidate_excluded_voxels.mif \
-ignorezero \
-output count \
-quiet


# ------------------------------------------------------------
# Step 5: Resample the co-registered 5TT image onto the exact
# upsampled DWI grid
#
# Linear interpolation is used because 5TT contains continuous
# tissue fractions rather than binary labels.
# ------------------------------------------------------------

mrgrid \
5tt_coreg.mif \
regrid \
-template \
dwi_eddy_BA_upsampled.mif \
-interp linear \
5tt_on_upsampled_grid.mif \
-force


# ------------------------------------------------------------
# Step 6: Extract WM and CSF tissue-fraction maps
#
# Standard MRtrix 5TT volume order:
# volume 2 = WM
# volume 3 = CSF
# ------------------------------------------------------------

mrconvert \
5tt_on_upsampled_grid.mif \
-coord 3 2 \
5tt_WM.mif \
-force

mrconvert \
5tt_on_upsampled_grid.mif \
-coord 3 3 \
5tt_CSF.mif \
-force


# ------------------------------------------------------------
# Step 7: Confirm matching spatial grids
#
# The first three dimensions, spatial voxel sizes, and transform
# matrices must match. The DWI additionally has a fourth volume
# dimension, which is expected and does not prevent voxel-wise QC.
# ------------------------------------------------------------

echo "=== Upsampled DWI ==="
mrinfo \
dwi_eddy_BA_upsampled.mif \
-size \
-spacing \
-transform

echo "=== Candidate excluded voxels ==="
mrinfo \
candidate_excluded_voxels.mif \
-size \
-spacing \
-transform

echo "=== 5TT WM ==="
mrinfo \
5tt_WM.mif \
-size \
-spacing \
-transform

echo "=== 5TT CSF ==="
mrinfo \
5tt_CSF.mif \
-size \
-spacing \
-transform


# ------------------------------------------------------------
# Step 8A: Identify valid 5TT coverage
#
# Sum all five tissue fractions at each voxel. A total fraction
# greater than 1e-6 indicates that 5TT information is available.
# ------------------------------------------------------------

mrmath \
5tt_on_upsampled_grid.mif \
sum \
5tt_total_fraction.mif \
-axis 3 \
-force

mrcalc \
5tt_total_fraction.mif \
0.000001 \
-gt \
valid_5tt_coverage.mif \
-datatype bit \
-force


# ------------------------------------------------------------
# Step 8B: Create the preliminary tissue-based rule
#
# WM > 0.25 OR CSF < 0.5
#
# This raw rule must not be used directly. Outside valid 5TT
# coverage, CSF may equal zero and incorrectly satisfy CSF < 0.5.
# ------------------------------------------------------------

mrcalc \
5tt_WM.mif \
0.25 \
-gt \
5tt_CSF.mif \
0.5 \
-lt \
-or \
5tt_restoration_rule_raw.mif \
-datatype bit \
-force


# ------------------------------------------------------------
# Step 8C: Restrict the rule to valid 5TT coverage
#
# valid rule = raw rule AND valid 5TT coverage
# ------------------------------------------------------------

mrcalc \
5tt_restoration_rule_raw.mif \
valid_5tt_coverage.mif \
-mult \
5tt_restoration_rule_valid.mif \
-datatype bit \
-force


# ------------------------------------------------------------
# Step 9A: Select candidate voxels supported for restoration
#
# selected = candidate excluded voxels AND valid 5TT rule
# ------------------------------------------------------------

mrcalc \
candidate_excluded_voxels.mif \
5tt_restoration_rule_valid.mif \
-mult \
voxels_to_restore_5tt_supported.mif \
-datatype bit \
-force


# ------------------------------------------------------------
# Step 9B: Identify all candidate voxels not selected by the
# valid 5TT rule
# ------------------------------------------------------------

mrcalc \
candidate_excluded_voxels.mif \
voxels_to_restore_5tt_supported.mif \
-sub \
0 \
-max \
candidate_not_selected_for_restoration.mif \
-datatype bit \
-force


# ------------------------------------------------------------
# Step 9C: Separate the reasons for non-selection
#
# Output 1: candidates without valid 5TT coverage
# ------------------------------------------------------------

mrcalc \
1 \
valid_5tt_coverage.mif \
-sub \
candidate_excluded_voxels.mif \
-mult \
candidate_without_valid_5tt.mif \
-datatype bit \
-force


# Output 2: candidates with valid 5TT coverage that did not meet
# the WM/CSF criterion

mrcalc \
1 \
5tt_restoration_rule_valid.mif \
-sub \
valid_5tt_coverage.mif \
-mult \
candidate_excluded_voxels.mif \
-mult \
candidate_valid_but_criterion_not_met.mif \
-datatype bit \
-force


# ------------------------------------------------------------
# Step 10: Report voxel counts
#
# For YTH001/BL, the expected results from the current test are:
# candidates = 616
# selected by valid 5TT rule = 305
# not selected = 311
# without valid 5TT = 311
# valid 5TT but criterion not met = 0
# ------------------------------------------------------------

echo "All candidate excluded voxels:"
mrstats \
candidate_excluded_voxels.mif \
-ignorezero \
-output count \
-quiet

echo "Selected by the valid 5TT restoration rule:"
mrstats \
voxels_to_restore_5tt_supported.mif \
-ignorezero \
-output count \
-quiet

echo "All candidates not selected:"
mrstats \
candidate_not_selected_for_restoration.mif \
-ignorezero \
-output count \
-quiet

echo "Candidates without valid 5TT coverage:"
mrstats \
candidate_without_valid_5tt.mif \
-ignorezero \
-output count \
-quiet

echo "Candidates with valid 5TT but criterion not met:"
mrstats \
candidate_valid_but_criterion_not_met.mif \
-ignorezero \
-output count \
-quiet


# ------------------------------------------------------------
# Step 11: Verify that the two complementary outputs reproduce
# the complete candidate map exactly
#
# The mismatch count and maximum must both equal zero.
# ------------------------------------------------------------

mrcalc \
voxels_to_restore_5tt_supported.mif \
candidate_not_selected_for_restoration.mif \
-or \
candidate_reconstructed_from_two_groups.mif \
-datatype bit \
-force

mrcalc \
candidate_reconstructed_from_two_groups.mif \
candidate_excluded_voxels.mif \
-sub \
-abs \
candidate_partition_mismatch.mif \
-datatype bit \
-force

echo "Candidate-partition mismatched voxels:"
mrstats \
candidate_partition_mismatch.mif \
-ignorezero \
-output count \
-quiet

echo "Maximum candidate-partition mismatch:"
mrstats \
candidate_partition_mismatch.mif \
-output max \
-quiet


# ------------------------------------------------------------
# Step 12: Visual QC
#
# First overlay: candidates selected by the valid 5TT rule
# Second overlay: remaining candidates not selected
#
# In mrview, assign different colour maps (for example, red and
# cyan), turn interpolation off for both binary overlays, inspect
# all three planes, and scroll through adjacent slices.
# ------------------------------------------------------------

mrview \
dwi_eddy_BA_upsampled.mif \
-overlay.load \
voxels_to_restore_5tt_supported.mif \
-overlay.opacity 0.8 \
-overlay.threshold_min 0.5 \
-overlay.threshold_max 1 \
-overlay.interpolation 0 \
-overlay.load \
candidate_not_selected_for_restoration.mif \
-overlay.opacity 0.8 \
-overlay.threshold_min 0.5 \
-overlay.threshold_max 1 \
-overlay.interpolation 0

# ------------------------------------------------------------
#final output-corrected mask
#Step 13:Restore only candidate voxels supported by valid
#5TT coverage and satisfying WM > 0.25 OR CSF < 0.5
# ------------------------------------------------------------

mrcalc \
dwi_mask_upsampled.mif \
voxels_to_restore_5tt_supported.mif \
-or \
dwi_mask_upsampled_5tt_corrected.mif \
-datatype bit \
-force
# ============================================================

# ------------------------------------------------------------
# Step 14: Verify the final mask correction
#
# Report the number of voxels in the original mask, the number
# selected for restoration, and the final corrected mask.
#
# The increase from the original to the corrected mask should
# equal the number of voxels selected for restoration.
# ------------------------------------------------------------

echo "Original mask voxels:"
mrstats \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-ignorezero \
-output count \
-quiet

echo "Voxels selected for restoration:"
mrstats \
voxels_to_restore_5tt_supported.mif \
-ignorezero \
-output count \
-quiet

echo "Corrected mask voxels:"
mrstats \
dwi_mask_upsampled_5tt_corrected.mif \
-ignorezero \
-output count \
-quiet
# ------------------------------------------------------------ # ------------------------------------------------------------ # ------------------------------------------

```
Output interpretation:


| Output file                                   | Interpretation                                                                                                                       
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `old_mask_on_upsampled_grid.mif`              | Original 2 mm mask resampled onto the exact 1.25 mm upsampled DWI grid.                                                                      |
| `all_mask_loss_without_erosion.mif`           | All old-versus-new mask differences before erosion, including superficial grid-boundary differences.                                         |
| `old_mask_on_upsampled_grid_eroded.mif`       | Resampled original mask after one-voxel erosion to reduce superficial boundary effects.                                                      |
| `candidate_excluded_voxels.mif`               | Voxels inside the eroded resampled original mask but absent from the new upsampled mask.                                                     |
| `5tt_on_upsampled_grid.mif`                   | Co-registered 5TT image resampled onto the exact upsampled DWI grid.                                                                         |
| `5tt_WM.mif`                                  | White-matter tissue-fraction map extracted from the resampled 5TT image.                                                                     |
| `5tt_CSF.mif`                                 | CSF tissue-fraction map extracted from the resampled 5TT image.                                                                              |
| `5tt_total_fraction.mif`                      | Sum of all five 5TT tissue fractions at each voxel, used to identify valid 5TT coverage.                                                     |
| `valid_5tt_coverage.mif`                      | Binary map of voxels containing non-zero 5TT tissue information.                                                                             |
| `5tt_restoration_rule_raw.mif`                | Raw `WM > 0.25 OR CSF < 0.5` rule. It is not safe for direct restoration because areas without 5TT information may also satisfy `CSF < 0.5`. |
| `5tt_restoration_rule_valid.mif`              | Raw restoration rule restricted to voxels with valid 5TT coverage.                                                                           |
| `voxels_to_restore_5tt_supported.mif`         | Candidate excluded voxels selected for restoration by the valid 5TT rule.                                                                    |
| `candidate_without_valid_5tt.mif`             | Candidate voxels for which 5TT cannot provide a valid tissue classification.                                                                 |
| `candidate_valid_but_criterion_not_met.mif`   | Candidate voxels with valid 5TT coverage that fail both restoration conditions.                                                              |
| `candidate_not_selected_for_restoration.mif`  | All candidate voxels not selected by the valid 5TT restoration rule.                                                                         |
| `candidate_reconstructed_from_two_groups.mif` | Reconstruction of the complete candidate map by combining selected and non-selected candidate voxels.                                        |
| `candidate_partition_mismatch.mif`            | Verification map comparing the reconstructed and original candidate maps; voxel count and maximum should both be zero.                       |
| `dwi_mask_upsampled_5tt_corrected.mif`        | Final corrected binary DWI mask after adding only candidate voxels supported by the valid 5TT restoration rule.                              |

# ------------------------------------------------------------ 
Important notes
Visually confirm the anatomical alignment of 5tt_coreg.mif with the DWI before relying on its tissue fractions.
The tissue rule is applied only to candidate_excluded_voxels.mif, not to the entire image.
Voxels outside valid 5TT coverage are never automatically selected through the CSF < 0.5 condition.
No automatic fillh, fillh26, dilation, or erosion-based closing is applied to the final mask.
The original dwi_mask_upsampled.mif is never overwritten.
Only dwi_mask_upsampled_5tt_corrected.mif should be treated as the final corrected mask; the other files are intermediate or QC outputs.

#find restore components
RESTORE=~/Desktop/YTH001_BL_voxels_to_restore.mif
COMPONENTS=~/Desktop/YTH001_BL_restored_components.mif

maskfilter \
"$RESTORE" \
connect \
"$COMPONENTS" \
-force
----
max_label=$(mrstats \
"$COMPONENTS" \
-output max \
-quiet | awk '{printf "%d",$1}')

for i in $(seq 1 "$max_label"); do
    count=$(mrcalc \
        "$COMPONENTS" \
        "$i" \
        -eq \
        -quiet - | \
        mrstats - \
        -ignorezero \
        -output count \
        -quiet)

    printf "%s voxels | component %s\n" "$count" "$i"
done | sort -nr
------
