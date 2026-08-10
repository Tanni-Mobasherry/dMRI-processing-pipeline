Quantitative QC for Tissue Exclusion in Upsampled DWI Masks

Purpose

This supplementary QC procedure evaluates whether a newly generated 1.25 mm DWI mask excludes brain tissue that was included in the original native-resolution 2 mm mask.

The procedure:

resamples the original mask onto the upsampled DWI grid using nearest-neighbour interpolation;

identifies regions present in the original mask but absent from the new mask;

erodes the resampled original mask by one voxel to reduce differences caused only by boundary discretisation;

identifies and measures connected components within the remaining differences; and

visually determines whether the differences correspond to brain tissue or to CSF, ventricles, fissures, or minor boundary effects.

This comparison is a supplementary QC check. A non-zero difference image does not by itself demonstrate tissue loss; anatomical location must be assessed on the upsampled DWI.

Required inputs

For each subject and session:

Original native-resolution mask (2 mm):dmri/preprocessing/dec-fa/mask.mif

Upsampled DWI (1.25 mm):dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif

Newly generated upsampled mask (1.25 mm):dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif

Example: YTH001/BL

1. Resample the original mask onto the upsampled DWI grid

Nearest-neighbour interpolation is used because the input is a binary mask.

mrgrid \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/mask.mif \
regrid \
-template \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-interp nearest \
~/Desktop/YTH001_BL_old_mask_on_upsampled_grid.mif

2. Calculate all regions present in the old mask but absent from the new mask

The operation is:

max(old_mask_on_new_grid - new_mask, 0)

mrcalc \
~/Desktop/YTH001_BL_old_mask_on_upsampled_grid.mif \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-sub \
0 \
-max \
~/Desktop/YTH001_BL_mask_lost_regions.mif

Count the difference voxels:

mrstats \
~/Desktop/YTH001_BL_mask_lost_regions.mif \
-ignorezero \
-output count

Differences at this stage commonly include boundary shifts introduced by changing the voxel grid and should not automatically be interpreted as excluded tissue.

3. Reduce boundary-related differences

Erode the resampled original mask by one voxel:

maskfilter \
~/Desktop/YTH001_BL_old_mask_on_upsampled_grid.mif \
erode \
~/Desktop/YTH001_BL_old_mask_eroded.mif \
-npass 1

Subtract the new mask from the eroded original mask:

mrcalc \
~/Desktop/YTH001_BL_old_mask_eroded.mif \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled.mif \
-sub \
0 \
-max \
~/Desktop/YTH001_BL_deep_mask_loss.mif

Count the remaining voxels:

mrstats \
~/Desktop/YTH001_BL_deep_mask_loss.mif \
-ignorezero \
-output count

The resulting image identifies differences that persist after removing one voxel from the old mask boundary. These are candidates for closer review, not confirmed tissue loss.

4. Label connected components

maskfilter \
~/Desktop/YTH001_BL_deep_mask_loss.mif \
connect \
~/Desktop/YTH001_BL_deep_loss_components.mif

List the ten largest connected components:

max_label=$(mrstats \
~/Desktop/YTH001_BL_deep_loss_components.mif \
-output max \
-quiet | awk '{printf "%d",$1}')

for i in $(seq 1 "$max_label"); do
    count=$(mrcalc \
    ~/Desktop/YTH001_BL_deep_loss_components.mif \
    "$i" \
    -eq \
    -quiet \
    - | mrstats \
    - \
    -ignorezero \
    -output count \
    -quiet)

    printf "%s voxels | component %s\n" "$count" "$i"
done | sort -nr | head

5. Extract the largest component

In this example, component 1 was the largest component.

mrcalc \
~/Desktop/YTH001_BL_deep_loss_components.mif \
1 \
-eq \
~/Desktop/YTH001_BL_deep_loss_component_1.mif

6. Perform anatomical visual QC

mrview \
/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/dwi_eddy_BA_upsampled.mif \
-overlay.load \
~/Desktop/YTH001_BL_deep_loss_component_1.mif \
-overlay.opacity 0.8

In mrview:

turn interpolate off for the binary overlay;

inspect axial, coronal, and sagittal views;

scroll through adjacent slices;

determine whether the component overlaps visible brain tissue; and

distinguish tissue exclusion from ventricular CSF, sulcal CSF, fissures, and small grid-boundary differences.

Interpretation

Observation

Interpretation

QC classification

Difference lies within visible brain tissue

Potential exclusion of tissue intended for analysis

REVIEW

Difference is confined to ventricles, CSF, or fissures

Not meaningful brain-tissue exclusion

PASS

Only thin peripheral differences remain

Likely grid or boundary discretisation effect

Usually PASS, after visual confirmation

Location remains uncertain

Further inspection or supervisor review required

REVIEW

An overly generous mask is generally preferable to a mask that excludes brain tissue intended for analysis, because subsequent FOD estimation is restricted to the mask and the longitudinal analysis mask may be based on the intersection of session masks.

YTH001/BL example result

All old-versus-new differences: 8,624 voxels

Differences remaining after one-voxel erosion: 616 voxels

Largest connected component: 557 voxels

Anatomical assessment: the largest component was predominantly confined to ventricular/CSF spaces rather than visible brain tissue

Final interpretation: no meaningful brain-tissue exclusion identified (PASS)

At 1.25 mm isotropic resolution, each voxel has a volume of 1.953125 mm³. Voxel count alone should not determine the QC result; anatomical location is the decisive criterion.
