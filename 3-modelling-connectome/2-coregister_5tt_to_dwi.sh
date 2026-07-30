# Coregister 5TT to Diffusion Space
# Purpose:
# Coregister the five-tissue-type (5TT) image from longitudinal
# FreeSurfer anatomical space into diffusion space by estimating a
# rigid-body transformation between the brain-masked mean b0 image
#and the longitudinal T1 brain image.
#
# Subject:
#   YTH001_BL
#
# Inputs:
#
#   mean_b0_BA.mif
#   Generated in step:
#    1-diffusion-preprocessing/08_bias_field_correction.sh 
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preproc/bias-field-correction/"
#
#   dwi_mask.mif
#   Generated in step 
#   2-anatomy-parcellation/03_parcellation_to_diffusion_space.sh step:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation/"
#
#   T1_brain.nii.gz generated from 
#   (brain.mgz Longitudinal FreeSurfer brain image) in step:
#    02-anatomy-parcellation/03_parcellation_to_diffusion_niftyreg.sh 
#    "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation/"
#
#   5tt.mif
#   Generated in step 
#   3-modelling-connectome/01_generate_5tt.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
#
# Intermediate outputs:
#
#   mean_b0_BA_brain.mif
#   mean_b0_BA_brain.nii.gz
#   d2s.mat
#   d2s.txt
#
# Final output:
#
#   5tt_coreg.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
#==============================================================================================
# ------------------------------------------------------------
# Step 1: Create a brain-only mean b0 image & convert to NIFTI
# ------------------------------------------------------------
mrcalc \
    mean_b0_BA.mif \
    dwi_mask.mif \
    -mult \
    mean_b0_BA_brain.mif \
    -force \
    -quiet
mrconvert \
    mean_b0_BA_brain.mif \
    mean_b0_BA_brain.nii.gz \
    -force \
    -quiet
# ------------------------------------------------------------
# Step 2:  Estimates diffusion-to-structural registration.
# ------------------------------------------------------------
flirt \
    -in mean_b0_BA_brain.nii.gz \
    -ref T1_brain.nii.gz \
    -dof 6 \
    -cost normmi \
    -omat d2s.mat
# ------------------------------------------------------------
# Step 3: Convert the FSL transform to MRtrix format
# ------------------------------------------------------------
transformconvert \
    d2s.mat \
    mean_b0_BA_brain.nii.gz
    T1_brain.nii.gz \
    flirt_import \
    d2s.txt

# -------------------------------------------------------------------------------------------------------
# Step 4: Applies the transformation in the opposite direction, moving the T1-derived 5TT into diffusion
# -------------------------------------------------------------------------------------------------------
mrtransform \
    5tt.mif \
    -linear d2s.txt \
    -inverse \
    -template mean_b0_BA.mif \
    5tt_coreg.mif
    -force
