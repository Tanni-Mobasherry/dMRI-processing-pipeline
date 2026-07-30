# Coregister 5TT to Diffusion Space
# Purpose:
#Coregister the 5TT image to diffusion space using the transformation estimated between
# the diffusion mean b0 image and the longitudinal FreeSurfer brain image.
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
#   brain.mgz
#   Longitudinal FreeSurfer brain image:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/"
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
#   T1_brain.nii.gz
#   d2s.mat
#   d2s.txt
#
# Final output:
#
#   5tt_coreg.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
#==============================================================================================
# Convert the longitudinal FreeSurfer brain image
mrconvert \
    /Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/brain.mgz \
    /Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/T1_brain.nii.gz
# Register diffusion mean b0 to the longitudinal T1 image
flirt \
    -in mean_b0_BA.nii.gz \
    -ref T1_brain.nii.gz \
    -dof 6 \
    -cost normmi \
    -omat d2s.mat
# Convert the FSL transform to MRtrix format
transformconvert \
    d2s.mat \
    mean_b0_BA.nii.gz \
    T1_brain.nii.gz \
    flirt_import \
    d2s.txt


# Transform the 5TT image into diffusion space
mrtransform \
    5tt.mif \
    -linear d2s.txt \
    -inverse \
    5tt_coreg.mif
