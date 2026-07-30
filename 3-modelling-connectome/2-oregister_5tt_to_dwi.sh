# YTH001_BL: Coregister 5TT to Diffusion Space
# Purpose:
#   Estimate a rigid-body transformation between the diffusion
#   mean b0 image and the longitudinal FreeSurfer T1 image,
#   convert the FSL transformation into MRtrix format, and apply
#   the inverse transformation to the 5TT image.
# Subject:
#   YTH001_BL
# Inputs
#   /Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/5tt.mif
#   /Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preproc/bias-field-correction/mean_b0_BA.nii.gz
#   /Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/brain.mgz
#
# Output:
#   d2s.mat
#   d2s.txt
#   5tt_coreg.mif
#
# Output Directory: /Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/
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
