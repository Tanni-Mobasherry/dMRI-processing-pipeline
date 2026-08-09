# 05. DWI Upsampling and Brain Mask

## Purpose

## Input
dwi_eddy_BA.mif

## 1. DWI upsampling
mrgrid -vox 1.25

## Output
dwi_eddy_BA_upsampled.mif

## Upsampling QC
2 mm → 1.25 mm
104×104×72×198 → 166×166×115×198


# 05. DWI Upsampling and Brain Mask Generation

## Purpose

#This step upsamples the bias-field-corrected diffusion-weighted image from 2 mm to 1.25 mm isotropic spatial resolution
#before fibre orientation distribution estimation.
#Upsampling improves anatomical contrast and supports subsequent within-subject template construction, registration, and quantitative tractography.
#A new brain mask is then generated directly from the upsampled DWI.
#This ensures that the DWI and mask have identical dimensions, voxel sizes, and spatial transformations.

#Subject: YTH001/BL

## Input

#Bias-field-corrected DWI:
#YTH001/BL/dmri/preprocessing/bias-field-correction/
#└── dwi_eddy_BA.mif

## output
#YTH001/BL/dmri/modelling-connectome/longitudinal/upsampling/
#└── dwi_eddy_BA_upsampled.mif

