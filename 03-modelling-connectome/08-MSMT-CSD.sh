
# ================================================================
# Multi-Shell Multi-Tissue Constrained Spherical
#            Deconvolution (MSMT-CSD)
#
# Purpose:
# Estimate the white-matter fibre orientation distribution (WM FOD)
# together with the grey-matter and CSF tissue compartments using
# multi-shell multi-tissue constrained spherical deconvolution
# (MSMT-CSD).
#
# Group-average WM, GM and CSF response functions are used for all
# subjects and time points to improve consistency across the study.
#
# Subject:
#   YTH001_BL
#
# Inputs:
#
# dwi_eddy_BA_upsampled.mif
#
#    Bias-field-corrected DWI upsampled from 2 mm to 1.25 mm
#    isotropic resolution.
#
#    Generated during longitudinal DWI upsampling:
#
#    dmri/modelling-connectome/longitudinal/upsampling/
#    dwi_eddy_BA_upsampled.mif
#
#
#   dmri/modelling-connectome/longitudinal/group-response-function/"
#
# dwi_mask_upsampled_5tt_corrected.mif
#
#    Final 1.25 mm DWI mask after comparison with the original
#    mask and restoration of excluded voxels supported by the
#    co-registered 5TT image.
#
#    Restoration criteria:
#
#      WM fraction > 0.25
#      OR
#      CSF fraction < 0.5
#
#    Stored in:
#
#    dmri/modelling-connectome/longitudinal/mask/dwi_mask_upsampled_5tt_corrected.mif
#
# Outputs:
#
#   wmfod.mif
#   White-matter fibre orientation distribution (FOD)
#
#   gm.mif
#   Grey-matter tissue compartment
#
#   csf.mif
#   Cerebrospinal fluid tissue compartment
#
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/group-MSMT-CSD/"
# ================================================================


# ------------------------------------------------------------
# Step 1: Estimate WM FOD, GM and CSF tissue compartments
#         using MSMT-CSD with the group-average response
#         functions
# ------------------------------------------------------------

dwi2fod msmt_csd \
    dwi_eddy_BA_upsampled.mif \
    group_wm.txt wmfod.mif \
    group_gm.txt gm.mif \
    group_csf.txt csf.mif \
    -mask dwi_mask_upsampled_5tt_corrected.mif \
    -force
