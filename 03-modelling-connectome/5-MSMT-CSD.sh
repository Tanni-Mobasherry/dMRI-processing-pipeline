
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
#   dwi_eddy_BA.mif
#   Bias-field-corrected diffusion image generated in:
#   1-diffusion-preprocessing/09_bias_field_correction.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/"
#
#   group_wm.txt
#   group_gm.txt
#   group_csf.txt
#   Group-average tissue response functions generated in:
#   3-modelling-connectome/04_group_response_function.sh
#   "/Volumes/Toshiba-Ext/raw-data/group-response-function/"
#
#   mask.mif
#   Diffusion brain mask generated from the bias-field-corrected
#   DWI in:
#   1-diffusion-preprocessing/10_dec_fa.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/"
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
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/MSMT-CSD/"
# ================================================================


# ------------------------------------------------------------
# Step 1: Estimate WM FOD, GM and CSF tissue compartments
#         using MSMT-CSD with the group-average response
#         functions
# ------------------------------------------------------------

dwi2fod msmt_csd \
    dwi_eddy_BA.mif \
    group_wm.txt wmfod.mif \
    group_gm.txt gm.mif \
    group_csf.txt csf.mif \
    -mask mask.mif \
    -force
