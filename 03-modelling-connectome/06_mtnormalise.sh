# ================================================================
# Multi-Tissue Log-Domain Intensity Normalisation
#
# Purpose:
# Intensity-normalise the white-matter fibre orientation
# distribution (WM FOD), grey-matter tissue compartment, and CSF
# tissue compartment using MRtrix3 mtnormalise.
#
# This step improves the consistency of tissue amplitudes across
# subjects and time points for subsequent group-level analysis.
#
# Subject:
#   YTH001_BL
#
# Inputs:
#
#   wmfod.mif
#   gm.mif
#   csf.mif
#   Tissue compartments generated in:
#   3-modelling-connectome/05_MSMT-CSD.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/MSMT-CSD/"
#
#   mask.mif
#   Diffusion brain mask generated from the bias-field-corrected
#   DWI in:
#   1-diffusion-preprocessing/10_dec_fa.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/"
#
# Outputs:
#
#   wmfod_norm.mif
#   Intensity-normalised white-matter FOD
#
#   gm_norm.mif
#   Intensity-normalised grey-matter tissue compartment
#
#   csf_norm.mif
#   Intensity-normalised CSF tissue compartment
#
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/mtnormalise/"
# ================================================================


# ------------------------------------------------------------
# Step 1: Intensity-normalise the WM FOD, GM and CSF
#         tissue compartments
# ------------------------------------------------------------

mtnormalise \
    wmfod.mif wmfod_norm.mif \
    gm.mif gm_norm.mif \
    csf.mif csf_norm.mif \
    -mask mask.mif \
    -force
