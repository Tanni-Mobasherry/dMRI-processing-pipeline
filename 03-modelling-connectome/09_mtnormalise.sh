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
# Subject/session:
#   YTH001/BL
#
# Inputs:
#
#   wmfod.mif
#   gm.mif
#   csf.mif
#
#   Multi-tissue compartments generated using the group-average
#   response functions and stored in:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/updated-method/group-MSMT-CSD/"
#
#   dwi_mask_upsampled_5tt_corrected.mif
#
#   Corrected 1.25 mm diffusion mask. The initial mask was generated
#   from the upsampled DWI and subsequently checked against the
#   co-registered 5TT tissue compartments. Excluded voxels were
#   restored according to the agreed 5TT tissue criteria:
#   WM > 0.25 or CSF < 0.5.
#
#   The corrected mask is stored in:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/updated-method/mask/dwi_mask_upsampled_5tt_corrected.mif"
#
# Output directory:
#
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/updated-method/mtnormalise/"
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
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/updated-method/mtnormalise"
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
