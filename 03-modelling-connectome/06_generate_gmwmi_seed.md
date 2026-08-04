

# ================================================================
# YTH001_BL: Generate Grey Matter–White Matter Interface Seed Mask
#
# Purpose:
# Generate the grey matter–white matter interface (GMWMI) seed
# mask from the five-tissue-type image registered to diffusion
# space. This mask is used to seed anatomically constrained
# tractography.
#
# Subject:
#   YTH001_BL
#
# Input:
#
#   5tt_coreg.mif
#   Five-tissue-type image registered to diffusion space generated in:
#   3-modelling-connectome/02_coregister_5tt_to_dwi.sh
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
#
# Output:
#
#   gmwmSeed.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
# ================================================================

# ------------------------------------------------------------
# Step 1: Generate the GM–WM interface seed mask from the
#         diffusion-space 5TT image
# ------------------------------------------------------------

5tt2gmwmi \
    5tt_coreg.mif \
    gmwmSeed.mif \
    -force
