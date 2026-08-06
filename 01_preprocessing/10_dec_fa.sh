# ============================================================
# DEC-FA validation
# Test subject: YTH001 / BL
#
# Input:
#   dwi_eddy_BA.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/bias-field-correction/dwi_eddy_BA.mif"
#
# Outputs:
#   mask.mif
#   tensor.mif
#   fa.mif
#   dec.mif
#   dec_fa.mif
#"/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/"
# ============================================================

# ------------------------------------------------------------
# Step 1: Generate a brain mask from the bias-corrected DWI
# ------------------------------------------------------------

dwi2mask \
dwi_eddy_BA.mif \
mask.mif \
-force


# ------------------------------------------------------------
# Step 2: Fit the diffusion tensor model
# ------------------------------------------------------------

dwi2tensor \
dwi_eddy_BA.mif \
tensor.mif \
-mask mask.mif \
-force


# ------------------------------------------------------------
# Step 3: Compute fractional anisotropy (FA) and the
#         principal diffusion direction (DEC)
# ------------------------------------------------------------

tensor2metric \
tensor.mif \
-fa fa.mif \
-vector dec.mif \
-force


# ------------------------------------------------------------
# Step 4: Generate the directionally encoded colour FA image
# ------------------------------------------------------------

mrcalc \
dec.mif \
-abs \
fa.mif \
-mult \
dec_fa.mif \
-force

