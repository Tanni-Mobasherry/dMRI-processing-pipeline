# ============================================================
# DTI metrics
# Test subject: YTH001 / BL
#
# Input:
#   tensor.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dec-fa/tensor.mif"
#
# Outputs:
#   md.mif
#   ad.mif
#   rd.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/preprocessing/dti-metrics/"
# ============================================================


# ------------------------------------------------------------
# Step 1: Compute mean diffusivity (MD), axial diffusivity (AD),
#         and radial diffusivity (RD) from the fitted tensor
# ------------------------------------------------------------

tensor2metric \
tensor.mif \
-adc md.mif \
-ad ad.mif \
-rd rd.mif \
-force
