# ================================================================
# Step 14: Generate Five-Tissue-Type (5TT) Image
# Test subject: YTH001 / BL
#
# Purpose:
# Generate a five-tissue-type (5TT) image from the longitudinal
# FreeSurfer aparc.a2009s+aseg segmentation for anatomically
# constrained tractography.
#
# Input:
#   aparc.a2009s+aseg.mgz
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation/"
#    Copied from Sjoerd's longitudinal FreeSurfer results:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/"
# Outputs:
#   5tt.mif
#   5ttcheck.log
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome/"
# ================================================================

# ------------------------------------------------------------
# Step 1: Generate the five-tissue-type image
# ------------------------------------------------------------
5ttgen freesurfer \
    aparc.a2009s+aseg.mgz \
    5tt.mif \
    -force

# ------------------------------------------------------------
# Step 2: Validate the generated 5TT image
# ------------------------------------------------------------
5ttcheck \
    5tt.mif \
    > 5ttcheck.log 2>&1
