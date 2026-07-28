# ================================================================
# YTH001_BL: Generate 5-Tissue-Type (5TT) Image
#
# Purpose:
# Generate a five-tissue-type (5TT) image from the
# FreeSurfer aparc.a2009s+aseg segmentation.
#
# Subject:
#   YTH001_BL
#
# INPUT="/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation/aparc.a2009s+aseg.mgz"
#
# Output:
#   5tt.mif
#OUTPUT_DIR="/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/modelling-connectome"
# ================================================================

5ttgen freesurfer aparc.a2009s+aseg.mgz 5tt.mif
# Validate the generated 5TT image
5ttcheck 5tt.mif
