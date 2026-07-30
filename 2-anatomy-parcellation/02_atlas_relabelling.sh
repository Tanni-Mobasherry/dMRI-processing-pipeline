
# ============================================================
# Test subject: YTH001 / BL
#
# Inputs:
#   aparc.a2009s+aseg.mgz
#   Copied from Sjoerd's longitudinal FreeSurfer results:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/FS_longi/FS_BL/mri/"
#   Working copy:
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation/"
#
#   FreeSurferColorLUT.txt
#   "/Applications/freesurfer/FreeSurferColorLUT.txt"
#
#   fs_a2009s.txt
#   "/usr/local/mrtrix3/share/mrtrix3/labelconvert/fs_a2009s.txt"
#
# Outputs:
#   parcels_a2009s.mif
#   "/Volumes/Toshiba-Ext/raw-data/YTH001/BL/dmri/anatomy-parcellation/"
# ============================================================


# ------------------------------------------------------------
# Step 1: Convert the FreeSurfer anatomical parcellation
#         to the MRtrix3 lookup table - Atlas relabelling
# ------------------------------------------------------------

labelconvert \
aparc.a2009s+aseg.mgz \
/Applications/freesurfer/FreeSurferColorLUT.txt \
/usr/local/mrtrix3/share/mrtrix3/labelconvert/fs_a2009s.txt \
parcels_a2009s.mif \
-force
