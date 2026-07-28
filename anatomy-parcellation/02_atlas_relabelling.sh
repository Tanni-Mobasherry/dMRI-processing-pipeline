#!/bin/bash

# Atlas Relabelling — MRtrix labelconvert
# Example dataset: YTH001/BL
#
# Input:
#   aparc.a2009s+aseg.mgz
#   (from Sjoerd's longitudinal FreeSurfer processing)
#
# Output:
#   parcels_a2009s.mif
#
# Run this script from inside the YTH001 directory.

DATA="BL"
AP="$DATA/dmri/anatomy-parcellation"

labelconvert \
  "$AP/aparc.a2009s+aseg.mgz" \
  /Applications/freesurfer/FreeSurferColorLUT.txt \
  /usr/local/mrtrix3/share/mrtrix3/labelconvert/fs_a2009s.txt \
  "$AP/parcels_a2009s.mif" \
  -force
