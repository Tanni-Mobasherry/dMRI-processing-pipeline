#!/bin/bash

# DEC-FA generation
# Example dataset: YTH001/BL
# Run this script from inside the YTH001 directory.

DATA="BL"

dwi2mask \
  "$DATA/dmri/preproc/bias-field-correction/dwi_eddy_BA.mif" \
  "$DATA/dmri/dec-fa/mask.mif"

dwi2tensor \
  "$DATA/dmri/preproc/bias-field-correction/dwi_eddy_BA.mif" \
  "$DATA/dmri/dec-fa/tensor.mif" \
  -mask "$DATA/dmri/dec-fa/mask.mif"

tensor2metric \
  "$DATA/dmri/dec-fa/tensor.mif" \
  -fa "$DATA/dmri/dec-fa/fa.mif" \
  -vector "$DATA/dmri/dec-fa/dec.mif"

mrcalc \
  "$DATA/dmri/dec-fa/dec.mif" -abs \
  "$DATA/dmri/dec-fa/fa.mif" -mult \
  "$DATA/dmri/dec-fa/dec_fa.mif"
