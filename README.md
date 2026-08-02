# dMRI-processing-pipeline

A complete diffusion MRI preprocessing, modelling and quality control (QC) pipeline developed for structural connectome analysis.

## Overview

This repository contains the scripts and documentation used for processing diffusion MRI (dMRI) data from raw images to connectome generation.

The pipeline includes:

- Diffusion MRI preprocessing
- Bias field correction
- Denoising
- Gibbs ringing removal
- Motion, eddy-current and susceptibility distortion correction
- Brain masking
- Quality control (QC)
- Anatomical parcellation
- Diffusion-to-T1 registration
- Tissue segmentation
- Fibre orientation distribution (FOD) estimation
- Whole-brain tractography
- SIFT2 filtering
- Connectome construction
- DTI metric calculation (FA, MD, AD, RD)

---

## Repository Structure

```
01_preprocessing/
    Diffusion preprocessing

02_anatomy-parcellation/
    Anatomical processing and registration

03_modelling-connectome/
    Fibre modelling, tractography and connectome generation
```

---

## Software

This pipeline uses:

- MRtrix3
- FSL
- ANTs
- FreeSurfer
- NiftyReg

---

## Author

Tanni Mobasherry

The University of Western Australia
