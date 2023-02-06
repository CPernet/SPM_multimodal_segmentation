# SPM multimodal segmentation

This repository containst the [BIDS](https://bids.neuroimaging.io/) compliant code for a 2*2 analysis of the multimodal segmentation: T1w vs. T1w & T2w images (i.e. unimodal vs multimodal) but also 1 Gaussian vs. 2 Gaussians per main tissue classes.

## Environement

[Matlab](https://se.mathworks.com/) & [SPM12 ](https://www.fil.ion.ucl.ac.uk/spm/)  

# Background

  
# Analysis

## Code

`multispectral_segmentation_analysis.m` is the analysis script.  

`segment_images.m` is a function performs the SPM12 segmentation, get tissue volume information, get voxel distributions of normalized images, compute the DARTEL template and returns the SPM batch jobs.  

`create_decile_images.m` is a function to compute a decile image which can then be used to mperorm a shif function analysis.

## Results

