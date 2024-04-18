# SPM multimodal segmentation

This repository containst the [BIDS](https://bids.neuroimaging.io/) compliant code for a 2*2 analysis of the multimodal segmentation: T1w vs. T1w & T2w images (i.e. unimodal vs multimodal) but also 1 Gaussian vs. 2 Gaussians per main tissue classes. Our analysis showed that multispectral segmentation is more accurate but 2 Gaussians per tissue class (GM, WM and CSF) must be used.

## Environement

[Matlab](https://se.mathworks.com/) & [SPM12 ](https://www.fil.ion.ucl.ac.uk/spm/) to run the image analysis.

# Background

Using multimodal segmentation is supposed to be more accurate, but some studies suggested otherwise. One issue is that the generative model used is often not updated, and it seems that the discrepancy between those results and expectations comes from there. We thus segmented data 4 times, using T1w images only or using T1w and T2w images, and using 1 or 2 Gaussians per brain tissue.
  
# Analysis

## Code

[multispectral_segmentation_analysis.m](/code/multispectral_segmentation_analysis.m) is the script used to compute all the metrics and thus obtain the results.  
[segment_images.m](code/segment_images.m) is a subfunction called to perform the SPM12 segmentation, get tissue volume information, get voxel distributions of normalized images, compute the DARTEL template and return the SPM batch jobs.  
[create_decile_images.m](code/create_decile_images.m) is a subfunction called to compute deciles, which can then be used to perform the shift function analysis.
[statistical_analysis](code/statistical_analysis.m) is the code reading the 'results' and doing the stats and plots.
Note: we developed [modified_DunnIndex](code/modified_DunnIndex.m) which allows to look at tissue separation as for discrete segmentation, but this was not used in the main analysis/report. All the data are, however, available.

## Results

This folder contains the results obtained and saved as csv files. Simply running [statistical_analysis](code/statistical_analysis.m) will redo all the analyses and figures from the publication (and more).



