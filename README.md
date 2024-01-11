# Is there a neural common factor for visual illusions?
Maya A. Jastrzębowska, Ayberk Ozkirli, Aline F. Cretenoud, Bogdan Draganski, Michael H. Herzog
https://www.biorxiv.org/content/10.1101/2023.12.27.573437v1

The following Matlab and JASP files can be used to recreate the analyses and figures from our preprint "Is there a neural common factor for visual illusions?", in which we tested the association between illusion strengths in a battery of 13 illusions with the size of visual areas and population receptive field (pRF) sizes in a group of 30 human participants. Please see the preprint for the description and acronym labels of the 13 illusions and additional details. The data required to carry out the analyses is contained in the following files:
- ill_fMRI_illusionMagnitudes.csv: illusion magnitudes for each of the 13 tested illusions, for each of the 30 participants
- surface_areas_0_8.mat: surface areas of visual regions of interest (ROIs) V1, V2, V3 and V4 in left and right hemispheres for each of the 30 participants
- pRF_measures.mat: table of pRF sizes (sigmas) and corresponding goodness-of-fit (R^2) values across eccentricities, for each of the 30 participants, ROI and hemisphere

## Abstract
It is tempting to map interindividual variability in human perception to variability in brain structure or neural activity. Indeed, it has been shown that susceptibility to size illusions correlates with the size of primary visual cortex V1. Yet contrary to common belief, illusions correlate only weakly at the perceptual level, raising the question of how they can correlate with a localized neural measure. In addition, mounting evidence suggests that there is substantial interindividual variability not only in neural function and anatomy but also in the mapping between the two, which further challenges the findings of a neural common factor for illusions. To better understand these questions, here, we re-evaluated previous studies by correlating illusion strengths in a battery of 13 illusions with the size of visual areas and population receptive field sizes. We did not find significant correlations either at the perceptual level or between illusion susceptibility and visual functional neuroanatomy.

## Cortical idiosyncrasies and illusion magnitude
### Visual surface area
correlate_illMagn_surfaceArea.m

illMagn_surfaceArea_0_8.jasp

### PRF size and illusion magnitude
correlate_illMagn_pRFsize.m

illMagn_pRF_size.jasp

### Slope and intercept of pRF size as a function of eccentricity 
correlate_illMagn_pRFsizeSlopeIntercept.m

illMagn_pRFsizeSlopeIntercept.jasp
