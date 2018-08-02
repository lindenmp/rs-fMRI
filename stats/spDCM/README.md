# spDCM
<!-- **spDCM** contains code used to perform spectral DCM analysis found in [link to preprint]. -->

Before running this code, `FactorialSPM.m` and `SeedBased_func.m` from the ***stats/seed*** subdirectory were run using the TriStri parcellation found on [Figshare](https://figshare.com/articles/TriStri_nii/4903118)

# Publications

<!-- See the following publications for examples of this code in use: -->
<!-- - **[name]** L. Parkes, ..... [*journal*](link) (2018). -->

# Data

<!-- Fully processed rs-fMRI EPI data used in the above preprint are available for download from [Figshare](link). -->

# Scripts

Below is a summary of what each of the scripts in this subdirectory do with regard to the analysis reported above in Parkes et al.

## Generate search volumes
- `spDCM_GenerateSearchVolumes.m`

Generates search volumes used to constrain the generation of subject-specific DCM nodes.
For each of two striatal seeds (dorsal, ventral) and hemispheres, this script performs the following steps:
1. Generates anatomical masks by combining specific parcels in the AAL parcellation.
2. Removes voxels from these masks that are within 20mm of the seed.
3. Find the peak t-value from the second-level main effect of seed within each mask.
4. Generates spherical search volumes for each mask centred on the peak t-value.
5. Contrains each search volume to within boundaries of corresponding anatomical mask.
6. For search volumes that are present in both dorsal and ventral striatal circuits, subtracts overlap

The outcome is a series of search volumes that are informed a priori by anatomy and refined using seed-based functional connectivity analyses.
These can be viewed for each strial subregion is Figures 1 and S2 in Parkes et al.

## Generate subject-specific nodes for DCM
- `spDCM_GenerateSubjectNodes.m`

Uses the search volumes created above to guide the generation of subject-specific regions of interest (ROI) for DCM.
For each subject, spherical ROIs are centred on the peak seed-based functional connectivity within each search volume.

## First Level DCM
- `spDCM_FirstLevel.m`

Takes processed EPI data and the ROIs generated above to perform first-level spectral DCM analysis, including specification and full inversion of a sparsely connected parent and generation of nested DCMs, wherein each connection in the parent is systemtically switched off.

## Second Level DCM
- `spDCM_SecondLevel.m`

Takes first-level DCMs across participants and performs second-level Bayesian analysis of individual differences in disinhibition, compulsivity, and impulsivity, as well as case-control/case-case comparisons.
