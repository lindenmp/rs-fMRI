# Resting-state functional magnetic resonance imaging
*rs-fMRI* is a repository of code for preprocessing, denoising, and running quality control on resting-state fMRI datasets, using Matlab.

Resting-state fMRI is highly sensitive to artefacts caused by in-scanner movement. These artefacts can cause spurious correlations in the time-series data that hinder functional connectivity analyses. This code applies a wide range of popular data denoising methods designed to correct for these artefacts and produces a range of quality control benchmarks in order to facilitate the evaluation of their efficacy.

This code is broken into four main subdirectories:
1. **prepro**: process and denoising rs-fMRI data on [MASSIVE](https://www.massive.org.au).
2. **qc**: calculates quality control benchmarks and reproduces figures found in the following [manuscript](https://www.sciencedirect.com/science/article/pii/S1053811917310972)
3. **stats**: various second-level statistics
4. **func**: functions for repository

The details of each subdirectory can be found below.

Linden Parkes, linden.parkes@monash.edu, Brain & Mental Health Laboratory, Monash University

# Publications

See the following publications for examples of this code in use:
- **An evaluation of the efficacy, reliability, and sensitivity of motion correction strategies for resting-state functional MRI.** L. Parkes, B. D. Fulcher, M. Yucel, & A. Fornito. [*NeuroImage*](https://www.sciencedirect.com/science/article/pii/S1053811917310972) (2017).

# Data

Fully processed rs-fMRI time-series data for two of the three datasets (*CNP* and *NYU*) used in the above preprint are available for download from [Figshare](https://doi.org/10.4225/03/595193482c03e).
Raw, unprocessed data is also available for [CNP](https://openfmri.org/dataset/ds000030/) and [NYU](http://fcon_1000.projects.nitrc.org/indi/CoRR/html/)

## prepro

Code in the ***prepro*** subdirectory was developed and tested on the Multi-modal Australian ScienceS Imaging and Visualisation Environment, [MASSIVE](https://www.massive.org.au).
This code uses the following software:
- MATLAB R2014a (8.3.0.532) 64-bit (glnxa64) February 11, 2014
- SPM8 r5236
- REST 1.8
- FSL 5.0.9
- ANTs 1.9.v4
- ICA-AROMA

These scripts are hard coded for use with my own projects but can be edited to run on new data.

The preprocessing pipeline is broken into three main scripts
- `prepro_base.m`
- `prepro_noise.m`
- `prepro_extractTS_FSL.m`

All three of these are called by the script `run_prepro.m`.
This is the script that the user should edit.

`run_prepro.m` allows the user to select from a handful of different noise correction strategies but each of these methods works off of the output file (except ICA-AROMA) of `prepro_base`.
The decision to split the scripts was so that users may run all the noise correction options available on their dataset without unnecessarily repeating the base processing steps.
This is desirable because it allows the user to examine how each noise correction strategy compares to the others (see the **qc**).

The following is a brief summary of what the two scripts do

#### prepro_base.m
1. Segment T1 using SPM and erode the generated tissue masks
2. Discard first 4 volumes (optional. set using input argument)
3. Slice-timing correction (optional. set using input argument)
4. Despiking (optional. set using input argument)
5. EPI Realignment
6. Co-registration of realigned images to native T1
7. Spatially normalize T1 to MNI template
8. Application of T1-spatial normalization parameters to coregistered EPI and tissue masks
9. Mask out non-brain voxels
10. Linear detrending of realigned EPI time series  (optional. set using input argument)
11. Modal intensity normalisation (to 1000) (optional. set using inputargument)
12. Spatial smoothing (for ICA-AROMA)

#### prepro_noise.m
1. Correct noise in EPI output from step 11 `prepro_base` (or output from step 12 in case of ICA-AROMA)
2. Bandpass filter
3. Spatial smoothing (not in case of ICA-AROMA)

#### Basic usage (Matlab)

```matlab
>> WhichProject = 'UCLA';
>> WhichSessScan = 'Sess1_Scan1';
>> subject = 'sub-10159'; % <-- string containing subject ID
>> run_prepro(WhichMASSIVE,WhichProject,WhichSessScan,subject)
```

#### Usage (MASSIVE)
- `$ ./MASSIVE_Job.sh`

This script submits an array slurm job for processing multiple participants at once on M3.

## Reproducibility

The code in the **qc** subdirectory can be used to reproduce the figures in **An evaluation of the efficacy, reliability, and sensitivity of motion correction strategies for resting-state functional MRI.** L. Parkes, B. D. Fulcher, M. Yucel, & A. Fornito. [*NeuroImage*](https://www.sciencedirect.com/science/article/pii/S1053811917310972) (2017).

#### QC.m

This script computes the quality control benchmarks outlined in the above manuscript. Below are step by step instructions to produce the primary plots from Figure 1 (QC-FC) and Figure 2 (QC-FC distance-dependence) using MAC OS X.

1. Download the full *rs-fMRI* repository onto your computer (e.g., ~/Desktop/rs-fMRI/)
2. Download the fully processed data (approx. 30GB) from [Figshare](https://doi.org/10.4225/03/595193482c03e) and unpack onto desktop (~/Desktop/)
3. Download the Gordon parcels by visiting [here](http://www.nil.wustl.edu/labs/petersen/Resources.html) and downloading the Parcels from [Gordon et al., Cerebral Cortex](https://www.ncbi.nlm.nih.gov/pubmed/25316338).
4. Place **Parcels_MNI_222.nii** (and **Gordon_Centroids.txt** and **CommunityModified.txt** from [Figshare](https://doi.org/10.4225/03/595193482c03e)) into ~/Desktop/ROIs/Gordon/.
5. Load Matlab (tested on version version 2014a, 2015b, and 2017b) and navigate to the **qc** subdirectory - ```cd ~/Desktop/rs-fMRI/qc```
6. Run `QC.m`.

By default, this will produce plots for the Beijing Zang dataset.
If you want to reproduce the CNP dataset then change line 40 to:
```WhichProject = Projects{2};```
This will produce the plots for Figure 1b and Figure 2b (i.e., CNP under lenient exclusion criteria).
If you want to select stringent criteria then you will need to change ```WhichExclude``` to ```2``` (see lines 479 and 481).

If you have any questions or run into issues please let me know by email (linden.parkes@monash.edu).

## QC_ThePlot.m

`QC_ThePlot.m` produces a *carpet plot* (see [Power 2016, NeuroImage](https://doi.org/10.1016/j.neuroimage.2016.08.009)) so that the user can perform subject-level quality control.
This script is also hard coded for use with my own projects but can be edited to run on new data.

## stats

The **stats** subdirectory contains the code used to compute whole-brain differences in functional connectivity using the Network Based Statistic ([Zalesky et al., 2014. NeuroImage](https://doi.org/10.1016/j.neuroimage.2010.06.041))
