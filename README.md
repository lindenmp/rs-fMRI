# rfMRI-PrePro
fMRI preprocessing scripts for use with Monash's distributed computing platform, MASSIVE.

The preprocessing pipeline is broken into three functions
- prepro_base.m
- prepro_noise.m
- prepro_extractTS_FSL.m

All three of these are called by the script
- run_prepro.m

run_prepro.m is the script that the user should edit.


prepro_noise allows the user to select from a handful of different noise correction
strategies but each of these methods works off of the output file (except see ICA-AROMA) of the prepro_base script.
The decision to split the scripts was so that users may run all the noise correction options offered
by prepro_noise on their dataset without unnecessarily repeating the base processing steps.
This is desirable because it allows the user to compare how each noise correction strategy compares to the others (see rfMRI-QC repo for more).

Each script requires the user to edit the 'Add Paths' and 'Basic Inputs' sections that occur in both prepro_ scripts

The following is a brief summary of what the two scripts do

## prepro_base.m
    1 - Segment T1 using SPM (Generate tissue masks)

    
    2 - Discard first 4 volumes
    
    3 - Slice-timing correction (optional)
    
    4 - Despiking (optional)
    
    5 - EPI Realignment
    
    6 - Co-registration of realigned images to native T1
    
    7 - Spatially normalize T1 to MNI template
    
    8 - Application of T1-spatial normalization parameters to coregistered EPI
    
    9 - Mask out non-brain voxels
    
    10 - Linear detrending of realigned EPI time series (uses REST)
    
    11 - Modal intensity normalisation (to 1000)
    
    12 - Spatial smoothing


## prepro_noise.m
    1 - Correct noise in EPI output from step 11 of prepro_base script (or output from step 12 in case of ICA-AROMA)
        
    2 - Bandpass filter (includes demeaning)
    
    3 - Spatial smoothing (not in case of ICA-AROMA)	
	
## Usage (matlab)

	subject = '1008.2.48.009'; % <-- string containing subject ID
    
    discard = 1;

    slicetime = 1;
    
    despike = 1;
    
    smoothing = 'after';

    run_prepro(subject,discard,slicetime,despike,smoothing)

## Usage (Monash/MASSIVE users)

    For those users with access to MASSIVE (https://www.massive.org.au), you will two additional .sh scripts

    - MASSIVE_Setup.sh
    - MASSIVE_Job.sh
    
    These are both simple shell scripts that facilitate submitting multiple jobs to MASSIVE that will run concurrently,
    where each job is responsible for processing a single participant in your dataset.

    MASSIVE_Setup.sh reads in a text file containing a list of subject identifiers (e.g., DARIS IDs) and loops over these identifiers,
    each time calling the MASSIVE_Job.sh script and inputting a given subject identifier.

    MASSIVE_Job.sh is where you input all your SBATCH settings (e.g., cpu time, email address, MASSIVE account, etc.) and load the software (i.e., MASSIVE modules).
    MASSIVE_Job.sh also runs matlab and executes the run_prepro.m script for you.

## Testing

This code was developed and tested on the Multi-modal Australian ScienceS Imaging and Visualisation Environment (MASSIVE: https://www.massive.org.au).

This code uses the following software

- MATLAB R2014a (8.3.0.532) 64-bit (glnxa64) February 11, 2014
- SPM8 r5236
- REST 1.8
- FSL 5.0.9
- ANTs 1.9.v4
- ICA-AROMA