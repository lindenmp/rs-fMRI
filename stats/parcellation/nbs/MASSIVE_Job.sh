#!/bin/env bash

#SBATCH --job-name=fMRI-WholeBrain
#SBATCH --account=monash076
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=16G

# Assign input args
WhichProject=$1
WhichSplit=$2
WhichParc=$3
P=$4

# MASSIVE modules
module load spm8/matlab2014a.r6685

matlab -nodisplay -r "ComputeWholeBrainDiff('${WhichProject}','${WhichSplit}','${WhichParc}',${P}); exit"
