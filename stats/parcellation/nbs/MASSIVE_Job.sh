#!/bin/env bash

#SBATCH --job-name=fMRI-WholeBrain
#SBATCH --account=monash076
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=16G

# Assign input args
WhichProject=$1
WhichSplit=$2
WhichParc=$3
WhichNoise=$4
P=$5
outDir=$6
runCensor=$7

# MASSIVE modules
module load matlab/r2014a

matlab -nodisplay -r "ComputeWholeBrainDiff('${WhichProject}','${WhichSplit}','${WhichParc}','${WhichNoise}',${P},'${outDir}','${runCensor}'); exit"
