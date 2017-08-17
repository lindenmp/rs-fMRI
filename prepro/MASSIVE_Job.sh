#!/bin/env bash

#SBATCH --job-name=fMRI-PrePro
#SBATCH --account=monash076
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-100%50

# Assign input args
WhichMASSIVE=$1
WhichProject=$2
WhichSessScan=$3

case $WhichProject in
	OCDPG) SUBJECT_LIST="/gpfs/M2Home/projects/Monash076/Linden/Sublists/OCDPG.txt" ;;
	UCLA) SUBJECT_LIST="/gpfs/M2Home/projects/Monash076/Linden/Sublists/UCLA.txt" ;;
	NYU_2) SUBJECT_LIST="/gpfs/M2Home/projects/Monash076/Linden/Sublists/NYU_2.txt" ;;
	GoC) SUBJECT_LIST="/projects/kg98/kristina/code/sublists/GoC.txt" ;;
esac

# MASSIVE modules
module load rest/1.8
module load fsl/5.0.9
module load python/2.7.11-gcc
module load ants/1.9.v4
module load afni/16.2.16
module load spm8/matlab2014a.r5236

subject=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SUBJECT_LIST})
echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"

# Run pipeline
cd /gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/prepro/

smoothing=before
discard=1
slicetime=1
despike=1
detr=1
intnorm=1

matlab -nodisplay -r "run_prepro('${WhichMASSIVE}','${WhichProject}','${WhichSessScan}','${subject}','$smoothing',$discard,$slicetime,$despike,$detr,$intnorm); exit"
