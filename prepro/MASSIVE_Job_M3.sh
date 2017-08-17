#!/bin/env bash

#SBATCH --job-name=fMRI-PrePro
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-100%25

# Assign input args
WhichMASSIVE=$1
WhichProject=$2
WhichSessScan=$3

case $WhichProject in
	M3_OCDPG) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/OCDPG/OCDPG_SubjectIDs.txt" ;;
	M3_COBRE) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/COBRE/COBRE_SubjectIDs.txt" ;;
	M3_UCLA) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/UCLA/UCLA_SubjectIDs.txt" ;;
	M3_NAMIC) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/NAMIC/NAMIC_SubjectIDs.txt" ;;
esac

# MASSIVE modules
module load rest/1.8 # works on M2 and M3
module load fsl/5.0.9 # works on M2 and M3. FSLDIR is different though.
module load python/2.7.11-gcc # works on M2 and M3.
module load ants/1.9.v4 # works on M2 and M3
module load afni/16.2.16
module load spm8/matlab2015b.r6685 # M3

subject=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SUBJECT_LIST})
echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"

# Run pipeline
cd /home/lindenmp/kg98/Linden/Scripts/rs-fMRI/prepro/

smoothing=before
discard=1
slicetime=1
despike=1
detr=1
intnorm=1

matlab -nodisplay -r "run_prepro('${WhichMASSIVE}','${WhichProject}','${WhichSessScan}','${subject}','$smoothing',$discard,$slicetime,$despike,$detr,$intnorm); exit"
