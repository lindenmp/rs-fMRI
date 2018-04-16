#!/bin/env bash

#SBATCH --job-name=fMRI-PrePro
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-100

# Assign input args
WhichProject=$1
WhichSessScan=$2

case $WhichProject in
	OCDPG_DCM) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_DCM/OCDPG/OCDPG_SubjectIDs.txt" ;;
	OCDPG) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/OCDPG/OCDPG_SubjectIDs.txt" ;;
	UCLA) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/UCLA/UCLA_SubjectIDs.txt" ;;
	NYU_2) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/NYU_2/NYU_2_SubjectIDs.txt" ;;
	COBRE) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/COBRE/COBRE_SubjectIDs.txt" ;;
	TBS_fMRI) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/TBS_fMRI/TBS_fMRI_SubjectIDs.txt" ;;
	Beijing_Zang) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/Beijing_Zang/Beijing_Zang_SubjectIDs.txt" ;;
	COBRE_HCTSA) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/COBRE/COBRE_SubjectIDs.txt" ;;
	UCLA_HCTSA) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/UCLA/UCLA_SubjectIDs.txt" ;;
	NAMIC_HCTSA) SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/NAMIC/NAMIC_SubjectIDs.txt" ;;
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

smoothing=after
matlab -nodisplay -r "run_prepro('${WhichProject}','${WhichSessScan}','${subject}','$smoothing'); exit"
