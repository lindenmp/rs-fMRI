#!/bin/env bash

#SBATCH --job-name=FreeSurf
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8G

# OCDPG
# SBATCH --array=1-100

# UCLA
#SBATCH --array=1-186

# NYU_2
# SBATCH --array=1-29

# COBRE
# SBATCH --array=1-147

# Beijing_Zang
# SBATCH --array=1-192

# Assign input args
WhichProject=UCLA

case $WhichProject in
	OCDPG) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/OCDPG/OCDPG_SubjectIDs.txt" ;;
	UCLA) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/UCLA/UCLA_SubjectIDs.txt" ;;
	NYU_2) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/NYU_2/NYU_2_SubjectIDs.txt" ;;
	COBRE) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/COBRE/COBRE_SubjectIDs.txt" ;;
	Beijing_Zang) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/Beijing_Zang/Beijing_Zang_SubjectIDs.txt" ;;
esac

SUBJECTID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SUBJECT_LIST})
echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} ${SUBJECTID} ----- "
echo -e "\t\t\t --------------------------- \n"

# ------------------------------------------------------------------------------
echo -e "\nSETTING UP FREESURFER\n"
module load freesurfer/6.0
module load fsl/5.0.9
module load ants/1.9.v4 # works on M2 and M3

export WHERESMYSCRIPT="/projects/kg98/Linden/Scripts/rs-fMRI/Freesurfer"
export WHEREERRORLOG="/projects/kg98/Linden/Scripts/rs-fMRI/Freesurfer/errorLog"
export PARENTDIR="/scratch/kg98/Linden/ResProjects/rfMRI_denoise/${WhichProject}/data"
export FS_SUBJECTS_DIR="/scratch/kg98/Linden/ResProjects/rfMRI_denoise/${WhichProject}/data"

export FREESURFERDIR="/usr/local/freesurfer/6.0"
export FSLDIR="/usr/local/fsl/5.0.9"
export FREESURFER_HOME="${FREESURFERDIR}"
export FSFAST_HOME="${FREESURFERDIR}/fsfast"
export FSF_OUTPUT_FORMAT="nii"
export SUBJECTS_DIR="${FS_SUBJECTS_DIR}"
export FUNCTIONALS_DIR="${FREESURFERDIR}/sessions"
export MNI_DIR="${FREESURFERDIR}/mni"
export FSL_DIR="${FSLDIR}"
export FSLOUTPUTTYPE="NIFTI"
source ${WHERESMYSCRIPT}/FreeSurferEnv.sh

module load spm8/matlab2015b.r6685 # M3
export SPMDIR="/usr/local/spm8/matlab2015b.r6685/"

# ------------------------------------------------------------------------------
SUBJECTDIR=${SUBJECTS_DIR}/${SUBJECTID}/
PREPRODIR=${SUBJECTDIR}func/prepro/
RIBBONDIR=${SUBJECTDIR}mri/
WORKDIR=${SUBJECTDIR}ribbon2GM/

if [ -d "$WORKDIR" ]; then
	echo "Cleaning and re-initialising outdir"
	rm -r ${WORKDIR}
	mkdir -p ${WORKDIR}
elif [ ! -d "$WORKDIR" ]; then
	echo "Initialising outdir"
	mkdir -p ${WORKDIR}
fi
cd ${WORKDIR}

# Convert to ribbon.nii.gz
mri_convert ${RIBBONDIR}ribbon.mgz ribbon.nii.gz

# Create GM mask
fslmaths ribbon.nii.gz -thr 3 -uthr 3 -bin gm_l.nii.gz
fslmaths ribbon.nii.gz -thr 42 -uthr 42 -bin gm_r.nii.gz
fslmaths gm_l.nii.gz -add gm_r.nii.gz ribbon.nii.gz

# Transform to MNI space
antsApplyTransforms -d 3 -e 0 -i ribbon.nii.gz -r ${SPMDIR}templates/T1.nii -o ribbon_mni.nii.gz -n NearestNeighbor -t ${PREPRODIR}t12MNI_1Warp.nii.gz -t ${PREPRODIR}t12MNI_0GenericAffine.mat

# Get GM time series
EPI=dibwrat${SUBJECTID}_task-rest_bold.nii.gz
fslmeants -i ${PREPRODIR}${EPI} -o ribbonTS.txt -m ribbon_mni.nii.gz
EPI=sdibwrat${SUBJECTID}_task-rest_bold.nii.gz
fslmeants -i ${PREPRODIR}${EPI} -o sribbonTS.txt -m ribbon_mni.nii.gz

# ------------------------------------------------------------------------------
