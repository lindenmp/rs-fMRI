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
# SBATCH --array=1-100%50

# Assign input args
WhichProject=$1
WhichSessScan=$2

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
cd /gpfs/M2Home/projects/Monash076/Linden/scripts/rfMRI-PrePro/

discard=1
slicetime=1
despike=1
smoothing=after

matlab -nodisplay -r "run_prepro('${WhichProject}','${WhichSessScan}','${subject}',$discard,$slicetime,$despike,'$smoothing'); exit"

# ------------------------------------------------------------------------------
# Compressing outputs
# ------------------------------------------------------------------------------
echo -e "\t\t ----- Compressing outputs ----- \n"

datadir=/gpfs/M2Home/projects/Monash076/Linden/${WhichProject}/data/

if [[ "$WhichProject" == "OCDPG" ]] || [[ "$WhichProject" == "UCLA" ]] || [[ "$WhichProject" == "GoC" ]] ; then t1str=/t1/; preprostr=/rfMRI/prepro/; fi
if [[ "$WhichProject" == "NYU_2" ]]; then if [[ "$WhichSessScan" == "Sess1_Scan1" ]]; then t1str=/session_1/anat_1/; preprostr=/session_1/rest_1/prepro/; elif [[ "$WhichSessScan" == "Sess1_Scan2" ]]; then t1str=/session_1/anat_1/; preprostr=/session_1/rest_2/prepro/; elif [[ "$WhichSessScan" == "Sess2_Scan1" ]]; then t1str=/session_2/anat_1/; preprostr=/session_2/rest_1/prepro/; fi; fi

# T1 dir
cd ${datadir}${subject}${t1str}
count=`ls -1 *.nii 2>/dev/null | wc -l`
if [ $count != 0 ]; then gzip *.nii; fi

# Prepro dir
cd ${datadir}${subject}${preprostr}
count=`ls -1 *.nii 2>/dev/null | wc -l`
if [ $count != 0 ]; then gzip *.nii; fi

# export list of prepro subdirs
find ./ -maxdepth 1 -mindepth 1 -type d -printf '%f\n' >> ${subject}.txt
folders=$(<${subject}.txt)
for i in $folders; do
	if [[ "$i" != "mot" ]]; then cd ${datadir}${subject}${preprostr}${i}; count=`ls -1 *.nii 2>/dev/null | wc -l`; if [ $count != 0 ]; then echo -e "\t\tCompressing ${i} data"; gzip *.nii; fi; fi
done

cd ${datadir}${subject}${preprostr}
rm -f ${subject}.txt

echo -e "\t\t ----- Finished. ----- \n"