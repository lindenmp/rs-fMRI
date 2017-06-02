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
# SBATCH --array=1-2

# Assign input args
WhichMASSIVE=$1
WhichProject=$2
WhichSessScan=$3

case $WhichProject in
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

discard=1
slicetime=1
despike=1
smoothing=before

matlab -nodisplay -r "run_prepro('${WhichMASSIVE}','${WhichProject}','${WhichSessScan}','${subject}',$discard,$slicetime,$despike,'$smoothing'); exit"

# ------------------------------------------------------------------------------
# Compressing outputs
# ------------------------------------------------------------------------------
echo -e "\t\t ----- Compressing outputs ----- \n"

datadir=/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/${WhichProject:3}/data/

if [[ "$WhichProject" == "M3_COBRE" ]]; then t1str=/session_1/anat_1/; preprostr=/session_1/rest_1/prepro/; fi
if [[ "$WhichProject" == "M3_UCLA" ]] || [[ "$WhichProject" == "M3_NAMIC" ]]; then t1str=/anat/; preprostr=/func/prepro/; fi

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
