#!/bin/env bash

#SBATCH --job-name=fMRI-HCP-PrePro
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=64G
#SBATCH --array=1-100%50

# SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/HCP_BF_TimeSeries/100_subs.txt"
SUBJECT_LIST="/home/lindenmp/kg98/Linden/ResProjects/HCP_BF_TimeSeries/HCP_100.txt"

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
cd /home/lindenmp/kg98/Linden/Scripts/rs-fMRI/prepro-hcp/

matlab -nodisplay -r "HCP_run_prepro('${subject}'); exit"
