#!/bin/env bash

#SBATCH --job-name=spDCM
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=4G
#SBATCH --qos=shortq
# SBATCH --qos=rtq
#SBATCH -A kg98

#SBATCH --array=1-90%30
# SBATCH --array=1-353%30

# Assign input args
ProjDir=/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_DCM/OCDPG/
DataDir=${ProjDir}data/
WhichNoise=ICA-AROMA+2P
# WhichNoise=ICA-AROMA+2P+GSR
data=epi_prepro.nii
units=scans
N=185
TR=2.5
TE=0.03

# ProjDir=/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_DCM/GenCog/
# DataDir=${ProjDir}data/
# WhichNoise=ICA-FIX
# # WhichNoise=ICA-FIX+GSR
# data=epi_prepro.nii
# units=scans
# N=616
# TR=0.754
# TE=0.021

WhichSeed=TriStri
# WhichSeed=DiMartino

SUBJECT_LIST="${ProjDir}SecondLevel/SPM/Factorial/"${WhichNoise}"/${WhichSeed}/ParticipantIDs.txt"
subject=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SUBJECT_LIST})
datadir=${DataDir}${subject}/func/prepro/${WhichNoise}/

spmdir=/home/lindenmp/kg98/Linden/Scripts/Tools/spm12/

# MASSIVE modules
module load matlab/r2017a

echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"

for i in L R; do
	##############
	# Generate nodes
	##############
	# DORSAL
	cd /home/lindenmp/kg98/Linden/Scripts/rs-fMRI/stats/spDCM/
	if [ "$WhichSeed" = "TriStri" ]; then WhichStri=2; elif [ "$WhichSeed" = "DiMartino" ]; then WhichStri=3; fi
	searchVolDir=${ProjDir}SecondLevel/spDCM/SearchVolumes/${WhichNoise}/${WhichSeed}_${WhichStri}_${i}/
	con1Name=${DataDir}${subject}/func/prepro/${WhichNoise}/FirstLevel/${WhichSeed}_${i}/con_000${WhichStri}.nii
	seedName=${WhichSeed}_${WhichStri}_${i}.nii
	nodeDir=${DataDir}${subject}/func/prepro/${WhichNoise}/Nodes/${WhichSeed}_${WhichStri}_${i}/
	matlab -nodisplay -r "spDCM_GenerateSubjectNodes('${searchVolDir}','${con1Name}','${seedName}','${nodeDir}'); exit"

	# VENTRAL
	cd /home/lindenmp/kg98/Linden/Scripts/rs-fMRI/stats/spDCM/
	if [ "$WhichSeed" = "TriStri" ]; then WhichStri=3; elif [ "$WhichSeed" = "DiMartino" ]; then WhichStri=1; fi
	searchVolDir=${ProjDir}SecondLevel/spDCM/SearchVolumes/${WhichNoise}/${WhichSeed}_${WhichStri}_${i}/
	con1Name=${DataDir}${subject}/func/prepro/${WhichNoise}/FirstLevel/${WhichSeed}_${i}/con_000${WhichStri}.nii
	seedName=${WhichSeed}_${WhichStri}_${i}.nii
	nodeDir=${DataDir}${subject}/func/prepro/${WhichNoise}/Nodes/${WhichSeed}_${WhichStri}_${i}/
	matlab -nodisplay -r "spDCM_GenerateSubjectNodes('${searchVolDir}','${con1Name}','${seedName}','${nodeDir}'); exit"

	##############
	# run spDCM
	##############
	cd /home/lindenmp/kg98/Linden/Scripts/rs-fMRI/stats/spDCM/
	if [ "$WhichSeed" = "TriStri" ]; then
		dorsalStr='TriStri_2' # TriStri
		ventralStr='TriStri_3' # TriStri
		nodedir1=${DataDir}${subject}/func/prepro/${WhichNoise}/Nodes/${WhichSeed}_2_${i}/ # TriStri
		nodedir2=${DataDir}${subject}/func/prepro/${WhichNoise}/Nodes/${WhichSeed}_3_${i}/ # TriStri

		out=spDCM/${WhichSeed}_23_${i}/ # TriStri
		# out=spDCM/${WhichSeed}_23_${i}_alt/ # TriStri
	elif [ "$WhichSeed" = "DiMartino" ]; then
		dorsalStr='DiMartino_3' # DiMartino
		ventralStr='DiMartino_1' # DiMartino
		nodedir1=${DataDir}${subject}/func/prepro/${WhichNoise}/Nodes/${WhichSeed}_3_${i}/ # DiMartino
		nodedir2=${DataDir}${subject}/func/prepro/${WhichNoise}/Nodes/${WhichSeed}_1_${i}/ # DiMartino

		out=spDCM/${WhichSeed}_31_${i}/ # DiMartino
		# out=spDCM/${WhichSeed}_31_${i}_alt/ # DiMartino
	fi

	matlab -nodisplay -r "spDCM_FirstLevel('${spmdir}','${datadir}','${data}','${units}',${N},${TR},${TE},{'${nodedir1}','${nodedir2}'},'${out}',{'${dorsalStr}','${ventralStr}'}); exit"
done
