#!/bin/bash

for WhichProject in Beijing_Zang UCLA; do
	SUBJIDS=$(</home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/${WhichProject}/${WhichProject}_SubjectIDs.txt)
	BASEINDIR=/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/${WhichProject}/
	BASEOUTDIR=/home/lindenmp/kg98_scratch/Linden/ResProjects/figshare/rfMRI_denoise/${WhichProject}/

	mkdir -p ${BASEOUTDIR}
	# copy CSV
	cp ${BASEINDIR}${WhichProject}.csv ${BASEOUTDIR}

	INDIR=${BASEINDIR}data/
	OUTDIR=${BASEOUTDIR}data/

	for SUBJ in $SUBJIDS; do
		echo -e "Copying $SUBJ \n"

		STR=/func/prepro/

		for i in 6P 6P+2P 6P+2P+GSR 24P 24P+8P 24P+8P+4GSR 24P+aCC 24P+aCC+4GSR 24P+aCC50 24P+aCC50+4GSR 12P+aCC 12P+aCC50 ICA-AROMA+2P ICA-AROMA+2P+GSR ICA-AROMA+8P ICA-AROMA+8P+4GSR 24P+8P+4GSR+SpikeReg; do
			mkdir -p ${OUTDIR}${SUBJ}${STR}${i}
			cp ${INDIR}${SUBJ}${STR}${i}/cfg*.mat ${OUTDIR}${SUBJ}${STR}${i}/
		done

		mkdir -p ${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub
		cp ${INDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub/cfg*.mat ${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub
		if test "$(ls -A "${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub")"; then
			echo "The dir isnt empty"
		else
			echo "The dir is empty. deleting."
		    rmdir ${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub
		fi

		# Copy .txt files from ICA-AROMA
		mkdir -p ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output
		cp ${INDIR}${SUBJ}${STR}ICA-AROMA_output/classification_overview.txt ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output/
		cp ${INDIR}${SUBJ}${STR}ICA-AROMA_output/classified_motion_ICs.txt ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output/
		cp ${INDIR}${SUBJ}${STR}ICA-AROMA_output/feature_scores.txt ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output/

		# Also copy the mot files
		cp ${INDIR}${SUBJ}${STR}rp_*.txt ${OUTDIR}${SUBJ}${STR}
		cp -r ${INDIR}${SUBJ}${STR}raw_mov ${OUTDIR}${SUBJ}${STR}raw_mov

		# Also copy the dvars.txt file
		cp ${INDIR}${SUBJ}${STR}dvars.txt ${OUTDIR}${SUBJ}${STR}

		# Copy power 2014 scrub mask
		cp ${INDIR}${SUBJ}${STR}JP14_ScrubMask.txt ${OUTDIR}${SUBJ}${STR}
		cp ${INDIR}${SUBJ}${STR}fdPower.txt ${OUTDIR}${SUBJ}${STR}

		# Get numcomponents for aCompCor50 pipelines. This will be the same across aCompCor50 pipelines.
		cp ${INDIR}${SUBJ}${STR}12P+aCC50/aCC_num_wm.txt ${OUTDIR}${SUBJ}${STR}12P+aCC50/
		cp ${INDIR}${SUBJ}${STR}12P+aCC50/aCC_num_csf.txt ${OUTDIR}${SUBJ}${STR}12P+aCC50/
	done

	if test "$(ls -A "${BASEINDIR}NBS_liberal")"; then
		cp ${BASEINDIR}NBS_liberal ${BASEOUTDIR}
	fi

	if test "$(ls -A "${BASEINDIR}NBS_stringent")"; then
		cp ${BASEINDIR}NBS_stringent ${BASEOUTDIR}
	fi
	
done

##############################################################################################################################
# Get NYU
##############################################################################################################################

WhichProject=NYU_2
SUBJIDS=$(</home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/${WhichProject}/${WhichProject}_SubjectIDs.txt)
BASEINDIR=/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/${WhichProject}/
BASEOUTDIR=/home/lindenmp/kg98_scratch/Linden/ResProjects/rfMRI_denoise/figshare/${WhichProject}/

mkdir -p ${BASEOUTDIR}
# copy CSV
cp ${BASEINDIR}${WhichProject}.csv ${BASEOUTDIR}

INDIR=${BASEINDIR}data/
OUTDIR=${BASEOUTDIR}data/

for SUBJ in $SUBJIDS; do
	echo -e "Copying $SUBJ \n"

	for STR in /session_1/func_1/prepro/ /session_1/func_2/prepro/ /session_2/func_1/prepro/; do

		for i in 6P 6P+2P 6P+2P+GSR 24P 24P+8P 24P+8P+4GSR 24P+aCC 24P+aCC+4GSR 24P+aCC50 24P+aCC50+4GSR 12P+aCC 12P+aCC50 ICA-AROMA+2P ICA-AROMA+2P+GSR ICA-AROMA+8P ICA-AROMA+8P+4GSR 24P+8P+4GSR+SpikeReg; do
			mkdir -p ${OUTDIR}${SUBJ}${STR}${i}
			cp ${INDIR}${SUBJ}${STR}${i}/cfg*.mat ${OUTDIR}${SUBJ}${STR}${i}/
		done

		mkdir -p ${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub
		cp ${INDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub/cfg*.mat ${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub
		if test "$(ls -A "${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub")"; then
			echo "The dir isnt empty"
		else
			echo "The dir is empty. deleting."
		    rmdir ${OUTDIR}${SUBJ}${STR}24P+4P+2GSR+JP14Scrub
		fi

		# Copy .txt files from ICA-AROMA
		mkdir -p ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output
		cp ${INDIR}${SUBJ}${STR}ICA-AROMA_output/classification_overview.txt ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output/
		cp ${INDIR}${SUBJ}${STR}ICA-AROMA_output/classified_motion_ICs.txt ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output/
		cp ${INDIR}${SUBJ}${STR}ICA-AROMA_output/feature_scores.txt ${OUTDIR}${SUBJ}${STR}ICA-AROMA_output/

		# Also copy the mot files
		cp ${INDIR}${SUBJ}${STR}rp_*.txt ${OUTDIR}${SUBJ}${STR}
		cp -r ${INDIR}${SUBJ}${STR}raw_mov ${OUTDIR}${SUBJ}${STR}raw_mov

		# Also copy the dvars.txt file
		cp ${INDIR}${SUBJ}${STR}dvars.txt ${OUTDIR}${SUBJ}${STR}

		# Copy power 2014 scrub mask
		cp ${INDIR}${SUBJ}${STR}JP14_ScrubMask.txt ${OUTDIR}${SUBJ}${STR}
		cp ${INDIR}${SUBJ}${STR}fdPower.txt ${OUTDIR}${SUBJ}${STR}

		# Get numcomponents for aCompCor50 pipelines. This will be the same across aCompCor50 pipelines.
		cp ${INDIR}${SUBJ}${STR}12P+aCC50/aCC_num_wm.txt ${OUTDIR}${SUBJ}${STR}12P+aCC50/
		cp ${INDIR}${SUBJ}${STR}12P+aCC50/aCC_num_csf.txt ${OUTDIR}${SUBJ}${STR}12P+aCC50/
	done
done

