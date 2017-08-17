#!/bin/bash
# This script downloads dicoms from XNAT and convert them to nifti using dcm2niix
# Nifti files will be named and organised according to BIDS

# ------------------------------------------------------------------------------
# Linden Parkes, James Morrow, Kristina Sabaroedin Brain & Mental Health Laboratory, 2017
# ------------------------------------------------------------------------------

# Your project paths on your local machine/MASSIVE
PROJDIR=/projects/kg98/Linden/ResProjects/OCDPG
DICOMDIR=$PROJDIR/dicomdir/
NIFTIDIR=$PROJDIR/data/

# These variables is based on the name of your project and scans on XNAT
# Please edit them accordingly
# Often, scans are named differently on XNAT
# After running this script, please check whether all scans are there
# If they are missing, go back on XNAT and see if the scans are named differently
ANATOMICAL=t1_mprage_sag_p2_iso_1_ADNI
FUNCTIONAL=Resting_ep2d_p2_3mm
STUDY=MRH034_
SESSION=_MR01

# Text file containing subject IDs
# These IDs just need to be the last three digits (zero padded) i.e. 007, 098, 231, etc
SUBJIDS=$(</projects/kg98/Linden/ResProjects/OCDPG/xnat_SubjectIDs.txt)
# SUBJIDS=009

# quirky subjects on xnat
# SUBJIDS=107
# SESSION=_MR03
# SUBJIDS=110
# SESSION=_MR02

# load modules
module purge;
module load xnat-utils;
# Load the dcm2niix software
module load mricrogl/1.0.20170207

# create for loop to loop over IDs
for ID in $SUBJIDS; do 
	
	# Dynamic directories
	SUBDICOMDIR=${DICOMDIR}sub-$ID/
	OUTDIR=${NIFTIDIR}sub-$ID/
	EPIOUTDIR=${OUTDIR}func/
	T1OUTDIR=${OUTDIR}anat/

	# Create subject's DICOMS folders 
	if [ ! -d ${SUBDICOMDIR} ]; then mkdir -p ${SUBDICOMDIR}; fi
	if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; fi
	if [ ! -d ${EPIOUTDIR} ]; then mkdir -p ${EPIOUTDIR}; fi
	if [ ! -d ${T1OUTDIR} ]; then mkdir -p ${T1OUTDIR}; fi

    # ------------------------------------------------------------------------------
	# T1
    # ------------------------------------------------------------------------------
	# Download structural scans from XNAT
	cd $SUBDICOMDIR
	xnat-get ${STUDY}${ID}${SESSION} --scans $ANATOMICAL
	
	# Delete intermediate folders
	mv ${SUBDICOMDIR}${STUDY}${ID}${SESSION}/*$ANATOMICAL ${SUBDICOMDIR}

	# Run dcm2niix
	dcm2niix -f "sub-"$ID"_T1w" -o $T1OUTDIR -b -m n -z y $SUBDICOMDIR*$ANATOMICAL/

    # ------------------------------------------------------------------------------
    # rs-fMRI EPI
    # ------------------------------------------------------------------------------
	# Download resting state scans from XNAT
	cd $SUBDICOMDIR
	xnat-get ${STUDY}${ID}${SESSION} --scans $FUNCTIONAL
	
	# Delete intermediate folders
	mv ${SUBDICOMDIR}${STUDY}${ID}${SESSION}/*$FUNCTIONAL ${SUBDICOMDIR}

	# Run dcm2niix
	dcm2niix -f "sub-"$ID"_task-rest_bold" -o $EPIOUTDIR -b -m n -z y $SUBDICOMDIR*$FUNCTIONAL/

    # ------------------------------------------------------------------------------
    # Clean up
	# ------------------------------------------------------------------------------
	rm -rf ${SUBDICOMDIR}${STUDY}${ID}${SESSION}

done
