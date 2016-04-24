#!/bin/sh

#  generateSPORKfasta.sh
#  
#
#  Created by Gillian Hsieh on 4/19/16.
#

KNIFEOUTPUTDIR=${1}
SPORKDIR=${2}
NUMBASES=${3}
SPORK_INSTALL_DIR=${4}

cd ${4}

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
#STEM=Fetal_Adrenal_360_CTTGTA_L006

MODE="complete"
DENOVOCIRC="1"
ORIGDIR=${1}orig/
NUMFLANKING="150"
OUTPUTDIR=${2}fasta/${STEM}/

mkdir -p ${OUTPUTDIR}/


if [ -e ${2}fasta/${STEM}/${STEM}_SPORK_Junctions.fa ]
then
echo "SPORK fasta has already been generated"
else
rm ${2}fasta/${STEM}/*
ml load python/2.7.5
python denovo_pipeline_GH.py ${ORIGDIR} ${STEM} ${MODE} ${NUMFLANKING} ${NUMBASES} ${DENOVOCIRC} ${OUTPUTDIR}
python SPORK_reformat_header.py -f ${OUTPUTDIR} -s ${STEM}
echo "successfully generated SPORK fastas in ${OUTPUTDIR}" >> ${2}MasterError.txt
fi
