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

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
#STEM=SRR1027188

cd ${4}
pwd

MODE="complete"
DENOVOCIRC="1"
ORIGDIR=${1}orig/
NUMFLANKING="150"
OUTPUTDIR=${2}fasta/${STEM}/

mkdir -p ${OUTPUTDIR}/
ml load python/2.7.5

python denovo_pipeline_GH.py ${ORIGDIR} ${STEM} ${MODE} ${NUMFLANKING} ${NUMBASES} ${DENOVOCIRC} ${OUTPUTDIR}

python SPORK_reformat_header.py -f ${OUTPUTDIR} -s ${STEM}

echo "successfully generated SPORK fastas in ${OUTPUTDIR}" >> ${2}MasterError.txt