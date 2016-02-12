#!/bin/sh

#  AlignUnalignedtoFJ.sh
#  
#
#  Created by Gillian Hsieh on 2/5/16.
#

FJDir=${1}
OrigDir=${2}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
STILL_UNALIGNED_DIR=${1}FarJuncSecondary/${STEM}/

BOWTIE_INDEX=${1}BowtieIndex/${STEM}/${STEM}_FJ_Index
OUTPUT_SAM_DIR=${1}FarJunctionAlignments/${STEM}/

mkdir -p ${OUTPUT_SAM_DIR}
mkdir -p ${STILL_UNALIGNED_DIR}


for file in ${2}unaligned/*${STEM}*.fq
do
FILENAME=$(basename "$file" .fq)
BOWTIEPARAM="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50 --un ${STILL_UNALIGNED_DIR}still_${FILENAME}.fq"
COMMAND="${BOWTIEPARAM} -x ${BOWTIE_INDEX} -U ${file} -S ${OUTPUT_SAM_DIR}${FILENAME}.sam"

bowtie2 ${COMMAND}


done;