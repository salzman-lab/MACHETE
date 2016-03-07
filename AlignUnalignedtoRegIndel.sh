#!/bin/sh

#  AlignUnalignedtoRegIndel.sh
#  
#
#  Created by Gillian Hsieh on 2/16/16.
#

CircPipeDir=${1}
IndelNumber=${2}
FarJuncDir=${3}
BOWTIEPARAM=${4}

STEMFILE=${3}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

glmReportsDir=${1}circReads/glmReports/
unalignedDir=${1}orig/unaligned/
AlignedIndels=${1}orig/RegIndelAlignments/${STEM}/
mkdir -p ${AlignedIndels}

IndelIndex="/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/hg19_junctions_reg_indels_${2}"

for file in ${unalignedDir}*${STEM}*.fq
do
FILENAME=$(basename $file .fq)
bowtie2 ${BOWTIEPARAM} -x ${IndelIndex} -U ${file} -S ${AlignedIndels}${FILENAME}_indel${2}.sam
done