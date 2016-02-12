#!/bin/sh

#  BowtieIndexFJIndels.sh
#  
#
#  Created by Gillian Hsieh on 2/9/16.
#

FastaDir=${1}
IndelNumber=${2}
BowtieIndexDir=${3}
FJDir=${4}

STEMFILE=${4}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${3}${STEM}

bowtie2-build ${1}${STEM}/${STEM}_FJ_Indels_${2}.fa ${3}${STEM}/${STEM}_Indels_${2}

# 1= fasta file
# 2 = index name