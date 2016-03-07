#!/bin/sh

#  PEfinder_genomeAndReg_ENCODE.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

OrigDir=${1}
FJDir=${2}
BP_distance=${3}
INSTALLDIR=${4}

STEMFILE=${2}StemList.txt

#module load python/2.7.9
ml load python/2.7.5

STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

python ${INSTALLDIR}PEfinder_genomeAndReg_ENCODE.py -o ${1} -s ${STEM} -f ${2} -w 10000 -n ${3}

