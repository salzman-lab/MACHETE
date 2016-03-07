#!/bin/sh

#  AddIndelstoGLM.sh
#  
#
#  Created by Gillian Hsieh on 2/17/16.
#

CirclePipeDir=${1}
FJDir=${2}
window=${3}

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

ml load python/2.7.5
python /scratch/PI/horence/gillian/MACHETE/AddIndelstoGLM.py -c ${1} -s ${STEM} -w ${window}