#!/bin/sh

#  FindAlignmentArtifact.sh
#  
#
#  Created by Gillian Hsieh on 1/6/16.
#

FJFile=${1}
NumBParoundJunc=${2}
NumIndels=${3}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${1}FarJuncSecondary/AlignedIndels/RemoveNonOverlap/${STEM}/

module load python/2.7.5

python /scratch/PI/horence/gillian/MACHETE/FindAlignmentArtifact.py -s ${STEM} -p ${1} -n ${2} -x ${3}