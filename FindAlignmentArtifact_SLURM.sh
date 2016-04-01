#!/bin/sh

#  FindAlignmentArtifact.sh
#  
#
#  Created by Gillian Hsieh on 1/6/16.
#

FJFile=${1}
NumBParoundJunc=${2}
NumIndels=${3}
INSTALLDIR=${4}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${1}FarJuncSecondary/AlignedIndels/RemoveNonOverlap/${STEM}/

module load python/2.7.5

python ${INSTALLDIR}MakeIndelsHisto.py -s ${STEM} -f ${1} -w ${2} -x ${3}