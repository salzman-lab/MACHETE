#!/bin/sh

#  AppendNaiveRept.sh
#  
#
#  Created by Gillian Hsieh on 2/1/16.
#

FarJuncDir=${1}
GLMReportDir=${2}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${1}reports/AppendedReports/

ml load python/2.7.5

python /scratch/PI/horence/gillian/MACHETE/AddIndels_BadFJ_LinearFreq_ToNaiveRept.py -f ${1} -g ${2} -s ${STEM}