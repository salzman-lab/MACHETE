#!/bin/sh

#  FarJuncNaiveReport.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

FJDir=${1}
OrigDir=${2}
Window=${3}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

# module load python/2.7.9
ml load python/2.7.5

# python /srv/gsfs0/projects/salzman/gillian/createFarJunctionsIndex/FarJuncNaiveReport.py -o ${1} -i ${2} -w ${3}

python /scratch/PI/horence/gillian/MACHETE/FarJuncNaiveReport.py -s ${STEM} -f ${1} -i ${2} -w ${3}
#
#Variables
#-s STEM
#-f output dir (FJ Dir)
#-i orig Dir
#-w window - num bases on each side of junction