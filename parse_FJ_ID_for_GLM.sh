#!/bin/sh

#  parse_FJ_ID_for_GLM.sh
#  
#
#  Created by Gillian Hsieh on 2/22/16.
#

################
#Current categories
#FJgood -- genome, reg, FJ
#FJbad -- genome anomaly, reg anomaly, junc, junc anomaly, FJ anomaly
#################

FJDir=${1} # FJ Dir
circreads=${2} # - circreads dir

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`


#Linda's outputfile
FJ_Outputfile=${2}ids/${STEM}_1__output_FJ.txt

for file in ${1}reports/IDs*${STEM}*.txt
do
sed '1d' ${file} | grep -v "Unmapped\|unaligned" > ${FJ_Outputfile}
done
