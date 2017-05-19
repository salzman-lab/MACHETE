#!/bin/sh

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 8/18/15.
#

SORTINGDIR=${1}
FJDir=${2}

StemFile=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${StemFile}`

for file in ${1}*${STEM}*
    do
    FILENAME=$(basename $file)
    SORTEDNAME=sorted_${FILENAME}
        if ! [[ "$FILENAME" = *sorted* ]]
        then
        head -n 2 ${file} > ${SORTINGDIR}${SORTEDNAME}
        tail -n +3 ${file} | sort -k 1 >> ${SORTINGDIR}${SORTEDNAME}
        fi
    done;


echo "completed Sorting step for all files in ${1}" >> ${2}MasterError.txt