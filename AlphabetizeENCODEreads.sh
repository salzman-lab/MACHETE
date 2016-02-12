#!/bin/sh

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 8/18/15.
#

SORTINGDIR=${1}
StemFile=${2}

STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${2}`

for file in ${1}/*${STEM}*
do
FILENAME=$(basename $file)
SORTEDNAME=sorted_${FILENAME}
    if ! [[ "$FILENAME" = *sorted* ]]
    then
    head -n 2 ${file} > ${SORTINGDIR}${SORTEDNAME}
    tail -n +3 ${file} | sort -k 1 >> ${SORTINGDIR}${SORTEDNAME}
    fi
done;