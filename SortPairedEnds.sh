#!/bin/sh

#  SortPairedEnds.sh
#  
#
#  Created by Gillian Hsieh on 1/7/16.
#


FJFile=${1}
STEMFILE=${1}StemList.txt

STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

DIR_TO_SORT=${1}DistantPEFiles/${STEM}/

for file in ${DIR_TO_SORT}* ; do
FILENAME=$(basename $file)
SORTEDNAME=sorted_${FILENAME}

if ! [[ "$FILENAME" = *sorted* ]]
then
#echo "sort -k1,1.1 -k1,1.2 -k2,2 ${file} > ${DIR_TO_SORT}${SORTEDNAME}"
sort -k1,1.1 -k1,1.2 -k2,2 ${file} > ${DIR_TO_SORT}${SORTEDNAME}
fi

done;


#for c in `seq 1 24`; do
#A=${c}
#if [ ${c} -eq 23 ]; then A=X; fi
#if [ ${c} -eq 24 ]; then A=Y; fi
#
#
#sort -k1,1.1 -k1,1.2 -k2,2 ${1} > ${2}
#
#done;