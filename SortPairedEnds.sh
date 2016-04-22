#!/bin/sh

#  SortPairedEnds.sh
#  
#
#  Created by Gillian Hsieh on 1/7/16.
#

## This is a simple shell script SortPairedEnds.sh to sort the chrA_Distant_PE_frequency.txt files into alphabetical order.  It takes FJDir/DistantPEFiles/chrA_Distant_PE_frequency.txt and outputs to same directory, sorted_chrA_Distant_PE_frequency.txt using the linux "sort" command.
## The reason for sorting is to increase the speed of the next step.

FJDir=${1} ## MACHETE output directory
STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

DIR_TO_SORT=${1}DistantPEFiles/${STEM}/


## loops through all the files in <FJ dir>/DistantPEFiles/<STEM/* and if it hasn't been previously sorted, then sorts it by numerical value so that discordant pairs with locations that are lower integers are first in the file, and higher integers are later in the file.
## In the next step, when unpickling, exons can be retrieved in order from the pickle. Accessing the pickle can be otherwise time consuming.

for file in ${DIR_TO_SORT}* ; do
FILENAME=$(basename $file)
SORTEDNAME=sorted_${FILENAME}

if ! [[ "$FILENAME" = *sorted* ]]
then
sort -k1,1.1 -k1,1.2 -k2,2 ${file} > ${DIR_TO_SORT}${SORTEDNAME}
fi

done;

echo "SortPairedEnds.sh complete for ${STEM}" >> ${1}MasterError.txt