#!/bin/sh

#  BowtieAlignFJIndels.sh
#  
#
#  Created by Gillian Hsieh on 2/10/16.
#

# aligns still unaligned files in far junc secondary to the indel indices.

FJDir=${1}
BOWTIEPARAMETERS=${2}
IndelNum=${3}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

InputFQDir=${1}FarJuncSecondary/${STEM}/
OutputDir=${1}FarJuncSecondary/AlignedIndels/${STEM}/
mkdir -p ${OutputDir}
Index=${1}BowtieIndels/${STEM}/${STEM}_Indels_${3}

# bowtie2 <options> -x <index> -U <fq files> -S <sam output>

for file in ${InputFQDir}*
do
FILENAME=$(basename "$file" .fq)
bowtie2 ${2} -x ${Index} -U ${file} -S ${OutputDir}${FILENAME}_indels${3}.sam
done

