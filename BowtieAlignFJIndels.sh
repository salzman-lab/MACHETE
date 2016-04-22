#!/bin/sh

#  BowtieAlignFJIndels.sh
#  
#
#  Created by Gillian Hsieh on 2/10/16.
#

# This section calls the shell BowtieAlignFJIndels.sh to align the fq files FJdir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_1/2.fq to the Bowtie2 indices of the far junction indels created in the previous step.  Aligned indels are stored in FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam where N is the number of indels in the Bowtie index where the reads aligned.  The bowtie parameters include a max of ~4 mismatches / 100 basepairs, max #N is the read length, and prohibits gapped alignments or gaps in the read.


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

echo "BowtieAlignFJIndels.sh complete- check for ${1}FarJunctionSecondary/AlignedIndels/${STEM}/still_unaligned_${STEM}_indels${3}.sam" >> ${1}MasterError.txt