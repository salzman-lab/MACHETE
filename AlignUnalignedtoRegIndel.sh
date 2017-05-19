#!/bin/sh

#  AlignUnalignedtoRegIndel.sh
#  
#
#  Created by Gillian Hsieh on 2/16/16.
#
#To train the GLM, indel alignments are also created for the linear junctions.  The reference index of indels to the linear junctions is static and has already been created and is referenced above as "REG_INDEL_INDICES" on line 44.  The script AlignRegIndels calls bowtie to align reads that were unaligned the the KNIFE indices (in KNIFEdir/orig/unaligned/*.fq) to the REG_INDEL_INDICES, with the parameters of 1) approx 4 mismatches / 100 bases, maximum number N's = readlength, and no gapped alignments or read gaps.


CircPipeDir=${1}
IndelNumber=${2}
FarJuncDir=${3}
BOWTIEPARAM=${4}
REG_INDEL_INDICES=${5}


STEMFILE=${3}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

glmReportsDir=${1}circReads/glmReports/
unalignedDir=${1}orig/unaligned/
AlignedIndels=${1}orig/RegIndelAlignments/${STEM}/
mkdir -p ${AlignedIndels}

IndelIndex="${5}hg19_junctions_reg_indels_${2}"

if [ "$(ls -A ${AlignedIndels}unaligned_${STEM}_*_indel${2}.sam)" ]
then
    echo "Reg Alignments to ${STEM} exist, skipping step" >> ${3}MasterError.txt
else
    for file in ${unalignedDir}*${STEM}*.fq
    do
    FILENAME=$(basename $file .fq)
    bowtie2 ${BOWTIEPARAM} -x ${IndelIndex} -U ${file} -S ${AlignedIndels}${FILENAME}_indel${2}.sam
    done
    echo "AlignUnalignedtoRegIndel.sh complete for ${STEM} - check ${AlignedIndels}unaligned_${STEM}_*_indel${2}.sam" >> ${3}MasterError.txt
fi