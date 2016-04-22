#!/bin/sh

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 4/20/16.
#


CircPipeDir=${1}
IndelNumber=${2}
FarJuncDir=${3}
BOWTIEPARAM=${4}
REG_INDEL_INDICES=${5}

STEM=SRR3192409


AlignedIndels=${1}orig/RegIndelAlignments/${STEM}/

if [ -e ${AlignedIndels}unaligned_${STEM}_1_indel${2}.sam ] && [ -e ${AlignedIndels}unaligned_${STEM}_2_indel${2}.sam ]
then
echo "Reg Alignments exist, skipping step"
else
echo "running some crap"
fi