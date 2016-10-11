#!/bin/sh

#  parseFJforGLM.sh
#  
#
#  Created by Gillian Hsieh on 2/20/16.
#

FJ_File=${1}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

FJ_AlignmentsFiles=${1}FarJunctionAlignments/${STEM}/*.sam
FJ_outdir=${1}FJ_out_for_GLM/
mkdir -p ${FJ_outdir}
FJout=${FJ_outdir}${STEM}_FJ_score_output.txt

rm ${FJout}

for file in ${FJ_AlignmentsFiles}
do
sed '1,2d' ${file} | grep "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $14 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' >> ${FJout}
sed '1,2d' ${file} | grep -v "XS:i:" | awk 'x=length($10) {print $1 "\t" $4 "\t" $5 "\t" $12 "\t" $13 "\t" x "\t" $3}' | sed 's/AS:i://' | sed 's/XN:i://' >> ${FJout}

done



