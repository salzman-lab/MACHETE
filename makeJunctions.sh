#!/bin/sh
#  makeJunctions.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#


PICKLEDIR=${1}
FJFile=${2}
Chromosome=${3}
INSTALLDIR=${4}

STEMFILE=${2}StemList.txt
FASTADIR=${2}fasta/

STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

OUTPUTDIR=${FASTADIR}${STEM}/
mkdir -p ${OUTPUTDIR}
INPUTDIR=${2}DistantPEFiles/${STEM}/

ml load python/2.7.5
for file in ${INPUTDIR}/sorted_chr${3}_*; do
python ${INSTALLDIR}makeJunctions.py -p ${1} -f ${file} -o ${OUTPUTDIR} -s ${STEM}
done;

rm ${OUTPUTDIR}${STEM}_chr${3}_*_duplicates.fa