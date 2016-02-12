#!/bin/sh

#  linkfastafiles.sh
#  
#
#  Created by Gillian Hsieh on 1/31/16.
#

FJDir=${1}


STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
#STEM=${2}


ChrFastaDir=${1}fasta/${STEM}/
BigFastaDir=${1}fasta/
BigFastaFile=${BigFastaDir}${STEM}_FarJunctions.fa
BowtieIndex=${1}BowtieIndex/${STEM}/${STEM}_FJ_Index
mkdir -p ${1}BowtieIndex/${STEM}/

for c in `seq 1 24`; do

if [ ${c} -eq 1 ]; then cat ${ChrFastaDir}${STEM}_chr${c}FarJunctions.fa > ${BigFastaFile}
else
    A=${c}
    if [ ${c} -eq 23 ]; then A=X; fi
    if [ ${c} -eq 24 ]; then A=Y; fi
    echo $A
    cat ${ChrFastaDir}${STEM}_chr${A}FarJunctions.fa >> ${BigFastaFile}

fi
done;

bowtie2-build ${BigFastaFile} ${BowtieIndex}


