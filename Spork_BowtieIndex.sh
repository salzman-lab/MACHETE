#!/bin/sh

#  bowtieindexer.batch.sh
#
#
#  Created by Gillian Hsieh on 9/22/15.
#

#bowtie2-build [options]* <reference_in> <bt2_base>

FJDir=${1}
IndexName

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
#STEM=${2}

BowtieIndex=${1}BowtieIndex/${STEM}/${STEM}_FJ_Index
BigFastaFile=${1}fasta/*${STEM}*FarJunctions.fa

mkdir -p ${1}BowtieIndex/${STEM}/


module load bowtie/2.2.4

bowtie2-build ${BigFastaFile} ${BowtieIndex}