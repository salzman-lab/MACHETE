#!/bin/sh

#  bowtieindexer.batch.sh
#
#
#  Created by Gillian Hsieh on 9/22/15.
#

#bowtie2-build [options]* <reference_in> <bt2_base>

SPORKdir=${1}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
mkdir -p ${1}BowtieIndex/${STEM}/

BowtieIndex=${1}BowtieIndex/${STEM}/${STEM}_FJ_Index
FastaFile=${1}fasta/${STEM}/${STEM}_SPORK_Junctions.fa

bowtie2-build ${FastaFile} ${BowtieIndex}

echo "completed bowtie index of ${FastaFile} to ${BowtieIndex} " >> ${1}MasterError.txt