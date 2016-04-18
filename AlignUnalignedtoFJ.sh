#!/bin/sh

#  AlignUnalignedtoFJ.sh
#  
#
#  Created by Gillian Hsieh on 2/5/16.
#

#  AlignUnalignedtoFJ takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.

FJDir=${1}  ## MACHETE output dir
OrigDir=${2} ## KNIFE alignment directories

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
STILL_UNALIGNED_DIR=${1}FarJuncSecondary/${STEM}/

BOWTIE_INDEX=${1}BowtieIndex/${STEM}/${STEM}_FJ_Index
OUTPUT_SAM_DIR=${1}FarJunctionAlignments/${STEM}/

mkdir -p ${OUTPUT_SAM_DIR}
mkdir -p ${STILL_UNALIGNED_DIR}


for file in ${2}unaligned/*${STEM}*.fq
do
FILENAME=$(basename "$file" .fq)
BOWTIEPARAM="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50 --n-ceil L,0,1 -p 4 --un ${STILL_UNALIGNED_DIR}still_${FILENAME}.fq"
COMMAND="${BOWTIEPARAM} -x ${BOWTIE_INDEX} -U ${file} -S ${OUTPUT_SAM_DIR}${FILENAME}.sam"

bowtie2 ${COMMAND}


done;