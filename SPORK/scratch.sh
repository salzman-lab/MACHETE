#!/bin/sh

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 4/20/16.
#


CIRCPIPE_DIR=${1} #Main circular RNA pipeline output directory -- contains subfolders "orig", "circReads", "logs", and "sample stats"
OUTPUT_DIR=${2} # Output directory for MACHETE - does not need to exist already

INSTALLDIR="/scratch/PI/horence/gillian/MACHETE/"
RESOURCE_FLAG="-p owners"
CIRCREF="/share/PI/horence/circularRNApipeline_Cluster/index/"
genomeIndex=${CIRCREF}hg19_genome

START=1
for (( c=$START; c<=16; c++ ))
do
STEM=`awk 'FNR == '${c}' {print $1}' ${StemFile}`

FarJuncFasta=${2}fasta/${STEM}*FarJunctions.fa
BadFJDir=${2}BadFJ/${STEM}/
BadFJver2Dir=${2}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}
mkdir -p ${BadFJver2Dir}
r1file=${BadFJver2Dir}${STEM}_FarJunctions_R1.fa
r2file=${BadFJver2Dir}${STEM}_FarJunctions_R2.fa

genomeBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${genomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam"

BadFJv2j1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 --nodes=4 ${RESOURCE_FLAG} --time=24:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${genomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to genome: ${BadFJv2j1_id}"
