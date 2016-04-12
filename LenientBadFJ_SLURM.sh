#!/bin/sh

#  GoodvsBadFJ.sh
#  
#
#  Created by Gillian Hsieh on 1/18/16.
#

# this shell generates far junction alignments to genome,transcriptome,reg, and junc.
# after this shell, should run python script AddIndelsandBadFJtoNaiveReports

FarJuncDir=${1}
REFGENOME=${2}
INSTALLDIR=${3}
CIRCREF=${4}
RESOURCE_FLAG=${5}

##bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align
STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`


BadFJDir=${1}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}

FastaDir=${1}fasta/${STEM}*FarJunctions.fa
for file in ${FastaDir}
do
FarJuncFasta=${file}
echo ${FarJuncFasta}
done



#indices for alignment
#if [[ "$REFGENOME" = *HG38* ]]
#then
#### INDICES NOT READY YET
#fi

if [[ "$REFGENOME" = *HG19* ]]
then
genomeIndex=${4}hg19_genome
transcriptomeIndex=${4}hg19_transcriptome
regIndex=${4}hg19_junctions_reg
juncIndex=${4}hg19_junctions_scrambled
fi


# Align FarJunc fasta file to the above indices:


BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-.5 --trim5 120 --trim5 120 --n-ceil 300 --np 0 --rdg 5,.00001"

BadFJj1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o out.txt -e err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoGenome.sam | awk '{print $4}'`

BadFJj2_id=`sbatch -J ${STEM}FJ_to_transcriptome --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o out.txt -e err.txt /${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtotranscriptome.sam | awk '{print $4}'`

BadFJj3_id=`sbatch -J ${STEM}FJ_to_reg --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o out.txt -e err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoReg.sam | awk '{print $4}'`

BadFJj4_id=`sbatch -J ${STEM}FJ_to_junc --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o out.txt -e err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoJunc.sam | awk '{print $4}'`