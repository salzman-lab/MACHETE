#!/bin/sh

#  linkfastafiles.sh
#  
#
#  Created by Gillian Hsieh on 1/31/16.
#

##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.


FJDir=${1} #MACHETE output directory


STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

ChrFastaDir=${1}fasta/${STEM}/
BigFastaDir=${1}fasta/
BigFastaFile=${BigFastaDir}${STEM}_FarJunctions.fa
BowtieIndex=${1}BowtieIndex/${STEM}/${STEM}_FJ_Index
mkdir -p ${1}BowtieIndex/${STEM}/

## The script linkfastafiles.sh uses linux to concatenate the <FJDir>/fasta/<STEM>/<STEM>_chr1,2,3,...,X,Y_FarJunctions.fa into a single large fasta <FJDir>/fasta/<STEM>_FarJunctions.fa.
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

## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>_FJ_Index
bowtie2-build --threads 4 ${BigFastaFile} ${BowtieIndex}

echo "linkfastafiles and Bowtie Index built for Far Junctions for sample ${STEM}" >> ${1}MasterError.txt
