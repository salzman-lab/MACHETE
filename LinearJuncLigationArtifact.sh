#!/bin/sh

#  LinearJuncLigationArtifact.sh
#  
#
#  Created by Gillian Hsieh on 1/13/16.
#

CircReadsAndOrigParentDir=${1}
FJDir=${2}
NumIndels=${3}
NumBPOverlapAtJunc=${4}
HiPPcutoff=${5}
LoPPcutoff=${6}

origDir=${1}orig/
circReadsDir=${1}circReads/

mkdir -p {2}LinearIndelsCDF/
mkdir -p ${origDir}unaligned/AlignToRegIndels/
mkdir -p ${origDir}unaligned/AlignToRegIndels/RemoveNonOverlap/



# align orig/unaligned reads to indels indices in orig/unaligned/AlignToRegIndels
echo "Aligning reads to indels indices"

START=1
depend_str="--depend=afterok"
for (( c=$START; c<=${3}; c++))
do
    BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
    BOWTIEINDEX=/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/hg19_junctions_reg_indels_${c}
    BOWTIEOUTPUTDIR=${origDir}unaligned/AlignToRegIndels/

    for file in ${origDir}unaligned/*.fq
    do
        FILENAME=$(basename "$file" .fq)
        OUTPUTFILE=${FILENAME}_indels${c}.sam
        j_id=`sbatch -J AlignIndels --mem=55000 --time=24:0:0 -o out.txt -e err.txt /scratch/PI/horence/gillian/MACHETE/BowtieAligner.batch.sh "${BOWTIEPARAMETERS}" ${BOWTIEINDEX} ${file} ${BOWTIEOUTPUTDIR}${OUTPUTFILE} | awk '{print $4}'`
        depend_str=${depend_str}:${j_id}
    done

done

echo ${depend_str}

# call python script to remove all non-overlapping reads, and then compile the indel cdf files from the glmreports and reports directories.

python FindAlignmentArtifact_LinearJunc.py -s ${BOWTIEOUTPUTDIR} -c ${circReadsDir} -o ${2}LinearIndelsCDF/ -H ${5} -L ${6} -x ${4}


