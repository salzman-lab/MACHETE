#!/bin/sh

#  DenovoJunctions.sh
#  
#
#  Created by Gillian Hsieh on 11/3/15.
#


TRIMMEDFQDIR=${1}
OUTPUTDIR=${2} # with bowtie index for alignment
READLENGTH=${3}

TRIMLEN=$((READLENGTH * 2 / 3))  # plan to trim 2/3 of the read from each side and align only 1/3
echo $TRIMLEN

# BOWTIE COMMANDS
# -5/--trim5 <int>
# Trim <int> bases from 5' (left) end of each read before alignment (default: 0).
# -3/--trim3 <int>
# Trim <int> bases from 3' (right) end of each read before alignment (default: 0).

DENOVODIR=${2}/FarJuncDenovo/
BOWTIEINDEX=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_genome


mkdir -p ${DENOVODIR}


#align 3' ends to the index
for file in ${1}/*.fq
do
echo ${file}
FILENAME=$(basename "$file" .fq)

bowtie2 --trim5 ${TRIMLEN} --rdg 50,50 --rfg 50,50 --score-min L,0,-0.24 --no-unal --no-sq -x ${BOWTIEINDEX} -U ${file} -S ${DENOVODIR}${FILENAME}_right.sam
done;

#align 5' ends to the index
for file in ${1}/*.fq
do
echo ${file}
FILENAME=$(basename "$file" .fq)

bowtie2 --trim3 ${TRIMLEN} --rdg 50,50 --rfg 50,50 --score-min L,0,-0.24 --no-unal --no-sq -x ${BOWTIEINDEX} -U ${file} -S ${DENOVODIR}${FILENAME}_left.sam
done;

sh AlphebetizeENCODEreads.sh ${DENOVODIR}