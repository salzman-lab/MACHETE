#!/bin/sh

#  GoodvsBadFJ.sh
#  
#
#  Created by Gillian Hsieh on 1/18/16.
#

# this shell generates far junction alignments to genome,transcriptome,reg, and junc.
# after this shell, should run python script AddIndelsandBadFJtoNaiveReports

FarJuncDir=${1}


FarJuncFasta=${1}FarJunctions.fa

#bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align
OutputDir=${1}BadFJ/
mkdir -p ${OutputDir}


#indices for alignment
genomeIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_genome
transcriptomeIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_transcriptome
regIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_junctions_reg
juncIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_junctions_scrambled

# Align FarJunc fasta file to the above indices:


BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --np 0 --rdg 50,50 --rfg 50,50"

j1_id=`sbatch -J FJ.fa_to_genome --mem=55000 -p owners --time=24:0:0 -o out.txt -e err.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${FarJuncFasta} ${OutputDir}BadFJtoGenome.sam | awk '{print $4}'`

j2_id=`sbatch -J FJ.fa_to_transcriptome --mem=55000 -p owners --time=24:0:0 -o out.txt -e err.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${FarJuncFasta} ${OutputDir}BadFJtotranscriptome.sam | awk '{print $4}'`

j3_id=`sbatch -J FJ.fa_to_reg --mem=55000 -p owners --time=24:0:0 -o out.txt -e err.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${FarJuncFasta} ${OutputDir}BadFJtoReg.sam | awk '{print $4}'`

j4_id=`sbatch -J FJ.fa_to_genome --mem=55000 -p owners --time=24:0:0 -o out.txt -e err.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${FarJuncFasta} ${OutputDir}BadFJtoJunc.sam | awk '{print $4}'`

depend_str=--depend=afterok:${j1_id}:${j2_id}:${j3_id}:${j4_id}


