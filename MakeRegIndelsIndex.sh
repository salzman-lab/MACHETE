#!/bin/sh

#  MakeRegIndelsIndex.sh
#  
#
#  Created by Gillian Hsieh on 2/9/16.
#
##  This shell takes the linear junctions fasta file from KNIFE, expands it by removing or adding up to n indels on each side of the junction interface, and uses Bowtie to make the linear junctions indels fastas into Bowtie indices.

RegIndelsFasta=${1}
OutputDir=${2}
IndelNumber=${3}
genome=${4}
ResourceFlag=${5} # For example, "-p owners" or "-p horence" or can be blank

mkdir -p ${OutputDir}
mkdir -p ${OutputDir}/IndelIndices/
ml load python/2.7.5
python MakeRegIndelsFasta.py -i ${RegIndelsFasta} -o ${OutputDir} -g ${genome} -n ${IndelNumber}


START=1
for (( c=$START; c<=${IndelNumber}; c++ ))
do

sbatch -J IndexRegIndels ${ResourceFlag} --mem=55000 --nodes=4 --time=18:0:0 BowtieIndexer.batch.sh ${OutputDir}/${genome}_reg_indels_${c}.fa ${OutputDir}/IndelIndices/${genome}_reg_indels_${c}

done

