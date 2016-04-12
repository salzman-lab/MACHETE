#!/bin/sh
#  MakeIndelFiles.sh
#  
#
#  Created by Gillian Hsieh on 1/10/16.
#


FJDir=${1}
NumIndels=${2}
INSTALLDIR=${3}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`


mkdir -p ${1}FarJuncIndels/${STEM}/

for file in ${1}fasta/${STEM}*_FarJunctions.fa
do
InputFile=${file}
done


echo "Making FarJunctions.fa into Indel.fa files"

##module load python/2.7.9
##python /srv/gsfs0/projects/salzman/gillian/MACHETE/AddIndelsToFasta.py -o ${1} -n ${2}

ml load python/2.7.5
python ${INSTALLDIR}AddIndelsToFasta.py -i ${InputFile} -o ${1}FarJuncIndels/${STEM}/ -s ${STEM} -n ${2}

# PYTHON INPUTS
#-i infile
# -o outDir
# -s stem
#-n = max # indels
#