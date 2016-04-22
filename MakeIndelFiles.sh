#!/bin/sh
#  MakeIndelFiles.sh
#  
#
#  Created by Gillian Hsieh on 1/10/16.
#

#MakeIndelFiles.sh is a shell script that calls the python script AddIndelsToFasta.py.  It takes the FarJunctions fasta files as inputs (FJDir/fasta/<STEM>_FarJunctions.fa) and outputs five files called FJDir/FarJuncIndels/<STEM>/<STEM>_FJ_Indels_1,2,3,4,5.fa where the numbers 1-5 indicate the N's inserted on each side of the breakpoint or deletions on each side of the breakpoint.  For example, the FJ_indels_3 file is the same as the FarJunctions.fa file except every sequence has 3 N's on each side of the breakpoint (total of 6 N's inserted at the breakpoint), or 3 bases are deleted from each exon on each side of the breakpoint (total of 6 deletions at the breakpoint).

FJDir=${1}
NumIndels=${2}
INSTALLDIR=${3}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

InputFile=${1}fasta/${STEM}/${STEM}_FarJunctions.fa
if [ ! -e ${InputFile} ]
then
InputFile=${1}fasta/${STEM}/${STEM}_SPORK_Junctions.fa
fi

mkdir -p ${1}FarJuncIndels/${STEM}/

ml load python/2.7.5
python ${INSTALLDIR}AddIndelsToFasta.py -i ${InputFile} -o ${1}FarJuncIndels/${STEM}/ -s ${STEM} -n ${2}

echo "MakeIndelFiles.sh done -- check for ${1}FarJuncIndels/${STEM}/${STEM}_FJ_indels_*.fa where * = the numbers 1-5" >> ${1}MasterError.txt

