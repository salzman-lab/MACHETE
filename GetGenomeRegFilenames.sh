#!/bin/sh

#  GetGenomeRegFilenames.sh
#  
#
#  Created by Gillian Hsieh on 1/28/16.
#

origDir=${1}
FJDir=${2}

PEFiles=${2}PE_ID_files.txt


#module load python/2.7.9
ml load python/2.7.5

#python /srv/gsfs0/projects/salzman/gillian/createFarJunctionsIndex/PEfinder_genomeAndReg_ENCODE.py -i ${1} -o ${2} -w 10000 -n ${3}
python /scratch/PI/horence/gillian/MACHETE/writePEMatchIDFiles.py -o ${1} -f ${2}