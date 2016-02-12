#!/bin/sh

#  DistantPE_Counter_genome_ENCODE.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

FJFile=${1}
STEMFILE=${1}StemList.txt

STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir ${1}DistantPEFiles/${STEM}

#module load python/2.7.9
ml python/2.7.5

echo "python /scratch/PI/horence/gillian/MACHETE/DistantPE_Counter_genome_ENCODE.py -d ${1} -s ${STEM}"
python /scratch/PI/horence/gillian/MACHETE/DistantPE_Counter_genome_ENCODE.py -d ${1} -s ${STEM}


#python /srv/gsfs0/projects/salzman/gillian/createFarJunctionsIndex/DistantPE_Counter_genome_ENCODE.py -d ${1} -s ${STEM}

