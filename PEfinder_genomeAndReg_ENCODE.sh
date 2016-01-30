#!/bin/sh

#  PEfinder_genomeAndReg_ENCODE.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

OrigDir=${1}
FJDir=${2}
BP_distance=${3}

PEFiles=${2}PE_ID_files.txt


#module load python/2.7.9
ml load python/2.7.5

echo ${SLURM_ARRAY_TASK_ID}

G1=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${PEFiles}`
G2=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $2}' ${PEFiles}`
R1=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $3}' ${PEFiles}`
R2=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $4}' ${PEFiles}`

echo ${G1}
echo ${G2}
echo ${R1}
echo ${R2}

python /scratch/PI/horence/gillian/createFarJunctionsIndex/PEfinder_genomeAndReg_ENCODE.py -o ${1} -g1 ${G1} -g2 ${G2} -r1 ${R1} -r2 ${R2} -f ${2} -w 10000 -n ${3}



