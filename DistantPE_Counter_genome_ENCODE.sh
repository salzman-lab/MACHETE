#!/bin/sh

#  DistantPE_Counter_genome_ENCODE.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

Var1=${1}

#module load python/2.7.9
ml python/2.7.5


#python /srv/gsfs0/projects/salzman/gillian/createFarJunctionsIndex/DistantPE_Counter_genome_ENCODE.py -d ${1}

python /scratch/PI/horence/gillian/createFarJunctionsIndex/DistantPE_Counter_genome_ENCODE.py -d ${1}