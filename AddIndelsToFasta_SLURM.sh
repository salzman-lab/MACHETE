#!/bin/sh

#  AddIndelsToFasta_SLURM.sh
#  
#
#  Created by Gillian Hsieh on 1/11/16.
#

ml load python/2.7.5

python /scratch/PI/horence/gillian/createFarJunctionsIndex/AddIndelsToFasta.py -o ${1} -n ${2}