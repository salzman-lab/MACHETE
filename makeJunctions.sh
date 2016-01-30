#!/bin/sh

#  makeJunctions.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

var1=${1}
var2=${2}

#module load python/2.7
ml load python/2.7.5

#python /srv/gsfs0/projects/salzman/gillian/createFarJunctionsIndex/makeJunctions.py -p ${1} -o ${2}

python /scratch/PI/horence/gillian/createFarJunctionsIndex/makeJunctions.py -p ${1} -o ${2}