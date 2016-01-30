#!/bin/sh

#  FarJuncNaiveReport.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

var1=${1}
var2=${2}
var3=${3}

# module load python/2.7.9
ml load python/2.7.5

# python /srv/gsfs0/projects/salzman/gillian/createFarJunctionsIndex/FarJuncNaiveReport.py -o ${1} -i ${2} -w ${3}

python /scratch/PI/horence/gillian/createFarJunctionsIndex/FarJuncNaiveReport.py -o ${1} -i ${2} -w ${3}