#!/bin/sh

#  run_GLM.sh
#  
#
#  Created by Gillian Hsieh on 3/3/16.
#
InputFolder=${1} #directory that contains all the class input files
OutputFolder=${2}
UniqueIDList=${3}
Version=${4} #UseIndel vs NoIndel vs old
FJFolder=${5}
RegGLMDir=${6}

mkdir -p ${2}


if [ $# -ge 7 ]
then
MODE=${7}
else
MODE=sam
fi

if [[ "$Version" = UseIndel ]]
then
RSCRIPT=GLM_script_Apr11_UseIndel.r
fi
if [[ "$Version" = NoIndel ]]
then
RSCRIPT=GLM_script_Mar21_NoIndel.r
fi
if [[ "$Version" = old ]]
then
RSCRIPT=GLM_script.r
fi


if [[ "$MODE" = *owners* ]]
then
RESOURCE_FLAG="-p owners"
fi

mkdir -p ${2}err_and_out/

INSTALLDIR="/scratch/PI/horence/gillian/MACHETE/"

NUM_FILES=$((`more ${UniqueIDList} | wc -l`))
echo "${NUM_FILES} unique IDs"

#testing mode

#STEM="H3-AD015"

j_id=`sbatch -J GLM.r ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=5:0:0 -o ${2}err_and_out/out_independentGLM_r.txt -e ${2}err_and_out/err_independentGLM_r.txt ${INSTALLDIR}run_GLM_independent_step2.sh ${1} ${2} ${UniqueIDList} ${RSCRIPT} ${5} ${6} | awk '{print $4}'`

echo "Run GLM: ${j_id}"
