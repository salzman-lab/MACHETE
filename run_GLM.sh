#!/bin/sh

#  run_GLM.sh
#  
#
#  Created by Gillian Hsieh on 3/3/16.
#
CircularPipelineDir=${1} #directory that contains circReads, orig, logs, etc
FJDir=${2}

#testing mode

#STEM="H3-AD015"

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

INPUTDIR=${1}circReads/ids/
OUTPUTDIR=${2}reports/glmReports/
mkdir -p ${OUTPUTDIR}

for file in ${INPUTDIR}*${STEM}*output.txt
do
class_input=${file}
done

for file in ${INPUTDIR}*${STEM}*output_FJ.txt
do
FJ_input=${file}
done

ml load R/3.0.2
Rscript GLM_script.r ${FJ_input} ${class_input} ${STEM} ${OUTPUTDIR}