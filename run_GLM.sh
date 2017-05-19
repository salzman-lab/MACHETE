#!/bin/sh

#  run_GLM.sh
#  
#
#  Created by Gillian Hsieh on 3/3/16.
#
CircularPipelineDir=${1} #directory that contains circReads, orig, logs, etc
FJDir=${2}
INSTALLDIR=${3}

if [ $# -ge 4 ]
then
OUTPUTDIR=${4}
else
OUTPUTDIR=${2}reports/glmReports/
fi
echo $OUTPUTDIR

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
#STEM=SRR1594021

REG_INPUTDIR=${1}circReads/ids/
FJ_INPUTDIR=${2}GLM_classInput/

mkdir -p ${OUTPUTDIR}

for file in ${REG_INPUTDIR}*${STEM}*output.txt
do
reg_class_input=${file}
done

for file in ${FJ_INPUTDIR}*${STEM}*output_FJ.txt
do
FJ_input=${file}
done

for file in ${REG_INPUTDIR}*${STEM}*output_RegIndel.txt
do
RegIndel_input=${file}
done

for file in ${FJ_INPUTDIR}*${STEM}*output_FJIndels.txt
do
FJIndel_input=${file}
done

ml load R/3.0.2
Rscript ${INSTALLDIR}GLM_script_UseIndel.r ${FJ_input} ${reg_class_input} ${STEM} ${OUTPUTDIR} ${RegIndel_input} ${FJIndel_input}
#GLM script most recently updated May 6 -- copied from /scratch/PI/horence/julia/MACHETE/GLM_script_May6_UseIndel.r

echo "run_GLM.sh complete for ${STEM} -- check ${OUTPUTDIR}/*${STEM}* for Far Junction GLM outputs" >> ${2}MasterError.txt
