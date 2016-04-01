#!/bin/sh

#  run_GLM.sh
#  
#
#  Created by Gillian Hsieh on 3/3/16.
#
INPUTDIR=${1} #directory that contains circReads, orig, logs, etc
OUTPUTDIR=${2}
STEMFILE=${3}
RSCRIPT=${4}

#echo ${4}

#testing mode
#RSCRIPT=GLM_script_Mar22_UseIndel.r
#STEM="SRR1027190"
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`



mkdir -p ${OUTPUTDIR}

for file in ${INPUTDIR}*${STEM}*output.txt
do
class_input=${file}
done

for file in ${INPUTDIR}*${STEM}*output_FJ.txt
do
FJ_input=${file}
done

for file in ${INPUTDIR}*${STEM}*output_RegIndel.txt
do
RegIndel_input=${file}
done

for file in ${INPUTDIR}*${STEM}*output_FJIndels.txt
do
FJIndel_input=${file}
done

ml load R/3.0.2

if [[ "$RSCRIPT" = *UseIndel* ]]
then
echo "Rscript /scratch/PI/horence/gillian/MACHETE/${RSCRIPT} ${FJ_input} ${class_input} ${STEM} ${OUTPUTDIR} ${RegIndel_input} ${FJIndel_input} "
Rscript /scratch/PI/horence/gillian/MACHETE/${RSCRIPT} ${FJ_input} ${class_input} ${STEM} ${OUTPUTDIR} ${RegIndel_input} ${FJIndel_input}
else
echo "Rscript /scratch/PI/horence/gillian/MACHETE/${RSCRIPT} ${FJ_input} ${class_input} ${STEM} ${OUTPUTDIR}"
Rscript /scratch/PI/horence/gillian/MACHETE/${RSCRIPT} ${FJ_input} ${class_input} ${STEM} ${OUTPUTDIR}
fi



##
#fusion_class_input=args[1]
#class_input=args[2]
#srr= args[3]
#output_dir=args[4]
#reg_indel_class_input = args[5]
#FJ_indel_class_input = args[6]
