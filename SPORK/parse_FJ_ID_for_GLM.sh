#!/bin/sh

#  parse_FJ_ID_for_GLM.sh
#  
#
#  Created by Gillian Hsieh on 2/22/16.
#

##Make fJ Class input files for GLM
##parse_FJ_ID_for_GLM.sh is a simple shell script that takes the ID files generated above in FJDir/reports/IDs_<STEM>.txt as inputs and removes any unmapped or unaligned read IDs using the "grep -v" command.  The outputs are fed into FJDir/GLM_classInput/<STEM>__output_FJ.txt with the other class input files for the GLM.
###

FJDir=${1} # MACHETE output dir

GLM_class_inputs=${1}GLM_classInput/
mkdir ${GLM_class_inputs}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

FJ_Outputfile=${GLM_class_inputs}${STEM}_1__output_FJ.txt

for file in ${1}reports/IDs*${STEM}*.txt
do
sed '1d' ${file} | grep -v "Unmapped\|unaligned" > ${FJ_Outputfile}
done

echo "parse_FJ_ID_for_GLM.sh complete for ${STEM} - check for ${GLM_class_inputs}${STEM}__output_FJ.txt" >> ${1}MasterError.txt