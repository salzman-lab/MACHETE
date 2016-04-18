#!/bin/sh

#  AppendNaiveRept.sh
#  
#
#  Created by Gillian Hsieh on 2/1/16.
#
### The AppendNaiveRept.sh shell calls the AppendNaiveRept.py script.  This reads in the IndelsHistogram, BadFJ and BadFJ_ver2 files, and GLM report results and outputs all the results into a single file in /FJDir/reports/AppendedReports/<STEM>_naive_report_Appended.txt
FarJuncDir=${1}
GLMReportDir=${2}
INSTALLDIR=${3}
FJGLMReportsDir=${4}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${1}reports/AppendedReports/

ml load python/2.7.5

python ${INSTALLDIR}AppendNaiveRept.py -f ${1} -g ${2} -s ${STEM} -G ${4}