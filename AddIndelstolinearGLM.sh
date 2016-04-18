#!/bin/sh

#  AddIndelstoGLM.sh
#  
#
#  Created by Gillian Hsieh on 2/17/16.
#
## Append linear junctions GLM report with anomalies, indels
#AddIndelstolinearGLM.sh calls the python script KNIFEglmReportsForMachete.py.  This script parses circular and linear glmReports.  For the linear glmReports from KNIFE, the script collects any junctions where 1) the two exons from the linear report are from different genes or 2) the posterior probability is >0.9.  It adds on the rate of anomaly reads and indels to the reports and feeds them into FJDir/reports/AppendedReports.  For ciruclar reports, the script collects any junctions where the posterior probability is <0.9, appends the "Decoy" rate, and feeds the reports into FJDir/reports/Appended reports.
## The purpose of this script is to place all reports in a single directory for the user.

CirclePipeDir=${1}
FJDir=${2}
INSTALLDIR=${3}

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

ml load python/2.7.5
python ${3}KNIFEglmReportsForMachete.py -c ${1} -f ${2} -s ${STEM}