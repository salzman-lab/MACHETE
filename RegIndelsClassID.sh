#!/bin/sh

#  RegIndelsClassID.sh
#  
#
#  Created by Gillian Hsieh on 3/17/16.
#
FJDir=${1} # FJ Dir
circpipedir=${2}
WINDOW=${3}
INSTALLDIR=${4}

origDir=${2}orig/
circReads=${2}circReads/

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

#### TESTING
#STEM=SRR1027188

ml load python/2.7.5

python ${INSTALLDIR}RegIndels_ClassIDFile.py -s ${STEM} -c ${circReads} -i ${origDir} -w ${WINDOW}
