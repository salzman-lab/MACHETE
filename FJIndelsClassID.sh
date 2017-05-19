#!/bin/sh

#  FJIndelsClassID.sh
#  
#
#  Created by Gillian Hsieh on 3/17/16.

# FJ indels class output file
## This script calls FJIndelsClassID.sh
## This takes the FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/All_<STEM>_1/2_FJindels.sam and identifies read partners.  The same criteria to identify read partners as FarJuncNaiveReport.sh are applied.
## Output files are placed into FJ Dir/ GLM_classInput/<STEM>_output_FJIndels.txt

#
FJDir=${1} # FJ Dir
circpipedir=${2}
WINDOW=${3}
INSTALLDIR=${4}

origDir=${2}orig/
circReads=${2}circReads/

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

ml load python/2.7.5
python ${INSTALLDIR}FJIndels_ClassIDFile.py -s ${STEM} -c ${circReads} -f ${FJDir} -i ${origDir} -w ${WINDOW}

rm ${FJDir}/GLM_classInput/${STEM}_temp_output_FJIndels.txt

echo "FJIndelsClassID.sh complete for ${STEM} - check ${FJDir}GLM_classInput/${STEM}_output_FJIndels.txt" >> ${1}MasterError.txt