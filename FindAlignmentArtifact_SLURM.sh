#!/bin/sh

#  FindAlignmentArtifact.sh
#  
#
#  Created by Gillian Hsieh on 1/6/16.
#

## Calls FindAlignmentArtifact_SLURM.sh which is a shell that calls MakeIndelsHisto.py.  The MakeIndelsHisto.py script reads in the aligned indels from the sam files FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam.  It concatenates all the indels_1,2,3,4,5.sam files into a larger file All_<STEM>_1/2_indels.sam in the same directory as the original sam files. In the All_<STEM>_1/2_indels.sam files, any indels that did not overlap the junction by the user specified # base pairs around the breakpoint are removed. Additionally, since the FarJuncSecondary files were aligned independently to the indels_1-5 bowtie indices, the same read could align to multiple indices.  Therefore for each read, the read with the best alignment score is placed in the All_<STEM>_1/2_indels sam file, and all other alignments to other indices are discarded.
#  Then the python script checks the FJdir/FarJunctionAlignments/<STEM>/ sam files and creates an array for each junction.  The array is of length 2*#indels+1. In the case of 5 indels, the length is 11 and represents the number reads that aligned to the junction with [5Del, 4Del, 3Del, 2Del, 1Del, aligned to junction exactly, 1Ins, 2Ins, 3Ins, 4Ins, 5Ins]
## the junction name and this array are output to FJDir/IndelsHistogram/indels_<STEM>_1/2.txt.  These outputs will be used to generate the Appended reports later.

FJFile=${1}
NumBParoundJunc=${2}
NumIndels=${3}
INSTALLDIR=${4}

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${1}FarJuncSecondary/AlignedIndels/RemoveNonOverlap/${STEM}/

module load python/2.7.5

python ${INSTALLDIR}MakeIndelsHisto.py -s ${STEM} -f ${1} -w ${2} -x ${3}