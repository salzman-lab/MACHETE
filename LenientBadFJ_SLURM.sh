#!/bin/sh

#
#  
#
#  Created by Gillian Hsieh on 1/18/16.
#

FarJuncFasta=${1}
BadFJDir=${2}
SPORKDir=${3}
INSTALLDIR=${4}

#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

ml load python/2.7.5
python ${4}SplitFastaforBadFJ.py -i ${FarJuncFasta} -l 40 -o ${BadFJDir}

echo "Fasta sequences split into 40bp reads - LenientBadFJ_SLURM.sh complete for ${STEM}" >> ${3}MasterError.txt
echo "Unable to tell if BadFJ and BadFJ ver2 alignments complete. Please check ${1}BadFJ/${STEM} and ${1}/BadFJ_ver2/${STEM} to see if empty if Appended Reports fail without other failures" >> ${3}MasterError.txt