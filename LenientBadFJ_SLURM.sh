#!/bin/sh

#
#  
#
#  Created by Gillian Hsieh on 1/18/16.
#

FarJuncFasta=${1}
BadFJDir=${2}
INSTALLDIR=${3}

#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

ml load python/2.7.5
python ${3}SplitFastaforBadFJ.py -i ${FarJuncFasta} -l 40 -o ${BadFJDir}
