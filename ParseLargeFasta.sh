#!/bin/sh

#  ParseLargeFasta.sh
#  
#
#  Created by Gillian Hsieh on 10/11/16.
#
## This script takes the large FarJunctions.fa file and calls a script to parse it, for far junctions that were identified in the KNIFE unaligned reads

FJDIR=${1}
STEM=${2}
INSTALLDIR=${3}

ml load python/2.7.5
python ${INSTALLDIR}/filter_large_fasta.py -d ${FJDIR} -s ${STEM}