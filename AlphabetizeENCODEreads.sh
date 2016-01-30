#!/bin/sh

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 8/18/15.
#

file=${1}
SORTEDFILE=${2}

head -n 2 ${file} > ${SORTEDFILE}
tail -n +3 ${file} | sort -k 1 >> ${SORTEDFILE}
