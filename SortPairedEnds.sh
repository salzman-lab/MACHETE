#!/bin/sh

#  SortPairedEnds.sh
#  
#
#  Created by Gillian Hsieh on 1/7/16.
#

DISTANTPEFILE=${1}
DISTANTPESORTED=${2}

sort -k1,1.1 -k1,1.2 -k2,2 ${1} > ${2}