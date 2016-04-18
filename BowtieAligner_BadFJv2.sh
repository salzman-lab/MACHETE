#!/bin/sh

#  BowtieAligner.batch.sh
#  
#
#  Created by Gillian Hsieh on 11/10/15.
#

BOWTIEPARAM=${1}

module load bowtie/2.2.4


bowtie2 ${1}