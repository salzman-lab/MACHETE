#!/bin/sh

#  bowtieindexer.batch.sh
#  
#
#  Created by Gillian Hsieh on 9/22/15.
#

module load bowtie/2.2.4

bowtie2-build ${1} ${2}