#!/bin/sh

#  alignAndMakeReports.sh
#  THIS IS THE BIG WRAPPER
#
#  Created by Gillian Hsieh on 9/16/15.
#


ORIGDIR=${1}
OUTPUTDIR=${2}
NUMBASESAROUNDJUNC=${3} #default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70



rm ${1}genome/sorted*
rm ${1}reg/sorted*


sbatch -J aligning --mem=55000 --time=24:0:0 -e %err alignMkRprtsBowtie.sh ${1} ${2} ${3}
