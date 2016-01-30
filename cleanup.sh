#!/bin/sh

#  cleanup.sh
#  
#
#  Created by Gillian Hsieh on 1/7/16.
#

rm ${1}
rm ${2}
rm ${3}
rm ${4}
mv ${5} ${6}
rm ${7}


# What do those commands above correlate to?
#rm ${1}genome/sorted*
#rm ${1}reg/sorted*
#rm ${DISTANTPEFILE}
#rm ${2}FarJunctions_duplicates.fa
#mv ${2}*distant_pairs.txt ${OUTPUT_DIR}distant_pairs/
#rm ${BOWTIEDIR}*.txt
