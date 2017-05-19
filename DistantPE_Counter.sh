#!/bin/sh

#  DistantPE_Counter_genome_ENCODE.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#


## Because there are lot of repeat locations in the _distant_pairs.txt file generated above, the DistantPE_Counter_genome_ENCODE.py script is called by the shell of the same name to eliminate duplicate locations.  Another early problem was that the fasta generated later was too enormous to run in a timely fashion so the distant_pairs.txt file is split by this script into 24 smaller files based on the chromosome # of the upstream partner.
## The shell DistantPE_Counter_genome_ENCODE.sh takes in the FarJunction output directory and MACHETE installation directory and outputs <FJDir>/DistantPEFiles/<STEM>/chr1,2,3,4,...,X,Y_Distant_PE_frequency.txt

#  The output files chr?_Distant_PE_frequency.txt files contain three columns: chrA:M-N, chrB:P-Q, and R, where R is the number of times that these two exact windows were matched together.  R could be used to cull the fasta file if it gets too large, but at this point we are still looking for junctions between exons if only one read pair aligned discordantly.

FJDir=${1}  ## MACHETE output dir
INSTALLDIR=${2} ## MACHETE installation file (so python script can run)

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

mkdir -p ${1}DistantPEFiles/${STEM}

#module load python/2.7.9
ml python/2.7.5
python ${INSTALLDIR}DistantPE_Counter.py -d ${1} -s ${STEM}

echo "DistantPE_counter completed for ${STEM}" >> ${1}MasterError.txt