#!/bin/sh

#  PEfinder_genomeAndReg_ENCODE.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#


## finding mismatched paired end reads
## The shell PEfinder_genomeAndReg_ENCODE.sh takes the KNIFE alignment directory (OrigDir -1 ), the output directory (FJDir -2 ), the distance beyond which the user would consider alignments to be "discordant" (BP_distance -3 ), and the MACHETE installation directory (4) and calls a python script PEfinder_genomeAndReg_ENCODE.py.
## The python script identifies paired R1 and R2 from the genome and reg alignments and if the alignment location is > user defined bases apart then records them in an output file within the output directory: FarJunctionDirectory/DistantPEFiles/<STEM>_distant_pairs.txt


OrigDir=${1}  # KNIFE alignment dir
FJDir=${2} ## MACHETE output dir
BP_distance=${3} # user defined # base pairs after which reads would be considred discordant - currently using 100Kb
INSTALLDIR=${4} ## MACHETE's installation directory - so the python file can be found.

STEMFILE=${2}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

#module load python/2.7.9
ml load python/2.7.5

## If, for example, a read pair was found to be discordant, and R1= chrA:X, R2=chrB:Y, then the distant_pairs.txt file would contain the readID and chrA:M-N, chrB:P-Q where M-N is a window of 10,000 bases on each side of X and P-Q is a window of 10,000 bases on each side of Y. The window can be set below with the -w flag.

python ${INSTALLDIR}PEfinder.py -o ${1} -s ${STEM} -f ${2} -w 10000 -n ${3}

echo "PEfinder.sh step completed for ${STEM}" >> ${2}MasterError.txt

