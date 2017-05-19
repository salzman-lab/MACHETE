#!/bin/sh
#  makeJunctions.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#


#make Junction fasta file by extracting sequence info from pickles
## Pickle files are binary files used for storage of large amounts of information in Python.  Here, GTF files have been repackaged into pickle files and store sequence, exon name, and exon location information.  Accessing pickles can be time consuming so target locations from the discordant reads have been broken down by chromosome and alphabetized.
## Loop goes through integers 1-24 where that corresponds to chromosome # (23 = X, 24 = Y). For each of those chromosomes, the shell makeJunctions.sh calls makeJunctions.py
## MakeJunctions.py takes in FJDir/DistantPEFiles/sorted__chrA_Distant_PE_frequency.txt, and outputs FJDir/fasta/<STEM>/<STEM>_chrAFarJunctions.fa.  It reads in the discordant windows, and searches the pickle file for all exons names/exon locations/ exon sequences on either the sense or antisense strands within the discordant windows.  Then it makes all pairs of possible exon junctions. All pairs include fusions between two exons that are +/+, +/-, -/+, and -/-.  In the case that a sense and antisense exon are fused, then the sequence listed is the exact sequence that would occur if the fusion was read from the 5'->3' direction.  If the generated sequence was BLATted, the correct "strands" would appear in BLAT.  Similarly, for (-)/(-) exon pairs, if the generated sequence was BLATted, the exons would appear on the (-) strand in BLAT.
## THIS IS DIFFERENT IN KNIFE.  In KNIFE, if a (-)/(-) pair was BLATted, it would look as if it were +/+ because the KNIFE reverse complements (-) sequences.  Additionally in KNIFE, there is no way to detect inversions because -/+ and +/- fasta sequences are not generated.


PICKLEDIR=${1}  ## Python pickle storing sequence/ exon name / exon location info
FJFile=${2}  ## MACHETE output directory
Chromosome=${3}  ## the chromosome number of the input file
INSTALLDIR=${4} ## the directory where MACHETE is installed.

STEMFILE=${2}StemList.txt
FASTADIR=${2}fasta/
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

OUTPUTDIR=${FASTADIR}${STEM}/
mkdir -p ${OUTPUTDIR}
INPUTDIR=${2}DistantPEFiles/${STEM}/

ml load python/2.7.5
for file in ${INPUTDIR}/sorted_chr${3}_*; do
python ${INSTALLDIR}makeJunctions.py -p ${1} -f ${file} -o ${OUTPUTDIR} -s ${STEM}
done;


## FJDir/<STEM>_ChrA_FarJunction_duplicates.fa exists because when extracting the sequence/name/location from the pickles, it's faster just to generate the output than to search a dictionary for the previous existence of the same exon pair and sequence.  The python file then parses the output files and removes duplicate fasta entries.  After duplicates are removed, the intermediate ChrA_FarJunction_duplicates.fa file is also removed by the shell.
rm ${OUTPUTDIR}${STEM}_chr${3}_*_duplicates.fa


echo "makeJunctions.sh completed for ${STEM}" >> ${2}MasterError.txt
