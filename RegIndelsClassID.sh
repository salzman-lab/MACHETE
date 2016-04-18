#!/bin/sh

#  RegIndelsClassID.sh
#  
#
#  Created by Gillian Hsieh on 3/17/16.
#


#Calls RegIndelsClassID.sh shell which calls RegIndels_ClassIDFile.py.
# Inputs are the regular indel alignments located at KNIFEdir/orig/RegIndelAlignments/<STEM>/unaligned_<STEM>_1/2_indel1,2,3,4,5.sam.  These are concatenated into a single file in the same directory called All_<STEM>_1/2_Regindels.sam.  In the concatenation step, like for the far junctions indels, any reads are omitted if they fail to overlie the junction by the user specified overlap, and also if a read aligns to multiple indel indices, the one with the best alignment score is put into the concatenated file.
## Then, the partner reads of all reads from All_<STEM>_1/2_Regindels.sam are identified and labeled as "good" or "bad".  If a read partner is found in genome, it has priority over transcriptome, which has priority over reg, and finally junc.  A far junction R2 cannot be found in another dictionary, as FJ reads are generated from previously unaligned reads.  If the read partner is in genome, it must be located on the same chromosome, on the opposite reference strand from R1.  If a read partner is in reg, then the downstream exon must be upstream of the uptstream reg indel exon, or the upstream read partner exon must be downstream of the downstream reg indel exon, on the same chromosome. Reference strands must be opposite.    If the partner is in junc, then the reg indel alignment must be located inside the circle created by the scrambled junction, and on the opposite reference strand.  In this manner, class input files are generated for the reg indels, which are located at KNIFE dir/circReads/<STEM>_output_RegIndel.txt

FJDir=${1} # FJ Dir
circpipedir=${2}
WINDOW=${3}
INSTALLDIR=${4}

origDir=${2}orig/
circReads=${2}circReads/

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

#### TESTING
#STEM=SRR1027188

ml load python/2.7.5

python ${INSTALLDIR}RegIndels_ClassIDFile.py -s ${STEM} -c ${circReads} -i ${origDir} -w ${WINDOW}