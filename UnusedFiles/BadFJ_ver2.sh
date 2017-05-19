#!/bin/sh

#  BadFJOnly.sh
#  
#
#  Created by Gillian Hsieh on 6/24/16.
#

## these are input parameters

CIRCPIPE_DIR=${1} #Main circular RNA pipeline output directory -- contains subfolders "orig", "circReads", "logs", and "sample stats"
OUTPUT_DIR=${2} # Output directory for MACHETE - does not need to exist already
USERBPDIST=${3} # min # base pairs between two exons that can be considered a "fusion".  In testing, we used 100000.
REFGENOME=${4} # HG19 currently.  HG38 functionality to be added.
NUMBASESAROUNDJUNC=${5} #Number of bases that need to surround the junction for it to be considered a true alignment.  Default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70
NumIndels=5  ## if we want to use more than 5 insertions or 5 deletions per side, this can be changed but I hardcoded it in on 4/15/16


# OPTIONAL 6th input field is mode - owners, horence, etc
if [ $# -ge 6 ]
then
MODE=${6}
else
MODE=sam
fi
if [[ "$MODE" = *owners* ]]
then
RESOURCE_FLAG="-p owners"
fi
if [[ "$MODE" = *horence* ]]
then
RESOURCE_FLAG="-p horence"
fi


## REPLACE THESE THREE REFERENCE FIELDS AFTER INSTALLATION
## these were hard coded but when Sweet-Cordero lab was trying to install our scripts, I thought they might want to install on a different server.  If my paths were still sherlock paths, then none of the shells would work.

## location of MACHETE installation.  This is important for sub-shell scripts so they know where to find python scripts that are called.
INSTALLDIR="/scratch/PI/horence/gillian/MACHETE/"
## location of reference indices for KNIFE - these would only change based on HG38 vs HG19
CIRCREF="/share/PI/horence/circularRNApipeline_Cluster/index/"
## location of reference index for the reg indels - these never change
REG_INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/"
####

## pickle files are binary python files that take up less space.  These particular pickles contain information about reference exon locations and sequences. These are packaged from a GTF file.

## please enter path to "reference pickle" to HG38.  Only necessary if using HG38
## OF NOTE - these may not actually exist.  Currently no HG38 functionality
#if [[ "$REFGENOME" = *HG38* ]]
#then
#PICKLEDIR="/scratch/PI/horence/grch38_junctions/"
#fi

## please enter path to "reference pickle" for HG19. Only necessary if using HG19
if [[ "$REFGENOME" = *HG19* ]]
then
PICKLEDIR="/scratch/PI/horence/gillian/HG19exons/"
fi



## define variables
ORIG_DIR=${1}orig/  # KNIFE alignments
UNALIGNEDDIR=${ORIG_DIR}unaligned/ # KNIFE unaligned reads
GLM_DIR=${1}circReads/glmReports/ # KNIFE GLM reports
DistantPEDir=${2}DistantPEFiles/ # new directory that will contain mismatched paired end reads
FASTADIR=${2}fasta/ # MACHETE far junctions fasta dir
BOWTIE_DIR=${2}BowtieIndex/ # MACHETE bowtie indices for far junctions fasta
FARJUNCDIR=${2}FarJunctionAlignments/ # unaligned reads that aligned to FJ Bowtie Index
SECONDFARJUNCDIR=${2}FarJuncSecondary/ # unaligned reads that did not align to FJ Bowtie Index
BadFJDir=${2}BadFJ/ # MACHETE far junction fasta entries that align to a KNIFE reference index
StemFile=${2}BadFJ2StemList.txt # file containing unique IDs of each experiment e.g. SRR#

# make temporary and output directories
mkdir -p ${FASTADIR}
mkdir -p ${2}reports/
mkdir -p ${2}err_and_out/
mkdir -p ${BOWTIE_DIR}
mkdir -p ${FARJUNCDIR}
mkdir -p ${SECONDFARJUNCDIR}
mkdir -p ${BadFJDir}
mkdir -p ${DistantPEDir}
mkdir -p ${2}BowtieIndels/
mkdir -p ${2}FarJuncIndels/
mkdir -p ${SECONDFARJUNCDIR}AlignedIndels/
mkdir -p ${2}IndelsHistogram/
mkdir -p ${2}reports/AppendedReports/
mkdir -p ${GLM_DIR}AppendGLM/
mkdir -p ${2}GLM_classInput/




## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
if [ -e ${2}StemList.txt ]
then
echo "using existing StemList.txt"
else
ml load python/2.7.5
echo "generating StemList.txt from KNIFE output directory filenames"
python ${INSTALLDIR}writeStemIDFiles.py -o ${ORIG_DIR} -f ${2}
fi

NUM_FILES=$((`more ${StemFile} | wc -l`))
echo ${NUM_FILES}

##
## If there is homology between a FarJunctions fasta sequence and the genome or transcriptome or a linear junction or circular junction, then the fusion read is less likely.  Alignments of the FarJunctions fasta sequences to the KNIFE reference indices, genome, transcriptome, linear junctions (reg), and scrambled junctions (junc) are created with two different bowtie parameters.  Bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align. These are just aligning the FJ Fasta to the bad juncs with various alignment parameters. Any junctions aligning to here will eventually be tagged as "BadFJ=1" in the final reports whereas if junctions don't align, they will receive a "BadFJ=0" in the final reports.

if [[ "$REFGENOME" = *HG19* ]]
then
genomeIndex=${CIRCREF}hg19_genome
transcriptomeIndex=${CIRCREF}hg19_transcriptome
regIndex=${CIRCREF}hg19_junctions_reg
juncIndex=${CIRCREF}hg19_junctions_scrambled
fi

depend_str7="--depend=afterok"

START=1
for (( c=$START; c<=${NUM_FILES}; c++ ))
do
STEM=`awk 'FNR == '${c}' {print $1}' ${StemFile}`

FarJuncFasta=${2}fasta/${STEM}*FarJunctions.fa
BadFJDir=${2}BadFJ/${STEM}/
BadFJver2Dir=${2}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}
mkdir -p ${BadFJver2Dir}
r1file=${BadFJver2Dir}${STEM}_FarJunctions_R1.fa
r2file=${BadFJver2Dir}${STEM}_FarJunctions_R2.fa


#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

j7_id=`sbatch -J BadFJ_Split ${RESOURCE_FLAG} --mem=10000 --nodes=4 --time=10:0:0 -o ${2}err_and_out/out_6BadJunc_split.txt -e ${2}err_and_out/err_6BadJunc_split.txt ${INSTALLDIR}LenientBadFJ_SLURM.sh ${FarJuncFasta} ${BadFJver2Dir} ${2} ${INSTALLDIR} | awk '{print $4}'`
echo "BadfJ ver 2 Split -- ${j7_id}"
# Read gaps are disallowed in the first version of BadJuncs.  A second version of BadJuncs was created to also find genome/reg/junc/transcriptome alignments with gapped alignments.
## For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.

genomeBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 --np 0 -I 0 -X 50000 -f --ff -x ${genomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam"
transcriptomeBOWTIEPARAM="--no-unal --no-mixed  --no-sq -p 8 --np 0 -I 0 -X 50000 -f --ff -x ${transcriptomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam"
regBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 --np 0 -I 0 -X 50000 -f --ff -x ${regIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoReg.sam"
juncBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 --np 0 -I 0 -X 50000 -f --ff -x ${juncIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam"


## submit SLURM jobs to do the bowtie alignment for each of the BadFJ Ver 2 indices above
if [ -e ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoGenome.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j1_id=`sbatch -J ${STEM}FJ_to_genome --mem=200000 --nodes=8 ${RESOURCE_FLAG} --time=48:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${genomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to genome: ${BadFJv2j1_id}"
depend_str7=${depend_str7}:${BadFJv2j1_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j2_id=`sbatch -J ${STEM}FJ_to_genome --mem=200000 ${RESOURCE_FLAG} --nodes=8 --time=48:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${transcriptomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to transcriptome: ${BadFJv2j2_id}"
depend_str7=${depend_str7}:${BadFJv2j2_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtoReg.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoReg.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j3_id=`sbatch -J ${STEM}FJ_to_genome --mem=200000 ${RESOURCE_FLAG} --nodes=8 --time=48:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${regBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to reg: ${BadFJv2j3_id}"
depend_str7=${depend_str7}:${BadFJv2j3_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoJunc.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j4_id=`sbatch -J ${STEM}FJ_to_genome --mem=200000 ${RESOURCE_FLAG} --nodes=8 --time=48:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${juncBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to junc: ${BadFJv2j4_id}"
depend_str7=${depend_str7}:${BadFJv2j4_id}
fi
done
