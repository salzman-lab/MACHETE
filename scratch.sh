#!/bin/sh

#  scratch.sh
#  
#
#  Created by Gillian Hsieh on 4/18/16.
#

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
StemFile=${2}StemList.txt # file containing unique IDs of each experiment e.g. SRR#

# make temporary and output directories
mkdir -p ${OUTPUT_DIR}
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
mkdir -p ${SECONDFARJUNCDIR}AlignedIndels/RemoveNonOverlap/
mkdir -p ${2}IndelsHistogram/
mkdir -p ${2}reports/AppendedReports/
mkdir -p ${GLM_DIR}AppendGLM/

## remove the old error files
rm ${2}err_and_out/*

## loads python for sherlock.  on SCG3 python/2.7.9 is installed
ml load python/2.7.5

## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
python ${INSTALLDIR}writeStemIDFiles.py -o ${ORIG_DIR} -f ${2}
NUM_FILES=$((`more ${StemFile} | wc -l`))
echo ${NUM_FILES}

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
BadFJver2Dir=${1}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}
mkdir -p ${BadFJver2Dir}
r1file=${BadFJver2Dir}${STEM}_FarJunctions_R1.fa
r2file=${BadFJver2Dir}${STEM}_FarJunctions_R2.fa


#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

j7_id=`sbatch -J BadFJ_Split ${RESOURCE_FLAG} --mem=55000 --nodes=4 --time=10:0:0 -o ${2}err_and_out/out_6BadJunc_split.txt -e ${2}err_and_out/err_6BadJunc_split.txt ${INSTALLDIR}LenientBadFJ_SLURM.sh ${FarJuncFasta} ${BadFJver2Dir} ${INSTALLDIR} | awk '{print $4}'`

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --np 0 --rdg 50,50 --rfg 50,50"


## submit SLURM jobs to do bowtie alignments for each of BadFJ indices above
BadFJj1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoGenome.sam | awk '{print $4}'`

BadFJj2_id=`sbatch -J ${STEM}FJ_to_transcriptome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtotranscriptome.sam | awk '{print $4}'`

BadFJj3_id=`sbatch -J ${STEM}FJ_to_reg --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoReg.sam | awk '{print $4}'`

BadFJj4_id=`sbatch -J ${STEM}FJ_to_junc --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoJunc.sam | awk '{print $4}'`
depend_str7=${depend_str7}:${BadFJj1_id}:${BadFJj2_id}:${BadFJj3_id}:${BadFJj4_id}

## Read gaps are disallowed in the first version of BadJuncs.  A second version of BadJuncs was created to also find genome/reg/junc/transcriptome alignments with gapped alignments.
## For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.

genomeBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${genomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam"
transcriptomeBOWTIEPARAM="--no-unal --no-mixed  --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${transcriptomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam"
regBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${regIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoReg.sam"
juncBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${juncIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam"


## submit SLURM jobs to do the bowtie alignment for each of the BadFJ Ver 2 indices above
BadFJv2j1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${genomeBOWTIEPARAM}" | awk '{print $4}'`

BadFJv2j2_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${transcriptomeBOWTIEPARAM}" | awk '{print $4}'`

BadFJv2j3_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${regBOWTIEPARAM}" | awk '{print $4}'`

BadFJv2j4_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=5:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${juncBOWTIEPARAM}" | awk '{print $4}'`

depend_str7=${depend_str7}:${BadFJv2j1_id}:${BadFJv2j2_id}:${BadFJv2j3_id}:${BadFJv2j4_id}

done


echo "BadFJ alignments:  ${depend_str7} -- check for FJDir/BadFJ/<STEM>/* and FJDir/BadFJ_ver2/<STEM>/ alignment sam files"
