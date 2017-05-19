#!/bin/sh

#  BadFJOnly.sh
#  
#
#  Created by Gillian Hsieh on 6/24/16.
#

## these are input parameters

OUTPUT_DIR=${1} # Output directory for MACHETE - does not need to exist already
REFGENOME=${2} # HG19 currently.  HG38 functionality to be added.


# OPTIONAL 6th input field is mode - owners, horence, etc
if [ $# -ge 3 ]
then
MODE=${3}
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
BadFJDir=${1}BadFJ/ # MACHETE far junction fasta entries that align to a KNIFE reference index
StemFile=${1}StemList.txt # file containing unique IDs of each experiment e.g. SRR#


# make temporary and output directories

mkdir -p ${BadFJDir}
mkdir -p ${1}err_and_out/

## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
#if [ -e ${2}StemList.txt ]
#then
#echo "using existing StemList.txt"
#else
#ml load python/2.7.5
#echo "generating StemList.txt from KNIFE output directory filenames"
#python ${INSTALLDIR}writeStemIDFiles.py -o ${ORIG_DIR} -f ${2}
#fi

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

FarJuncFasta=${1}fasta/${STEM}*_filtered_FarJunctions.fa
BadFJDir=${1}BadFJ/${STEM}/
BadFJver2Dir=${1}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}
mkdir -p ${BadFJver2Dir}
r1file=${BadFJver2Dir}${STEM}_filtered_FarJunctions_R1.fa
r2file=${BadFJver2Dir}${STEM}_filtered_FarJunctions_R2.fa


#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

j7_id=`sbatch -J BadFJ_Split ${RESOURCE_FLAG} --mem=55000 --nodes=4 --time=10:0:0 -o ${1}err_and_out/out_6BadJunc_split.txt -e ${1}err_and_out/err_6BadJunc_split.txt ${INSTALLDIR}LenientBadFJ_SLURM.sh ${FarJuncFasta} ${BadFJver2Dir} ${1} ${INSTALLDIR} | awk '{print $4}'`
echo "BadfJ ver 2 Split -- ${j7_id}"

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 --np 0 -p 4 --rdg 50,50 --rfg 50,50"


## submit SLURM jobs to do bowtie alignments for each of BadFJ indices above
if [ -e ${BadFJDir}${STEM}_BadFJtoGenome.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtoGenome.sam exists.  To realign, please manually delete this file first"
else
BadFJj1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 --nodes=4 ${RESOURCE_FLAG} --time=36:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoGenome.sam | awk '{print $4}'`
echo "BadFJ to genome: ${BadFJj1_id}"
depend_str7=${depend_str7}:${BadFJj1_id}
fi

if [ -e ${BadFJDir}${STEM}_BadFJtotranscriptome.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtotranscriptome.sam exists.  To realign, please manually delete this file first"
else
BadFJj2_id=`sbatch -J ${STEM}FJ_to_transcriptome --mem=55000 --nodes=4  ${RESOURCE_FLAG} --time=36:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtotranscriptome.sam | awk '{print $4}'`
echo "BadFJ to transcriptome: ${BadFJj2_id}"
depend_str7=${depend_str7}:${BadFJj2_id}
fi

if [ -e ${BadFJDir}${STEM}_BadFJtoReg.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtoReg.sam exists.  To realign, please manually delete this file first"
else
BadFJj3_id=`sbatch -J ${STEM}FJ_to_reg --mem=55000 --nodes=4  ${RESOURCE_FLAG} --time=36:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoReg.sam | awk '{print $4}'`
echo "BadFJ to reg: ${BadFJj3_id}"
depend_str7=${depend_str7}:${BadFJj3_id}
fi

if [ -e ${BadFJDir}${STEM}_BadFJtoJunc.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtoJunc.sam exists.  To realign, please manually delete this file first"
else
BadFJj4_id=`sbatch -J ${STEM}FJ_to_junc --mem=55000 --nodes=4 ${RESOURCE_FLAG} --time=36:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoJunc.sam | awk '{print $4}'`
echo "BadFJ to junc: ${BadFJj4_id}"
depend_str7=${depend_str7}:${BadFJj4_id}
fi

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
BadFJv2j1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 --nodes=4 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${genomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to genome: ${BadFJv2j1_id}"
depend_str7=${depend_str7}:${BadFJv2j1_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j2_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --nodes=4 --time=36:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${transcriptomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to transcriptome: ${BadFJv2j2_id}"
depend_str7=${depend_str7}:${BadFJv2j2_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtoReg.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoReg.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j3_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --nodes=4 --time=36:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${regBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to reg: ${BadFJv2j3_id}"
depend_str7=${depend_str7}:${BadFJv2j3_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoJunc.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j4_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --nodes=4 --time=36:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${juncBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to junc: ${BadFJv2j4_id}"
depend_str7=${depend_str7}:${BadFJv2j4_id}
fi
done
