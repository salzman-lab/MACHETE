#!/bin/sh

# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for


CIRCPIPE_DIR=${1} #Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)
OUTPUT_DIR=${2} # Far Junc Dir
USERBPDIST=${3} #  using 100000 (100Kb) so far for testing purposes
REFGENOME=${4} # HG19 vs HG38  currently using /scratch/PI/horence/gillian/HG19exons, although could upgrade to HG38
NUMBASESAROUNDJUNC=${5} #default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70
NumIndels=${6} # current default = 5, arbitrary
NumBPOverlapAtJunc=${7} # this is different than input #5 because it is for after the indel alignment.  Current default 5+8 (8 more than the max # indels per side if default # indels is 5) = 13 ; this is arbitrary
# OPTIONAL 8 is mode - owners, horence, etc

if [ $# -ge 8 ]
then
MODE=${8}
else
MODE=sam
fi

## REPLACE THESE THREE FIELDS AFTER INSTALLATION
INSTALLDIR="/scratch/PI/horence/gillian/MACHETE/"
CIRCREF="/share/PI/horence/circularRNApipeline_Cluster/index/"
REG_INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/"
####


if [[ "$REFGENOME" = *HG38* ]]
then
PICKLEDIR="/scratch/PI/horence/grch38_junctions/"
fi

if [[ "$REFGENOME" = *HG19* ]]
then
PICKLEDIR="/scratch/PI/horence/gillian/HG19exons/"
## REPLACE WITH PYTHON PICKLE DIRECTORY
fi


if [[ "$MODE" = *owners* ]]
then
RESOURCE_FLAG="-p owners"
fi

#if [[ "$MODE" = *horence* ]]
#then
#RESOURCE_FLAG="-p horence"
#fi

ORIG_DIR=${1}orig/
GLM_DIR=${1}circReads/glmReports/
DistantPEDir=${2}DistantPEFiles/
FASTADIR=${2}fasta/
BOWTIE_DIR=${2}BowtieIndex/
UNALIGNEDDIR=${ORIG_DIR}unaligned/
FARJUNCDIR=${2}FarJunctionAlignments/
SECONDFARJUNCDIR=${2}FarJuncSecondary/
BadFJDir=${2}BadFJ/
StemFile=${2}StemList.txt


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



rm ${2}err_and_out/*


ml load python/2.7.5
python ${INSTALLDIR}writeStemIDFiles.py -o ${ORIG_DIR} -f ${2}

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
NUM_FILES=$((`more ${StemFile} | wc -l`))
echo ${NUM_FILES}
#
#
rm ${ORIG_DIR}reg/sorted*.sam
rm ${ORIG_DIR}genome/sorted*.sam
## sorting reg files

j1_id=`sbatch -J sortingReg ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_1sortReg.txt -e ${2}err_and_out/err_1sortReg.txt ${INSTALLDIR}AlphabetizeENCODEreads.sh ${ORIG_DIR}reg/ ${StemFile} | awk '{print $4}'`

echo "sorting reg files: job# ${j1_id}"


# sorting genome files
j2_id=`sbatch -J sortingGenome ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_1sortGenome.txt -e ${2}err_and_out/err_1sortGenome.txt ${INSTALLDIR}AlphabetizeENCODEreads.sh ${ORIG_DIR}genome/ ${StemFile} | awk '{print $4}'`

echo "sorting genome files: job# ${j2_id}"
#
#
## finding mismatched paired end reads


j3_id=`sbatch -J PEMismatch ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_2PEfinder.txt -e ${2}err_and_out/err_2PEfinder.txt --depend=afterok:${j1_id}:${j2_id} ${INSTALLDIR}PEfinder_genomeAndReg_ENCODE.sh ${ORIG_DIR} ${2} ${3} ${INSTALLDIR} | awk '{print $4}'`

echo "Outputting Mismatched paired ends: job ${j3_id}"

### window around position in genome is 10,000

### counting rates of mismatched PE

#
j4_id=`sbatch -J CountMismatch ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_3PEcounter.txt -e ${2}err_and_out/err_3PEcounter.txt --depend=afterok:${j3_id} ${INSTALLDIR}DistantPE_Counter_genome_ENCODE.sh ${2} ${INSTALLDIR} | awk '{print $4}'`
echo "counting mismatch rates: ${j4_id}"

# sort mismatched PE by chromosome

j5_id=`sbatch -J SortPE ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_4PEsort.txt -e ${2}err_and_out/err_4PEsort.txt --depend=afterok:${j4_id} ${INSTALLDIR}SortPairedEnds.sh ${2} | awk '{print $4}'`

echo "sorting mismatched PE files: ${j5_id}"

#make Junction fasta file by extracting sequence info from pickles

depend_str6="--depend=afterok"
for c in `seq 1 24`; do
A=${c}
    if [ ${c} -eq 23 ]; then A=X; fi
    if [ ${c} -eq 24 ]; then A=Y; fi

    j6_id=`sbatch -J MakeFJFasta ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_5makefasta.txt -e ${2}err_and_out/err_5makefasta.txt --depend=afterok:${j5_id} ${INSTALLDIR}makeJunctions.sh ${PICKLEDIR} ${2} ${A} ${INSTALLDIR} | awk '{print $4}'`
    depend_str6=${depend_str6}:${j6_id}
done

echo "make fusion fasta files: ${depend_str6}"


##
##
##make single FJ fasta from all the fastas and then call bowtie indexer
##
#
j6a_id=`sbatch -J FJ_Index ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_5FJIndexing.txt -e ${2}err_and_out/err_5FJIndexing.txt ${depend_str6} ${INSTALLDIR}linkfastafiles.sh ${2} | awk '{print $4}'`
echo "make FJ bowtie indices for each experiment: ${j6a_id}"

##
##
# make BadJunc directory --  bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align
#

j7_id=`sbatch -J BadJuncs ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_6BadJunc.txt -e ${2}err_and_out/err_6BadJunc.txt --depend=afterok:${j6a_id} ${INSTALLDIR}GoodvsBadFJ_SLURM.sh ${2} ${4} ${INSTALLDIR} ${CIRCREF} ${RESOURCE_FLAG} | awk '{print $4}'`
echo "Identify Bad FJ's: ${j7_id}"


# align unaligned files to the FJ bowtie index
#
j8_id=`sbatch -J AlignFJ ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_7AlignFJ.txt -e ${2}err_and_out/err_7AlignFJ.txt --depend=afterok:${j6a_id} ${INSTALLDIR}AlignUnalignedtoFJ.sh ${2} ${ORIG_DIR} | awk '{print $4}'`

echo "align unaligned reads to FJ index: ${j8_id}"
##
#
#
###make FJ naive report
#
j9_id=`sbatch -J FJNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_8NaiveRpt.txt -e ${2}err_and_out/err_8NaiveRpt.txt --depend=afterok:${j8_id} ${INSTALLDIR}FarJuncNaiveReport.sh ${2} ${ORIG_DIR} ${5} ${INSTALLDIR} | awk '{print $4}'`

echo "make naive rpt - ${j9_id}"
###
##
####### GLM #######
##Make Class input files for GLM
##
###
j15a_id=`sbatch -J AddFJtoIDfile ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=5000 --nodes=1 --time=1:0:0 -o ${2}err_and_out/out_15FJforGLM.txt -e ${2}err_and_out/err_15FJforGLM.txt --depend=afterok:${j9_id} ${INSTALLDIR}parse_FJ_ID_for_GLM.sh ${2} "${1}circReads/" | awk '{print $4}'`

echo "make FJ class input files: ${j15a_id}"

#
##
###
###### ESTIMATE LIGATION ARTIFACT
#
###
###make FarJunctions.fa => Indel.fa files
####
j10_id=`sbatch -J MakeFJIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_10FJIndels.txt -e ${2}err_and_out/err_10FJIndels.txt --depend=afterok:${j6a_id} ${INSTALLDIR}MakeIndelFiles.sh ${2} ${6} ${INSTALLDIR} | awk '{print $4}'`

echo "make indel files: ${j10_id}"
#
#

##
#
#
# make Bowtie Indices for Far Junc Indel files
#
depend_str11="--depend=afterok"

START=1
for (( c=$START; c<=${6}; c++ ))
do
j11_id=`sbatch -J IndexFJIndel ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_11indexindels.txt -e ${2}err_and_out/err_11indexindels.txt --depend=afterok:${j10_id} ${INSTALLDIR}BowtieIndexFJIndels.sh ${2}FarJuncIndels/ ${c} ${2}BowtieIndels/ ${2} | awk '{print $4}'`
depend_str11=${depend_str11}:${j11_id}

done
echo "index indels ${depend_str11}"



##Align FarJuncSecondary (unaligned to FJ index) to FJ indels


#
##
depend_str13="--depend=afterok"
START=1
for (( c=$START; c<=${6}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
j13_id=`sbatch -J AlignIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_12alignindels.txt -e ${2}err_and_out/err_12alignindels.txt ${depend_str11} ${INSTALLDIR}BowtieAlignFJIndels.sh ${2} "${BOWTIEPARAMETERS}" ${c} | awk '{print $4}'`
    depend_str13=${depend_str13}:${j13_id}

done
echo "align to indels: ${depend_str13}"
#

#
## loop through AlignedIndels directory
## things that don't overlap indels junction by ${7} are removed.
## for every junction, a string is created.  For example, if 3 indel files exist, the string [0, 0, 2, 1, 5, 0, 1] represents 0*3 deletions, 0* 2 deletions, 2 * 1 deletion, 1 * no indels, 5 * 1 insertion, 0 * 2 insertions, 1 * 3 insertions.
## strings are output into IndelsHistogram folder
#
##
j14_id=`sbatch -J FilterIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_13filterIndels.txt -e ${2}err_and_out/err_13filterIndels.txt ${depend_str13} ${INSTALLDIR}FindAlignmentArtifact_SLURM.sh ${2} ${7} ${6} ${INSTALLDIR}| awk '{print $4}'`

echo "make indels histo: ${j14_id}"

###
####
##### REG INDELS ################
##
# Align unaligned files to the expanded reg junctions with indels
AlignedIndels=${1}/orig/RegIndelAlignments/
mkdir -p ${AlignedIndels}

depend_str16="--depend=afterok"
START=1
for (( c=$START; c<=${6}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
j16_id=`sbatch -J AlignRegIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_15AlignRegIndels.txt -e ${2}err_and_out/err_15AlignRegIndels.txt ${INSTALLDIR}AlignUnalignedtoRegIndel.sh ${1} ${c} ${2} "${BOWTIEPARAMETERS}" ${REG_INDEL_INDICES} | awk '{print $4}'`
    depend_str16=${depend_str16}:${j16_id}
done
echo "Aligning unaligned files to linear junc indels: ${depend_str16}"


#

####
####
###
####
######## MAKE REG AND FJ INDELS CLASS OUTPUT FILES ###########
####
###
### reg indels class output file
#
j18_id=`sbatch -J RegIndelClass ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_18RegIndelsClassOutput.txt -e ${2}err_and_out/err_18RegIndelsClassOutput.txt ${depend_str16} ${INSTALLDIR}RegIndelsClassID.sh ${2} ${1} ${7} ${INSTALLDIR} | awk '{print $4}'`
echo " Reg Indels Class Output: ${j18_id}"

# FJ indels class output file

#
j19_id=`sbatch -J FJIndelClass ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_19FJIndelsClassOutput.txt -e ${2}err_and_out/err_19FJIndelsClassOutput.txt --depend=afterok:${j14_id} ${INSTALLDIR}FJIndelsClassID.sh ${2} ${1} ${7} ${INSTALLDIR} | awk '{print $4}'`
echo "FJ Indels Class Output: ${j19_id}"


###### RUN GLM ###########################
#
## Run GLM
##
j15b_id=`sbatch -J GLM.r ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=1:0:0 -o ${2}err_and_out/out_15GLM_r.txt -e ${2}err_and_out/err_15GLM_r.txt --depend=afterok:${j15a_id}:${j18_id}:${j19_id} ${INSTALLDIR}run_GLM.sh ${1} ${2} ${INSTALLDIR} | awk '{print $4}'`

echo "Run GLM: ${j15b_id}"


###
###

## Append linear junctions GLM report with anomalies, indels
#
j17_id=`sbatch -J AppendRegGLM ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=2 --time=24:0:0 -o ${2}err_and_out/out_17AppendRegGLM.txt -e ${2}err_and_out/err_17AppendGLM.txt ${depend_str16} ${INSTALLDIR}AddIndelstoGLM.sh ${1} ${2} ${7} | awk '{print $4}'`
echo " Appending linearJuncs GLM report: ${j17_id}"
##
####

### Append Naive report:  Add Indels, quality of junctions (good/bad), frequency of junction participation in linear junctions, and GLM report fields to the Naive report
###
j15_id=`sbatch -J AppendNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_14AppendRpt.txt -e ${2}err_and_out/err_14AppendRpt.txt --depend=afterok:${j7_id}:${j9_id}:${j15b_id} ${INSTALLDIR}AppendNaiveRept.sh ${2} ${GLM_DIR} ${INSTALLDIR} | awk '{print $4}'`

echo "append naive rpt: ${j15_id}"
##
