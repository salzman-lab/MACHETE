#!/bin/sh



## Put rob's fasta in a new directory fasta/<STEM>_SPORK_FarJunctions.fa
## then run MACHETE


# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for


CIRCPIPE_DIR=${1} #Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)
OUTPUT_DIR=${2} # SPORK directory
USERBPDIST=${3} #  using 100000 (100Kb) so far for testing purposes, NOT used in SPORK
REFGENOME=${4} # HG19 vs HG38  currently using /scratch/PI/horence/gillian/HG19exons, although could upgrade to HG38
NUMBASESAROUNDJUNC=${5} #default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70
NumIndels=5 # current default = 5, arbitrary


if [ $# -ge 6 ]
then
MODE=${6}
else
MODE=sam
fi

## REPLACE THESE THREE FIELDS AFTER INSTALLATION
INSTALLDIR="/scratch/PI/horence/gillian/MACHETE/"
CIRCREF="/share/PI/horence/circularRNApipeline_Cluster/index/"
REG_INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/"
SPORKINSTALLDIR=${INSTALLDIR}SPORK/
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
mkdir -p ${2}BowtieIndels/
mkdir -p ${2}FarJuncIndels/
mkdir -p ${SECONDFARJUNCDIR}AlignedIndels/
mkdir -p ${2}IndelsHistogram/
mkdir -p ${2}reports/AppendedReports/
mkdir -p ${GLM_DIR}AppendGLM/
mkdir -p ${2}GLM_classInput/


if [ "$(ls -A ${2}/err_and_out/*)" ]
then
echo "old error files exist -- removing"
rm ${2}err_and_out/*
fi

if [ -e ${2}MasterError.txt ]
then
rm ${2}MasterError.txt
fi

if [ -e ${2}StemList.txt ]
then
echo "using existing StemList.txt"
else
ml load python/2.7.5
echo "generating StemList.txt from KNIFE output directory filenames"
python ${INSTALLDIR}writeStemIDFiles.py -o ${ORIG_DIR} -f ${2}
fi

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
NUM_FILES=$((`more ${StemFile} | wc -l`))
echo ${NUM_FILES}
#
#
## GENERATE SPORK FASTAS
j1_id=`sbatch -J SporkFa ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=12:0:0 -o ${2}err_and_out/out_1SporkFa.txt -e ${2}err_and_out/err_1SporkFa.txt ${INSTALLDIR}generateSPORKfasta.sh ${1} ${2} ${5} ${SPORKINSTALLDIR} | awk '{print $4}'`
echo "make SPORK fastas: ${j1_id}"
###
##


###call bowtie indexer
##
#
j6a_id=`sbatch -J FJ_Index ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_5FJIndexing.txt -e ${2}err_and_out/err_5FJIndexing.txt --depend=afterok:${j1_id} ${INSTALLDIR}Spork_BowtieIndex.sh ${2} | awk '{print $4}'`
echo "make SPORK bowtie indices for each experiment: ${j6a_id}"

##
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

SPORKFasta=${2}fasta/${STEM}/*SPORK_Junctions.fa
BadFJDir=${2}BadFJ/${STEM}/
BadFJver2Dir=${2}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}
mkdir -p ${BadFJver2Dir}
r1file=${BadFJver2Dir}${STEM}*_R1.fa
r2file=${BadFJver2Dir}${STEM}*_R2.fa


#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

j7_id=`sbatch -J BadFJ_Split ${RESOURCE_FLAG} --mem=55000 --nodes=4 --time=10:0:0 -o ${2}err_and_out/out_6BadJunc_split.txt -e ${2}err_and_out/err_6BadJunc_split.txt --depend=afterok:${j6a_id} ${INSTALLDIR}LenientBadFJ_SLURM.sh ${SPORKFasta} ${BadFJver2Dir} ${2} ${INSTALLDIR} | awk '{print $4}'`
echo "BadfJ ver 2 Split -- ${j7_id}"

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --np 0 --rdg 50,50 --rfg 50,50"


## submit SLURM jobs to do bowtie alignments for each of BadFJ indices above
BadFJj1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${SPORKFasta} ${BadFJDir}${STEM}_BadFJtoGenome.sam | awk '{print $4}'`
echo "BadFJ to genome: ${BadFJj1_id}"

BadFJj2_id=`sbatch -J ${STEM}FJ_to_transcriptome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${SPORKFasta} ${BadFJDir}${STEM}_BadFJtotranscriptome.sam | awk '{print $4}'`
echo "BadFJ to transcriptome: ${BadFJj2_id}"


BadFJj3_id=`sbatch -J ${STEM}FJ_to_reg --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${SPORKFasta} ${BadFJDir}${STEM}_BadFJtoReg.sam | awk '{print $4}'`
echo "BadFJ to reg: ${BadFJj3_id}"


BadFJj4_id=`sbatch -J ${STEM}FJ_to_junc --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${SPORKFasta} ${BadFJDir}${STEM}_BadFJtoJunc.sam | awk '{print $4}'`
echo "BadFJ to junc: ${BadFJj4_id}"
depend_str7=${depend_str7}:${BadFJj1_id}:${BadFJj2_id}:${BadFJj3_id}:${BadFJj4_id}

## Read gaps are disallowed in the first version of BadJuncs.  A second version of BadJuncs was created to also find genome/reg/junc/transcriptome alignments with gapped alignments.
## For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.

genomeBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${genomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam"
transcriptomeBOWTIEPARAM="--no-unal --no-mixed  --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${transcriptomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam"
regBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${regIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoReg.sam"
juncBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${juncIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam"


## submit SLURM jobs to do the bowtie alignment for each of the BadFJ Ver 2 indices above
BadFJv2j1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${genomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to genome: ${BadFJv2j1_id}"

BadFJv2j2_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${transcriptomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to transcriptome: ${BadFJv2j2_id}"

BadFJv2j3_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${regBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to reg: ${BadFJv2j3_id}"

BadFJv2j4_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${juncBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to junc: ${BadFJv2j4_id}"

depend_str7=${depend_str7}:${BadFJv2j1_id}:${BadFJv2j2_id}:${BadFJv2j3_id}:${BadFJv2j4_id}

done




# align unaligned files to the FJ bowtie index
#
j8_id=`sbatch -J AlignFJ ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_7AlignFJ.txt -e ${2}err_and_out/err_7AlignFJ.txt --depend=afterok:${j6a_id} ${INSTALLDIR}AlignUnalignedtoFJ.sh ${2} ${ORIG_DIR} | awk '{print $4}'`

echo "align unaligned reads to FJ index: ${j8_id}"
##
#
#
###make FJ naive report
#
j9_id=`sbatch -J FJNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=10:0:0 -o ${2}err_and_out/out_8NaiveRpt.txt -e ${2}err_and_out/err_8NaiveRpt.txt --depend=afterok:${j8_id} ${INSTALLDIR}FarJuncNaiveReport.sh ${2} ${ORIG_DIR} ${5} ${INSTALLDIR} | awk '{print $4}'`

echo "make naive rpt - ${j9_id}"
####

##Make Class input files for GLM
##
###
j15a_id=`sbatch -J AddFJtoIDfile ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=5000 --nodes=1 --time=1:0:0 -o ${2}err_and_out/out_15FJforGLM.txt -e ${2}err_and_out/err_15FJforGLM.txt --depend=afterok:${j9_id} ${INSTALLDIR}parse_FJ_ID_for_GLM.sh ${2} | awk '{print $4}'`

echo "make FJ class input files: ${j15a_id}"

#
##
###
###### ESTIMATE LIGATION ARTIFACT
#
###
###make FarJunctions.fa => Indel.fa files
####
j10_id=`sbatch -J MakeFJIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_10FJIndels.txt -e ${2}err_and_out/err_10FJIndels.txt --depend=afterok:${j6a_id} ${INSTALLDIR}MakeIndelFiles.sh ${2} ${NumIndels} ${INSTALLDIR} | awk '{print $4}'`

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
for (( c=$START; c<=${NumIndels}; c++ ))
do
j11_id=`sbatch -J IndexFJIndel ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=60000 --nodes=6 --time=24:0:0 -o ${2}err_and_out/out_11indexindels.txt -e ${2}err_and_out/err_11indexindels.txt --depend=afterok:${j10_id} ${INSTALLDIR}BowtieIndexFJIndels.sh ${2}FarJuncIndels/ ${c} ${2}BowtieIndels/ ${2} | awk '{print $4}'`
depend_str11=${depend_str11}:${j11_id}

done
echo "index indels ${depend_str11}"



##Align FarJuncSecondary (unaligned to FJ index) to FJ indels


#
##
depend_str13="--depend=afterok"
START=1
for (( c=$START; c<=${NumIndels}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,1 --rdg 50,50 --rfg 50,50"
j13_id=`sbatch -J AlignIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=60000 --nodes=6 --time=24:0:0 -o ${2}err_and_out/out_12alignindels.txt -e ${2}err_and_out/err_12alignindels.txt ${depend_str11} ${INSTALLDIR}BowtieAlignFJIndels.sh ${2} "${BOWTIEPARAMETERS}" ${c} | awk '{print $4}'`
depend_str13=${depend_str13}:${j13_id}

done
echo "align to indels: ${depend_str13}"
#

#
## loop through AlignedIndels directory
## things that don't overlap indels junction by ${5} are removed.
## for every junction, a string is created.  For example, if 3 indel files exist, the string [0, 0, 2, 1, 5, 0, 1] represents 0*3 deletions, 0* 2 deletions, 2 * 1 deletion, 1 * no indels, 5 * 1 insertion, 0 * 2 insertions, 1 * 3 insertions.
## strings are output into IndelsHistogram folder
#
##
j14_id=`sbatch -J FilterIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=5:0:0 -o ${2}err_and_out/out_13filterIndels.txt -e ${2}err_and_out/err_13filterIndels.txt ${depend_str13} ${INSTALLDIR}FindAlignmentArtifact_SLURM.sh ${2} ${5} ${NumIndels} ${INSTALLDIR}| awk '{print $4}'`

echo "make indels histo: ${j14_id}"

##
####
##### REG INDELS ################
##
# Align unaligned files to the expanded reg junctions with indels
AlignedIndels=${1}/orig/RegIndelAlignments/
mkdir -p ${AlignedIndels}

depend_str16="--depend=afterok"
START=1
for (( c=$START; c<=${NumIndels}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --n-ceil L,0,1 --rfg 50,50"
j16_id=`sbatch -J AlignRegIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_15AlignRegIndels.txt -e ${2}err_and_out/err_15AlignRegIndels.txt ${INSTALLDIR}AlignUnalignedtoRegIndel.sh ${1} ${c} ${2} "${BOWTIEPARAMETERS}" ${REG_INDEL_INDICES} | awk '{print $4}'`
depend_str16=${depend_str16}:${j16_id}
done
echo "Aligning unaligned files to linear junc indels: ${depend_str16}"


#
### reg indels class output file
#
j18_id=`sbatch -J RegIndelClass ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_18RegIndelsClassOutput.txt -e ${2}err_and_out/err_18RegIndelsClassOutput.txt ${depend_str16} ${INSTALLDIR}RegIndelsClassID.sh ${2} ${1} ${5} ${INSTALLDIR} | awk '{print $4}'`
echo " Reg Indels Class Output: ${j18_id}"

#
#
#####
####
###
####
######## FJ INDELS CLASS OUTPUT FILES ###########
####
####

#
j19_id=`sbatch -J FJIndelClass ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=60000 --nodes=6 --time=24:0:0 -o ${2}err_and_out/out_19FJIndelsClassOutput.txt -e ${2}err_and_out/err_19FJIndelsClassOutput.txt --depend=afterok:${j14_id} ${INSTALLDIR}FJIndelsClassID.sh ${2} ${1} ${5} ${INSTALLDIR} | awk '{print $4}'`
echo "FJ Indels Class Input: ${j19_id}"


###### RUN GLM ###########################
#
## Run GLM
##
j15b_id=`sbatch -J GLM.r ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=1:0:0 -o ${2}err_and_out/out_15GLM_r.txt -e ${2}err_and_out/err_15GLM_r.txt --depend=afterok:${j15a_id}:${j18_id}:${j19_id} ${INSTALLDIR}run_GLM.sh ${1} ${2} ${INSTALLDIR} | awk '{print $4}'`

echo "Run GLM: ${j15b_id}"


###
####
#
### Append linear junctions GLM report with anomalies, indels
##
j17_id=`sbatch -J AppendRegGLM ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=2 --time=24:0:0 -o ${2}err_and_out/out_17AppendRegGLM.txt -e ${2}err_and_out/err_17AppendGLM.txt ${depend_str16} ${INSTALLDIR}AddIndelstolinearGLM.sh ${1} ${2} ${INSTALLDIR} | awk '{print $4}'`
echo " Appending linearJuncs GLM report: ${j17_id}"
###
####

### Append Naive report:  Add Indels, quality of junctions (good/bad), frequency of junction participation in linear junctions, and GLM report fields to the Naive report
###:${j7_id}:${j9_id}
j15_id=`sbatch -J AppendNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_14AppendRpt.txt -e ${2}err_and_out/err_14AppendRpt.txt ${depend_str7}:${j15b_id} ${INSTALLDIR}AppendNaiveRept.sh ${2} ${GLM_DIR} ${INSTALLDIR} ${2}reports/glmReports/ | awk '{print $4}'`

echo "append naive rpt: ${j15_id}"
##
