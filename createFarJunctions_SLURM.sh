#!/bin/sh

## This shell is the entirety of the MACHETE.
## it takes KNIFE input, a reference genome, and some user parameters, and outputs junctions
## greater than "USERBPDIST" below.  In testing we have only used 100Kb because KNIFE looks
## for junctions at less than 100Kb.


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
StemFile=${2}StemList.txt # file containing unique IDs of each experiment e.g. SRR#

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


## remove the old error files
if [ "$(ls -A ${2}/err_and_out/*)" ]
then
echo "old error files exist -- removing"
rm ${2}err_and_out/*
fi

if [ -e ${2}MasterError.txt ]
then
rm ${2}MasterError.txt
fi


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

## if the program has been run before, there will be "sorted" reg and genome files.
## these are removed if they are present.
## All files from the original KNIFE alignments are sorted into alphabetical order because it is faster for python to identify read partners in two alphabetically sorted files than it is to find read pairs in two huge files where read IDs are scrambled.


## the shell AlphabetizeKNIFEreads.sh takes directories reg and genome, where we plan to search for mismatched paired ends, and sorts them alphabetically using the linux sort function

## sorting reg files
j1_id=`sbatch -J sortingReg ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_1sortReg.txt -e ${2}err_and_out/err_1sortReg.txt ${INSTALLDIR}AlphabetizeKNIFEreads.sh ${ORIG_DIR}reg/ ${2} | awk '{print $4}'`
echo "sorting reg files: job# ${j1_id}"

# sorting genome files
j2_id=`sbatch -J sortingGenome ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_1sortGenome.txt -e ${2}err_and_out/err_1sortGenome.txt ${INSTALLDIR}AlphabetizeKNIFEreads.sh ${ORIG_DIR}genome/ ${2} | awk '{print $4}'`
echo "sorting genome files: job# ${j2_id}"


## finding mismatched paired end reads
## The shell PEfinder.sh takes the KNIFE alignment directory (OrigDir -1 ), the output directory (FJDir -2 ), the distance beyond which the user would consider alignments to be "discordant" (BP_distance -3 ), and the MACHETE installation directory (4) and calls a python script PEfinder_genomeAndReg_ENCODE.py.
## The python script identifies paired R1 and R2 from the genome and reg alignments and if the alignment location is > user defined bases apart then records them in an output file within the output directory: FarJunctionDirectory/DistantPEFiles/<STEM>_distant_pairs.txt
## If, for example, a read pair was found to be discordant, and R1= chrA:X, R2=chrB:Y, then the distant_pairs.txt file would contain the readID and chrA:M-N, chrB:P-Q where M-N is a window of 10,000 bases on each side of X and P-Q is a window of 10,000 bases on each side of Y.
## The window of 10,000 bases can be set by the user in the shell PEfinder.sh

j3_id=`sbatch -J PEMismatch ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_2PEfinder.txt -e ${2}err_and_out/err_2PEfinder.txt --depend=afterok:${j1_id}:${j2_id} ${INSTALLDIR}PEfinder.sh ${ORIG_DIR} ${2} ${3} ${INSTALLDIR} | awk '{print $4}'`

echo "Outputting Mismatched paired ends: job ${j3_id} - check for FJDir/DistantPEFiles/*_distant_pairs.txt"


### counting rates of mismatched PE
## Because there are lot of repeat locations in the _distant_pairs.txt file generated above, the DistantPE_Counter.py script is called by the shell of the same name to 1) eliminate duplicate locations.  Another early problem was that the fasta generated later was too enormous to run in a timely fashion so the distant_pairs.txt file is split by this script into 24 smaller files based on the chromosome # of the upstream partner.
## The shell DistantPE_Counter_genome_ENCODE.sh takes in the FarJunction output directory and MACHETE installation directory and outputs <FJDir>/DistantPEFiles/<STEM>/chr1,2,3,4,...,X,Y_Distant_PE_frequency.txt
#  The chrA_Distant_PE_frequency.txt files contain three columns: chrA:M-N, chrB:P-Q, and R, where R is the number of times that these two exact windows were matched together.  R could be used to cull the fasta file if it gets too large, but at this point we are still looking for junctions between exons if only one read pair aligned discordantly.


j4_id=`sbatch -J CountMismatch ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_3PEcounter.txt -e ${2}err_and_out/err_3PEcounter.txt --depend=afterok:${j3_id} ${INSTALLDIR}DistantPE_Counter.sh ${2} ${INSTALLDIR} | awk '{print $4}'`
echo "counting mismatch rates: ${j4_id} - check for FJdir/DistantPEFiles/<STEM>/chr1,2,3,..,X,Y_Distant_PE_frequency.txt"



# sort mismatched PE by chromosome

## This is a simple shell script SortPairedEnds.sh to sort the chrA_Distant_PE_frequency.txt files into alphabetical order.  It takes FJDir/DistantPEFiles/chrA_Distant_PE_frequency.txt and outputs to same directory, sorted_chrA_Distant_PE_frequency.txt using the linux "sort" command.
## The reason for sorting is to increase the speed of the next step.

j5_id=`sbatch -J SortPE ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_4PEsort.txt -e ${2}err_and_out/err_4PEsort.txt --depend=afterok:${j4_id} ${INSTALLDIR}SortPairedEnds.sh ${2} | awk '{print $4}'`

echo "sorting mismatched PE files: ${j5_id} - check for FJDir/DistantPEFiles/sorted_*.txt"

#make Junction fasta file by extracting sequence info from pickles
## Pickle files are binary files used for storage of large amounts of information in Python.  Here, GTF files have been repackaged into pickle files and store sequence, exon name, and exon location information.  Accessing pickles can be time consuming so target locations from the discordant reads have been broken down by chromosome and alphabetized.
## Loop goes through integers 1-24 where that corresponds to chromosome # (23 = X, 24 = Y). For each of those chromosomes, the shell makeJunctions.sh calls makeJunctions.py
## MakeJunctions.py takes in FJDir/DistantPEFiles/sorted__chrA_Distant_PE_frequency.txt, and outputs FJDir/fasta/<STEM>/<STEM>_chrAFarJunctions.fa.  It reads in the discordant windows, and searches the pickle file for all exons names/exon locations/ exon sequences on either the sense or antisense strands within the discordant windows.  Then it makes all pairs of possible exon junctions. All sequences are 300 base pairs long - 150 bases on each side of the breakpoint.  For exons that are fewer than 150 bases, the remainder of the 150 bases is padded with N's.  All pairs include fusions between two exons that are +/+, +/-, -/+, and -/-.  In the case that a sense and antisense exon are fused, then the sequence listed is the exact sequence that would occur if the fusion was read from the 5'->3' direction.  If the generated sequence was BLATted, the correct "strands" would appear in BLAT.  Similarly, for (-)/(-) exon pairs, if the generated sequence was BLATted, the exons would appear on the (-) strand in BLAT.
## THIS IS DIFFERENT IN KNIFE.  In KNIFE, if a (-)/(-) pair was BLATted, it would look as if it were +/+ because the KNIFE reverse complements (-) sequences.  Additionally in KNIFE, there is no way to detect inversions because -/+ and +/- fasta sequences are not generated.


depend_str6="--depend=afterok"
for c in `seq 1 24`; do
A=${c}
    if [ ${c} -eq 23 ]; then A=X; fi
    if [ ${c} -eq 24 ]; then A=Y; fi

    j6_id=`sbatch -J MakeFJFasta ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_5makefasta.txt -e ${2}err_and_out/err_5makefasta.txt --depend=afterok:${j5_id} ${INSTALLDIR}makeJunctions.sh ${PICKLEDIR} ${2} ${A} ${INSTALLDIR} | awk '{print $4}'`
    depend_str6=${depend_str6}:${j6_id}
done

echo "make fusion fasta files: ${depend_str6} - check for FJDir/fasta/<STEM>/<STEM>_chrAFarJunctions.fa"


##
##
##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.
## The script linkfastafiles.sh uses linux to concatenate the <FJDir>/fasta/<STEM>/<STEM>_chr1,2,3,...,X,Y_FarJunctions.fa into a single large fasta <FJDir>/fasta/<STEM>_FarJunctions.fa.
## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>/<STEM>_FJ_Index
#
j6a_id=`sbatch -J FJ_Index ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_5FJIndexing.txt -e ${2}err_and_out/err_5FJIndexing.txt ${depend_str6} ${INSTALLDIR}linkfastafiles.sh ${2} | awk '{print $4}'`
echo "make FJ bowtie indices for each experiment: ${j6a_id} - check for FJDir/BowtieIndex/<STEM>_FJ_Index"

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

FarJuncFasta=${2}fasta/${STEM}*FarJunctions.fa
BadFJDir=${2}BadFJ/${STEM}/
BadFJver2Dir=${2}BadFJ_ver2/${STEM}/
mkdir -p ${BadFJDir}
mkdir -p ${BadFJver2Dir}
r1file=${BadFJver2Dir}${STEM}_FarJunctions_R1.fa
r2file=${BadFJver2Dir}${STEM}_FarJunctions_R2.fa


#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the last 40 non-N reads from the sequence.

j7_id=`sbatch -J BadFJ_Split ${RESOURCE_FLAG} --mem=55000 --nodes=4 --time=10:0:0 -o ${2}err_and_out/out_6BadJunc_split.txt -e ${2}err_and_out/err_6BadJunc_split.txt --depend=afterok:${j6a_id} ${INSTALLDIR}LenientBadFJ_SLURM.sh ${FarJuncFasta} ${BadFJver2Dir} ${2} ${INSTALLDIR} | awk '{print $4}'`
echo "BadfJ ver 2 Split -- ${j7_id}"

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --rdg 50,50 --rfg 50,50"


## submit SLURM jobs to do bowtie alignments for each of BadFJ indices above
if [ -e ${BadFJDir}${STEM}_BadFJtoGenome.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtoGenome.sam exists.  To realign, please manually delete this file first"
else
BadFJj1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoGenome.sam | awk '{print $4}'`
echo "BadFJ to genome: ${BadFJj1_id}"
depend_str7=${depend_str7}:${BadFJj1_id}
fi

if [ -e ${BadFJDir}${STEM}_BadFJtotranscriptome.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtotranscriptome.sam exists.  To realign, please manually delete this file first"
else
BadFJj2_id=`sbatch -J ${STEM}FJ_to_transcriptome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtotranscriptome.sam | awk '{print $4}'`
echo "BadFJ to transcriptome: ${BadFJj2_id}"
depend_str7=${depend_str7}:${BadFJj2_id}
fi

if [ -e ${BadFJDir}${STEM}_BadFJtoReg.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtoReg.sam exists.  To realign, please manually delete this file first"
else
BadFJj3_id=`sbatch -J ${STEM}FJ_to_reg --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoReg.sam | awk '{print $4}'`
echo "BadFJ to reg: ${BadFJj3_id}"
depend_str7=${depend_str7}:${BadFJj3_id}

fi

if [ -e ${BadFJDir}${STEM}_BadFJtoJunc.sam ]
then
echo "${BadFJDir}${STEM}_BadFJtoJunc.sam exists.  To realign, please manually delete this file first"
else
BadFJj4_id=`sbatch -J ${STEM}FJ_to_junc --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJDir}out.txt -e ${BadFJDir}err.txt --depend=afterok:${j6a_id} ${INSTALLDIR}BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${FarJuncFasta} ${BadFJDir}${STEM}_BadFJtoJunc.sam | awk '{print $4}'`
echo "BadFJ to junc: ${BadFJj4_id}"
depend_str7=${depend_str7}:${BadFJj4_id}
fi

## Read gaps are disallowed in the first version of BadJuncs.  A second version of BadJuncs was created to also find genome/reg/junc/transcriptome alignments with gapped alignments.
## For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.

genomeBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${genomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam"
transcriptomeBOWTIEPARAM="--no-unal --no-mixed  --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${transcriptomeIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam"
regBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${regIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoReg.sam"
juncBOWTIEPARAM="--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff -x ${juncIndex} -1 ${r1file} -2 ${r2file} -S ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam"


## submit SLURM jobs to do the bowtie alignment for each of the BadFJ Ver 2 indices above
if [ -e ${BadFJver2Dir}${STEM}_BadFJtoGenome.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoGenome.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j1_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${genomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to genome: ${BadFJv2j1_id}"
depend_str7=${depend_str7}:${BadFJv2j1_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtotranscriptome.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j2_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${transcriptomeBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to transcriptome: ${BadFJv2j2_id}"
depend_str7=${depend_str7}:${BadFJv2j2_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtoReg.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoReg.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j3_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${regBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to reg: ${BadFJv2j3_id}"
depend_str7=${depend_str7}:${BadFJv2j3_id}
fi

if [ -e ${BadFJver2Dir}${STEM}_BadFJtoJunc.sam ]
then
echo "${BadFJver2Dir}${STEM}_BadFJtoJunc.sam exists.  To realign, please manually delete this file first"
else
BadFJv2j4_id=`sbatch -J ${STEM}FJ_to_genome --mem=55000 ${RESOURCE_FLAG} --time=12:0:0 -o ${BadFJver2Dir}out.txt -e ${BadFJver2Dir}err.txt --depend=afterok:${j7_id} ${INSTALLDIR}BowtieAligner_BadFJv2.sh "${juncBOWTIEPARAM}" | awk '{print $4}'`
echo "BadFJ_ver2 to junc: ${BadFJv2j4_id}"
depend_str7=${depend_str7}:${BadFJv2j4_id}
fi
done




# align unaligned files to the FJ bowtie index
# This calls the shell AlignUnalignedtoFJ.  It takes the inputs of the MACHETEoutput directory and the KNIFE unaligned reads (KNIFEdir/orig/unaligned/).  It calls on Bowtie2 to align the unaligned reads for each <STEM> to the Far Junctions bowtie indices located at FJDir/BowtieIndex/<STEM>/<STEM>_FJ_Index.   Bowtie2 parameters include alignment score with mismatch rate of ~4/100 bases, prohibiting read gaps in the reference or given sequence, and N ceiling = read length (e.g. a read consisting of 100% N's would be discarded).  The aligned reads are output to /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam.  Reads that continue to fail to align are output to /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq.
#
j8_id=`sbatch -J AlignFJ ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_7AlignFJ.txt -e ${2}err_and_out/err_7AlignFJ.txt --depend=afterok:${j6a_id} ${INSTALLDIR}AlignUnalignedtoFJ.sh ${2} ${ORIG_DIR} | awk '{print $4}'`

echo "align unaligned reads to FJ index: ${j8_id} - check for /FJDir/FarJunctionAlignments/<STEM>/unaligned_<STEM>_R1/2.sam and /FJDir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_R1/2.fq"
##
#
#
###make FJ naive report
## FarJuncNaiveReport.sh is a shell script that calls the python script FarJuncNaiveReport.py to generate the "Naive Reports".  Inputs include the MACHETE output directory, paths to the KNIFE alignment files, the amount a read should overlap the junction in order to be considered a "true" junctional alignment, and the MACHETE installation directory.
## see the FarJuncNaiveReport.sh for more info on FarJuncNaiveReport.py and details about how alignments are selected as "true" or "false", and how the a p value is calculated.
## The rate of true or anomaly alignments and p values are output to FJDir/reports/<STEM>_naive_report.txt.  Specific read ID's are also tracked and information on them can be found in FJDir/reports/IDs_<STEM>.txt.
#
j9_id=`sbatch -J FJNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_8NaiveRpt.txt -e ${2}err_and_out/err_8NaiveRpt.txt --depend=afterok:${j8_id} ${INSTALLDIR}FarJuncNaiveReport.sh ${2} ${ORIG_DIR} ${5} ${INSTALLDIR} | awk '{print $4}'`

echo "make naive rpt - ${j9_id}"
###
##
##Make fJ Class input files for GLM
##parse_FJ_ID_for_GLM.sh is a simple shell script that takes the ID files generated above in FJDir/reports/IDs_<STEM>.txt as inputs and removes any unmapped or unaligned read IDs using the "grep -v" command.  The outputs are fed into the  FJDir/GLM_classInput/<STEM>__output_FJ.txt with the other class input files for the GLM.
###
j15a_id=`sbatch -J AddFJtoIDfile ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=5000 --nodes=1 --time=1:0:0 -o ${2}err_and_out/out_15FJforGLM.txt -e ${2}err_and_out/err_15FJforGLM.txt --depend=afterok:${j9_id} ${INSTALLDIR}parse_FJ_ID_for_GLM.sh ${2} | awk '{print $4}'`

echo "make FJ class input files: ${j15a_id}"

#
##
###
###### ESTIMATE LIGATION ARTIFACT
## Ligation artifact refers to the rate at which cDNA are artificially ligated, producing what appears to be a true junctional alignment.  In this step, reads that remained unaligned to the Far Junctions Bowtie index that are located in FarJuncSecondary are then aligned to new indel fasta indices under the assumption that in any local area, the rate of false ligation of a junction is similar to the rate of false ligation of any two cDNAs that do not end at exon boundaries.
#
###
### MakeIndelFiles.sh is a shell script that calls the python script AddIndelsToFasta.py.  It takes the FarJunctions fasta files as inputs (FJDir/fasta/<STEM>_FarJunctions.fa) and outputs five files called FJDir/FarJuncIndels/<STEM>/<STEM>_FJ_Indels_1,2,3,4,5.fa where the numbers 1-5 indicate the N's inserted on each side of the breakpoint or deletions on each side of the breakpoint.  For example, the FJ_indels_3 file is the same as the FarJunctions.fa file except every sequence has 3 N's on each side of the breakpoint (total of 6 N's inserted at the breakpoint), or 3 bases are deleted from each exon on each side of the breakpoint (total of 6 deletions at the breakpoint).
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
## The shell script BowtieIndexFJIndels.sh calls bowtie to index the indels_N.fa files that were created in the previous step.  The indels are output to the directory FJDir/BowtieIndels/<STEM>/<STEM>_Indels_N where N is the number of indels in that index.

depend_str11="--depend=afterok"

START=1
for (( c=$START; c<=${NumIndels}; c++ ))
do
j11_id=`sbatch -J IndexFJIndel ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_11indexindels.txt -e ${2}err_and_out/err_11indexindels.txt --depend=afterok:${j10_id} ${INSTALLDIR}BowtieIndexFJIndels.sh ${2}FarJuncIndels/ ${c} ${2}BowtieIndels/ ${2} | awk '{print $4}'`
depend_str11=${depend_str11}:${j11_id}

done
echo "index indels ${depend_str11}"



##Align FarJuncSecondary (unaligned to FJ index) to FJ indels

# This section calls the shell BowtieAlignFJIndels.sh to align the fq files FJdir/FarJuncSecondary/<STEM>/still_unaligned_<STEM>_1/2.fq to the Bowtie2 indices of the far junction indels created in the previous step.  Aligned indels are stored in FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam where N is the number of indels in the Bowtie index where the reads aligned.  The bowtie parameters include a max of ~4 mismatches / 100 basepairs, max #N is the read length, and prohibits gapped alignments or gaps in the read.
#
##
depend_str13="--depend=afterok"
START=1
for (( c=$START; c<=${NumIndels}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,1 -p 4 --rdg 50,50 --rfg 50,50"
j13_id=`sbatch -J AlignIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_12alignindels.txt -e ${2}err_and_out/err_12alignindels.txt ${depend_str11} ${INSTALLDIR}BowtieAlignFJIndels.sh ${2} "${BOWTIEPARAMETERS}" ${c} | awk '{print $4}'`
    depend_str13=${depend_str13}:${j13_id}

done
echo "align to indels: ${depend_str13}"
#

#
## Calls FindAlignmentArtifact_SLURM.sh which is a shell that calls MakeIndelsHisto.py.  The MakeIndelsHisto.py script reads in the aligned indels from the sam files FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/still_unaligned_<STEM>_indels<N>.sam.  It concatenates all the indels_1,2,3,4,5.sam files into a larger file All_<STEM>_1/2_indels.sam in the same directory as the original sam files. In the All_<STEM>_1/2_indels.sam files, any indels that did not overlap the junction by the user specified # base pairs around the breakpoint are removed. Additionally, since the FarJuncSecondary files were aligned independently to the indels_1-5 bowtie indices, the same read could align to multiple indices.  Therefore for each read, the read with the best alignment score is placed in the All_<STEM>_1/2_indels sam file, and all other alignments to other indices are discarded.
#  Then the python script checks the FJdir/FarJunctionAlignments/<STEM>/ sam files and creates an array for each junction.  The array is of length 2*#indels+1. In the case of 5 indels, the length is 11 and represents the number reads that aligned to the junction with [5Del, 4Del, 3Del, 2Del, 1Del, aligned to junction exactly, 1Ins, 2Ins, 3Ins, 4Ins, 5Ins]
## the junction name and this array are output to FJDir/IndelsHistogram/indels_<STEM>_1/2.txt.  These outputs will be used to generate the Appended reports later.

##
j14_id=`sbatch -J FilterIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=5:0:0 -o ${2}err_and_out/out_13filterIndels.txt -e ${2}err_and_out/err_13filterIndels.txt ${depend_str13} ${INSTALLDIR}FindAlignmentArtifact_SLURM.sh ${2} ${5} ${NumIndels} ${INSTALLDIR}| awk '{print $4}'`

echo "make indels histo/FilterIndels: ${j14_id}"

###
####
##### REG INDELS ################
##
#  To train the GLM, indel alignments are also created for the linear junctions.  The reference index of indels to the linear junctions is static and has already been created and is referenced above as "REG_INDEL_INDICES" on line 44.  The script AlignRegIndels calls bowtie to align reads that were unaligned the the KNIFE indices (in KNIFEdir/orig/unaligned/*.fq) to the REG_INDEL_INDICES, with the parameters of 1) approx 4 mismatches / 100 bases, maximum number N's = readlength, and no gapped alignments or read gaps.


AlignedIndels=${1}/orig/RegIndelAlignments/
mkdir -p ${AlignedIndels}

depend_str16="--depend=afterok"
START=1
for (( c=$START; c<=${NumIndels}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,1 -p 4 --rdg 50,50 --rfg 50,50"
j16_id=`sbatch -J AlignRegIndels ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_15AlignRegIndels.txt -e ${2}err_and_out/err_15AlignRegIndels.txt ${INSTALLDIR}AlignUnalignedtoRegIndel.sh ${1} ${c} ${2} "${BOWTIEPARAMETERS}" ${REG_INDEL_INDICES} | awk '{print $4}'`
    depend_str16=${depend_str16}:${j16_id}
done
echo "Aligning unaligned files to linear junc indels: ${depend_str16}"


####
######## MAKE REG AND FJ INDELS CLASS OUTPUT FILES ###########
####
###
### reg indels class output file
#Calls RegIndelsClassID.sh shell which calls RegIndels_ClassIDFile.py.
# Inputs are the regular indel alignments located at KNIFEdir/orig/RegIndelAlignments/<STEM>/unaligned_<STEM>_1/2_indel1,2,3,4,5.sam.  These are concatenated into a single file in the same directory called All_<STEM>_1/2_Regindels.sam.  In the concatenation step, like for the far junctions indels, any reads are omitted if they fail to overlie the junction by the user specified overlap, and also if a read aligns to multiple indel indices, the one with the best alignment score is put into the concatenated file.
## Then, the partner reads of all reads from All_<STEM>_1/2_Regindels.sam are identified and labeled as "good" or "bad".  If a read partner is found in genome, it has priority over transcriptome, which has priority over reg, and finally junc.  A far junction R2 cannot be found in another dictionary, as FJ reads are generated from previously unaligned reads.  If the read partner is in genome, it must be located on the same chromosome, on the opposite reference strand from R1.  If a read partner is in reg, then the downstream exon must be upstream of the uptstream reg indel exon, or the upstream read partner exon must be downstream of the downstream reg indel exon, on the same chromosome. Reference strands must be opposite.    If the partner is in junc, then the reg indel alignment must be located inside the circle created by the scrambled junction, and on the opposite reference strand.  In this manner, class input files are generated for the reg indels, which are located at KNIFE dir/circReads/ids/<STEM>_output_RegIndel.txt
#
j18_id=`sbatch -J RegIndelClass ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=4 --time=24:0:0 -o ${2}err_and_out/out_18RegIndelsClassInput.txt -e ${2}err_and_out/err_18RegIndelsClassInput.txt ${depend_str16} ${INSTALLDIR}RegIndelsClassID.sh ${2} ${1} ${5} ${INSTALLDIR} | awk '{print $4}'`
echo " Reg Indels Class Output: ${j18_id}"

# FJ indels class output file
## This script calls FJIndelsClassID.sh
## This takes the FJdir/FarJunctionSecondary/AlignedIndels/<STEM>/All_<STEM>_1/2_FJindels.sam and identifies read partners.  The same criteria to identify read partners as FarJuncNaiveReport.sh are applied (see above).
## Output files are placed into FJDir/GLM_classInput/<STEM>_output_FJIndels.txt

#
j19_id=`sbatch -J FJIndelClass ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=80000 --nodes=8 --time=24:0:0 -o ${2}err_and_out/out_19FJIndelsClassInput.txt -e ${2}err_and_out/err_19FJIndelsClassInput.txt --depend=afterok:${j14_id} ${INSTALLDIR}FJIndelsClassID.sh ${2} ${1} ${5} ${INSTALLDIR} | awk '{print $4}'`
echo "FJ Indels Class Input: ${j19_id}"


###### RUN GLM ###########################
#
## Run GLM
##  This calls the GLM script.  Class input files from KNIFE dir/circReads/ids/ are fed into the GLM and GLM reports are generated in FJDir/reports/glmReports.  Please see GLM documentation for additional information.

j15b_id=`sbatch -J GLM.r ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=5:0:0 -o ${2}err_and_out/out_15GLM_r.txt -e ${2}err_and_out/err_15GLM_r.txt --depend=afterok:${j15a_id}:${j18_id}:${j19_id} ${INSTALLDIR}run_GLM.sh ${1} ${2} ${INSTALLDIR} | awk '{print $4}'`

echo "Run GLM: ${j15b_id}"



## Append linear junctions GLM report with anomalies, indels
#AddIndelstolinearGLM.sh calls the python script KNIFEglmReportsForMachete.py.  This script parses circular and linear glmReports.  For the linear glmReports from KNIFE, the script collects any junctions where 1) the two exons from the linear report are from different genes or 2) the posterior probability is >0.9.  It adds on the rate of anomaly reads and indels to the reports and feeds them into FJDir/reports/AppendedReports.  For ciruclar reports, the script collects any junctions where the posterior probability is <0.9, appends the "Decoy" rate, and feeds the reports into FJDir/reports/Appended reports.
## The purpose of this script is to place all reports in a single directory for the user.
#
j17_id=`sbatch -J AppendRegGLM ${RESOURCE_FLAG} --array=1-${NUM_FILES}  --mem=55000 --nodes=2 --time=24:0:0 -o ${2}err_and_out/out_17AppendRegGLM.txt -e ${2}err_and_out/err_17AppendGLM.txt ${depend_str16} ${INSTALLDIR}AddIndelstolinearGLM.sh ${1} ${2} ${INSTALLDIR} | awk '{print $4}'`
echo " Appending linearJuncs GLM report: ${j17_id}"
##
####

### The AppendNaiveRept.sh shell calls the AppendNaiveRept.py script.  This reads in the IndelsHistogram, BadFJ and BadFJ_ver2 files, and GLM report results and outputs all the results into a single file in /FJDir/reports/AppendedReports/<STEM>_naive_report_Appended.txt
###
j15_id=`sbatch -J AppendNaiveRpt ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --nodes=4 --time=2:0:0 -o ${2}err_and_out/out_14AppendRpt.txt -e ${2}err_and_out/err_14AppendRpt.txt ${depend_str7}:${j15b_id} ${INSTALLDIR}AppendNaiveRept.sh ${2} ${GLM_DIR} ${INSTALLDIR} ${2}reports/glmReports/ | awk '{print $4}'`

echo "append naive rpt: ${j15_id}"
##
