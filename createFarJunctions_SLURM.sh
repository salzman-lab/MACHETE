#!/bin/sh

# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for


ORIG_DIR=${1}
OUTPUT_DIR=${2}
USERBPDIST=${3}
PICKLEFILE=${4}
NUMBASESAROUNDJUNC=${5} #default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70

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

DISTANTPEFILE=${2}Distant_PE_frequency.txt
DISTANTPESORTED=${2}Distant_PE_frequency_sorted.txt
FASTAFILE=${2}FarJunctions.fa
BOWTIE_DIR=${2}BowtieIndex/
BOWTIEINDEX=${2}BowtieIndex/Index
UNALIGNEDDIR=${1}unaligned/
FARJUNCDIR=${2}FarJunctionAlignments/
SECONDFARJUNCDIR=${2}FarJuncSecondary/
FarJuncFasta=${2}FarJunctions.fa
BadFJDir=${2}BadFJ/
genomeIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_genome
transcriptomeIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_transcriptome
regIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_junctions_reg
juncIndex=/share/PI/horence/circularRNApipeline_SLURM/index/hg19_junctions_scrambled

mkdir -p ${2}reports/
mkdir -p ${2}err_and_out/
mkdir -p ${BOWTIE_DIR}
mkdir -p ${FARJUNCDIR}
mkdir -p ${SECONDFARJUNCDIR}
mkdir -p ${BadFJDir}

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
NUM_FILES=$((`ls ${1}genome/*.sam | wc -l`/2))
echo ${NUM_FILES}

echo "${1}"
echo "${2}"
echo "${3}"
echo "${4}"
echo "${5}"

rm ${1}reg/sorted*.sam
rm ${1}genome/sorted*.sam

# sorting reg and genome files


depend_str="--depend=afterok"

for file in ${1}reg/*
do
#echo ${file}
FILENAME=$(basename $file)
if ! [[ "$FILENAME" = *sorted* ]]
    then
    SORTEDFILE=${1}reg/sorted_${FILENAME}
    j1_id=`sbatch -J sortingReg${FILENAME} ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_1sort1${FILENAME}.txt -e ${2}err_and_out/err_1sort1${FILENAME}.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/AlphabetizeENCODEreads.sh ${file} ${SORTEDFILE} | awk '{print $4}'`
   depend_str=${depend_str}:${j1_id}
fi
done;


for file in ${1}genome/*
do
#echo ${file}
FILENAME=$(basename $file)
if ! [[ "$FILENAME" = *sorted* ]]
then
    SORTEDFILE=${1}genome/sorted_${FILENAME}
#   echo ${SORTEDFILE}
    j2_id=`sbatch -J sortingGenome${FILENAME} ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_1sort2${FILENAME}.txt -e ${2}err_and_out/err_1sort2${FILENAME}.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/AlphabetizeENCODEreads.sh ${file} ${SORTEDFILE} | awk '{print $4}'`
    depend_str=${depend_str}:${j2_id}

fi
done;

echo ${depend_str}

# write PE match ID files.
echo "writing PE match ID file"

j3_id=`sbatch -J MakeFileList ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_2PEfinder.txt -e ${2}err_and_out/err_2PEfinder.txt ${depend_str} /scratch/PI/horence/gillian/createFarJunctionsIndex/GetGenomeRegFilenames.sh ${1} ${2} ${3} | awk '{print $4}'`

# finding mismatched paired end reads

j3a_id=`sbatch -J PEMismatch ${RESOURCE_FLAG} --array=1-${NUM_FILES} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_2PEfinder.txt -e ${2}err_and_out/err_2PEfinder.txt --depend=afterok:${j3_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/PEfinder_genomeAndReg_ENCODE.sh ${1} ${2} ${3} | awk '{print $4}'`

 #window around position in genome is 10,000


# counting rates of mismatched PE
j4_id=`sbatch -J CountMismatch ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_3PEcounter.txt -e ${2}err_and_out/err_3PEcounter.txt --depend=afterok:${j3a_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/DistantPE_Counter_genome_ENCODE.sh ${2} | awk '{print $4}'`

# sort mismatched PE by chromosome
j5_id=`sbatch -J SortPE ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_4PEsort.txt -e ${2}err_and_out/err_4PEsort.txt --depend=afterok:${j4_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/SortPairedEnds.sh ${DISTANTPEFILE} ${DISTANTPESORTED} | awk '{print $4}'`

# make Junction fasta file by extracting sequence info from pickles
j6_id=`sbatch -J MakeFJFasta ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_5makefasta.txt -e ${2}err_and_out/err_5makefasta.txt --depend=afterok:${j5_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/makeJunctions.sh ${4} ${2} | awk '{print $4}'`


# make FJ bowtie index from FJ Fasta file
j7_id=`sbatch -J FJ_Indexing ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_6makeIndex.txt -e ${2}err_and_out/err_6makeIndex.txt --depend=afterok:${j6_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieIndexer.batch.sh ${FASTAFILE} ${BOWTIEINDEX} | awk '{print $4}'`

# make BadJunc directory --  bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align

BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --np 0 --rdg 50,50 --rfg 50,50"

BadFJj1_id=`sbatch -J FJ.fa_to_genome --mem=55000 ${RESOURCE_FLAG} --time=24:0:0 -o ${2}err_and_out/out_6BadFJ_genome.txt -e ${2}err_and_out/err_6BadFJ_genome.txt --depend=afterok:${j7_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${genomeIndex} ${FarJuncFasta} ${BadFJDir}BadFJtoGenome.sam | awk '{print $4}'`

BadFJj2_id=`sbatch -J FJ.fa_to_transcriptome --mem=55000 ${RESOURCE_FLAG} --time=24:0:0 -o ${2}err_and_out/out_6BadFJ_transcriptome.txt -e ${2}err_and_out/err_6BadFJ_transcriptome.txt --depend=afterok:${j7_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${transcriptomeIndex} ${FarJuncFasta} ${BadFJDir}BadFJtotranscriptome.sam | awk '{print $4}'`

BadFJj3_id=`sbatch -J FJ.fa_to_reg --mem=55000 ${RESOURCE_FLAG} --time=24:0:0 -o ${2}err_and_out/out_6BadFJ_reg.txt -e ${2}err_and_out/err_6BadFJ_reg.txt --depend=afterok:${j7_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${regIndex} ${FarJuncFasta} ${BadFJDir}BadFJtoReg.sam | awk '{print $4}'`

BadFJj4_id=`sbatch -J FJ.fa_to_genome --mem=55000 ${RESOURCE_FLAG} --time=24:0:0 -o ${2}err_and_out/out_6BadFJ_junc.txt -e ${2}err_and_out/err_6BadFJ_junc.txt --depend=afterok:${j7_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${juncIndex} ${FarJuncFasta} ${BadFJDir}BadFJtoJunc.sam | awk '{print $4}'`

depend_str_BadFJ=--depend=afterok:${BadFJj1_id}:${BadFJj2_id}:${BadFJj3_id}:${BadFJj4_id}


# align unaligned files to the FJ bowtie index

depend_str="--depend=afterok"

for file in ${UNALIGNEDDIR}*.fq
do
    depend_str=${depend_str}:
    FILENAME=$(basename "$file" .fq)
    echo "creating ${FILENAME}"
    # if i wanted to feed unaligned files into a SecondFarJuncDir --un ${SECONDFARJUNCDIR}still_${FILENAME}.fq tag would be applied
    BOWTIEPARAM="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50 --un ${SECONDFARJUNCDIR}still_${FILENAME}.fq"

    j8_id=`sbatch -J Align${FILENAME} ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_7Align${FILENAME}.txt -e ${2}err_and_out/err_7Align${FILENAME}.txt --depend=afterok:${j7_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh "${BOWTIEPARAM}" ${2}BowtieIndex/Index ${file} ${FARJUNCDIR}${FILENAME}.sam | awk '{print $4}'`
    depend_str=${depend_str}${j8_id}
done;

echo ${depend_str}

# make FJ naive report

j9_id=`sbatch -J FJNaiveRpt ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_8NaiveRpt.txt -e ${2}err_and_out/err_8NaiveRpt.txt ${depend_str} /scratch/PI/horence/gillian/createFarJunctionsIndex/FarJuncNaiveReport.sh ${2} ${1} ${5}| awk '{print $4}'`

# clean up old files

j10_id=`sbatch -J CleanUp ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_9clean.txt -e ${2}err_and_out/err_9clean.txt --depend=afterok:${j9_id} /scratch/PI/horence/gillian/createFarJunctionsIndex/cleanup.sh ${1}genome/sorted* ${1}reg/sorted* ${DISTANTPEFILE} ${2}FarJunctions_duplicates.fa ${2}*distant_pairs.txt ${OUTPUT_DIR}distant_pairs/ ${BOWTIEDIR}*.txt| awk '{print $4}'`





