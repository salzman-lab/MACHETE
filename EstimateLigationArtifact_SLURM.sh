#!/bin/sh

#  EstimateLigationArtifact.sh
#  
#
#  Created by Gillian Hsieh on 10/26/15.
#

BaseDirectory=${1}
NumIndels=${2}
NumBPOverlapAtJunc=${3}


echo "Base Dir ${1}"


FarJuncSecondaryfile=${1}FarJuncSecondary/

mkdir -p ${FarJuncSecondaryfile}AlignedIndels/
mkdir -p ${1}BowtieIndels/
mkdir -p ${1}FarJuncIndels/
mkdir -p ${FarJuncSecondaryfile}AlignedIndels/RemoveNonOverlap/
mkdir -p ${1}IndelsHistogram/
mkdir -p ${1}reports/withIndels/


echo "Making FarJunctions.fa into Indel.fa files"

python /scratch/PI/horence/gillian/createFarJunctionsIndex/AddIndelsToFasta.py -o ${1} -n ${2}

echo "Making Indel.fa files into indices"

mv ${1}*FarJunctions_indels_*.fa ${1}FarJuncIndels/

depend_str1="--depend=afterok"

for file in ${1}FarJuncIndels/FarJunctions_indels_*.fa
do
    echo ${file}
    depend_str1=${depend_str1}:
    FILENAME=$(basename "$file" .fa)
    INDEXNAME=${1}BowtieIndels/${FILENAME}
    j_id=`sbatch -J Indexing --mem=55000 --time=24:0:0 -o out.txt -e err.txt /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieIndexer.batch.sh ${file} ${1}BowtieIndels/${FILENAME} | awk '{print $4}'`
    depend_str1=${depend_str1}${j_id}

done

echo ${depend_str1}

## FarJuncSecondary is the directory containing unaligned reads that did not align to the far junc index.
## This step aligns all files in FarJuncsecondary to each of the Indels1, Indels2, Indels3, etc indices
## Indel alignments output to FarJuncSecondary/AlignedIndels

depend_str2="--depend=afterok"
echo "Aligning reads to indels indices"

START=1
for (( c=$START; c<=${2}; c++ ))
do
    BOWTIEINDEX=${1}BowtieIndels/Indels${c}
    BOWTIEOUTPUT=${1}FarJuncSecondary/AlignedIndels/
    for file in ${1}FarJuncSecondary/*.fq
    do
        FILENAME=$(basename "$file" .fq)
        OUTPUTFILE=${FILENAME}_indels${c}.sam
        j_id=`sbatch -J AlignIndels --mem=55000 --time=24:0:0 -o out.txt -e err.txt ${depend_str1} /scratch/PI/horence/gillian/createFarJunctionsIndex/BowtieAligner.batch.sh ${BOWTIEINDEX} ${file} ${BOWTIEOUTPUT}${OUTPUTFILE} | awk '{print $4}'`
        depend_str2=${depend_str2}:${j_id}

    done

done

echo ${depend_str2}


## This section loops through AlignedIndels directory.
## alignments that don't overlap the junction by at least ${3} (user set BP len) are discarded.
## for every junction, a string is created.  For example, if 3 indel files exist, the string [0, 0, 2, 1, 5, 0, 1] represents 0*3 deletions, 0* 2 deletions, 2 * 1 deletion, 1 * no indels, 5 * 1 insertion, 0 * 2 insertions, 1 * 3 insertions.
## strings are output into IndelsHistogram folder

echo "Creating indel distributions in IndelsHistogram folder"


j_id=`sbatch -J CountIndels --mem=55000 --time=24:0:0 -o out.txt -e err.txt ${depend_str2} /scratch/PI/horence/gillian/createFarJunctionsIndex/FindAlignmentArtifact.sh ${1} ${3} ${2}| awk '{print $4}'`

depend_str3=--depend=afterok:${j_id}

##  the IndelsHistogram strings are added to the naive reports.

echo "Adding Indel distributions to Naive reports"

j_id=`sbatch -J AddIndelToRept --mem=55000 --time=24:0:0 -o out.txt -e err.txt ${depend_str3} /scratch/PI/horence/gillian/createFarJunctionsIndex/AddIndelsToNaiveRept.sh ${1}| awk '{print $4}'`

depend_str4=--depend=afterok:${j_id}


                                                                                           