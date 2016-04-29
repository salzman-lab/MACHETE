#!/bin/sh
cd /share/PI/horence/automatedForACClist/
# ./downloadAndGenerateReports.sh /scratch/PI/horence/ACC_FILES/SRR_BCancer_Acc_List.txt /scratch/PI/horence/reads/AutomatedACC /scratch/PI/horence/alignments/AutomatedACC large

ACC_LIST=${1}
READ_PARDIR=${2}
ALIGN_PARDIR=${3}

if [ $# -ge 4 ]
then
  MODE=${4}
else
  MODE=sam_large
fi

# select correct prefix name to use for bowtie index files
if [[ $MODE = *mouse* ]]
then
  bt_prefix="mm10"
elif [[ $MODE = *rat* ]]
then
  bt_prefix="rn5"
elif [[ $MODE = *fly* ]]
then
  bt_prefix="dm3"
elif [[ $MODE = *pombe* ]]
then
  bt_prefix="ASM294v2_23"
elif [[ $MODE = *crypto* ]]
then
  bt_prefix="cryptococcus_neoformans_grubii_h99"
elif [[ $MODE = *cerevisiae* ]]
then
  bt_prefix="Scer"
elif [[ $MODE = *mikatae* ]]
then
  bt_prefix="Smik"
elif [[ $MODE = *bayanus* ]]
then
  bt_prefix="Sban"
elif [[ $MODE = *HSV* ]]
then
  bt_prefix="KOS"
elif [[ $MODE = *zebrafish* ]]
then
  bt_prefix="danRer10"
elif [[ $MODE = *elegan* ]]
then
  bt_prefix="celegans"
elif [[ $MODE = *capsas* ]]
then
  bt_prefix="capsaspora_atcc_30864_2"
elif [[ $MODE = *rosetta* ]]
then
  bt_prefix="salpingoeca_rosetta_1"
else
  bt_prefix="hg19"
fi

# things we don't allow passing as params now
REPORTDIR_NAME=circReads
JUNCTION_MIDPOINT=150
DENOVOCIRC=1
NUM_SAMPLES=1  # since we are treating each sample as its own dataset
NUM_FILES=2  # max number of files possible if this is PE
TRIM_DIR=/share/PI/horence/scripts/trimmingScripts/
CODE_DIR="../circularRNApipeline_SLURM"  # where the analysis pipeline scripts will be run from

# resource allocations 
EXTRACT_MAX_RT="12:0:0"
EXTRACT_VMEM="15000"
WPIF_MAX_RT="1:0:0"
WPIF_VMEM="15000"
WTIF_MAX_RT="1:0:0"
WTIF_VMEM="15000"
TRIM_MAX_RT="12:0:0"
TRIM_VMEM="35000"

# change resource allocations based on job size
if [[ "$MODE" = *large* ]]
then
  JUNC_VMEM="35000"
  GENOME_VMEM="35000"
  TRANSC_VMEM="35000"
  REG_VMEM="35000"
  RIBO_VMEM="25000"
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_VMEM="45000"
  FILTER_VMEM="59000"
  PREPROCESS_MAX_RT="23:0:0"
  FILTER_MAX_RT="23:0:0"
  PFA_VMEM="40000"
elif [[ "$MODE" = *bigmem* ]]
then
  JUNC_VMEM="35000"
  GENOME_VMEM="35000"
  TRANSC_VMEM="35000"
  REG_VMEM="35000"
  RIBO_VMEM="25000"
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_VMEM="55000"
  FILTER_VMEM="150000"  # this is used for the naive reports
  PREPROCESS_MAX_RT="16:0:0"  # this is a limit put on bigmem
  FILTER_MAX_RT="16:0:0"  # this is a limit put on bigmem
  PFA_VMEM="85000"  # this is used for glm reports
else
  JUNC_VMEM="20000"
  GENOME_VMEM="6000"
  TRANSC_VMEM="6000"
  REG_VMEM="20000"
  RIBO_VMEM="5000"
  ALIGN_MAX_RT="12:0:00"
  PREPROCESS_VMEM="15000"
  FILTER_VMEM="15000"
  PREPROCESS_MAX_RT="4:0:0"
  FILTER_MAX_RT="12:0:0"
  PFA_VMEM="40000"
fi

if [[ "$MODE" = *bam* ]]
then
  ALIGN_MAX_RT="23:0:0"
  PREPROCESS_MAX_RT="23:0:0"
  FILTER_MAX_RT="23:0:0"
fi

# we can ask specifically to run on our own machine
if [[ "$MODE" = *horence* ]]
then
  RESOURCE_FLAG="-p horence"
elif [[ "$MODE" = *owners* ]]
then
  RESOURCE_FLAG="-p owners"
fi

# sometimes we need more memory for analysis part only
if [[ "$MODE" = *bigmem* ]]
then
  ANALYSIS_RESOURCE_FLAG="-p bigmem --qos=bigmem"
else
  ANALYSIS_RESOURCE_FLAG=$RESOURCE_FLAG
fi

# start processing samples 1 at a time
while read line    
do
  DATASET_NAME=`echo $line | tr -d '\n\r'`
  
  #### download and extract fastq files into directory under the parent read directory  ####

  RAW_READ_DIR=${READ_PARDIR}/${DATASET_NAME}
  mkdir -p ${RAW_READ_DIR}

  e_id=`sbatch -J Extract${DATASET_NAME} ${RESOURCE_FLAG} --time=${EXTRACT_MAX_RT} --mem=${EXTRACT_VMEM} -D ${RAW_READ_DIR} -o ${DATASET_NAME}Extract_%A_%a.out -e ${DATASET_NAME}Extract_%A_%a.err extractFastq.sh ${DATASET_NAME} | awk '{print $4}'`
  
  #### trim extracted files ####
  TRIMMED_DATASET_NAME=${DATASET_NAME}_cutAdapt  # this is the name that will be used for the trimmed read directory and the alignment output
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${TRIMMED_DATASET_NAME}.txt  # file pipeline will read sample info from
  TRIMMED_READ_DIR=${READ_PARDIR}/${TRIMMED_DATASET_NAME}
  mkdir -p ${TRIMMED_READ_DIR}
  
  # write some intermediate files that help us figure out the parameters for trim galore 
  wpif_id=`sbatch -J WPIF${DATASET_NAME} ${RESOURCE_FLAG} --time=${WPIF_MAX_RT} --mem=${WPIF_VMEM} -D ${TRIM_DIR} -o ${TRIMMED_READ_DIR}/${DATASET_NAME}WPIF_%A.out -e ${TRIMMED_READ_DIR}/${DATASET_NAME}WPIF_%A.err --depend=afterok:${e_id} writePairedIdFiles.sh ${RAW_READ_DIR} ${TRIMMED_READ_DIR} Trim${DATASET_NAME} | awk '{print $4}'`
  
  # now trim
  trim_id=`sbatch -J Trim${DATASET_NAME} ${RESOURCE_FLAG} --time=${TRIM_MAX_RT} --mem=${TRIM_VMEM} -o ${TRIMMED_READ_DIR}/Trim${DATASET_NAME}_%A_%a.out -e ${TRIMMED_READ_DIR}/Trim${DATASET_NAME}_%A_%a.err --depend=afterok:${wpif_id} trimGaloreOnePair.sh ${TRIMMED_READ_DIR}/pairedFileIdList.txt ${TRIMMED_READ_DIR}/metaData.txt ${TRIMMED_READ_DIR} | awk '{print $4}'`
  

done < ${ACC_LIST}
