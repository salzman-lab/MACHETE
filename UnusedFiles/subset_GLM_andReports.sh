#!/bin/sh

#  All_Samples_GLMandReports.sh
#  
#
#  Created by Gillian Hsieh on 4/23/16.
#
## This shell will loop through the following data sets on Sherlock to regenerate the GLM reports and appended naive reports
## 1. Ewing Sarcoma KD
## 2. Normal Breast
## 3. bladder cancer
## 4. CML_test (CML from caltech)
## 5. prostate cancer
## 6. normal_fetal
## 7. Engstrom simulated data set
## 8. SEQC - normal human and brain sets
## 9. Ov_RNaseR_Qatar
## 10. ovarian data from salzman lab !!!!! -RErun FJ indels class input?
## 11. CML_UConn

OUTPUT_DIR=${1} ## This is the output directory for ALL above samples.  Directories for each sample will be created.
if [ $# -ge 2 ]
then
MODE=${2}
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


INSTALLDIR="/scratch/PI/horence/gillian/MACHETE/"

mkdir -p ${1}



############# NORMAL BREAST ########################


CircPipe_normal_breast=/scratch/PI/horence/gillian/normal_breast/circpipe/
FJDir_normal_breast=/scratch/PI/horence/gillian/normal_breast/FarJunc_Mar9/

Output_normal_breast=${1}normal_breast/
mkdir -p ${Output_normal_breast}
mkdir -p ${Output_normal_breast}err/
mkdir -p ${Output_normal_breast}glmReports/
mkdir -p ${Output_normal_breast}AppendedReports/

NUM_FILES_normal_breast=$((`more ${FJDir_normal_breast}StemList.txt | wc -l`))
echo ${NUM_FILES_normal_breast}

j2_id=`sbatch -J AppendNaiveRpt_normal_breast ${RESOURCE_FLAG} --array=1-${NUM_FILES_normal_breast} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_normal_breast}err/out_AppendNaive.txt -e ${Output_normal_breast}err/err_AppendNaive.txt ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_normal_breast} ${CircPipe_normal_breast}circReads/glmReports/ ${INSTALLDIR} ${Output_normal_breast}glmReports/ ${Output_normal_breast}AppendedReports/ | awk '{print $4}'`

echo "normal_breast Append naive rpt: ${j2_id}"



##################### ov_salzman ################################

CircPipe_ov_salzman=/scratch/PI/horence/alignments/OvarianCancer2014_cutAdapt/
FJDir_ov_salzman=/scratch/PI/horence/gillian/ov_salzman/FarJunc/

Output_ov_salzman=${1}ov_salzman/
mkdir -p ${Output_ov_salzman}
mkdir -p ${Output_ov_salzman}err/
mkdir -p ${Output_ov_salzman}glmReports/
mkdir -p ${Output_ov_salzman}AppendedReports/

NUM_FILES_ov_salzman=$((`more ${FJDir_ov_salzman}StemList.txt | wc -l`))
echo ${NUM_FILES_ov_salzman}


j2_id=`sbatch -J AppendNaiveRpt_ov_salzman ${RESOURCE_FLAG} --array=1-${NUM_FILES_ov_salzman} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_ov_salzman}err/out_AppendNaive.txt -e ${Output_ov_salzman}err/err_AppendNaive.txt ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_ov_salzman} ${CircPipe_ov_salzman}circReads/glmReports/ ${INSTALLDIR} ${Output_ov_salzman}glmReports/ ${Output_ov_salzman}AppendedReports/ | awk '{print $4}'`

echo "ov_salzman Append naive rpt: ${j2_id}"

##################### CML_UConn ################################

CircPipe_CML_UConn=/scratch/PI/horence/gillian/CML_UConn/circpipe_K562/
FJDir_CML_UConn=/scratch/PI/horence/gillian/CML_UConn/FarJunc/

Output_CML_UConn=${1}CML_UConn/
mkdir -p ${Output_CML_UConn}
mkdir -p ${Output_CML_UConn}err/
mkdir -p ${Output_CML_UConn}glmReports/
mkdir -p ${Output_CML_UConn}AppendedReports/

NUM_FILES_CML_UConn=$((`more ${FJDir_CML_UConn}StemList.txt | wc -l`))
echo ${NUM_FILES_CML_UConn}


j2_id=`sbatch -J AppendNaiveRpt_CML_UConn ${RESOURCE_FLAG} --array=1-${NUM_FILES_CML_UConn} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_CML_UConn}err/out_AppendNaive.txt -e ${Output_CML_UConn}err/err_AppendNaive.txt ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_CML_UConn} ${CircPipe_CML_UConn}circReads/glmReports/ ${INSTALLDIR} ${Output_CML_UConn}glmReports/ ${Output_CML_UConn}AppendedReports/ | awk '{print $4}'`

echo "CML_UConn Append naive rpt: ${j2_id}"
