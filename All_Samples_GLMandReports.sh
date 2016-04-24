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



#######EWING SARCOMA KD###################
CircPipe_Ewing=/scratch/PI/horence/gillian/Ewing/circpipe/
FJDir_Ewing=/scratch/PI/horence/gillian/Ewing/FarJunc/

Output_Ewing=${1}EwingKD/
mkdir -p ${Output_Ewing}
mkdir -p ${Output_Ewing}err/
mkdir -p ${Output_Ewing}glmReports/
mkdir -p ${Output_Ewing}AppendedReports/

NUM_FILES_Ewing=$((`more ${FJDir_Ewing}StemList.txt | wc -l`))
echo ${NUM_FILES_Ewing}

j1_id=`sbatch -J GLM_ewing.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_Ewing} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_Ewing}err/out_GLM_r.txt -e ${Output_Ewing}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_Ewing} ${FJDir_Ewing} ${INSTALLDIR} ${Output_Ewing}glmReports/ | awk '{print $4}'`
echo "Ewing GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_ewing ${RESOURCE_FLAG} --array=1-${NUM_FILES_Ewing} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_Ewing}err/out_AppendNaive.txt -e ${Output_Ewing}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_Ewing} ${CircPipe_Ewing}circReads/glmReports/ ${INSTALLDIR} ${Output_Ewing}glmReports/ ${Output_Ewing}AppendedReports/ | awk '{print $4}'`

echo "Ewing Append naive rpt: ${j2_id}"
#
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

j1_id=`sbatch -J GLM_normal_breast.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_normal_breast} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_normal_breast}err/out_GLM_r.txt -e ${Output_normal_breast}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_normal_breast} ${FJDir_normal_breast} ${INSTALLDIR} ${Output_normal_breast}glmReports/ | awk '{print $4}'`
echo "normal_breast GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_normal_breast ${RESOURCE_FLAG} --array=1-${NUM_FILES_normal_breast} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_normal_breast}err/out_AppendNaive.txt -e ${Output_normal_breast}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_normal_breast} ${CircPipe_normal_breast}circReads/glmReports/ ${INSTALLDIR} ${Output_normal_breast}glmReports/ ${Output_normal_breast}AppendedReports/ | awk '{print $4}'`

echo "normal_breast Append naive rpt: ${j2_id}"


################### BLADDER CANCER #########################

CircPipe_bladder=/scratch/PI/horence/gillian/bladder/circpipe/
FJDir_bladder=/scratch/PI/horence/gillian/bladder/FarJunc/

Output_bladder=${1}bladder/
mkdir -p ${Output_bladder}
mkdir -p ${Output_bladder}err/
mkdir -p ${Output_bladder}glmReports/
mkdir -p ${Output_bladder}AppendedReports/

NUM_FILES_bladder=$((`more ${FJDir_bladder}StemList.txt | wc -l`))
echo ${NUM_FILES_bladder}

j1_id=`sbatch -J GLM_bladder.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_bladder} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_bladder}err/out_GLM_r.txt -e ${Output_bladder}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_bladder} ${FJDir_bladder} ${INSTALLDIR} ${Output_bladder}glmReports/ | awk '{print $4}'`
echo "bladder GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_bladder ${RESOURCE_FLAG} --array=1-${NUM_FILES_bladder} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_bladder}err/out_AppendNaive.txt -e ${Output_bladder}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_bladder} ${CircPipe_bladder}circReads/glmReports/ ${INSTALLDIR} ${Output_bladder}glmReports/ ${Output_bladder}AppendedReports/ | awk '{print $4}'`

echo "bladder Append naive rpt: ${j2_id}"


################# CML_TEST (CALTECH) ###########################

CircPipe_CML_test=/scratch/PI/horence/gillian/CML_test/aligned/CML/
FJDir_CML_test=/scratch/PI/horence/gillian/CML_test/FarJunc_Feb22/

Output_CML_test=${1}CML_test/
mkdir -p ${Output_CML_test}
mkdir -p ${Output_CML_test}err/
mkdir -p ${Output_CML_test}glmReports/
mkdir -p ${Output_CML_test}AppendedReports/

NUM_FILES_CML_test=$((`more ${FJDir_CML_test}StemList.txt | wc -l`))
echo ${NUM_FILES_CML_test}

j1_id=`sbatch -J GLM_CML_test.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_CML_test} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_CML_test}err/out_GLM_r.txt -e ${Output_CML_test}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_CML_test} ${FJDir_CML_test} ${INSTALLDIR} ${Output_CML_test}glmReports/ | awk '{print $4}'`
echo "CML_test GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_CML_test ${RESOURCE_FLAG} --array=1-${NUM_FILES_CML_test} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_CML_test}err/out_AppendNaive.txt -e ${Output_CML_test}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_CML_test} ${CircPipe_CML_test}circReads/glmReports/ ${INSTALLDIR} ${Output_CML_test}glmReports/ ${Output_CML_test}AppendedReports/ | awk '{print $4}'`

echo "CML_test Append naive rpt: ${j2_id}"

################### prostate ########################################

CircPipe_prostate=/scratch/PI/horence/gillian/prostate/aligned/circpipe/
FJDir_prostate=/scratch/PI/horence/gillian/prostate/aligned/FarJunc_Feb11/

Output_prostate=${1}prostate/
mkdir -p ${Output_prostate}
mkdir -p ${Output_prostate}err/
mkdir -p ${Output_prostate}glmReports/
mkdir -p ${Output_prostate}AppendedReports/

NUM_FILES_prostate=$((`more ${FJDir_prostate}StemList.txt | wc -l`))
echo ${NUM_FILES_prostate}

j1_id=`sbatch -J GLM_prostate.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_prostate} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_prostate}err/out_GLM_r.txt -e ${Output_prostate}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_prostate} ${FJDir_prostate} ${INSTALLDIR} ${Output_prostate}glmReports/ | awk '{print $4}'`
echo "prostate GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_prostate ${RESOURCE_FLAG} --array=1-${NUM_FILES_prostate} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_prostate}err/out_AppendNaive.txt -e ${Output_prostate}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_prostate} ${CircPipe_prostate}circReads/glmReports/ ${INSTALLDIR} ${Output_prostate}glmReports/ ${Output_prostate}AppendedReports/ | awk '{print $4}'`

echo "prostate Append naive rpt: ${j2_id}"

############################# normal_fetal #####################

CircPipe_normal_fetal=/scratch/PI/horence/gillian/normal_fetal/fetalcircpipe/
FJDir_normal_fetal=/scratch/PI/horence/gillian/normal_fetal/FarJunc/

Output_normal_fetal=${1}normal_fetal/
mkdir -p ${Output_normal_fetal}
mkdir -p ${Output_normal_fetal}err/
mkdir -p ${Output_normal_fetal}glmReports/
mkdir -p ${Output_normal_fetal}AppendedReports/

NUM_FILES_normal_fetal=$((`more ${FJDir_normal_fetal}StemList.txt | wc -l`))
echo ${NUM_FILES_normal_fetal}

j1_id=`sbatch -J GLM_normal_fetal.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_normal_fetal} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_normal_fetal}err/out_GLM_r.txt -e ${Output_normal_fetal}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_normal_fetal} ${FJDir_normal_fetal} ${INSTALLDIR} ${Output_normal_fetal}glmReports/ | awk '{print $4}'`
echo "normal_fetal GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_normal_fetal ${RESOURCE_FLAG} --array=1-${NUM_FILES_normal_fetal} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_normal_fetal}err/out_AppendNaive.txt -e ${Output_normal_fetal}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_normal_fetal} ${CircPipe_normal_fetal}circReads/glmReports/ ${INSTALLDIR} ${Output_normal_fetal}glmReports/ ${Output_normal_fetal}AppendedReports/ | awk '{print $4}'`

echo "normal_fetal Append naive rpt: ${j2_id}"

####################### Engstrom ###############################

CircPipe_Engstrom=/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom/
FJDir_Engstrom=/scratch/PI/horence/gillian/Engstrom/FarJunc/

Output_Engstrom=${1}Engstrom/
mkdir -p ${Output_Engstrom}
mkdir -p ${Output_Engstrom}err/
mkdir -p ${Output_Engstrom}glmReports/
mkdir -p ${Output_Engstrom}AppendedReports/

NUM_FILES_Engstrom=$((`more ${FJDir_Engstrom}StemList.txt | wc -l`))
echo ${NUM_FILES_Engstrom}

j1_id=`sbatch -J GLM_Engstrom.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_Engstrom} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_Engstrom}err/out_GLM_r.txt -e ${Output_Engstrom}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_Engstrom} ${FJDir_Engstrom} ${INSTALLDIR} ${Output_Engstrom}glmReports/ | awk '{print $4}'`
echo "Engstrom GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_Engstrom ${RESOURCE_FLAG} --array=1-${NUM_FILES_Engstrom} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_Engstrom}err/out_AppendNaive.txt -e ${Output_Engstrom}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_Engstrom} ${CircPipe_Engstrom}circReads/glmReports/ ${INSTALLDIR} ${Output_Engstrom}glmReports/ ${Output_Engstrom}AppendedReports/ | awk '{print $4}'`

echo "Engstrom Append naive rpt: ${j2_id}"


###################### SEQ-C ####################################

CircPipe_SEQC=/scratch/PI/horence/gillian/SEQC_study_set/circpipe_SEQC/
FJDir_SEQC=/scratch/PI/horence/gillian/SEQC_study_set/FarJunc/

Output_SEQC=${1}SEQC/
mkdir -p ${Output_SEQC}
mkdir -p ${Output_SEQC}err/
mkdir -p ${Output_SEQC}glmReports/
mkdir -p ${Output_SEQC}AppendedReports/

NUM_FILES_SEQC=$((`more ${FJDir_SEQC}StemList.txt | wc -l`))
echo ${NUM_FILES_SEQC}

j1_id=`sbatch -J GLM_SEQC.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_SEQC} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_SEQC}err/out_GLM_r.txt -e ${Output_SEQC}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_SEQC} ${FJDir_SEQC} ${INSTALLDIR} ${Output_SEQC}glmReports/ | awk '{print $4}'`
echo "SEQC GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_SEQC ${RESOURCE_FLAG} --array=1-${NUM_FILES_SEQC} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_SEQC}err/out_AppendNaive.txt -e ${Output_SEQC}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_SEQC} ${CircPipe_SEQC}circReads/glmReports/ ${INSTALLDIR} ${Output_SEQC}glmReports/ ${Output_SEQC}AppendedReports/ | awk '{print $4}'`

echo "SEQC Append naive rpt: ${j2_id}"


#################### Ov_RNaseR_Qatar ############################

CircPipe_ov_RNase_Qatar=/scratch/PI/horence/gillian/ov_RNaseR_Qatar/ovCircPipe/
FJDir_ov_RNase_Qatar=/scratch/PI/horence/gillian/ov_RNaseR_Qatar/FarJunc/

Output_ov_RNase_Qatar=${1}ov_RNase_Qatar/
mkdir -p ${Output_ov_RNase_Qatar}
mkdir -p ${Output_ov_RNase_Qatar}err/
mkdir -p ${Output_ov_RNase_Qatar}glmReports/
mkdir -p ${Output_ov_RNase_Qatar}AppendedReports/

NUM_FILES_ov_RNase_Qatar=$((`more ${FJDir_ov_RNase_Qatar}StemList.txt | wc -l`))
echo ${NUM_FILES_ov_RNase_Qatar}

j1_id=`sbatch -J GLM_ov_RNase_Qatar.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_ov_RNase_Qatar} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_ov_RNase_Qatar}err/out_GLM_r.txt -e ${Output_ov_RNase_Qatar}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_ov_RNase_Qatar} ${FJDir_ov_RNase_Qatar} ${INSTALLDIR} ${Output_ov_RNase_Qatar}glmReports/ | awk '{print $4}'`
echo "ov_RNase_Qatar GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_ov_RNase_Qatar ${RESOURCE_FLAG} --array=1-${NUM_FILES_ov_RNase_Qatar} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_ov_RNase_Qatar}err/out_AppendNaive.txt -e ${Output_ov_RNase_Qatar}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_ov_RNase_Qatar} ${CircPipe_ov_RNase_Qatar}circReads/glmReports/ ${INSTALLDIR} ${Output_ov_RNase_Qatar}glmReports/ ${Output_ov_RNase_Qatar}AppendedReports/ | awk '{print $4}'`

echo "ov_RNase_Qatar Append naive rpt: ${j2_id}"

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

j1_id=`sbatch -J GLM_ov_salzman.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_ov_salzman} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_ov_salzman}err/out_GLM_r.txt -e ${Output_ov_salzman}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_ov_salzman} ${FJDir_ov_salzman} ${INSTALLDIR} ${Output_ov_salzman}glmReports/ | awk '{print $4}'`
echo "ov_salzman GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_ov_salzman ${RESOURCE_FLAG} --array=1-${NUM_FILES_ov_salzman} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_ov_salzman}err/out_AppendNaive.txt -e ${Output_ov_salzman}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_ov_salzman} ${CircPipe_ov_salzman}circReads/glmReports/ ${INSTALLDIR} ${Output_ov_salzman}glmReports/ ${Output_ov_salzman}AppendedReports/ | awk '{print $4}'`

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

j1_id=`sbatch -J GLM_CML_UConn.r ${RESOURCE_FLAG} --array=1-${NUM_FILES_CML_UConn} --mem=55000 --nodes=4 --time=5:0:0 -o ${Output_CML_UConn}err/out_GLM_r.txt -e ${Output_CML_UConn}err/out_GLM_r.txt ${INSTALLDIR}run_GLM.sh ${CircPipe_CML_UConn} ${FJDir_CML_UConn} ${INSTALLDIR} ${Output_CML_UConn}glmReports/ | awk '{print $4}'`
echo "CML_UConn GLM: ${j1_id}"


j2_id=`sbatch -J AppendNaiveRpt_CML_UConn ${RESOURCE_FLAG} --array=1-${NUM_FILES_CML_UConn} --mem=55000 --nodes=4 --time=2:0:0 -o ${Output_CML_UConn}err/out_AppendNaive.txt -e ${Output_CML_UConn}err/out_AppendNaive.txt --depend=afterok:${j1_id} ${INSTALLDIR}AppendNaiveRept.sh ${FJDir_CML_UConn} ${CircPipe_CML_UConn}circReads/glmReports/ ${INSTALLDIR} ${Output_CML_UConn}glmReports/ ${Output_CML_UConn}AppendedReports/ | awk '{print $4}'`

echo "CML_UConn Append naive rpt: ${j2_id}"
