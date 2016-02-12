FJDIR=${1}
NumIndels=${2}
ORIG_DIR=${3}

UNALIGNEDDIR=${3}unaligned/
FARJUNCDIR=${1}FarJunctionAlignments/
SECONDFARJUNCDIR=${1}FarJuncSecondary/


START=1
for (( c=$START; c<=${2}; c++ ))
do
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
BOWTIEINDEX=${1}BowtieIndels/Indels${c}
BOWTIEOUTPUT=${1}FarJuncSecondary/AlignedIndels/

for file in ${UNALIGNEDDIR}*.fq
do
FILENAME=$(basename "$file" .fq)
FarJuncSecondaryFile=${SECONDFARJUNCDIR}still_${FILENAME}.fq
OUTPUTFILE=${FILENAME}_indels${c}.sam
j13_id=`sbatch -J AlignIndels ${RESOURCE_FLAG} --mem=55000 --time=24:0:0 -o ${2}err_and_out/out_12alignindels.txt -e ${2}err_and_out/err_12alignindels.txt ${depend_str12} /scratch/PI/horence/gillian/MACHETE/BowtieAligner.batch.sh "${BOWTIEPARAMETERS}" ${BOWTIEINDEX} ${FarJuncSecondaryFile} ${BOWTIEOUTPUT}${OUTPUTFILE} | awk '{print $4}'`
depend_str13=${depend_str13}:${j13_id}
done
done

echo "align indels ${depend_str13}"