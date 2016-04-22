#!/bin/sh

#  FarJuncNaiveReport.sh
#  
#
#  Created by Gillian Hsieh on 1/8/16.
#

## FarJuncNaiveReport.sh is a shell script that calls the python script FarJuncNaiveReport.py to generate the "Naive Reports".
## files are output to FJDir/reports/<STEM>_naive_report.txt and FJDir/reports/IDs_<STEM>.txt.



FJDir=${1} ## MACHETE output dir
OrigDir=${2} ## KNIFE alignment files
Window=${3} ## Num bases that read must overlap junction to be considered
INSTALLDIR=${4} ## MACHETE installation directory

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`

## FarJuncNaiveReport.py sequentially opens the far junctions read 1 and 2 alignments, genome read 1 and 2 alignments, reg (linear junction) 1 and 2 alignments, junc (scrambled junction) 1 and 2 alignments, and the unaligned files.
## Alignments are designated as "true" or "false".  All true alignments are used to calculate naive p values using the read length and the alignment score.  Using an estimated mismatch rate of 0.01 (per Illumina specifications), the mismatch rate is compared to the estimated mismatch rate using a Poisson cdf.
## First FJ_1 is opened and compared to FJ2 to look for read partners.  If a read pair R1 and R2 in these files align to the exact same junction, those are considered a "true" FJ alignment.  If a read pair R1 and R2 do not align to the same junction, these are a "false" alignment.
## If read partners R1 or R2 are not found in the FJ Directory, then the genome alignment files are opened.  If R1 is in FJ and R2 in genome, then the following criteria are used to designate those partners as "true" -- R2 must occur on the same chromosome as either upstream of the 5' gene or downstream of the 3' gene, and the two reads must have aligned to opposite reference strands.  If these criteria are not all met, then the read is considered "false". The same is done for FJ R2 and genome R1.
## If the read partner is not in genome, the next library opened is the linear (reg) alignments.  For R1 in FJ and R2 in reg, the following criteria are used to designate read partners as "true".  For the location criterion, the downstream exon of the reg R2 must occur on the same chromosome and upstream of the upstream exon of the FJ R1 or the upstream exon of the reg R2 must occur on the same chromosome and downstream of the downstream exon of FJ R1.  Whichever of the R1 and R2 exons that meet the location criterion must also be on the same strand (+ or -).  Because the KNIFE sequences for (-) stranded exons are reverse complemented whereas the MACHETE sequences are not reverse complemented, for (-) R1 and R2 exons the reference strands must match (0 and 0 or 16 and 16, indicating that alignments were made to opposite sides of the cDNA), but for (+) R1 and R2 exons, the reference strands must be opposite (0 and 16). The same is done for FJ R2 and reg R1.
## if the read partner R2 is not in FJ, genome, or reg, the partner is searched for in scrambled junction alignments.  For R1 in FJ and R2 in junc, all reads are considered "false".  The same is done for FJ R2 and junc R1.
## if a FJ R1's read partner R2 is in none of the above, the unaligned reads are searched for R2. The same is done for FJ R2 and unaligned R1.
## all reads where R1 or R2 was in FJ and the partner is in none of the above are tagged as "unmapped".  One caveat - if a read was found in any of the junctional alignment files (FJ, reg, or junc) but did not overlap the junction by the below specified overlap (-w ${3}) then that read is not considered.  For example, for FJ R1 that did overlap the junction by the desired number of bases and genome R2 that did NOT overlap the junction by the desired number of bases, the script would continue to search for R2 in reg, junc, and unaligned.  If the R2 does not occur in any of those other alignment files, then R2 would be labeled "unmapped" even though it technically did map to a genome alignment.
## output files are in FJDir/reports/<STEM>_naive_report.txt for P values and total numbers of true and false reads (Designated by genome, genome anomaly, reg, reg anomaly, etc), and individual read IDs and partners and their information are output into FJDir/reports/IDs_<STEM>.txt

ml load python/2.7.5
python ${INSTALLDIR}FarJuncNaiveReport.py -s ${STEM} -f ${1} -i ${2} -w ${3}

echo "FarJuncNaiveRept.sh complete for ${STEM} - check for ${1}reports/${STEM}_naive_report.txt and ${1}reports/IDs_${STEM}.txt" >> ${1}MasterError.txt
