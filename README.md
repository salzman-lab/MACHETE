# MACHETE
More Accurate vs Mismatched Alignment CHimEra Tracking Engine

MACHETE is a fusion detection software that is in development.

REQUIREMENTS: 

1. KNIFE - https://github.com/lindaszabo/KNIFE
   NOTE --  before running KNIFE, please make sure your fastq files are in the following format:

"UNIQUE STRING"_"read1 or 2 file notation".fq

some examples of acceptable things are pairs called ABCD_EFG_HIJ_KLM_1.fq & ABCD_EFG_HIJ_KLM_2.fq 

or  ABCD_EFG_HIJ_KLM_R1.fq & ABCD_EFG_HIJ_KLM_R2.fq 

or  ABCD_EFG_HIJ_KLM_001.fq & ABCD_EFG_HIJ_KLM_002.fq 

The last field  before ".fq" has to be the part that indicates whether the file contains R1 or R2 partners.

2. R version 3.0.2 - or a later version but with the package "data.table" installed
3. Bowtie2 
4. python version 2.7.5

Before you use:
Please download "PICKLES" -- the HG19 annotated exon information created for the purpose of this script

For Sherlock users, this is the directory "/scratch/PI/horence/gillian/HG19exons/"
For Sherlock users, please download the REG INDELS INDICES from - /scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/

For other users, please email glhsieh@stanford.edu for this directory.

Requires paired end read data.


Please open and change the script createFarJunctions_SLURM.sh

line31 - change your INSTALLDIR to the full path to the MACHETE script.  Last char must be "/"
line32 - change CIRCREF to the path to the reference libraries generated for KNIFE e.g. directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices.
line41 - PICKLEDIR - please change this path to the path that you used to store the HG19exons directory above

Running MACHETE:
First run KNIFE script completely to generate linear and scrambled junction reports and alignments.

sh createFarJunctions_SLURM.sh <1. KNIFE parent directory> <2. output directory> <3. discordant read distance> <4. ref genome> <5. junction overlapping reads> <6. #indels to use> <7. indel junction overlapping reads> 

NOTE -- ALL directory inputs must end in "/"

1. KNIFE parent directory - contains output from the KNIFE algorithm.  This directory is the path to the directory that contains "circReads", "orig", "logs", "sampleStats"
2. output directory - must already be in existence.
3. discordant read distance -- In testing mode, using 100000 base pairs to identify paired end reads that aligned discordantly.  
4. "HG19" is the only option currently
5. Currently using KNIFE convention of "8" for files with read lengths < 70 and "13" for files with read lengths > 70
6. <optional for sherlock use> -- "owners" if you want to run in owners queue, otherwise leave #6 blank


Example command:
sh createFarJunctions_SLURM.sh /scratch/PI/horence/alignments/EWS_FLI_bigmem/ /scratch/PI/horence/alignments/EWS_FLI_bigmem/FarJunc/ 100000 HG19 13 owners 


======================================================

DEBUGGING THE CODE
once you run the MACHETE, you will notice a file in the MACHETE output directory called "MasterError.txt"
Below is the workflow for the MACHETE and there are several parallel steps and several sequential steps.  All of these must complete for each unique pair of files for the MACHETE to be complete.


SEQUENCE1 and SEQUENCE 5 are started in parallel


SEQUENCE 1 
j1 - sort reg
j2 - sort genome
j3 - PE finder.sh
j4 - DistantPE_counter
j5 - sortPairedEnds
j6 - make Junctions
j6a - linkFastaFiles

From j6, three separate sequences (SEQUENCE 2, 3, 4) are called in parallel

SEQUENCE 2
j7 - LenientBadFJ.sh
dependstr7 - BadFJ & BadFJver2


SEQUENCE 3
j8 - alignUnalignedFJ
j9 - FJNaiveRept
j15a - parse_FJ_ID_for_GLM


SEQUENCE 4
j10 - makeIndelFiles
j11 - BowtieIndexFJindels
J13 - BowtieAlignFJIndels
j14 - FindAlignmentArtifact
j19 - FJIndelsClassOutput


SEQUENCE 5 (concurrently with SEQUENCE 1)
j16 - AlignUnalignedtoRegIndel-> j17 - AddIndelstoLinearGLM

j18 - regIndelsClassOutput


When SEQUENCE 3-5 are complete then the script j15b - run_GLM.sh is called


When SEQUENCE2 and j15b - run_GLM.sh are complete -- j15 - AppendNaiveRept.sh is called. This completes the MACHETE.  

=========================================
Output files

 MACHETE output dir/reports/AppendedReports - these are the final output files of MACHETE containing p values and posterior probabilities for reads.
1. naive_report_Appended.txt -- The first part of the report sequentially opens the far junctions read 1 and 2 alignments, genome read 1 and 2 alignments, reg (linear junction) 1 and 2 alignments, junc (scrambled junction) 1 and 2 alignments, and the unaligned files.
Alignments are designated as "true" or "false". Each column - genome, genome anomaly, genome p value, reg, reg anomaly, reg p value, etc...) refers to the count of alignments of a read partner of FJ to align to one of those libraries, or the p value. All true alignments from each category (genome, reg, junc) are used to calculate p values using the read length and the alignment score.  Using an estimated mismatch rate of 0.01 (per Illumina specifications), the mismatch rate is compared to the estimated mismatch rate using a Poisson cdf.
First FJ_1 is opened and compared to FJ2 to look for read partners.  If a read pair R1 and R2 in these files align to the exact same junction, those are considered a "true" FJ alignment.  If a read pair R1 and R2 do not align to the same junction, these are a "false" alignment.
If read partners R1 or R2 are not found in the FJ Directory, then the genome alignment files are opened.  If R1 is in FJ and R2 in genome, then the following criteria are used to designate those partners as "true" -- R2 must occur on the same chromosome as either upstream of the 5' gene or downstream of the 3' gene, and the two reads must have aligned to opposite reference strands.  If these criteria are not all met, then the read is considered "false". The same is done for FJ R2 and genome R1.
If the read partner is not in genome, the next library opened is the linear (reg) alignments.  For R1 in FJ and R2 in reg, the following criteria are used to designate read partners as "true".  For the location criterion, the downstream exon of the reg R2 must occur on the same chromosome and upstream of the upstream exon of the FJ R1 or the upstream exon of the reg R2 must occur on the same chromosome and downstream of the downstream exon of FJ R1.  Whichever of the R1 and R2 exons that meet the location criterion must also be on the same strand (+ or -).  Because the KNIFE sequences for (-) stranded exons are reverse complemented whereas the MACHETE sequences are not reverse complemented, for (-) R1 and R2 exons the reference strands must match (0 and 0 or 16 and 16, indicating that alignments were made to opposite sides of the cDNA), but for (+) R1 and R2 exons, the reference strands must be opposite (0 and 16). The same is done for FJ R2 and reg R1.
If the read partner R2 is not in FJ, genome, or reg, the partner is searched for in scrambled junction alignments.  For R1 in FJ and R2 in junc, all reads are considered "false".  The same is done for FJ R2 and junc R1.
If a FJ R1's read partner R2 is in none of the above, the unaligned reads are searched for R2. The same is done for FJ R2 and unaligned R1.
All reads where R1 or R2 was in FJ and the partner is in none of the above are tagged as "unmapped".  One caveat - if a read was found in any of the junctional alignment files (FJ, reg, or junc) but did not overlap the junction by the below specified overlap (-w ${3}) then that read is not considered.  For example, for FJ R1 that did overlap the junction by the desired number of bases and genome R2 that did NOT overlap the junction by the desired number of bases, the script would continue to search for R2 in reg, junc, and unaligned.  If the R2 does not occur in any of those other alignment files, then R2 would be labeled "unmapped" even though it technically did map to a genome alignment.

The Indels columns refer to the number of reads that aligned to a far junction : the number of times a read aligned to one of the indel indices. These are used for the GLM.

If there is homology between a FarJunctions fasta sequence and the genome or transcriptome or a linear junction or circular junction, then the fusion read is less likely.  Alignments of the FarJunctions fasta sequences to the KNIFE reference indices, genome, transcriptome, linear junctions (reg), and scrambled junctions (junc) are created with two different bowtie parameters.  Bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align. These are just aligning the FJ Fasta to the bad juncs with various alignment parameters. Any junctions aligning to here will eventually be tagged as "BadFJ=1" in the final reports whereas if junctions don't align, they will receive a "BadFJ=0" in the final reports.
For BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:  A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  For BadFJ ver2 we use bowtie to align the reads1 and 2 as if they were paired end reads from the same strand.  We impose a minimum gap of 0 between the two and a maximum gap of 50,000 bases to allow up to a 50,000 base gapped alignment.
ExonL and ExonR refer to the number of times that the left or right exon of that particular fusion participated in a linear junction.
The last several columns refer to the posterior probability of the junction, as calculated by the GLM.  Currently we are using both a 3 coefficient posterior probability (p_predicted.x) and a 5 coefficient posterior probability (p_predicted.y). 
2.  Other files in the AppendedReports directory include a circJunc reports and linearJunc reports.  These are parsed from the KNIFE.  The purpose of this is to place all reports in a single directory for the user.  KNIFE reports junctions up to 100Kb apart and MACHETE reports junctions beyond the user setting, which in our default, has been 100Kb. For the linear glmReports from KNIFE, the script collects any junctions where 1) the two exons from the linear report are from different genes or 2) the posterior probability is >0.9.  It adds on the rate of anomaly reads and indels to the reports and feeds them into FJDir/reports/AppendedReports.  For circular reports, the script collects any junctions where the posterior probability is <0.9, appends the "Decoy" rate, and feeds the reports into FJDir/reports/Appended reports. 
