# MACHETE
More Accurate vs Mismatched Alignment CHimEra Tracking Engine

MACHETE is a fusion detection software that is in development.

REQUIREMENTS: 

1. KNIFE - https://github.com/lindaszabo/KNIFE
   NOTE *** - before running KNIFE, please make sure your fastq files are in the following format:
<UNIQUE STRING>_<read1 or 2 file notation>.fastq
some examples of acceptable things are pairs called ABCD_EFG_HIJ_KLM_1.fq & ABCD_EFG_HIJ_KLM_2.fq 
or  ABCD_EFG_HIJ_KLM_R1.fq & ABCD_EFG_HIJ_KLM_R2.fq 
or  ABCD_EFG_HIJ_KLM_001.fq & ABCD_EFG_HIJ_KLM_002.fq 
The last field  before ".fq" has to be the part that indicates whether the file contains R1 or R2 partners.

2. R version 3.0.2
3. Bowtie2 
4. python version 2.7.5

Before you use:
Please download "PICKLES" -- the HG19 annotated exon information created for the purpose of this script

For Sherlock users, this is the directory "/scratch/PI/horence/gillian/HG19exons/"

For other users, please email glhsieh@stanford.edu for this directory.

Requires paired end read data.
Matching paired files must have identical names except for that which denotes them as read 1 or read 2. For example "abcdef_1.fq" and "abcdef_2.fq" OR "abcedf_R1.fq" and "abcdef_R2.fq"


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
6. Currently using "5" but this does not affect our stats at this time
7. Currently using "13" as an empirically derived convention.


Example command:
sh createFarJunctions_SLURM.sh /scratch/PI/horence/alignments/EWS_FLI_bigmem/ /scratch/PI/horence/alignments/EWS_FLI_bigmem/FarJunc/ 100000 HG19 13 5 13
