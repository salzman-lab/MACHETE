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

Requires paired end reads but aligns in a single end fashion.


Please open and change the script createFarJunctions_SLURM.sh

line31 - change your INSTALLDIR to the full path to the MACHETE script.  Last char must be "/"
line32 - change CIRCREF to the path to the reference libraries generated for KNIFE e.g. directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices.
line41 - PICKLEDIR - please change this path to the path that you used to store the HG19exons directory above

Running MACHETE:
First run KNIFE script completely to generate linear and scrambled junction reports and alignments.

sh createFarJunctions_SLURM.sh <1. KNIFE parent directory> <2. output directory> <3. discordant read distance> <4. ref genome> <5. #indels to use> <6. special queue> 

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
Output directories

1. DistantPEfiles - 
<STEM>_distant_pairs.txt -- list of three columns - readID / position1 / position2.  Position1 and Position2 are the windows surrounding "discordant" bases from either the genome or reg alignment files.
    Subdirectory <STEM>:
    A. file chr(1,2,3,..,X,Y)_Distant_PE_frequency.txt.  If, for the pair of windows, the upstream location is on chrA, then it will be binned into the chrA_Distant_PE_frequency.txt file.  Three columns appear -- Window1, Window2, and the number of times that these two windows were matched. 
    B. file sorted_chr(1,2,3,..)_Distant_PE_frequency.txt -- File A, sorted on the numerical value of the first column followed by the numerical values of the second column, to reduce parsing time later.
2. fasta -
    Subdirectory <STEM>:
    A: <STEM>_chr(1,2,3,..)FarJunctions.fa - a fusion reference fasta file is created from the corresponding sorted_chr(1,2,3,...)Distant_PE_Frequency. 
    <STEM>_FarJunctions.fa -- the concatenated list of all fasta entries from fasta/<STEM>/<STEM>_chr*_FarJunctions.fa
3. BadFJ/<STEM> -
    A. err/out.txt- error and output files from bowtie alignment of fasta to reference dictionaries
    B. <STEM>_BadFJto(Genome/Reg/Junc/transcriptome).sam - aligns fasta/<STEM>_FarJunctions.fa to the KNIFE reference indices Genome, Reg, Junc, and transcriptome with bowtie parameters "-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --np 0 --rdg 50,50 --rfg 50,50"
4. BadFJ_ver2/<STEM> - 
    A. err/out.txt- error and output files from bowtie alignment of fasta to reference dictionaries
    B. <STEM>_FarJunctions_R1/2.fa - contains reads from fasta/<STEM>_FarJunctions.fa where all N's are removed from reads and the first 40 nt of the original fasta read is fed into the R1 file and the last 40 nt of the original fasta read is fed into the R2 file. 
    C. <STEM>_BadFJto(Genome/Reg/Junc/transcriptome).sam - aligns BadFJ_ver2/<STEM>_FarJunctions_R1/2.fa to the KNIFE reference indices Genome, Reg, Junc, and transcriptome as if they were paired end alignments to allow a read gap.  Bowtie parameters are "--no-unal --no-mixed --no-sq -p 8 -I 0 -X 50000 -f --ff" 
5. BowtieIndex - 
    Subdirectory <STEM>
    A: <STEM>_FJ_Index*bt2 -- corresponding bowtie indices are created for each of the fasta files from the fasta/<STEM>_chr(1,2,3...)_FarJunctions.fa
6. FarJunctionAlignments - 
    Subdirectory <STEM>
    A. unaligned_<STEM>_R1/2.sam - contains reads that were unaligned to the KNIFE reference indices that have been aligned to the bowtie indices in the BowtieIndex/<STEM>/ directory
7. FarJuncSecondary - 
    Subdirectory <STEM>
    A. still_unaligned_<STEM>_R1/2.fq - contains reads that were unaligned to the KNIFE reference indices that did not align to the FarJunctions bowtie indices in the BowtieIndex/<STEM>/ directory
    Subdirectory AlignedIndels/<STEM>/
    A. still_unaligned_<STEM>_R1/2_indels(1-5).sam - contains reads that were unaligned to the KNIFE reference indices, did not align to the FarJunctions bowtie indices, but DID align to an the expanded FarJunctions indels bowtie indices.
    B. All_<STEM>_R1/2_indels.sam -- concatenates all still_unaligned_<STEM>_R1/2_indels(1-5).sam files but removes reads that aligned to multiple indel indices and keeps the one with the best alignment score and removes reads that aligned to the indel indices but did not overlap the junction by the user specified number of nt.
8. FarJuncIndels

< "IF script doesn't have "INDELS" in glmReports dir means that indels were omitted " >"'




