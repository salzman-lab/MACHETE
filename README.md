# MACHETE
More Accurate vs Mismatched Alignment CHimEra Tracking Engine

MACHETE is a fusion detection software that is in development.

PREREQUISITE SOFTWARE: 
1. KNIFE - https://github.com/lindaszabo/KNIFE
2. R version 3.0.2 - or a later version but with the package "data.table" installed
3. Bowtie2 
4. python version 2.7.5
5. SLURM job scheduler

INSTRUCTIONS BEFORE YOU RUN
1. Paired fastq files must named identically and end in _1.fq and _2.fq
Example:
ACCEPTABLE -- MySample_1.fq, MySample_2.fq
UNACCEPTABLE -- MySample1_001.fq, MySample2_001.fq

2.Generate or download a directory of necessary pickles
Pickles are a method of storing serialized data in Python.  We use pickles to store annotated exon information.  

A: Generating the pickles directory from a gtf and fasta file
Use makeExonDB.py to create the pickles directory.  The OutputDirectory is the pickles directory.
usage: python makeExonDB.py -f <genome.fasta> -a <genome.gtf> -o <OutputDirectory>

This can be adapted to any build of any genome.  The downloadable version is the hg19 genome
Genome.fasta is the name of the fasta file used, ex: hg19_genome.fasta
Genome.gtf is the name of the gtf file used, ex: hg19_genome.gtf
OutputDirectory is a path to the pickle directory which is chosen by the user

B: Downloading the pickles directory
At Stanford -- For Sherlock users who are using the HG19 genome, copy or point MACHETE at this directory "/scratch/PI/horence/gillian/HG19exons/"
Outside of Stanford - please email glhsieh@stanford.edu.  The Pickles directory is too large to fit on GitHub.

3. Generate or download an index of indels made from the KNIFE linear junction index. 
Generate the index.  The RegIndelsIndices will be created in a subfolder of the OutputDirectory called “IndelIndices”.
 Run MakeRegIndelsIndex.sh using the following command

Sh MakeRegIndelsIndex.sh <linear junctions fasta> <output directory> <# indels desired> <genome name> <Resource flag (optional)>

Linear junctions fasta is the path to the fasta file containing all linear junctions that is created for KNIFE
Output directory is the path to a directory where the user plans to store linear junction indels
# indels desired is the integer number of indels that will be used in searching for alignment artifact. We have chosen 5
Genome is the name of the genome that is being used. This is only used to name output files.  For example, if “HG19” is entered, output files will be named HG19_reg_indels_1.fa, etc.
Resource flag is an optional field for users of the Stanford SLURM network to specify which queue should be requested, eg “-p owners” or “-p horence” or can be left blank.

Download the index
At Stanford -- for Sherlock users who are using the HG19 genome, copy or point MACHETE at this directory “/scratch/PI/horence/gillian/HG19_reg_indels/IndelIndices/”
Outside of Stanford - please email glhsieh@stanford.edu.  The RegIndels directory is too large to fit on GitHub

PREPARING THE SHELL SCRIPT
Open createFarJunctions_SLURM.sh

line40 - change your INSTALLDIR to the full path to the MACHETE script. 
line42 - change CIRCREF to the path to the reference libraries generated for KNIFE e.g. directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices.
line44  - change REGINDELINDICES to the directory above under “INSTRUCTIONS BEFORE YOU RUN”, bullet #3.  This path should end with “IndelIndices”
line59 - change PICKLEDIR to the directory above under “INSTRUCTIONS BEFORE YOU RUN”, bullet #2


USING A DIFFERENT GENOME OTHER THAN HG19
Change line 44 REGINDELINDICES to the chosen new genome linear junction indices that were generated in “INSTRUCTIONS BEFORE YOU RUN”, bullet #3.
either change or duplicate lines 57 to lines 60.  Choose a genome name to replace “HG19”, and change the PICKLEDIR as needed.  
Either change or duplicate lines 200-206.  Choose a genome name to replace “HG19” and change the paths to the KNIFE generated reference indices to reflect the genome change.

RUNNING MACHETE:

First run KNIFE script completely to generate linear and scrambled junction reports and alignments.

sh createFarJunctions_SLURM.sh <1. KNIFE parent directory> <2. output directory> <3. discordant read distance> <4. ref genome> <5. #indels to use> <6. special queue> 

1. KNIFE parent directory - contains output from the KNIFE algorithm.  This directory is the path to the directory that contains "circReads", "orig", "logs", "sampleStats"
2. Output directory - if not already existing, will be created.
3. discordant read distance -- For testing, we have used 100000 base pairs to identify paired end reads that aligned discordantly.  
4. "HG19" is the only option currently.  If you are have created HG38 or another organism, pickles and a linear junction index must be generated above, instead of downloaded.  Additionally, lines 200-204 can be changed to point at KNIFE reference indices for that organism.
5. Currently using KNIFE convention of "8" for files with read lengths < 70 and "13" for files with read lengths > 70
6. <optional for sherlock use> -- "owners" if you want to run in owners queue, otherwise leave #6 blank


Example command:
sh createFarJunctions_SLURM.sh /scratch/PI/horence/alignments/EWS_FLI_bigmem/ /scratch/PI/horence/alignments/EWS_FLI_bigmem/FarJunc/ 100000 HG19 13 owners 


For the an explanation of outputs of MACHETE, please reference our paper:

Statistical algorithms improve accuracy of gene fusion detection
Gillian Hsieh, Rob Bierman, Linda Szabo, Alex Gia Lee, Donald E. Freeman, Nathaniel Watson, E. Alejandro Sweet-Cordero, Julia Salzman
Nucleic Acids Res. 2017 Jul 27; 45(13): e126. Published online 2017 May 24. doi: 10.1093/nar/gkx453 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5737606/.
