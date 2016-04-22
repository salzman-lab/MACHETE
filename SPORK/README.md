# SPORK
Acronym: **S**mall **P**owerful **O**rthogonal **R**ead mapper from **K**[NIFE](https://github.com/lindaszabo/KNIFE "KNIFE")

Based on the implementation of the denovo pipeline as discussed in [Szabo 2015](http://www.genomebiology.com/2015/16/1/126 "Szabo 2015 Genome Bio")

## Overview
SPORK takes as input either single or paired read data in [fastq format](http://maq.sourceforge.net/fastq.shtml "fastq format reference") 
and creates statistically based spliced alignments using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml "Bowtie2") 
as the underlying read mapper. SPORK is currently reliant on KNIFE to preprocess RNASeq data and remove reads that map completely to 
a given reference or previously annotated splice sites. SPORK will likely report a huge amount of "*de novo*" junctions if it is run on 
and unfiltered fastq file and actual hits will be obscured.

The algorithm splits each read into thirds. The 1st and 3rd third are attempted to be mapped to the reference without allowing gaps. 
If 

## Installation
The following installation steps assume the KNIFE has been successfully installed and has been used
to produce the necessary unaligned inputs.

## Expected Input Directory Tree
Since SPORK operates on the unaligned output of KNIFE it requires a specific directory structure of the input.

## Parameters
There are six mandatory parameters to SPORK:

1. Parent directory
2. Run Name

## Adjustable constants
1. Bin size
2. Junction number spanning bases

## Script specifics
SPORK is a small package containing two main python files and class util files.

The code is run from denovo_pipeline.py which calls denovo_utils.py and denovo_consensus_utils.py.

There are five very small class files in the denovo_class_utils directory:

1. Junction.py
2. BinPair.py
3. SAMEntry.py
4. FastQEntry.py
5. GTFEntry.py

#### denovo_pipeline.py
This is the main script that gets called to start the program.

#### denovo_utils.py
This is the largest script and provides functionality to the main program.

## Contact
This code was developed and is maintained by Rob Bierman: *rbierman[at]stanford.edu*

Please contact me with any questions or suggestions.

