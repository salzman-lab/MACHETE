# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:32:18 2016

@author: Gillian
"""

import argparse
import os


parser=argparse.ArgumentParser()
parser.add_argument("-i", "--inputFasta", required=True, help="fasta file for processing")
parser.add_argument("-l", "--ReadLength", required=True, help= "length of R1 and R2")
parser.add_argument("-o", "--outputDir", required=True, help = "directory for output")
args=parser.parse_args()

if args.outputDir[-1]!="/":
    args.outputDir+="/"


basename= os.path.basename(args.inputFasta)

if basename[-3:]==".fa":
    basename=basename[:-3]
if basename[-6:]==".fasta":
    basename=basename[:-6]
    
    
FastaOut_R1= open(args.outputDir+basename+"_R1.fa", mode = "w")
FastaOut_R2= open(args.outputDir+basename+"_R2.fa", mode = "w")

FastaFile = open(args.inputFasta, mode="rU")
for line in FastaFile:
    if line[0]==">":
        FastaOut_R1.write(line.strip()+"\n")
        FastaOut_R2.write(line.strip()+"\n")
        continue

    line=line.strip().replace("N","")
    FastaOut_R1.write(line[0:int(args.ReadLength)]+"\n")
    FastaOut_R2.write(line[-int(args.ReadLength):]+"\n")

FastaFile.close()
FastaOut_R1.close
FastaOut_R2.close()
