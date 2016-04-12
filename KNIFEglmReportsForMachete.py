# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 21:16:37 2016

@author: Gillian
"""


import argparse
import glob

class juncInfo:
    def __init__(self, line_raw):
        line=line_raw.strip().split("\t")
        self.junction=line[0]
        self.leftgene=line[0].replace(":","|").split("|")[1]
        self.rightgene=line[0].replace(":","|").split("|")[3]
        self.posterior=line[2]

parser=argparse.ArgumentParser()
parser.add_argument("-l", "--glmReportsDir", required=True, help = "KNIFE glmReports dir")
parser.add_argument("-f", "--FJDir", required=True, help = "FarJunc directory")
parser.add_argument("-s", "--stem", required=True, help = "unique stem in file names" )
args=parser.parse_args()

if args.glmReportsDir[-1]!="/":
    args.glmDir+="/"
if args.FJDir[-1]!="/":
    args.FJDir+="/"

KNIFEfiles=[]
for filename in glob.glob(args.glmReportsDir+"*"):
    if args.stem in filename:
        KNIFEfiles.append(filename)

fout= open(args.FJDir+"reports/AppendedReports/"+args.stem+"_circ_and_linear_Juncs.txt", mode="w")
fout.write("junction\tnumReads\tp_predicted\tp_value\n")
for filename in KNIFEfiles:
    f1= open(filename, mode="rU")
    for line in f1:
        if line[0:8]=="junction":
            continue            
        JuncLine=juncInfo(line)
        if JuncLine.leftgene==JuncLine.rightgene:
            continue            
        fout.write(line.strip()+"\n")            
    f1.close()

fout.close()