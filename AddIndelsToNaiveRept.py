# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:13:58 2015

@author: Gillian
"""


# this function takes a [a,b,c,x,d,e,f] and sums a-f vs x

def addindels(bracketedstring):
    IndelArray = bracketedstring.replace("[","").replace("]","").split(",")
    midpoint= len(IndelArray)/2
    NumJuncAlignments=int(IndelArray[midpoint])
    NumIndels=0
    for entry in IndelArray:
        NumIndels+=int(entry)
    NumIndels=NumIndels-NumJuncAlignments
    return NumJuncAlignments, NumIndels


import glob
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--FarJuncDir", required=True, help = "Far Junction Directory")
args=parser.parse_args()

indelsDir = args.FarJuncDir + "IndelsHistogram/"
reportsDir = args.FarJuncDir + "reports/"
outputDir = args.FarJuncDir + "reports/withIndels/"

for reportfile in glob.glob(reportsDir+"*_naive_report.txt"):
    stemlist = []
    path,filename=os.path.split(reportfile)
    if filename[:-17] not in stemlist:
        stemlist+=filename[:-17]

    f1 = open(reportfile, mode = "rU")
    f2 = open(indelsDir + "unaligned_" + filename[:-17] + "_1.txt")
    f3 = open(indelsDir + "unaligned_" + filename[:-17] + "_2.txt")
    fout=open(outputDir+ filename[:-17] + "_naive_report_withIndels.txt", mode="w")
    
    
    #populate a junction dictionary so each key is the junction name
    # each value is a string, where entry 0 is the _1 file indels
    # and entry 1 is the _2 file indels
    
    JunctionIndelsDict={}

    for line_raw in f2:
        line = line_raw.strip().split("\t")
        if line[0] not in JunctionIndelsDict:
            JunctionIndelsDict[line[0]]=["0","0"]
        numNoIndels, numIndels= addindels(line[1])
        JunctionIndelsDict[line[0]][0]= str(numNoIndels)+":"+str(numIndels)
    
    for line_raw in f3:
        line = line_raw.strip().split("\t")
        if line[0] not in JunctionIndelsDict:
            JunctionIndelsDict[line[0]]= ["0","0"]
        numNoIndels, numIndels= addindels(line[1])
        JunctionIndelsDict[line[0]][1] = str(numNoIndels)+":"+str(numIndels)
    
    for line_raw in f1:
        if line_raw[0] =="@":
            fout.write(line_raw.strip()+"\t"+"_1 NoIndels:Indels \t_2 NoIndels:Indels\n")
        junc= line_raw.strip().split("\t")[0]
        
        if junc in JunctionIndelsDict:
            fout.write(line_raw.strip()+"\t"+JunctionIndelsDict[junc][0]+"\t"+JunctionIndelsDict[junc][1]+"\n")
        else:
            fout.write(line_raw)
    