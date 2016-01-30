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





## START PROGRAM


import glob
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--FarJuncDir", required=True, help = "Far Junction Directory")
parser.add_argument("-g", "--glmReportsDir", required= True, help = "glmReports Directory")

args=parser.parse_args()

if args.FarJuncDir[-1] != "/":
    args.FarJuncDir+= "/"

indelsDir = args.FarJuncDir + "IndelsHistogram/"
reportsDir = args.FarJuncDir + "reports/"
outputDir = args.FarJuncDir + "reports/withIndels/"
BadFJDir = args.FarJuncDir + "BadFJ/"


# parses through bad FJ alignments.  for every FJ that aligned to genome/transcriptome/reg/junc
# the FJ is stored in the Bad FJ dictionary, for tagging in the report file.
BadFJDictionary={}

for BadFJfile in glob.glob(BadFJDir+"*.sam"):
    f1=open(BadFJfile, mode="rU")


    for line_raw in f1:
        if line_raw[0]=="@":
            continue
    
        badjunction=line_raw.strip().split("\t")[0]
        BadFJDictionary[badjunction]=1
        
    f1.close()


## parses through the report file and adds # indels.

for reportfile in glob.glob(reportsDir+"*_naive_report.txt"):
    stemlist = []
    path,filename=os.path.split(reportfile)
    if filename[:-17] not in stemlist:
        stemlist+=filename[:-17]

    f1 = open(reportfile, mode = "rU")
    f2 = open(indelsDir + "unaligned_" + filename[:-17] + "_1.txt", mode="rU")
    f3 = open(indelsDir + "unaligned_" + filename[:-17] + "_2.txt", mode="rU")

    glmReportfile = glob.glob(args.glmReportsDir + "*" + filename[:-17] + "*linear*.txt")

    f4 = open(glmReportfile[0], mode="rU")
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
    
    f2.close()
    
    
    for line_raw in f3:
        line = line_raw.strip().split("\t")
        if line[0] not in JunctionIndelsDict:
            JunctionIndelsDict[line[0]]= ["0","0"]
        numNoIndels, numIndels= addindels(line[1])
        JunctionIndelsDict[line[0]][1] = str(numNoIndels)+":"+str(numIndels)

    f3.close()

## parse through reports and get all exons that participated.
    exonDict = {}
    for line_raw in f1:
        if line_raw[0]=="@":
            continue
        
        exonL=line_raw.strip().split("\t")[0].split(":")[1]       
        exonR=line_raw.strip().split("\t")[0].split(":")[4]
        
        if exonL not in exonDict:
            exonDict[exonL]=0
        if exonR not in exonDict:
            exonDict[exonR]=0
    f1.close()

## open GLMreports file and see what is the max # of times exon appeared.
    for line_raw in f4:
        if line_raw.strip().split("\t")[0]=="junction":
            continue


        # if posterior prob >0.8 then if L or R exon count is the greatest in exonDict, 
## then replace exonDict entry with new count    
   
        if float(line_raw.strip().split("\t")[2])>= 0.8: 
            exonL = line_raw.strip().split("\t")[0].replace("|",":").split(":")[1]
            exonR = line_raw.strip().split("\t")[0].replace("|",":").split(":")[3]
            
            if exonL in exonDict:
                if int(line_raw.strip().split("\t")[1]) > exonDict[exonL]:
                    exonDict[exonL]= int(line_raw.strip().split("\t")[1])
                    
            if exonR in exonDict:
                if int(line_raw.strip().split("\t")[1]) > exonDict[exonR]:
                    exonDict[exonR]= int(line_raw.strip().split("\t")[1])                    

    f4.close()


## OUTPUTTING- includes report file, no indel: indel ratio, quality of FJ (bad = 1, good = 0)
## and max # times that an exon appeared with good quality in GLM reports file.
    f1 = open(reportfile, mode = "rU")
    for line_raw in f1:
        if line_raw[0] =="@":
            fout.write(line_raw.strip()+"\t"+"_1 NoIndels:Indels \t_2 NoIndels:Indels\tBadFJ=1\tExonL\tExonR\n")
            continue

        junc= line_raw.strip().split("\t")[0]
            

        if junc in BadFJDictionary:
            BadFJ="1"
        else:
            BadFJ="0"
        
#        try: 
        exonL= str(exonDict[line_raw.strip().split("\t")[0].split(":")[1]])
#        except: 
 #           print line_raw            
            #print line_raw.strip().split("\t")[0].split(":")[1]
            #print exonDict[line_raw.strip().split("\t")[0].split(":")[1]]
            
        exonR=str(exonDict[line_raw.strip().split("\t")[0].split(":")[4]])      
        
        
        if junc in JunctionIndelsDict:
            fout.write(line_raw.strip()+"\t"+JunctionIndelsDict[junc][0]+"\t"+JunctionIndelsDict[junc][1]+"\t"+ BadFJ + "\t" + exonL + "\t" + exonR + "\n")
        else:
            fout.write(line_raw.strip() +"\t-\t-\t"+ BadFJ + "\t" + exonL + "\t" + exonR + "\n" ) 
    