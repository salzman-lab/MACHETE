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
parser.add_argument("-g", "--glmReportsDir", required= True, help = "linear junc glmReports Directory")
parser.add_argument("-s", "--stem", required=True, help = "file stem")

args=parser.parse_args()

if args.FarJuncDir[-1] != "/":
    args.FarJuncDir+= "/"

if args.glmReportsDir[-1] != "/":
    args.glmReportsDir+="/"

indelsDir = args.FarJuncDir + "IndelsHistogram/"
reportsDir = args.FarJuncDir + "reports/"
outputDir = args.FarJuncDir + "reports/AppendedReports/"
BadFJDir = args.FarJuncDir + "BadFJ/" + args.stem + "/"
NoGLM=False

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


## Parse thru FJ GLM and load results into a dictionary for appending onto the naive report
FJ_GLM_Dict={}

for FJ_GLM_file in glob.glob(reportsDir+"glmReports/"+args.stem+"*FUSION_W_ANOM*"):
    GLMfile=open(FJ_GLM_file, mode="rU")
    
    for line_raw in GLMfile:
        if line_raw[0:4]=="junc":
            FJ_GLM_Dict["header"]=line_raw.strip()[9:]
            continue
        
        junc=line_raw.strip().split("\t")[0]
        value=line_raw.strip()[len(junc)+1:]           
        
        FJ_GLM_Dict[junc]=value        
       
    GLMfile.close()

## parses through the report file and adds # indels.

for reportfile in glob.glob(reportsDir+ "*" + args.stem + "*naive*.txt"):

    f1 = open(reportfile, mode = "rU")
    indelfiles = sorted(glob.glob(indelsDir+ "*" + args.stem + "*.txt"))
    
    f2 = open(indelfiles[0], mode="rU")  # the indels  histo _1 file
    f3 = open(indelfiles[1], mode="rU")  # the indels  histo _2 file

    glmReportfile = glob.glob(args.glmReportsDir + "*" + args.stem + "*linear*.txt")

    if glmReportfile==[]:
        NoGLM=True
    else:
        f4 = open(glmReportfile[0], mode="rU")
    
    fout=open(outputDir+ args.stem + "_naive_report_Appended.txt", mode="w")
    
    
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

    if NoGLM==False:
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
            fout.write(line_raw.strip()+"\t_1NoIndels:Indels\t_2NoIndels:Indels\tBadFJ=1\tExonL\tExonR\t"+ FJ_GLM_Dict["header"]+"\n")
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
        if NoGLM==True:
            exonL="No glmReport"
            exonR="No glmReport"
        
        if junc in JunctionIndelsDict:
            fout.write(line_raw.strip()+"\t"+JunctionIndelsDict[junc][0]+"\t"+JunctionIndelsDict[junc][1]+"\t"+ BadFJ + "\t" + exonL + "\t" + exonR + "\t")
        else:
            fout.write(line_raw.strip() +"\t-\t-\t"+ BadFJ + "\t" + exonL + "\t" + exonR + "\t" ) 
            
        
        if junc in FJ_GLM_Dict:
            fout.write(FJ_GLM_Dict[junc]+"\n")
        else:
            fout.write("-\t-\t-\t-\t-\t-\t-\n")

    fout.close()
    