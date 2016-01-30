# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:51:44 2015

@author: Gillian
"""
#import argparse
import os
import glob
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--FarJuncDir", required=True, help = "path to FarJunc directory")
parser.add_argument("-H", "--hiPPcutoff", required=True, help = " collects reads above this high posterior probability cutoff")
parser.add_argument("-L", "--lowPPcutoff", required=True, help = "collects reads below this low posterior probability cutoff")

args = parser.parse_args()


if args.FarJuncDir[-1] != "/":
    args.FarJuncDir += "/"


HiPPCutoff= float(args.hiPPcutoff)
LoPPCutoff= float(args.lowPPcutoff)
BinSize = 0.01

naivereportdir = args.FarJuncDir + "reports/"
indelsdir = args.FarJuncDir + "/FarJuncSecondary/AlignedIndels/RemoveNonOverlap/"

StemDict = {}

for filename in glob.glob(naivereportdir+"*_naive_report*.txt"):
    head, tail = os.path.split(filename)
    stem = tail[:-17] # THIS IS WITHOUT THE "with indels" tag that probably adds more onto the stem
    if stem not in StemDict:
        StemDict[stem]=0
    StemDict[stem]+=1

print StemDict


for stem in StemDict:
    naivereportfile=glob.glob(naivereportdir + "*" + stem + "*_naive_report*.txt")[0]
    print naivereportfile
    indelfile = indelsdir + "*" + stem + "*.sam"
    
        
    #TEST CASE 
    #HiPPCutoff = 0.95
    #LoPPCutoff = 0.1
    #BinSize = 0.01
    #
    #glmrptfile = "/Users/Gillian/Desktop/sherlock/regJuncIndels/C088AACX2CACGAT_1__linearJuncProbs.txt"
    #naivereportfile= "/Users/Gillian/Desktop/sherlock/regJuncIndels/C088AACX2CACGAT_1_report.txt"
        ## these are the indel files with the nonoverlapping reads removed.
    #indelfile = "/Users/Gillian/Desktop/sherlock/regJuncIndels/unaligned_C088AACX2CACGAT*indel*.sam"
    
    # open  file and feed all junctions into a junction dictionary

    f_naive = open(naivereportfile, mode = "rU")
    HiProbJunctionDict={}
    LoProbJunctionDict={}

    counter=0
    
    for line_raw in f_naive:
        counter+=1

        while line_raw[0]=="@":
     #       print line_raw[0]
            try: line_raw = f_naive.next()
            except: break
        
        
     #   print line_raw
        if line_raw[0]=="@": break
        
        junction=line_raw.strip().split("\t")[0]
        NetP=line_raw.strip().split("\t")[15]
        
        if NetP=="-":
            continue
        else:        
            NetP_float = float(NetP)
        
        if NetP_float >= HiPPCutoff:
            if junction not in HiProbJunctionDict:
                HiProbJunctionDict[junction]=[0.0]*7
            HiProbJunctionDict[junction][0]+= NetP_float
        elif NetP_float <= LoPPCutoff:
            if junction not in LoProbJunctionDict:
                LoProbJunctionDict[junction]=[0.0]*7
            LoProbJunctionDict[junction][0]+= NetP_float
            
#  "true occurrence" currently means partner aligns to genome or linear junc 
#  "anomaly occurrence" currently means partner aligns to any anomaly bin or unaligned or unmapped -- 
            #(genome anom, reg anom, junc, junc anom, far junc anom, unaligned bins combined)            
            
        true_occurrence = int(line_raw.strip().split("\t")[1])+int(line_raw.strip().split("\t")[4])
        anomaly_occurrence = int(line_raw.strip().split("\t")[2]) + int(line_raw.strip().split("\t")[5])+ int(line_raw.strip().split("\t")[7])+ int(line_raw.strip().split("\t")[8])+ int(line_raw.strip().split("\t")[11])+ int(line_raw.strip().split("\t")[13])
    
    
        if junction in HiProbJunctionDict:
            
            HiProbJunctionDict[junction][1] += true_occurrence
            HiProbJunctionDict[junction][2] += anomaly_occurrence
        elif junction in LoProbJunctionDict:
            LoProbJunctionDict[junction][1] += true_occurrence
            LoProbJunctionDict[junction][2] += anomaly_occurrence
    
    f_naive.close()


#if junction dictionary is empty, skip to the next naive file
    if HiProbJunctionDict=={} and LoProbJunctionDict=={}:
        continue


    for IndelFile in glob.glob(indelfile):
        #print IndelFile
        f_indel = open(IndelFile, mode = "rU")
        for line_raw in f_indel:
            while line_raw[0] == "@":
                line_raw = f_indel.next()
                #print line_raw
            junction_indel = line_raw.strip().split("\t")[2]
            junction = junction_indel[:-5]
    
            if junction in HiProbJunctionDict:
                if junction_indel[-4:-1]=="INS":
                    HiProbJunctionDict[junction][3] += 1
                elif junction_indel[-4:-1]=="DEL":
                    HiProbJunctionDict[junction][4] += 1
                                
            elif junction in LoProbJunctionDict:
                if junction_indel[-4:-1]=="INS":
                    LoProbJunctionDict[junction][3] += 1
                elif junction_indel[-4:-1]=="DEL":
                    LoProbJunctionDict[junction][4] += 1
    
    f_indel.close()
    
    
    
        # KEY - junction name - FROM GLM (CURR NAIVE) FILE
        # value column 0 - posterior probability - FROM GLM (CURR NAIVE) FILE
        # column 1 - number of occurrences of junction - NAIVE RPT
        # column 2 - number of anomaly reads - NAIVE RPT
        # column 3 - number of insertions - FROM INDELS HISTOGRAM
        # column 4 - number of deletions - FROM INDELS HISTOGRAM
        # column 5 - number of insertions + deletions - CALCULATE
        # column 6 - number of insertion/deletion/anomalies - CALCULATE
    
    
    
    ## removes junctions where anom + ins + del +  occurences <= 10
    ## for all other junctions, calculates ins+del & ins+del+anom categories
    
    BadHiProb={}
    BadLoProb={}
    
    for junction in HiProbJunctionDict:
        string = HiProbJunctionDict[junction]
        if string[1]+string[2]+string[3]+string[4] <= 2:
            BadHiProb[junction]=0
        elif string[1] == 0:
            BadHiProb[junction]=0
        else:
            HiProbJunctionDict[junction][5]= HiProbJunctionDict[junction][3]+HiProbJunctionDict[junction][4]
            HiProbJunctionDict[junction][6]= HiProbJunctionDict[junction][2]+HiProbJunctionDict[junction][3]+HiProbJunctionDict[junction][4]
        
    for junction in LoProbJunctionDict:
        string = LoProbJunctionDict[junction]
        if string[1]+string[2]+string[3]+string[4] <= 2:
            BadLoProb[junction]=0
        elif string[1] == 0:
            BadLoProb[junction]=0
        else:
            LoProbJunctionDict[junction][5]= LoProbJunctionDict[junction][3]+LoProbJunctionDict[junction][4]
            LoProbJunctionDict[junction][6]= LoProbJunctionDict[junction][2]+LoProbJunctionDict[junction][3]+LoProbJunctionDict[junction][4]
    
    for junction in BadHiProb:
        del HiProbJunctionDict[junction]    
    for junction in BadLoProb:
        del LoProbJunctionDict[junction]    
    
    ## print results
    fout = open(args.FarJuncDir + "FJIndelsCDF/"+stem+"_FJAlignmentArtifact.txt", mode = "w")
    #fout = open("/Users/Gillian/Desktop/sherlock/regJuncIndels/LinearAlignmentArtifact.txt", mode ="w")
    
    
    HiProbAnomalyArray=[] 
    HiProbInsertionArray=[]
    HiProbDeletionArray=[]
    HiProbIns_Del_Array=[]
    HiProbIns_Del_Anom_Array = []
    
    LoProbAnomalyArray=[] 
    LoProbInsertionArray=[]
    LoProbDeletionArray=[]
    LoProbIns_Del_Array=[]
    LoProbIns_Del_Anom_Array = []
    

    
        # KEY - junction name - FROM GLM (CURR NAIVE) FILE
        # value column 0 - posterior probability - FROM GLM (CURR NAIVE) FILE
        # column 1 - number of occurrences of junction - NAIVE RPT
        # column 2 - number of anomaly reads - NAIVE RPT
        # column 3 - number of insertions - FROM INDELS HISTOGRAM
        # column 4 - number of deletions - FROM INDELS HISTOGRAM
        # column 5 - number of insertions + deletions - CALCULATE
        # column 6 - number of insertion/deletion/anomalies - CALCULATE
    
        
    
    fout.write("Junction\t[NetPP, #true alignment, #anomaly, #insert, #del, #indels, #Ins/DelAnom]\n")
    fout.write("Junctions with Posterior Probability > " + str(HiPPCutoff) + "\n")
    for junction in HiProbJunctionDict:
        fout.write(junction+"\t"+str(HiProbJunctionDict[junction])+"\n")
        HiProbAnomalyArray.append( HiProbJunctionDict[junction][2]/HiProbJunctionDict[junction][1] )
        HiProbInsertionArray.append( HiProbJunctionDict[junction][3]/HiProbJunctionDict[junction][1])
        HiProbDeletionArray.append(HiProbJunctionDict[junction][4]/HiProbJunctionDict[junction][1])
        HiProbIns_Del_Array.append(HiProbJunctionDict[junction][5]/HiProbJunctionDict[junction][1])
        HiProbIns_Del_Anom_Array.append(HiProbJunctionDict[junction][6]/HiProbJunctionDict[junction][1])
    
    
    fout.write("\nJunctions with Posterior Probability <" + str(LoPPCutoff) + "\n")
    for junction in LoProbJunctionDict:    
        fout.write(junction+"\t"+str(LoProbJunctionDict[junction])+"\n")
        LoProbAnomalyArray.append( LoProbJunctionDict[junction][2]/LoProbJunctionDict[junction][1] )
        LoProbInsertionArray.append( LoProbJunctionDict[junction][3]/LoProbJunctionDict[junction][1])
        LoProbDeletionArray.append(LoProbJunctionDict[junction][4]/LoProbJunctionDict[junction][1])
        LoProbIns_Del_Array.append(LoProbJunctionDict[junction][5]/LoProbJunctionDict[junction][1])
        LoProbIns_Del_Anom_Array.append(LoProbJunctionDict[junction][6]/LoProbJunctionDict[junction][1])
    
    fout.close()
    # KEY - junction name - FROM GLM FILE
    # value column 0 - posterior probability - FROM GLM FILE
    # column 1 - number of occurrences of junction - NAIVE RPT
    # column 2 - number of anomaly reads - NAIVE RPT
    # column 3 - number of insertions - FROM INDELS HISTOGRAM
    # column 4 - number of deletions - FROM INDELS HISTOGRAM
    # column 5 - number of insertions + deletions - CALCULATE
    # column 6 - number of insertion/deletion/anomalies - CALCULATE
    
    
    def MakeTextHisto(HiOutfile, LoOutfile, HiProbArray, LoProbArray):
        HiProbHistogram = [0]*100
        for entry in HiProbArray:
            HistoBin = int(entry/BinSize)
            if HistoBin >= len(HiProbHistogram):
                for i in range(0,HistoBin-len(HiProbHistogram)+1):
                    HiProbHistogram.append(0)
            HiProbHistogram[HistoBin]+=1
            
        for i in range(0, len(HiProbHistogram)):
            HiOutfile.write(str(HiProbHistogram[i])+"\t")
        
        LoProbHistogram = [0]*100
        for entry in LoProbArray:
            HistoBin = int(entry/BinSize)
            if HistoBin >= len(LoProbHistogram):
                for i in range(0,HistoBin-len(LoProbHistogram)+1):
                    LoProbHistogram.append(0)
            LoProbHistogram[HistoBin]+=1
    
        for i in range(0, len(LoProbHistogram)):
            LoOutfile.write(str(LoProbHistogram[i])+"\t")
    
   
    Hifout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_HiProb_Anomaly.txt", mode = "w")
    Lofout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_LoProb_Anomaly.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbAnomalyArray, LoProbAnomalyArray)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    
    Hifout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_HiProb_Insertion.txt", mode = "w")
    Lofout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_LoProb_Insertion.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbInsertionArray, LoProbInsertionArray)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    Hifout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_HiProb_Deletion.txt", mode = "w")
    Lofout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_LoProb_Deletion.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbDeletionArray, LoProbDeletionArray)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    
    Hifout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_HiProb_InsDel.txt", mode = "w")
    Lofout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_LoProb_InsDel.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbIns_Del_Array, LoProbIns_Del_Array)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    
    Hifout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_HiProb_InsDelAnom.txt", mode = "w")
    Lofout=open(args.FarJuncDir + "FJIndelsCDF/" + stem + "_LoProb_InsDelAnom.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbIns_Del_Anom_Array, LoProbIns_Del_Anom_Array)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()



