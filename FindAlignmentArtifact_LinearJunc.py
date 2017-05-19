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
parser.add_argument("-s", "--samDir", required=True, help = "path to sam files of unaligned reads after alignment to linear junction indels directory; dir that you want to remove overlapping reads from")
parser.add_argument("-c", "--circReadsDir", required=True,help="path to circReads directory containing glmReports and reports directories")
parser.add_argument("-o", "--outDir", required=True, help = "output directory where you want to put linear junc indel CDF files")
parser.add_argument("-H", "--hiPPcutoff", required=True, help = " collects reads above this high posterior probability cutoff")
parser.add_argument("-L", "--lowPPcutoff", required=True, help = "collects reads below this low posterior probability cutoff")
parser.add_argument("-x", "--overlap", required=True, help = "alignment maps to at least -o reads on each side of the junction breakpoint")

args = parser.parse_args()


if args.samDir[-1] != "/":
    args.samDir += "/"
    
if args.circReadsDir[-1] != "/": 
    args.circReadsDir += "/"


HiPPCutoff= float(args.hiPPcutoff)
LoPPCutoff= float(args.lowPPcutoff)
BinSize = 0.01

glmReportsdir = args.circReadsDir+"glmReports/"
naivereportdir = args.circReadsDir + "reports/"

NonOverlapDir = args.samDir + "RemoveNonOverlap/"


# remove non-junctional reads and place in nonoverlap


for filename in glob.glob(args.samDir+"*.sam"):
    (inpath,infile)=os.path.split(filename)
    
    f1 = open(filename, mode ="rU")
    fout = open(NonOverlapDir + infile, mode ="w")
    print "opened file: " + filename
    
    
    for line in f1:
        while line[0] == "@":
            fout.write(line)
            line = f1.next()
        
        line_list = line.strip().split("\t")
        if int(line_list[3]) > (150-int(args.overlap)) or (int(line_list[3])+len(line_list[9])) <  (150+int(args.overlap)):
            #print "range = " + line_list[3] + "-" + str(int(line_list[3])+len(line_list[9]))
            pass
        else:
            #print "overlap contained"
            fout.write(line)


f1.close()
fout.close()




# Get info from glm file, naive file, and indels
# get ratio of # insertions, deletions, and/or anomaly reads to # true linear reads


StemDict = {}

for filename in glob.glob(NonOverlapDir+"*.sam"):
    head, tail = os.path.split(filename)
    stem = tail.split("_")[1]
    if stem not in StemDict:
        StemDict[stem]=0
    StemDict[stem]+=1

print StemDict


for stem in StemDict:

    for name in glob.glob(glmReportsdir + "*" + stem + "*linear*.txt"):
        glmrptfile = name
    for name in glob.glob(naivereportdir + "*" + stem + "*.txt"):        
        if "denovo" not in name:
            naivereportfile = name
    indelfile = NonOverlapDir + "*" + stem + "*.sam"

        # KEY - junction name - FROM GLM FILE
        # value column 0 - posterior probability - FROM GLM FILE
        # column 1 - number of occurrences of junction - NAIVE RPT
        # column 2 - number of anomaly reads - NAIVE RPT
        # column 3 - number of insertions - FROM INDELS HISTOGRAM
        # column 4 - number of deletions - FROM INDELS HISTOGRAM
        # column 5 - number of insertions + deletions - CALCULATE
        # column 6 - number of insertion/deletion/anomalies - CALCULATE
    
    # open  file and feed all junctions into a junction dictionary  
    f_glm = open(glmrptfile, mode ="rU")
    f_naive = open(naivereportfile, mode = "rU")
    HiProbJunctionDict={}
    LoProbJunctionDict={}
    
    for line_raw in f_glm:
        if line_raw.strip().split("\t")[0]=="junction":
            line_raw = f_glm.next()
            
        junction=line_raw.strip().split("\t")[0]
        postprob=line_raw.strip().split("\t")[2]
        postprob_float = float(postprob)
        if postprob_float >= HiPPCutoff:
            if junction not in HiProbJunctionDict:
                HiProbJunctionDict[junction]=[0.0]*7
            HiProbJunctionDict[junction][0]+= postprob_float
        elif postprob_float <= LoPPCutoff:
            if junction not in LoProbJunctionDict:
                LoProbJunctionDict[junction]=[0.0]*7
            LoProbJunctionDict[junction][0]+= postprob_float
           
    
    f_glm.close()
    
    for line_raw in f_naive:
        while line_raw[0] =="@":
            line_raw = f_naive.next()
    
            
        junction=line_raw.strip().split("\t")[0]
        linear_occurrence = int(line_raw.strip().split("\t")[1])
        anomaly_occurrence = int(line_raw.strip().split("\t")[2]) + int(line_raw.strip().split("\t")[3])
    
    
        if junction in HiProbJunctionDict:
            
            HiProbJunctionDict[junction][1] += linear_occurrence
            HiProbJunctionDict[junction][2] += anomaly_occurrence
        elif junction in LoProbJunctionDict:
            LoProbJunctionDict[junction][1] += linear_occurrence
            LoProbJunctionDict[junction][2] += anomaly_occurrence
    
    f_naive.close()
    
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
    
    
    
        # KEY - junction name - FROM GLM FILE
        # value column 0 - posterior probability - FROM GLM FILE
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
        if string[1]+string[2]+string[3]+string[4] <= 10:
            BadHiProb[junction]=0
        else:
            HiProbJunctionDict[junction][5]= HiProbJunctionDict[junction][3]+HiProbJunctionDict[junction][4]
            HiProbJunctionDict[junction][6]= HiProbJunctionDict[junction][2]+HiProbJunctionDict[junction][3]+HiProbJunctionDict[junction][4]
        
    for junction in LoProbJunctionDict:
        string = LoProbJunctionDict[junction]
        if string[1]+string[2]+string[3]+string[4] <= 10:
            BadLoProb[junction]=0
        else:
            LoProbJunctionDict[junction][5]= LoProbJunctionDict[junction][3]+LoProbJunctionDict[junction][4]
            LoProbJunctionDict[junction][6]= LoProbJunctionDict[junction][2]+LoProbJunctionDict[junction][3]+LoProbJunctionDict[junction][4]
    
    for junction in BadHiProb:
        del HiProbJunctionDict[junction]    
    for junction in BadLoProb:
        del LoProbJunctionDict[junction]    
    
    ## print results
    fout = open(args.outDir + stem+ "_LinearAlignmentArtifact.txt", mode = "w")
    
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
    

    
    Hifout=open(args.outDir + stem + "_HiProb_Anomaly.txt", mode = "w")
    Lofout=open(args.outDir + stem + "_LoProb_Anomaly.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbAnomalyArray, LoProbAnomalyArray)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    
    Hifout=open(args.outDir + stem + "_HiProb_Insertion.txt", mode = "w")
    Lofout=open(args.outDir + stem + "_LoProb_Insertion.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbInsertionArray, LoProbInsertionArray)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    Hifout=open(args.outDir + stem + "_HiProb_Deletion.txt", mode = "w")
    Lofout=open(args.outDir + stem + "_LoProb_Deletion.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbDeletionArray, LoProbDeletionArray)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    
    Hifout=open(args.outDir + stem + "_HiProb_InsDel.txt", mode = "w")
    Lofout=open(args.outDir + stem + "_LoProb_InsDel.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbIns_Del_Array, LoProbIns_Del_Array)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()
    
    
    Hifout=open(args.outDir + stem + "_HiProb_InsDelAnom.txt", mode = "w")
    Lofout=open(args.outDir + stem + "_LoProb_InsDelAnom.txt", mode ="w")
    MakeTextHisto(Hifout, Lofout, HiProbIns_Del_Anom_Array, LoProbIns_Del_Anom_Array)
    Hifout.write("\n")
    Lofout.write("\n")
    Hifout.close()
    Lofout.close()



