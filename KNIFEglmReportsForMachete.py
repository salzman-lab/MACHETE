# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 21:16:37 2016

@author: Gillian
"""


import argparse
import glob


parser=argparse.ArgumentParser()
parser.add_argument("-c", "--circparentDir", required=True, help = "KNIFE output dir with circReads, orig, etc")
parser.add_argument("-f", "--FJDir", required=True, help = "FarJunc directory")
parser.add_argument("-s", "--stem", required=True, help = "unique stem in file names" )
args=parser.parse_args()


## if input directories don't end in "/", then add it
if args.circparentDir[-1]!="/":
    args.circparentDir+="/"
if args.FJDir[-1]!="/":
    args.FJDir+="/"

## defines various input directories
glmDir = args.circparentDir+"circReads/glmReports/"
naiveRptDir= args.circparentDir + "circReads/reports/"
RegIndelsDir = args.circparentDir + "orig/RegIndelAlignments/" + args.stem + "/"


## identify glmReports from KNIFE that need to be copied into FJ Appended reports
KNIFEfiles=[]
for filename in glob.glob(glmDir+"*"):
    if args.stem in filename:
        KNIFEfiles.append(filename)


## loop through Reg Indels Files and get all the junctions and indels.  These files
## have already previously been parsed to remove duplicate read IDs and alignments that 
## do not overlie the junction.

## open the RegIndels directory and find number of Reg Indels for each junction
## OF NOTE, Linda only uses the R1 file to calculate her # linear reads so 
## at this point I am omitting indels from R2 files so that I can report an accurate
## noIndel:Indel ratio.
RegIndels_1={}  #key = junction name, value = # indels in R1 file
#RegIndels_2={}  # key = junction name, value = # indels in R2 file

RegIndelsFile_1= open( sorted(glob.glob(RegIndelsDir+"All_*.sam"))[0], mode="rU")
#RegIndelsFile_2= open( sorted(glob.glob(RegIndelsDir+"All_*.sam"))[1], mode="rU")

for line in RegIndelsFile_1:
    line=line.strip().split("\t")
    if line[2][:-5] not in RegIndels_1:  
        RegIndels_1[line[2][:-5]]=0  
    RegIndels_1[line[2][:-5]]+=1  
    
#for line in RegIndelsFile_2:
#    line=line.strip().split("\t")
#    if line[2][:-5] not in RegIndels_2:
#        RegIndels_2[line[2][:-5]]=0
#    RegIndels_2[line[2][:-5]]+=1
    
RegIndelsFile_1.close()
#RegIndelsFile_2.close()

## loop through naive reports to extract rates of linear anomalies and circular decoys
AnomalyandDecoyDict={} # key = junction, value = [x,y] where x = # anomaly reads, y = # decoy reads

NaiveRptFile=open(glob.glob(naiveRptDir+args.stem+"_*_report.txt")[0], mode="rU")
for line in NaiveRptFile:
    if line[0]=="@":
        continue
    line=line.strip().split("\t")
    junction=line[0]
    AnomalyandDecoyDict[junction]=[line[2], line[6]]

NaiveRptFile.close()

## Open output directories
## take reg and circ GLM reports and append them before placing them into the FJ Appended
## reports directory




## criteria for reads from circular GLM report to be added:
##   1. posterior probability > 0.9 only
## CAVEAT -- circ reads entries could be either a circular RNA or an internal tandem duplication
## CAVEAT#2 -- inversions at <100KB are missed because Linda's index doesn't 
##              include inversions 
circ_output= open(args.FJDir+"reports/AppendedReports/"+args.stem+"_circJunc_reports.txt", mode="w")
circInputFile= open( glob.glob(glmDir+"*"+args.stem+"*circJuncProbs.txt")[0], mode="rU")

for line in circInputFile:
    if line[0:8]=="junction":
        circ_output.write(line.strip()+"\tdecoy\n")
        continue
    junction=line.strip().split("\t")[0]
    posterior=float(line.strip().split("\t")[2])
    if posterior > 0.90:
        circ_output.write(line.strip()+"\t")
        circ_output.write(AnomalyandDecoyDict[junction][1]+"\n")
circInputFile.close()
circ_output.close()


## criteria for reads from linear GLM report to be added:
##   1. the two exons of the linear report are from different genes
##   2. posterior probability > 0.9

linear_output=open(args.FJDir+"reports/AppendedReports/"+args.stem+"_linearJunc_reports.txt", mode="w")
linearInputFile= open( glob.glob(glmDir+"*"+args.stem+"*linearJuncProbs.txt")[0], mode="rU")
for line in linearInputFile:
    if line[0:8]=="junction":
        linear_output.write(line.strip()+"\tanomaly\tnoIndel:Indel1\n")
        continue
    junction=line.strip().split("\t")[0]
    numReads=line.strip().split("\t")[1]
    posterior=float(line.strip().split("\t")[2])
    leftExon=junction.replace(":","|").split("|")[1]
    rightExon=junction.replace(":","|").split("|")[3]
    if leftExon!=rightExon or posterior > 0.90:
        linear_output.write(line.strip()+"\t")
        linear_output.write(AnomalyandDecoyDict[junction][0]+"\t")
        if junction not in RegIndels_1:
            RegIndels_1[junction]=0
        linear_output.write(numReads+":"+str(RegIndels_1[junction])+"\n")


linear_output.close()
linearInputFile.close()