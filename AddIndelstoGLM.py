# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:48:00 2016

@author: Gillian
"""


import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--circleDir", required=True, help= "directory to circular RNA pipeline, subdirectories include orig and circReads")
parser.add_argument("-s", "--stem", required=True, help = "file stem")
parser.add_argument("-w", "--window", required=True, help = "# bases needed on each side of indel to call it an alignment" )
args=parser.parse_args()

if args.circleDir[-1] != "/":
    args.circleDir+= "/"


glmReportsfile = glob.glob(args.circleDir+"circReads/glmReports/*" + args.stem + "*linear*")
indelfiles = glob.glob(args.circleDir+"orig/RegIndelAlignments/" + args.stem + "/*.sam")
reportsDir = args.circleDir + "circReads/reports/"
outputfile= args.circleDir+"circReads/glmReports/AppendGLM/" + args.stem + "_linear_report.txt"

junctionDict={}
indels_1Dict={} # [ins,del]
indels_2Dict={} # [ins, del]

#cycle thorugh glmReports file and identify all junctions in glmReports
f1 = open( glmReportsfile[0], mode = "rU")
for line_raw in f1:
    if line_raw[0:3]!="chr": 
        continue
        
    line=line_raw.strip().split("\t")
    if line[0] not in junctionDict:
        junctionDict[line[0]]=[0]*5 # will be populated with reports tab items
f1.close()
        
# identify the report file
for entry in glob.glob(reportsDir + "*" + args.stem + "*"):
    if "denovo" not in entry:
        reportsfile = entry

# loop thru report file.  If junction is in report file AND glm file, 
#then append the linear/anom/unmap/multimap/circle onto the glm report.

f2 = open(reportsfile, mode = "rU")
for line_raw in f2:
    if line_raw[0:3]!="chr": 
        continue
    line = line_raw.strip().split("\t")
    
    if line[0] in junctionDict:
        junctionDict[line[0]][0]+=int(line[1]) #linear 
        junctionDict[line[0]][1]+=int(line[2]) # anomaly
        junctionDict[line[0]][2]+=int(line[3]) # unampped
        junctionDict[line[0]][3]+=int(line[4]) # multimapped
        junctionDict[line[0]][4]+=int(line[5]) # circle
f2.close()

#print "indelfiles:"
#print sorted(indelfiles)
       

for i in range(0,len(indelfiles)/2):
    f3 = open(sorted(indelfiles)[i], mode = "rU")
    
    for line_raw in f3:
        if line_raw[0]=="@":
            continue
        line=line_raw.strip().split("\t")
        junction=line[2][:-5]
        
        SeqLen = len(line[9])
#        print SeqLen
        SeqStart = int(line[3])     
#        print SeqStart
        SeqEnd = SeqStart+SeqLen
#        print SeqEnd
        NumIndels = int(line[2][-1])
#        print NumIndels
        if line[2][-4:-1]=="INS":
            IndexSeq=300+NumIndels
        if line[2][-4:-1]=="DEL":
            IndexSeq=300-NumIndels
            
#        print IndexSeq
# ONLY COUNT the indel if it overlies the junction by user defined window.            
        if SeqStart < (IndexSeq/2- int(args.window)) and SeqEnd > (IndexSeq/2+ int(args.window)):
            if line[2][:-5] not in indels_1Dict:
                indels_1Dict[junction]=[0]*2
            if line[2][-4:-1]=="INS":
                indels_1Dict[junction][0]+=1
            if line[2][-4:-1]=="DEL":
                indels_1Dict[junction][1]+=1
                    
    f3.close()
    
for i in range (len(indelfiles)/2, len(indelfiles)):
    f4 = open(sorted(indelfiles)[i], mode = "rU")

    for line_raw in f4:
        if line_raw[0]=="@":
            continue
        line=line_raw.strip().split("\t")
        junction=line[2][:-5]
        if line[2][:-5] not in indels_2Dict:
            indels_2Dict[junction]=[0]*2
        if line[2][-4:-1]=="INS":
            indels_2Dict[junction][0]+=1
        if line[2][-4:-1]=="DEL":
            indels_2Dict[junction][1]+=1                
    f4.close()
    

fout= open(outputfile, mode = "w")
fout.write("junction\tnumReads\tp_predicted\tp_value\tlinear\tanomaly\tunmapped\tmultimapped\tcirc\t_1 NoIndels:Indels\t_2 NoIndels:Indels\n")

f1 = open( glmReportsfile[0], mode = "rU")
for line_raw in f1:
    if line_raw[0:3]!="chr":
        continue
    junction = line_raw.strip().split("\t")[0]
    fout.write(line_raw.strip() + "\t")

    for i in range(0,5):
        fout.write(str(junctionDict[junction][i])+"\t")
    
    NumLinear = str(junctionDict[junction][0])
    
    if junction not in indels_1Dict:
        indels_1Dict[junction]=[0]*2
    if junction not in indels_2Dict:
        indels_2Dict[junction]=[0]*2       

    NumIndels_1 = str(indels_1Dict[junction][0]+ indels_1Dict[junction][1])
    NumIndels_2 = str(indels_2Dict[junction][0]+ indels_2Dict[junction][1])
    
    fout.write(NumLinear+":"+NumIndels_1+"\t"+NumLinear+":"+NumIndels_2+"\n")

f1.close()
fout.close()


