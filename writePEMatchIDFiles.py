# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 12:51:03 2016

@author: Gillian
"""
import os
import glob
import argparse


parser=argparse.ArgumentParser()
parser.add_argument("-o", "--origDir", required=True, help = "path to orig dir" )
parser.add_argument("-f", "--FJDir", required=True, help = "path to FJ dir")
args=parser.parse_args()

if args.origDir[-1] != "/":
    args.origDir += "/"

genomeDir = args.origDir+"genome/"
regDir= args.origDir+"reg/"

fout = open(args.FJDir+"PE_ID_files.txt", mode="w")

FileList={}
FileStem={}

for name in glob.glob(os.path.join(genomeDir+ "*sorted*.sam")):
#for name in glob.glob(os.path.join(genomeDir+ "*.sam")):
    (genomepath,filename)=os.path.split(name)
    if not filename in FileList:
        FileList[filename] =0
    FileList[filename]+=1
print FileList

i=1
if len(FileList)>2:
    for i in range(0,len(FileList.keys()[0].split("_"))):
        if FileList.keys()[0].split("_")[i] in ["1","2"]:
            continue
        if FileList.keys()[0].split("_")[i] != FileList.keys()[1].split("_")[i] or FileList.keys()[0].split("_")[i] != FileList.keys()[2].split("_")[i]:
            break

for name in FileList:
    stem = name.split("_")[i]
    if not stem in FileStem:
        FileStem[stem] = 0
    FileStem[stem] +=1
print FileStem

counter =0

for stem in FileStem:
    SamFiles = []
    for name in glob.glob(os.path.join(genomeDir+"sorted*"+stem+"*.sam")):
#    for name in glob.glob(os.path.join(genomeDir+"*"+stem+"*.sam")):
        SamFiles.append(name)
    for name in glob.glob(os.path.join(regDir+"sorted*"+stem+"*.sam")):
#    for name in glob.glob(os.path.join(regDir+"*"+stem+"*.sam")):
        SamFiles.append(name)
    
    sortedSamFiles = sorted(SamFiles)
    
    print sortedSamFiles
    fout.write(sortedSamFiles[0]+"\t"+sortedSamFiles[1]+"\t"+sortedSamFiles[2]+"\t"+sortedSamFiles[3]+"\n")    

fout.close()