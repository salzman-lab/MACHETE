# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:59:42 2016

@author: Gillian
"""

import argparse
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--origDir", required=True, help = "orig directory")
parser.add_argument("-f", "--FJDir", required=True, help = "Far Junctions Directory")
args = parser.parse_args()


if args.origDir[-1] != "/":
    args.origDir += "/"
    
if args.FJDir[-1] != "/":
    args.FJDir += "/"
    
# change the input path to the path where your file exists
os.chdir(args.origDir)
FileList={}
FileStem={}

#create file list of paired files
for name in glob.glob(args.origDir+"genome/*.sam"):
    if "sorted" not in name:    
        (inDir,filename)=os.path.split(name)
        if filename not in FileList:
            FileList[filename] =0
        FileList[filename]+=1
#print FileList

# ran into problems with that every filename had a different unique stem, sometimes occurring
# SAMESAME_SAME_SAME_DIFFERENT_SAME_SAME.sam and if you split this string by _ then
# the "different" part occurs at a variable index each time.
i=1
FoutStems=open(args.FJDir + "/StemList.txt", mode="w")
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

else:
    stem=FileList.keys()[0].split("_")[0]
    FileStem[stem]=0

for stem in FileStem:
    FoutStems.write(stem+"\n")
    
    
FoutStems.close()

#print FileStem
