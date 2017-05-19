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

FoutStems=open(args.FJDir + "/StemList.txt", mode="w")

FileList=[]
FileListPairs={}
FileListEntries=[]
#
##create file list of paired files
#for name in glob.glob(args.origDir+"genome/*.sam"):
#    if "sorted" not in name:    
#        (inDir,filename)=os.path.split(name)
#        if filename not in FileList:
#            FileList[filename] =0
#        FileList[filename]+=1
##print FileList

for name in glob.glob(args.origDir+"genome/*.sam"):
    if "sorted" not in name:
        (inDir,filename)=os.path.split(name)
        FileList.append(filename)

#print FileList
#print len(FileList[0].split("_"))-3

Counter=len(FileList)
#print Counter

while Counter!=0:
    firstEntry=FileList[0].split("_")[:-3]
    for i in range(1,Counter):
        compareEntry=FileList[i].split("_")[:-3]      
        if firstEntry==compareEntry:
            Stem="_".join(firstEntry)
#            print Stem
            FoutStems.write(Stem+"\n")

            del FileList[0]
            del FileList[i-1]
            Counter = Counter-2
            break
            
    
FoutStems.close()
            

