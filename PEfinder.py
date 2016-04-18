# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:35:23 2015

@author: Gillian
"""

import os
import glob
import argparse
import sys
import utils_os

class ReadInfo:
    def __init__ (self,ID,Q,MAPQ,C,X,side):
        self.id = ID
        self.chromosome = C        
        if side == 0:
            self.StartPosition = max(1,int(X)-int(args.window))
            self.StopPosition = int(X)+int(args.window)
        if side == "R":
            self.StartPosition = int(X)
            self.StopPosition = int(X)+int(args.window)
        if side =="L":
            self.StartPosition = max(1,int(X)-int(args.window))
            self.StopPosition = int(X)
                    
        if Q in ["0","16"] and int(MAPQ) >= 10 and len(C)<3 and C!="M":
#        if Q in ["0","16"] and len(C)<3 and C!="M":  #removed MAPQ as a quality restriction
            self.quality=1
        else:
            self.quality=0

def ID(string):
    if string[-2:] == "/1" or string[-2:] == "/2":
        return string[:-2]
    else:
        return string


def comparereads(read1, read2):
    output = str(read1.id)+"\t"+str(read1.chromosome)+":"+str(read1.StartPosition)+"-"+str(read1.StopPosition)+"\t"+str(read2.chromosome)+":"+str(read2.StartPosition)+"-"+str(read2.StopPosition)+"\n"
            
    if read1.quality == 1 and read2.quality ==1:
        if read1.chromosome != read2.chromosome:
            fout.write(output)
        elif read2.StartPosition < read1.StartPosition-UserBPdistance and read2.StopPosition<read1.StartPosition-UserBPdistance:
            fout.write(output)
        elif read2.StartPosition > read1.StopPosition+UserBPdistance and read2.StopPosition>read1.StopPosition+UserBPdistance:
            fout.write(output)


##  RUN MODE
#parser=argparse.ArgumentParser()
#parser.add_argument("-i","--inputDir", required = True, help = "path to orig directory containing genome & reg files")
#parser.add_argument("-o", "--outputDir", required = True, help = "path to output directory")
#parser.add_argument("-w", "--window", required=True, help = "size = w of window where if read occurs at X, then window starts at X-w and ends at X+w")
#parser.add_argument("-n","--UserBPdistance", required = True, help = "looking for PE > X base pairs apart. Linda's default window is 100K.")
#args=parser.parse_args()
#
#
#
#if args.inputDir[-1] != "/":
#    inpath = args.inputDir + "/"
#else:
#    inpath = args.inputDir
#
#if args.outputDir[-1] != "/":
#    outpath = args.outputDir + "/"
#else:
#    outpath = args.outputDir
#    
#utils_os.createDirectory(outpath)
#
#UserBPdistance = int(args.UserBPdistance)
#
#print "Analyzing" + str(args.inputDir)
#
#
##TESTING MODE
##inpath = "/Users/Gillian/Desktop/transcriptometest/"
##outpath = "/Users/Gillian/Desktop/transcriptometest/"
##UserBPdistance = 100000
#
## change the input path to the path where your file exists
#os.chdir(inpath)
#FileList={}
#FileStem={}
#
##create file list of paired files
#for name in glob.glob(os.path.join(inpath,"genome/sorted*.sam")):
#    (inpath1,filename)=os.path.split(name)
#    if not filename in FileList:
#        FileList[filename] =0
#    FileList[filename]+=1
#print FileList
#
## ran into problems with that every filename had a different unique stem, sometimes occurring
## SAMESAME_SAME_SAME_DIFFERENT_SAME_SAME.sam and if you split this string by _ then
## the "different" part occurs at a variable index each time.
#i=1
#if len(FileList)>2:
#    for i in range(0,len(FileList.keys()[0].split("_"))):
#        if FileList.keys()[0].split("_")[i] in ["1","2"]:
#            continue
#        if FileList.keys()[0].split("_")[i] != FileList.keys()[1].split("_")[i] or FileList.keys()[0].split("_")[i] != FileList.keys()[2].split("_")[i]:
#            break
#
#for name in FileList:
#    stem = name.split("_")[i]
#    if not stem in FileStem:
#        FileStem[stem] = 0
#    FileStem[stem] +=1
#print FileStem
#
#
## for each pair of read1 and read2 files -- 
#counter =0
#for stem in FileStem:
#    SamFiles = []
#    for name in glob.glob(os.path.join(inpath+"genome/", "sorted*" + stem + "*.sam")):
#        SamFiles.append(name)
#    for name in glob.glob(os.path.join(inpath+"reg/", "sorted*"+ stem+"*.sam")):
#        SamFiles.append(name)
#
#    print SamFiles
   
   
parser=argparse.ArgumentParser()
parser.add_argument("-o", "--origDir", required=True, help = "orig directory")
parser.add_argument("-s","--stem", required = True, help = "stem name for file comparison")
parser.add_argument("-f", "--outputDir", required = True, help = "Far Junctions directory (output dir)")
parser.add_argument("-w", "--window", required=True, help = "size = w of window where if read occurs at X, then window starts at X-w and ends at X+w")
parser.add_argument("-n","--UserBPdistance", required = True, help = "looking for PE > X base pairs apart. Linda's default window is 100K.")
args=parser.parse_args()


if args.origDir[-1] != "/":
    inpath = args.origDir + "/"
else:
    inpath = args.origDir

if args.outputDir[-1] != "/":
    outpath = args.outputDir + "/DistantPEFiles/"
else:
    outpath = args.outputDir + "DistantPEFiles/"
    
utils_os.createDirectory(outpath)

UserBPdistance = int(args.UserBPdistance)



# change the input path to the path where your file exists
os.chdir(inpath)


genomefiles = sorted(glob.glob(inpath + "genome/sorted*" + args.stem + "*.sam"))
regfiles = sorted(glob.glob(inpath + "reg/sorted*" + args.stem + "*.sam"))

print genomefiles
print regfiles



# print SamFiles
#opening paired files
g1 = open(genomefiles[0], mode ="rU")
g2 = open(genomefiles[1], mode="rU")
r1 = open(regfiles[0], mode ="rU")
r2 = open(regfiles[1], mode = "rU")
            
#creating output files
fout = open(outpath + args.stem +"_distant_pairs.txt", mode="w")
print"writing to output file: " + args.stem +"_distant_pairs.txt"
fout.write("Distant pairs more than "+str(UserBPdistance)+"BP apart in the genome: \n\nID\tposition1\tposition2\n")
sys.stdout.flush() #force print command above to output to screen (gives user an idea of progress)

# get rid of headers

line_raw2 = g2.next()
line_raw4 = r2.next()

while line_raw2[0] == "@":
    line_raw2=g2.next()
while line_raw4[0] == "@":
    line_raw4=r2.next()
    
lineG2 = line_raw2.strip().split("\t")
lineR2 = line_raw4.strip().split("\t")
# looping through reads in genome file 1



for line_raw1 in g1:
#        print "in genome 1"
    if line_raw1[0]=="@":
        continue
    lineG1 = line_raw1.strip().split("\t")

#comparing to genome file 2 

    while ID(lineG2[0])<ID(lineG1[0]) and lineG2!=["FILE_END"]:
#            print "genome 2 < genome 1"
        try: lineG2 = g2.next().strip().split("\t")
        except: lineG2 = ["FILE_END"]
    
    if ID(lineG1[0])==ID(lineG2[0]):
#            print "genome 1 and genome 2 match"
        readG1 = ReadInfo(ID(lineG1[0]), lineG1[1], lineG1[4], lineG1[2][3:], lineG1[3], 0)
        readG2 = ReadInfo(ID(lineG1[0]), lineG2[1], lineG2[4], lineG2[2][3:], lineG2[3], 0)
        
#            fout.write("genome1 - genome 2:\n")
        comparereads(readG1, readG2)            
        
        try: lineG2 = g2.next().strip().split("\t")
        except: lineG2 = ["FILE_END"]

    while ID(lineR2[0])<ID(lineG1[0]) and lineR2!=["FILE_END"]:
#            print "reg2 < genome 1"
        try: lineR2 = r2.next().strip().split("\t")
        except: lineR2 = ["FILE_END"]
        
        
    if ID(lineG1[0])==ID(lineR2[0]):
        readG1 = ReadInfo(ID(lineG1[0]), lineG1[1], lineG1[4], lineG1[2][3:], lineG1[3], 0)
        R2loc = lineR2[2].replace("|", " ").replace(":"," ").split(" ")
        if int(R2loc[2])>=int(R2loc[4]):
            read1R2 = ReadInfo(ID(lineR2[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[2], "R")
            read2R2 = ReadInfo(ID(lineR2[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[4], "L")
        else:
            read1R2 = ReadInfo(ID(lineR2[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[2], "L")
            read2R2 = ReadInfo(ID(lineR2[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[4], "R")                
        
#            fout.write("genome1 - reg 2: \n")
        comparereads(readG1, read1R2)
        comparereads(readG1, read2R2)            
        
        try: lineR2 = r2.next().strip().split("\t")
        except: lineR2 = ["FILE_END"]
r2.close()
g2.close()

#===========LOOPING THROUGH REG 1 =============

g2 = open(genomefiles[1], mode="rU")
r2 = open(regfiles[1], mode = "rU")

lineG2 = g2.next().strip().split("\t")
lineR2 = r2.next().strip().split("\t")

while lineG2[0][0] == "@":
    lineG2=g2.next().strip().split("\t")
while lineR2[0][0] == "@":
    lineR2=r2.next().strip().split("\t")
    
for line_raw3 in r1:
   
    if line_raw3[0]=="@":
        continue
    
    lineR1 = line_raw3.strip().split("\t")

    while ID(lineG2[0])<ID(lineR1[0]) and lineG2 != ["FILE_END"]:
 #           print lineG2[0][:-3] + " < " + lineR1[0][:-3]
#            print "genome2 read < reg 1 read"
        try: lineG2 = g2.next().strip().split("\t")
        except: lineG2 = ["FILE_END"]

    
    if ID(lineR1[0])==ID(lineG2[0]):
#            print "reg1 and genome2 match"
        R1loc = lineR1[2].replace("|", " ").replace(":"," ").split(" ")
        #print R1loc
        
        if int(R1loc[2])>=int(R1loc[4]):
            read1R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[2], "R")
            read2R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[4], "L")
        else:
            read1R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[2], "L")
            read2R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[4], "R")                

        readG2 = ReadInfo(ID(lineG2[0]), lineG2[1], lineG2[4], lineG2[2][3:], lineG2[3], 0)
#            fout.write("reg1 - genome 2: \n")            
        comparereads(read1R1, readG2)
        comparereads(read2R1, readG2)
        
        try: lineG2 = g2.next().strip().split("\t")
        except: lineG2 = ["FILE_END"]
    
    while ID(lineR2[0])<ID(lineR1[0]) and lineR2!=["FILE_END"]:
#            print "reg2 read < reg 1 read"
        try: lineR2 = r2.next().strip().split("\t")
        except: lineR2 = ["FILE_END"]
        
    if ID(lineR1[0])==ID(lineR2[0]):
        R1loc = lineR1[2].replace("|", " ").replace(":"," ").split(" ")
 #           print "reg 1 and reg 2 match"
        if int(R1loc[2])>=int(R1loc[4]):
            read1R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[2], "R")
            read2R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[4], "L")
        else:
            read1R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[2], "L")
            read2R1 = ReadInfo(ID(lineR1[0]), lineR1[1], lineR1[4], R1loc[0][3:], R1loc[4], "R")                
 
        
        R2loc = lineR2[2].replace("|", " ").replace(":"," ").split(" ")
        if int(R2loc[2])>=int(R2loc[4]):
            read1R2 = ReadInfo(ID(lineR1[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[2], "R")
            read2R2 = ReadInfo(ID(lineR1[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[4], "L")
        else:
            read1R2 = ReadInfo(ID(lineR1[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[2], "L")
            read2R2 = ReadInfo(ID(lineR1[0]), lineR2[1], lineR2[4], R2loc[0][3:], R2loc[4], "R")                
 
#            fout.write("reg1 - reg 2: \n")             
        comparereads(read1R1, read1R2)
        comparereads(read1R1, read2R2)
        comparereads(read2R1, read1R2)
        comparereads(read2R1, read2R2)            
        
        try: lineR2 = r2.next().strip().split("\t")
        except: lineR2= ["FILE_END"]
  
  
         
fout.close()            
g1.close()
g2.close()
r1.close()
r2.close()

   