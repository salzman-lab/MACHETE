# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:51:44 2015

@author: Gillian
"""
import argparse
import os
import glob


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--stem", required=True, help = "file stem name")
parser.add_argument("-p", "--parentDir", required=True, help = "path to parent directory of input and output directories; FJ directory")
parser.add_argument("-n", "--overlap", required=True, help = "required overlapping nt on each side of junction, including # indels")
parser.add_argument("-x", "--NumIndels", required=True, help = "number of indels allowed from previous file")

args = parser.parse_args()


if args.parentDir[-1] != "/": 
    args.inputDir += "/"
    
## WITH USE ON FAR JUNCTIONS
FarJunctionAlignments = args.parentDir + "FarJunctionAlignments/" + args.stem + "/"
AlignedIndels = args.parentDir+"FarJuncSecondary/AlignedIndels/"+args.stem + "/"
RemoveNonOverlap = args.parentDir+"FarJuncSecondary/AlignedIndels/RemoveNonOverlap/"+args.stem + "/"
IndelsHistogram = args.parentDir + "IndelsHistogram/"


## TEMPORARY, FOR USE ALIGNING A RANDOM SET OF FILES, CHANGE AS NECESSARY
#FarJunctionAlignments= args.parentDir
#AlignedIndels=args.parentDir+"AlignToIndels/"
#RemoveNonOverlap=AlignedIndels+"RemoveNonOverlap/"
#IndelsHistogram=args.parentDir+"IndelsHistogram"


for filename in glob.glob(AlignedIndels+"*.sam"):
    (inpath,infile)=os.path.split(filename)

    f1 = open(filename, mode ="rU")
    fout = open(RemoveNonOverlap + infile, mode ="w")
    print "opened file: " + filename
    
    
    for line in f1:
        if line[0] == "@":
            fout.write(line)
            continue
                
        line_list = line.strip().split("\t")
        if int(line_list[3]) > (150-int(args.overlap)) or (int(line_list[3])+len(line_list[9])) <  (150+int(args.overlap)):
            #print "range = " + line_list[3] + "-" + str(int(line_list[3])+len(line_list[9]))
            pass
        else:
            #print "overlap contained"
            fout.write(line)

        
    f1.close()
    fout.close()
    

#    ==============Now nonoverlapping are removed. time to make histogram. ==================


# open  file and feed all junctions into a junction dictionary

for filename in glob.glob(FarJunctionAlignments+ "*.sam"):
    (path, infile) = os.path.split(filename)
    stem = infile[:-4] ##stem includes file stem and the _1 or _2 info.
   
    f1 = open(filename, mode="rU")
    JunctionDict = {}
    print "opening " + args.stem
    for line_raw in f1:
        if line_raw[0]=="@":
            continue
        
        line=line_raw.strip().split("\t")

        if int(line[3]) > (150-int(args.overlap)) or (int(line[3])+len(line[9])) <  (150+int(args.overlap)):

            if line[2] not in JunctionDict:
                JunctionDict[line[2]]= [0]*(2*int(args.NumIndels)+1)
        
            JunctionDict[line[2]][int(args.NumIndels)]+=1
        
    f1.close()


# sequentially open all indel files and compare them to the original juction dictionary    
    for file2 in glob.glob(RemoveNonOverlap+"*"):
        if stem not in file2:
            continue
        print file2
        f2 = open(file2, mode = "rU")
        for line_raw in f2:
            if line_raw[0] == "@":
                continue
            
        
            line=line_raw.strip().split("\t")
            position = int(line[2][-1])
            if line[2][-4:-1]=="DEL":
                position = position * -1
                
                # populates dictionary "JunctionDict" where each 
                # key is junction name
                # value is an array representing a histogram of the number of times
                # each position appeared
                # [x deletions, x+1 deletions....0 indels..., x insertions]
            
            if line[2][:-5] in JunctionDict:
                JunctionDict[line[2][:-5]][position+int(args.NumIndels)]+=1
                        
        f2.close()
    
    # print results 
    fout = open(IndelsHistogram+stem+".txt", mode = "w")
        
    for key in JunctionDict:
        fout.write(key+"\t"+str(JunctionDict[key])+"\n")
        
    fout.close()
    
        
