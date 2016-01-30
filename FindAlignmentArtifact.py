# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:51:44 2015

@author: Gillian
"""
import argparse
import os
import glob


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--parentDir", required=True, help = "path to parent directory of input and output directories")
parser.add_argument("-n", "--overlap", required=True, help = "required overlapping nt on each side of junction")
parser.add_argument("-x", "--NumIndels", required=True, help = "number of indels allowed from previous file")

args = parser.parse_args()


if args.parentDir[-1] != "/": 
    args.inputDir += "/"
    
## WITH USE ON FAR JUNCTIONS
FarJunctionAlignments = args.parentDir + "FarJunctionAlignments/"
AlignedIndels = args.parentDir+"FarJuncSecondary/AlignedIndels/"
RemoveNonOverlap = AlignedIndels+"RemoveNonOverlap/"
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
    

#    ==============Now nonoverlapping are removed. time to make histogram. ==================    


for filename in glob.glob(FarJunctionAlignments+"*.sam"):
    # for each junction, there is going to be a string with len 2* NumIndels + 1 
    # each bin from left to right will represent x deletions....0....x insertion
    # parse thru FarJunction file to grab all the junctions
    # then parse thru AlignedIndels file to grab junctions from the AlignedIndels file with that same 
    # file name on it.
  #  *** args.NumIndels  
    (path, infile) = os.path.split(filename)
    stem = infile[:-4]


# open  file and feed all junctions into a junction dictionary  
    f1 = open(filename, mode="rU")
    JunctionDict = {}
    print "opening " + stem
    for line_raw in f1:
        while line_raw[0]=="@":
            line_raw=f1.next()
        
        line=line_raw.strip().split("\t")
        
        if line[2] not in JunctionDict:
            JunctionDict[line[2]]= [0]*(2*int(args.NumIndels)+1)
        
        JunctionDict[line[2]][int(args.NumIndels)]+=1
        
    f1.close()


# sequentially open all indel files and compare them to the original juction dictionary    
    for file2 in glob.glob(RemoveNonOverlap+"*"+stem+"*"):
        print file2
        f2 = open(file2, mode = "rU")
        for line_raw in f2:
            if line_raw[0] == "@":
                line_raw=f2.next()
        
    
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
    