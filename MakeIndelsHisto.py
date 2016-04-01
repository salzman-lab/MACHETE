# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 15:25:42 2016

@author: Gillian
"""

import argparse
import glob


def ID(string):
    if string[-2:] == "/1" or string[-2:] == "/2":
        return string[:-2]
    else:
        return string


class ReadInfoFJ:
    def __init__ (self,line_raw):
        line = line_raw.strip().split("\t")
        self.ID = ID(line[0])       
        self.refstrand = line[1]
        self.junction = line[2]
        if self.junction[-4:-1]=="DEL":
            self.indel=-int(self.junction[-1:])
        elif self.junction[-4:-1]=="INS":
            self.indel=int(self.junction[-1:])
        self.MAPQ = int(line[4])
        self.AS= int(line[11].split(":")[2])
        self.NumOfBases = len(line[9])
        self.offset = int(line[3])
        if "XS:i:" in line[12]:
            self.NumN=line[13][5:]
        else:
            self.NumN=line[12][5:]


        JuncInfo = line[2].replace(":"," ").replace("|"," ").split(" ")
        self.chr_left=JuncInfo[0]
        self.loc_left=JuncInfo[2]
        self.strand_left = JuncInfo[3]
        self.chr_right = JuncInfo[4]
        self.loc_right = JuncInfo[6]
        self.strand_right = JuncInfo[7]



parser = argparse.ArgumentParser()
parser.add_argument("-s", "--stem", required=True, help = "file stem name")
parser.add_argument("-f", "--FJDir", required=True, help = "FJ directory")
parser.add_argument("-w", "--overlap", required=True, help = "required overlapping nt on each side of junction, including # indels")
parser.add_argument("-x", "--NumIndels", required=True, help = "number of indels allowed from previous file")

args = parser.parse_args()


if args.FJDir[-1] != "/": 
    args.FJDir += "/"
    
#
### make a concatenated list of all the far junction alignment files

#
#input_lines=fileinput.input(FJ1_list)
#fout_FJ1.writelines(input_lines)
#input_lines=fileinput.input(FJ2_list)
#fout_FJ2.writelines(input_lines)
#
#fout_FJ1.close()
#fout_FJ2.close()


## now go through all primary alignments and make a dictionary of possible far juncs:

FJDict_1={} ## key = junc name, Val = # [0,0,0,0,X,0,0,0,0]  0's # of indels
FJDict_2={}

FJFile1=sorted(glob.glob(args.FJDir+"FarJunctionAlignments/"+args.stem+"/*.sam"))[0]
FJFile2=sorted(glob.glob(args.FJDir+"FarJunctionAlignments/"+args.stem+"/*.sam"))[1]

f1= open(FJFile1, mode="rU")
print "opening" + FJFile1

linecount=0
goodlinecount=0
newjunc=0
for line_raw in f1:
    if line_raw[0]=="@":
        continue
    linecount+=1
    FJread=ReadInfoFJ(line_raw)
    if FJread.offset <= (150-int(args.overlap)) and FJread.offset+FJread.NumOfBases >= (150+int(args.overlap)):
        if FJread.junction not in FJDict_1:
            FJDict_1[FJread.junction]=[0]*(2*int(args.NumIndels)+1)
            newjunc+=1
        FJDict_1[FJread.junction][int(args.NumIndels)]+=1
        goodlinecount+=1
f1.close()
print "new junc added" + str(newjunc)
print linecount
print goodlinecount


linecount=0
goodlinecount=0
newjunc=0

f2= open(FJFile2, mode="rU")
print "opening" + FJFile2
for line_raw in f2:
    if line_raw[0]=="@":
        continue
    FJread=ReadInfoFJ(line_raw)
    linecount+=1
    if FJread.offset <= (150-int(args.overlap)) and FJread.offset+FJread.NumOfBases >= (150+int(args.overlap)):
        if FJread.junction not in FJDict_2:
            FJDict_2[FJread.junction]=[0]*(2*int(args.NumIndels)+1)
            newjunc+=1

        FJDict_2[FJread.junction][int(args.NumIndels)]+=1
        goodlinecount+=1
f2.close()
print "new junc added" + str(newjunc)
print linecount
print goodlinecount

## make a dictionary of readIDs that aligned to any FJ indel



AlignedFJFiles=[]

for name in glob.glob(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ args.stem + "/*.sam"):
    if "All_" not in name:
        AlignedFJFiles.append(name)

FJ1_list= sorted(AlignedFJFiles)[0:len(AlignedFJFiles)/2]
FJ2_list= sorted(AlignedFJFiles)[len(AlignedFJFiles)/2:]
#
#print "FJ1 indels list"
#print FJ1_list
#print "FJ2 indels list"
#print FJ2_list

IndelsReadIDs={}

for name in FJ1_list:
    print "FJ1 indels"
    print name
    f1= open(name, mode="rU")
    
    for line in f1:
        if line[0]=="@":
            continue
        
        read = ReadInfoFJ(line)
        # if the read overlaps the junction
        if read.offset <= (150-int(args.overlap)+read.indel) and read.offset+read.NumOfBases >= (150+int(args.overlap)+read.indel):

            # if the read isn't in dictionary then add it 
            if read.ID not in IndelsReadIDs:             
                IndelsReadIDs[read.ID]= line
            # if the read is in the dictionary, compare it to existing read. If AS is better, then replace existing read
            else:
                compareRead= ReadInfoFJ(IndelsReadIDs[read.ID])
                if int(compareRead.AS) >= int(read.AS):
                    pass
                else:
                    IndelsReadIDs[read.ID]= line                            
    f1.close()
  
# write all distinct readIDs to an All_1_indels file  
fout_FJ1= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ args.stem + "/All_" + args.stem + "_1_indels.sam", mode="w")

for key in IndelsReadIDs:
    fout_FJ1.write(IndelsReadIDs[key].strip()+"\n")
fout_FJ1.close()
    

## CLEAR Read IDs dictionary and do the same with FJ2 list
IndelsReadIDs={}

for name in FJ2_list:
    print "FJ2 indels"
    print name    
    
    f1= open(name, mode="rU")
    
    for line in f1:
        if line[0]=="@":
            continue
        
        read = ReadInfoFJ(line)
        # if the read overlaps the junction
        if read.offset <= (150-int(args.overlap)+read.indel) and read.offset+read.NumOfBases >= (150+int(args.overlap)+read.indel):

            # if the read isn't in dictionary then add it 
            if read.ID not in IndelsReadIDs:             
                IndelsReadIDs[read.ID]= line
            # if the read is in the dictionary, compare it to existing read. If AS is better, then replace existing read
            else:
                compareRead= ReadInfoFJ(IndelsReadIDs[read.ID])
                if int(compareRead.AS) >= int(read.AS):
                    pass
                else:
                    IndelsReadIDs[read.ID]= line                            
    f1.close()
  

# write all distinct readIDs to an All_1_indels file  
fout_FJ2= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ args.stem + "/All_" + args.stem + "_2_indels.sam", mode="w")

for key in IndelsReadIDs:
    fout_FJ2.write(IndelsReadIDs[key].strip()+"\n")
fout_FJ2.close()

## parse AllIndels_1 and AllIndels_2 files to see if they aligned to the same juncs with indels as an FJ
## if yes, then add to junction "dictionary"

Indels1=open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ args.stem + "/All_" + args.stem + "_1_indels.sam", mode="rU")
for line in Indels1:
    read=ReadInfoFJ(line) 
    if read.junction[:-5] in FJDict_1:
        FJDict_1[read.junction[:-5]][int(args.NumIndels) + read.indel]+=1
Indels1.close()

Indels2=open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ args.stem + "/All_" + args.stem + "_2_indels.sam", mode="rU")
for line in Indels2:
    read=ReadInfoFJ(line) 
    if read.junction[:-5] in FJDict_2:
        FJDict_2[read.junction[:-5]][int(args.NumIndels) + read.indel]+=1
Indels2.close()

## output indels histo

Outfile1 = open(args.FJDir + "IndelsHistogram/indels_" + args.stem + "_1.txt", mode="w")
for key in FJDict_1:
  Outfile1.write(key+"\t"+str(FJDict_1[key])+"\n")
Outfile1.close()


Outfile2 = open(args.FJDir + "IndelsHistogram/indels_" + args.stem + "_2.txt", mode="w")
for key in FJDict_2:
  Outfile2.write(key+"\t"+str(FJDict_2[key])+"\n")
Outfile2.close()
