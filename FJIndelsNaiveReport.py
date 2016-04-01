# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 16:12:27 2015

@author: Gillian
"""
# takes alignments from Indels and finds read partner. 
# all Far Junction and Scrambled junction reads where the alignment does not overlap the 
# Junction by args.window bp will be thrown out.

# This program then tells if read partners "makes sense" or not
# final Output categories -- [0] R1 junc name
            #   [1] genome - R2 location < 100mill bp away
            #   [2] genome anomaly - R2 location > 100mill bp away
            #   [3] genome p value
            #   [4] reg - R2 closest location <100mill bp away
            #   [5] reg anomaly - R2 closest location > 100 mill bp away
            #   [6] reg p value
            #   [7] junc / scrambled - R2 closest location <100mill bp away
            #   [8] junc anomaly - R2 closest location > 100mill bp away
            #   [9] junc p value
            #   [10] FarJunc - R2 aligned to same FarJunc
            #   [11] FarJunc anomaly - R2 aligned to diff FarJunc
            #   [12] FarJunc p value
            #   [13] unaligned - R2 didn't align
            #   [14] unmapped - R2 missing in action
            #   [15] P val for all non-anomaly classes



################
#Current categories
#FJgood -- genome, reg, FJ
#FJbad -- genome anomaly, reg anomaly, junc, junc anomaly, FJ anomaly
#################



import argparse
import os
import glob
import fileinput
from math import sqrt
from scipy.stats import norm
from collections import Counter

def AddToDict(inputtype, TargetDict, line_raw_comparison, line_raw_FJ):

    lineFJ = ReadInfoFJ(line_raw_FJ)
#    if lineFJ.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
#        print "ERROR AT LINE 42"
#        print inputtype
#        print line_raw_comparison
#        print line_raw_FJ
#    

    if lineFJ.junction not in TargetDict:  # add junction to target dictionary if it doesn't exist
        TargetDict[lineFJ.junction] = [0,0,0.0,0.0] 
    
        
    if inputtype=="FJ": # if comparing Far Junc to Far Junc, they have to be identical
        line2= ReadInfoFJ(line_raw_comparison)
        
        ## output R1 -  offset, MAPQ, AS, #N, readlen, junc name, strand
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.AS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.AS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand

        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

        if lineFJ.junction == line2.junction and lineFJ.refstrand in ["0","16"] and line2.refstrand in ["0","16"] and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
            IDfiletype = "FJgood,FarJunction"
        else:
            TargetDict[lineFJ.junction][1]+=1
            IDfiletype = "FJbad,FarJuncAnom"
            addAS = 0.0
            addNumofBases = 0.0
            
        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases

        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype=="reg" or inputtype =="junc": #if reg or junc read, then one side has to be within 100KB, and meets refstrand criteria below
        line2 = ReadInfoJunc(line_raw_comparison)

        ## output R1 -  offset, MAPQ, AS, #N, readlen, junc name, strand
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.AS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.AS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand


        if inputtype =="reg": IDfiletype = "FJgood,Regular"
        if inputtype == "junc": IDfiletype = "FJbad,Junction"

        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

        
        if lineFJ.strand_left==line2.strand and lineFJ.strand_left=="-" and lineFJ.chr_left == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_left), int(lineFJ.loc_left)+100000000) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand: 
            TargetDict[lineFJ.junction][0]+=1 
        elif lineFJ.strand_left==line2.strand and lineFJ.strand_left=="-" and lineFJ.chr_left == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_left), int(lineFJ.loc_left)+100000000) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand: 
            TargetDict[lineFJ.junction][0]+=1  
        elif lineFJ.strand_left==line2.strand and lineFJ.strand_left=="+" and lineFJ.chr_left == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_left)-100000000,int(lineFJ.loc_left)) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_left==line2.strand and lineFJ.strand_left=="+" and lineFJ.chr_left == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_left)-100000000,int(lineFJ.loc_left)) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
                     
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "-" and lineFJ.chr_right == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_right)-100000000,int(lineFJ.loc_right)) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "-" and lineFJ.chr_right == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_right)-100000000,int(lineFJ.loc_right)) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "+" and lineFJ.chr_right == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_right),int(lineFJ.loc_right)+100000000) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "+" and lineFJ.chr_right == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_right),int(lineFJ.loc_right)+100000000) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
            
        else:
            TargetDict[lineFJ.junction][1]+=1
            if inputtype =="reg": IDfiletype = "FJbad,RegAnomaly"
            if inputtype == "junc": IDfiletype = "FJbad,JuncAnomaly"
            addAS = 0.0
            addNumofBases = 0.0


        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
       
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype == "genome": #comparing FJ to genome, has to be within 100Kbp, meet ref strand criteria (opp refstrand if + read, same refstrand if - read)

        line2 = ReadInfoGenome(line_raw_comparison)

        ## output R1 -  offset, MAPQ, AS, #N, readlen, junc name, strand
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.AS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.loc) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.AS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.chr + "\t" + line2.refstrand




        IDfiletype = "FJgood,genome"
        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

 
        if lineFJ.strand_left=="-" and lineFJ.chr_left == line2.chr and int(line2.loc) in range(int(lineFJ.loc_left), int(lineFJ.loc_left)+100000000) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand: 
            TargetDict[lineFJ.junction][0]+=1  
        elif lineFJ.strand_left=="+" and lineFJ.chr_left == line2.chr and int(line2.loc) in range(int(lineFJ.loc_left)-100000000,int(lineFJ.loc_left)) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right == "-" and lineFJ.chr_right == line2.chr and int(line2.loc) in range(int(lineFJ.loc_right)-100000000,int(lineFJ.loc_right)) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right == "+" and lineFJ.chr_right == line2.chr and int(line2.loc) in range(int(lineFJ.loc_right),int(lineFJ.loc_right)+100000000) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        else:
            TargetDict[lineFJ.junction][1] +=1
            IDfiletype = "FJbad,genomAnomaly"
            addAS = 0.0
            addNumofBases = 0.0           
            
        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")

    if inputtype == "unaligned":
        if lineFJ.ID not in FJDict:
            if lineFJ.junction not in TargetDict:
                TargetDict[lineFJ.junction] = 0 
            TargetDict[lineFJ.junction]+=1
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.AS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
        IDfiletype = "unaligned"
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\n")
        
     

    return TargetDict

def ID(string):
    if string[-2:] == "/1" or string[-2:] == "/2":
        return string[:-2]
    else:
        return string


def Pvalue(AS, NumBases):
    ExpectedMMrate = 0.01
    ExpectedMM = ExpectedMMrate*float(NumBases)
    if NumBases == 0.0:
        return "-"
    ActualMM = float(AS)/(-6.0)
    Z = (ActualMM - ExpectedMM ) / sqrt(float(NumBases)*(ExpectedMMrate)*(1.0-ExpectedMMrate))
    return norm.cdf(-Z)


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
        
class ReadInfoGenome:
    def __init__ (self, line_raw):
        line = line_raw.strip().split("\t")
        self.ID= ID(line[0])
        self.refstrand = line[1]
        self.chr = line[2]
        self.loc = line[3]
        self.MAPQ = int(line[4])
        self.AS = int(line[11].split(":")[2])
        self.NumOfBases=len(line[9])
        self.offset = int(line[3])
        if "XS:i:" in line[12]:
            self.NumN=line[13][5:]
        else:
            self.NumN=line[12][5:]

        
class ReadInfoJunc:
    def __init__ (self, line_raw):
        line = line_raw.strip().split("\t")
        self.ID = ID(line[0])
        self.refstrand = line[1]
        self.junction = line[2]
        self.chr = line[2].split("|")[0]
        self.loc_1 = line[2].replace(":","|").split("|")[2]
        self.loc_2 = line[2].replace(":","|").split("|")[4]
        self.strand = line[2][-1]
        self.MAPQ = int(line[4])
        self.AS = int(line[11].split(":")[2])
        self.NumOfBases = len(line[9])      
        self.offset = int(line[3])
        if "XS:i:" in line[12]:
            self.NumN=line[13][5:]
        else:
            self.NumN=line[12][5:]

        
#=========================================
#start here

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--stem", required=True, help = "stem name of file to generate report")
parser.add_argument("-f", "--FJDir", required = True, help = "path to aligned junction reads")
parser.add_argument("-i", "--origDir", required=True, help = "path to orig dir containing genome reads")
parser.add_argument("-w", "--window", required=True, help = "# of bases needed on each side of the junction")
args = parser.parse_args()
window= int(args.window)

# f1 = open("/Users/Gillian/Desktop/sherlock/unaligned_ENCFF000HOC1_1.sam", mode ="rU")
# f2 = open("/Users/Gillian/Desktop/sherlock/20000_ENCFF000HOC2_1_genome_output.sam", mode ="rU")

if args.FJDir[-1]!= "/":
    args.FJDir+="/"
if args.origDir[-1]!="/":
    args.origDir+="/"

stem = args.stem


FarJunctionfiles=[]
genomefiles=[]
regfiles=[]
junctionfiles=[]
unalignedfiles=[]

for name in glob.glob(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/*.sam"):
    print name
    FarJunctionfiles.append(name)
    # FarJunctionFiles contains indel alignments for _1 and _2 files to indels 1-5


for name in glob.glob(os.path.join(args.origDir,"genome/*" + stem + "*.sam")):        
    print name
    if "sorted" not in name:       
        genomefiles.append(name)

for name in glob.glob(os.path.join(args.origDir,"reg/*" + stem + "*.sam")):
    print name
    if "sorted" not in name:       
        regfiles.append(name) 
for name in glob.glob(os.path.join(args.origDir,"junction/*" + stem + "*.sam")):
    print name
    if "sorted" not in name:       
        junctionfiles.append(name) 
for name in glob.glob(os.path.join(args.origDir,"unaligned/*" + stem + "*.fq")):
    print name
    if "sorted" not in name:       
        unalignedfiles.append(name)         
    

# opening all files for a particular stem
print sorted(FarJunctionfiles)
print sorted(genomefiles)
print sorted(regfiles)
print sorted(junctionfiles)
print sorted(unalignedfiles)


## concatenate all fJ indels files into a single big indels file
fout_FJ1= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/All_" + stem + "_1_indels.sam", mode="w")
fout_FJ2= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/All_" + stem + "_2_indels.sam", mode="w")

FJ1_list= sorted(FarJunctionfiles)[0:len(FarJunctionfiles)/2]
FJ2_list= sorted(FarJunctionfiles)[len(FarJunctionfiles)/2:]
input_lines=fileinput.input(FJ1_list)
fout_FJ1.writelines(input_lines)
input_lines=fileinput.input(FJ2_list)
fout_FJ2.writelines(input_lines)

fout_FJ1.close()
fout_FJ2.close()

## open big indels files

f1_FarJunc= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/All_" + stem + "_1_indels.sam", mode ="rB")
f2_FarJunc= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/All_" + stem + "_2_indels.sam", mode="rB")

f1_genome = open(sorted(genomefiles)[0], mode="rB")
f2_genome = open(sorted(genomefiles)[1], mode="rB")
f1_reg = open(sorted(regfiles)[0], mode="rB")
f2_reg = open(sorted(regfiles)[1], mode="rB")
f1_junc= open(sorted(junctionfiles)[0], mode="rB")
f2_junc= open(sorted(junctionfiles)[1], mode="rB")
f1_unaligned=open(sorted(unalignedfiles)[0], mode="rB")
f2_unaligned=open(sorted(unalignedfiles)[1], mode="rB")


# ID file ReadID and different buckets.
    # [0] = readID
    # [1] = R2 in genome        
    # [2] = R2 in genome anomaly  
    # [3] = reg
    # [4] = reg anom
    # [5] = junc
    # [6] = junc anom
    # [7] = FarJunc
    # [8] = FarJunc anom
    # [9] = unaligned
    # [10] = unmapped


IDfile = open(args.FJDir+"FJIndelsreports//IDs_"+stem+".txt", mode= "w")
IDfile.write("ID\tclass\tR1_offset\tR1_MAPQ\tR1_AS\tR1_NumN\tR1_Readlength\tR1_JuncName\tR1_strand\tR2_offset\tR2_MAPQ\tR2_AS\tR2_NumN\tR2_Readlength\tR2_JuncName\tR2_strand\n")

#populate all reads and junctions into separate dictionaries
AllFJRead1= {}
AllFJRead2= {}
AllJunctions = {}
genomeDict = {}  # for all these dictionaries, [0] = reg, [1] = anom
regDict = {}        # [2] = sum of AS, [3] = read length
juncDict = {}
FJDict = {}
unalignedDict = {}
unmappedDict= {} # start with all readIDs.  if a partner is seen, then remove from list.


# populate AllFJRead1 dictionary - all read 1's from FarJunction alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening FarJunc _1 file"

for line_raw in f1_FarJunc:
    if line_raw[0] =="@":
        continue
    
    FJ1read = ReadInfoFJ(line_raw)
    
    if FJ1read.offset <= (150+FJ1read.indel-window) and (FJ1read.offset+FJ1read.NumOfBases)>= 150+FJ1read.indel+window:    
        AllFJRead1[FJ1read.ID] = [line_raw, 0]
        if FJ1read.junction not in AllJunctions:
            AllJunctions[FJ1read.junction]=0
        AllJunctions[FJ1read.junction] +=1
        unmappedDict[FJ1read.ID] = FJ1read.junction
        
f1_FarJunc.close()
IDfile.flush()

# populate AllFJRead2 dictionary - all read 2's from FarJunc alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening farJunc _2 file"
for line_raw in f2_FarJunc:
    if line_raw[0] =="@":
        continue
    
    FJ2read = ReadInfoFJ(line_raw)

#    if FJ1read.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
#        print "ERROR AT LINE 409"

    # if R1 and R2 both in Far Junc, then add to FJ-FJ list
    if FJ2read.ID in AllFJRead1:
        #print "found FJ read"
        #AllFJRead1[FJ2read.ID][1]="FJ"
        if FJ2read.offset <= (150+FJ2read.indel-window) and (FJ2read.offset+FJ2read.NumOfBases)>= 150+FJ2read.indel+window and AllFJRead1[FJ2read.ID][1]==0:    
            FJDict = AddToDict("FJ",FJDict,line_raw,AllFJRead1[FJ2read.ID][0])
            AllFJRead1[FJ2read.ID][1]="FJ"
            if FJ2read.ID in unmappedDict:
                del unmappedDict[FJ2read.ID]
        # otherwise add to F2 read
    else:
        AllFJRead2[FJ2read.ID]= [line_raw, 0]
        unmappedDict[FJ2read.ID] = FJ2read.junction
        
    if FJ2read.junction not in AllJunctions:
        AllJunctions[FJ2read.junction]=0
        
    AllJunctions[FJ2read.junction]+=1
f2_FarJunc.close()
IDfile.flush()



# compare FJ read 1 to genome read 2
for line_raw in f2_genome:
    if line_raw[0] =="@":
        continue
    g2read = ReadInfoGenome(line_raw)
    
    if g2read.ID in AllFJRead1 and AllFJRead1[g2read.ID][1]==0:
        #print "found genome R2"+ g2read.ID     
        if g2read.ID in unmappedDict:
            del unmappedDict[g2read.ID]
        genomeDict = AddToDict("genome", genomeDict, line_raw, AllFJRead1[g2read.ID][0])
        AllFJRead1[g2read.ID][1]="genome"

f2_genome.close()    
IDfile.flush()


# compare FJ read 2 to genome read 1

for line_raw in f1_genome:
    if line_raw[0] =="@":
        continue
    g1read = ReadInfoGenome(line_raw)
    
    if g1read.ID in AllFJRead2 and AllFJRead2[g1read.ID][1]==0:
        #print "found genome R1"+g1read.ID      
        if g1read.ID in unmappedDict:
            del unmappedDict[g1read.ID]
        genomeDict = AddToDict("genome", genomeDict, line_raw, AllFJRead2[g1read.ID][0])
        AllFJRead2[g1read.ID][1]="genome"
f1_genome.close()    
IDfile.flush()


    
# compare FJ read 1 to reg read 2
for line_raw in f2_reg:
    if line_raw[0] =="@":
        continue
    reg2read = ReadInfoJunc(line_raw)

    if reg2read.offset <= (150-window) and (reg2read.offset+reg2read.NumOfBases)>=( 150+window):    
        if reg2read.ID in AllFJRead1 and AllFJRead1[reg2read.ID][1]==0:
            #print "found reg R2" + reg2read.ID            
            if reg2read.ID in unmappedDict:
                del unmappedDict[reg2read.ID]
            regDict = AddToDict("reg", regDict, line_raw, AllFJRead1[reg2read.ID][0])
            AllFJRead1[reg2read.ID][1]="reg"
f2_reg.close()
IDfile.flush()

   

# compare FJ read 2 to reg read 1

for line_raw in f1_reg:
    if line_raw[0] =="@":
        continue
    reg1read = ReadInfoJunc(line_raw)
    
    if reg1read.offset <= (150-window) and (reg1read.offset+reg1read.NumOfBases)>=( 150+window):    
        if reg1read.ID in AllFJRead2 and AllFJRead2[reg1read.ID][1]==0:
            #print "found reg R1: " + reg1read.ID            
            if reg1read.ID in unmappedDict:
                del unmappedDict[reg1read.ID]
            regDict = AddToDict("reg", regDict, line_raw, AllFJRead2[reg1read.ID][0])
            AllFJRead2[reg1read.ID][1]="reg"
f1_reg.close()
IDfile.flush()

     
# compare FJ read 1 to junc read 2
for line_raw in f2_junc:
    if line_raw[0] =="@":
        continue
    junc2read = ReadInfoJunc(line_raw)
    
    if junc2read.offset <= (150-window) and (junc2read.offset+junc2read.NumOfBases)>=( 150+window):
        if junc2read.ID in AllFJRead1 and AllFJRead1[junc2read.ID][1]==0:
            #print "found junc R2 " + junc2read.ID            
            if junc2read.ID in unmappedDict:
                del unmappedDict[junc2read.ID]
            juncDict = AddToDict("junc", juncDict, line_raw, AllFJRead1[junc2read.ID][0])
            AllFJRead1[junc2read.ID][1]="junc"
f2_junc.close()
IDfile.flush()



# compare FJ read 2 to junc read 1
for line_raw in f1_junc:
    if line_raw[0] =="@":
        continue
    junc1read = ReadInfoJunc(line_raw)

    if junc1read.offset <= (150-window) and (junc1read.offset+junc1read.NumOfBases)>= (150+window):
        if junc1read.ID in AllFJRead2 and AllFJRead2[junc1read.ID][1]==0:
            #print "found junc R1: " + junc1read.ID            
            if junc1read.ID in unmappedDict:
                del unmappedDict[junc1read.ID]            
            juncDict = AddToDict("junc", juncDict, line_raw, AllFJRead2[junc1read.ID][0])
            AllFJRead2[junc1read.ID][1]="junc"
f1_junc.close()
IDfile.flush()


  
# compare FJ read 1 to unaligned read 2

for line_raw in f2_unaligned:
    if line_raw[0]=="@":
        if line_raw[1:][:-3] in AllFJRead1 and AllFJRead1[line_raw[1:][:-3]][1]==0:
            if line_raw[1:][:-3] in unmappedDict:
                del unmappedDict[line_raw[1:][:-3]]
            unalignedDict = AddToDict("unaligned", unalignedDict, line_raw, AllFJRead1[line_raw[1:][:-3]][0])
            AllFJRead1[line_raw[1:][:-3]][1]="unaligned"
f2_unaligned.close()
IDfile.flush()



# compare FJ read 2 to unaligned read 1
       
for line_raw in f1_unaligned:
    if line_raw[0]=="@":
        if line_raw[1:][:-3] in AllFJRead2 and AllFJRead2[line_raw[1:][:-3]][1]==0:
            if line_raw[1:][:-3] in unmappedDict:
                del unmappedDict[line_raw[1:][:-3]]                
            unalignedDict = AddToDict("unaligned", unalignedDict, line_raw, AllFJRead2[line_raw[1:][:-3]][0])
            AllFJRead2[line_raw[1:][:-3]][1]="unaligned"
f1_unaligned.close()
IDfile.flush()

# output header
outputfile = "FJIndelsreports/" + stem + "indels_naive_report.txt"  
fout = open(args.FJDir+outputfile, mode = "w")
print "fout: " + args.FJDir+outputfile

fout.write("@Junction\tgenome\tgenome-anomaly\tgenome-pval\treg\treg-anomaly\treg-pval\tjunc\tjunc-anom\tjunc-pval\tFarJunc\tFarJunc-anom\tFarJunc-pval\tunaligned\tNoPartner\tNetPValue\n")

for key in unmappedDict:
    IDfile.write(key+"\t"+unmappedDict[key]+"\tUnmapped\n")
IDfile.close()
#
#
### TESTING MODE - SEE WHAT ALLFJREAD1 and 2 FILE SHOW
#AllReadOutfile=open("AllJuncDict.txt", mode="w")
#
#AllReadOutfile.write("AllFJRead1:\n")
#
#for key in AllFJRead1:
#    AllReadOutfile.write(key+"\t"+str(AllFJRead1[key][1])+"\n")
#
#
#AllReadOutfile.write("AllFJRead2:\n")
#
#for key in AllFJRead2:
#    AllReadOutfile.write(key+"\t"+str(AllFJRead2[key][1])+"\n")
#
#AllReadOutfile.close()
#


## WRITE ALL JUNCTIONS
for key in AllJunctions:

    for dict in [genomeDict, regDict, juncDict, FJDict, unalignedDict]:
        if key not in dict:
            dict[key] = [0,0,0.0,0.0]

    
    NumUnmapped = Counter(unmappedDict.values())[key]
# calculates P value for all genome/reg/junc/FJ dicts combined, excluding all anomalous reads
    NetAS = genomeDict[key][2]+regDict[key][2]+juncDict[key][2]+FJDict[key][2]
    NetNumBases = genomeDict[key][3]+regDict[key][3]+juncDict[key][3]+FJDict[key][3]
    NetP = Pvalue(NetAS,NetNumBases)

# writing to output file
    fout.write(key+"\t") # write junction [0]
    fout.write(str(genomeDict[key][0])+"\t") # [1]  number of make-sense maps in genome
    fout.write(str(genomeDict[key][1])+"\t") # [2]  nonsense maps in genome
    fout.write(str(Pvalue(genomeDict[key][2], genomeDict[key][3]))+"\t") #[3] p value for genome reads
    fout.write(str(regDict[key][0])+"\t") # [4]  number of make-sense maps in reg
    fout.write(str(regDict[key][1])+"\t") # [5]  nonsense maps in reg
    fout.write(str(Pvalue(regDict[key][2], regDict[key][3]))+"\t") #[6] p value for reg reads
    fout.write(str(juncDict[key][0])+"\t") # [7]  number of make-sense maps in junc
    fout.write(str(juncDict[key][1])+"\t") # [8]  nonsense maps in junc
    fout.write(str(Pvalue(juncDict[key][2], juncDict[key][3]))+"\t") #[9] p value for junc reads
    fout.write(str(FJDict[key][0])+"\t") # [10]  number of make-sense maps in FJ
    fout.write(str(FJDict[key][1])+"\t") # [11]  nonsense maps in FJ
    fout.write(str(Pvalue(FJDict[key][2], FJDict[key][3]))+"\t") #[12] p value for FJ reads
    fout.write(str(unalignedDict[key][0])+"\t") # [13]  number junctions whose partner was unaligned
    fout.write(str(NumUnmapped)+"\t") # [14] - no read partner in any of the files.
    fout.write(str(NetP)+"\n") #[15] net P value of all non-anomalous reads
fout.close()



# takes alignments from Far Junctions and finds read partner. 
# tells if read partner makes sense or not
# final Output categories -- [0] R1 junc name
        #   [1] genome - R2 location < 100mill bp away
        #   [2] genome anomaly - R2 location > 100mill bp away
        #   [3] genome p value
        #   [4] reg - R2 closest location <100mill bp away
        #   [5] reg anomaly - R2 closest location > 100mill bp away
        #   [6] reg p value
        #   [7] junc / scrambled - R2 closest location <100mill bp away
        #   [8] junc anomaly - R2 closest location > 100mill bp away
        #   [9] junc p value
        #   [10] FarJunc - R2 aligned to same FarJunc
        #   [11] FarJunc anomaly - R2 aligned to diff FarJunc
        #   [12] FarJunc p value
        #   [13] unaligned - R2 didn't align
        #   [14] unmapped - R2 missing 
        #  [15] P value for all non-anomaly classes

