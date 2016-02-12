# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 16:12:27 2015

@author: Gillian
"""
# takes alignments from Far Junctions and finds read partner. 
# all Far Junction and Scrambled junction reads where the alignment does not overlap the 
# Junction by args.window bp will be thrown out.

# This program then tells if read partners "makes sense" or not
# final Output categories -- [0] R1 junc name
            #   [1] genome - R2 location < 100K bp away
            #   [2] genome anomaly - R2 location > 100K bp away
            #   [3] genome p value
            #   [4] reg - R2 closest location <100Kbp away
            #   [5] reg anomaly - R2 closest location > 100Kbp away
            #   [6] reg p value
            #   [7] junc / scrambled - R2 closest location <100Kbp away
            #   [8] junc anomaly - R2 closest location > 100Kbp away
            #   [9] junc p value
            #   [10] FarJunc - R2 aligned to same FarJunc
            #   [11] FarJunc anomaly - R2 aligned to diff FarJunc
            #   [12] FarJunc p value
            #   [13] unaligned - R2 didn't align
            #   [14] unmapped - R2 missing in action
            #   [15] P val for all non-anomaly classes



import argparse
import os
import glob
from math import sqrt
from scipy.stats import norm
from collections import Counter

def AddToDict(inputtype, TargetDict, line_raw_comparison, line_raw_FJ):

    lineFJ = ReadInfoFJ(line_raw_FJ)
    if lineFJ.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
        print "ERROR AT LINE 42"
        print inputtype
        print line_raw_comparison
        print line_raw_FJ
    
    
    if lineFJ.junction not in TargetDict:  # add junction to target dictionary if it doesn't exist
        TargetDict[lineFJ.junction] = [0,0,0.0,0.0] 
    
        
    if inputtype=="FJ": # if comparing Far Junc to Far Junc, they have to be identical
        line2= ReadInfoFJ(line_raw_comparison)    
        IDfileoutputR1 = lineFJ.ID+"\t"+lineFJ.junction+"\t"+lineFJ.refstrand+"\t"
        IDfileoutputR2 = line2.refstrand + "\t" + line2.junction + "\t" + str(line2.AS) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.offset) + "\t" + str(line2.NumOfBases)
        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

        if lineFJ.junction == line2.junction and lineFJ.refstrand in ["0","16"] and line2.refstrand in ["0","16"] and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
            IDfiletype = "FarJunction\t"
        else:
            TargetDict[lineFJ.junction][1]+=1
            IDfiletype = "FarJuncAnom\t"
            addAS = 0.0
            addNumofBases = 0.0
            
        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases

        IDfile.write(IDfileoutputR1+IDfiletype+IDfileoutputR2+"\n")


    if inputtype=="reg" or inputtype =="junc": #if reg or junc read, then one side has to be within 100KB, and meets refstrand criteria below
        line2 = ReadInfoJunc(line_raw_comparison)

        IDfileoutputR1 = lineFJ.ID+"\t"+lineFJ.junction+"\t"+lineFJ.refstrand+"\t"
        IDfileoutputR2 = line2.refstrand + "\t" + line2.junction + "\t" + str(line2.AS) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.offset) + "\t" + str(line2.NumOfBases)

        if inputtype =="reg": IDfiletype = "Regular\t"
        if inputtype == "junc": IDfiletype = "Junction\t"

        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

        
        if lineFJ.strand_left==line2.strand and lineFJ.strand_left=="-" and lineFJ.chr_left == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_left), int(lineFJ.loc_left)+100001) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand: 
            TargetDict[lineFJ.junction][0]+=1 
        elif lineFJ.strand_left==line2.strand and lineFJ.strand_left=="-" and lineFJ.chr_left == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_left), int(lineFJ.loc_left)+100001) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand: 
            TargetDict[lineFJ.junction][0]+=1  
        elif lineFJ.strand_left==line2.strand and lineFJ.strand_left=="+" and lineFJ.chr_left == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_left)-100000,int(lineFJ.loc_left)) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_left==line2.strand and lineFJ.strand_left=="+" and lineFJ.chr_left == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_left)-100000,int(lineFJ.loc_left)) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
                     
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "-" and lineFJ.chr_right == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_right)-100000,int(lineFJ.loc_right)) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "-" and lineFJ.chr_right == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_right)-100000,int(lineFJ.loc_right)) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "+" and lineFJ.chr_right == line2.chr and int(line2.loc_1) in range(int(lineFJ.loc_right),int(lineFJ.loc_right)+100001) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right==line2.strand and lineFJ.strand_right == "+" and lineFJ.chr_right == line2.chr and int(line2.loc_2) in range(int(lineFJ.loc_right),int(lineFJ.loc_right)+100001) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
            
        else:
            TargetDict[lineFJ.junction][1]+=1
            if inputtype =="reg": IDfiletype = "RegAnomaly\t"
            if inputtype == "junc": IDfiletype = "JuncAnomaly\t"
            addAS = 0.0
            addNumofBases = 0.0


        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
       
        IDfile.write(IDfileoutputR1+IDfiletype+IDfileoutputR2+"\n")


    if inputtype == "genome": #comparing FJ to genome, has to be within 100Kbp, meet ref strand criteria (opp refstrand if + read, same refstrand if - read)

        line2 = ReadInfoGenome(line_raw_comparison)

        IDfileoutputR1 = lineFJ.ID+"\t"+lineFJ.junction+"\t"+lineFJ.refstrand+"\t"
        IDfileoutputR2 = line2.refstrand + "\t" + line2.chr+":"+line2.loc + "\t" + str(line2.AS) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.offset) + "\t" + str(line2.NumOfBases)
        IDfiletype = "genome\t"
        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

 
        if lineFJ.strand_left=="-" and lineFJ.chr_left == line2.chr and int(line2.loc) in range(int(lineFJ.loc_left), int(lineFJ.loc_left)+100001) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand: 
            TargetDict[lineFJ.junction][0]+=1  
        elif lineFJ.strand_left=="+" and lineFJ.chr_left == line2.chr and int(line2.loc) in range(int(lineFJ.loc_left)-100000,int(lineFJ.loc_left)) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right == "-" and lineFJ.chr_right == line2.chr and int(line2.loc) in range(int(lineFJ.loc_right)-100000,int(lineFJ.loc_right)) and lineFJ.refstrand in ("0","16") and lineFJ.refstrand ==line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        elif lineFJ.strand_right == "+" and lineFJ.chr_right == line2.chr and int(line2.loc) in range(int(lineFJ.loc_right),int(lineFJ.loc_right)+100001) and lineFJ.refstrand in ("0","16") and line2.refstrand in ("0","16") and lineFJ.refstrand!=line2.refstrand:
            TargetDict[lineFJ.junction][0] +=1
        else:
            TargetDict[lineFJ.junction][1] +=1
            IDfiletype = "genomAnomaly\t"
            addAS = 0.0
            addNumofBases = 0.0           
            
        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
        IDfile.write(IDfileoutputR1+IDfiletype+IDfileoutputR2+"\n")

    if inputtype == "unaligned":
        if lineFJ.ID not in FJDict:
            if lineFJ.junction not in TargetDict:
                TargetDict[lineFJ.junction] = 0 
            TargetDict[lineFJ.junction]+=1
        IDfileoutputR1 = lineFJ.ID+"\t"+lineFJ.junction+"\t"+lineFJ.refstrand+"\t"
        IDfiletype = "unaligned\t"
        IDfile.write(IDfileoutputR1 + IDfiletype + "\n")
        
     

    return TargetDict

def ID(string):
    if string[-2:] == "/1" or string[-2:] == "/2":
        return string[:-3]
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
        self.MAPQ = int(line[4])
        self.AS= int(line[11].split(":")[2])
        self.NumOfBases = len(line[9])
        self.offset = int(line[3])        
        
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

for name in glob.glob(args.FJDir + "FarJunctionAlignments/"+ stem + "/*.sam"):
    print name
    FarJunctionfiles.append(name)

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

f1_FarJunc= open(sorted(FarJunctionfiles)[0], mode ="rB")
f2_FarJunc= open(sorted(FarJunctionfiles)[1], mode="rB")

f1_genome = open(sorted(genomefiles)[0], mode="rB")
f2_genome = open(sorted(genomefiles)[1], mode="rB")
f1_reg = open(sorted(regfiles)[0], mode="rB")
f2_reg = open(sorted(regfiles)[1], mode="rB")
f1_junc= open(sorted(junctionfiles)[0], mode="rB")
f2_junc= open(sorted(junctionfiles)[1], mode="rB")
f1_unaligned=open(sorted(unalignedfiles)[0], mode="rB")
f2_unaligned=open(sorted(unalignedfiles)[1], mode="rB")

#
#
#for file1 in FarJunctionfiles:
##        print file1
#    sorted(FarJunctionfiles[0])    
#    if "_1." in file1 or "_1_" in file1 or "001_" in file1:
#        f1_FarJunc = open(file1, mode="rB")
#        print "FarJunc_1: " + file1
#    else:
#        f2_FarJunc = open(file1, mode="rB")
#        print "FarJunc_2: " + file1
#
#for file1 in genomefiles:
#    if "_1." in file1 or "_1_" in file1 or "001_" in file1:
#        f1_genome = open(file1, mode="rB")
#        print "genome_1: " + file1
#    else:
#        f2_genome = open(file1, mode="rB")
#        print "genome_2: " + file1
#
#for file1 in regfiles:
#    if "_1." in file1 or "_1_" in file1 or "001_" in file1:
#        f1_reg = open(file1, mode="rB")
#        print "reg_1: " + file1
#    else:
#        f2_reg = open(file1, mode="rB")
#        print "reg_2: " + file1
#
#for file1 in junctionfiles:
#    if "_1." in file1 or "_1_" in file1 or "001_" in file1:
#        f1_junc = open(file1, mode="rB")
#        print "scramble_1: " + file1
#
#    else:
#        f2_junc = open(file1, mode="rB")
#        print "scramble_2: " + file1
#
#for file1 in unalignedfiles:
#    if "_1." in file1 or "_1_" in file1 or "001_" in file1:
#        f1_unaligned = open(file1, mode="rB")
#        print "unalign_1: " + file1
#    else:
#        f2_unaligned = open(file1, mode="rB")
#        print "unalign_2: " + file1

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


IDfile = open(args.FJDir+"reports/IDs_"+stem+".txt", mode= "w")
IDfile.write("ID\tFarJunction\tFarJuncStrand\tPartnerBin\tPartnerStrand\tPartnerlocation\tPartnerAS\tPartnerMAPQ\tPartnerOffset\tPartnerReadLength\n\n")
 

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
    if FJ1read.offset < (150-window) and (FJ1read.offset+FJ1read.NumOfBases)> 150+window:    
        AllFJRead1[FJ1read.ID] = line_raw
        if FJ1read.junction not in AllJunctions:
            AllJunctions[FJ1read.junction]=0
        AllJunctions[FJ1read.junction] +=1
        unmappedDict[FJ1read.ID] = FJ1read.junction
    if FJ1read.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
        print "ERROR AT LINE 395"
f1_FarJunc.close()
IDfile.flush()

# populate AllFJRead2 dictionary - all read 2's from FarJunc alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening farJunc _2 file"
for line_raw in f2_FarJunc:
    if line_raw[0] =="@":
        continue
    
    FJ2read = ReadInfoFJ(line_raw)


    if FJ1read.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
        print "ERROR AT LINE 409"

    # if R1 and R2 both in Far Junc, then add to FJ-FJ list
    if FJ2read.ID in AllFJRead1:
        if FJ2read.offset < (150-window) and (FJ2read.offset+FJ2read.NumOfBases)> 150+window:    
            FJDict = AddToDict("FJ",FJDict,line_raw,AllFJRead1[FJ2read.ID])
            if FJ2read.ID in unmappedDict:
                del unmappedDict[FJ2read.ID]
        # otherwise add to F2 read
    else:
        AllFJRead2[FJ2read.ID]= line_raw
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
    
    if g2read.ID in AllFJRead1:
        if g2read.ID in unmappedDict:
            del unmappedDict[g2read.ID]
        genomeDict = AddToDict("genome", genomeDict, line_raw, AllFJRead1[g2read.ID])
f2_genome.close()    
IDfile.flush()

    
# compare FJ read 1 to reg read 2
for line_raw in f2_reg:
    if line_raw[0] =="@":
        continue
    reg2read = ReadInfoJunc(line_raw)
    
    if reg2read.ID in AllFJRead1:
        if reg2read.ID in unmappedDict:
            del unmappedDict[reg2read.ID]
        regDict = AddToDict("reg", regDict, line_raw, AllFJRead1[reg2read.ID])
f2_reg.close()
IDfile.flush()

        
# compare FJ read 1 to junc read 2
for line_raw in f2_junc:
    if line_raw[0] =="@":
        continue
    junc2read = ReadInfoJunc(line_raw)
    
    if junc2read.offset < (150-window) and (junc2read.offset+junc2read.NumOfBases)>( 150+window):
        if junc2read.ID in AllFJRead1:
            if junc2read.ID in unmappedDict:
                del unmappedDict[junc2read.ID]
            juncDict = AddToDict("junc", juncDict, line_raw, AllFJRead1[junc2read.ID])
f2_junc.close()
IDfile.flush()

  
# compare FJ read 1 to unaligned read 2

for line_raw in f2_unaligned:
    if line_raw[0]=="@":
        if line_raw[1:][:-3] in AllFJRead1:
            if line_raw[1:][:-3] in unmappedDict:
                del unmappedDict[line_raw[1:][:-3]]
            unalignedDict = AddToDict("unaligned", unalignedDict, line_raw, AllFJRead1[line_raw[1:][:-3]])
f2_unaligned.close()
IDfile.flush()


# compare FJ read 2 to genome read 1

for line_raw in f1_genome:
    if line_raw[0] =="@":
        continue
    g1read = ReadInfoGenome(line_raw)
    
    if g1read.ID in AllFJRead2:
        if g1read.ID in unmappedDict:
            del unmappedDict[g1read.ID]
        genomeDict = AddToDict("genome", genomeDict, line_raw, AllFJRead2[g1read.ID])
f1_genome.close()    
IDfile.flush()

# compare FJ read 2 to reg read 1

for line_raw in f1_reg:
    if line_raw[0] =="@":
        continue
    reg1read = ReadInfoJunc(line_raw)
    
    if reg1read.ID in AllFJRead2:
        if reg1read.ID in unmappedDict:
            del unmappedDict[reg1read.ID]
        regDict = AddToDict("reg", regDict, line_raw, AllFJRead2[reg1read.ID])
f1_reg.close()
IDfile.flush()

# compare FJ read 2 to junc read 1
for line_raw in f1_junc:
    if line_raw[0] =="@":
        continue
    junc1read = ReadInfoJunc(line_raw)

    if junc1read.offset < (150-window) and (junc1read.offset+junc1read.NumOfBases)> (150+window):
        if junc1read.ID in AllFJRead2:
            if junc1read.ID in unmappedDict:
                del unmappedDict[junc1read.ID]            
            juncDict = AddToDict("junc", juncDict, line_raw, AllFJRead2[junc1read.ID])
f1_junc.close()
IDfile.flush()

# compare FJ read 2 to unaligned read 1
       
for line_raw in f1_unaligned:
    if line_raw[0]=="@":
        if line_raw[1:][:-3] in AllFJRead2:
            if line_raw[1:][:-3] in unmappedDict:
                del unmappedDict[line_raw[1:][:-3]]                
            unalignedDict = AddToDict("unaligned", unalignedDict, line_raw, AllFJRead2[line_raw[1:][:-3]])
f1_unaligned.close()
IDfile.flush()

# output header
outputfile = "reports/" + stem + "_naive_report.txt"  
fout = open(args.FJDir+outputfile, mode = "w")
print "fout: " + args.FJDir+outputfile

fout.write("@Junction\tgenome\tgenome-anomaly\tgenome-pval\treg\treg-anomaly\treg-pval\tjunc\tjunc-anom\tjunc-pval\tFarJunc\tFarJunc-anom\tFarJunc-pval\tunaligned\tNoPartner\tNetPValue\n")

for key in unmappedDict:
    IDfile.write(key+"\tUnmapped\n")
IDfile.close()

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
        #   [1] genome - R2 location < 100K bp away
        #   [2] genome anomaly - R2 location > 100K bp away
        #   [3] genome p value
        #   [4] reg - R2 closest location <100Kbp away
        #   [5] reg anomaly - R2 closest location > 100Kbp away
        #   [6] reg p value
        #   [7] junc / scrambled - R2 closest location <100Kbp away
        #   [8] junc anomaly - R2 closest location > 100Kbp away
        #   [9] junc p value
        #   [10] FarJunc - R2 aligned to same FarJunc
        #   [11] FarJunc anomaly - R2 aligned to diff FarJunc
        #   [12] FarJunc p value
        #   [13] unaligned - R2 didn't align
        #   [14] unmapped - R2 missing 
        #  [15] P value for all non-anomaly classes

