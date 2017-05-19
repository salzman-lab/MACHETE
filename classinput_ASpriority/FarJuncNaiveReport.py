# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 16:12:27 2015

@author: Gillian
"""
# takes alignments from Far Junctions and finds read partner. 
# all Far Junction and Scrambled junction reads where the alignment does not overlap the 
# Junction by args.window bp will be thrown out.

## Also creates class input files for Far Junctions (R1 = FJ, R2 = something else)

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
from math import ceil
from scipy.stats import poisson
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
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand

        addAS = lineFJ.AS+line2.AS
        addNumofBases = lineFJ.NumOfBases + line2.NumOfBases

        if lineFJ.junction == line2.junction and lineFJ.refstrand!=line2.refstrand:
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
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand


        if inputtype == "junc": 
            IDfiletype = "FJbad,Junction"
            addAS = 0.0
            addNumofBases = 0.0
       
        if inputtype =="reg": 
            IDfiletype = "FJbad,RegAnomaly"
            addAS = 0.0
            addNumofBases = 0.0
            
            if lineFJ.chr_left == line2.chr:
                if lineFJ.strand_left==line2.strand:
                    if lineFJ.strand_left=="-":
                        if lineFJ.refstrand ==line2.refstrand:
                            if int(lineFJ.loc_left) <= int(line2.loc_2):
                                IDfiletype = "FJgood,Regular"
                    elif lineFJ.strand_left=="+":
                        if lineFJ.refstrand !=line2.refstrand:
                            if int(lineFJ.loc_left) >= int (line2.loc_2):
                                IDfiletype = "FJgood,Regular"
            if IDfiletype == "FJbad,RegAnomaly":
                if lineFJ.chr_right == line2.chr:
                    if lineFJ.strand_right==line2.strand:
                        if lineFJ.strand_right=="-":
                            if lineFJ.refstrand ==line2.refstrand:
                                if int(lineFJ.loc_right) >= int(line2.loc_1):
                                    IDfiletype = "FJgood,Regular"
                        elif lineFJ.strand_right=="+":
                            if lineFJ.refstrand !=line2.refstrand:
                                 if int(lineFJ.loc_right) <= int(line2.loc_1):
                                    IDfiletype = "FJgood,Regular"
                               
        if IDfiletype == "FJgood,Regular":
            TargetDict[lineFJ.junction][0] +=1
            addAS = lineFJ.AS+line2.AS
            addNumofBases = lineFJ.NumOfBases + line2.NumOfBases
        else:
            TargetDict[lineFJ.junction][1] +=1

        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
       
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype == "genome": #comparing FJ to genome, meet ref strand criteria (opp refstrand if + read, same refstrand if - read)

        line2 = ReadInfoGenome(line_raw_comparison)

        ## output R1 -  offset, MAPQ, AS, #N, readlen, junc name, strand
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.loc) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.chr + "\t" + line2.refstrand


        IDfiletype = "FJbad,genomAnomaly"
        addAS = 0.0
        addNumofBases = 0.0         

        if lineFJ.chr_left == line2.chr:
            #print "checking left"
            if lineFJ.strand_left=="-":
                #print "left exon neg"
                if lineFJ.refstrand ==line2.refstrand:
                    #print "left strand reference good"                    
                    if int(lineFJ.loc_left) <= int(line2.loc):
                        #print "left exon in expected location"                      
                        IDfiletype = "FJgood,genome"           
            elif lineFJ.strand_left=="+":
                #print "left exon pos"
                if lineFJ.refstrand != line2.refstrand:
                    #print "left strand good"
                    if int(lineFJ.loc_left) >= int(line2.loc):      
                        #print "left strand location good"
                        IDfiletype = "FJgood,genome"
        ## if left not the same, then compare right        
        if IDfiletype == "FJbad,genomAnomaly": 
            if lineFJ.chr_right == line2.chr:  
                #print "right chromosome correct"
                if lineFJ.strand_right=="-":
                    #print "right exon neg"
                    if lineFJ.refstrand ==line2.refstrand:
                        #print "right reference strand correct"
                        if int(lineFJ.loc_right) >= int(line2.loc):
                            #print "location correct"
                            IDfiletype = "FJgood,genome"           
                elif lineFJ.strand_right=="+":
                    #print " right strand pos"
                    if lineFJ.refstrand !=line2.refstrand:
                        #print "right ref strand correct"
                        if int(lineFJ.loc_right) <= int(line2.loc):
                            #print "right exon location correct"
                            IDfiletype = "FJgood,genome"
                   
  
        if IDfiletype == "FJgood,genome":
            TargetDict[lineFJ.junction][0] +=1
            addAS = lineFJ.AS+line2.AS
            addNumofBases = lineFJ.NumOfBases + line2.NumOfBases
        else:
            TargetDict[lineFJ.junction][1] +=1

        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
          
        TargetDict[lineFJ.junction][2] += addAS
        TargetDict[lineFJ.junction][3] += addNumofBases 
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")

    if inputtype == "unaligned":
        if lineFJ.ID not in FJDict:
            try: TargetDict[lineFJ.junction][0]+=1
            except: 
                print lineFJ.junction
                print TargetDict[lineFJ.junction]
                print line_raw_comparison
            
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction+"\t"+lineFJ.refstrand        
        IDfiletype = "unaligned"
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\n")
        
     

    return TargetDict

def ID(string):
    if string[-2:] == "/1" or string[-2:] == "/2":
        return string[:-2]
    else:
        return string

## actual MM - round up to nearest integer = X
## expected MM  - no need to round = lambda
## return 1- poisson.cdf(X, lambda)

def Pvalue(AS, NumBases):
    ExpectedMMrate = 0.01
    ExpectedMM = ExpectedMMrate*float(NumBases)
    if NumBases == 0.0:
        return "-"
    ActualMM = int(ceil(float(AS)/(-6.0)))
    prob = 1 - poisson.cdf(ActualMM, ExpectedMM)
    return prob


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
        if "XS:i:" in line[12]:
            self.NumN=line[13][5:]
            self.adjAS=int(line[11].split(":")[2])+int(line[13][5:])
        else:
            self.NumN=line[12][5:]
            self.adjAS=int(line[11].split(":")[2])+int(line[12][5:])
        



        JuncInfo = line[2].replace(":"," ").replace("|"," ").split(" ")
        self.chr_left=JuncInfo[0]
        self.loc_left= JuncInfo[2]
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
            self.adjAS=int(line[11].split(":")[2])+int(line[13][5:])
        else:
            self.NumN=line[12][5:]
            self.adjAS=int(line[11].split(":")[2])+int(line[12][5:])

        
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
            self.adjAS=int(line[11].split(":")[2])+int(line[13][5:])
        else:
            self.NumN=line[12][5:]
            self.adjAS=int(line[11].split(":")[2])+int(line[12][5:])

        
#=========================================
#start here

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--stem", required=True, help = "stem name of file to generate report")
parser.add_argument("-f", "--FJDir", required = True, help = "path to aligned junction reads")
parser.add_argument("-i", "--origDir", required=True, help = "path to orig dir containing genome reads")
parser.add_argument("-w", "--window", required=True, help = "# of bases needed on each side of the junction")
args = parser.parse_args()
window= int(args.window)


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
    FarJunctionfiles.append(name)

for name in glob.glob(os.path.join(args.origDir,"genome/*" + stem + "*.sam")):        
    if "sorted" not in name:       
        genomefiles.append(name)

for name in glob.glob(os.path.join(args.origDir,"reg/*" + stem + "*.sam")):
    if "sorted" not in name:       
        regfiles.append(name) 
for name in glob.glob(os.path.join(args.origDir,"junction/*" + stem + "*.sam")):
    if "sorted" not in name:       
        junctionfiles.append(name) 
for name in glob.glob(os.path.join(args.origDir,"unaligned/*" + stem + "*.fq")):
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


IDfile = open(args.FJDir+"reports/temp_IDs_"+stem+".txt", mode= "w")
IDfile.write("ID\tclass\tR1_offset\tR1_MAPQ\tR1_adjAS\tR1_NumN\tR1_Readlength\tR1_JuncName\tR1_strand\tR2_offset\tR2_MAPQ\tR2_adjAS\tR2_NumN\tR2_Readlength\tR2_JuncName\tR2_strand\n")

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

linecounter=0
newjunccounter=0
goodlinecounter=0

for line_raw in f1_FarJunc:
    if line_raw[0] =="@":
        continue
    linecounter+=1
    
    FJ1read = ReadInfoFJ(line_raw)
    if FJ1read.offset <= (150-window) and (FJ1read.offset+FJ1read.NumOfBases)>= 150+window:  
        goodlinecounter+=1
        # Dict AllFJRead1 contains key = read ID of all FJ R1
        # Value = [ FJ read info, indicator of which library the R2 is in]
        # indicator = 0 if no R2 detected
        # indicator = FJ if R2 in FJ, genome if R2 in genome, reg if r2 in reg, etc...
        AllFJRead1[FJ1read.ID] = [line_raw, 0]
        if FJ1read.junction not in AllJunctions:
            AllJunctions[FJ1read.junction]=0
            newjunccounter+=1
        AllJunctions[FJ1read.junction] +=1
        unmappedDict[FJ1read.ID] = FJ1read.junction
        
f1_FarJunc.close()
print "lines "+ str(linecounter)
print "good lines " + str(goodlinecounter)
print "independent juncs" + str(newjunccounter)

IDfile.flush()

# populate AllFJRead2 dictionary - all read 2's from FarJunc alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening farJunc _2 file"


linecounter=0
newjunccounter=0
goodlinecounter=0
overlapwithFJ1=0

for line_raw in f2_FarJunc:
    if line_raw[0] =="@":
        continue
    linecounter+=1
    FJ2read = ReadInfoFJ(line_raw)

#    if FJ1read.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
#        print "ERROR AT LINE 409"

    # if R1 and R2 both in Far Junc, then add to FJ-FJ list
    if FJ2read.ID in AllFJRead1:
        overlapwithFJ1+=1
        #print "found FJ read"
        #AllFJRead1[FJ2read.ID][1]="FJ"
        if FJ2read.offset <= (150-window) and (FJ2read.offset+FJ2read.NumOfBases)>= 150+window and AllFJRead1[FJ2read.ID][1]==0:    
            FJDict = AddToDict("FJ",FJDict,line_raw,AllFJRead1[FJ2read.ID][0])

            AllFJRead1[FJ2read.ID][1]="FJ"
            if FJ2read.ID in unmappedDict:
                del unmappedDict[FJ2read.ID]
        # otherwise add to F2 read
    else:
        if FJ2read.offset <= (150-window) and (FJ2read.offset+FJ2read.NumOfBases)>= 150+window:    
            goodlinecounter+=1
            AllFJRead2[FJ2read.ID]= [line_raw, 0]
            unmappedDict[FJ2read.ID] = FJ2read.junction        
            if FJ2read.junction not in AllJunctions:
                newjunccounter+=1
                AllJunctions[FJ2read.junction]=0
        
            AllJunctions[FJ2read.junction]+=1
f2_FarJunc.close()

print "lines "+ str(linecounter)
print "good lines " + str(goodlinecounter)
print "independent juncs" + str(newjunccounter)
print "overlapping with FJ1 " + str(overlapwithFJ1)
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

    if reg2read.offset <= (150-window) and (reg2read.offset+reg2read.NumOfBases)>=(150+window):    
        if reg2read.ID in AllFJRead1:
            if AllFJRead1[reg2read.ID][1]==0 or AllFJRead1[reg2read.ID][1]=="genome":
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
        if reg1read.ID in AllFJRead2:
            if AllFJRead2[reg1read.ID][1]==0 or AllFJRead2[reg1read.ID][1]=="genome":
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
        readID = line_raw.strip().split(" ")[0][1:]
        readID = ID(readID)
        if readID in AllFJRead1 and AllFJRead1[readID][1]==0:
            if readID in unmappedDict:
                del unmappedDict[readID]
            unalignedDict = AddToDict("unaligned", unalignedDict, line_raw, AllFJRead1[readID][0])
            AllFJRead1[readID][1]="unaligned"
f2_unaligned.close()
IDfile.flush()



# compare FJ read 2 to unaligned read 1
       
for line_raw in f1_unaligned:
    if line_raw[0]=="@":
        readID = line_raw.strip().split(" ")[0][1:]
        readID = ID(readID)
        if readID in AllFJRead2 and AllFJRead2[readID][1]==0:
            if readID in unmappedDict:
                del unmappedDict[readID]                
            unalignedDict = AddToDict("unaligned", unalignedDict, line_raw, AllFJRead2[readID][0])
            AllFJRead2[readID][1]="unaligned"
f1_unaligned.close()
IDfile.flush()

# output header
outputfile = "reports/" + stem + "_naive_report.txt"  
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




#############################################################################
## This section of code takes the written ID file above (temp_IDs_STEM.txt) and
## removes duplicate entries of genome and reg. The same readID may be found 
## in both libraries and would both be in the ID file.
## The new ID file removes duplicates and only keeps the readID with the 
## best alignment score.
tempIDfile = open(args.FJDir+"reports/temp_IDs_"+stem+".txt", mode= "rU")
newIDfile = open(args.FJDir+"reports/IDs_"+stem+".txt", mode= "w")

##grep col 2 for "genom", "Regular" or "RegAnomaly". if not found ,write to
## new file immediately.
## if found, feed into dictionary (key=readID, value= entire line from temp file)
## if duplicate entry, then compare R2 AS. if AS larger, then replace
## value with new value from new R2
## at completion of file, write entire dictionary into new ID file.

GenomeAndRegReadIDs={}

for line in tempIDfile:
    line=line.strip()
    if "unaligned" in line:
        continue

    if "Unmapped" in line:
        continue
    
    readID = line.split("\t")[0]
    classID = line.split("\t")[1]
    AS_new=line.split("\t")[11]
    

    
    if "genom" in classID:
        ## if readID has been seen previously, then replace value in dictionary
    ## only if AS is greater.
        if readID in GenomeAndRegReadIDs:
            AS_old=GenomeAndRegReadIDs[readID].split("\t")[11]
            if int(AS_new)>int(AS_old):
                GenomeAndRegReadIDs[readID]=line
        else:
            GenomeAndRegReadIDs[readID]=line
    elif "Regular" in classID:
        ## do the same if reg 
        if readID in GenomeAndRegReadIDs:
            AS_old=GenomeAndRegReadIDs[readID].split("\t")[11]
            if int(AS_new)>int(AS_old):
                GenomeAndRegReadIDs[readID]=line
        else:
            GenomeAndRegReadIDs[readID]=line
    elif "RegAnomaly" in classID:
        ## do the same if reg anomaly
        if readID in GenomeAndRegReadIDs:
            AS_old=GenomeAndRegReadIDs[readID].split("\t")[11]
            if int(AS_new)>int(AS_old):
                GenomeAndRegReadIDs[readID]=line
        else:
            GenomeAndRegReadIDs[readID]=line
    else:
        ## if not genome/genome anomaly/ reg/ reg anomaly then 
    ## write line directly in new file.
        newIDfile.write(line)

for entry in GenomeAndRegReadIDs:
    newIDfile.write(GenomeAndRegReadIDs[entry]+"\n")

tempIDfile.close()
newIDfile.close()

