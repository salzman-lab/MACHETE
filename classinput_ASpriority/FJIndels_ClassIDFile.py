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


def AddToDict(inputtype, line_raw_comparison, line_raw_FJ):

    lineFJ = ReadInfoFJ(line_raw_FJ)
        
    if inputtype=="FJ": # if comparing Far Junc to Far Junc, they have to be identical
        line2= ReadInfoFJ(line_raw_comparison)
        
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction[:-5]+"\t"+lineFJ.refstrand        
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand

        if lineFJ.junction == line2.junction and lineFJ.refstrand in ["0","16"] and line2.refstrand in ["0","16"] and lineFJ.refstrand!=line2.refstrand:
#            TargetDict[lineFJ.junction][0] +=1
            IDfiletype = "FJgood,FarJunction,"+lineFJ.junction[-4:]
        else:
#            TargetDict[lineFJ.junction][1]+=1
            IDfiletype = "FJbad,FarJuncAnom,"+lineFJ.junction[-4:]
#            addAS = 0.0
#            addNumofBases = 0.0
#            
#        TargetDict[lineFJ.junction][2] += addAS
#        TargetDict[lineFJ.junction][3] += addNumofBases

        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype=="reg" or inputtype =="junc": #if reg or junc read and meets refstrand criteria below
        line2 = ReadInfoJunc(line_raw_comparison)

        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction[:-5]+"\t"+lineFJ.refstrand        
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand

        if inputtype == "junc": 
            IDfiletype = "FJbad,Junction,"+lineFJ.junction[-4:]
        
        if inputtype =="reg": 
            IDfiletype = "FJbad,RegAnomaly,"+lineFJ.junction[-4:]

            if lineFJ.chr_left == line2.chr:
                if lineFJ.strand_left==line2.strand:
                    if lineFJ.strand_left=="-":
                        if lineFJ.refstrand ==line2.refstrand:
                            if int(lineFJ.loc_left) <= int(line2.loc_2):
                                IDfiletype = "FJgood,Regular,"+lineFJ.junction[-4:]
                    elif lineFJ.strand_left=="+":
                        if lineFJ.refstrand !=line2.refstrand:
                            if int(lineFJ.loc_left) >= int (line2.loc_2):
                                IDfiletype = "FJgood,Regular,"+lineFJ.junction[-4:]
            if IDfiletype == "FJbad,RegAnomaly,"+lineFJ.junction[-4:]:
                if lineFJ.chr_right == line2.chr:
                    if lineFJ.strand_right==line2.strand:
                        if lineFJ.strand_right=="-":
                            if lineFJ.refstrand ==line2.refstrand:
                                if int(lineFJ.loc_right) >= int(line2.loc_1):
                                    IDfiletype = "FJgood,Regular,"+lineFJ.junction[-4:]
                        elif lineFJ.strand_right=="+":
                            if lineFJ.refstrand !=line2.refstrand:
                                 if int(lineFJ.loc_right) <= int(line2.loc_1):
                                    IDfiletype = "FJgood,Regular,"+lineFJ.junction[-4:]
 
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")
        IDfile.flush()

    if inputtype == "genome": #comparing FJ to genome, has to be within 100Kbp, meet ref strand criteria (opp refstrand if + read, same refstrand if - read)

        line2 = ReadInfoGenome(line_raw_comparison)

        ## output R1 -  offset, MAPQ, AS, #N, readlen, junc name, strand
        IDfileoutputR1 =  str(lineFJ.offset) +"\t" + str(lineFJ.MAPQ) +"\t" + str(lineFJ.adjAS) + "\t" + lineFJ.NumN + "\t"+ str(lineFJ.NumOfBases) + "\t" +lineFJ.junction[:-5]+"\t"+lineFJ.refstrand        
         ## output R1 - offset, MAPQ, AS, #N, readlen, junc name, strand   
        IDfileoutputR2 = str(line2.loc) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.chr + "\t" + line2.refstrand


        IDfiletype = "FJbad,genomAnomaly,"+lineFJ.junction[-4:]   
        #compare left
        if lineFJ.chr_left == line2.chr:
            if lineFJ.strand_left=="-":
                if lineFJ.refstrand ==line2.refstrand:
                    if int(lineFJ.loc_left) <= int(line2.loc):
                        IDfiletype = "FJgood,genome,"+lineFJ.junction[-4:]                
            elif lineFJ.strand_left=="+":
                if lineFJ.refstrand != line2.refstrand:
                    if int(lineFJ.loc_left) >= int(line2.loc):                    
                        IDfiletype = "FJgood,genome,"+lineFJ.junction[-4:] 
        ## if left not the same, then compare right        
        if IDfiletype == "FJbad,genomAnomaly,"+lineFJ.junction[-4:]: 
            if lineFJ.chr_right == line2.chr:  
                if lineFJ.strand_right=="-":
                    if lineFJ.refstrand ==line2.refstrand:
                        if int(lineFJ.loc_right) >= int(line2.loc):
                            IDfiletype = "FJgood,genome,"+lineFJ.junction[-4:]                
                elif lineFJ.strand_right=="+":
                    if lineFJ.refstrand !=line2.refstrand:
                        if int(lineFJ.loc_right) <= int(line2.loc):
                            IDfiletype = "FJgood,genome,"+lineFJ.junction[-4:] 
            
    
        IDfile.write(line_raw_FJ.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")

#   
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
            self.adjAS=int(line[11].split(":")[2])+int(line[13][5:])
        else:
            self.NumN=line[12][5:]
            self.adjAS=int(line[11].split(":")[2])+int(line[12][5:])


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
parser.add_argument("-c", "--circReads", required= True, help = "path to circReads Dir")
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
if args.circReads[-1]!="/":
    args.circReads+="/"

stem = args.stem


FarJunctionfiles=[]
FarJunction_noIndelfiles=[]
genomefiles=[]
regfiles=[]
junctionfiles=[]

for name in glob.glob(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/*.sam"):
    print name
    if "All_" not in name:
        FarJunctionfiles.append(name)
    # FarJunctionFiles contains indel alignments for _1 and _2 files to indels 1-5


for name in glob.glob(args.FJDir+ "FarJunctionAlignments/"+ stem + "/*.sam" ):
    FarJunction_noIndelfiles.append(name)

for name in glob.glob(os.path.join(args.origDir,"genome/*" + stem + "*.sam")):        
#    print name
    if "sorted" not in name:       
        genomefiles.append(name)

for name in glob.glob(os.path.join(args.origDir,"reg/*" + stem + "*.sam")):
#    print name
    if "sorted" not in name:       
        regfiles.append(name) 
for name in glob.glob(os.path.join(args.origDir,"junction/*" + stem + "*.sam")):
#    print name
    if "sorted" not in name:       
        junctionfiles.append(name) 
#for name in glob.glob(os.path.join(args.origDir,"unaligned/*" + stem + "*.fq")):
##    print name
#    if "sorted" not in name:       
#        unalignedfiles.append(name)         


# opening all files for a particular stem
print sorted(FarJunction_noIndelfiles)
print sorted(genomefiles)
print sorted(regfiles)
print sorted(junctionfiles)

## open big indels files

f1_FarJunc= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/All_" + stem + "_1_indels.sam", mode ="rB")
f2_FarJunc= open(args.FJDir + "FarJuncSecondary/AlignedIndels/"+ stem + "/All_" + stem + "_2_indels.sam", mode="rB")

f1_FJ_noIndel=open(sorted(FarJunction_noIndelfiles)[0], mode="rB")
f2_FJ_noIndel=open(sorted(FarJunction_noIndelfiles)[1], mode="rB")


IDfile = open(args.FJDir+"GLM_classInput/"+ args.stem + "_temp_output_FJIndels.txt", mode= "w")
IDfile.write("ID\tclass\tR1_offset\tR1_MAPQ\tR1_AS\tR1_NumN\tR1_Readlength\tR1_JuncName\tR1_strand\tR2_offset\tR2_MAPQ\tR2_AS\tR2_NumN\tR2_Readlength\tR2_JuncName\tR2_strand\n")

#populate all reads and junctions into separate dictionaries
AllFJRead1= {}
AllFJRead2= {}


# populate AllFJRead1 dictionary - all read 1's from FarJunction alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening FarJunc _1 file"

for line_raw in f1_FarJunc:
    if line_raw[0] =="@":
        continue
    
    FJ1read = ReadInfoFJ(line_raw)
    
    if FJ1read.offset <= (150+FJ1read.indel-window) and (FJ1read.offset+FJ1read.NumOfBases)>= 150+FJ1read.indel+window:    
        AllFJRead1[FJ1read.ID] = [line_raw, 0]
        
f1_FarJunc.close()
IDfile.flush()

# populate AllFJRead2 dictionary - all read 2's from FarJunc alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening farJunc _2 file"
for line_raw in f2_FarJunc:
    if line_raw[0] =="@":
        continue
    
    FJ2read = ReadInfoFJ(line_raw)

    # if R1 and R2 both in Far Junc, then add to FJ-FJ list
    if FJ2read.ID in AllFJRead1:
        if FJ2read.offset <= (150+FJ2read.indel-window) and (FJ2read.offset+FJ2read.NumOfBases)>= 150+FJ2read.indel+window and AllFJRead1[FJ2read.ID][1]==0:    
            AddToDict("FJ",line_raw,AllFJRead1[FJ2read.ID][0])
            AllFJRead1[FJ2read.ID][1]="FJ"
    else:
        AllFJRead2[FJ2read.ID]= [line_raw, 0]
        
#    if FJ2read.junction not in AllJunctions:
#        AllJunctions[FJ2read.junction]=0
#        
#    AllJunctions[FJ2read.junction]+=1
f2_FarJunc.close()
IDfile.flush()

# compare FJ with indels_1 to FJ with no indels _ 2
print "comparing indels with FJ _2"
for line_raw in f2_FJ_noIndel:
    if line_raw[0] =="@":
        continue
    FJ2read = ReadInfoFJ(line_raw)
    
    if FJ2read.ID in AllFJRead1 and AllFJRead1[FJ2read.ID][1]==0:
        AddToDict("FJ", line_raw, AllFJRead1[FJ2read.ID][0])
        AllFJRead1[FJ2read.ID][1]="FJ"

f2_FJ_noIndel.close()    
IDfile.flush()


# compare FJ with indels _2 to FJ with no indels _1 

print "comparing indels with FJ _1"

for line_raw in f1_FJ_noIndel:
    if line_raw[0] =="@":
        continue
    FJ1read = ReadInfoFJ(line_raw)
    
    if FJ1read.ID in AllFJRead2 and AllFJRead2[FJ1read.ID][1]==0:
        #print "found genome R1"+g1read.ID      
#        if g1read.ID in unmappedDict:
#            del unmappedDict[g1read.ID]
        AddToDict("FJ",line_raw, AllFJRead2[FJ1read.ID][0])
        AllFJRead2[FJ1read.ID][1]="FJ"
f1_FJ_noIndel.close()    
IDfile.flush()




f2_genome = open(sorted(genomefiles)[1], mode="rB")

# compare FJ read 1 to genome read 2
print "comparing indels with genome_2"

for line_raw in f2_genome:
    if line_raw[0] =="@":
        continue
    g2read = ReadInfoGenome(line_raw)
    
    if g2read.ID in AllFJRead1 and AllFJRead1[g2read.ID][1]==0:
        AddToDict("genome", line_raw, AllFJRead1[g2read.ID][0])
        AllFJRead1[g2read.ID][1]="genome"

f2_genome.close()    

IDfile.flush()

f1_genome = open(sorted(genomefiles)[0], mode="rB")

# compare FJ read 2 to genome read 1
print "comparing indels with genome _1"
for line_raw in f1_genome:
    if line_raw[0] =="@":
        continue
    g1read = ReadInfoGenome(line_raw)
    
    if g1read.ID in AllFJRead2 and AllFJRead2[g1read.ID][1]==0:
        #print "found genome R1"+g1read.ID      
#        if g1read.ID in unmappedDict:
#            del unmappedDict[g1read.ID]
        AddToDict("genome", line_raw, AllFJRead2[g1read.ID][0])
        AllFJRead2[g1read.ID][1]="genome"
f1_genome.close()

IDfile.flush()

f2_reg = open(sorted(regfiles)[1], mode="rB")


# compare FJ read 1 to reg read 2
print "comparing indels with reg _2"

for line_raw in f2_reg:
    if line_raw[0] =="@":
        continue
    reg2read = ReadInfoJunc(line_raw)

    if reg2read.offset <= (150-window) and (reg2read.offset+reg2read.NumOfBases)>=( 150+window):    
        if reg2read.ID in AllFJRead1 and AllFJRead1[reg2read.ID][1]==0:
#            print "found reg R2:" + reg2read.ID
#            if reg2read.ID in unmappedDict:
#                del unmappedDict[reg2read.ID]
            AddToDict("reg", line_raw, AllFJRead1[reg2read.ID][0])
            AllFJRead1[reg2read.ID][1]="reg"
f2_reg.close()
IDfile.flush()



f1_reg = open(sorted(regfiles)[0], mode="rB")

# compare FJ read 2 to reg read 1
print "comparing indels with reg _1"

for line_raw in f1_reg:
    if line_raw[0] =="@":
        continue
    reg1read = ReadInfoJunc(line_raw)
    
    if reg1read.offset <= (150-window) and (reg1read.offset+reg1read.NumOfBases)>=( 150+window):    
        if reg1read.ID in AllFJRead2 and AllFJRead2[reg1read.ID][1]==0:
#            print "found reg R1: " + reg1read.ID
#            if reg1read.ID in unmappedDict:
#                del unmappedDict[reg1read.ID]
            AddToDict("reg", line_raw, AllFJRead2[reg1read.ID][0])
            AllFJRead2[reg1read.ID][1]="reg"
f1_reg.close()
IDfile.flush()


f2_junc= open(sorted(junctionfiles)[1], mode="rB")


# compare FJ read 1 to junc read 2
print "comparing indels with junc _2"

for line_raw in f2_junc:
    if line_raw[0] =="@":
        continue
    junc2read = ReadInfoJunc(line_raw)
    
    if junc2read.offset <= (150-window) and (junc2read.offset+junc2read.NumOfBases)>=( 150+window):
        if junc2read.ID in AllFJRead1 and AllFJRead1[junc2read.ID][1]==0:
            #print "found junc R2 " + junc2read.ID
#            if junc2read.ID in unmappedDict:
#                del unmappedDict[junc2read.ID]
            AddToDict("junc", line_raw, AllFJRead1[junc2read.ID][0])
            AllFJRead1[junc2read.ID][1]="junc"
f2_junc.close()
IDfile.flush()


f1_junc= open(sorted(junctionfiles)[0], mode="rB")

# compare FJ read 2 to junc read 1
print "comparing indels with junc _1"

for line_raw in f1_junc:
    if line_raw[0] =="@":
        continue
    junc1read = ReadInfoJunc(line_raw)

    if junc1read.offset <= (150-window) and (junc1read.offset+junc1read.NumOfBases)>= (150+window):
        if junc1read.ID in AllFJRead2 and AllFJRead2[junc1read.ID][1]==0:
            #print "found junc R1: " + junc1read.ID
#            if junc1read.ID in unmappedDict:
#                del unmappedDict[junc1read.ID]            
            AddToDict("junc", line_raw, AllFJRead2[junc1read.ID][0])
            AllFJRead2[junc1read.ID][1]="junc"
f1_junc.close()
IDfile.flush()

IDfile.close()


#############################################################################
## This section of code takes the written ID file above (temp_IDs_STEM.txt) and
## removes duplicate entries of genome and reg. The same readID may be found 
## in both libraries and would both be in the ID file.
## The new ID file removes duplicates and only keeps the readID with the 
## best alignment score.
tempIDfile = open(args.FJDir+"GLM_classInput/"+ args.stem + "_temp_output_FJIndels.txt", mode= "rU")
newIDfile = open(args.FJDir+"GLM_classInput/"+ args.stem + "_output_FJIndels.txt", mode= "w")

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
