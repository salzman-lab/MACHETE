# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 16:12:27 2015

@author: Gillian
"""
# takes alignments from Indels and finds read partner. 
# all Far Junction and Scrambled junction reads where the alignment does not overlap the 
# Junction by args.window bp will be thrown out.


################
#Current categories
#linear -- genomeGood, RegGood, JuncGood, RegIndelGood 
#anomaly -- genomeBad, RegBad, JuncBad, RegIndelBad
##################



import argparse
import os
import glob

def AddToDict(inputtype, line_raw_comparison, line_raw_RI):

    lineRI = ReadInfoRI(line_raw_RI)

 
    if inputtype=="RI": # if comparing Far Junc to Far Junc, they have to be identical
        line2= ReadInfoRI(line_raw_comparison)
        
        IDfileoutputR1 =  str(lineRI.offset) +"\t" + str(lineRI.MAPQ) +"\t" + str(lineRI.adjAS) + "\t" + lineRI.NumN + "\t"+ str(lineRI.NumOfBases) + "\t" +lineRI.junction[:-5]+"\t"+lineRI.refstrand        
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand

        if lineRI.junction == line2.junction and lineRI.refstrand in ["0","16"] and line2.refstrand in ["0","16"] and lineRI.refstrand!=line2.refstrand:
            IDfiletype = "linear,RegIndelGood,"+lineRI.junction[-4:]
        else:
            IDfiletype = "anomaly,RegIndelBad,"+lineRI.junction[-4:]
 
        IDfile.write(line_raw_RI.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype=="reg":
        line2= ReadInfoJunc(line_raw_comparison)
        
        IDfileoutputR1 =  str(lineRI.offset) +"\t" + str(lineRI.MAPQ) +"\t" + str(lineRI.adjAS) + "\t" + lineRI.NumN + "\t"+ str(lineRI.NumOfBases) + "\t" +lineRI.junction[:-5]+"\t"+lineRI.refstrand        
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand
        
        IDfiletype = "linear,RegGood,"+lineRI.junction[-4:]
        
        if lineRI.chr==line2.chr and lineRI.refstrand!=line2.refstrand and lineRI.strand==line2.strand:
            if lineRI.strand=="+":
                if line2.loc_right < lineRI.loc_left or line2.loc_left > lineRI.loc_right:
                    pass
                else:
                    IDfiletype="anomaly,RegBad,"+lineRI.junction[-4:]
            elif lineRI.strand=="-": 
                if line2.loc_right > lineRI.loc_left or line2.loc_left < lineRI.loc_right:
                    pass
                else:
                    IDfiletype="anomaly,RegBad,"+lineRI.junction[-4:]
        else:
            IDfiletype="anomaly,RegBad,"+lineRI.junction[-4:]

        IDfile.write(line_raw_RI.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype =="junc": #if reg or junc read, then one side has to be within 100KB, and meets refstrand criteria below
        
        line2 = ReadInfoJunc(line_raw_comparison)

        IDfileoutputR1 =  str(lineRI.offset) +"\t" + str(lineRI.MAPQ) +"\t" + str(lineRI.adjAS) + "\t" + lineRI.NumN + "\t"+ str(lineRI.NumOfBases) + "\t" +lineRI.junction[:-5]+"\t"+lineRI.refstrand        
        IDfileoutputR2 = str(line2.offset) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.junction + "\t" + line2.refstrand

        IDfiletype = "linear,JuncGood,"+lineRI.junction[-4:]
        
        if lineRI.chr==line2.chr and lineRI.refstrand!=line2.refstrand and lineRI.strand==line2.strand and min(line2.loc_left, line2.loc_right)< lineRI.loc_left and min(line2.loc_left, line2.loc_right)< lineRI.loc_right and max (line2.loc_left, line2.loc_right) > lineRI.loc_left and max(line2.loc_left, line2.loc_right)> lineRI.loc_right:
            pass
        else:
            IDfiletype="anomaly,JuncBad,"+lineRI.junction[-4:]
    
        IDfile.write(line_raw_RI.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")


    if inputtype == "genome": #comparing FJ to genome, has to be within 100Kbp, meet ref strand criteria (opp refstrand if + read, same refstrand if - read)

        line2 = ReadInfoGenome(line_raw_comparison)

        IDfileoutputR1 =  str(lineRI.offset) +"\t" + str(lineRI.MAPQ) +"\t" + str(lineRI.adjAS) + "\t" + lineRI.NumN + "\t"+ str(lineRI.NumOfBases) + "\t" +lineRI.junction[:-5]+"\t"+lineRI.refstrand        
        IDfileoutputR2 = str(line2.loc) + "\t" + str(line2.MAPQ) + "\t"+ str(line2.adjAS) + "\t"+ line2.NumN + "\t" + str(line2.NumOfBases) + "\t" + line2.chr + "\t" + line2.refstrand

        IDfiletype = "linear,genomeGood,"+lineRI.junction[-4:]

        if lineRI.chr==line2.chr and lineRI.refstrand!=line2.refstrand and lineRI.strand=="+":
            if line2.loc>lineRI.loc_right or line2.loc<lineRI.loc_left:
                pass
            else: 
                IDfiletype = "anomaly,genomeBad,"+lineRI.junction[-4:]
        elif lineRI.chr==line2.chr and lineRI.refstrand!=line2.refstrand and lineRI.strand=="-":    
            if line2.loc>lineRI.loc_left or line2.loc > lineRI.loc_right:
                pass
            else:
                IDfiletype = "anomaly,genomeBad,"+lineRI.junction[-4:]
        else:
            IDfiletype = "anomaly,genomeBad,"+lineRI.junction[-4:]

        IDfile.write(line_raw_RI.split("\t")[0]+"\t"+IDfiletype+"\t"+IDfileoutputR1+"\t"+IDfileoutputR2+"\n")
    
     
def ID(string):
    if string[-2:] == "/1" or string[-2:] == "/2":
        return string[:-2]
    else:
        return string



class ReadInfoRI:
    def __init__ (self,line_raw):
        line = line_raw.strip().split("\t")
        self.ID = ID(line[0])       
        self.refstrand = line[1]
        self.junction = line[2]
        if line[2][-4:-1]=="DEL":
#            print "is a deletion"
#            print line[2][-1:]
            self.indel=-int(line[2][-1:])
        
        if line[2][-4:-1]=="INS":
            self.indel=int(line[2][-1:])
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
        self.chr=JuncInfo[0]
        self.loc_left=int(JuncInfo[2])
        self.loc_right = int(JuncInfo[4])
        self.strand = JuncInfo[6]
        
class ReadInfoGenome:
    def __init__ (self, line_raw):
        line = line_raw.strip().split("\t")
        self.ID= ID(line[0])
        self.refstrand = line[1]
        self.chr = line[2]
        self.loc = int(line[3])
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
        self.loc_left = int(line[2].replace(":","|").split("|")[2])
        self.loc_right = int(line[2].replace(":","|").split("|")[4])
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
parser.add_argument("-i", "--origDir", required=True, help = "path to orig dir containing genome reads")
parser.add_argument("-w", "--window", required=True, help = "# of bases needed on each side of the junction")
args = parser.parse_args()
window= int(args.window)

# f1 = open("/Users/Gillian/Desktop/sherlock/unaligned_ENCFF000HOC1_1.sam", mode ="rU")
# f2 = open("/Users/Gillian/Desktop/sherlock/20000_ENCFF000HOC2_1_genome_output.sam", mode ="rU")



if args.origDir[-1]!="/":
    args.origDir+="/"
if args.circReads[-1]!="/":
    args.circReads+="/"


stem = args.stem

regIndelfiles=[]
regfiles=[]
genomefiles=[]
junctionfiles=[]

for name in glob.glob(args.origDir + "RegIndelAlignments/"+ stem + "/*.sam"):
    print name
    if "All_" not in name:       
        regIndelfiles.append(name) 
    # Regfiles contains indel alignments for _1 and _2 files to indels 1-5


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


# opening all files for a particular stem
#print sorted(regIndelfiles)
#print sorted(genomefiles)
#print sorted(regfiles)
#print sorted(junctionfiles)


## concatenate all reg indels files into a single big indels file

regIndel1_list = sorted(regIndelfiles)[0:len(regIndelfiles)/2]
regIndel2_list = sorted(regIndelfiles)[len(regIndelfiles)/2:]


IndelsReadIDs={}

for name in regIndel1_list:
    print "reg1 indels"
    print name
    f1= open(name, mode="rU")
    
    for line in f1:
        if line[0]=="@":
            continue
        
        read = ReadInfoRI(line)
        # if the read overlaps the junction
        if read.offset <= (150-int(args.window)+read.indel) and read.offset+read.NumOfBases >= (150+int(args.window)+read.indel):

            # if the read isn't in dictionary then add it 
            if read.ID not in IndelsReadIDs:             
                IndelsReadIDs[read.ID]= line
            # if the read is in the dictionary, compare it to existing read. If AS is better, then replace existing read
            else:
                compareRead= ReadInfoRI(IndelsReadIDs[read.ID])
                if int(compareRead.AS) >= int(read.AS):
                    pass
                else:
                    IndelsReadIDs[read.ID]= line                            
    f1.close()
  
# write all distinct readIDs to an All_1_indels file  
fout_regIndel1= open(args.origDir + "RegIndelAlignments/"+ stem + "/All_" + stem + "_1_Regindels.sam", mode="w")

for key in IndelsReadIDs:
    fout_regIndel1.write(IndelsReadIDs[key].strip()+"\n")
    
fout_regIndel1.close()
    

## CLEAR Read IDs dictionary and do the same with FJ2 list
IndelsReadIDs={}

for name in regIndel2_list:
    print "regIndel2 indels"
    print name    
    
    f1= open(name, mode="rU")
    
    for line in f1:
        if line[0]=="@":
            continue
        
        read = ReadInfoRI(line)
        # if the read overlaps the junction
        if read.offset <= (150-int(args.window)+read.indel) and read.offset+read.NumOfBases >= (150+int(args.window)+read.indel):

            # if the read isn't in dictionary then add it 
            if read.ID not in IndelsReadIDs:             
                IndelsReadIDs[read.ID]= line
            # if the read is in the dictionary, compare it to existing read. If AS is better, then replace existing read
            else:
                compareRead= ReadInfoRI(IndelsReadIDs[read.ID])
                if int(compareRead.AS) >= int(read.AS):
                    pass
                else:
                    IndelsReadIDs[read.ID]= line                            
    f1.close()
  

# write all distinct readIDs to an All_1_indels file  
fout_regIndel2= open(args.origDir + "RegIndelAlignments/"+ stem + "/All_" + stem + "_2_Regindels.sam", mode="w")

for key in IndelsReadIDs:
    fout_regIndel2.write(IndelsReadIDs[key].strip()+"\n")
fout_regIndel2.close()


## open big indels files

f1_regIndel= open(args.origDir + "RegIndelAlignments/"+ stem + "/All_" + stem + "_1_Regindels.sam", mode ="rB")
f2_regIndel= open(args.origDir + "RegIndelAlignments/"+ stem + "/All_" + stem + "_2_Regindels.sam", mode="rB")



# ID file ReadID and different buckets.
    # [0] = readID
    # [1] = R2 in genome        
    # [2] = R2 in genome anomaly  
    # [3] = reg
    # [4] = reg anom
    # [5] = junc
    # [6] = junc anom



IDfile = open(args.circReads+"ids/"+stem+"_temp_output_RegIndel.txt", mode= "w")
IDfile.write("ID\tclass\tR1_offset\tR1_MAPQ\tR1_adjAS\tR1_NumN\tR1_Readlength\tR1_JuncName\tR1_strand\tR2_offset\tR2_MAPQ\tR2_adjAS\tR2_NumN\tR2_Readlength\tR2_JuncName\tR2_strand\n")

#populate all reads and junctions into separate dictionaries
AllRegIndelRead1= {}
AllRegIndelRead2= {}
#AllJunctions = {}
#genomeDict = {}  # for all these dictionaries, [0] = reg, [1] = anom
#regDict = {}        # [2] = sum of AS, [3] = read length
#juncDict = {}
#unmappedDict= {} # start with all readIDs.  if a partner is seen, then remove from list.


# populate AllFJRead1 dictionary - all read 1's from FarJunction alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening All_regIndels_1 file"

for line_raw in f1_regIndel:
    if line_raw[0] =="@":
        continue
    
    RegIndel1read = ReadInfoRI(line_raw)
    
    if RegIndel1read.offset <= (150+RegIndel1read.indel-window) and (RegIndel1read.offset+RegIndel1read.NumOfBases)>= 150+RegIndel1read.indel+window:    
        AllRegIndelRead1[RegIndel1read.ID] = [line_raw, 0]
#        if FJ1read.junction not in AllJunctions:
#            AllJunctions[FJ1read.junction]=0
#        AllJunctions[FJ1read.junction] +=1
#        unmappedDict[AllRegRead1.ID] = AllRegRead1.junction
        
f1_regIndel.close()
IDfile.flush()

# populate AllFJRead2 dictionary - all read 2's from FarJunc alignments
# in order for R1 to feed into dictionary, must overlap entire offset (userspecified)
print "opening All_regIndels _2 file"
for line_raw in f2_regIndel:
    if line_raw[0] =="@":
        continue
    
    RegIndel2read = ReadInfoRI(line_raw)

#    if FJ1read.junction=="chr1:S100A4:153516097:-|chr1:IFI16:158985661:+|strandcross":
#        print "ERROR AT LINE 409"

    # if R1 and R2 both in Far Junc, then add to FJ-FJ list
    if RegIndel2read.ID in AllRegIndelRead1:
        #print "found FJ read"
        #AllFJRead1[FJ2read.ID][1]="FJ"
        if RegIndel2read.offset <= (150+RegIndel2read.indel-window) and (RegIndel2read.offset+RegIndel2read.NumOfBases)>= 150+RegIndel2read.indel+window and AllRegIndelRead1[RegIndel2read.ID][1]==0:    
            AddToDict("RI",line_raw,AllRegIndelRead1[RegIndel2read.ID][0])
            AllRegIndelRead1[RegIndel2read.ID][1]="RI"
#            if FJ2read.ID in unmappedDict:
#               del unmappedDict[FJ2read.ID]
        # otherwise add to F2 read
    else:
        AllRegIndelRead2[RegIndel2read.ID]= [line_raw, 0]
#        unmappedDict[FJ2read.ID] = FJ2read.junction
        
#    if RegIndel2read.junction not in AllJunctions:
#        AllJunctions[FJ2read.junction]=0
        
#    AllJunctions[FJ2read.junction]+=1
f2_regIndel.close()
IDfile.flush()



f2_genome = open(sorted(genomefiles)[1], mode="rB")


# compare FJ read 1 to genome read 2
print "comparing Indels to genome_2"
for line_raw in f2_genome:
    if line_raw[0] =="@":
        continue
    g2read = ReadInfoGenome(line_raw)
    
    if g2read.ID in AllRegIndelRead1 and AllRegIndelRead1[g2read.ID][1]==0:
        #print "found genome R2"+ g2read.ID     
#        if g2read.ID in unmappedDict:
#            del unmappedDict[g2read.ID]
        AddToDict("genome", line_raw, AllRegIndelRead1[g2read.ID][0])
        AllRegIndelRead1[g2read.ID][1]="genome"

f2_genome.close()    
IDfile.flush()

f1_genome = open(sorted(genomefiles)[0], mode="rB")

# compare FJ read 2 to genome read 1
print "comparing Indels to genome_1"

for line_raw in f1_genome:
    if line_raw[0] =="@":
        continue
    g1read = ReadInfoGenome(line_raw)
    
    if g1read.ID in AllRegIndelRead2 and AllRegIndelRead2[g1read.ID][1]==0:
        #print "found genome R1"+g1read.ID      
#        if g1read.ID in unmappedDict:
#            del unmappedDict[g1read.ID]
        AddToDict("genome", line_raw, AllRegIndelRead2[g1read.ID][0])
        AllRegIndelRead2[g1read.ID][1]="genome"
f1_genome.close()    
IDfile.flush()
#


f2_reg = open(sorted(regfiles)[1], mode="rB")

# compare FJ read 1 to reg read 2
print "comparing Indels to reg_2"

for line_raw in f2_reg:
    if line_raw[0] =="@":
        continue
    reg2read = ReadInfoJunc(line_raw)

    if reg2read.offset <= (150-window) and (reg2read.offset+reg2read.NumOfBases)>=( 150+window):    
        if reg2read.ID in AllRegIndelRead1 and AllRegIndelRead1[reg2read.ID][1]==0:
            #print "found reg R2: " + reg2read.ID
#            if reg2read.ID in unmappedDict:
#                del unmappedDict[reg2read.ID]
            AddToDict("reg", line_raw, AllRegIndelRead1[reg2read.ID][0])
            AllRegIndelRead1[reg2read.ID][1]="reg"
f2_reg.close()
IDfile.flush()

f1_reg = open(sorted(regfiles)[0], mode="rB")


# compare FJ read 2 to reg read 1
print "comparing Indels to reg_1"

for line_raw in f1_reg:
    if line_raw[0] =="@":
        continue
    reg1read = ReadInfoJunc(line_raw)
    
    if reg1read.offset <= (150-window) and (reg1read.offset+reg1read.NumOfBases)>=( 150+window):    
        if reg1read.ID in AllRegIndelRead2 and AllRegIndelRead2[reg1read.ID][1]==0:
            #print "found reg R1: " + reg1read.ID
#            if reg1read.ID in unmappedDict:
#                del unmappedDict[reg1read.ID]
            AddToDict("reg", line_raw, AllRegIndelRead2[reg1read.ID][0])
            AllRegIndelRead2[reg1read.ID][1]="reg"
f1_reg.close()
IDfile.flush()



f2_junc= open(sorted(junctionfiles)[1], mode="rB")
     
# compare FJ read 1 to junc read 2
print "comparing Indels to junc_2"

for line_raw in f2_junc:
    if line_raw[0] =="@":
        continue
    junc2read = ReadInfoJunc(line_raw)
    
    if junc2read.offset <= (150-window) and (junc2read.offset+junc2read.NumOfBases)>=( 150+window):
        if junc2read.ID in AllRegIndelRead1 and AllRegIndelRead1[junc2read.ID][1]==0:
#            print "found junc R2: " + junc2read.ID
#            if junc2read.ID in unmappedDict:
#                del unmappedDict[junc2read.ID]
            AddToDict("junc", line_raw, AllRegIndelRead1[junc2read.ID][0])
            AllRegIndelRead1[junc2read.ID][1]="junc"
f2_junc.close()
IDfile.flush()


f1_junc= open(sorted(junctionfiles)[0], mode="rB")

# compare FJ read 2 to junc read 1
print "comparing Indels to junc_1"

for line_raw in f1_junc:
    if line_raw[0] =="@":
        continue
    junc1read = ReadInfoJunc(line_raw)

    if junc1read.offset <= (150-window) and (junc1read.offset+junc1read.NumOfBases)>= (150+window):
        if junc1read.ID in AllRegIndelRead2 and AllRegIndelRead2[junc1read.ID][1]==0:
#            print "found junc R1: " + junc1read.ID
#            if junc1read.ID in unmappedDict:
#                del unmappedDict[junc1read.ID]            
            AddToDict("junc", line_raw, AllRegIndelRead2[junc1read.ID][0])
            AllRegIndelRead2[junc1read.ID][1]="junc"
f1_junc.close()
IDfile.flush()

IDfile.close()

#


#############################################################################
## This section of code takes the written ID file above (temp_IDs_STEM.txt) and
## removes duplicate entries of genome and reg. The same readID may be found 
## in both libraries and would both be in the ID file.
## The new ID file removes duplicates and only keeps the readID with the 
## best alignment score.
tempIDfile = open(args.circReads+"ids/"+stem+"_temp_output_RegIndel.txt", mode= "rU")
newIDfile = open(args.circReads+"ids/"+stem+"_output_RegIndel.txt", mode= "w")



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
