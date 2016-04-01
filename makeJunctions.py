# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 10:24:15 2015

@author: Gillian
"""
from itertools import islice
import os
import glob
import cPickle as pickle
import re
from Bio import SeqIO
import time
import argparse
import sys

class locationparameters:
    def __init__(self,read):
        if read[-1] == "-":
            self.strand="-"
        else:
            self.strand="+"
        
        read = read.replace(":"," ").replace("-"," ").replace("strand=", " ").split(" ")
        self.chr = read[0]
        self.start = int(read[1])
        self.stop = int(read[2])
        
        
# make a list of 26 dictionaries with each bin representing a different chromosome
    # each dictionary contains key, value with UI, start/stop tuple respectively

def SortedLocation(UI, Chr, Start, Stop):
    try: ChrBin = int(Chr)-1
    except: 
        x = {"X":23, "Y":24, "M":25}
        ChrBin = x[Chr]
    # each chromosome has a dictionary. Add location of read to chromosome dictionary
#    print ChrBin
    try: 
        LocationTuples_AllChr[ChrBin][UI]=(Start,Stop)
    except:
        LocationTuples_AllChr[ChrBin]={UI:(Start,Stop)}
#    print LocationTuples_AllChr[ChrBin]
        

def ChrAdjuster(chrname):
    return 200000*(int(chrname.replace("chr","").replace("X","23").replace("Y","24").replace("M","25")) -1)

def reverseChrAdjuster(input):
    x = {22:"X", 23:"Y", 24:"M" }
    if int(input)>=23:
        chrname = x[input]
    else:
        chrname = str(int+1)
    return "chr"+chrname
    
def StrandAdjuster(strand):
    if strand == 1:
        return 0
    if strand ==-1:
        return 5000000
        
def UnpickleExons(LocationDict,Picklefile,chrID): 
    AllExonsPos = []
    AllExonsNeg = []
    
    if LocationDict==0:
        return AllExonsPos,AllExonsNeg

    ExonPickle = pickle.load(Picklefile)

    for read in ExonPickle:
        chrID, strand, exons = read
        ExonsInRange = exons.keys()
        if strand == 1:
#            print "unpickling pos strand"
            AllExonsPos = exons
        else:
#            print "unpickling neg strand"
            AllExonsNeg = exons
        
        if len(ExonsInRange)>0:
            for index in LocationDict:
                start, stop = LocationDict[index]
#                print index + ":" + str(start)+ "-" + str(stop)
                if index[-1]=="A":
                    ExonsListA[int(index[:-1])+ChrAdjuster(chrID)+StrandAdjuster(strand)]= [x for x in ExonsInRange if x[0]>= start and x[1]<= stop]
#                    print "# exons going into list A"
#                    print len(ExonsListA[int(index[:-1])+ChrAdjuster(chrID)+StrandAdjuster(strand)])
#                    print int(index[:-1])+ChrAdjuster(chrID)+StrandAdjuster(strand)
                if index[-1]=="B":
                    ExonsListB[int(index[:-1])+ChrAdjuster(chrID)+StrandAdjuster(strand)]= [x for x in ExonsInRange if x[0] >= start and x[1] <= stop]
#                    print "# exons going into list B"
#                    print len(ExonsListB[int(index[:-1])+ChrAdjuster(chrID)+StrandAdjuster(strand)])
#                    print int(index[:-1])+ChrAdjuster(chrID)+StrandAdjuster(strand)
#    print "A plus strand" 
#    print ExonsListA[ChrAdjuster(chrID):UI_counter+ChrAdjuster(chrID)]
#    print "A minus strand"
#    print ExonsListA[ChrAdjuster(chrID)+StrandAdjuster(-1):UI_counter+ChrAdjuster(chrID)+StrandAdjuster(-1)]
#    print "B plus strand" 
#    print ExonsListB[ChrAdjuster(chrID):UI_counter+ChrAdjuster(chrID)]
#    print "B minus strand"
#    print ExonsListB[ChrAdjuster(chrID)+StrandAdjuster(-1):UI_counter+ChrAdjuster(chrID)+StrandAdjuster(-1)]
    Picklefile.close()    
    
#    print "All Exons Pos:" 
#    print AllExonsPos
#    print "All Exons Neg:" 
#    print AllExonsNeg
    return AllExonsPos, AllExonsNeg
#
    
        
def UnpickleSequence(chrnum, PickleRecfile, AllExonsPos, AllExonsNeg):

    exonSeqRec = pickle.load(PickleRecfile)    
    
    # Exons A, + strand
    counter = 0
    for x in ExonsListA[ChrAdjuster(chrnum):ChrAdjuster(chrnum)+UI_counter]:
        if x!=0 and len(x)>0:       
            for i in x:
                #print i
                if NameListA[ChrAdjuster(chrnum)+counter]==0:
                    NameListA[ChrAdjuster(chrnum)+counter] = []
                if SeqListA[ChrAdjuster(chrnum)+counter]==0:
                    SeqListA[ChrAdjuster(chrnum)+counter] = []
                
                NameListA[ChrAdjuster(chrnum)+counter]+= [AllExonsPos[i].qualifiers["gene_id"][0]]
                SeqListA[ChrAdjuster(chrnum)+counter]+= [AllExonsPos[i].extract(exonSeqRec.seq)]
                                 
                #print SeqListA[ChrAdjuster(chrnum)+counter]
        counter +=1
    # Exons A, - strand
    counter = 0
    for x in ExonsListA[ChrAdjuster(chrnum)+StrandAdjuster(-1):ChrAdjuster(chrnum)+UI_counter+StrandAdjuster(-1)]:
        if x!=0 and len(x)>0:       
            for i in x:
                #print i
                if NameListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]==0:
                    NameListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)] = []
                if SeqListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]==0:
                    SeqListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)] = []
            
                NameListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]+= [AllExonsNeg[i].qualifiers["gene_id"][0]]
                SeqListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]+= [AllExonsNeg[i].extract(exonSeqRec.seq)]
                       
                #print SeqListA[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]
        counter +=1
    
#    # Exons B, + strand    
    counter = 0
    for x in ExonsListB[ChrAdjuster(chrnum):ChrAdjuster(chrnum)+UI_counter]:
        if x!=0 and len(x)>0:       
            for i in x:

                if NameListB[ChrAdjuster(chrnum)+counter]==0:
                    NameListB[ChrAdjuster(chrnum)+counter] = []
                if SeqListB[ChrAdjuster(chrnum)+counter]==0:
                    SeqListB[ChrAdjuster(chrnum)+counter] = []
                
                NameListB[ChrAdjuster(chrnum)+counter]+= [AllExonsPos[i].qualifiers["gene_id"][0]]
                SeqListB[ChrAdjuster(chrnum)+counter]+= [AllExonsPos[i].extract(exonSeqRec.seq)]
        counter +=1
#
#    # exons B, - strand    
    counter = 0
    for x in ExonsListB[ChrAdjuster(chrnum)+StrandAdjuster(-1):ChrAdjuster(chrnum)+UI_counter+StrandAdjuster(-1)]:
        if x!=0 and len(x)>0:       
            for i in x:
                if NameListB[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]==0:
                    NameListB[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)] = []
                if SeqListB[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]==0:
                    SeqListB[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)] = []
            
                NameListB[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]+= [AllExonsNeg[i].qualifiers["gene_id"][0]]
                SeqListB[ChrAdjuster(chrnum)+counter+StrandAdjuster(-1)]+= [AllExonsNeg[i].extract(exonSeqRec.seq)]
        counter +=1
    PickleRecfile.close()
        
        
def JunctionEmpty(ExonsA, ExonsB):
    if ExonsA in [0,[]] or ExonsB in [0,[]]:
        return True
    else:
        return False


def MakeJunc(chrA, strandA, NameA, LocA, chrB, strandB, NameB, LocB):
    x= {1:"+", -1:"-"}
    if chrA!=chrB:
        JuncType = "fusion"
    elif strandA != strandB:
        JuncType = "strandcross"
    elif strandA==1 and strandB==1:
        if int(LocA) < int(LocB):
            JuncType = "reg"
        else:
            JuncType = "rev"
    else:
        if int(LocA) > int(LocB):
            JuncType = "reg"
        else: 
            JuncType = "rev"

    trueLocA=int(LocA)
    trueLocB=int(LocB)

    if strandA==-1:
        trueLocA+=1
    if strandB==1:
        trueLocB+=1

    return "chr"+str(chrA)+":"+str(NameA)+":"+str(trueLocA)+":"+x[strandA]+"|chr"+str(chrB)+":"+str(NameB)+":"+str(trueLocB)+":"+x[strandB]+"|"+JuncType
            
def MakeSeq(seqA, seqB):
    pad = "N"*150
    boundary = 150+len(seqA)
    Sequence = pad + seqA + seqB + pad
    return Sequence[boundary-150:boundary+150]


#def AddToLib(AllJunctions, Junction, Sequence, outputfile):
#    if Junction not in AllJunctions:
#        AllJunctions.append(Junction)
#        outputfile.write(">"+str(Junction)+"\n"+str(Sequence)+"\n")
#        return AllJunctions
#    else:
#        return AllJunctions


def MakePairs(outputfile, AllJunctions, UI, chrA, chrB, ExonsApos, ExonsAneg, ExonsBpos, ExonsBneg, SeqApos, SeqAneg, SeqBpos, SeqBneg, NameApos, NameAneg, NameBpos, NameBneg): 
    
#    print "MakePairs"
#    print ExonsApos
#    print ExonsAneg
#    print ExonsBpos
#    print ExonsBneg
#    print NameApos
#    print NameAneg
#    print NameBpos
#    print NameBneg
#    
    # case A pos, B pos 
    if not JunctionEmpty(ExonsApos, ExonsBpos):
#        print "Making all pairs for A pos, B pos"
        for i in range(0, len(ExonsApos)):
            for j in range(0,len(ExonsBpos)):
               ABJunc= MakeJunc(chrA, 1, NameApos[i], ExonsApos[i][1], chrB, 1, NameBpos[j], ExonsBpos[j][0])
               BAJunc = MakeJunc(chrB, 1, NameBpos[j], ExonsBpos[j][1], chrA, 1, NameApos[i], ExonsApos[i][0])
               ABSeq= MakeSeq(SeqApos[i], SeqBpos[j])
               BASeq = MakeSeq(SeqBpos[j], SeqApos[i])
#               AllJunctions = AddToLib(AllJunctions, ABJunc, ABSeq, outputfile)
#               AllJunctions = AddToLib(AllJunctions, BAJunc, BASeq, outputfile)
               outputfile.write(">"+str(ABJunc)+"\n"+str(ABSeq)+"\n")
               outputfile.write(">"+str(BAJunc)+"\n"+str(BASeq)+"\n")

  
    #CASE A pos and B neg   
    if not JunctionEmpty(ExonsApos, ExonsBneg):
#        print "Making all pairs for A pos, B neg"
        for i in range(0, len(ExonsApos)):
            for j in range(0,len(ExonsBneg)):
               ABJunc= MakeJunc(chrA, 1, NameApos[i], ExonsApos[i][1], chrB, -1, NameBneg[j], ExonsBneg[j][1])
               BAJunc = MakeJunc(chrB, -1, NameBneg[j], ExonsBneg[j][0], chrA, 1, NameApos[i], ExonsApos[i][0])
               ABSeq= MakeSeq(SeqApos[i], SeqBneg[j])
               BASeq = MakeSeq(SeqBneg[j], SeqApos[i])
#               AllJunctions = AddToLib(AllJunctions, ABJunc, ABSeq, outputfile)
#               AllJunctions = AddToLib(AllJunctions, BAJunc, BASeq, outputfile)
               outputfile.write(">"+str(ABJunc)+"\n"+str(ABSeq)+"\n")
               outputfile.write(">"+str(BAJunc)+"\n"+str(BASeq)+"\n")             

    #CASE A neg and B pos   
    if not JunctionEmpty(ExonsAneg, ExonsBpos):
#        print "Making all pairs for A neg, B pos"
        for i in range(0, len(ExonsAneg)):
            for j in range(0,len(ExonsBpos)):
               ABJunc= MakeJunc(chrA, -1, NameAneg[i], ExonsAneg[i][0], chrB, 1, NameBpos[j], ExonsBpos[j][0])
               BAJunc = MakeJunc(chrB, 1, NameBpos[j], ExonsBpos[j][1], chrA, -1, NameAneg[i], ExonsAneg[i][1])
               ABSeq= MakeSeq(SeqAneg[i], SeqBpos[j])
               BASeq = MakeSeq(SeqBpos[j], SeqAneg[i])
#               AllJunctions = AddToLib(AllJunctions, ABJunc, ABSeq, outputfile)
#               AllJunctions = AddToLib(AllJunctions, BAJunc, BASeq, outputfile)
               outputfile.write(">"+str(ABJunc)+"\n"+str(ABSeq)+"\n")
               outputfile.write(">"+str(BAJunc)+"\n"+str(BASeq)+"\n")             

    #CASE A neg and B neg   
    if not JunctionEmpty(ExonsAneg, ExonsBneg):
#        print "Making all pairs for A neg, B neg"
        for i in range(0, len(ExonsAneg)):
            for j in range(0,len(ExonsBneg)):
               ABJunc= MakeJunc(chrA, -1, NameAneg[i], ExonsAneg[i][0], chrB, -1, NameBneg[j], ExonsBneg[j][1])
               BAJunc = MakeJunc(chrB, -1, NameBneg[j], ExonsBneg[j][0], chrA, -1, NameAneg[i], ExonsAneg[i][1])
               ABSeq= MakeSeq(SeqAneg[i], SeqBneg[j])
               BASeq = MakeSeq(SeqBneg[j], SeqAneg[i])
#               AllJunctions = AddToLib(AllJunctions, ABJunc, ABSeq, outputfile)
#               AllJunctions = AddToLib(AllJunctions, BAJunc, BASeq, outputfile)
               outputfile.write(">"+str(ABJunc)+"\n"+str(ABSeq)+"\n")
               outputfile.write(">"+str(BAJunc)+"\n"+str(BASeq)+"\n")             



                
# ==================PROGRAM START ==============



parser=argparse.ArgumentParser()
parser.add_argument("-p", "--pickle", required = True, help = "path to pickle directory")
parser.add_argument("-f", "--infile", required= True, help = "file to unpickle")
parser.add_argument("-o", "--outDir", required=True, help = " directory to output files")
parser.add_argument("-s", "--stem", required=True, help = "unique identifying stem of file name")
args=parser.parse_args()



if args.pickle[-1] != "/":
    path = args.pickle+"/"
else:
    path = args.pickle
    
        
if args.outDir[-1] != "/":
    outDir = args.outDir + "/"
else:
    outDir = args.outDir
    

##TEST ENVIRONMENT
#path = "/Users/Gillian/Desktop/pickles/"
#fastaDir = "/Users/Gillian/Desktop/ERP000710output/"    
    

infilepath, infilename = os.path.split(args.infile)
chromosome=infilename.split("_")[1]



os.chdir(path)
exonDir = path + "exons"
recordDir = path+"records"
EMPTY_SEQUENCE = "NOSEQUENCE"
patt_exonfilename = re.compile(".+?exonsByStrand_(.+?)\.pkl")
patt_exonfile = re.compile("exonsByStrand_.+\.pkl")


fout = open(outDir+args.stem+"_"+chromosome+"_FarJunctions_duplicates.fa", mode = "wb")

#****SLICE FILE

with open(args.infile, mode ="rU") as f1:
    while True:
        next_n_lines = list(islice(f1,50000))
        print "next " + str(len(next_n_lines)) + "lines"
        if not next_n_lines: break
#        print next_n_lines
#        sys.stdout.flush()
    
        UI_counter = 0 #UI stands for unique identifier
        AllChrA = [0] * 10000000 # 10 million possible read differences
        AllChrB = [0] * 10000000
        LocationTuples_AllChr = [0]*26
        ExonsListA = [0]*10000000 # 5 mill possible tuples, add batch of additional 1 mill slots if necessary
        ExonsListB = [0]*10000000
        SeqListA = [0]*10000000
        SeqListB = [0]*10000000
        NameListA = [0]*10000000
        NameListB = [0]*10000000
        AllJunctions = [0]
        runtimestart = time.clock()

#        print "reading in lines"
        for line_raw in next_n_lines:
            
            line = line_raw.strip().split("\t")
 #           print line
 #           sys.stdout.flush()
            AllChrA[UI_counter] = line[0]
            AllChrB[UI_counter] = line[1]
            RangeA = locationparameters(line[0])
            RangeB = locationparameters(line[1])
            SortedLocation(str(UI_counter)+"A", RangeA.chr, RangeA.start, RangeA.stop)
            SortedLocation(str(UI_counter)+"B", RangeB.chr, RangeB.start, RangeB.stop)
            UI_counter+=1
        #    if UI_counter==20: break
         
        for i in range(0,25):
            if i<22: chrnum = "chr"+str(i+1)
            if i==22: chrnum ="chrX"
            if i==23: chrnum = "chrY"
            if i==24: chrnum = "chrM"
            
            exonfile = exonDir + "/exonsByStrand_" + chrnum + ".pkl"
            recordfile = recordDir + "/rec_" +chrnum + ".pkl"
            
         
            if patt_exonfile.search(exonfile): # only parse if this is an exon pickled file
                Picklefile = open(exonfile, 'rb')       
                PickleRecfile = open(recordfile, 'rb')

#                print "getting Exon pickle " + str(i)                
                ExonsOfChrPos, ExonsOfChrNeg = UnpickleExons(LocationTuples_AllChr[i], Picklefile, chrnum)
                
                
                if len(ExonsOfChrPos)>0 or len(ExonsOfChrNeg)>0:
#                    print "getting sequence from pickle " + str(i)
                    UnpickleSequence(chrnum, PickleRecfile, ExonsOfChrPos, ExonsOfChrNeg)
  
#        print "outputting junctions"          
        for i in range(0,UI_counter):
            RangeA = locationparameters(AllChrA[i])
            RangeB = locationparameters(AllChrB[i])
            APosStart = ChrAdjuster(RangeA.chr)
            APosStop= ChrAdjuster(RangeA.chr)+UI_counter
            ANegStart= ChrAdjuster(RangeA.chr)+StrandAdjuster(-1)
            ANegStop =ChrAdjuster(RangeA.chr)+UI_counter+StrandAdjuster(-1)
            BPosStart = ChrAdjuster(RangeB.chr)
            BPosStop= ChrAdjuster(RangeB.chr)+UI_counter
            BNegStart= ChrAdjuster(RangeB.chr)+StrandAdjuster(-1)
            BNegStop =ChrAdjuster(RangeB.chr)+UI_counter+StrandAdjuster(-1)
            MakePairs(fout, AllJunctions, i, RangeA.chr, RangeB.chr, ExonsListA[APosStart:APosStop][i], ExonsListA[ANegStart:ANegStop][i], ExonsListB[BPosStart:BPosStop][i], ExonsListB[BNegStart:BNegStop][i], SeqListA[APosStart:APosStop][i], SeqListA[ANegStart:ANegStop][i], SeqListB[BPosStart:BPosStop][i], SeqListB[BNegStart:BNegStop][i], NameListA[APosStart:APosStop][i], NameListA[ANegStart:ANegStop][i], NameListB[BPosStart:BPosStop][i], NameListB[BNegStart:BNegStop][i])
#            print NameListA[APosStart:APosStop][i]
#            print NameListA[ANegStart:ANegStop][i]
#            print NameListB[BPosStart:BPosStop][i]
#            print NameListB[BNegStart:BNegStop][i]
#            
#            print ExonsListA[APosStart:APosStop][i]
#            print ExonsListA[ANegStart:ANegStop][i]
#            print ExonsListB[BPosStart:BPosStop][i]
#            print ExonsListB[BNegStart:BNegStop][i]
            
             
        fout.flush()
        runtimestop = time.clock()   
        print "Unpickler run time: " + str(runtimestop-runtimestart)
        
        del UI_counter #UI stands for unique identifier
        del AllChrA # 10 million possible read differences
        del AllChrB
        del LocationTuples_AllChr
        del ExonsListA  # 5 mill possible tuples, add batch of additional 1 mill slots if necessary
        del ExonsListB 
        del SeqListA
        del SeqListB 
        del NameListA 
        del NameListB 
        del AllJunctions        
        
fout.close()
f1.close()

## ===========REMOVE DUPLICATES====================

f1 = open(outDir+args.stem+"_"+chromosome+"_FarJunctions_duplicates.fa", mode = "rb")
fout = open(outDir+args.stem+"_"+chromosome+"FarJunctions.fa", mode = "wb")


junctiondict={}

counter=0
            
for line in f1:
    start= time.clock()
    if line not in junctiondict and line[0]==">":
        junctiondict[line]=1
        fout.write(line+f1.next())
        counter+=1
        if counter%50000 ==0:
            fout.flush()
#            print "flushing, " + str(time.clock()-start)
f1.close()
fout.close()
        