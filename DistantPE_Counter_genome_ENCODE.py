# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:36:05 2015

@author: Gillian
"""

# this class assigns two paired reads to "ChrA:X1-X2, ChrB:Y1-Y2" and forces the chrA read to be "smaller"
class Paired_Reads:
    def __init__ (self, S1, S2):
        PE1 = S1.replace(":"," ").replace("-"," ").split(" ")
        PE2 = S2.replace(":"," ").replace("-"," ").split(" ") 
        ChrA = PE1[0]
        X1 = int(PE1[1])
        X2 = int(PE1[2])
        ChrB = PE2[0]
        Y1 = int(PE2[1])
        Y2 = int(PE2[2])
        PE1first = True

        
        if ChrB<ChrA:
            PE1first = False
        if ChrA==ChrB and Y1 < X1:
            PE1first = False
        try:
            if int(ChrA) <= int(ChrB):
                PE1first = True
            elif int(ChrB) < int(ChrA):
                PE1first = False
        except: 
            pass
        if PE1first == True:
            self.chrA = ChrA
            self.chrB = ChrB
            self.x1 = X1
            self.x2 = X2
            self.xlen = X2-X1
            self.y1 = Y1
            self.y2 = Y2
            self.ylen = Y2-Y1
            self.PE1 = str(ChrA+":"+str(X1)+"-"+str(X2))
            self.PE2 = str(ChrB+":"+str(Y1)+"-"+str(Y2))
            self.bothPE = str(ChrA+":"+str(X1)+"-"+str(X2)+"\t"+ChrB+":"+str(Y1)+"-"+str(Y2))
        if PE1first == False:
            self.chrA = ChrB
            self.chrB = ChrA
            self.x1 = Y1
            self.x2 = Y2
            self.xlen = Y2-Y1
            self.y1 = X1
            self.y2 = X2
            self.ylen = X2-X1
            self.PE2 = str(ChrA+":"+str(X1)+"-"+str(X2))
            self.PE1 = str(ChrB+":"+str(Y1)+"-"+str(Y2))                
            self.bothPE = str(ChrB+":"+str(Y1)+"-"+str(Y2)+"\t"+ChrA+":"+str(X1)+"-"+str(X2))
            

def writetofile(outfile,key,value):
    outfile.write("%s\t %s" % (key, value) + "\n")

import os
import glob
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--Dir", required = True, help = "path to Far Junc parent directory")
parser.add_argument("-s", "--stem", required=True, help="unique Stem that identifies the file to count distant PE from")
args=parser.parse_args()

PairLibrary={}  #documents frequency of that pair occurring

if args.Dir[-1] != "/":
    path = args.Dir + "/DistantPEFiles/"
else:
    path = args.Dir + "DistantPEFiles/"

outpath= path + args.stem+ "/"
    
os.chdir(path)

for name in glob.glob(os.path.join(path + "*" + args.stem  + "*_distant_pairs.txt")):
    (path,filename)=os.path.split(name)
    f1= open(name, mode ="rU")
    print name
    f1.next()
    f1.next()
    f1.next()
    
    for line_raw in f1:
        line = line_raw.strip().split("\t")
        read = Paired_Reads(line[1], line[2])
        
        if not read.bothPE in PairLibrary:
            PairLibrary[read.bothPE]=1
        else:
            PairLibrary[read.bothPE]+=1
        
    f1.close()

outfile1 = open(outpath+"chr1_Distant_PE_frequency.txt", mode="w")
outfile2 = open(outpath+"chr2_Distant_PE_frequency.txt", mode="w")
outfile3 = open(outpath+"chr3_Distant_PE_frequency.txt", mode="w")
outfile4 = open(outpath+"chr4_Distant_PE_frequency.txt", mode="w")
outfile5 = open(outpath+"chr5_Distant_PE_frequency.txt", mode="w")
outfile6 = open(outpath+"chr6_Distant_PE_frequency.txt", mode="w")
outfile7 = open(outpath+"chr7_Distant_PE_frequency.txt", mode="w")
outfile8 = open(outpath+"chr8_Distant_PE_frequency.txt", mode="w")
outfile9 = open(outpath+"chr9_Distant_PE_frequency.txt", mode="w")
outfile10 = open(outpath+"chr10_Distant_PE_frequency.txt", mode="w")
outfile11 = open(outpath+"chr11_Distant_PE_frequency.txt", mode="w")
outfile12= open(outpath+"chr12_Distant_PE_frequency.txt", mode="w")
outfile13 = open(outpath+"chr13_Distant_PE_frequency.txt", mode="w")
outfile14 = open(outpath+"chr14_Distant_PE_frequency.txt", mode="w")
outfile15 = open(outpath+"chr15_Distant_PE_frequency.txt", mode="w")
outfile16 = open(outpath+"chr16_Distant_PE_frequency.txt", mode="w")
outfile17 = open(outpath+"chr17_Distant_PE_frequency.txt", mode="w")
outfile18 = open(outpath+"chr18_Distant_PE_frequency.txt", mode="w")
outfile19 = open(outpath+"chr19_Distant_PE_frequency.txt", mode="w")
outfile20 = open(outpath+"chr20_Distant_PE_frequency.txt", mode="w")
outfile21 = open(outpath+"chr21_Distant_PE_frequency.txt", mode="w")
outfile22= open(outpath+"chr22_Distant_PE_frequency.txt", mode="w")
outfileX = open(outpath+"chrX_Distant_PE_frequency.txt", mode="w")
outfileY = open(outpath+"chrY_Distant_PE_frequency.txt", mode="w")

FileDict={"1":outfile1, "2":outfile2,"3":outfile3,"4":outfile4, "5":outfile5, "6":outfile6,"7":outfile7,"8":outfile8,"9":outfile9,"10":outfile10, "11":outfile11,"12":outfile12,"13":outfile13, "14":outfile14,"15":outfile15, "16":outfile16,"17":outfile17,"18":outfile18,"19":outfile19,"20":outfile20, "21":outfile21,"22":outfile22,"X":outfileX, "Y":outfileY }

for key, value in sorted(PairLibrary.iteritems(), key=lambda (k,v): (v,k), reverse=True ):
    chromosome= key.split(":")[0]
    writetofile(FileDict[chromosome], key, value)


outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()
outfile5.close()
outfile6.close()
outfile7.close()
outfile8.close()
outfile9.close()
outfile10.close()
outfile11.close()
outfile12.close()
outfile13.close()
outfile14.close()
outfile15.close()
outfile16.close()
outfile17.close()
outfile18.close()
outfile19.close()
outfile20.close()
outfile21.close()
outfile22.close()
outfileX.close()
outfileY.close()

    
