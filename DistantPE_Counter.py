# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:36:05 2015

@author: Gillian
"""

# this class assigns two paired reads to "ChrA:X1-X2, ChrB:Y1-Y2" and forces the chrA read to be "smaller"
class Paired_Reads:
    def __init__ (self, S1, S2):
        PE1 = S1.replace(":"," ").replace("-"," ").replace("strand="," ").split(" ")
        PE2 = S2.replace(":"," ").replace("-"," ").replace("strand="," ").split(" ") 
        ChrA = PE1[0]
        X1 = int(PE1[1])
        X2 = int(PE1[2])
        StrandA = PE1[3]
        ChrB = PE2[0]
        Y1 = int(PE2[1])
        Y2 = int(PE2[2])
        StrandB = PE2[3]
        PE1first = True
        
        if StrandA=="":
            StrandA= "-"
        if StrandB=="":
            StrandB="-"
        
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
            self.strandA = StrandA
            self.y1 = Y1
            self.y2 = Y2
            self.ylen = Y2-Y1
            self.strandB = StrandB
            self.PE1 = str(ChrA+":"+str(X1)+"-"+str(X2))
            self.PE2 = str(ChrB+":"+str(Y1)+"-"+str(Y2))
            self.bothPE = str(ChrA+":"+str(X1)+"-"+str(X2)+"strand="+StrandA+"\t"+ChrB+":"+str(Y1)+"-"+str(Y2)+"strand="+StrandB)
        if PE1first == False:
            self.chrA = ChrB
            self.chrB = ChrA
            self.x1 = Y1
            self.x2 = Y2
            self.strandA = StrandB
            self.xlen = Y2-Y1
            self.y1 = X1
            self.y2 = X2
            self.strandb= StrandA
            self.ylen = X2-X1
            self.PE2 = str(ChrA+":"+str(X1)+"-"+str(X2))
            self.PE1 = str(ChrB+":"+str(Y1)+"-"+str(Y2))                
            self.bothPE = str(ChrB+":"+str(Y1)+"-"+str(Y2)+"strand="+StrandB+"\t"+ChrA+":"+str(X1)+"-"+str(X2)+"strand="+StrandA)    
            
import os
import glob
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--Dir", required = True, help = "path to working directory")
args=parser.parse_args()

PairLibrary={}  #documents frequency of that pair occurring

if args.Dir[-1] != "/":
    path = args.Dir + "/"
else:
    path = args.Dir
    
    
os.chdir(path)

for name in glob.glob(os.path.join(path,"*_distant_pairs.txt")):
    (path,filename)=os.path.split(name)
    f1= open(name, mode ="rU")  
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
    
fout = open("Distant_PE_frequency.txt", mode="w")

for key, value in sorted(PairLibrary.iteritems(), key=lambda (k,v): (v,k), reverse=True ):
    fout.write("%s\t %s" % (key, value) + "\n")

fout.close()


    
