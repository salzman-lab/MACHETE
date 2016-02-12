# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 15:17:27 2015

@author: Gillian
"""


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", required=True, help = "input file")
parser.add_argument("-o", "--outDir", required=True, help = "output directory")
parser.add_argument("-s", "--stem", required = True, help = "stem name" )
parser.add_argument("-n", "--MaxInDel", required=True, help = "# of indels on each side to test")

args=parser.parse_args()

if args.outDir[-1]!="/":
    args.outDir+="/"
    

# indels are inserted into the FarJunction.fa file
# insertions are SEQUENCEANNNNNSEQUENCEB where 2N is the number of indels given in the argument parser
# deletions are SEQA-AAAAA[AAABBB]BBBSEQB where N*AAA and N*BBB are removed from the junction interface.
# x,..,5, 4, 3, 2, 1 deletions on each side, then 2, 4, 6, 8, 10,.., 2X N's inserted into junction


counter=0
f1 = open(args.infile, mode = "rU")


for i in range(1,int(args.MaxInDel)+1):  
    f1.seek(0)

    fout = open(args.outDir + args.stem + "_FJ_Indels_" +str(i)+".fa", mode ="w")
    print "writing indels"+str(i)+".fa"
    
    for line_raw in f1:
        
        counter+=1
        

        if line_raw[0] == ">" :
            JunctionName=line_raw.strip()
            JunctionSeq=""
        
        else:
            JunctionSeq += line_raw.strip()
            
        if len(JunctionSeq)>290:
            LeftExon = JunctionSeq[0:150]
            RightExon = JunctionSeq[150:300]
            #        print JunctionName
            #print JunctionSeq
        
            fout.write(JunctionName + "|DEL"+str(i)+"\n")
            fout.write(LeftExon[0:-i]+RightExon[i:]+"\n")
            InsertN = "N"*(i*2)
            fout.write(JunctionName + "|INS"+str(i)+"\n")
            fout.write(LeftExon+InsertN+RightExon+"\n")
        
        if counter==5000:
            fout.flush()
    fout.close()

f1.close()
