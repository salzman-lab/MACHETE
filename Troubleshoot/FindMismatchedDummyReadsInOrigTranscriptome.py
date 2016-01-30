# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:35:23 2015

@author: Gillian
"""

import argparse

###  RUN MODE
parser=argparse.ArgumentParser()
parser.add_argument("-m","--Mfile", required = True, help = "mismatch reads file")
parser.add_argument("-f1","--file1", required=True, help = "file1")
parser.add_argument("-f2","--file2", required=True, help = "file2")
parser.add_argument("-n", "--name", required=True, help = "name")

args=parser.parse_args()
mismatchfile = args.Mfile
file1 = args.file1
file2 = args.file2


##TESTING MODE
#mismatchfile = "/Users/Gillian/Desktop/transcriptometest/ENCFF000HOC1_mismatched_reads.txt"
#file1 = "blah"
#file2 = "blah"
# change the input path to the path where your file exists

# for each pair of read1 and read2 files -- 

m1 = open(mismatchfile, mode="rU")
f1 = open(file1, mode = "rU")
f2 = open(file2, mode = "rU")
fout = open(args.name+"_MismatchedDummyReadsinOrigTranscriptome.txt", mode="w")

m1.next()
m1.next()


for line_rawM in m1:
    lineM = m1.next().strip().split("\t")
    line1 = f1.next().strip().split("\t")
    line2 = f2.next().strip().split("\t")
     
    if lineM[0] < line1[0][:-3] or lineM[0]< line2[0][:-3]:
        continue

    while line1[0][:-3]<lineM[0]:
        try:
            line1 = f1.next().strip().split("\t")
        except: break
    
    while line2[0][:-3]<lineM[0]:
        try:
            line2 = f2.next().strip().split("\t")
        except: break

    if lineM[0] == line1[0][:-3] and lineM[0] == line2[0][:-3]:        
        output = lineM[0] + "\t"+ line1[2]+ "\t" + line2[2]+"\n"
#            print output 
        
        if line1[4]>=10 and line2[4]>=10:
            fout.write(output)
   
            # if at the end of file2, then break
            # if at the end of file1, for loop will automatically discontinue
        try: 
            line1 = f1.next().strip().split("\t")
            line2 = f2.next().strip().split("\t")
        except: break    
            
fout.close()            
f1.close()
f2.close()
   