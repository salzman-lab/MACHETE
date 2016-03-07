# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:38:44 2016

@author: Gillian
"""

IDFile=open("/scratch/PI/horence/gillian/CML_test/aligned/CML/circReads/ids/ENCFF000HOC1_1__output_withFJ.txt", mode="rU")
AlignmentFile1=open("/scratch/PI/horence/gillian/CML_test/FarJunc_Feb22/FarJunctionAlignments/ENCFF000HOC1/unaligned_ENCFF000HOC1_1.sam", mode ="rU")
AlignmentFile2=open("/scratch/PI/horence/gillian/CML_test/FarJunc_Feb22/FarJunctionAlignments/ENCFF000HOC1/unaligned_ENCFF000HOC1_2.sam", mode ="rU")

ID_Dict={}

for line in IDFile:
    if line[0]=="@":
        continue
    if "FJ" not in line.strip().split("\t")[3]:
        continue    
    ID_Dict[line.strip().split("\t")[0]]=0

IDFile.close()


for line in AlignmentFile1:
    if line[0]=="@":
        continue
    
    if line.strip().split("\t")[0] in ID_Dict:
        ID_Dict[line.strip().split("\t")[0]]+=1
AlignmentFile1.close()  

        
for line in AlignmentFile2:
    if line[0]=="@":
        continue
    if line.strip().split("\t")[0] in ID_Dict:
        ID_Dict[line.strip().split("\t")[0]]+=1
AlignmentFile2.close()    

for key in ID_Dict:
    if ID_Dict[key]==0:
        print key
    
        