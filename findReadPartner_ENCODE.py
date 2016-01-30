# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:35:23 2015

@author: Gillian
"""


import argparse
import sys

###  RUN MODE
parser=argparse.ArgumentParser()
parser.add_argument("-1","--R1file", required = True, help = "path to file with reads of interest")
parser.add_argument("-2", "--R2file", required = True, help = "path to file where you want to find buddy reads")
parser.add_argument("-o", "--outfile", required=True, help = "name of output file x.txt (need to type txt)")
args=parser.parse_args()

#TESTING MODE
#path = "/Users/Gillian/Desktop/sherlock/"

# change the input path to the path where your file exists

f1 = open(args.R1file, mode = "rU")
f2 = open(args.R2file, mode = "rU")

print args.R1file
print args.R2file

#create file list of paired files
fout = open(args.outfile, mode="w")  
fout.write("All read IDs from " + args.R1file + " and their partners in " + args.R2file+"\n")
sys.stdout.flush() #force print command above to output to screen (gives user an idea of progress)
    
for line_raw1 in f1:
 
    if line_raw1[0] == "@":       
        continue
    
    line1 = line_raw1.strip().split("\t")    
    try:
        line2=f2.next().strip().split("\t")
    except: break

    while line2[0][0] == "@":
        line2=f2.next().strip().split("\t")
        
#    print line1
#    print line2
#        # if readID in file1 is less than readID in file2, move to next line in file1
    if line1[0][:-3]<line2[0][:-3]:
#            print "line1<line2"    
        continue
    #if readID in file2 is less than readID in file1, move to next line in file2
    while line2[0][:-3]<line1[0][:-3]:
        try:
            line_raw2 = f2.next()
            line2 = line_raw2.strip().split("\t")
        except: break
#            print "line2<line1"
        
    # if readID in file1 and file2 are equivalent, then these reads can be compared 
    if line1[0][:-3]==line2[0][:-3]:
#             print "found matching Read IDs"
        output = line1[0][:-3] + "\t"+ line1[2]+ "\t" + line2[2]+":"+line2[3]+"\n"
#            print output 
        fout.write(output)
   
            # if at the end of file2, then break
            # if at the end of file1, for loop will automatically discontinue
        try: 
            line_raw2 = f2.next()
            line2 = line_raw2.strip().split("\t")
        except: break    
        
fout.close()            
f1.close()
f2.close()
       
