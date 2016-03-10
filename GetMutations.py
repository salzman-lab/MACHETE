# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 18:26:59 2016

@author: Gillian
"""

class MAF:
    def __init__(self, line_raw, RefGenome):
        line=line_raw.strip().split("\t")
        
        self.build=line[3]
        self.chr=line[4]
        self.loc=line[5]
        self.type=line[9] # options: SNP, DNP, TNP, ONP, INS, DEL
        self.RefSeq=line[10].strip()
        self.Seq1=line[11].strip()
        self.Seq2=line[12].strip()
        self.validation=line[24] # options "Valid" or "Untested" -- should we use untested?
        self.verification=line[23] # options "Verified" or "Unknown"     
#A "Verified" Verification_Status confirms a mutation as germline or somatic using the same technology as the one that discovered the mutation.
#A "Valid" Validation_Status confirms a mutation as germline or somatic using a different technology than the one that discovered the mutation.


        #ForgiveSeq1 = the first mutation to forgive
        #ForgiveSeq2 = either integer 0 or the second mutation to forgive

        self.ForgiveSeq2=0
        if self.Seq1!=self.RefSeq:
            self.ForgiveSeq1=self.Seq1
            if self.Seq2!=self.RefSeq and self.Seq2!=self.Seq1:
                self.ForgiveSeq2=self.Seq2
        if self.Seq1==self.RefSeq and self.Seq2!=self.RefSeq:
            self.ForgiveSeq1=self.Seq2
            
        #ForgiveSeq2 doesn't exist if there is only one sequence to forgive at that location
        
        ## REQUIREMENTS TO CONTINUE
        
        if RefGenome=="HG19" and self.build!="37":
            self.goodMAF=False
        elif RefGenome=="HG38" and self.build!="38":
            self.goodMAF=False
        elif self.type=="INS" or self.type=="DEL":
            self.goodMAF=False
        elif self.ForgiveSeq1==self.RefSeq and self.ForgiveSeq2==self.RefSeq:
            self.goodMAF=False

        #### could force script to only allow validated, but currently we will use all
        #### mutations including untested/unknown

#        elif self.validation=="Untested" and self.verification=="Unknown":
#            self.goodMAF=False
        else:
            self.goodMAF=True
        
      
########################################################
### END FUNCTIONS, BEGIN WORK ##########################        
########################################################
        

import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-v", "--VariantFile", required = True, help = "Variant calling .maf file")
#parser.add_argument("-o", "--origDir", required = True, help = "orig dir containing genome/reg/junc/transcriptome alignments")
parser.add_argument("-f", "--FJDir", required=True, help = "Far Junction base directory")
#parser.add_argument("-s", "--stem", required=True, help="unique Stem that identifies the file to count distant PE from")
parser.add_argument("-b", "--build", required=True, help = "genome build that you are using, HG19 or HG38")
args=parser.parse_args()

#if args.origDir[-1]!="/":
#    args.origDir += "/"
if args.FJDir[-1]!="/":
    args.FJDir += "/"


VariantDict={} 
# key of variant dict: (chr, absolute location) tuple
# value of variant dict: (ref base(s), forgiven seq 1, forgiven seq 2 ) tuple
# only accounts for SNP, DNP, TNP, ONP
# does not account for INS or DEL
# MAF file only allows "+" strand notation so I will only
# keep track of "lowest" location 
# chr5, 5000, and is a SNP, then I will note loc (5,5000)
# if DNP, then loc (5, 5000) means chr5: 5000-5001
# if TNP then loc (5,5000) means chr 5: 5000-5002
# if ONP then loc (5, 5000) means chr5: 5000-5000+length of ONP

## read in Mutation Annotation Format file
MAFfile=open(args.VariantFile, mode= "rU")

for line_raw in MAFfile:
    if line_raw[0:11]=="Hugo_Symbol":
        continue   
    MAFinfo=MAF(line_raw, args.build)
    if MAFinfo.goodMAF==True:
        VariantDict[(MAFinfo.chr, MAFinfo.loc)]=(MAFinfo.RefSeq, MAFinfo.ForgiveSeq1, MAFinfo.ForgiveSeq2)
MAFfile.close()

## write out a quick text file with the variant dictionary per julia's request
# in the future, if they have their own variant input file then I can read it in and skip
# reading in the MAF file.

MutationLog=open(args.FJDir+"ForgiveMutations/MutationList.txt", mode="w")

MutationLog.write("chr\tloc\tRefSeq\tMutation1\tMutation2\t BUILD:"+ args.build +"\n")

for key in VariantDict:
    (chromosome,location)=key
    (RefSeq, Mutation1, Mutation2) = VariantDict[key] 

    MutationLog.write("chr"+chromosome+"\t"+location+"\t"+RefSeq+"\t"+Mutation1)
    if Mutation2==0:
        MutationLog.write("\n")
    else:
        print "Mutation1!=Mutation2"
        print Mutation2
        MutationLog.write(Mutation2+"\n")

MutationLog.close()


# go through genome, reg, junc, transcriptome, FJ
# if same mutation, record read ID, # mismatches to forgive, which file from which it is forgiven.

