#Identified junction class
#TODO give better description of members and methods
#TODO this class has too many members. They should be moved to the underlying SAMEntries

#Imports
import sys
import subprocess
import time
from SAMEntry import SAMEntry

class Junction(object):
    __slots__ = ["consensus","upstream_seq","downstream_seq","score","upstream_chromosome","downstream_chromosome",
                "upstream_start","downstream_start","upstream_stop","downstream_stop",
                "splice_site","bin_pair","bin_pair_group","linear","upstream_strand",
                "downstream_strand","span","upstream_genes","downstream_genes","splice_type","splice_gap",
                "constants_dict","took_reverse_compliment","splice_ind","constitutive_jcts","fusion"]
    def __init__(self,consensus,score,bin_pair_group,took_reverse_compliment,constants_dict):
        #Read in arguments
        self.consensus = consensus
        self.score = score
        self.bin_pair_group = bin_pair_group
        self.took_reverse_compliment = took_reverse_compliment
        self.constants_dict = constants_dict

        #Find chromosome, bin_pair and strand info from the first mapped read
        rep_bin_pair = self.bin_pair_group[0]
        self.upstream_chromosome = rep_bin_pair.five_prime_SAM.chromosome
        self.downstream_chromosome = rep_bin_pair.three_prime_SAM.chromosome
        self.bin_pair = rep_bin_pair.bin_pair
        sys.stdout.flush()
        five_prime_bin,three_prime_bin,strand_info = self.bin_pair.split("_")
        five_prime_chr = five_prime_bin.split(":")[0]
        five_prime_bin = five_prime_bin.split(":")[1]
        three_prime_chr = three_prime_bin.split(":")[0]
        three_prime_bin = three_prime_bin.split(":")[1]
        linear = True if int(five_prime_bin) <= int(three_prime_bin) else False
        self.linear = not linear if self.took_reverse_compliment else linear
        self.upstream_strand = rep_bin_pair.five_prime_strand
        self.downstream_strand = rep_bin_pair.three_prime_strand

        #Find splice ind, splice site, and upstream and downstream locations
        self.splice_ind = 0
        self.upstream_seq = ""
        self.downstream_seq = ""
        self.upstream_genes = "no_genes"
        self.downstream_genes = "no_genes"
        self.upstream_seq = ""
        self.downstream_seq = ""
        self.upstream_start = -1
        self.downstream_start = -1
        self.upstream_stop = -1
        self.downstream_stop = -1
        self.splice_site = -1
        self.splice_ind = -1
        self.span = 0
        self.splice_type = "None"
        self.splice_gap = 0
        self.fusion = False

        #Add self to the list of junctions that constitute this one
        self.constitutive_jcts = [self]

    def write_preliminary(self,fasta_name):
        with open(fasta_name,"w") as fasta_file:
            for junction in junctions:
                fasta_file.write(">jct_bin_pair:"+junction.bin_pair+"\n")
                fasta_file.write(junction.consensus+"\n")

    #TODO I think that this should rely more on the consensus
    def find_start_stop(self):
        #Use the first bin_pair_group to get guesses for positions
        rep_bin_pair = self.bin_pair_group[0]
        self.upstream_start = rep_bin_pair.five_prime_SAM.upstream_pos
        self.upstream_stop = rep_bin_pair.five_prime_SAM.downstream_pos
        self.downstream_start = rep_bin_pair.three_prime_SAM.upstream_pos
        self.downstream_stop = rep_bin_pair.three_prime_SAM.downstream_pos

        #Loop through the other bin_pair_groups to find widest range
        for bin_pair in self.bin_pair_group[1:]:
            five_prime_upstream_pos = bin_pair.five_prime_SAM.upstream_pos
            five_prime_downstream_pos = bin_pair.five_prime_SAM.downstream_pos
            three_prime_upstream_pos = bin_pair.three_prime_SAM.upstream_pos
            three_prime_downstream_pos = bin_pair.three_prime_SAM.downstream_pos
            if five_prime_upstream_pos  < self.upstream_start: self.upstream_start = five_prime_upstream_pos
            if five_prime_downstream_pos  > self.upstream_stop:  self.upstream_stop  = five_prime_downstream_pos

            if three_prime_upstream_pos < self.downstream_start: self.downstream_start = three_prime_upstream_pos
            if three_prime_downstream_pos > self.downstream_stop:  self.downstream_stop  = three_prime_downstream_pos

        self.upstream_stop = self.upstream_start+self.splice_ind
        self.downstream_start = self.downstream_stop-(len(self.consensus)-self.splice_ind)
        self.span = self.downstream_start-self.upstream_stop
        self.splice_site = self.upstream_start+self.splice_ind

    #Just adds the junction to the list of jcts that compose this one
    #Waits until a call to "combine junctions" to process adding them together
    def add_constitutive_junction(self,other):
        self.constitutive_jcts.append(other)

    #Want to combine all the bin pairs and see how good the "consensus" consensus looks
    #The "consensus" consensus will define the upstream and downstream stop/starts
    #First find the min and max splice inds to know how much left and right padding to add
    def combine_constitutive_junctions(self):
        if len(self.constitutive_jcts) == 1: #Only composed of itself
            return True
        print "Num constitutive jcts",len(self.constitutive_jcts)
        group_consensuses = []
        group_bin_pairs = []
        splice_inds = [jct.splice_ind for jct in self.constitutive_jcts]
        min_splice_ind = min(splice_inds)
        max_splice_ind = max(splice_inds)
        for constitutive_jct in self.constitutive_jcts:
            consensus = constitutive_jct.consensus
            splice_ind = constitutive_jct.splice_ind
            group_bin_pairs += constitutive_jct.bin_pair_group
            left_padding = " "*(max_splice_ind-splice_ind)
            print "linear: ",constitutive_jct.linear," splice type: ",constitutive_jct.splice_type
            print "splice ind:",splice_ind,"left padding:",len(left_padding)
            boxed_consensus = consensus[:splice_ind-1]+"|"+consensus[splice_ind]+"|"+consensus[splice_ind:]
            left_padded_consensus = left_padding+boxed_consensus
            group_consensuses.append(left_padded_consensus)

        max_len_consensus = max([len(consensus) for consensus in group_consensuses])
        group_consensuses = [consensus+" "*(max_len_consensus-len(consensus)) for consensus in group_consensuses]

        for consensus in group_consensuses:
            print len(consensus),consensus

        consensus_of_consensus = ""
        for base in range(max_len_consensus):
            prev_base = " "
            found_mismatch = False
            for consensus in group_consensuses:
                if prev_base == " " and consensus[base] != prev_base:
                    prev_base = consensus[base]
                elif consensus[base] != " " and consensus[base] != prev_base:
                    found_mismatch = True
                    break
            if found_mismatch:
                consensus_of_consensus += "N"
            elif prev_base != " ":
                consensus_of_consensus += prev_base
        
        print "====================="
        print consensus_of_consensus
        sys.stdout.flush()

    #Check to see if this jct represents a fusion
    def check_fusion(self):
#        if self.upstream_genes != "no_genes" and self.downstream_genes != "no_genes":
        upstream_gene_list = self.upstream_genes.split(",")
        downstream_gene_list = self.downstream_genes.split(",")
        found_shared_gene = False
        for upstream_gene in upstream_gene_list:
            if upstream_gene in downstream_gene_list:
                found_shared_gene = True
                break
        
        self.fusion = not found_shared_gene



#        else:
#            self.fusion = False


    #Format the junction to print in fasta form
    def fasta_string(self):
        fasta_str = ""
        fasta_str += ">|chromosome1:"+self.upstream_chromosome+"|"
        fasta_str += "genes1:"+self.upstream_genes+"|"
        fasta_str += "start1:"+str(self.upstream_start)+"|"
        fasta_str += "stop1:"+str(self.upstream_stop)+"|"
        fasta_str += "strand1:"+self.upstream_strand+"|_"

        fasta_str += "|chromosome2:"+self.downstream_chromosome+"|"
        fasta_str += "genes2:"+self.downstream_genes+"|"
        fasta_str += "start2:"+str(self.downstream_start)+"|"
        fasta_str += "stop2:"+str(self.downstream_stop)+"|"
        fasta_str += "strand2:"+self.downstream_strand+"|"

        fasta_str += "splice:"+str(self.splice_ind)+"|"
        fasta_str += "span:"+str(self.span)+"|"
        fasta_str += "score:"+str(self.score)+"|"
        fasta_str += "fusion:"+str(self.fusion)+"|"
        fasta_str += "num:"+str(len(self.bin_pair_group))+"|"
        fasta_str += "splice-gap:"+str(self.splice_gap)+"|"
        fasta_str += "splice-type:"+self.splice_type+"|\n"

        # Add N padding to the consensus to get a uniform len
        splice_flank_len = int(self.constants_dict["splice_flank_len"])
        left_padding = "N"*(splice_flank_len-self.splice_ind)
        right_padding = "N"*(splice_flank_len-(len(self.consensus)-self.splice_ind))
        full_consensus = left_padding+self.consensus+right_padding

        #NOTE Helpful to verify that all the splice inds wind up at the same location
        #full_consensus = left_padding+self.consensus[:self.splice_ind]+"|"+self.consensus[self.splice_ind:]+right_padding

        fasta_str += full_consensus+"\n"
        return fasta_str

    #More human readable format
    def __str__(self):
        out_str = ""
        out_str += "Junction with bin pair ["+self.bin_pair+"] with ["+str(len(self.bin_pair_group))+"] reads mapped\n"
        out_str += "Linear " if self.linear else "Non-Linear "
        out_str += "Upstream on the "+self.upstream_strand+" strand and downstream on the "+self.downstream_strand+"\n"
        out_str += "5' map position ["+str(self.upstream_start)+"-"+str(self.upstream_stop)+"]\n"
        out_str += "3' map position ["+str(self.downstream_start)+"-"+str(self.downstream_stop)+"]\n"
        out_str += "This causes a span of ["+str(self.span)+"]\n"
        out_str += "Consensus with score ["+str(self.score)+"] and splice site ["+str(self.splice_site)+"]:\n"
        out_str += self.consensus+"\n"
        out_str += self.upstream_seq+"\n"
        out_str += " "*len(self.upstream_seq)+self.downstream_seq+"\n"
        out_str += "Upstream genes ["+self.upstream_genes+"]\n"
        out_str += "Downstream genes ["+self.downstream_genes+"]\n"
        return out_str
 
    def __lt__(self,other):
        return self.bin_pair < other.bin_pair


    def save(obj):
        return (obj.__class__, obj.__dict__)

    def load(cls, attributes):
        obj = cls.__new__(cls)
        obj.__dict__.update(attributes)
        return obj

    """
    def __getstate__(self):
        state = (self.consensus,self.upstream_seq,self.downstream_seq,self.score,
                self.upstream_chromosome,self.downstream_chromosome,self.upstream_start,self.downstream_start,
                self.upstream_stop,self.downstream_stop,self.splice_site,self.bin_pair,self.bin_pair_group,self.linear
                self.upstream_strand,self.downstream_strand,self.span,self.upstream_genes,self.downstream_genes,self.found_splice,
                self.constants_dict,self.took_reverse_compliment,self.splice_ind,self.constitutive_jcts,self.fusion)
        return state

    def __setstate__(self,state):
        obj = cls.__new__(cls)
        obj.__dict__.update(attributes)
        self.consensus,self.upstream_seq,self.downstream_seq,self.score,
        self.upstream_chromosome,self.downstream_chromosome,self.upstream_start,self.downstream_start,
        self.upstream_stop,self.downstream_stop,self.splice_site,self.bin_pair,self.bin_pair_group,self.linear
        self.upstream_strand,self.downstream_strand,self.span,self.upstream_genes,self.downstream_genes,self.found_splice,
        self.constants_dict,self.took_reverse_compliment,self.splice_ind,self.constitutive_jcts,self.fusion = state
    """
