# General Imports
import subprocess
import time
import sys
import re

# Specific Imports
sys.path.append('./denovo_class_utils/')
from denovo_consensus_utils import *
from Junction import Junction
from BinPair import BinPair
from GTFEntry import GTFEntry
from SAMEntry import SAMEntry
from FastQEntry import FastQEntry

#######################################
#   Get Reference and GTF from Mode   #
#######################################
# Finding the correct reference index based on the mode
# Human defaults are used and currently the only supported mode
# Could imagine lots of if else statements with other supported references though
# Putting this in the utils file should allow for easy reference addition
def get_reference_and_gtf_from_mode(mode):
    index_path = "/scratch/PI/horence/rob/index/"
    reference = index_path+"hg19"
    gtf = "/scratch/PI/horence/rob/index/hg19_genes.gtf"

    return reference,gtf


#################################
#   Find Bin Pair Group Ranges  #
#################################
# Run through the bin pair groupings to find the end index of each group
# Example bin pair list:
# [0] chr10:45_chr10:55
# [2] chr10:45_chr10:55
# [3] chr10:98_chr10:99
# [4] chr10:145_chr10:155
# [5] chr10:145_chr10:155
# [6] chr10:145_chr10:155
# 
# Output: [2,3,6]
def find_bin_pair_group_ranges(bin_pairs,constants_dict):
    bin_pair_group_ends = []
    prev_bin_pair = ""
    for bin_pair_ind in range(len(bin_pairs)):
        curr_bin_pair = bin_pairs[bin_pair_ind].bin_pair
        if curr_bin_pair != prev_bin_pair:
            bin_pair_group_ends.append(bin_pair_ind)
            prev_bin_pair = curr_bin_pair

    # Remove the 0 at the front of the bin_pair list
    bin_pair_group_ends = bin_pair_group_ends[1:]

    # Add on the last bin_pair end which will necessarily be the end of the bin_pair list
    bin_pair_group_ends.append(len(bin_pairs))

    # Find the bin pair ranges
    # TODO this whole function can definitely be streamlined
    bin_pair_group_ranges = []
    for bin_pair_ind in range(len(bin_pair_group_ends)-1):
        start_ind = bin_pair_group_ends[bin_pair_ind]
        stop_ind = bin_pair_group_ends[bin_pair_ind+1]
        if stop_ind-start_ind > constants_dict["group_member_cutoff"]:
            bin_pair_group_ranges.append([start_ind,stop_ind])
    print bin_pair_group_ends
    print bin_pair_group_ranges
    sys.stdout.flush()
    return bin_pair_group_ranges


################################
#   Build Junction Sequences   #
################################
# Run through the bin_pair_group_ends to perform:
# (1) checking if there are enough in each group
# (2) padding the sequences in the first bin
# (3) creating a consensus in the first bin
# (4) scoring the consensus in the first bin
# (5) repeat (2)-(4) for the second bin
# (6) average the two consensus scores and see if they are below a cutoff
def build_junction_sequences(bin_pairs,bin_pair_group_ranges,full_path_name,constants_dict):
    group_member_cutoff = constants_dict["group_member_cutoff"]
    consensus_score_cutoff = constants_dict["consensus_score_cutoff"]
    bin_size = constants_dict["bin_size"]
    reference = constants_dict["reference"]
    denovo_junctions = []

    # Look back at the original full path to get seq lines
    unaligned_file = open(full_path_name,"r")
    unaligned_reads = unaligned_file.readlines() 
    unaligned_file.close()
    unaligned_reads = [unaligned_reads[ind] for ind in range(len(unaligned_reads)) if ind%4 == 0 or ind%4 == 1]
    id_to_seq = {}

    # Build the dictionary of read_id to the full read
    for ind in range(0,len(unaligned_reads),2):
        key = unaligned_reads[ind].replace("\n","")
        value = unaligned_reads[ind+1].replace("\n","")
        id_to_seq[key] = value


    for bin_pair_group_range in bin_pair_group_ranges:
        junction_num = "("+str(bin_pair_group_ranges.index(bin_pair_group_range)+1)+"/"+str(len(bin_pair_group_ranges))+")"
        #print junction_num
        sys.stdout.flush()

        start_build_junction = time.time()
        start_ind = bin_pair_group_range[0]
        stop_ind = bin_pair_group_range[1]
        group_members = bin_pairs[start_ind:stop_ind]

        # If there are not enough group members then skip this group
        if len(group_members) < group_member_cutoff: continue

        # Otherwise start thinking about getting strandedness
        five_prime_strand = group_members[0].five_prime_strand
        three_prime_strand = group_members[0].three_prime_strand

        # Find the consensus sequence and score
        # Takes just the 5' ends to get the pos
        # the full original unaligned seq
        mapped_reads = [member.five_prime_SAM for member in group_members]
        start_build_consensus = time.time()
        bin_consensus,bin_score = build_and_score_consensus(mapped_reads,five_prime_strand,id_to_seq,bin_size)
        write_time("--Time to build a single consensus n="+str(len(mapped_reads))+": "+junction_num,start_build_consensus,constants_dict["timer_file_path"])

        # TODO what do I do if the five and three prime strands are not the same?
        # If the reverse compliment was taken above then take the rev compliment of the consensus too
        took_reverse_compliment = False
        if five_prime_strand == "-" and three_prime_strand == "-":
            group_members = [member.take_reverse_compliment() for member in group_members]
            bin_consensus = reverse_compliment(bin_consensus)
            took_reverse_compliment = True

        # If the bin score is high enough then add it
        if bin_score < consensus_score_cutoff:
            denovo_junction = Junction(bin_consensus,bin_score,group_members,took_reverse_compliment,constants_dict)
            denovo_junctions.append(denovo_junction)

        write_time("-Time to completely process a single junction: "+junction_num,start_build_junction,constants_dict["timer_file_path"])
    return denovo_junctions


########################
#   Find Splice Inds   #
########################
# Runs bowtie on all of the possible splice sites of all possible junctions
# Returns a dict keyed by jct_id and valued by a list of cut sites
def find_splice_inds(denovo_junctions,constants_dict):
    # Gather info from the constants dictionary
    splice_finder_temp_name = constants_dict["out_dir"]+"splice_finder_temp_"
    min_score = constants_dict["splice_finding_min_score"]
    max_mismatches = int(constants_dict["splice_finding_allowed_mismatches"])
    read_gap_score = constants_dict["read_gap_score"]
    ref_gap_score = constants_dict["ref_gap_score"]
    allowed_mappings = constants_dict["splice_finding_allowed_mappings"]
    num_threads = constants_dict["num_threads"]
    reference = constants_dict["reference"]

    five_prime_fa_file = splice_finder_temp_name+"5_prime.fa"
    three_prime_fa_file = splice_finder_temp_name+"3_prime.fa"
    five_temp_file = open(five_prime_fa_file,"w")
    three_temp_file = open(three_prime_fa_file,"w")

    # Write out all the possible splice sites for every jct out to a 5' and 3' file
    for jct_ind in range(len(denovo_junctions)):
        junction = denovo_junctions[jct_ind]
        first_third = len(junction.consensus)/3
        second_third = first_third*2
        split_boundaries = range(first_third,second_third)
        five_prime_list = [junction.consensus[:split_boundary] for split_boundary in split_boundaries]
        five_prime_fa_list = [">jct_"+str(jct_ind)+"_ind_"+str(ind)+"\n"+five_prime_list[ind] for ind in range(len(five_prime_list))]
        three_prime_list = [junction.consensus[split_boundary:] for split_boundary in split_boundaries]
        three_prime_fa_list = [">jct_"+str(jct_ind)+"_ind_"+str(ind)+"\n"+three_prime_list[ind] for ind in range(len(three_prime_list))]

        five_temp_file.write("\n".join(five_prime_fa_list)+"\n")
        three_temp_file.write("\n".join(three_prime_fa_list)+"\n")

    # Don't forget to close the files
    five_temp_file.close()
    three_temp_file.close()
   
    # Map the temp files above to the reference and save in temp sam files
    # Need to specify the -f flag because the inputs are fasta files
    five_prime_mapped_name = splice_finder_temp_name+"5_prime.sam"
    three_prime_mapped_name = splice_finder_temp_name+"3_prime.sam"
    with open(five_prime_mapped_name,"w") as five_prime_mapped:
        subprocess.call(
            ["bowtie2", "-f", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
             "-x", reference, five_prime_fa_file], stdout=five_prime_mapped)
    with open(three_prime_mapped_name,"w") as three_prime_mapped:
        subprocess.call(
            ["bowtie2", "-f", "--no-sq", "--no-unal", min_score, read_gap_score, ref_gap_score, "-p", num_threads,
             "-x", reference, three_prime_fa_file], stdout=three_prime_mapped)

    # Walk through the 5' sam output and build up a dictionary by seq number
    five_prime_sam_file = open(five_prime_mapped_name,"r")
    five_sam_line = five_prime_sam_file.readline()
    while five_sam_line and "@" == five_sam_line[0]: #Read past the header lines
        five_sam_line = five_prime_sam_file.readline()
    
    three_prime_sam_file = open(three_prime_mapped_name,"r")
    three_sam_line = three_prime_sam_file.readline()
    while three_sam_line and "@" == three_sam_line[0]: #Read past the header lines
        three_sam_line = three_prime_sam_file.readline()

    best_splices = {}
    sam_five_list = []
    sam_three_list = []
    prev_jct_ind = 0
    while five_sam_line and three_sam_line:
        five_sam_entry = SAMEntry(five_sam_line)
        three_sam_entry = SAMEntry(three_sam_line)
        five_jct_ind = int(five_sam_entry.read_id.split("_")[1])
        three_jct_ind = int(three_sam_entry.read_id.split("_")[1])

        # Both the 5' and 3' sam are the prev_jct_ind
        if five_jct_ind == prev_jct_ind and three_jct_ind == prev_jct_ind:
            if five_sam_entry.num_mismatches <= max_mismatches:
                sam_five_list.append(five_sam_entry)
            if three_sam_entry.num_mismatches <= max_mismatches:
                sam_three_list.append(three_sam_entry)
            five_sam_line = five_prime_sam_file.readline()
            three_sam_line = three_prime_sam_file.readline()

        # Only the 5' sam is at the prev_jct_ind
        elif five_jct_ind == prev_jct_ind:
            if five_sam_entry.num_mismatches <= max_mismatches:
                sam_five_list.append(five_sam_entry)
            five_sam_line = five_prime_sam_file.readline()

        # Only the 3' sam is at the prev_jct_ind
        elif three_jct_ind == prev_jct_ind:
            if three_sam_entry.num_mismatches <= max_mismatches:
                sam_three_list.append(three_sam_entry)
            three_sam_line = three_prime_sam_file.readline()

        # Niether the 5' nor 3' sam is at the prev_jct_ind
        else:
            prev_consensus = denovo_junctions[prev_jct_ind].consensus
            best_splices[prev_jct_ind] = get_best_splice(sam_five_list,sam_three_list,prev_consensus,max_mismatches)
            sam_five_list = []
            sam_three_list = []
            prev_jct_ind = five_jct_ind

    # Have to push the last jct lists into the shared_dict
    prev_consensus = denovo_junctions[prev_jct_ind].consensus
    best_splices[prev_jct_ind] = get_best_splice(sam_five_list,sam_three_list,prev_consensus,max_mismatches)
    
    # Close the 5' and 3' sam files
    five_prime_sam_file.close()
    three_prime_sam_file.close()

    # Loop through the jcts assigning the splice site info
    for jct_ind in range(len(denovo_junctions)):
        jct = denovo_junctions[jct_ind]

        # Default splice options if splice sites not found
        splice_ind = len(jct.consensus)/2
        splice_type = "None"
        splice_gap = len(jct.consensus)

        # Use calculted best splice if was previously found
        if jct_ind in best_splices:
            splice_ind,splice_type,splice_gap = best_splices[jct_ind]

        # Assign the type of jct found Full is most confident, Half is 2nd most, None is 3rd most
        jct.splice_type = splice_type
        jct.splice_ind = splice_ind
        jct.splice_gap = splice_gap
        jct.upstream_seq = jct.consensus[:jct.splice_ind]
        jct.downstream_seq = jct.consensus[jct.splice_ind:]
        jct.find_start_stop()


#######################
#   Get Best Splice   #
#######################
# This is a helper function for the splice site finder
# If there are multiple found splice, return the best from the shared sams dict from the two sam lists
def get_best_splice(sam_five_list,sam_three_list,consensus,max_mismatches):
    shared_dict = {}
    id_dict = {}
    best_sam_five = None
    best_sam_five_len = 0
    best_sam_three = None
    best_sam_three_len = 0
    for sam_five in sam_five_list:
        id_dict[sam_five.read_id] = sam_five
        if len(sam_five.seq) > best_sam_five_len:
            best_sam_five_len = len(sam_five.seq)
            best_sam_five = sam_five
    for sam_three in sam_three_list:
        if len(sam_three.seq) > best_sam_three_len:
            best_sam_three_len = len(sam_three.seq)
            best_sam_three = sam_three
        if sam_three.read_id in id_dict:
            sam_five = id_dict[sam_three.read_id]
            shared_dict[sam_three.read_id] = [sam_five,sam_three]

    # This is the gap between the 5' and 3' mapping pieces
    # Will be 0 if there is a perfect map:
    #       0 1 2 3
    #               4 5 6 7
    # Or the number of missing bases in a non-perfect map:
    #       0 1 2 3
    #                   6 7
    gap_len = len(consensus)-best_sam_five_len-best_sam_three_len

    # If there is at least one shared perfect splice ind find the one w/ least mismatches
    if len(shared_dict) > 0:
        best_key = ""
        min_mismatches = max_mismatches+1
        for key in shared_dict:
            sam1,sam2 = shared_dict[key]
            num_mismatches = sam1.num_mismatches+sam2.num_mismatches
            if num_mismatches < min_mismatches:
                best_key = key
                min_mismatches = num_mismatches
        best_sam1,best_sam2 = shared_dict[best_key]
        return len(best_sam1.seq),"Full",gap_len
    else:
        # There is a mapping for the left and right pieces, although there is space in between
        if best_sam_five and best_sam_three:
            five_ind = best_sam_five_len
            three_ind = len(consensus)-best_sam_three_len
            splice_ind = (five_ind+three_ind)/2
            return splice_ind,"Gapped",gap_len
        # If there are no mappings for the three-side piece
        elif best_sam_five:
            splice_ind = best_sam_five_len
            return splice_ind,"Five_Only",gap_len
        # If there are no mappings for the five-side piece
        elif best_sam_three:
            splice_ind = len(consensus)-best_sam_three_len
            return splice_ind,"Three_Only",gap_len
        # There are no 5' and 3' mappings (I think this case shouldn't happen...)
        else:
            sys.stderr.write("ERROR: Somehow got to jct w/ no mappings past filter. Gap len:"+str(gap_len)+" cons len: "+str(len(consensus)))
            splice_ind = len(consensus)/2
            return splice_ind,"None",gap_len


#####################
#   Generate GTFS   #
#####################
# Helper function to generate a list of gtf objects from a gtf_file
def generate_gtfs(gtf_file_name):
    # Group the gtfs by name using a dictionary
    gtf_dict_by_name = {}
    with open(gtf_file_name,"r") as gtf_file:
        for gtf_line in gtf_file.readlines():
            gtf = GTFEntry(gtf_line)
            key = gtf.chromosome+"_"+gtf.gene_name
            if key not in gtf_dict_by_name:
                gtf_dict_by_name[key] = [gtf]
            else:
                gtf_dict_by_name[key].append(gtf)

    # For each key, convert all the gtfs into a single gtf w/ min and max start and stop pos
    gtfs = []
    for key in gtf_dict_by_name:
        initial_gtf = gtf_dict_by_name[key][0]
        start_pos = initial_gtf.start
        stop_pos = initial_gtf.stop
        for gtf in gtf_dict_by_name[key][1:]:
            if gtf.start < start_pos: start_pos = gtf.start
            if gtf.stop > stop_pos: stop_pos = gtf.stop
        initial_gtf.start = start_pos
        initial_gtf.stop = stop_pos
        gtfs.append(initial_gtf)
    return gtfs


########################
#   Get JCT GTF info   #
########################
# Reads the entire gtf line by line and checks against every junction
# Reading the genes in as all exons to try and get start and stop site
def get_jct_gtf_info(junctions,gtfs):
    # Loop through the collapsed gtfs to see if the junction is in the range
    for gtf in gtfs:
        for junction in junctions:
            if gtf.chromosome == junction.upstream_chromosome:
                if gtf.start <= junction.upstream_start <= gtf.stop or gtf.start <= junction.upstream_stop <= gtf.stop:
                    if junction.upstream_genes == "no_genes": junction.upstream_genes = gtf.gene_name
                    else: junction.upstream_genes += ","+gtf.gene_name
            if gtf.chromosome == junction.downstream_chromosome:
                if gtf.start <= junction.downstream_start <= gtf.stop or gtf.start <= junction.downstream_stop <= gtf.stop:
                    if junction.downstream_genes == "no_genes": junction.downstream_genes = gtf.gene_name
                    else: junction.downstream_genes += ","+gtf.gene_name
            junction.check_fusion()


########################
#   Get SAM GTF Info   #
########################
# Assigns gene names to a single sam_entry. Doesn't return, just updates the argument
def get_sam_gtf_info(gtfs,sam_entry):
    for gtf in gtfs:
        if gtf.chromosome == sam_entry.chromosome:
            if gtf.start <= sam_entry.upstream_pos <= gtf.stop:
                sam_entry.genes.append(gtf.gene_name)
            elif gtf.start <= sam_entry.downstream_pos <= gtf.stop:
                sam_entry.genes.append(gtf.gene_name)

# Collapse junctions by splice site
# Currently keeps linear and non-linear separated even if they share a splice site
def collapse_junctions(junctions):
    splice_to_jct_dict = {}
    for junction in junctions:
        splice_key = junction.upstream_chromosome+"-"+junction.downstream_chromosome+":"+str(junction.splice_site)+":"+str(junction.linear)
        if splice_key not in splice_to_jct_dict:
            splice_to_jct_dict[splice_key] = [junction]
        else:
            splice_to_jct_dict[splice_key].append(junction)
            
    collapsed_junctions = []
    for jcts_by_splice in splice_to_jct_dict.itervalues():
        first_junction = jcts_by_splice[0]
        for shared_jct in jcts_by_splice[1:]:
            first_junction.add_constitutive_junction(shared_jct)
        first_junction.combine_constitutive_junctions()
        collapsed_junctions.append(first_junction)
    return collapsed_junctions


####################
#   Assign Class   #
####################
# Assigns a pair of R1 and R2 to a class based on certain factors. These are both SAMEntry objects.
# Possible classes. R1 always maps to a denovo jct, and R2 somewhere else.
# TODO update the logic to allow 'Fusion' classification
# TODO use regex to parse the chromosome rather than all this messy splitting (it looks terrible)
# [1] Linear
# [2] Linear Anomally
# [3] Circular
# [4] Circular Anomally
# [5] None <-- kind of in the place of fusions for now
def assign_class(sam_R1,sam_R2):
    jct_chrom_1 = sam_R1.chromosome.split("|_|")[0].split("|")[1]
    jct_chrom_2 = sam_R1.chromosome.split("|_|")[1].split("|")[0]

    # If the jct splices 2 chromosomes together just skip it for now
    if jct_chrom_1 != jct_chrom_2:
        return "None"
    # If the R1 and R2 are on different chromosomes just skip it for now
    if jct_chrom_1 != sam_R2.chromosome:
        return "None"

    span = int(sam_R1.chromosome.split("|_|")[1].split("|")[5].split(":")[1])

    #Linear case
    if span > 0:
        if sam_R1.strand != sam_R2.strand:
            return "Linear"
        else:
            return "Linear_Anomaly"

    #Circular case
    else:
        start_1 = int(sam_R1.chromosome.split("|_|")[0].split("|")[3].split(":")[1])
        stop_1 = int(sam_R1.chromosome.split("|_|")[0].split("|")[4].split(":")[1])
        start_2 = int(sam_R1.chromosome.split("|_|")[1].split("|")[2].split(":")[1])
        stop_2 = int(sam_R1.chromosome.split("|_|")[1].split("|")[3].split(":")[1])
        if sam_R1.strand == sam_R2.strand:
            return "Circular_Anomaly"
        elif start_2 <= sam_R2.upstream_pos <= stop_1:
            return "Circular"
        else:
            return "Circular_Anomaly"


############################
#   Write GLM Class File   #
############################
# Simple function to print out GLM class file in the right format
def write_glm_class_file(class_file_name,sam_list):
    header = ""
    header += "id\t"
    header += "class\t"
    header += "pos\t"
    header += "qual\t"
    header += "aScore\t"
    header += "numN\t"
    header += "readLen\t"
    header += "junction\t"
    header += "strand\t"
    header += "posR2\t"
    header += "qualR2\t"
    header += "aScoreR2\t"
    header += "numNR2\t"
    header += "readLenR2\t"
    header += "junctionR2\t"
    header += "strandR2\n"
    with open(class_file_name,"w") as class_file:
        class_file.write(header)
        for read_pair in sam_list:
            sam_R1,sam_R2,pair_class,sam_R1_genes,sam_R2_genes = read_pair

            #Add general info to the out_line
            out_line = ""
            out_line += str(sam_R1.read_id.split("/")[0])+"\t"
            out_line += str(pair_class)+"\t"

            #Add R1 info to the out_line
            out_line += str(sam_R1.upstream_pos)+"\t"
            out_line += str(sam_R1.mapping_quality)+"\t"
            out_line += str(sam_R1.alignment_score)+"\t"
            out_line += str(sam_R1.num_Ns)+"\t"
            out_line += str(len(sam_R1.seq))+"\t"
            out_line += str(sam_R1.junction())+"|"+sam_R1_genes+"\t"
            out_line += str(sam_R1.strand)+"\t"

            #Add R2 info to the out_line
            out_line += str(sam_R2.upstream_pos)+"\t"
            out_line += str(sam_R2.mapping_quality)+"\t"
            out_line += str(sam_R2.alignment_score)+"\t"
            out_line += str(sam_R2.num_Ns)+"\t"
            out_line += str(len(sam_R2.seq))+"\t"
            out_line += str(sam_R2.junction())+"|"+sam_R2_genes+"\t"
            out_line += str(sam_R2.strand)+"\n"

            #Write the built up out_line to the glm input class file
            class_file.write(out_line)


##################
#   Write Time   #
##################
# Helper function to write out the timing of something
# Takes in the message, a start time in seconds, and a timer_file_path
# Appends to the timer file by default, first call should overwrite
def write_time(message,start_time,timer_file_path,append=True,uniform_len=70):
    seconds_duration = float(time.time()-start_time)
    minutes_duration = seconds_duration/60
    hours_duration = seconds_duration/3600
    seconds_str= ("{:3.2f}".format(seconds_duration)).rjust(5," ")
    minutes_str= ("{:3.2f}".format(minutes_duration)).rjust(5," ")
    hours_str= ("{:3.2f}".format(hours_duration)).rjust(5," ")

    if len(message) < uniform_len:
        message += " "*(uniform_len-len(message))
    else:
        message = message[:uniform_len]
    time_out_str  = message+":    "
    time_out_str += seconds_str+":seconds    "
    time_out_str += minutes_str+":minutes    "
    time_out_str += hours_str+":hours\n"
    timer_file = open(timer_file_path,"a") if append else open(timer_file_path,"w")
    timer_file.write(time_out_str)
    timer_file.close()


##########################
#   Reverse Compliment   #
##########################
# Just a little helper function to give the reverse compliment of a sequence
def reverse_compliment(seq):
    comp_dict = {"A":"T","a":"t","T":"A","t":"a",
                     "C":"G","c":"g","G":"C","g":"c",
                     "N":"N","n":"n"}
    rev_comp_seq = "".join([comp_dict[base] for base in seq])[::-1]
    return rev_comp_seq



