#!/usr/bin/python2.6
#Set python version

#General Imports
import matplotlib.pyplot as plt
import sys
               
def plot_score_vs_num(input_jct_name):
    nums = []
    scores = []
    with open(input_jct_name) as input_jct:
        for junction_line in input_jct.readlines():
            if ">" not in junction_line: continue
            jct_dict = parse_fasta_header(junction_line)
            nums.append(jct_dict["num"])
            scores.append(jct_dict["score"])
    plt.plot(nums,scores,"ob")
    plt.xlabel("Number in Bin Pair Group")
    plt.ylabel("Consensus Mismatch Score")
    plt.axis([0,400,0,1])
    plt.show()
    print "Should be showing"

#>|chr1|HFM1|start:91852903|stop:91852938|strand:-|_|chr1|HFM1|start:91852874|stop:91852926|strand:+|span:-64|score:0.464935064935|num:15| 
def hist_found_site_vs_score(input_jct_name):
    scores_with_splices = []
    with open(input_jct_name) as input_jct:
        for junction_line in input_jct.readlines():
            if ">" not in junction_line: continue
            jct_dict = parse_fasta_header(junction_line)
            scores_with_splices.append(float(jct_dict["score"]))
    plt.hist(scores_with_splices,20)
    plt.xlabel("Consensus Mismatch Score")
    plt.ylabel("Number Found Splice Sites")
    plt.show()

def parse_fasta_header(fasta_line):
    jct_dict = {}
    split_fasta = fasta_line.split("|")
    chr1 = split_fasta[1]
    jct_dict["chr1"] = chr1
    genes1 = split_fasta[2].split(",")
    jct_dict["genes1"] = genes1
    start1 = int(split_fasta[3].split(":")[1])
    jct_dict["start1"] = start1
    stop1 = int(split_fasta[4].split(":")[1])
    jct_dict["stop1"] = stop1
    strand1 = split_fasta[5].split(":")[1]
    jct_dict["strand1"] = strand1

    chr2 = split_fasta[7]
    jct_dict["chr2"] = chr2
    genes2 = split_fasta[8].split(",")
    jct_dict["genes2"] = genes2
    start2 = int(split_fasta[9].split(":")[1])
    jct_dict["start2"] = start2
    stop2 = int(split_fasta[10].split(":")[1])
    jct_dict["stop2"] = stop2
    strand2 = split_fasta[11].split(":")[1]
    jct_dict["strand2"] = strand2

    span = int(split_fasta[12].split(":")[1])
    jct_dict["span"] = span
    score = float(split_fasta[13].split(":")[1])
    jct_dict["score"] = score
    num = int(split_fasta[14].split(":")[1])
    jct_dict["num"] = num
    return jct_dict

#def num_denovo_jcts_per_chr(input_jct)name):
    

#If running from command line then take in the junction file name
#Has a default which might be easier to change
if __name__ == "__main__":
    path = "../test_python_perl/compare/orig/unaligned/forDenovoIndex/test_temp_out_dir/"
    fq_name = "unaligned_MS3_28_S28_L003_R1"
    input_jct_name = path+fq_name+"/novel_junctions_"+fq_name+".fasta"
    if len(sys.argv) == 2:
        print "Using input jct name"
        input_jct_name = sys.argv[2]
    else:
        print "Using defualt jct name"

    #Make different plots
    plot_score_vs_num(input_jct_name)
    #hist_found_site_vs_score(input_jct_name)
    print "Got here"
