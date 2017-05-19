#General Imports
import sys

#Currently gives information regarding the counts of:
# (1) fastas representing bin pairs from diff. chroms
# (2) fastas representing bin pairs separated by > 1000nt
# (3) fastas representing bin pairs separated by <= 1000n
#Example fasta header to parse
#>|chr1|HFM1|start:91852903|stop:91852938|strand:-|_|chr1|HFM1|start:91852874|stop:91852926|strand:+|span:-64|score:0.464935064935|num:15| 
def find_shared_chrs(input_fasta_name):
    path = "/".join(input_fasta_name.split("/")[:-1])+"/"
    fasta_name = input_fasta_name.split("/")[-1]
    shared_chr_file = open(path+"shared_chr_"+fasta_name,"w")
    not_shared_chr_file = open(path+"not_shared_chr_"+fasta_name,"w")
    shared_chrs = []
    not_shared_chrs = []
    with open(input_fasta_name,"r") as input_fasta:
        fasta_lines = input_fasta.readlines()
        for fasta_line_ind in range(len(fasta_lines)):
            fasta_line = fasta_lines[fasta_line_ind]
            #Only interested in the header lines not sequences
            if ">" in fasta_line:
                first_bin = fasta_line.split("|_")[0]
                first_chr = first_bin.split("|")[1]
                first_genes = first_bin.split("|")[2].split(",")
                first_start = int(first_bin.split("|")[3].split("start:")[1])
                first_stop = int(first_bin.split("|")[4].split("stop:")[1])

                second_bin = fasta_line.split("|_")[1]
                second_chr = second_bin.split("|")[1]
                second_genes = second_bin.split("|")[2].split(",")
                second_start = int(second_bin.split("|")[3].split("start:")[1])
                second_stop = int(second_bin.split("|")[4].split("stop:")[1])

                span = int(second_bin.split("|")[6].split("span:")[1])

                if first_genes == second_genes:
                    shared_chrs.append(fasta_line)
                else:
                    not_shared_chrs.append(fasta_line)

        shared_chrs.sort()
        not_shared_chrs.sort()
        for shared_chr in shared_chrs:
            shared_chr_file.write(shared_chr)
        for not_shared_chr in not_shared_chrs:
            not_shared_chr_file.write(not_shared_chr)

    shared_chr_file.close()
    not_shared_chr_file.close()

#If running from command line then take in the junction fasta file name
#Has a default which might be easier to change
if __name__ == "__main__":
    input_fasta_name = "../checking_superset/novel_junctions_still_unaligned_140701_HAVERS_0463_AC4GKDACXX_L8_CGATGT_R1.fasta"
    find_shared_chrs(input_fasta_name)
