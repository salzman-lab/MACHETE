#Goal is to take a large fasta file and a naive report file
#then make a new fasta file of just the naive reports (filtering)
import time
import sys
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--directory", required=True, help="FJ directory")
parser.add_argument("-s", "--stem", required=True, help="unique file stem")
args=parser.parse_args()


#Setup the paths
base_path = args.directory
name = args.stem


large_fasta_name = base_path+"fasta/"+name+"_FarJunctions.fa"
naive_report_name = base_path+"reports/"+name+"_naive_report.txt"
output_fasta_name = base_path+"fasta/"+name+"_filtered_FarJunctions.fa"

#Read in the desired headers
start_time = time.time()
desired_headers = []
with open(naive_report_name,"r") as naive_report:
    for naive_line in naive_report:
        if "@" not in naive_line:
            desired_headers.append(naive_line.split("\t")[0])

#Walk through the large fasta file keeping only the ones
#that have headers in the desired headers

output_file = open(output_fasta_name,"w")
with open(large_fasta_name,"r") as large_fasta:
    header_line = large_fasta.readline()
    while header_line:
        check_header = header_line.strip()[1:]
        seq_line = large_fasta.readline()
        if check_header in desired_headers:
            output_file.write(header_line)
            output_file.write(seq_line)
        header_line = large_fasta.readline()

output_file.close()
print (time.time()-start_time)/60,"mins"
