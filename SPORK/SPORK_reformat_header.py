import re
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-f", "--fastadir", required=True, help="directory of fasta file to convert")
parser.add_argument("-s", "--stem", required=True, help ="unique ID for file")
args=parser.parse_args()

if args.fastadir[-1]!="/":
    args.fastadir += "/"

fasta_file_name = args.fastadir+args.stem+"_SPORK_preprocess_Junctions.fa"
fusion_file_name = args.fastadir+args.stem+"_novel_fusions.fa"
cleaned_file_name = args.fastadir+args.stem+"_SPORK_Junctions.fa"

with open(fasta_file_name,"r") as fasta_file:
    with open(fusion_file_name,"w") as fusion_file:
        fasta_line = fasta_file.readline()
        while fasta_line:
            if "True" in fasta_line:
                fusion_file.write(fasta_line)
                fasta_line = fasta_file.readline()
                fusion_file.write(fasta_line)
            fasta_line = fasta_file.readline()




chrom1_pattern = re.compile("chromosome1:(.*?)\|")
chrom2_pattern = re.compile("chromosome2:(.*?)\|")
strand1_pattern = re.compile("strand1:(.*?)\|")
strand2_pattern = re.compile("strand2:(.*?)\|")
genes1_pattern = re.compile("genes1:(.*?)\|")
genes2_pattern = re.compile("genes2:(.*?)\|")
start1_pattern = re.compile("start1:(.*?)\|")
start2_pattern = re.compile("start2:(.*?)\|")

with open(fusion_file_name,"r") as fusion_file:
    with open(cleaned_file_name,"w") as cleaned_file:
        fusion_line = fusion_file.readline()
        while fusion_line:
            if ">" in fusion_line:
                chrom1 = re.findall(chrom1_pattern,fusion_line)[0]
                chrom2 = re.findall(chrom2_pattern,fusion_line)[0]
                strand1 = re.findall(strand1_pattern,fusion_line)[0]
                strand2 = re.findall(strand2_pattern,fusion_line)[0]
                genes1 = re.findall(genes1_pattern,fusion_line)[0]
                genes2 = re.findall(genes2_pattern,fusion_line)[0]
                start1 = re.findall(start1_pattern,fusion_line)[0]
                start2 = re.findall(start2_pattern,fusion_line)[0]

                new_header = ">"
                new_header += chrom1+":"+genes1+":"+start1+":"+strand1+"|"
                new_header += chrom2+":"+genes2+":"+start2+":"+strand2+"|"
                new_header += "fusion\n"
                cleaned_file.write(new_header)
                
            else:
                cleaned_file.write(fusion_line)
            fusion_line = fusion_file.readline()

                #|chromosome1:chr11|genes1:EWSR1|start1:29683023|stop1:29683123|strand1:+|_|chromosome2:chr22|genes2:FLI1|start2:128675263|stop2:128675328|strand2:+|splice:100|span:98992140|score:0.0733915636983|fusion:True|num:63|splice-gap:65|splice-type:Five_Only|
                #chr11:EWSR1:29683123:+|chr11:FLI:128675263:+|fusion
