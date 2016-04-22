fasta_file_name = "novel_junctions.fasta"
fusion_file_name = "novel_fusions.fasta"

with open(fasta_file_name,"r") as fasta_file:
    with open(fusion_file_name,"w") as fusion_file:
        fasta_line = fasta_file.readline()
        while fasta_line:
            if "True" in fasta_line:
                fusion_file.write(fasta_line)
                fasta_line = fasta_file.readline()
                fusion_file.write(fasta_line)
            fasta_line = fasta_file.readline()
