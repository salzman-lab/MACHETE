
fusion_file_name = "gillian_headers_fusions.fasta"
padded_file_name = "padded_gillian.fasta"

with open(fusion_file_name,"r") as fusion_file:
    with open(padded_file_name,"w") as padded_file:
        fusion_line = fusion_file.readline()
        while fusion_line:
            if ">" not in fusion_line:
                padded_file.write("N"*50+fusion_line.strip()+"N"*50+"\n")
            else:
                padded_file.write(fusion_line)
            fusion_line = fusion_file.readline()
