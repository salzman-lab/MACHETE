import sys
input_file_name = "unaligned_MS2_20_S20_L002_R1_denovo_report.txt"

#Example lines
#chr5|KIAA1191:175779624|KIAA1191:175777740|reg|-    1   0   0   0   0   0   0.958759573772  (0,0),
#chr8|ST3GAL1:134478136|ST3GAL1:134477200|reg|-  2   0   0   0   0   0   0.988514264524  (0,0),(-3,-9),

with open(input_file_name,"r") as input_file:
    all_lines = input_file.readlines()
    all_lines = all_lines[2:]
    for line in all_lines:
        split_by_pipe = line.split("|")
        chromosome = split_by_pipe[0]
        pos_1 = split_by_pipe[1].split(":")[1]
        pos_2 = split_by_pipe[2].split(":")[1]
        print chromosome+"\t"+pos_1+"\t"+pos_2
