fastq_R1 = "still_unaligned_140701_HAVERS_0463_AC4GKDACXX_L8_GATCAG_R1.fq"
fastq_R2 = "still_unaligned_140701_HAVERS_0463_AC4GKDACXX_L8_GATCAG_R2.fq"

R1_read_ids = []
with open(fastq_R1,"r") as fastq_file:
    fastq_lines = fastq_file.readlines()
    for fastq_line_ind in range(0,len(fastq_lines),4):
        R1_read_ids.append(fastq_lines[fastq_line_ind])
R1_read_ids.sort()

R2_read_ids = []
with open(fastq_R2,"r") as fastq_file:
    fastq_lines = fastq_file.readlines()
    for fastq_line_ind in range(0,len(fastq_lines),4):
        R2_read_ids.append(fastq_lines[fastq_line_ind])
R2_read_ids.sort()

shared_read_ids = []
R1_ind = 0
R2_ind = 0
while R1_ind < len(R1_read_ids) and R2_ind < len(R2_read_ids):
    print str(R1_ind)+"/"+str(len(R1_read_ids))+"------"+str(R2_ind)+"/"+str(len(R2_read_ids))
    if R1_read_ids[R1_ind] < R2_read_ids[R2_ind]:
        R1_ind += 1
    elif R1_read_ids[R1_ind] > R2_read_ids[R2_ind]:
        R2_ind += 1
    else:
        shared_read_ids.append(R1_read_ids[R1_ind])
        R1_ind += 1
        R2_ind += 1

for shared_read_id in shared_read_ids:
    print shared_read_id
print len(shared_read_ids)
