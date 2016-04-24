#Mapped Read class
class SAMEntry(object):
    __slots__ = ["read_id","strand","chromosome","upstream_pos","downstream_pos","seq","seq_quality",
                 "num_gaps","num_mismatches","num_Ns","mapping_quality","alignment_score","genes"]
    def __init__(self,full_line):
        split_line = full_line.split("\t")
        #Example bowtie2 SAM line w/ annotation:
        #[0:read_id               ] K00180:68:H5CF7BBXX:3:1122:4422:11442/5_prime
        #[1 :strand (0/16)        ] 0
        #[2 :chromosome           ] chr21
        #[3 :position             ] 9827122
        #[4 :map quality          ] 42  
        #[5 :CIGAR string         ] 33M 
        #[6 :name of mate         ] *
        #[7 :pos of mate          ] 0
        #[8 :template len         ] 0
        #[9 :sequence             ] CTTTGGTCGCTCGCTCCTCTCCTACTTGGATAA
        #[10:quality string       ] <AAAFJJJFJJJJJFJFJJFJJJJJAJJFAFFJ
        #[11:Alignment score      ] AS:i:0
        #[12:number of N's        ] XN:i:0
        #[13:number of mismatches ] XM:i:0
        #[14:number gap opens     ] XO:i:0
        #[15:number gap extensions] XG:i:0
        #[16:edit distance        ] NM:i:0
        #[17:string for mismatches] MD:Z:33
        #[18:whether or not paired] YT:Z:UU
        self.read_id         = split_line[0]
        try:
            self.strand          = "+" if split_line[1] == "0" else "-"
        except:
            print full_line
        self.chromosome      = split_line[2]
        self.upstream_pos    = int(split_line[3])
        self.mapping_quality = int(split_line[4])
        self.seq             = split_line[9]
        self.seq_quality     = split_line[10]
        self.alignment_score = int(split_line[11].split(":")[-1])
        self.num_Ns          = int(split_line[12].split(":")[-1])
        self.num_mismatches  = int(split_line[13].split(":")[-1])
        self.num_gaps        = int(split_line[14].split(":")[-1])

        self.downstream_pos = self.upstream_pos+len(self.seq)
        self.genes = []

    def junction(self):
        out_str = ""
        out_str += "jct|"+self.chromosome.split("|")[1] if "|" in self.chromosome else self.chromosome
        out_str += "|"+str(self.upstream_pos)
        out_str += "|"+str(self.downstream_pos)
        out_str += "|"+str(self.strand)
        return out_str

    def __str__(self):
        return self.read_id+"\t"+self.seq+"\n"

    def __lt__(self,other):
        same_chr = self.chromosome == other.chromosome
        if not same_chr: return self.chromosome < other.chromosome
        else: return self.upstream_pos < other.upstream_pos

    def save(obj):
        return (obj.__class__, obj.__dict__)

    def load(cls, attributes):
        obj = cls.__new__(cls)
        obj.__dict__.update(attributes)
        return obj

