#Fastq entry class
#TODO better document this class

#Imports
import sys

class FastQEntry(object):
    __slots__ = ["read_id","seq","plus_line","quality"]
    def __init__(self,read_id,seq,plus_line,quality):
        self.read_id = read_id
        self.seq = seq
        self.plus_line = plus_line
        self.quality = quality
        self.clean()

    def get_edge_thirds(self):
        third_len = len(self.seq)/3
        if third_len*2 > len(self.seq):
            print "Read is too short to split"
            sys.exit(1)
        five_prime_seq = self.seq[:third_len]
        three_prime_seq = self.seq[2*third_len:]
        five_prime_read = FastQEntry(self.read_id+"/5_prime",five_prime_seq,self.plus_line,self.quality[:third_len])
        three_prime_read = FastQEntry(self.read_id+"/3_prime",three_prime_seq,self.plus_line,self.quality[2*third_len:])
        return five_prime_read,three_prime_read
        
    def clean(self):
        self.read_id = self.read_id.replace(" ","_").replace("\t","_").replace("\n","")
        self.seq = self.seq.replace("U","T").replace("\n","")
        self.plus_line = self.plus_line.replace(" ","_").replace("\n","")
        self.quality = self.quality.replace("\n","")
    
    def __str__(self):
        ret_str = ""
        ret_str += self.read_id+"\n"
        ret_str += self.seq+"\n"
        ret_str += self.plus_line+"\n"
        ret_str += self.quality+"\n"
        return ret_str

    def __lt__(self,other):
        return self.read_id < other.read_id

    def save(obj):
        return (obj.__class__, obj.__dict__)

    def load(cls, attributes):
        obj = cls.__new__(cls)
        obj.__dict__.update(attributes)
        return obj

