#Class to store GTF entries
#TODO give this better documentation

#Imports
import re
import sys

class GTFEntry(object):
    __slots__ = ["chromosome","source","feature","start","stop","score","strand","frame","gene_name"]
    def __init__(self,gtf_line='chr\tsrc\tfeat\t-1\t-1\t-1\t+\t0\tgene_name "default";'):
        split_gtf_line = gtf_line.split("\t")
        self.chromosome = split_gtf_line[0]
        self.source = split_gtf_line[1]
        self.feature = split_gtf_line[2]
        self.start = int(split_gtf_line[3])
        self.stop = int(split_gtf_line[4])
        self.score = split_gtf_line[5]
        self.strand = split_gtf_line[6]
        self.frame = split_gtf_line[7]
        group_info = split_gtf_line[8]
        gene_name_pattern = re.compile('gene_name "(.*?)";')
        gene_name = gene_name_pattern.findall(group_info)
        if len(gene_name) == 0:
            gene_id_pattern = re.compile('gene_id "(.*?)";')
            gene_name = gene_id_pattern.findall(group_info)
            if len(gene_name) == 0:
                print "Error in gtf init. No found gene_name or gene_id"
                sys.stderr.write("Error in gtf init. No found gene_name or gene_id")
                sys.exit(1)
        self.gene_name = gene_name[0]

    def __str__(self):
        ret_str = ""
        ret_str += "chromosome: "+self.chromosome+"\t"
        ret_str += "name: "+self.gene_name+"\t"
        ret_str += "start: "+str(self.start)+"\t"
        ret_str += "stop: "+str(self.stop)+"\t"
        ret_str += "strand: "+self.strand+"\t"
        return ret_str

    def __lt__(self,other):
        return self.start < other.start

    def save(obj):
        return (obj.__class__, obj.__dict__)

    def load(cls, attributes):
        obj = cls.__new__(cls)
        obj.__dict__.update(attributes)
        return obj

