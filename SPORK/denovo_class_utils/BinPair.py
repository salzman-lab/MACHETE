#BinPair class
class BinPair(object):
    __slots__ = ["five_prime_SAM","three_prime_SAM","five_prime_bin","three_prime_bin",
                 "five_prime_strand","three_prime_strand","bin_pair","five_prime_chr","three_prime_chr"]
    def __init__(self,five_prime_SAM,three_prime_SAM,five_prime_bin,three_prime_bin):
        self.five_prime_SAM = five_prime_SAM
        self.three_prime_SAM = three_prime_SAM
        self.five_prime_bin = int(five_prime_bin)
        self.three_prime_bin = int(three_prime_bin)
        self.five_prime_strand = five_prime_SAM.strand
        self.three_prime_strand = three_prime_SAM.strand
        self.five_prime_chr = five_prime_SAM.chromosome
        self.three_prime_chr = three_prime_SAM.chromosome
        self.bin_pair = self.five_prime_chr+":"+str(self.five_prime_bin)+"_"+self.three_prime_chr+":"+str(self.three_prime_bin)+"_("+self.five_prime_strand+","+self.three_prime_strand+")"

    #Schematic of what the reverse compliment aims to do:
    #================
    #Before flipping:
    #================
    #
    #(+ strand)-----------------------------------------
    #
    #             (*)  (#)                (*)  (#)
    #              ======                  ======
    #(- strand)----| 3' |------------------| 5' |-------
    #              ======                  ======
    #
    #(*)'s are the upstream positions and (#)'s are the downstream positions
    #
    #===============
    #After flipping:
    #===============
    #
    #             (*)  (#)                (*)  (#)
    #              ======                  ======
    #(+ strand)----| 5' |------------------| 3' |-------
    #              ======                  ======
    #
    #(- strand)-----------------------------------------
    #
    #So the 5' box got the orig 3' positions and should have the rev comp seq of the orig 3' seq
    #the same thing happens to the 3' box
    #TODO I should really test that this is working properly w/ unittests
    def take_reverse_compliment(self):
        #Swap the strand both in the underlying SAMs and the other members
        #TODO make the members actually just call the underlying SAMs each time
        self.five_prime_SAM.strand = "-" if self.five_prime_SAM.strand == "+" else "+"
        self.three_prime_SAM.strand = "-" if self.three_prime_SAM.strand == "+" else "+"
        self.five_prime_strand = "-" if self.five_prime_strand == "+" else "+"
        self.three_prime_strand = "-" if self.three_prime_strand == "+" else "+"

        #Swap the positions of the five and three prime SAMs
        hold_five_prime_upstream_pos = self.five_prime_SAM.upstream_pos
        hold_five_prime_downstream_pos = self.five_prime_SAM.downstream_pos
        self.five_prime_SAM.upstream_pos = self.three_prime_SAM.upstream_pos
        self.five_prime_SAM.downstream_pos = self.three_prime_SAM.downstream_pos
        self.three_prime_SAM.upstream_pos = hold_five_prime_upstream_pos
        self.three_prime_SAM.downstream_pos = hold_five_prime_downstream_pos

        #Take the reverse compliment of the 5' and 3' seqs
        #The list comprehensions are complicated but I like doing it in one line:
        #   The [rev_comp_dict[base] for base in self.five_prime_SAM.seq] builds a list of complimentary bases
        #   The "".join takes that list and turns it into a string
        #   The [::-1] at the very end reverses the string to turn the compliment string into the rev comp string
        rev_comp_dict = {"A":"T","a":"t","T":"A","t":"a",
                         "C":"G","c":"g","G":"C","g":"c",
                         "N":"N","n":"n"}
        rev_comp_5_prime_seq = "".join([rev_comp_dict[base] for base in self.five_prime_SAM.seq])[::-1]
        rev_comp_3_prime_seq = "".join([rev_comp_dict[base] for base in self.three_prime_SAM.seq])[::-1]

        #Then put the orig 3' rev comp into the new 5' and vice versa
        self.five_prime_SAM.seq = rev_comp_3_prime_seq
        self.three_prime_SAM.seq = rev_comp_5_prime_seq

        #Switch the chromosomes too in case this is a fusion
        hold_chromosome = self.five_prime_SAM.chromosome
        self.five_prime_SAM.chromosome = self.three_prime_SAM.chromosome
        self.three_prime_SAM.chromosome = hold_chromosome

        #Using this in a list comprehension, so I want it to return itself
        return self

    def __str__(self):
        return "Bin Pair: "+self.bin_pair+" Left: "+self.five_prime_SAM.read_id+" Right: "+self.three_prime_SAM.read_id+"\n"

    def __lt__(self,other):
        return self.bin_pair < other.bin_pair

    def save(obj):
        return (obj.__class__, obj.__dict__)

    def load(cls, attributes):
        obj = cls.__new__(cls)
        obj.__dict__.update(attributes)
        return obj

