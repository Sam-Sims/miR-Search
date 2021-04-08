from Bio import SeqIO


class MicroRNA:
    def __init__(self, seq):
        self.mir_seq = seq
        print("Creating miRNA object with sequence: " + self.mir_seq)

    def trim_mir(self):
        print("Trimming miRNA...")
        self.mir_seq = self.mir_seq[:8]
        print("Done!")

    def reverse_complement(self):
        print("Finding the reverse complement of ", self.mir_seq)
        self.mir_seq = self.mir_seq.reverse_complement()
        print("Done!")

    def back_transcribe(self):
        print("Back transcribing RNA sequence...")
        self.mir_seq = self.mir_seq.back_transcribe()
        print("Done!")

    def auto_process(self):
        print("Auto processing miRNA sequence...")
        self.trim_mir()
        self.reverse_complement()
        self.back_transcribe()
        print("Processing complete!")
        print("Final sequence: ", self.mir_seq)

    def return_seq(self):
        return self.mir_seq

    def find_6mer(self):
        return self.mir_seq[1:7]

    def find_7mera1(self):
        return self.mir_seq[1:7] + "A"

    def find_7merm8(self):
        return self.mir_seq[:7]

    def find_8mer(self):
        return self.mir_seq[:7] + "A"


class MicroRNASearch:
    def __init__(self, utr):
        self.utr = utr


class FASTAParse:
    def read_mir(self, path):
        for fasta_sequence in SeqIO.parse(open(path), 'fasta'):
            print("Reading miRNA at " + path + "...")
            print("miRNA: ",  fasta_sequence.seq)
            return fasta_sequence.seq

    def read_3utr(self, path):
        for fasta_sequence in SeqIO.parse(path, "fasta"):
            print("Reading 3' UTR region... ")
            print(fasta_sequence.seq)
            return fasta_sequence.seq
