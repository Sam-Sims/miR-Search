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
        print("Final search sequence: ", self.mir_seq)

    def return_seq(self):
        return str(self.mir_seq)

    def find_6mer(self):
        return str(self.mir_seq[1:7])

    def find_7mera1(self):
        return str(self.mir_seq[1:7] + "A")

    def find_7merm8(self):
        return str(self.mir_seq[:7])

    def find_8mer(self):
        return str(self.mir_seq[:7] + "A")


class FASTAParse:
    def __init__(self, v_print):
        self.v_print = v_print

    def read_mir(self, path):
        for fasta_sequence in SeqIO.parse(open(path), 'fasta'):
            print("Reading miRNA at " + path + "...")
            print("miRNA: ", fasta_sequence.seq)
            return fasta_sequence.seq

    def read_multi_3utr(self, path):
        utr_list = []
        for fasta_sequence in SeqIO.parse(path, "fasta"):
            self.v_print(1, "Reading 3' UTR region... ")
            self.v_print(1, fasta_sequence.seq)
            utr_list.append(fasta_sequence)
        return utr_list

    def read_mirbase(self, input_mir):
        print("Reading database...")
        for fasta_sequence in SeqIO.parse(open("database/mirs.fasta"), 'fasta'):
            if input_mir in fasta_sequence.id:
                print(fasta_sequence.id)
                return fasta_sequence.seq
        return "Error"
