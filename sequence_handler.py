from Bio import SeqIO
import re


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
        return self.mir_seq

    def find_6mer(self):
        print(str(self.mir_seq[1:7]))
        return str(self.mir_seq[1:7])

    def find_7mera1(self):
        return self.mir_seq[1:7] + "A"

    def find_7merm8(self):
        return self.mir_seq[:7]

    def find_8mer(self):
        return self.mir_seq[:7] + "A"


class MicroRNASearch:
    def __init__(self, input_utr):
        self.utr = input_utr

    def search_6mer(self, mir):
        print("Searching for 6mer: ", mir)
        seed_locations = []
        for i in re.finditer(mir, self.utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
        if not seed_locations:
            return False, seed_locations  # still need to return list although empty
        else:
            return True, seed_locations  # return tuple - first always bool indicating if seed site exsists, second a list of locations within the utr - if no locations list will be empty

    def search_7merm8(self, mir):  # 7mer-m8 = posistions 2-7 + match at posistion 8
        print("Searching for 7mer-m8: ", mir)
        seed_locations = []
        is_7mera1 = False
        for i in re.finditer(mir, self.utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
            _end_check_a = _end + 1
            if self.utr[_start:_end_check_a].endswith('A'):
                is_7mera1 = True
            else:
                is_7mera1 = False
        if not seed_locations:  # if list empty - no 7mers return false
            return False, seed_locations  # still need to return list although empty
        else:  # if list has items - i.e 7mers found
            if is_7mera1:  # check if A at pos1 in miRNA - if false return false as this would be 7mer-a1
                return False, seed_locations  # still need to return list although empty
            else:  # Else must be a 7mer-m8
                return True, seed_locations

    def search_7mera1(self, mir):  # probably ineffecient searching again when already searched for 7mers - LOOK IN TO
        print("Searching for 7mer-a1: ", mir)
        seed_locations = []
        is_7mera1 = False
        for i in re.finditer(mir, self.utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
        if not seed_locations:
            return False, seed_locations
        else:
            return True, seed_locations

    def search_8mers(self, mir):
        print("Searching for 8mer: ", mir)
        seed_locations = []
        for i in re.finditer(mir, self.utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
        if not seed_locations:
            return False, seed_locations  # still need to return list although empty
        else:
            return True, seed_locations


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

    def read_multi_3utr(self, path):
        utr_list = []
        for fasta_sequence in SeqIO.parse(path, "fasta"):
            print("Reading 3' UTR region... ")
            print(fasta_sequence.seq)
            utr_list.append(fasta_sequence)
        return utr_list


