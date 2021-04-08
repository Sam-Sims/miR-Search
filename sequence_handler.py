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

    def search_6mer(mir, utr):
        _mir = mir[1:7]
        print("Searching for 6mer: ", _mir)
        seed_locations = []
        for i in re.finditer(_mir, utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
        if not seed_locations:
            return False, seed_locations  # still need to return list although empty
        else:
            return True, seed_locations  # return tuple - first always bool indicating if seed site exsists, second a list of locations within the utr - if no locations list will be empty

    def search_7merm8(mir, utr):  # 7mer-m8 = posistions 2-7 + match at posistion 8
        _mir = mir[:7]  # remove pos 1
        print("Searching for 7mer-m8: ", _mir)
        seed_locations = []
        is_7mera1 = False
        for i in re.finditer(_mir, utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
            _end_check_a = _end + 1
            if utr[_start:_end_check_a].endswith('A'):
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

    def search_7mera1(mir, utr):  # probably ineffecient searching again when already searched for 7mers - LOOK IN TO
        _mir = mir[1:7] + "A"  # pos 1 always an A then look at sites 2-7
        print("Searching for 7mer-a1: ", _mir)
        seed_locations = []
        is_7mera1 = False
        for i in re.finditer(_mir, utr):
            _start = i.start()
            _end = i.end()
            _seed_locations = [_start, _end]
            seed_locations.append(_seed_locations)
        if not seed_locations:
            return False, seed_locations
        else:
            return True, seed_locations

    def search_8mers(mir, utr):
        _mir = mir[:7] + "A"  # look at sites 1-7
        print("Searching for 8mer: ", _mir)
        seed_locations = []
        for i in re.finditer(_mir, utr):
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
