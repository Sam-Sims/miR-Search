class Microrna:
    def __init__(self, seq):
        self.seq = seq

    def my_reverse_complement(self):
        print("Finding reverse complement...")
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        calculated_complement = "".join(complement.get(base, base) for base in self.seq)
        reverse_complement = calculated_complement[::-1]
        self.seq = reverse_complement
        print(reverse_complement)
        print("Done!")
        return self.seq

    def clean_mir(self):
        print("Cleaned microRNA - ", self.seq[:8])
        self.seq = self.seq[:8]
        return self.seq # Returns the first 1-8 posistions - allows us to look for 8mers then progress down into 7, 6

    def back_trans(self):
        print("Back translating...")
        complement = {'U': 'A'}
        calculated_complement = "".join(complement.get(base, base) for base in self.seq)
        self.seq = calculated_complement
        return self.seq
