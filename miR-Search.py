from Bio import SeqIO
import sequence
import re

def read_mir_testing(test):
    for fasta_sequence in SeqIO.parse(open('test-data/site-testing/' + test), 'fasta'):
        print("Reading", fasta_sequence.seq)
        return fasta_sequence.seq

def read_3utr_testing(test):
    print("TEST")
    for fasta_sequence in SeqIO.parse("test-data/site-testing/" + test, "fasta"):
        print("Reading 3' UTR region... ")
        print(fasta_sequence.seq)
        return fasta_sequence.seq


def read_mir():
    for fasta_sequence in SeqIO.parse(open('test-data/hsa-miR-451a.fasta'), 'fasta'):
        print("Reading", fasta_sequence.seq)
        return fasta_sequence.seq


def read_3utr():
    for fasta_sequence in SeqIO.parse("test-data/OSR1.fasta", "fasta"):
        print("Reading 3' UTR region... ")
        print(fasta_sequence.seq)
        return fasta_sequence.seq


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
        return False, seed_locations # still need to return list although empty
    else:
        return True, seed_locations # return tuple - first always bool indicating if seed site exsists, second a list of locations within the utr - if no locations list will be empty


def search_7merm8(mir, utr): # 7mer-m8 = posistions 2-7 + match at posistion 8
    _mir = mir[:7] # remove pos 1
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
    if not seed_locations: # if list empty - no 7mers return false
        return False, seed_locations  # still need to return list although empty
    else: # if list has items - i.e 7mers found
        if is_7mera1: # check if A at pos1 in miRNA - if false return false as this would be 7mer-a1
            return False, seed_locations  # still need to return list although empty
        else: # Else must be a 7mer-m8
            return True, seed_locations


def search_7mera1(mir, utr): # probably ineffecient searching again when already searched for 7mers - LOOK IN TO
    _mir = mir[1:7] + "A" # pos 1 always an A then look at sites 2-7
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
    _mir = mir[:7] + "A" # look at sites 1-7
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


def main():
    mir_to_search = read_mir()
    seq = sequence.microrna(mir_to_search)
    search_mir = seq.clean_mir()
    search_mir = seq.my_reverse_complement()
    search_mir = seq.back_trans() # back translate after reverse complement or A are swapped to T
    print(search_mir)
    print("Done!")
    print("Final target site search string is: ", search_mir)

    utr = read_3utr_testing("test-utr-8mer.fasta")

    sixmer = search_6mer(str(search_mir), str(utr))
    sevenmerm8 = search_7merm8(str(search_mir), str(utr))
    sevenmera1 = search_7mera1(str(search_mir), str(utr))
    eightmer = search_8mers(str(search_mir), str(utr))

    is6mer = False
    is7merm8 = False
    is7mera1 = False
    is8mer = False


    if sixmer[0]:
        print("6mers found at: ")
        print(sixmer[1])
        is6mer = True
        if sevenmerm8[0]:
            print("7mer-m8s found at: ")
            print(sevenmerm8[1])
            is7merm8 = True
        else:
            print("No 7mer-m8 sites found")
            if sevenmera1[0]:
                print("7mer-a1s found at: ")
                print(sevenmera1[1])
                is7mera1 = True
                if eightmer[0]:
                    print("8mers also found at: ")
                    print(eightmer[1])
                    is8mer = True
                else:
                    print("7mer no 8mer")
            else:
                print("No 7mer-a1 sites found")
                print("Checking for 8mers")
                if eightmer[0]:
                    print("8mers found at: ")
                    print(eightmer[1])
                    is8mer = True
                else:
                    print("No targets found")

    else:
        print("No seed sites found!")

    print("6mers: " + str(is6mer))
    print("7mer-a1: " + str(is7mera1))
    print("7mer-m8: " +str(is7merm8))
    print("8mer: " +str(is8mer))


if __name__ == "__main__":
    main()