import re

class TargetSearch:
    def __init__(self, utr_input, mir_input):
        self.utr = utr_input # list of biopython seq objects
        self.mir = mir_input # search string


    def search_6mer(self):
        print("Searching for 6mer:", self.mir)
        master_list = []
        for seq in self.utr:
            print("Searching ", seq.id, "...")
            seed_locations = []
            seq_string = str(seq.seq) # convert to string for regex
            for i in re.finditer(self.mir, seq_string):
                _start = i.start()
                _end = i.end()
                _seed_locations = [_start, _end]
                seed_locations.append(seq.id)
                seed_locations.append(_seed_locations)
            if not seed_locations:
                print("False")
            else:
                print("True")
                master_list.append(seed_locations)
        print(master_list)
        print(len(master_list))
        print(self.mir)