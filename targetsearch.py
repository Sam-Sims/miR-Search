import re


class TargetSearch:
    def __init__(self, utr_input):
        self.utr = utr_input # list of biopython seq objects

    def search_6mer(self, mir):
        print("Searching for 6mer:", mir)
        result = self._search(mir)
        print(result)
        print(len(result))

    def search_7mera1(self, mir):
        print("Searching for 7mer-a1:", mir)
        result = self._search(mir)
        print(result)
        print(len(result))

    def search_7merm8(self, mir):
        print("Searching for 7mer-m8: ", mir)
        result = self._search(mir)
        print(result)
        print(len(result))

    def search_8mer(self, mir):
        print("Searching for 8mer: ", mir)
        result = self._search(mir)
        print(result)
        print(len(result))

    def _search(self, mir):
        master_list = []
        for seq in self.utr:
            #print("Searching ", seq.id, "...")
            seed_locations = []
            seq_string = str(seq.seq)  # convert to string for regex
            for i in re.finditer(mir, seq_string):
                _start = i.start()
                _end = i.end()
                _seed_locations = [_start, _end]
                seed_locations.append(seq.id)
                seed_locations.append(_seed_locations)
            if not seed_locations:
                ans = 1+1
            else:
                master_list.append(seed_locations)
        return master_list
