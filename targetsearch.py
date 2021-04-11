import re
import pickle, sys


class TargetSearch:
    def __init__(self, utr_input):
        self.utr = utr_input  # list of biopython seq objects

    def search_6mer(self, mir):
        print("Searching for 6mer:", mir)
        result = self._search(mir)
        # print(result)
        print(len(result))
        return result

    def search_7mera1(self, mir):
        print("Searching for 7mer-a1:", mir)
        result = self._search(mir)
        # print(result)
        print(len(result))
        return result

    def search_7merm8(self, mir):
        print("Searching for 7mer-m8: ", mir)
        result = self._search(mir)
        # print(result)
        print(len(result))
        return result

    def search_8mer(self, mir):
        print("Searching for 8mer: ", mir)
        result = self._search(mir)
        # print(result)
        print(len(result))
        return result

    def _search(self, mir):
        master_list = []
        for seq in self.utr:
            # print("Searching ", seq.id, "...")
            seed_locations = []
            seq_string = str(seq.seq)  # convert to string for regex
            for i in re.finditer(mir, seq_string):
                _start = i.start()
                _end = i.end()
                _seed_locations = [_start, _end]
                seed_locations.append(seq.id)
                seed_locations.append(_seed_locations)
            if not seed_locations:
                ans = 1 + 1  # just to get if for now so console doesnt fill up with print
            else:
                master_list.append(seed_locations)
        return master_list

    def generate_gene_value_dict(self, six, sevena1, sevenm8,
                                 eight):  # produces dict of genes as keys and list of t/f if they are 6,7,8mer as values
        gene_id_dict = {}
        print(six)
        for seq in self.utr:  # extract every gene and create a dict - I think I do this somewhere else look at code repition?
            # this is probobaly horrible but works
            targets = []
            res6 = any(seq.id in sublist for sublist in six)  # this is also time consuming - optimise?
            if res6:
                print(seq.id, " in 6mer list")
                targets.append(res6)
            else:
                print(seq.id, " not in 6mer list")
                targets.append(res6)
            res7 = any(seq.id in sublist for sublist in sevena1)  # this is also time consuming - optimise?
            if res7:
                print(seq.id, " in 7a1mer list")
                targets.append(res7)
            else:
                print(seq.id, " not in 7a1mer list")
                targets.append(res7)
            res7m = any(seq.id in sublist for sublist in sevenm8)  # this is also time consuming - optimise?
            if res7m:
                print(seq.id, " in 7m8mer list")
                targets.append(res7m)
            else:
                print(seq.id, " not in 7m8mer list")
                targets.append(res7m)
            res8 = any(seq.id in sublist for sublist in eight)  # this is also time consuming - optimise?
            if res8:
                print(seq.id, " in 8mer list")
                targets.append(res8)
            else:
                print(seq.id, " not in 8mer list")
                targets.append(res8)
            print(targets)
            gene_id_dict[str(seq.id)] = targets
        #print(gene_id_dict)  # dict containg gene id then boolean value for 6mer 7mera1, 7merm8 8mer in that order
        with open('gene_id_dict.pickle', 'wb') as handle:
            pickle.dump(gene_id_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return gene_id_dict

    def calc_gene_targets(self):
        with open('gene_id_dict.pickle', 'rb') as handle:
            gene_dict = pickle.load(handle)
        '''
        sys.stdout = open("gene_dict.txt", "a")
        print(gene_dict)
        sys.stdout.close()
        '''
        s = []
        se1 = []
        se8 = []
        e = []
        dup = []
        for key, value in gene_dict.items():
            #print(key)
            if value[3]: # if 8mer true (everything else will true since subsets)
                #print("8mer: ", key, value)
                e.append(key) # store target as 8mer
            elif value[0] and not value[1] and not value[2] and not value[3]:  # if 6mer and nothing else
                #print("6mer: ", key, value)
                s.append(key)
            elif value[1]:
                if value[1] and value[2]:
                    print("Duplicates: ", key, value)
                    dup.append(key)
                else:
                    se1.append(key)
            elif value[2]:
                se8.append(key)
        print(s)
        print(se1)
        print(se8)
        print(e)
        print(dup)



