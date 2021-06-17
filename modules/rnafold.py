from subprocess import Popen, PIPE, STDOUT
import re
from operator import itemgetter


def calc_location(value, flank):  # messy parsing of the location as it is read as a string not list
    target_start_utr = int(value.split(",")[0][1:])
    target_start_utr_flank = target_start_utr - int(flank)  # flank region downstream
    if target_start_utr_flank < 0:
        target_start_utr_flank = 0
    target_end_utr_temp = value.split(",")[1][1:]
    target_end_utr = int(target_end_utr_temp[:-1])
    target_end_utr_flank = target_end_utr + int(flank)  # add flank region upstream
    return target_start_utr_flank, target_end_utr_flank, target_start_utr, target_end_utr


def call_rna_cofold(stdin):
    str_stdin = str(stdin)
    p = Popen(['RNAcofold'], stdout=PIPE, stdin=PIPE, stderr=PIPE)  # run RNA cofold and pipe stdout
    stdout_data = p.communicate(input=str_stdin.encode())[0]
    pattern = "\(([^)]*)\)[^(]*$"  # regex to match the number in the last bracket of a string
    result = re.findall(pattern, str(stdout_data))
    return result[0]


def call_rna_fold(stdin):
    str_stdin = str(stdin)
    p = Popen(['RNAfold', '-d2', '--noLP'], stdout=PIPE, stdin=PIPE, stderr=PIPE)  # run RNA fold and pipe stdout
    stdout_data = p.communicate(input=str_stdin.encode())[0]
    #print(stdout_data)

    # Regex to extract dot bracket notation from stdout
    dot_bracket_pattern = "n(.*)"
    dot_bracket = re.search(dot_bracket_pattern, str(stdout_data))
    dot_bracket_pattern = "^\S*"
    dot_bracket = re.search(dot_bracket_pattern, dot_bracket.group())
    final_dot_bracket = dot_bracket.group()[1:]

    pattern = "\(([^)]*)\)[^(]*$"  # regex to match the number in the last bracket of a string
    result = re.findall(pattern, str(stdout_data))

    # pattern = "(\d+(?:\.\d+)?)" # regex patern to find all decimal numbers
    # result = re.findall(pattern, str(stdout_data))
    return result[0], final_dot_bracket


def return_percentage(dict_to_process):
    combined_dict = {}
    for key, value in dict_to_process.items():
        temp_dict = {}
        print("Sorting: ", key)
        sorted_dict = dict(sorted(value.items(), key=lambda item: item[1]))  # sort dict low to high
        n = int(len(sorted_dict) * 0.20)  # calc 20% length of dict
        bot_20 = dict(sorted(sorted_dict.items(), key=itemgetter(1))[:n])
        top_20 = dict(sorted(sorted_dict.items(), key=itemgetter(1), reverse=True)[:n])
        temp_dict["top"] = top_20
        temp_dict["bot"] = bot_20
        combined_dict[key] = temp_dict
    return combined_dict


def return_percentage_subset(dict_to_process):
    sorted_dict = dict(sorted(dict_to_process.items(), key=lambda item: item[1]))  # sort dict low to high
    n = int(len(sorted_dict) * 0.20)  # calc 20% length of dict
    bot_20 = dict(sorted(sorted_dict.items(), key=itemgetter(1))[:n])
    top_20 = dict(sorted(sorted_dict.items(), key=itemgetter(1), reverse=True)[:n])
    combined_dict = {"top": top_20, "bot": bot_20}
    return combined_dict


def get_nested_keys(d, keys):  # recursive function to retrieve all keys in a nested dictionary
    for key, value in d.items():
        if isinstance(value, dict):
            get_nested_keys(value, keys)
        else:
            keys.append(key)
    keys_list = []
    return keys_list


def return_shape_transcripts(shape_data):
    keys_list = []
    get_nested_keys(shape_data, keys_list)
    return keys_list


def subset_target_data(target_data, shape_transcripts):
    flat_dict = {}
    master_dict = {}
    for i in target_data:
        for key, value in i.items():
            flat_dict[key] = value

    for key, value in flat_dict.items():
        for i in shape_transcripts:
            if i in key:
                master_dict[key] = value
    return master_dict


def run_fold_subset(dict_of_targets, utr, flank, cofold):
    energy_score_dict = {}
    energy_score_cofold = {}
    print("Flanking region set to: ", str(flank))
    counter = 0
    for key, value in dict_of_targets.items():
        print("Processing: ", key)
        split_header = key.split("|")
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        seed = utr.get(key)[utr_locations[2]:utr_locations[3]]
        #print("SEED: ", seed)
        seed = str(seed)
        # rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        # SeqIO.write(rec, "temp.fasta", "fasta")

        free_energy = call_rna_fold(target_seq)
        # print("Free energy: " + free_energy)
        if cofold:
            cofold_string = prep_cofold(free_energy[1], target_seq, seed)
            print("COFOLD STRING", cofold_string)
            cofold_energy = call_rna_cofold(cofold_string)
            energy_score_cofold[split_header[2]] = cofold_energy
        else:
            energy_score_dict[split_header[2]] = free_energy[0]
    #print("COUNTER: ", counter)
    if cofold:
        return energy_score_cofold
    else:
        return energy_score_dict


def process_dot_bracket(dot_bracket):
    stack = []
    pairs = {}
    for count, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(count)
        elif char == ')':
            popped = stack.pop()
            pairs[popped] = count
    print(pairs)
    return pairs

def match_sequence(pairing_table, sequence):
    #print(sequence)
    paired_dict = {}
    for count, char in enumerate(sequence):
        if count in pairing_table.keys():
            matched_pair = pairing_table.get(count)
            key = str(count) + char
            paired_dict[key] = str(matched_pair) + sequence[matched_pair]
        else:
            key = char + str(count)
            paired_dict[key] = "N"
    return paired_dict

def find_seed(seed, target_seq):
    for i in re.finditer(seed, target_seq):
        start = int(i.start())
        end = int(i.end())
        seed_locations = [start, end]
    return seed_locations


def match_unpaired(match_list, dot_bracket, length):
    '''
    index_start = 0
    forward = False
    if "(" in dot_bracket:
        forward = True
    if forward:
        for count, char in enumerate(match_list):
            if char != "N":
                index_start = count
    '''
    #print("DOT", dot_bracket)
    #print(match_list)
    n = len(match_list)
    if all(elem == match_list[0] for elem in match_list): # if all elements match aka all bases have N aka all bases unpaired
        print("NO PAIRS") # do nothing
    else: # do a kind of bubble sort comparing each value to the next until solved
        if "(" in dot_bracket: # IF seed site downstream - ( indicates that  the closing ) will be further along rather than closer
            while "N" in match_list: # iterate over the list until solved
                for j in range (n-1):
                    if isinstance(match_list[j], int) and match_list[j + 1] == "N":
                        if match_list[j] - 1 in match_list:
                            #print("Bulge?")
                            match_list[j + 1] = match_list[j] - 1
                            match_list[match_list.index(match_list[j] - 1)] = 'N'
                        else:
                            if match_list[j] - 1 < 0:
                                match_list[j + 1] = 0
                            else:
                                match_list[j + 1] = match_list[j] - 1
                    elif match_list[j] == "N" and isinstance(match_list[j + 1], int):
                        if match_list[j + 1] + 1 in match_list:
                            #print("bulge?")
                            match_list[j] = match_list[j + 1] + 1
                            match_list[match_list.index(match_list[j + 1] + 1)] = 'N'
                        else:
                            if match_list[j + 1] + 1 > length - 1:
                                match_list[j] = length - 1
                            else:
                                match_list[j] = match_list[j + 1] + 1
        elif ")" in dot_bracket:
            match_list.reverse()
            #print("REV: ", match_list)
            while "N" in match_list:
                for j in range (n-1):
                    if isinstance(match_list[j], int) and match_list[j + 1] == "N":
                        if match_list[j] + 1 in match_list:
                            #print("Bulge?")
                            match_list[j + 1] = match_list[j] + 1
                            match_list[match_list.index(match_list[j] + 1)] = 'N'
                        else:
                            if match_list[j] + 1 > length - 1:
                                match_list[j + 1] = length - 1
                            else:
                                match_list[j + 1] = match_list[j] + 1
                    elif match_list[j] == "N" and isinstance(match_list[j + 1], int):
                        if match_list[j + 1] - 1 in match_list:
                            #print("bulge?")
                            match_list[j] = match_list[j + 1] - 1
                            match_list[match_list.index(match_list[j + 1] - 1)] = 'N'
                        else:
                            if match_list[j + 1] - 1 < 0:
                                match_list[j] = 0
                            else:
                                match_list[j] = match_list[j + 1] - 1

    #print("FIXED: ", match_list)
    return match_list



        


def prep_cofold(dot_bracket, target_seq, seed):
    #print(target_seq)
    #print(dot_bracket)
    target_seq = str(target_seq)
    seed_loc = find_seed(seed, target_seq)
    seed_dot_bracket = dot_bracket[seed_loc[0]:seed_loc[1]]
    #print(seed_dot_bracket)
    pairing_table = process_dot_bracket(dot_bracket)
    seq_match = match_sequence(pairing_table, target_seq)
    inverse_pairing_table = {v: k for k, v in pairing_table.items()}
    seed_span = list(range(seed_loc[0], seed_loc[1]))
    if "(" in seed_dot_bracket and ")" in seed_dot_bracket and "." in seed_dot_bracket:
        return 0
    else:
        match_list =[]
        for count, char in enumerate(seed_dot_bracket):
            if char == "(":
                global_pos = seed_span[count]  # posistion in whole dot bracket struct
                if global_pos in pairing_table.keys():
                    match = pairing_table.get(global_pos)
                    match_list.append(match)
            elif char == ")":
                global_pos = seed_span[count]  # posistion in whole dot bracket struct
                if global_pos in inverse_pairing_table.keys():
                    match = inverse_pairing_table.get(global_pos)  # use the inverse of the pairing table as looking downstream
                    match_list.append(match)
            elif char == ".":
                match_list.append("N")
        paired_list = match_unpaired(match_list, seed_dot_bracket, len(target_seq))
        #print("SEED", seed)
        input1 = seed
        paired_seq =[]
        for i in match_list:
            if i == "N":
                break
            else:
                temp = target_seq[i]
                paired_seq.append(temp)
        input2 = ''.join(paired_seq)
        final_string = input1 + "&" + input2
        print(final_string)
        return final_string



def run_fold(dict_of_tagets, utr, flank):
    energy_score_dict_6mer = {}
    energy_score_dict_7mera1 = {}
    energy_score_dict_7merm8 = {}
    energy_score_dict_8mer = {}
    print("Flanking region set to: ", str(flank))
    for key, value in dict_of_tagets[0].items():
        split_header = key.split("|")
        # print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        # rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        # SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        # print("Free energy: " + free_energy)
        energy_score_dict_6mer[split_header[2]] = free_energy[0]
    for key, value in dict_of_tagets[1].items():
        split_header = key.split("|")
        # print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        # rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        # SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        # print("Free energy: " + free_energy)
        energy_score_dict_7mera1[split_header[2]] = free_energy[0]
    for key, value in dict_of_tagets[2].items():
        split_header = key.split("|")
        # print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        # rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        # SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        # print("Free energy: " + free_energy)
        energy_score_dict_7merm8[split_header[2]] = free_energy[0]
    for key, value in dict_of_tagets[3].items():
        split_header = key.split("|")
        # print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        # rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        # SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        # print("Free energy: " + free_energy)
        energy_score_dict_8mer[split_header[2]] = free_energy[0]
    combined_dict = {"6mer": energy_score_dict_6mer, "7mera1": energy_score_dict_7mera1,
                     "7merm8": energy_score_dict_7merm8,
                     "8mer": energy_score_dict_8mer}  # creates a dict of dict first
    return combined_dict
