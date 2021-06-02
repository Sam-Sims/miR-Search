
from subprocess import Popen, PIPE, STDOUT
import re
from operator import itemgetter



def calc_location(value, flank): # messy parsing of the location as it is read as a string not list
    target_start_utr = int(value.split(",")[0][1:])
    target_start_utr = target_start_utr - int(flank)  # flank region downstream
    if target_start_utr < 0:
        target_start_utr = 0
    target_end_utr_temp = value.split(",")[1][1:]
    target_end_utr = int(target_end_utr_temp[:-1])
    target_end_utr = target_end_utr + int(flank)  # add flank region upstream
    return target_start_utr, target_end_utr

def call_rna_fold(stdin):
    str_stdin = str(stdin)
    p = Popen(['RNAfold', '-d2' , '--noLP'], stdout=PIPE, stdin=PIPE, stderr=PIPE) # run RNA fold and pipe stdout
    stdout_data = p.communicate(input=str_stdin.encode())[0]
    #print(stdout_data)
    pattern = "\(([^)]*)\)[^(]*$" # regex to match the number in the last bracket of a string
    result = re.findall(pattern, str(stdout_data))

    #pattern = "(\d+(?:\.\d+)?)" # regex patern to find all decimal numbers
    #result = re.findall(pattern, str(stdout_data))
    return result[0]

def return_percentage(dict_to_process):
    combined_dict = {}
    for key, value in dict_to_process.items():
        temp_dict = {}
        print("Sorting: ", key)
        sorted_dict = dict(sorted(value.items(), key=lambda item: item[1])) # sort dict low to high
        n = int(len(sorted_dict) * 0.20) # calc 20% length of dict
        bot_20 = dict(sorted(sorted_dict.items(), key=itemgetter(1))[:n])
        top_20 = dict(sorted(sorted_dict.items(), key=itemgetter(1), reverse=True)[:n])
        temp_dict["top"] = top_20
        temp_dict["bot"] = bot_20
        combined_dict[key] = temp_dict
    return combined_dict



def run_fold(dict_of_tagets, utr, flank):
    energy_score_dict_6mer = {}
    energy_score_dict_7mera1 = {}
    energy_score_dict_7merm8 = {}
    energy_score_dict_8mer = {}
    print("Flanking region set to: ", str(flank))
    for key, value in dict_of_tagets[0].items():
        split_header = key.split("|")
        #print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        #rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        #SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        #print("Free energy: " + free_energy)
        energy_score_dict_6mer[split_header[2]] = free_energy
    for key, value in dict_of_tagets[1].items():
        split_header = key.split("|")
        #print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        #rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        #SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        #print("Free energy: " + free_energy)
        energy_score_dict_7mera1[split_header[2]] = free_energy
    for key, value in dict_of_tagets[2].items():
        split_header = key.split("|")
        #print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        #rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        #SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        #print("Free energy: " + free_energy)
        energy_score_dict_7merm8[split_header[2]] = free_energy
    for key, value in dict_of_tagets[3].items():
        split_header = key.split("|")
        #print("Processing: ", key)
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        #rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        #SeqIO.write(rec, "temp.fasta", "fasta")
        free_energy = call_rna_fold(target_seq)
        #print("Free energy: " + free_energy)
        energy_score_dict_8mer[split_header[2]] = free_energy
    combined_dict = {"6mer": energy_score_dict_6mer, "7mera1": energy_score_dict_7mera1,
                     "7merm8": energy_score_dict_7merm8,
                     "8mer": energy_score_dict_8mer}  # creates a dict of dict first
    return combined_dict


