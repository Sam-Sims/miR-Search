from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import Popen, PIPE, STDOUT
import re



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
    print(stdout_data)
    pattern = "\(([^)]*)\)[^(]*$" # regex to match last bracket
    result = re.findall(pattern, str(stdout_data))

    #pattern = "(\d+(?:\.\d+)?)" # regex patern to find all decimal numbers
    #result = re.findall(pattern, str(stdout_data))
    return result[0]



def run_fold(dict_of_tagets, utr, flank):
    for key, value in dict_of_tagets[0].items():
        print("Processing: ", key)
        print("Flanking region set to: ", str(flank))
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        #rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        #SeqIO.write(rec, "temp.fasta", "fasta")
        print("Free energy: " + call_rna_fold(target_seq))

