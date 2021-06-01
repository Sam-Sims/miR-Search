from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import Popen, PIPE, STDOUT



def calc_location(value, flank): # messy parsing of the location as it is read as a string not list
    target_start_utr = int(value.split(",")[0][1:])
    target_start_utr = target_start_utr - int(flank)  # flank region downstream
    target_end_utr_temp = value.split(",")[1][1:]
    target_end_utr = int(target_end_utr_temp[:-1])
    target_end_utr = target_end_utr + int(flank)  # add flank region upstream
    return target_start_utr, target_end_utr

def call_rna_fold(stdin):
    print(stdin)
    str_stdin = str(stdin)
    p = Popen(['RNAfold', '-p', '-d2' , '--noLP'], stdout=PIPE, stdin=PIPE, stderr=PIPE) # run RNA fold and pipe stdout
    stdout_data = p.communicate(input=str_stdin.encode())[0]
    print(stdout_data)



def run_fold(dict_of_tagets, utr, flank):
    for key, value in dict_of_tagets[0].items():
        print("Processing: ", key)
        print("Flanking region set to: ", str(flank))
        utr_locations = calc_location(value, flank)
        target_seq = utr.get(key)[utr_locations[0]:utr_locations[1]]
        print(target_seq)
        rec = SeqRecord(Seq(str(target_seq)), id=key, description="")
        SeqIO.write(rec, "temp.fasta", "fasta")
        call_rna_fold(target_seq)

