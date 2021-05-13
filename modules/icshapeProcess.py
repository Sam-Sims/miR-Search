import pandas as pd
import csv

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def return_transcripts(input):
    df = pd.read_csv(input, sep='	', header=None, names=list(range(10000)))
    print(df)
    # df.to_csv("test.csv", index=False)
    valid_transcripts = df[0].tolist()
    print(valid_transcripts)
    with open('transcipts_with_shape_data.txt', 'a') as f:
        for i in valid_transcripts:
            f.write("%s\n" % i)
    return valid_transcripts


def match_transcript_id(utr_input):
    recs = []
    seqs = {}
    with open('transcipts_with_shape_data.txt') as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    for record in SeqIO.parse(utr_input, "fasta"):
        res = any(ele in str(record.name) for ele in content)
        if res:
            seqs[record.name] = record.seq
    for name, seq in seqs.items():
        print(name)
        rec = SeqRecord(Seq(str(seq)), id=name, description="")
        recs.append(rec)

    SeqIO.write(recs, "UTR_matched_to_shapeseq.fasta", "fasta")


def count():
    i = 0
    for record in SeqIO.parse("Transcripts/final_trans.fasta", "fasta"):
        i = i+1
    print(i)


def _extract_data(header, attribute):
    # LIST OF ATTRIBUTES
    # 0: Gene ID
    # 1: Gene name
    # 2: Transcript ID
    # 3: UTR Start
    # 4: UTR End
    # 5: Transcript Start
    # 6: Transcript End
    split_header = header.split("|")
    return split_header[attribute]


def prepare_targets(target_file):
    df = pd.read_csv(target_file)
    return df


def _return_dict(df, target_type):
    location_string = target_type + "_location"
    list_to_zip = df[target_type].tolist()
    cleaned_list = [x for x in list_to_zip if str(x) != 'nan']
    ret_dict = dict(zip(cleaned_list, df[location_string]))
    return ret_dict


def gen_utr_hash_table(utr_file):
    utr_dict = {}
    for seq_record in SeqIO.parse(utr_file, "fasta"):
        utr_dict[seq_record.name] = seq_record.seq
    return utr_dict


def align(hashed_utr, key, curr_trans_length, value, row):
    utr_length = len(hashed_utr.get(key))
    utr_start = curr_trans_length - utr_length
    target_start_utr = int(value.split(",")[0][1:])
    target_start_trans = utr_start + target_start_utr
    target_end_utr_temp = value.split(",")[1][1:]
    target_end_utr = int(target_end_utr_temp[:-1])
    target_end_trans = utr_start + target_end_utr
    shape_data = row[4:]
    target_shape_score = shape_data[target_start_trans:target_end_trans]
    return target_shape_score


def process_target_shape_data(input, target_df):
    shape_6mer = {}
    dict_6mer = _return_dict(target_df, "6mer")

    shape_7mera1 = {}
    dict_7mera1 = _return_dict(target_df, "7mera1")

    shape_7merm8 = {}
    dict_7merm8 = _return_dict(target_df, "7merm8")

    shape_8mer = {}
    dict_8mer = _return_dict(target_df, "8mer")

    hashed_utr = gen_utr_hash_table("Transcripts/final_trans.fasta") # hash utr file for easy look up - dont have to iterate it each time
    shape_score_dict_6mer = {}
    shape_score_dict_7mera1 = {}
    shape_score_dict_7merm8 = {}
    shape_score_dict_8mer = {}

    with open(input) as f:
        read = csv.reader(f, delimiter="\t")
        for row in read:
            curr_trans = row[0]
            curr_trans_length = int(row[1])
            for key, value in dict_6mer.items(): # maybe rework header info here to allow for a direct look up in dict (dict.get())
                if curr_trans == _extract_data(key, 2):
                    print("6mer target data found for " + curr_trans)
                    target_shape_score = align(hashed_utr, key, curr_trans_length, value, row)
                    shape_score_dict_6mer[curr_trans] = target_shape_score
            for key, value in dict_7mera1.items():
                if curr_trans == _extract_data(key, 2):
                    print("7mera1 target data found for " + curr_trans)
                    target_shape_score = align(hashed_utr, key, curr_trans_length, value, row)
                    #key_for_dict = curr_trans + " "
                    shape_score_dict_7mera1[curr_trans] = target_shape_score
            for key, value in dict_7merm8.items():
                if curr_trans == _extract_data(key, 2):
                    print("7merm8 target data found for " + curr_trans)
                    target_shape_score = align(hashed_utr, key, curr_trans_length, value, row)
                    #key_for_dict = curr_trans + " "
                    shape_score_dict_7merm8[curr_trans] = target_shape_score
            for key, value in dict_8mer.items():
                if curr_trans == _extract_data(key, 2):
                    print("8mer target data found for " + curr_trans)
                    target_shape_score = align(hashed_utr, key, curr_trans_length, value, row)
                    #key_for_dict = curr_trans + " "
                    shape_score_dict_8mer[curr_trans] = target_shape_score
    f.close()
    return shape_score_dict_6mer, shape_score_dict_7mera1, shape_score_dict_7merm8, shape_score_dict_8mer
