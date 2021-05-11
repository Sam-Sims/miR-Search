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


def match_transcript_id():
    recs = []
    seqs = {}
    with open('transcipts_with_shape_data.txt') as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    for record in SeqIO.parse("Transcripts/UTR_w_transcripts_removed_non_seq.fasta", "fasta"):
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
    for record in SeqIO.parse("final.fasta", "fasta"):
        i = i+1
    print(i)

def align_scores(input):
    with open(input) as fd:
        rd = csv.reader(fd, delimiter="\t")
        for row in rd:
            curr_trans = row[0]
            print(len(row))
            for record in SeqIO.parse("Transcripts/UTR_matched_SHAPE_UNIQUE.fasta", "fasta"):
                if curr_trans in record.name:
                    print("Record in")
                    curr_seq = str(record.seq)
                    print(len(curr_seq))
