import multiprocessing as mp
import os
import re
from multiprocessing.dummy import Pool as ThreadPool
from xml.dom import minidom

import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class PYMart:
    def __init__(self, server, genelist, output_file, split, v_print):
        manager = mp.Manager()
        self._q = manager.Queue()
        self._server = server
        self.gene_list = genelist  # mart export needs to have ensembl IDs in first column.
        self._completed_path = "utr_check/completed.csv"
        self.output_file = output_file
        self.split = split
        self.vprint = v_print

    def set_split(self):
        self.split = True

    def download_mart(self, thread_num):
        df = pd.read_csv(self.gene_list)
        ids = df[df.columns[0]]
        num_threads = thread_num
        pool = ThreadPool(num_threads)
        pool.apply_async(self._listener)  # spawn listener to watch for queue activity
        pool.map(self._worker, ids)
        self._q.put('kill')  # listener watches for kill to exit
        pool.close()
        pool.join()

    def _listener(self):
        print("Listening...")
        with open(self.output_file, 'a') as f:
            while 1:
                m = self._q.get()
                if m == 'kill':
                    f.write('killed')
                    break
                f.write(str(m) + '\n')
                f.flush()

    def _worker(self, id):
        mydoc = minidom.parse('templates/template_utr.xml')
        filter = mydoc.getElementsByTagName('Filter')
        for elem in filter:
            elem.attributes['value'].value = id
            break
        string = mydoc.toxml()
        r = requests.get(self._server + string)
        print("Processing: ", id)
        print("Sending request...")

        if not r.ok:
            with open("error.txt", "a") as f:
                f.write(id + "\n")
            f.close()
        else:
            print("Request status code: ", r.status_code)

            if not self.split:  # if not spliting output add the utr to a queue so listener can write to file
                self._q.put(r.text)
                print("Added 3utr to queue \n")
            else:  # split the output and have each thread write to a seperate file
                _id = str(id) + ".fasta"
                with open(_id, 'a') as f:
                    f.write(str(r.text) + '\n')
                f.close()

    def run_verify(self, input):
        print("Running comparison check...")
        if not os.path.exists('utr_check'):
            os.makedirs('utr_check')
        self._create_check_file(input)  # Extracts all fasta headers from the completed file including gene ID
        self._check()  # Compares them agaisnt the original list of gene IDs
        self._check_completion_stamp(
            input)  # biomart returns own completion stamp with the response - check for any failures

    def run_check(self):  # funct for auto mode
        print("Running comparison check...")
        if not os.path.exists('utr_check'):
            os.makedirs('utr_check')
        self._create_check_file(
            self.output_file)  # Extracts all fasta headers from the completed file including gene ID
        self._check()  # Compares them agaisnt the original list of gene IDs
        self._check_completion_stamp(
            self.output_file)  # biomart returns own completion stamp with the response - check for any failures

    def _create_check_file(self, input):
        print("Generating check file...")
        with open(self._completed_path, "a") as f:
            f.write("gene, \n")
        f.close()
        for record in SeqIO.parse(input, "fasta"):
            with open(self._completed_path, "a") as f:
                f.write(record.id + "\n")
            f.close()
        print("Completed sequences stored in utr_check/completed.csv")

    def _check(self):
        print("Comparing...")
        human_genes = pd.read_csv('master_gene_list.csv')
        completed = pd.read_csv(self._completed_path)
        completed['gene'] = completed['gene'].str[:15]
        completed.drop(completed.columns[1], axis=1, inplace=True)
        unique_completed = completed.gene.unique()
        hgl = human_genes["Gene stable ID"].tolist()
        compare = list(set(hgl) - set(unique_completed))
        with open('utr_check/failed.csv', 'a') as f:
            for i in compare:
                f.write(i + "\n")
            f.close()
        print("Failed gene IDs stored in utr_check/failed.csv")

    def _check_completion_stamp(self, input):
        for seq_record in SeqIO.parse(input, "fasta"):
            # print(seq_record.seq)
            if "failed" in str(seq_record.seq):
                print(seq_record.name, " may of failed")
                with open('utr_check/failed_completion_stamp.csv', 'a') as f:
                    f.write(seq_record.name + "\n")
                f.close()

    def clean_utr(self, file):  # This entire function is probably bad - look at the file writing
        unique = False  # removes all duplicate IDs keeping those with longest seq only - shouldnt be true if unique_header true
        remove_dup = False  # removes duplicates by comparing sequences
        unique_header = True  # IF true will remove dup and keep longest if transcript ID present
        rm_list = []

        print("Cleaning UTR file...")
        print("Removing failed records...")
        for seq_record in SeqIO.parse(file, "fasta"):
            # print(seq_record.seq)
            if "failed" in str(seq_record.seq):
                print(seq_record.name, " may of failed - removing")
                rm_list.append(seq_record.name)
            else:
                with open("temp.fasta", "a") as f:
                    SeqIO.write(seq_record, f, "fasta")
        print("Done!")
        print("Removing records with no 3'UTR sequence...")
        for seq_record in SeqIO.parse("temp.fasta", "fasta"):
            if seq_record.seq == "Sequenceunavailable[success]" or seq_record.seq == "Sequenceunavailable" or seq_record.seq == "Sequenceunavailablekilled":
                print("No 3'UTR for", seq_record.name, " - removing")
                rm_list.append(str(seq_record.name))
            else:
                print("3'UTR sequence found for ", seq_record.name)
                with open("temp_2.fasta", "a") as f:
                    SeqIO.write(seq_record, f, "fasta")
        print("Done!")
        print("Removing all success markers...")
        f = open("temp_2.fasta", "r")
        file_string = str(f.read())
        f.close()
        # pat = "\[(?:[^\[\]]++|(?0))*+]"
        temp_new_file_string = re.sub("\[(?:[^\[\]]+|)]", "",
                                      file_string)  # regex to find any [success] markers and remove them
        print("Done!")
        print("Removing non sequence strings...")
        new_file_string = temp_new_file_string.replace("killed", "")
        print("Done!")
        print("Removing duplicate sequences...")

        with open("temp3.fasta", "a") as f:
            f.write(new_file_string)
        f.close()

        if remove_dup:
            # REMOVES DUPLICATES - takes long time
            # Python sets might be better here - to look into
            print("Removing duplicate sequences")
            seen = []
            records = []
            for record in SeqIO.parse("temp3.fasta", "fasta"):
                if str(record.seq) not in seen:
                    self.vprint(1, record.id + " not seen")
                    seen.append(str(record.seq))
                    records.append(record)
                else:
                    self.vprint(1, "Record seen: ", record)
            SeqIO.write(records, "temp4.fasta", "fasta")

            self.vprint(1, "Writing cleaned fasta file as temp4.fasta")
            print("Done!")

        if unique:
            print("Generating unique UTR file")  # this removes all duplicate gene headers leaving the longest sequence
            seqs = {}
            recs = []
            for seq_record in SeqIO.parse("temp4.fasta", "fasta"):
                if seq_record.name not in seqs:
                    self.vprint(1, "Unique found: " + seq_record.name)
                    seqs[seq_record.name] = seq_record.seq
                else:
                    if len(seqs[seq_record.name]) <= len(seq_record.seq):
                        self.vprint(1, "Longer record found")
                        seqs[seq_record.name] = seq_record.seq
                    print("Duplicate found, ", seq_record.name)
            for name, seq in seqs.items():
                rec = SeqRecord(Seq(str(seq)), id=name, description="")
                recs.append(rec)
            SeqIO.write(recs, "final_no_trans.fasta", "fasta")
            print("Done!")

        if unique_header:  # same as above function but if headers are not unique (i.e containing transcript ids)
            print("Generating unique UTR file (unique header mode)")
            seqs = {}
            recs = []
            transcript = {}
            for seq_record in SeqIO.parse("temp4.fasta", "fasta"):
                if seq_record.name[0:15] not in seqs:
                    self.vprint(1, "Unique found: " + seq_record.name[0:15])
                    seqs[seq_record.name[0:15]] = seq_record.seq
                else:
                    if len(seqs[seq_record.name[0:15]]) <= len(seq_record.seq):
                        self.vprint(1, "longer record found")
                        seqs[seq_record.name[0:15]] = seq_record.seq
                    print("Duplicate found, ", seq_record.name[0:15])
                transcript[seq_record.name] = seq_record.seq
            for name, seq in seqs.items():
                print(name)
                rec = SeqRecord(Seq(str(seq)), id=name, description="")
                recs.append(rec)

            SeqIO.write(recs, "temp5.fasta", "fasta")
            new_seqs = {}
            new_recs = []
            for seq in SeqIO.parse("temp5.fasta", "fasta"):
                for key, value in transcript.items():
                    if seq.seq == str(value):
                        new_seqs[key] = seq.seq
            for name1, seq1 in new_seqs.items():
                print(name1)
                new_rec = SeqRecord(Seq(str(seq1)), id=name1, description="")
                new_recs.append(new_rec)
            SeqIO.write(new_recs, "final_trans.fasta", "fasta")

        print("Removing temp files...")
        os.remove("temp.fasta")
        os.remove("temp_2.fasta")
        os.remove("temp3.fasta")
        os.remove("temp4.fasta")
        print("FASTA file cleaned!")
        with open("removed_during_cleaning.txt", "a") as f:
            for i in rm_list:
                f.write("%s\n" % i)
        f.close()
