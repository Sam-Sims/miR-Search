import requests, sys
from xml.dom import minidom
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp
from Bio import SeqIO
import os


class PYMart:
    def __init__(self):
        manager = mp.Manager()
        self.q = manager.Queue()
        self.server = "http://www.ensembl.org/biomart/martservice?query="
        self.biomart_master_path = "mart_export_test.csv" # mart export needs to have ensembl IDs in first column.

    def download_mart(self):
        df = pd.read_csv(self.biomart_master_path)
        ids = df[df.columns[0]]
        num_threads = 50
        pool = ThreadPool(num_threads)
        pool.apply_async(self.listener) # spawn listener to watch for queue activity
        pool.map(self.get_utr, ids)
        self.q.put('kill') # listener watches for kill to exit
        pool.close()
        pool.join()

    def listener(self):
        print("Listening...")
        with open("utr_test.fasta", 'a') as f:
            while 1:
                m = self.q.get()
                if m == 'kill':
                    f.write('killed')
                    break
                f.write(str(m) + '\n')
                f.flush()

    def get_utr(self, id):
        mydoc = minidom.parse('template.xml')
        filter = mydoc.getElementsByTagName('Filter')
        for elem in filter:
            elem.attributes['value'].value = id
        string = mydoc.toxml()
        r = requests.get(self.server + string)
        print("Processing: ", id)
        print("Sending request...")

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        print("Request status code: ", r.status_code)

        '''
        file_string = "temp/temp" + id + ".fasta"
        print(file_string)
        f = open(file_string, "a")
        f.write(r.text)
        f.close()
        for seq_record in SeqIO.parse("temp.fasta", "fasta"):
            print(seq_record.seq)
        os.remove(file_string)

        '''
        self.q.put(r.text)
        print("Added 3utr to queue \n")

    def run_check(self):
        print("Running comparison check...")
        self.create_check_file() # Extracts all fasta headers from the completed file including gene ID
        self.check() # Compares them agaisnt the original list of gene IDs

    def create_check_file(self):
        print("Generating check file...")
        with open("completed.csv", "a") as f:
            f.write("gene, \n")
        f.close()
        for record in SeqIO.parse("utr_test.fasta", "fasta"):
            with open("completed.csv", "a") as f:
                f.write(record.id + "\n")
            f.close()
        print("Completed sequences stored in completed.csv")

    def check(self):
        print("Comparing...")
        human_genes = pd.read_csv(self.biomart_master_path)
        completed = pd.read_csv("completed.csv")
        completed['gene'] = completed['gene'].str[:15]
        completed.drop(completed.columns[1], axis=1, inplace=True)
        unique_completed = completed.gene.unique()
        hgl = human_genes["Gene stable ID"].tolist()
        compare = list(set(hgl) - set(unique_completed))
        with open('failed.csv', 'a') as f:
            for i in compare:
                f.write(i + "\n")
            f.close()
        print("Failed gene IDs stored in failed.csv")
