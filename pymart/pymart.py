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

    def download_mart(self):
        df = pd.read_csv("mart_export_test.csv")
        ids = df[df.columns[0]]
        num_threads = 50
        pool = ThreadPool(num_threads)
        pool.apply_async(self.listener) # spawn listener to watch for queue activity
        pool.map(self.get_utr, ids)
        self.q.put('kill') # listener watches for kill to exit
        pool.close()
        pool.join()

    def run_check(self):
        self.create_check_file()
        self.check()

    def create_check_file(self):
        for record in SeqIO.parse("utr.fasta", "fasta"):
            with open("../completed.csv", "a") as f:
                f.write(record.id + "\n")
            f.close()

    def check(self):
        human_genes = pd.read_csv("../mart_export.csv")
        completed = pd.read_csv("../completed.csv")
        completed['gene'] = completed['gene'].str[:15]
        completed.drop(completed.columns[1], axis=1, inplace=True)
        print(completed)
        unique_completed = completed.gene.unique()
        print(unique_completed)
        print(len(unique_completed))
        print(human_genes)
        hgl = human_genes["Gene stable ID"].tolist()
        compare = list(set(hgl) - set(unique_completed))
        print(compare)
        with open('../failed.csv', 'a') as f:
            for i in compare:
                f.write(i + "\n")
            f.close()

    def listener(self):
        print("Listening...")
        with open("utr_test.fasta", 'a') as f:
            while 1:
                m = self.q.get()
                if m == 'kill':
                    f.write('killed')
                    break
                f.write(str(m) + '\n')
                print("Wrote to file")
                f.flush()

    def get_utr(self, id):
        mydoc = minidom.parse('template.xml')
        print("Read XML template")
        filter = mydoc.getElementsByTagName('Filter')
        for elem in filter:
            print("Previous: ", elem.attributes['value'].value)
            elem.attributes['value'].value = id
            print("New: ", elem.attributes['value'].value)
        string = mydoc.toxml()
        r = requests.get(self.server + string)
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


