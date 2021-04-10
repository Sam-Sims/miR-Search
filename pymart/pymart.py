import requests, sys, os, re
from xml.dom import minidom
import pandas as pd
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp
from Bio import SeqIO


class PYMart:
    def __init__(self, server, genelist, output_file):
        manager = mp.Manager()
        self._q = manager.Queue()
        self._server = server
        self.gene_list = genelist # mart export needs to have ensembl IDs in first column.
        self._completed_path = "utr_check/completed.csv"
        self.output_file = output_file

    def download_mart(self, thread_num):
        df = pd.read_csv(self.gene_list)
        ids = df[df.columns[0]]
        num_threads = thread_num
        pool = ThreadPool(num_threads)
        pool.apply_async(self._listener) # spawn listener to watch for queue activity
        pool.map(self._get_utr, ids)
        self._q.put('kill') # listener watches for kill to exit
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

    def _get_utr(self, id):
        mydoc = minidom.parse('template.xml')
        filter = mydoc.getElementsByTagName('Filter')
        for elem in filter:
            elem.attributes['value'].value = id
        string = mydoc.toxml()
        r = requests.get(self._server + string)
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
        self._q.put(r.text)
        print("Added 3utr to queue \n")

    def run_check(self):
        print("Running comparison check...")
        if not os.path.exists('utr_check'):
            os.makedirs('utr_check')
        self._create_check_file() # Extracts all fasta headers from the completed file including gene ID
        self._check() # Compares them agaisnt the original list of gene IDs
        self._check_completion_stamp()

    def _create_check_file(self):
        print("Generating check file...")
        with open(self._completed_path, "a") as f:
            f.write("gene, \n")
        f.close()
        for record in SeqIO.parse(self.output_file, "fasta"):
            with open(self._completed_path, "a") as f:
                f.write(record.id + "\n")
            f.close()
        print("Completed sequences stored in utr_check/completed.csv")

    def _check(self):
        print("Comparing...")
        human_genes = pd.read_csv(self.gene_list)
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

    def _check_completion_stamp(self):
        for seq_record in SeqIO.parse(self.output_file, "fasta"):
            #print(seq_record.seq)
            if "failed" in str(seq_record.seq):
                print(seq_record.name, " may of failed")
                with open('utr_check/failed_completion_stamp.csv', 'a') as f:
                    f.write(seq_record.name + "\n")
                f.close()

    def clean_utr(self, file):
        rm_list = []
        print("Cleaning UTR file...")
        print("Removing failed records...")
        for seq_record in SeqIO.parse(file, "fasta"):
            #print(seq_record.seq)
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
        #pat = "\[(?:[^\[\]]++|(?0))*+]"
        temp_new_file_string = re.sub("\[(?:[^\[\]]+|)]", "", file_string) # regex to find any [success] markers and remove them
        print("Done!")
        print("Removing non sequence strings...")
        new_file_string = temp_new_file_string.replace("killed", "")
        print("Done!")
        print("Removing duplicate sequences...")

        with open("temp3.fasta", "a") as f:
            f.write(new_file_string)
        f.close()

        '''
        sequences = {}
        for seq_record in SeqIO.parse("pymart_out.fasta", "fasta"):
            sequence = str(seq_record.seq).upper()
            if sequence not in sequences:
                sequences[sequence] = seq_record
            else:
                print("No dupe")
        print(sequences)
        with open("cleaned_utr_test.fasta", "w+") as output_file:
            # Just read the hash table and write on the file as a fasta format
            for sequence in sequences:
                output_file.write(">" + str(sequences[sequence].seq))
        '''
        # REMOVES DUPLICATES
        seen = []
        records = []
        for record in SeqIO.parse("temp3.fasta", "fasta"):
            if str(record.seq) not in seen:
                print(record.id, " not seen")
                seen.append(str(record.seq))
                records.append(record)
            else:
                print("Record seen")
            print("done checking")
        print("all done - writing")
        SeqIO.write(records, "cleaned_utr_test.fasta", "fasta")
        print("done write")

        print("Writing cleaned fasta file as cleaned_utr_test.fasta")
        print("Done!")
        print("Removing temp files...")
        os.remove("temp.fasta")
        os.remove("temp_2.fasta")
        print("FASTA file cleaned!")
        with open("removed_during_cleaning.txt", "a") as f:
            for i in rm_list:
                f.write("%s\n" % i)
        f.close()


