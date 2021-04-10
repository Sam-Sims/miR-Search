from Bio import SeqIO
import sequence_handler, targetsearch
from pymart import pymart
import argparse


def read_mir_testing(test):
    for fasta_sequence in SeqIO.parse(open('test-data/site-testing/' + test), 'fasta'):
        print("Reading", fasta_sequence.seq)
        return fasta_sequence.seq


def read_3utr_testing(test):
    print("TEST")
    for fasta_sequence in SeqIO.parse("test-data/site-testing/" + test, "fasta"):
        print("Reading 3' UTR region... ")
        print(fasta_sequence.seq)
        return fasta_sequence.seq


def print_menu():
    print('''\
 ---------------------------------------------------------
            _ ____      ____                      _     
  _ __ ___ (_)  _ \    / ___|  ___  __ _ _ __ ___| |__  
 | '_ ` _ \| | |_) |___\___ \ / _ \/ _` | '__/ __| '_ \ 
 | | | | | | |  _ <_____|__) |  __/ (_| | | | (__| | | |
 |_| |_| |_|_|_| \_\   |____/ \___|\__,_|_|  \___|_| |_|
 
 ---------------------------------------------------------                                                
    ''')
    print("Please make a selection:")
    print("1. Download 3'UTR regions")
    print("2. Run Search")

def check_menu_choice(ans):
    try:
        if ans == '1':
            run_pymart()
        elif ans == '2':
            run_search()
        else:
            print('Error, a valid answer was not supplied!')
    except Exception as e:
        print('An error occurred' + str(e))


def run_search(input_mir, input_utr):
    fp = sequence_handler.FASTAParse()
    input_mir = fp.read_mir(input_mir)

    mir = sequence_handler.MicroRNA(input_mir)
    mir.auto_process()

    utr_seq_ob_list = fp.read_multi_3utr(input_utr) # list of biopython seq objects for each record in master utr file

    ts = targetsearch.TargetSearch(utr_seq_ob_list)
    ts.search_6mer(mir.find_6mer())
    ts.search_7mera1(mir.find_7mera1())
    ts.search_7merm8(mir.find_7merm8())
    ts.search_8mer(mir.find_8mer())
    print(mir.return_seq())



    # SINGLE UTR SEARCH
    #utr = fp.read_3utr(input_utr)
    #so = sequence_handler.MicroRNASearch(utr)
    #is6mer = so.search_6mer(mir.find_6mer())
    #print(is6mer)

    #is6mer = False
    is7merm8 = False
    is7mera1 = False
    is8mer = False


'''
    if sixmer[0]:
        print("6mers found at: ")
        print(sixmer[1])
        is6mer = True
        if sevenmerm8[0]:
            print("7mer-m8s found at: ")
            print(sevenmerm8[1])
            is7merm8 = True
        else:
            print("No 7mer-m8 sites found")
            if sevenmera1[0]:
                print("7mer-a1s found at: ")
                print(sevenmera1[1])
                is7mera1 = True
                if eightmer[0]:
                    print("8mers also found at: ")
                    print(eightmer[1])
                    is8mer = True
                else:
                    print("7mer no 8mer")
            else:
                print("No 7mer-a1 sites found")
                print("Checking for 8mers")
                if eightmer[0]:
                    print("8mers found at: ")
                    print(eightmer[1])
                    is8mer = True
                else:
                    print("No targets found")

    else:
        print("No seed sites found!")

    print("6mers: " + str(is6mer))
    print("7mer-a1: " + str(is7mera1))
    print("7mer-m8: " +str(is7merm8))
    print("8mer: " +str(is8mer))
'''

def init_argparse():
    parser = argparse.ArgumentParser()
    # PYMART ARGS
    parser.add_argument("-d", "--download", help="Runs the pymart downloader. Requires --input.", action="store_true")
    parser.add_argument("-i", "--input", help="Specifies the location of gene list in csv format. NOTE the first column must contain ensembl gene IDs")
    parser.add_argument("-o", "--output", help="Specifies the location of the FASTA output", default="pymart_out.fasta")
    parser.add_argument("-c", "--check", help="Will run a comparison check between the output gene sequence and the input gene list. Used with --download", action="store_true")
    parser.add_argument("-a", "--auto", help="Runs pymart in automode, performing the comparison check and file tidy (same as running -c -s) (recommended)", action="store_true")
    parser.add_argument("-t", "--threads", help="Specifies the number of threads pymart will use when downloading. Default = 50", default=50)
    parser.add_argument("-s", "--scrub", help="Will clean the specified FASTA file. Recomended to be run after downloading.", metavar="file")
    parser.add_argument("--url", help="Specifies the url to append xml query to. Defaults to http://www.ensembl.org/biomart/martservice?query= (recomended not to change this", default="http://www.ensembl.org/biomart/martservice?query=")
    #MIRSEARCH ARGS
    parser.add_argument("-m", "--search", help="Runs the miR-Search programme", nargs=2)

    args = parser.parse_args()
    if args.download and (args.input is None):
        parser.error("--download requires --file")
    return args


def run_pymart(args):
    pm = pymart.PYMart(args.url, args.input, args.output)
    if args.scrub:
        print(args.scrub)
        pm.clean_utr(args.scrub)
    else:
        pm.download_mart(args.threads)
        if args.check:
            pm.run_check()
        if args.auto:
            pm.run_check()
            pm.clean_utr(args.output)


def main():
    #print_menu()
    #ans = input()
    #check_menu_choice(ans)
    run_search("test-data/hsa-miR-451a.fasta", "cleaned_utr_test.fasta")
    '''
    args = init_argparse()
    if args.download or args.scrub:
        run_pymart(args)
    else:
        print("Please specify a task using --download or --search")
    '''



if __name__ == "__main__":
    main()