import argparse

from Bio import SeqIO

import sequence_handler
import targetsearch
from pymart import pymart


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

    utr_seq_ob_list = fp.read_multi_3utr(input_utr)  # list of biopython seq objects for each record in master utr file

    ts = targetsearch.TargetSearch(utr_seq_ob_list)
    print(mir.find_6mer())
    print(mir.find_7mera1())
    print(mir.find_7merm8())
    print(mir.find_8mer())
    # sixmer_target_list = ts.search_6mer(mir.find_6mer())
    # sevenmera1_target_list = ts.search_7mera1(mir.find_7mera1())
    # sevenmerm8_target_list = ts.search_7merm8(mir.find_7merm8())
    # eightmer_target_list = ts.search_8mer(mir.find_8mer())
    # gene_dict = ts.generate_gene_value_dict(sixmer_target_list, sevenmera1_target_list, sevenmerm8_target_list, eightmer_target_list)
    ts.calc_gene_targets()

    # SINGLE UTR SEARCH
    # utr = fp.read_3utr(input_utr)
    # so = sequence_handler.MicroRNASearch(utr)
    # is6mer = so.search_6mer(mir.find_6mer())
    # print(is6mer)


def init_argparse():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    pymart_parser = subparsers.add_parser('pymart')
    search_parser = subparsers.add_parser('search')

    # PYMART ARGS
    # parser.add_argument("-d", "--download", help="Runs the pymart downloader. Requires --input.", action="store_true")
    pymart_parser.add_argument("-i", "--input",
                               help="Specifies a list of ensembl gene IDs to download in CSV format. NOTE the first column must contain the IDs.",
                               required=True)
    pymart_parser.add_argument("-o", "--output", help="Specifies the location of the output file",
                               default="pymart_out.fasta", required=False)
    pymart_parser.add_argument("-c", "--check",
                               help="Will run a comparison check between the output gene sequence and the input gene list.",
                               action="store_true", required=False)
    pymart_parser.add_argument("-a", "--auto",
                               help="Runs pymart in automode, performing the comparison check and file tidy (same as running -c -s) (recommended)",
                               action="store_true", required=False)
    pymart_parser.add_argument("-t", "--threads",
                               help="Specifies the number of threads pymart will use when downloading. Default = 50",
                               default=50, required=False)
    pymart_parser.add_argument("-s", "--scrub",
                               help="Will clean the specified FASTA file. Recomended to be run after downloading.",
                               metavar="file", required=False)
    pymart_parser.add_argument("--url",
                               help="Specifies the url to append xml query to. Defaults to http://www.ensembl.org/biomart/martservice?query= (recomended not to change this",
                               default="http://www.ensembl.org/biomart/martservice?query=", required=False)
    pymart_parser.add_argument("--split", help="Splits the output into seperate files. Default FALSE. Ignores --output",
                               default=False, required=False, action="store_true")

    # MIRSEARCH ARGS

    args = parser.parse_args()
    return args


def run_pymart(args):
    pm = pymart.PYMart(args.url, args.input, args.output, args.split)
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


global v_print


def main():
    # print_menu()
    # ans = input()
    # check_menu_choice(ans)
    # run_search("test-data/hsa-miR-451a.fasta", "cleaned_utr_test.fasta")
    args = init_argparse()
    if args.command == "pymart":
        run_pymart(args)
    elif args.command == "search":
        print("Run search")
    else:
        print("Please specify a task using pymart or search. Usage miR-Search pymart; miR-Search search")


if __name__ == "__main__":
    main()
