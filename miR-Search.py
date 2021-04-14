import argparse
import os
import pickle

import sequence_handler
import targetsearch
from pymart import pymart

v_print = None


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



def run_search(args):
    fp = sequence_handler.FASTAParse()
    input_mir = args.input[0]
    input_utr = args.input[1]
    mir_name = input_mir
    input_mir = fp.read_mir(input_mir)

    mir = sequence_handler.MicroRNA(input_mir)
    mir.auto_process()

    utr_seq_ob_list = fp.read_multi_3utr(input_utr, v_print)  # list of biopython seq objects for each record in master utr file

    ts = targetsearch.TargetSearch(utr_seq_ob_list, mir_name, v_print)

    sixmer_target_list = ts.search_6mer(mir.find_6mer())
    sevenmera1_target_list = ts.search_7mera1(mir.find_7mera1())
    sevenmerm8_target_list = ts.search_7merm8(mir.find_7merm8())
    eightmer_target_list = ts.search_8mer(mir.find_8mer())
    if args.cache:
        if not os.path.isfile("gene_id_dict.pickle"):
            v_print(2, "No cache found despite using --cache") # ignore pycharm error here
            gene_dict = ts.generate_gene_value_dict(sixmer_target_list, sevenmera1_target_list, sevenmerm8_target_list,
                                                    eightmer_target_list)
        else:
            print("Cache found!")
            print("Reading from cache...")
            with open('gene_id_dict.pickle', 'rb') as handle:
                gene_dict = pickle.load(handle)

            handle.close()
    else:
        gene_dict = ts.generate_gene_value_dict(sixmer_target_list, sevenmera1_target_list, sevenmerm8_target_list,
                                                eightmer_target_list)

    targets = ts.calc_gene_targets(gene_dict)
    ts.print_targets(targets)


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
                               help="Specifies the url to append xml query to. Defaults to -t 10 (recomended not to change this",
                               default="http://www.ensembl.org/biomart/martservice?query=", required=False)
    pymart_parser.add_argument("--split", help="Splits the output into seperate files. Default FALSE. Ignores --output",
                               default=False, required=False, action="store_true")
    pymart_parser.add_argument('-v', '--verbosity', action="count",
                               help="increase output verbosity (e.g., -vv is more than -v)")

    # MIRSEARCH ARGS
    search_parser.add_argument("-i", "--input",
                               help="Specifies the input for miR-Search. First: miRNA file; Second: UTR file. Files in FASTA format.",
                               required=True, nargs=2, metavar=("miRNA", "UTR"))
    search_parser.add_argument('-v', '--verbosity', action="count", help="increase output verbosity (e.g., -vv is more than -v)")
    search_parser.add_argument('--cache', action="store_true",
                               help="Will use the gene_dict cache if there is one")

    args = parser.parse_args()
    return args


def run_pymart(args):
    pm = pymart.PYMart(args.url, args.input, args.output, args.split)
    if args.scrub:
        print(args.scrub)
        pm.clean_utr(args.scrub)
    else:
        pm.download_mart(int(args.threads))
        if args.check:
            pm.run_check()
        if args.auto:
            pm.run_check()


def main():
    args = init_argparse()

    if args.verbosity: # if verbose the define print function use v_print(level, string) level =1 info, 2=warn, 3=error
        def _v_print(*verb_args):
            if verb_args[0] > (3 - args.verbosity):
                print(verb_args[1])
    else:
        _v_print = lambda *a: None  # do nothinh

    global v_print
    v_print = _v_print

    if args.command == "pymart":
        v_print(1, "egg")
        run_pymart(args)
    elif args.command == "search":
        print("Run search")
        run_search(args)
    else:
        print("Please specify a task using pymart or search. Usage miR-Search pymart; miR-Search search")


if __name__ == "__main__":
    main()
