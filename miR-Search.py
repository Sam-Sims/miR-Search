#!/usr/bin/env python3
import argparse
import os, sys
import pickle

from modules import pymart, ggplot_builder as pb, sequence_handler, targetsearch, icshapeProcess as sp

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



def init_argparse():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    pymart_parser = subparsers.add_parser('pymart')
    search_parser = subparsers.add_parser('mir-search')
    format_parser = subparsers.add_parser('ggplot-format')
    process_parser = subparsers.add_parser('icshape-align')

    # PYMART ARGS
    # parser.add_argument("-d", "--download", help="Runs the pymart downloader. Requires --input.", action="store_true")
    pymart_parser.add_argument("-i", "--input",
                               help="Specifies a list of ensembl gene IDs to download in CSV format. NOTE the first column must contain the IDs.",
                               required=False)
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
    pymart_parser.add_argument("--verify",
                               help="Will clean the specified FASTA file. Recomended to be run after downloading.",
                               metavar="file", required=False)
    pymart_parser.add_argument("--url",
                               help="Specifies the url to append xml query to. Defaults to -t 10 (recomended not to change this",
                               default="http://www.ensembl.org/biomart/martservice?query=", required=False)
    pymart_parser.add_argument("--split", help="Splits the output into seperate files. Default FALSE. Ignores --output",
                               default=False, required=False, action="store_true")
    pymart_parser.add_argument('-m', '--mir', action="store_true",
                               help="Download miRNA sequences.", required=False)
    pymart_parser.add_argument('-v', '--verbosity', action="count",
                               help="Increase output verbosity (e.g., -vv is more than -v)")

    # MIRSEARCH ARGS
    search_parser.add_argument("-i", "--input",
                               help="Specifies the input for miR-Search. First: miRNA file; Second: UTR file. Files in FASTA format.",
                               required=True, nargs=2, metavar=("miRNA", "UTR"))
    search_parser.add_argument("-m", "--mir",
                               help="Use if input mir is mirbase ID",
                               required=False, action="store_true")
    search_parser.add_argument('-v', '--verbosity', action="count", help="increase output verbosity (e.g., -vv is more than -v)")
    search_parser.add_argument('--cache', action="store_true",
                               help="Will use the gene_dict cache if there is one")

    # GGPLOT-FORMAT ARGS
    format_parser.add_argument("-t", "--target",
                               help="Runs the module in target mode - to process the output from the mir-search module. First: sleuth results; Second: mirR-Search output. Files in CSV format.",
                               required=False, nargs=2, metavar=("sleuth", "miR-Search"))
    format_parser.add_argument('-v', '--verbosity', action="count",
                               help="increase output verbosity (e.g., -vv is more than -v)")
    format_parser.add_argument("-o", "--output", help="Specifies the location of the output file",
                               default="out.csv", required=False)
    format_parser.add_argument("-s", "--shape", help="Runs the module in shape mode - to process the output from the icshape-align module. First: sleuth results; Second: icshape-align output. Files in CSV format.", required=False, nargs=2, metavar=("sleuth", "icshape-align"))

    # icSHAPE-align ARGS
    process_parser.add_argument("-s", "--shape",
                               help="Specifies icSHAPE input file.",
                               required=False)
    process_parser.add_argument("-u", "--utr",
                                help="Specifies utr input file.",
                                required=False)
    process_parser.add_argument("-t", "--target",
                                help="Specifies target input file.",
                                required=False)
    process_parser.add_argument('-v', '--verbosity', action="count",
                               help="increase output verbosity (e.g., -vv is more than -v)")
    process_parser.add_argument("-r", "--ret",
                                help="Returns a list of transcripts from the icSHAPE output",
                                required=False)
    process_parser.add_argument("-m", "--match",
                                help="Generates a new FASTA file containg only UTR sequences for that have shape coverage for. Takes UTR file from pymart as input - run icshape-align --ret first.",
                                required=False)


    args = parser.parse_args()
    return args

def run_search(args):
    input_mir = args.input[0]
    input_utr = args.input[1]

    # Load FASTA files
    fp = sequence_handler.FASTAParse(v_print)
    utr_seq_ob_list = fp.read_multi_3utr(input_utr)  # returns list of biopython seq objects for each record in utr file
    mir_name = input_mir # store name for later
    if args.mir:
        print("test")
        input_mir = fp.read_mirbase(input_mir)
    else:
        input_mir = fp.read_mir(input_mir)

    # process mir
    mir = sequence_handler.MicroRNA(input_mir) # create a sequence handler object of input mir
    mir.auto_process() # trims, finds reverse complement and back transcribes mir object

    # create target search object - might want to rework this feels weird
    ts = targetsearch.TargetSearch(utr_seq_ob_list, mir_name, v_print)

    # returns dict where key = gene id/name, value = location for all targets
    sixmer_target_list = ts.search_6mer(mir.find_6mer())
    sevenmera1_target_list = ts.search_7mera1(mir.find_7mera1())
    sevenmerm8_target_list = ts.search_7merm8(mir.find_7merm8())
    eightmer_target_list = ts.search_8mer(mir.find_8mer())

    # returns a dict containg each gene as key and a list of 4 bools as the value - bools indicate if gene has 6,7a1,7m8 or 8mer target site (in that order) for the input mir
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

    # logic to determine single target site based on gene dict - this is because sites can be subsets of eachother. i.e 6mer subset of 7mer and 8mer
    targets = ts.calc_gene_targets(gene_dict)
    ts.print_targets(targets)


def run_pymart(args):
    pm = pymart.PYMart(args.url, args.input, args.output, args.split, args.mir)
    #pm.run_check()
    if args.scrub:
        print(args.scrub)
        pm.clean_utr(args.scrub)
    elif args.verify:
        pm.run_verify(args.verify)
    elif args.mir:
        pm.set_split()
        pm.download_mart(int(args.threads))
    else:
        pm.download_mart(int(args.threads))
        if args.check:
            pm.run_check()
        if args.auto:
            pm.run_check()
            pm.clean_utr(args.output)


def run_process(args):
    if args.ret:
        sp.return_transcripts(args.ret)
    elif args.match:
        sp.match_transcript_id(args.match)
    elif args.shape and args.target and args.utr:
        sp.auto_process(args)
    else:
        print("Please speicify a mode")



def run_format(args):
    if args.shape:
        sleuth = pb.prepare_sleuth_results_icshape(args.shape[0])
        with open(args.shape[1], 'rb') as handle:
            percent_dict = pickle.load(handle)
        handle.close()
        shape = pb.prepare_shape_data(percent_dict)
        pb.merge_shape(sleuth, shape, args.output)
    elif args.target:
        sleuth = pb.prepare_sleuth_results(args.input[0])
        targets = pb.prepare_targets(args.input[1])
        pb.merge(sleuth, targets, args.output)
    else:
        print("Error - mode not supplied. Use -t or -s to specify a mode. -h for more information")
        sys.exit()


def main():
    args = init_argparse()
    if args.command == "pymart":
        run_pymart(args)
    elif args.command == "mir-search":
        run_search(args)
    elif args.command == "ggplot-format":
        run_format(args)
    elif args.command == "icshape-align":
        run_process(args)
    else:
        print("Please specify a module to run. miR-Search -h for more details")
        sys.exit()

    if args.verbosity:  # if verbose the define print function use v_print(level, string) level =1 info, 2=warn, 3=error
        def _v_print(*verb_args):
            if verb_args[0] > (3 - args.verbosity):
                print(verb_args[1])
    else:
        _v_print = lambda *a: None  # do nothinh

    global v_print
    v_print = _v_print


if __name__ == "__main__":
    main()
