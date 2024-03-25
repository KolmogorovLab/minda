#! /usr/bin/env python
import sys
import os
import argparse
import logging
import pandas as pd

from minda.__version__ import __version__
from minda.decompose import get_caller_name, get_df, get_intersected_df, get_decomposed_dfs
from minda.ensemble import get_support_df as get_ensemble_support_df
from minda.truthset import get_support_df as get_truthset_support_df
from minda.truthset import get_base_df
from minda.stats import get_results

logger = logging.getLogger()


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    if not debug:
        level = logging.INFO
    
    console_log.setLevel(level)
    file_handler.setLevel(level)
    
    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def _version():
    return __version__

def run(args):
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    log_file = os.path.join(args.out_dir, "minda.log")
    _enable_logging(log_file, debug=False, overwrite=True)

    version =  _version()
    logger.info("MindaV" + version)

    if args.command == 'truthset':
        base = args.base

    # check whether input is tsv or list of vcfs
    if args.tsv != None:
        vcf_df = pd.read_csv(args.tsv, sep='\t', header=None)
        
        vcf_list = vcf_df.iloc[:,0].to_list()
        tsv_directory = os.path.abspath(args.tsv)
        vcf_list = [os.path.abspath(path) if not os.path.isabs(path) else path for path in vcf_list]
        
        column_count = vcf_df.shape[1]
        if column_count >= 2:
            caller_names = vcf_df.iloc[:,1].fillna("unknown").to_list()
        else:
            caller_names = []
        
        if column_count == 3:
            prefixes = vcf_df.iloc[:,2].fillna("unk").to_list()
            caller_names = [prefixes[i] + "_" + caller_names[i] for i in range(len(caller_names))]
 
    else:
        column_count = 1
        vcf_list = args.vcfs
        caller_names = []

    
    if args.command == "ensemble" and len(vcf_list) < 2:
        sys.exit("Provide a minimum of 2 VCF files.")
    elif args.command == "ensemble" and args.min_support != None and len(vcf_list) < args.min_support:
        sys.exit("Number of VCF files should be less than or equal minimum number of support.")
    elif len(vcf_list) < 1 and args.command == "truthset":
        sys.exit("Provide a minimum of 1 comparison VCF file.") 

    if args.command == 'truthset':
        vcf_list.append(base)
        if len(caller_names) > 0:
            caller_names.append("base")
    
    if caller_names == []:
        for i in range(len(vcf_list)):
            vcf = vcf_list[i]
            caller_name = get_caller_name(vcf)
            caller_names.append(caller_name)

    if len(caller_names) > len(set(caller_names)):
        caller_names = [chr(ord('A') + i) + "_" + caller_names[i] for i in range(len(caller_names))]
        prefixed = True
    elif len(caller_names) == len(set(caller_names)) and column_count == 3:
        prefixed = True
    else:
        prefixed = False    

    max_svlengths =[]
    decomposed_dfs_list = []    
    for i in range(len(vcf_list)):
        caller_name = caller_names[i]
        vcf = vcf_list[i]
        if args.bed == None:
            df = get_df(vcf)
        else:
            df = get_intersected_df(vcf, args.bed)
        decomposed_dfs = get_decomposed_dfs(caller_name, df, args.filter, args.min_size, prefixed, args.vaf, args.sample_name, args.out_dir)
        decomposed_dfs_list.append(decomposed_dfs[:2])
        max_svlengths.append(decomposed_dfs[2])

    max_len = max(max_svlengths)

    if args.command == 'ensemble' and args.conditions != None:
        conditions = eval(args.conditions)
        support_df = get_ensemble_support_df(decomposed_dfs_list, caller_names, args.tolerance, conditions, args.vaf, args.command, args.out_dir, args.sample_name, args, version, args.multimatch)
        results = get_results(decomposed_dfs_list, support_df, caller_names, args.out_dir, args.sample_name, max_len, args.tolerance, args.vaf, args.command, args)
        logger.info(f"\n{results[0]}")
    
    elif args.command == 'ensemble' and args.min_support != None:
        conditions = eval(f"[[caller_names,'>=', {args.min_support}]]")
        support_df = get_ensemble_support_df(decomposed_dfs_list, caller_names, args.tolerance, conditions, args.vaf, args.command, args.out_dir, args.sample_name, args, version, args.multimatch)
        results = get_results(decomposed_dfs_list, support_df, caller_names, args.out_dir, args.sample_name, max_len, args.tolerance, args.vaf, args.command, args)
        logger.info(f"\n{results[0]}")

    else:
        base_df = get_base_df(decomposed_dfs_list, args.tolerance, args.multimatch) 
        support_df = get_truthset_support_df(base_df, caller_names, args.vaf, args.out_dir, args.sample_name)
        results = get_results(decomposed_dfs_list, support_df, caller_names, args.out_dir, args.sample_name, max_len, args.tolerance, args.vaf,args.command, args)
        logger.info(f"\n{results[0]}")


def main():
    parser=argparse.ArgumentParser(description="Minda - VCF evaluation tool for germline and somatic structural variant callers")
    subparser=parser.add_subparsers(dest="command")

     #defaults  ------------------------------------------------
    FILTER = ["PASS"]
    TOLERANCE = 500
    
    # TRUTHSET ------------------------------------------------
    truthset = subparser.add_parser("truthset", help='benchmark VCF(s) against a base VCF')

    # required arguements
    truthset.add_argument("--out_dir", help='path to out directory', dest="out_dir", type=str, required=True)
    truthset.add_argument("--base", help='path of base VCF', dest="base", type=str, required=True)
    
    # mutally exclusive arguments
    truthset_input = truthset.add_mutually_exclusive_group(required=True)
    truthset_input.add_argument('--tsv', action="store", dest="tsv", help="tsv file path")
    truthset_input.add_argument('--vcfs', action="store", dest="vcfs", nargs="+", help="vcf file path(s)")

    # # optional arguments 
    truthset.add_argument("--bed", help=f'path to bed file for filtering records with BedTool intersect', dest="bed", type=str)
    truthset.add_argument("--filter", help=f'filter records by FILTER column; default="{FILTER}"', dest="filter", type=str, nargs="*", default=FILTER) 
    truthset.add_argument("--min_size", help=f'filter records by SVSIZE in INFO column', dest="min_size", type=int)
    truthset.add_argument("--tolerance", help=f'maximum allowable bp distance between base and caller breakpoint; default={TOLERANCE}', dest="tolerance", type=int, default=TOLERANCE)
    truthset.add_argument("--sample_name", help=f'name of sample', dest="sample_name", type=str)
    truthset.add_argument("--vaf", help=f'filter out records below a given VAF treshold', dest="vaf", type=float)
    truthset.add_argument("--multimatch", help=f'allow more than one record from the same caller to match a single truthset record', dest="multimatch", action='store_true')

    # ENSEMBLE ------------------------------------------------
    ensemble = subparser.add_parser("ensemble", help='create an ensemble call list from multiple VCF and, optionally, benchmark each VCF against')

    # required arguements
    ensemble.add_argument("--out_dir", help='path to out directory', dest="out_dir", type=str, required=True)
    
    # mutally exclusive arguments
    ensemble_input = ensemble.add_mutually_exclusive_group(required=True)
    ensemble_input.add_argument('--tsv', action="store", dest="tsv", help="tsv file path")
    ensemble_input.add_argument('--vcfs', action="store", dest="vcfs", nargs="+", help="vcf file path(s)")

    ensemble_support = ensemble.add_mutually_exclusive_group(required=True)
    ensemble_support.add_argument("--conditions", help=f'specific conditions to support a call', dest="conditions", type=str)
    ensemble_support.add_argument("--min_support", help=f'minimumn number of callers to support a call', dest="min_support", type=int)

    # optional arguments 
    ensemble.add_argument("--bed", help=f'path to bed file for filtering records with BedTool intersect', dest="bed", type=str)
    ensemble.add_argument("--filter", help=f'filter records by FILTER column; default="{FILTER}"', dest="filter", type=str, nargs="*", default=FILTER) 
    ensemble.add_argument("--min_size", help=f'filter records by SVSIZE in INFO column', dest="min_size", type=int)
    ensemble.add_argument("--tolerance", help=f'maximum allowable bp distance between base and caller breakpoint; default={TOLERANCE}', dest="tolerance", type=int, default=TOLERANCE)
    ensemble.add_argument("--sample_name", help=f'name of sample', dest="sample_name", type=str)
    ensemble.add_argument("--vaf", help=f'filter out records below a given VAF treshold', dest="vaf", type=float)
    ensemble.add_argument("--multimatch", help=f'allow more than one record from the same caller to match a single ensemble record', dest="multimatch", action='store_true')

    # ------------------------------------------------
    args, remaining_args = parser.parse_known_args()
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)
