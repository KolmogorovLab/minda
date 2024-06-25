#!/usr/bin/env python3

import sys
import os
import json

from intervaltree import Interval, IntervalTree
from collections import defaultdict, namedtuple


def parse_repeatmasker(filename):
    """
    Parses repeatmakser 'out' file and returns index of annotated repeats
    """
    chr_trees = defaultdict(IntervalTree)
    num_rec = 0

    for line in open(filename, "r"):
        fields = line.split()
        if len(fields) != 15:
            continue

        divergence, chrom, start, end, repeat_id, family = float(fields[1]), fields[4], int(fields[5]), int(fields[6]), fields[9], fields[10]
        if end - start < MIN_REPEAT:
            continue

        chr_trees[chrom][start:end] = (repeat_id, family, start)
        num_rec += 1
    print(num_rec)

    return chr_trees


def get_bed_intervals(filename):
    """
    Parses bed file with intervals and returns index
    """
    chr_trees = defaultdict(IntervalTree)

    for line in open(filename, "r"):
        if line.startswith("#"):
            continue

        fields = line.strip().split()
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        chr_trees[chrom][start:end] = (start, end)

    return chr_trees


def get_vcf_breakpoints(filename):
    """
    Extracts breakpoint coordinates from a vcf
    """
    chr_trees = defaultdict(IntervalTree)
    for line in open(filename, "r"):
        if line.startswith("#"):
            continue

        fields = line.strip().split()
        chr_1, pos_1, info = fields[0], int(fields[1]), fields[7]
        chr_2, pos_2 = None, None
        tags = info.split(";")
        for tag in tags:
            if tag.startswith("CHR2"):
                chr_2 = tag[5:]
            if tag.startswith("END"):
                pos_2 = int(tag[4:])

        if chr_2 is None:
            chr_2 = chr_1
        chr_trees[chr_1][pos_1 : pos_1 + 1] = (chr_1, pos_1)
        if pos_2 is not None:
            chr_trees[chr_2][pos_2 : pos_2 + 1] = (chr_2, pos_2)

    return chr_trees


MindaEntry = namedtuple("MindaEntry", ["minda_num", "minda_str", "chr_x", "pos_x", "list_x", "chr_y", "pos_y",
                                       "list_y", "sv_type", "sv_len", "vaf", "support", "original_line", "is_ensemble"])
def parse_minda_csv(filename, num_callsets):
    """
    Parses Minda output file
    """
    callset_list = None
    minda_num = 0
    minda_entries = {}

    TOOLS_BEGIN = 13

    for line in open(filename, "r"):
        if line.startswith("#"):
            callset_list = line.strip().split("\t")[TOOLS_BEGIN : TOOLS_BEGIN + num_callsets]
            print("Callsets:", callset_list)
            continue

        fields = line.strip().split("\t")

        support_dict = {}
        for (caller, supp) in zip(callset_list, fields[TOOLS_BEGIN : TOOLS_BEGIN + num_callsets]):
            support_dict[caller] = True if supp != "False" else False
        is_ensemble = (fields[12] == "True")

        #if not len(fields[10]):
        #    print(line)
        #    continue

        #VAF not always available, default to 1.0
        try:
            vaf = float(fields[10])
        except ValueError:
            vaf = 1.0

        entry = MindaEntry(minda_num=minda_num, minda_str=fields[10],
                           chr_x=fields[0], pos_x=int(fields[1]), list_x=fields[3],
                           chr_y=fields[4], pos_y=int(fields[5]), list_y=fields[7],
                           sv_type=fields[8], sv_len=int(fields[9]), vaf=vaf,
                           support=support_dict, original_line=line.strip(), is_ensemble=is_ensemble)
        minda_entries[minda_num] = entry
        minda_num += 1

    return minda_entries


def get_confident_calls(minda_records):
    """
    Extracts confident calls based on given minimum support
    """
    confident_calls = set()
    for rec in minda_records.values():
        if rec.is_ensemble:
            confident_calls.add(rec.minda_num)

    """
    callset_list = list(next(iter(minda_records.values())).support.keys())
    for rec in minda_records.values():
        support_tools = set(x for x in rec.support if rec.support[x])
        against_tools = set(callset_list) - support_tools
        techs = set([x.split("_")[-1] for x in support_tools])
        if len(support_tools) >= min_tools and len(techs) >= min_tech:
            confident_calls.add(rec.minda_num)
    """

    return confident_calls


def filter_calls(minda_records, min_vaf, min_sv_len, remove_ins):
    filtered_minda = {}
    for rec in minda_records.values():
        if rec.vaf < min_vaf:
            continue
        #if rec.pos_y == 0:
        #    continue
        if rec.chr_x == rec.chr_y and abs(rec.pos_y - rec.pos_x) < min_sv_len:
            continue
        if remove_ins and rec.sv_type == "INS":
            continue

        filtered_minda[rec.minda_num] = rec
    return filtered_minda


def stratify_breakends(minda_records, confident_calls, annotation_dir, remove_insertions):
    """
    Performs various types of stratification
    """
    repeatmasker_file = os.path.join(annotation_dir, REPEAT_MASKER)
    segdup_file = os.path.join(annotation_dir, SEGDUPS_BED)
    vntr_file = os.path.join(annotation_dir, VNTR_BED)
    chr_sizes_file = os.path.join(annotation_dir, CHR_LEN_BED)

    index_repeatmasker = parse_repeatmasker(repeatmasker_file)
    index_segdup = get_bed_intervals(segdup_file)
    index_vntr = get_bed_intervals(vntr_file)
    #index_germline = get_vcf_breakpoints(germline_vcf)

    #chr sizes
    chr_sizes = defaultdict(int)
    for line in open(chr_sizes_file, "r"):
        fields = line.split()
        chr_sizes[fields[0]] = int(fields[1])

    #cluster indexing
    index_ensemble = defaultdict(IntervalTree)
    for entry in minda_records.values():
        index_ensemble[entry.chr_x][entry.pos_x : entry.pos_x + 1] = entry
        index_ensemble[entry.chr_y][entry.pos_y : entry.pos_y + 1] = entry

    strat_categories = ["hom_repeat", "segdup", "vntr", "low_vaf", "low_len",
                        "bnd_dup", "bnd_chain"]
    if not remove_insertions:
        strat_categories = ["insertion"] + strat_categories
    strat_entries = defaultdict(set)

    def _get_intervals(entry, index, threshold):
        ovlp_1 = index[entry.chr_x][entry.pos_x - threshold : entry.pos_x + threshold]
        ovlp_2 = index[entry.chr_y][entry.pos_y - threshold : entry.pos_y + threshold]
        return [o[2] for o in ovlp_1], [o[2] for o in ovlp_2]

    def _support_tools(rec):
        return set(x for x in rec.support if rec.support[x])

    def _coords_number(entries):
        clusters = []
        for e in entries:
            for (c, p) in [(e.chr_x, e.pos_x), (e.chr_y, e.pos_y)]:
                match = False
                for cl in clusters:
                    if cl[0] == c and abs(cl[1] - p) <= BND_AREA:
                        match = True
                if not match:
                    clusters.append((c, p))
        return len(clusters)

    def _is_duplicate(e1, e2):
        match_x = e1.chr_x == e2.chr_x and abs(e1.pos_x - e2.pos_x) <= BND_AREA
        match_y = e1.chr_y == e2.chr_y and abs(e1.pos_y - e2.pos_y) <= BND_AREA
        cross_x = e1.chr_x == e2.chr_y and abs(e1.pos_x - e2.pos_y) <= BND_AREA
        cross_y = e1.chr_y == e2.chr_x and abs(e1.pos_y - e2.pos_x) <= BND_AREA
        return (match_x and match_y) or (cross_x and cross_y)

    for entry in minda_records.values():
        #Homologous repeats at breakends
        ovlp_1, ovlp_2 = _get_intervals(entry, index_repeatmasker, REPEAT_AREA)
        for x in ovlp_1:
            for y in ovlp_2:
                if x[1] == y[1] and x[2] != y[2]:   #same repeat family, but different repeat
                    strat_entries[entry.minda_num].add("hom_repeat")

        #same segdup section
        segdup_1, segdup_2 = _get_intervals(entry, index_segdup, REPEAT_AREA)
        if len(set(s[0] for s in segdup_1) & set(s[0] for s in segdup_2)) > 0:
        #if len(segdup_1) > 0 and len(segdup_2) > 0:
            strat_entries[entry.minda_num].add("segdup")

        #same vntr section
        vntr_1, vntr_2 = _get_intervals(entry, index_vntr, REPEAT_AREA)
        if len(set(s[0] for s in vntr_1) & set(s[0] for s in vntr_2)) > 0:
        #if len(vntr_1) > 0 and len(vntr_2) > 0:
            strat_entries[entry.minda_num].add("vntr")

        #low-ish vaf
        if entry.vaf < LOW_VAF:
            strat_entries[entry.minda_num].add("low_vaf")

        if entry.chr_x == entry.chr_y and abs(entry.pos_y - entry.pos_x) < LOW_LEN:
            strat_entries[entry.minda_num].add("low_len")

        #is insertion
        if entry.sv_type == "INS":
            strat_entries[entry.minda_num].add("insertion")

        """
        #near telomere
        if min(entry.pos_x, chr_sizes[entry.chr_x] - entry.pos_x) < TELOMERE_LEN or \
                min(entry.pos_y, chr_sizes[entry.chr_y] - entry.pos_y) < TELOMERE_LEN:
            strat_entries[entry.minda_num].add("telomere")
        """

        """
        #near germline SV breakpoints
        germ_1, germ_2 = _get_intervals(entry, index_germline, BND_AREA)
        if len(germ_1) > 0 and len(germ_2) > 0:
            strat_entries[entry.minda_num].add("germline")

        """

        """
        #near multiple condifent breakpoints
        if len(set(o.minda_num for o in ens_bnds_1 if o.minda_num in confident_calls)) > 1 or \
                len(set(o.minda_num for o in ens_bnds_2 if o.minda_num in confident_calls)) > 1:
            strat_entries[entry.minda_num].add("truth_cluster")
        """

        ens_bnds_1, ens_bnds_2 = _get_intervals(entry, index_ensemble, BND_AREA)

        #duplication
        for r in ens_bnds_1:
            if r != entry and _is_duplicate(r, entry):
                strat_entries[entry.minda_num].add("bnd_dup")
                strat_entries[r.minda_num].add("bnd_dup")
        for r in ens_bnds_2:
            if r != entry and _is_duplicate(r, entry):
                strat_entries[entry.minda_num].add("bnd_dup")
                strat_entries[r.minda_num].add("bnd_dup")

        #3+ chain of breakends
        left_chain, right_chain = None, None
        extra_chain_entries = set()
        for r in ens_bnds_1:
            if r != entry and not _is_duplicate(r, entry):
                left_chain = r
                extra_chain_entries.add(r.minda_num)
        for r in ens_bnds_2:
            if r != entry and not _is_duplicate(r, entry):
                right_chain = r
                extra_chain_entries.add(r.minda_num)

        if None not in [left_chain, right_chain] and left_chain != right_chain:
            strat_entries[entry.minda_num].add("bnd_chain")
            for e in extra_chain_entries:
                strat_entries[e].add("bnd_chain")

        #print(entry.original_line + "\t" + ",".join(list(strat_entries[entry.minda_num])))

    return strat_categories, strat_entries


def compute_fp_fn(minda_records, confident_calls, strat_category, strat_entries, print_table):
    callset_list = list(next(iter(minda_records.values())).support.keys())
    callset_list.sort(key=lambda x: x.split("_")[0])
    tools_tp, tools_fp, tools_fn, tools_f1 = defaultdict(set), defaultdict(set), defaultdict(set), defaultdict(float)
    tools_recall, tools_precision = defaultdict(set), defaultdict(set)

    for rec in minda_records.values():
        support_tools = set(x for x in rec.support if rec.support[x])
        against_tools = set(callset_list) - support_tools

        if strat_category is not None:
            estrat = strat_entries[rec.minda_num]
            if len(estrat) == 0:
                estrat = set(["Unclassified"])
            if strat_category not in estrat:
                continue
            #if strat_category == "Unclassified" and rec.minda_num in confident_calls:
            #    print("FN", rec.chr_x, rec.pos_x, rec.chr_y, rec.pos_y, against_tools)

        if rec.minda_num in confident_calls:
            for tool in support_tools:
                tools_tp[tool].add(rec.minda_num)
            for tool in against_tools:
                tools_fn[tool].add(rec.minda_num)
        else:
            for tool in support_tools:
                tools_fp[tool].add(rec.minda_num)

    if print_table:
        print(f"\n=== Stratifying by: {strat_category} ===\n")
        print("#Tool\tTP\tFP\tFN\tprecision\trecall\tF1")

    for tool in callset_list:
        if len(tools_tp[tool]) > 0:
            precision = len(tools_tp[tool]) / (len(tools_tp[tool]) + len(tools_fp[tool]))
            recall = len(tools_tp[tool]) / (len(tools_tp[tool]) + len(tools_fn[tool]))
            f1_score = 2 * precision * recall / (precision + recall)
        else:
            precision, recall, f1_score = 0, 0, 0
        tp, fp, fn, = len(tools_tp[tool]), len(tools_fp[tool]), len(tools_fn[tool])
        tools_f1[tool] = f1_score
        tools_recall[tool] = recall
        tools_precision[tool] = precision

        if print_table:
            if PRETTY_PRINT:
                print(f"{tool:20s}\t{tp}\t{fp}\t{fn}\t{precision:.4f}\t{recall:.4f}\t{f1_score:.4f}")
            else:
                print(f"{tool}\t{tp}\t{fp}\t{fn}\t{precision:.4f}\t{recall:.4f}\t{f1_score:.4f}")

    return tools_tp, tools_fp, tools_fn, tools_recall, tools_precision, tools_f1


def summary_errors(minda_records, confident_calls, categories, strat_entries):
    callset_list = list(next(iter(minda_records.values())).support.keys())
    callset_list.sort(key=lambda x: x.split("_")[0])

    by_tool_fp = defaultdict(dict)
    by_tool_fn = defaultdict(dict)
    by_tool_tp = defaultdict(dict)
    by_tool_f1_score = defaultdict(dict)
    by_tool_recall = defaultdict(dict)
    by_tool_precision = defaultdict(dict)
    for cat in categories:
        tp, fp, fn, recall, precision, f1_score = \
                compute_fp_fn(minda_records, confident_calls, cat, strat_entries, print_table=False)
        for tool in tp:
            by_tool_tp[tool][cat] = len(tp[tool])
            by_tool_fp[tool][cat] = len(fp[tool])
            by_tool_fn[tool][cat] = len(fn[tool])
            by_tool_f1_score[tool][cat] = "{:.4f}".format(f1_score[tool])
            by_tool_recall[tool][cat] = "{:.4f}".format(recall[tool])
            by_tool_precision[tool][cat] = "{:.4f}".format(precision[tool])

    def print_with(stats, title):
        print(f"\n\t == {title} == \n")
        #print("{:20s}\t".format("#Tool") + "\t".join(categories))
        print("{}\t".format("#Tool") + "\t".join(categories))
        for tool in callset_list:
            numbers = [str(stats[tool][cat]) for cat in categories]
            if PRETTY_PRINT:
                print(f"{tool:20s}\t" + "\t".join(numbers))
            else:
                print(f"{tool}\t" + "\t".join(numbers))

    print_with(by_tool_tp, "True positives")
    print_with(by_tool_fp, "False positives")
    print_with(by_tool_fn, "False negatives")
    print_with(by_tool_recall, "Recall")
    print_with(by_tool_precision, "Precision")
    print_with(by_tool_f1_score, "F1 scores")


#annotation filenames 
REPEAT_MASKER = "repeatmasker.out"
SEGDUPS_BED = "segdups.bed"
VNTR_BED = "trf.bed"
CHR_LEN_BED = "chr.fasta.fai"

#stratification
MIN_REPEAT = 100
REPEAT_AREA = 5
BND_AREA = 500
LOW_VAF = 0.10
LOW_LEN = 100
#TELOMERE_LEN = 10000

#evaluation
#MIN_TECH = 2
#MIN_TOOLS = 4
MIN_VAF = 0.00
MIN_SV_LEN = 50
REMOVE_INS = False
PRETTY_PRINT = True


def minda_stratification(annotation_dir, minda_support_tsv, num_callsets):
    print(f"Params: min_vaf:{MIN_VAF} insertions_removed:{REMOVE_INS}")

    minda_records = parse_minda_csv(minda_support_tsv, num_callsets)
    minda_records = filter_calls(minda_records, MIN_VAF, MIN_SV_LEN, remove_ins=REMOVE_INS)
    confident_calls = get_confident_calls(minda_records)
    compute_fp_fn(minda_records, confident_calls, None, None, print_table=True)

    strat_categories, strat_calls = stratify_breakends(minda_records, confident_calls, annotation_dir,
                                                       REMOVE_INS)
    categories_unk = strat_categories + ["Unclassified"]

    for category in categories_unk:
        compute_fp_fn(minda_records, confident_calls, category, strat_calls, print_table=True)
    summary_errors(minda_records, confident_calls, categories_unk, strat_calls)


def main():
    if len(sys.argv) != 4:
        print("Usage: minda_stratify.py annotation_dir minda_support_tsv num_callsets")
        return 1

    annotation_dir = sys.argv[1]
    #germline_vcf = sys.argv[2]
    minda_csv = sys.argv[2]
    num_callsets = int(sys.argv[3])
    minda_stratification(annotation_dir, minda_csv, num_callsets)


if __name__ == "__main__":
    main()
