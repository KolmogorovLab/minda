import sys
import os
from datetime import datetime
import pandas as pd 
from collections import Counter

def _get_tp_fn_fp(support_df, decomposed_dfs, caller_name, vaf, command):

    paired_df = decomposed_dfs[0].merge(decomposed_dfs[1], on='Minda_ID')

    if command == "ensemble":
        base_df = support_df[support_df['ensemble'] == True]
    else:
        base_df = support_df[support_df.iloc[:,13] == True]
            
    tp_ids = [id for ids in base_df["Minda_IDs"].to_list() for id in ids]
    
    # create tp, fn, fp dfs
    if command == "ensemble":
        fn_columns = ['#CHROM_x', 'POS_x', 'locus_group_x', 'ID_list_x', \
                          '#CHROM_y', 'POS_y', 'locus_group_y', 'ID_list_y', \
                          'SVTYPE', 'SVLEN', 'VAF', 'Minda_IDs']

    if command == "truthset":
        fn_columns = ['#CHROM_x', 'POS_x', 'ID_x', 'INFO_x', \
                            '#CHROM_y', 'POS_y', 'ID_y', 'INFO_y', \
                            'SVTYPE', 'SVLEN', 'VAF', 'Minda_IDs'] 
        
    fn_df = base_df[base_df[f'{caller_name}'] == False][fn_columns]
    tp_df = paired_df[paired_df['Minda_ID'].isin(tp_ids)]
    fp_df = paired_df[~paired_df['Minda_ID'].isin(tp_ids)]
    
    
    return tp_df, fn_df, fp_df, base_df, paired_df
    
def _get_stats_df(tp_df, fn_df, fp_df, paired_df, base_df, caller_name, max_len, out_dir, sample_name, command, vaf, version):
    
    tp = tp_df.shape[0]
    fn = fn_df.shape[0]
    fp = fp_df.shape[0]
    
    # make tsv
    tp_df.to_csv(f'{out_dir}/{sample_name}_{caller_name}_tp.tsv', sep='\t', index=False)
    fn_df.to_csv(f'{out_dir}/{sample_name}_{caller_name}_fn.tsv', sep='\t', index=False)
    fp_df.to_csv(f'{out_dir}/{sample_name}_{caller_name}_fp.tsv', sep='\t', index=False)
    # dfs = [tp_df, fn_df, fp_df]
    # df_names = ["tp", "fn", "fp"]
    # date = datetime.today().strftime('%Y-%m-%d')
    # for i in range(len(dfs)):
    #     df = dfs[i]
    #     df_name = df_names[i]
    #     with open(f'{out_dir}/{sample_name}_{caller_name}_{df_name}.vcf', 'w') as file:
    #         file.write(f'##fileformat=VCFv4.2\n##fileDate={date}\n##source=MindaV{version}\n')
    #         file.write('##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=INS,Description="Insertion">\n##ALT=<ID=DUP,Description="Duplication">\n##ALT=<ID=INV,Description="Inversion">\n')
    #         file.write('##FILTER=<ID=PASS,Description="Default">\n')
    #         file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">\n##INFO=<ID=SUPP_VEC,Number=.,Type=String,Description="IDs of support records">\n')
    #         if vaf != None:
    #             file.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">\n')
    #         command_str = " ".join(sys.argv)
    #         file.write(f"cmd: {command_str}\n")
    #         df.to_csv(file, sep="\t", index=False)

    # caluluate stats
    if tp+fp == 0:
        sys.exit(f"{caller_name} has no TP or FP records. Please double check input files.")
    precision = tp/(tp+fp) 
    if tp+fn == 0:
        sys.exit(f"{caller_name} has no TP or FN records. Please double check input files.")
    recall = tp/(tp+fn)
    f1 = (2*precision*recall)/(precision+recall)

    caller_len = len(paired_df)
    base_len = len(base_df)

    # overall df
    columns = ['True Positives', 'False Negatives', 'False Positives', 'Precision', 'Recall', 'F1 Score', 'Caller Records', 'Ensemble Records']
    data = [[tp, fn, fp, precision, recall, f1, caller_len, base_len]]
    overall_df = pd.DataFrame(data, columns=columns, index=[caller_name])

    # SV type dfs
    tp_type_df = tp_df['SVTYPE_y'].value_counts().to_frame(name=caller_name).rename_axis("SVTYPE").T.sort_index(axis=1)
    fn_type_df = fn_df['SVTYPE'].value_counts().to_frame(name=caller_name).T.sort_index(axis=1) 
    fp_type_df = fp_df['SVTYPE_y'].value_counts().to_frame(name=caller_name).rename_axis("SVTYPE").T.sort_index(axis=1)
    

    # SV len dfs 
    ranges = [ -1, 0, 50, 100, 1000, 10000]#, max_len]
    # ensure bins must increase monotonically
    ranges = [x for x in ranges if x < max_len] + [max_len]
    tp_len_df = tp_df['SVLEN'].value_counts(bins=ranges, sort=False).to_frame(name=caller_name).rename_axis("SVLEN").T 
    fn_len_df = fn_df['SVLEN'].value_counts(bins=ranges, sort=False).to_frame(name=caller_name).T 
    fp_len_df = fp_df['SVLEN'].value_counts(bins=ranges, sort=False).to_frame(name=caller_name).rename_axis("SVLEN").T 
        
    return overall_df, tp_type_df, fn_type_df, fp_type_df, tp_len_df, fn_len_df, fp_len_df


def get_results(decomposed_dfs_list, base_dfs, caller_names, out_dir, sample_name, max_len, tolerance, vaf, command, args, version):
    
    # tp, fn, fp dfs for each caller
    stats_dfs_list = []
    for i in range(len(decomposed_dfs_list)):
        
        decomposed_dfs = decomposed_dfs_list[i]
        caller_name = caller_names[i]
        tp_df, fn_df, fp_df, base_df, paired_df = _get_tp_fn_fp(base_dfs, decomposed_dfs, caller_name, vaf, command)
        stats_dfs = _get_stats_df(tp_df, fn_df, fp_df, paired_df, base_df, caller_name, max_len, out_dir, sample_name, command, vaf, version)            
        stats_dfs_list.append(stats_dfs)

    overall_results_df = pd.concat([df[0] for df in stats_dfs_list])
    tp_type_results_df = pd.concat([df[1] for df in stats_dfs_list]).fillna(0).astype(int)
    fn_type_results_df = pd.concat([df[2] for df in stats_dfs_list]).fillna(0).astype(int)
    fp_type_results_df = pd.concat([df[3] for df in stats_dfs_list]).fillna(0).astype(int)
    tp_len_results_df = pd.concat([df[4] for df in stats_dfs_list]).fillna(0).astype(int)
    fn_len_results_df = pd.concat([df[5] for df in stats_dfs_list]).fillna(0).astype(int)
    fp_len_results_df = pd.concat([df[6] for df in stats_dfs_list]).fillna(0).astype(int)

    results_dfs = [overall_results_df, tp_type_results_df, fn_type_results_df, fp_type_results_df, \
                    tp_len_results_df, fn_len_results_df, fp_len_results_df] 
    
    #headings = ['OVERALL\n\n', '\n\nSV TYPE RESULTS\nTrue Positives\n\n', 'False Negatives\n', 'False Positives\n',\
                #'\n\nSV LENGTH RESULTS\nTrue Positives\n', 'False Negatives\n', 'False Positives\n']
    #user_input = ", ".join([f"{key}={value}" for key, value in vars(args).items() if value is not None and key!= "func"])
    # with open(f'{out_dir}/{sample_name}_minda_results.txt', 'w') as file:
    #     file.write("MINDA ENSEMBLE RESULTS\n\n")
    #     for i in range(len(results_dfs)):
    #         heading = headings[i]
    #         df = results_dfs[i]
    #         file.write(headings[i])
    #         if df.isna().all().all():
    #             file.write("None" + '\n\n')
    #         else:
    #             file.write(df.to_string() + '\n\n')
    #     file.write(f'##minda_args: {user_input}\n')

    file_names = ['overall', 'SV_type_TP', 'SV_type_FN', 'SV_type_FP',\
                  'SV_length_TP', 'SV_length_FN', 'SV_length_FP']
    
    if not os.path.isdir(args.out_dir + "/results"):
        os.makedirs(args.out_dir + "/results")
   
    for i in range(len(results_dfs)):
        file_name = file_names[i]
        df = results_dfs[i].copy()
        df.insert(0, "Caller", caller_names)
        df.to_csv(args.out_dir + f"/results/{file_name}.tsv", sep='\t', index=False)
        

    return overall_results_df, tp_type_results_df, fn_type_results_df, fp_type_results_df, tp_len_results_df, fn_len_results_df, fp_len_results_df, paired_df