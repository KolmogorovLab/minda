import pandas as pd 

def _get_tp_fn_fp(support_df, decomposed_dfs, caller_name, vaf, command):

    paired_df = decomposed_dfs[0].merge(decomposed_dfs[1], on='Minda_ID')
    #base_df = support_df[support_df.iloc[:,-1] == True]
    if command == "ensemble":
        base_df = support_df[support_df['ensemble'] == True]
    else:
        base_df = support_df[support_df.iloc[:,-1] == True]
    tp_ids = [id for ids in base_df["Minda_IDs"].to_list() for id in ids]

    # create tp, fn, fp dfs
    if command == "ensemble":
        if vaf != None:
            fn_columns = ['#CHROM_x', 'POS_x', 'locus_group_x', 'ID_list_x', \
                          '#CHROM_y', 'POS_y', 'locus_group_y', 'ID_list_y', \
                          'SVTYPE', 'SVLEN', 'VAF', 'Minda_IDs']
        else:
            fn_columns = ['#CHROM_x', 'POS_x', 'locus_group_x', 'ID_list_x', \
                          '#CHROM_y', 'POS_y', 'locus_group_y', 'ID_list_y', \
                          'SVTYPE', 'SVLEN', 'Minda_IDs']
    if command == "truthset":
        if vaf != None:
            fn_columns = ['#CHROM_x', 'POS_x', 'ID_x', 'INFO_x', \
                            '#CHROM_y', 'POS_y', 'ID_y', 'INFO_y', \
                            'SVTYPE', 'SVLEN', 'VAF', 'Minda_IDs']
        
        else:
            fn_columns = ['#CHROM_x', 'POS_x', 'ID_x', 'INFO_x', \
                            '#CHROM_y', 'POS_y', 'ID_y', 'INFO_y', \
                            'SVTYPE', 'SVLEN', 'Minda_IDs'] 
        
    fn_df = base_df[base_df[f'{caller_name}'] == False][fn_columns]
    tp_df = paired_df[paired_df['Minda_ID'].isin(tp_ids)]
    fp_df = paired_df[~paired_df['Minda_ID'].isin(tp_ids)]
    
    return tp_df, fn_df, fp_df, base_df, paired_df
    
def _get_stats_df(tp_df, fn_df, fp_df, paired_df, base_df, caller_name, max_len, out_dir, sample_name):

    # overall calls dfs
    tp = tp_df.shape[0]
    fn = fn_df.shape[0]
    fp = fp_df.shape[0]
    
    # make csvs
    tp_df.to_csv(f'{out_dir}/{sample_name}_{caller_name}_tp.tsv', sep='\t', index=False)
    fn_df.to_csv(f'{out_dir}/{sample_name}_{caller_name}_fn.tsv', sep='\t', index=False)
    fp_df.to_csv(f'{out_dir}/{sample_name}_{caller_name}_fp.tsv', sep='\t', index=False)

    # caluluate stats
    precision = tp/(tp+fp) 
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


def get_results(decomposed_dfs_list, base_dfs, caller_names, out_dir, sample_name, max_len, tolerance, vaf, command, args):
    
    # tp, fn, fp dfs for each caller
    stats_dfs_list = []
    for i in range(len(decomposed_dfs_list)):
        
        decomposed_dfs = decomposed_dfs_list[i]
        caller_name = caller_names[i]
        tp_df, fn_df, fp_df, base_df, paired_df = _get_tp_fn_fp(base_dfs, decomposed_dfs, caller_name, vaf, command)
        stats_dfs = _get_stats_df(tp_df, fn_df, fp_df, paired_df, base_df, caller_name, max_len, out_dir, sample_name)            
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
    headings = ['OVERALL\n\n', '\n\nSV TYPE RESULTS\nTrue Positives\n\n', 'False Negatives\n', 'False Positives\n',\
                '\n\nSV LENGTH RESULTS\nTrue Positives\n', 'False Negatives\n', 'False Positives\n']
    user_input = ", ".join([f"{key}={value}" for key, value in vars(args).items() if value is not None and key!= "func"])
    
    with open(f'{out_dir}/{sample_name}_minda_results.txt', 'w') as file:
        file.write("MINDA ENSEMBLE RESULTS\n\n")
        for i in range(len(results_dfs)):
            heading = headings[i]
            df = results_dfs[i]
            file.write(headings[i])
            file.write(df.to_string() + '\n\n')
        file.write(f'##minda_args: {user_input}\n')

    return overall_results_df, tp_type_results_df, fn_type_results_df, fp_type_results_df, tp_len_results_df, fn_len_results_df, fp_len_results_df, paired_df