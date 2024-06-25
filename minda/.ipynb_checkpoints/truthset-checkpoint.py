import pandas as pd
import numpy as np

def get_base_df(decomposed_dfs_list, tolerance, multimatch):
    dfs_1 = [dfs_list[0] for dfs_list in decomposed_dfs_list]
    dfs_2 = [dfs_list[1] for dfs_list in decomposed_dfs_list]

    # create collective comparison df
    start_dfs = pd.concat(dfs_1)
    end_dfs = pd.concat(dfs_2)
    
    # create base df
    base_1_df = dfs_1[-1]
    base_2_df = dfs_2[-1]

    # find which comparison start loci are within tolerance range of base start loci
    base_1_loci = list(zip(base_1_df['#CHROM'], base_1_df['POS']))
    start_loci = list(zip(start_dfs['#CHROM'], start_dfs['POS']))
    base_2_loci = list(zip(base_2_df['#CHROM'], base_2_df['POS']))
    end_loci = list(zip(end_dfs['#CHROM'], end_dfs['POS']))
    
    start_index_lists = []
    for i in range(len(base_1_loci)):
        base_locus = base_1_loci[i]
        index_list = []
        
        for j in range(len(start_loci)): # in order to get the correct index cannot use "for start_locus in start_loci"
            
            start_locus = start_loci[j]
            if base_locus[0] == start_locus[0]:
                
                distance = abs(base_locus[1] - start_locus[1])
                if distance <= tolerance:
                    start_index = j
                    index_list.append(start_index)
        
        start_index_lists.append(index_list)
        if len(start_index_lists) != (i+1): # ensure each base record has a list even if no comp calls within tolerance range
            start_index_lists.append([])

    # if start loci within tolerance range, check that end also is
    running_list = []
    comp_minda_ids = start_dfs.Minda_ID.to_list()
    minda_id_lists = []
    for i in range(len(start_index_lists)):
        index_list = start_index_lists[i]
        base_locus = base_2_loci[i]

        
        minda_id_list = []
        for index in index_list:
            end_locus = end_loci[index]
            if base_locus[0] == end_locus[0]:
                #print(base_locus, end_locus)
                distance = abs(base_locus[1] - end_locus[1])
                if distance <= tolerance:
                    minda_id = comp_minda_ids[index]
                    
                    if multimatch == False:
                        caller = minda_id.rsplit('_', 1)[0]
                        if any(id.startswith(caller) for id in minda_id_list) == False and minda_id not in running_list:   
                            minda_id_list.append(minda_id)
                            running_list.append(minda_id)
                    else:
                        minda_id_list.append(minda_id)
        minda_id_lists.append(minda_id_list)
        if len( minda_id_lists) != (i+1): # ensure each base record has a list even if no comp calls within tolerance range
             minda_id_lists.append([])

    # merge start & end base dfs & create column of Minda IDs for calls within tolerance range
    base_df = base_1_df.merge(base_2_df, left_index=True, right_index=True)
    base_df["Minda_IDs"] = minda_id_lists
    
    return base_df


def get_support_df(base_df, caller_names, vaf, out_dir, sample_name):
    
    minda_id_lists = base_df.Minda_IDs.to_list()
    # create call columns for each caller
    for caller_name in caller_names:
        caller_column = []
        for minda_id_list in minda_id_lists:
            
            call_boolean  = any(value.startswith(caller_name) for value in minda_id_list)
            caller_column.append(call_boolean)
        base_df[f'{caller_name}'] = caller_column
    
    # if vaf == None:
    #     base_df['VAF_x'] = np.nan
    
    column_names = ['#CHROM_x', 'POS_x', 'ID_x', 'INFO_x', \
                    '#CHROM_y', 'POS_y', 'ID_y', 'INFO_y', \
                    'SVTYPE_x', 'SVLEN', 'VAF_x', 'Minda_ID_x','Minda_IDs'] + [caller_names[-1]] + caller_names[:-1]
        
    support_df = base_df[column_names].rename(columns={'SVTYPE_x':'SVTYPE', 'VAF_x':'VAF', 'Minda_ID_x': 'Minda_ID'}).copy()

    support_df.to_csv(f'{out_dir}/{sample_name}_support.tsv', sep='\t', index=False)
    
    return support_df

