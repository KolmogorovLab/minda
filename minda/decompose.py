import os
import logging
import pandas as pd
import numpy as np
from pybedtools import BedTool

logger = logging.getLogger()

def get_caller_name(vcf):
    """
    Extracts the name of the caller from vcf.

    """
    found_source = False
    with open(vcf, 'r') as file:
        for line in file:
            if line.startswith("##source"):
                source_line = line.strip()
                line_length = len(source_line)
                start = source_line.find("=") + 1
                if line_length > 50:
                    stop = source_line.find(" ")
                    if stop == -1:
                        stop = 51
                else:
                    stop = line_length
    
                caller_name = source_line[start:stop]
                found_source = True

    if not found_source:
        caller_name = "Unknown source"
    
    return caller_name

def _get_column_names(vcf):
    """
    Extracts the column names from vcf

    """
    with open(vcf, 'r') as file:
        for line in file:
            if line.startswith("#CHROM"):
                column_names = line.strip().split('\t')
                return column_names
    raise ValueError("No line found starting with #CHROM.")


def get_df(vcf):
    """
    Create a df from vcf.
    
    """
    columns = _get_column_names(vcf)
    df = pd.read_csv(vcf, comment='#', sep='\t', names=columns, dtype={'#CHROM': 'str', 'POS':'Int64'})
    usecols = ['#CHROM', 'POS', 'ID', 'ALT', 'FILTER', 'INFO']
    df = df[usecols]
    return df


def get_intersected_df(vcf, bed):
    """
    Create a df that only includes records that interesect intervals of the bed file.

    """
    bed_to_bt = BedTool(bed)
    columns = _get_column_names(vcf)
    vcf_to_bt = BedTool(vcf)
    intersect_obj = vcf_to_bt.intersect(bed_to_bt, u=True)
    df = BedTool.to_dataframe(intersect_obj, header=None, names=columns, dtype={'#CHROM': 'str', 'POS':'int'})
    return df


def _get_sorted_df(df):
    """
    Sorts dataframe by #CHROM and POS
    
    """
    chrom_value = df.iloc[0,0]
    if chrom_value.startswith("chr"):
        chrom_set = set(df["#CHROM"].str.slice(start=3).to_list())
    else:
        chrom_set = set(df["#CHROM"].to_list())
        
    str_chrom_list = []
    int_chrom_list = []
    for chrom_str in chrom_set:
        try:
            chrom = int(chrom_str)
            int_chrom_list.append(chrom)
        except ValueError:
            chrom = chrom_str
            str_chrom_list.append(chrom)
    if chrom_value.startswith("chr"):
        chrom_sort = sorted(int_chrom_list) + sorted(str_chrom_list)
        chrom_sort = ['chr' + str(chrom_sort[i]) for i in range(len(chrom_sort))]
    else:
        chrom_sort = sorted(int_chrom_list) + sorted(str_chrom_list)
        chrom_sort = [str(chrom_sort[i]) for i in range(len(chrom_sort))]
        
    df = df.sort_values(by=['#CHROM', 'POS'], key=lambda x: x.map({v: i for i, v in enumerate(chrom_sort)})).reset_index(drop=True)
    
    return df


def _get_alt_mate_index(df):
    
    #create mate df & create list of df values
    alt_df = df.ALT.str.extract(r'((chr)?\w+):(\d+)').rename(columns={0: "MATE_#CHROM", 1: "regex_noncap_group ", 2:"MATE_POS"})
    alt_df.MATE_POS = alt_df.MATE_POS.astype(pd.Int64Dtype())
    alt_df.drop(columns="regex_noncap_group ", inplace=True) 
    mate_df = df[['#CHROM','POS']].merge(alt_df, left_index=True, right_index=True)
    mate_df_list = mate_df.values.tolist()

    # find the index of mate record
    mate_indices = []
    for i in range(len(mate_df_list)):
        record = mate_df_list[i]
        chrom = record[0]
        pos = record[1]
        mate_chrom = record[2]
        mate_pos = record[3]

        matching_rows = mate_df[(mate_df['#CHROM'] == mate_chrom) & \
        (mate_df['POS'] == mate_pos) & \
        (mate_df['MATE_#CHROM'] == chrom) & \
        (mate_df['MATE_POS'] == pos)]

        matching_row_indices = matching_rows.index.to_list()
        matching_row_indices = [index for index in matching_row_indices if index != i]
        if len(matching_row_indices) == 1:
            mate_index = matching_row_indices[0] 
        else:
            mate_index = -1
        mate_indices.append(mate_index)

    df["MATE_INDEX"] = mate_indices

    # get mate results
    unique_indices_count = df['MATE_INDEX'].nunique()
    first_unique_index = df["MATE_INDEX"].value_counts().to_frame().index[0]
    unpaired_recrods_count = len(df[df["MATE_INDEX"] == -1])
    paired_records_count = len(df[df["MATE_INDEX"] != -1])

    logger.info(f"Number of unique indices: {unique_indices_count}")
    if unique_indices_count == len(df):
        logger.info(f"{paired_records_count} paired records found...")
    elif unique_indices_count == 1 and first_unique_index == -1:
        logger.info("No paired records found...")
    else:
        logger.info(f"{paired_records_count} paired records and {unpaired_recrods_count} unpaired records found...")

    return df


def _get_paired_alt_dfs(alt_df):
    
    # check if BNDS are a single record or two
    mask = alt_df['MATE_INDEX'] == -1
    singleton_df = alt_df[mask]
    singleton_count = len(singleton_df)
    if singleton_count == alt_df.shape[0]: # ALT records are singletons
        logger.debug(f"(1) Number of singleton records: {singleton_count} {alt_df.shape[0]}")
        alt_df_1 = alt_df.copy()
        alt_df_2 = alt_df.copy()  
        alt_df_2['#CHROM'] = alt_df_2.ALT.str.extract(r'(chr\w+|\w+):')[0].to_list()
        alt_df_2['POS'] = alt_df_2.ALT.str.extract(r':(\d+)')[0].astype(pd.Int64Dtype()).to_list() 
        paired_alt_dfs = [alt_df_1, alt_df_2]
        logger.debug(f"(1) Number of alt/alt_1/alt_2 records: {alt_df.shape[0]} {alt_df_1.shape[0]} {alt_df_2.shape[0]}")
        logger.info(f"Number of paired records paired by ALT column: {alt_df_1.shape[0]} {alt_df_2.shape[0]}")
        logger.info(f"Number of unpaired records paired by MATE_ID: 0 0")
    
    # alt_df pairs based on mate index
    else:
        logger.debug(f"(2) Number of singleton/alt records: {singleton_count} {alt_df.shape[0]}")
        alt_df_1 = alt_df[(alt_df.index < alt_df.MATE_INDEX) & (alt_df.MATE_INDEX != -1)]
        alt_df_2 = alt_df[(alt_df.index > alt_df.MATE_INDEX) & (alt_df.MATE_INDEX != -1)]
        alt_df_2.index = alt_df_2["MATE_INDEX"].to_list()
        paired_alt_dfs = [alt_df_1, alt_df_2]

        mateless_alt_df = alt_df[alt_df.MATE_INDEX == -1]
        logger.debug(f"(2a) Number of alt/alt_1/alt_2/mateless records: {len(alt_df)} {len(alt_df_1)} {len(alt_df_2)} {len(mateless_alt_df)}")
        logger.info(f"Number of paired records paired by ALT column: {alt_df_1.shape[0]} {alt_df_2.shape[0]}")

        # alt_df pairs based on MATEID in INFO
        if mateless_alt_df.shape[0] > 0:
            
            mate_id_alt_df = _get_mate_id_df(mateless_alt_df)
            mate_id_alt_df_1 = mate_id_alt_df[(mate_id_alt_df.index < mate_id_alt_df.MATE_INDEX) & (mate_id_alt_df.MATE_INDEX != -1)]
            mate_id_alt_df_2 = mate_id_alt_df[(mate_id_alt_df.index > mate_id_alt_df.MATE_INDEX) & (mate_id_alt_df.MATE_INDEX != -1)]
            mate_id_alt_df_2.index = mate_id_alt_df_2["MATE_INDEX"].to_list()
            alt_to_info_df = mate_id_alt_df[mate_id_alt_df.MATE_INDEX == -1] # for SEVERUS INS
            logger.debug(f"(2b) Number of mate_id/mate_id/mate_id/alt_to_info records: {mate_id_alt_df.shape[0]} {mate_id_alt_df_1.shape[0]}  {mate_id_alt_df_2.shape[0]} {alt_to_info_df.shape[0]}")
            #alt_df_1 = pd.concat([alt_df_1, mate_id_alt_df_1])
            #alt_df_2 = pd.concat([alt_df_2, mate_id_alt_df_2])
            non_empty_1 = [df for df in [alt_df_1, mate_id_alt_df_1] if not df.empty]
            alt_df_1 = pd.concat(non_empty_1).sort_index()
            non_empty_2 = [df for df in [alt_df_2, mate_id_alt_df_2] if not df.empty]
            alt_df_2 = pd.concat(non_empty_2).sort_index()
            
            paired_alt_dfs = [alt_df_1, alt_df_2]
            logger.debug(f"(2b) Total number of alt_1/alt_2 records: {alt_df_1.shape[0]} {alt_df_2.shape[0]}")
            if alt_to_info_df.shape[0] > 0:
                paired_alt_dfs.append(alt_to_info_df)
            logger.info(f"Number of unpaired records paired by MATE_ID: {len(mate_id_alt_df_1)} {len(mate_id_alt_df_2)}")
        else:
            logger.info(f"Number of unpaired records paired by MATE_ID: 0 0")
    return paired_alt_dfs


def _get_mate_id_df(df):
    """
    Finds the index of the mate ID listed in INFO in the ID column. If not found, index of -1 assigned. 
    
    """
    mate_id_pattern = r'MATEID=([^;]+)(?=;|$)'
    mate_id_list = df.INFO.str.extract(mate_id_pattern)[0].to_list()
    
    mate_indices = []
    for i in range(len(mate_id_list)):
        mate_id = mate_id_list[i]
        matching_rows = df[(df['ID'] == mate_id)]
        matching_row_indices = matching_rows.index
        
        if len(matching_row_indices) == 1:
            mate_index = matching_row_indices[0]    
        else:
            mate_index = -1
        mate_indices.append(mate_index)

    #df['MATE_INDEX'] = mate_indices
    df.loc[:, 'MATE_INDEX'] = mate_indices
        
    return df
    

def _get_paired_info_dfs(info_df):
    
    mate_pos_list = info_df.INFO.str.extract(r'SVLEN=(-?\d+)')[0].astype(pd.Int64Dtype()).abs().to_list()
    info_df_1 = info_df.copy()
    info_df_2 = info_df.copy()
    info_df_2['POS'] = info_df_2['POS'] + mate_pos_list
    info_df_2['END'] = info_df_2.INFO.str.extract(r'END=(-?\d+)')[0].astype(pd.Int64Dtype()).abs().to_list()
    #info_df_2['POS'].fillna(info_df_2['END'], inplace=True)
    info_df_2.fillna({'POS':info_df_2['END']}, inplace=True)
    nan_indices = info_df_2[info_df_2['POS'].isna()].index
    info_df_1 =info_df_1.drop(index=nan_indices, errors='ignore')
    info_df_2 =info_df_2.drop(index=nan_indices, errors='ignore')
    dropped_singleton_count  = info_df.shape[0] - info_df_1.shape[0]
    
    logger.info(f"Number of unpaired records paired by INFO column: {info_df_1.shape[0]} {info_df_2.shape[0]}")
    logger.info(f"Number of singleton records dropped: {dropped_singleton_count}")
    
    return info_df_1, info_df_2


def _check_df_order(df_1, df_2):
    
    # row by row for start and end df, check that the order by sorting
    for i in range(len(df_1)):
        
        df_number = [0]
        row_1 = df_1.iloc[i].to_frame().T
        row_1['row_number'] = df_number
        
        df_number = [1]
        row_2 = df_2.iloc[i].to_frame().T
        row_2['row_number'] = df_number
    
        order_df = pd.concat([row_1, row_2]).reset_index(drop=True)
        sorted_order_df =  _get_sorted_df(order_df)
    
        # if sort is out of order, what the chrom & pos values of the start & end dfs
        if order_df.equals(sorted_order_df) == False: 
            df_1.at[i,'#CHROM'] = sorted_order_df.iloc[0]['#CHROM']
            df_1.at[i, 'POS'] = sorted_order_df.iloc[0]['POS']
            df_2.at[i,'#CHROM'] = sorted_order_df.iloc[1]['#CHROM']
            df_2.at[i, 'POS'] = sorted_order_df.iloc[1]['POS']
                       
    return df_1, df_2


def _write_removed_records(id_set_difference, caller_name, step, written_count, sample_name, out_dir):
    id_difference_list = sorted(list(id_set_difference))
    path = f'{out_dir}/{sample_name}_removed_records.txt'
    file_check = os.path.isfile(path)
    with open(path,'a') as file:
        if file_check == False:
            file.write('REMOVED RECORDS\n')  
        if written_count == 0:
            step =  "***** " + caller_name + " *****" + "\n" + step
        file.write(f'\n{step}\n')
        file.write('\n'.join(id_difference_list))
        file.write('\n')

def get_decomposed_dfs(caller_name, df, filter, min_size, prefixed, vaf, sample_name, out_dir):
    """
    Decomposes df records into start and end dfs.

    """
    logger.info(f"DECOMPOSING {caller_name} RECORDS...")
    logger.info(f"Original number of records: {len(df)}")
    
    written_count = 0
    original_id_set = set(df.ID.to_list())
    
    # filter_df
    if filter != None:
        df = df[df['FILTER'].isin(filter)]
    logger.info(f"Number of after filtering by FILTER column: {len(df)}")

    # write removed ids to txt
    filter_id_set = set(df.ID.to_list())
    id_set_difference = original_id_set.difference(filter_id_set)
    if len(list(id_set_difference)) > 0:
        step = "FILTER"
        _write_removed_records(id_set_difference, caller_name, step, written_count, sample_name, out_dir)
        written_count +=1
    
    # sort df
    df = _get_sorted_df(df)

    # change EVENTTYPE to SVTYPE (for GRIDSS/GRIPSS)
    # create SVTYPE column
    df.loc[:, 'INFO'] = df.INFO.str.replace('EVENTTYPE', 'SVTYPE')
    svtype_pattern = r'SVTYPE=([^;]+)(?=;|$)'
    svtype_column = df.INFO.str.extract(svtype_pattern)[0].to_list()
    df["SVTYPE"] = svtype_column

    # create SVLEN column
    # change SVINSLEN to SVLEN only if SVLEN not in info (for nanomonsv)
    df.loc[~df['INFO'].str.contains('SVLEN', na=False), 'INFO'] = df['INFO'].str.replace('SVINSLEN', 'SVLEN')
    
    # add VAF column
    if vaf != None:
        df['VAF'] = df.INFO.str.extract(r';VAF=([\d.]+)')[0].astype('float').to_list()
    
    # get indices of mate rows
    df = _get_alt_mate_index(df)

    # create paired ALT dfs
    alt_df = df[df['ALT'].str.contains(r'(?:chr)?\w+:\d+', na=False)].copy()

    #create paired INFO dfs
    info_df = df.drop(index=alt_df.index, errors='ignore') 
    logger.debug(f"Number of INFO records: {info_df.shape[0]}")

    # get ALT paired dfs
    paired_alt_dfs = _get_paired_alt_dfs(alt_df)
    if len(paired_alt_dfs) == 3:
        alt_df_1 = paired_alt_dfs[0]
        alt_df_2 = paired_alt_dfs[1]
        alt_to_info_df = paired_alt_dfs[2]

        info_df = pd.concat([info_df, alt_to_info_df])
        logger.debug(f"Total number of INFO records: {info_df.shape[0]}")
    else:
        alt_df_1 = paired_alt_dfs[0]
        alt_df_2 = paired_alt_dfs[1]
        

    info_df_1, info_df_2 = _get_paired_info_dfs(info_df)

    non_empty_1 = [df for df in [info_df_1, alt_df_1] if not df.empty]
    decomposed_df_1 = pd.concat(non_empty_1).sort_index()
    non_empty_2 = [df for df in [info_df_2, alt_df_2] if not df.empty]
    decomposed_df_2 = pd.concat(non_empty_2).sort_index()

    # write removed ids to txt
    singleton_id_set = set(decomposed_df_1.ID.to_list() + decomposed_df_2.ID.to_list())
    id_set_difference = filter_id_set.difference(singleton_id_set)
    if len(id_set_difference) > 0:
        step = "SINGLETON"
        _write_removed_records(id_set_difference, caller_name, step, written_count, sample_name, out_dir)
        written_count += 1

    # check that start and end record are in correct df
    decomposed_df_1, decomposed_df_2 = _check_df_order(decomposed_df_1, decomposed_df_2)

    # write removed ids to txt
    order_id_set = set(decomposed_df_1.ID.to_list() + decomposed_df_2.ID.to_list())
    id_set_difference = singleton_id_set.difference(order_id_set)
    step = "END/START ORDER"
    if len(id_set_difference) > 0:
            step = "END/START ORDER"
            _write_removed_records(id_set_difference, caller_name, step, written_count, sample_name, out_dir)
            written_count += 1

    decomposed_df_1['Minda_ID'] = f'{caller_name}_' + (decomposed_df_1.index + 1).astype(str)
    decomposed_df_2['Minda_ID'] = f'{caller_name}_' + (decomposed_df_2.index + 1).astype(str)
    logger.info(f"Number of decomposed records after pairing: {decomposed_df_1.shape[0]} {decomposed_df_2.shape[0]}")

    # create SVLEN column determined on start & end df (not all vcfs have SVLEN in INFO)
    decomposed_df_1['SVLEN'] = decomposed_df_1.apply(lambda row: -1 if row['#CHROM'] != decomposed_df_2.loc[row.name, '#CHROM'] else int(abs(row['POS'] - decomposed_df_2.loc[row.name, 'POS'])), axis=1)
    max_svlen = decomposed_df_1['SVLEN'].max()
    
    if min_size != None:
        decomposed_df_1 = decomposed_df_1[(decomposed_df_1['SVLEN'] >= min_size) | (decomposed_df_1['SVLEN'] == -1)]
        decomposed_df_2 = decomposed_df_2[decomposed_df_2.index.isin(decomposed_df_1.index)]
        logger.info(f"Number of records after size filtering: {len(decomposed_df_1)} {len(decomposed_df_2)}")

    # write removed ids to txt
    svlen_id_set = set(decomposed_df_1.ID.to_list() + decomposed_df_2.ID.to_list())
    id_set_difference = order_id_set.difference(svlen_id_set)
    if len(id_set_difference) > 0:
        step = "SVLEN"
        _write_removed_records(id_set_difference, caller_name, step, written_count, sample_name, out_dir)
        written_count += 1

    # filter low VAFs such that if either the start or end VAF is too low, records from both dfs are removed
    if vaf != None:
        decomposed_df_1 = decomposed_df_1[decomposed_df_1['VAF'] >= vaf]
        decomposed_df_2 = decomposed_df_2[decomposed_df_2['VAF'] >= vaf]
        minda_ids_list = pd.merge(decomposed_df_1, decomposed_df_2, on='Minda_ID')['Minda_ID'].to_list()
        decomposed_df_1 = decomposed_df_1[decomposed_df_1['Minda_ID'].isin(minda_ids_list)]
        decomposed_df_2 = decomposed_df_2[decomposed_df_2['Minda_ID'].isin(minda_ids_list)]
        logger.info(f"Number of records after VAF filtering: {len(decomposed_df_1)} {len(decomposed_df_2)}")

    # write removed ids to txt
    vaf_id_set = set(decomposed_df_1.ID.to_list() + decomposed_df_2.ID.to_list())
    id_set_difference = svlen_id_set.difference(vaf_id_set)
    step = "VAF"
    if len(id_set_difference) > 0:
        step = "VAF"
        _write_removed_records(id_set_difference, caller_name, step, written_count, sample_name, out_dir)
        written_count += 1
          
    logger.info(f"Total number of decomposed records: {decomposed_df_1.shape[0]} {decomposed_df_2.shape[0]}")

    if prefixed == True:
        prefix = caller_name.split('_', 1)[0]
        decomposed_df_1.ID = prefix + "_" + decomposed_df_1['ID'].astype(str)
        decomposed_df_2.ID = prefix + "_" + decomposed_df_2['ID'].astype(str)

    return decomposed_df_1, decomposed_df_2, max_svlen