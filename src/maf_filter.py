############################################################################################
# FileName     [ maf_filter.py ]
# PackageName  [ src ]
# Synopsis     [ Sifting data from MAF files ]
# Author       [ Cheng-Hua (Anita) Lu ]
# Copyright    [ 2020 7]
############################################################################################

import csv
from termcolor import colored
import pandas as pd
from tqdm import tqdm, trange
import os

# Read src/auxiliary_file/rna_tissue_consensus.tsv
def read_TSV(tsv_file):
    print(colored("Reading TSV auxiliary file....\nIt may take a while....\n", "yellow"))
    ALL_DICT = {}
    df = pd.read_csv(tsv_file, sep="\t", header = 0)
    pbar = tqdm(total = df.shape[0])
    for i in range(df.shape[0]):
        pbar.update(1)
        if df.iloc[i]['Gene name'] not in ALL_DICT:
            ALL_DICT[df.iloc[i]['Gene name']] = {}
        ALL_DICT[df.iloc[i]['Gene name']][df.iloc[i]['Tissue']] = df.iloc[i]['NX']
    pbar.close()
    return ALL_DICT

# Read MAF(Mutation Annotation Format) file
# Output: head:"" or "#version...." and dataframe
def fast_read_maf(maf_file):
    file = open(maf_file, 'r')
    head = file.readline()
    if head.startswith('#'):
        df = pd.read_csv(maf_file, sep="\t", skiprows = 1, header = 0,encoding = 'gb2312', low_memory=False)
    else:
        df = pd.read_csv(maf_file, sep="\t", header = 0, encoding = 'gb2312', low_memory=False)
        head = ""
    return head, df


# GI: Genome Interval
# Type of "inteval": List or Dictionary
# Examples:
## List: [1,2,3,4] or [1, 3, 5:7]
## Dict: {"1":[1, 1000000], "3":[30000, 20000000]}
def Genome_Interval(data, interval):
    if interval == False:
        return True
    if isinstance(interval, list):  # Only select Genome Number
        if data['Chromosome'] in interval:
            return True
        else:
            return False
    elif isinstance(interval, dict):
        chromo = str(data['Chromosome'])
        start, end = data['Start_Position'], data['End_Position']
        if chromo in interval:
            i_start, i_end = interval[chromo][0], interval[chromo][1]
            if (start >= i_start and start <= i_end) or (end >= i_start and end <= i_end):
                return True
            else:
                return False
        else:
            return False
    else:
        return False
    

# CI: Caller Information
# Type of "info": List[DP_N,DP_T,AD_T]
# Examples:
## List: [15,15,6]
def Call_Information(data, info):
    if info == False:
        return True
    DP_N, DP_T, AD_T = info[0], info[1], info[2]
    if data['t_depth'] > DP_T and data['t_alt_count'] >= AD_T and data['n_depth'] > DP_N and data['n_alt_count'] == '.':
        return True
    else:
        return False
    

# TE: Tissue Expression
# Type of "tissue": List
# Examples:
## List: ['stomach','B-cells','thymus']
def Tissue_Expression(data, tissue, ALL_DICT):
    if tissue == False:
        return True
    Hugo = data['Hugo_Symbol']
    PASS = True
    if Hugo in ALL_DICT:
        for item in tissue:
            if item in ALL_DICT[Hugo]:
                if float(ALL_DICT[Hugo][item]) <= 0:
                    PASS = False
            else:
                PASS = False
    else:
        PASS = False
    return PASS


# PF: Population Frequency
# Type of "info": Bool
# Examples:
## info = 1 or 0
def Population_Frequency(data, info):
    if info == False:
        return True
    return data['FILTER'] == 'PASS'
    
# H: Hypermutator or sample exclusion
# Type of "COUNT": Integer
# Examples:
## COUNT = 500
def Hypermutator(df, COUNT):
    sample_dict = {}
    rm_list = []
    for i in range(df.shape[0]):
        data = df.iloc[i]
        if data['Tumor_Sample_Barcode'] not in sample_dict:
            sample_dict[data['Tumor_Sample_Barcode']] = [i]
        else:
            sample_dict[data['Tumor_Sample_Barcode']].append(i)
    for item in sample_dict:
        if len(sample_dict[item]) > COUNT:
            rm_list += sample_dict[item]
    
    df = df.drop(df.index([rm_list]))
    return df

# MAF filtering
# Write a new MAF file to output_file
def maf_filter(maf_file, flt_list, ALL_DICT, output_file):
    head, df = fast_read_maf(maf_file)
    rm_list = []
    pbar = tqdm(total = df.shape[0])
    for i in range(df.shape[0]):
        pbar.update(1)
        GI = Genome_Interval(df.iloc[i], flt_list[0])
        CI = Call_Information(df.iloc[i], flt_list[1])
        TE = Tissue_Expression(df.iloc[i], flt_list[2], ALL_DICT)
        PF = Population_Frequency(df.iloc[i], flt_list[3])
        if not (GI and CI and TE and PF):
            rm_list.append(i)
    pbar.close()
    df = df.drop(df.index[rm_list])
    if flt_list[4] != False:
        df  = Hypermutator(df, flt_list[4])
    df.to_csv(output_file, sep="\t", index= False, header = list(df.columns.values))
    if head != "":
        with open(output_file, "r+") as filtered_file:
            lines = filtered_file.readlines()
            filtered_file.seek(0)
            filtered_file.write(head)
            filtered_file.writelines(lines)         
