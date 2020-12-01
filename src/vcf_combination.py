############################################################################################
# FileName     [ vcf_combination.py ]
# PackageName  [ src ]
# Synopsis     [ Combine 2 or more vcf files for a sample with different variants ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020. 7 ]
############################################################################################

from .vcf_tool import *
import vcf
import os, sys, operator
from tqdm import tqdm, trange
from termcolor import colored
from datetime import datetime
from tabulate import tabulate
import numpy as np

# Combine metadata
def combine_metadata(read_list, caller_list, NT):
    if len(read_list) != len(caller_list):
        print(colored("Error when reading combination files!!", "red"))
        return
    Combine_metadata = read_list[0].metadata.copy()
    Combine_metadata.clear()
    # filedate
    Combine_metadata['filedate'] = []
    Combine_metadata['filedate'].append(datetime.now().strftime("%Y%m%d"))
    # Fileformat
    fileformat =[i.metadata['fileformat'] for i in read_list if 'fileformat' in i.metadata]
    if not any([i != fileformat[0] for i in fileformat]):
        Combine_metadata['fileformat'] = [fileformat[0]]
    # Source
    Combine_metadata['source'] = [caller_list]
    # Source version
    Combine_metadata['source_version'] = [[]]
    for i in range(len(caller_list)):
        if caller_list[i] == "MuSe":
            Combine_metadata['source_version'][0].append(read_list[i].metadata['MuSE_version'][0])
        elif caller_list[i] == "Mutect2":
            Combine_metadata['source_version'][0].append(read_list[i].metadata['MutectVersion'][0])
        elif caller_list[i] == "Strelka2":
            Combine_metadata['source_version'][0].append(read_list[i].metadata['source_version'][0])
        else:
            Combine_metadata['source_version'][0].append('.')
    # Reference
    Combine_metadata['reference'] = []
    for idx, file in enumerate(read_list):
        if 'reference' in file.metadata:
            Combine_metadata['reference'].append(file.metadata['reference'])
        else:
            Combine_metadata['reference'].append('.')
    # NORMAL and TUMOR
    Combine_metadata['NORMAL'] = [NT[0]]
    Combine_metadata['TUMOR'] = [NT[1]]
    return Combine_metadata

# Combine infos, filters, formats
def combine_infos_filters_formats(read_list):
    combine_info = read_list[[len(i.infos) for i in read_list].index(max([len(i.infos) for i in read_list]))].infos.copy()
    combine_filter = read_list[[len(i.filters) for i in read_list].index(max([len(i.filters) for i in read_list]))].filters.copy()
    combine_format = read_list[[len(i.formats) for i in read_list].index(max([len(i.formats) for i in read_list]))].formats.copy()
    for idx, file in enumerate(read_list):
        for info in file.infos:
            if info not in combine_info:
                combine_info[info] = file.infos[info] 
        for flt in file.filters:
            if flt not in combine_filter:
                combine_filter[flt] = file.filters[flt]  
        for fmt in file.formats:
            if fmt not in combine_format:
                combine_format[fmt] = file.formats[fmt]
    combine_info['CALLS'] = vcf.parser._Info(id= "CALLS", num = 1, type= "String", desc="Sign caller names that contain this variant", source = None, version = None)
    combine_info['REJECT'] = vcf.parser._Info(id= "REJECT", num = 1, type= "String", desc="Sign caller names that reject this variant", source = None, version = None)
    return combine_info, combine_filter, combine_format

# Combine the above five functions in one
def generate_header(read_list, caller_list, NT):
    read_header = read_list[0]
    read_header.metadata = combine_metadata(read_list, caller_list, NT)
    read_header.infos, read_header.filters, read_header.formats = combine_infos_filters_formats(read_list)
    return read_header

# Get Somatic list
# A few dictionaries in this list, one dictionary for a file
# If don't want to remove chrX, chrY, set Somatic to False
def get_somatic_list(read_list, caller_list, Somatic = True):
    Somatic_list = []
    chromo_num_list =[str(i) for i in range(1,23)]#+["M","X","Y"]
    for read in read_list:
        new_dict = dict(zip(chromo_num_list, [[] for i in range(len(chromo_num_list))]))
        for record in read:
            if str(record.CHROM) in chromo_num_list:
                new_dict[str(record.CHROM)].append(record)
        Somatic_list.append(new_dict)
    somatic_print = [([i]+[len(file[i]) for file in Somatic_list])for i in chromo_num_list]
    num_list = [sum(somatic_print[ic][1:]) for ic, i in enumerate(somatic_print)]
    somatic_print = [chrom for chrom in somatic_print if any(item != 0 for item in chrom[1:])]
    print(colored("The number of each CHROM in each input vcf file:\n", 'green'))
    print(tabulate(somatic_print, headers=['#']+caller_list, tablefmt='orgtbl'))
    print("\n")
    return Somatic_list, num_list

# Merge similar record in different vcf files, and sort them in a correct order
def merge_and_sort(somatic_list, NEW_FORMAT_STR, INFO_list, caller_list, FORMAT_list, Sample_list, num_list):
    # 22 lists in before list
    # Each of the list is related to 1~22 chromosome respectively.
    list_order = [str(i) for i in range(1,23)]
    before_list, no = [[] for i in range(len(list_order))], len(somatic_list)
    print(colored("Merging....\n",'yellow'))
    for idx, CHROM in enumerate(list_order):
        if num_list[idx] != 0:
            pbar = tqdm(total = num_list[idx], desc = "chr"+CHROM+" ")
            candidate, candidate_POS, empty = [], [], False
            for file in somatic_list:
                if len(file[CHROM]) != 0:
                    candidate.append(file[CHROM][0])
                    candidate_POS.append(file[CHROM][0].POS)
                else:
                    candidate.append(sys.maxsize)
                    candidate_POS.append(sys.maxsize)
            while not empty:
                iMin, Min = min(enumerate([i for i in candidate_POS]), key=operator.itemgetter(1))
                if Min != sys.maxsize:
                    similar = list(set(np.where([(i.POS == Min and i.REF == candidate[iMin].REF and i.ALT == candidate[iMin].ALT) for i in candidate if i != sys.maxsize])[0])-set([iMin])-set([ic for ic, i in enumerate(candidate) if i == sys.maxsize]))
                    if len(similar) != 0:
                        #Deal with similar record
                        somatic_list[iMin][CHROM].pop(0)
                        for item in similar:
                            somatic_list[item][CHROM].pop(0)
                        candidate[iMin].INFO = merge_info(candidate, INFO_list, caller_list, iMin, similar)
                        candidate[iMin].samples = merge_samples(candidate, FORMAT_list, Sample_list, iMin, similar)
                        candidate[iMin].FORMAT = NEW_FORMAT_STR
                        before_list[idx].append(candidate[iMin])
                        candidate[iMin] = somatic_list[iMin][CHROM][0] if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize
                        candidate_POS[iMin] = somatic_list[iMin][CHROM][0].POS if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize
                        for item in similar:
                            candidate[item] = somatic_list[item][CHROM][0] if len(somatic_list[item][CHROM]) > 0 else sys.maxsize
                            candidate_POS[item] = somatic_list[item][CHROM][0].POS if len(somatic_list[item][CHROM]) > 0 else sys.maxsize
                        pbar.update(1+len(similar))
                    else:   # move this chosen record into before_list
                        somatic_list[iMin][CHROM].pop(0)
                        candidate[iMin].INFO = get_new_info(candidate[iMin].INFO, INFO_list, caller_list[iMin], len(candidate[iMin].FILTER) == 0 )
                        candidate[iMin].samples = get_new_samples(candidate[iMin], FORMAT_list, Sample_list)
                        candidate[iMin].FORMAT = NEW_FORMAT_STR
                        before_list[idx].append(candidate[iMin])
                        candidate[iMin] = somatic_list[iMin][CHROM][0] if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize
                        candidate_POS[iMin] = somatic_list[iMin][CHROM][0].POS if len(somatic_list[iMin][CHROM]) > 0 else sys.maxsize   
                        pbar.update(1)   
                else:
                    empty = True
    return before_list

# Part2 in vcf_combination
# change info for each record
def get_new_info(record_info, info_list, caller, PASS):
    for info in info_list:
        if info not in record_info:
            record_info[info] = None
    record_info['CALLS'] = caller
    record_info['REJECT'] = None if PASS else caller
    return record_info

# Part2 in vcf_combination
# change samples' data for each record
def get_new_samples(record, format_list, sample_list):
    origin_format = record.FORMAT.split(':')
    for sample in record.samples:
        new_list= []
        for item in format_list:
            if item not in origin_format:
                new_list.append(None)
            else:
                new_list.append(sample[item])
        sample.data = tuple(new_list)
    return record.samples

def merge_info(candidate, INFO_list, caller_list, iMin, similar):
    # merge similar record to one
    merge_INFO = candidate[iMin].INFO.copy()
    for s in similar:
        for info in candidate[s].INFO:
            if info not in merge_INFO:
                merge_INFO[info] =  candidate[s].INFO[info]   
            else:
                if merge_INFO[info] != candidate[s].INFO[info]:
                    merge_INFO[info] = None
    for info in INFO_list:
        if info not in merge_INFO:
            merge_INFO[info] = None
    c = [iMin]+similar
    merge_INFO['CALLS'] = "_".join(list(np.array(caller_list)[c]))
    reject = "_".join(list(np.array(caller_list)[[i for i in c if len(candidate[i].FILTER) != 0]]))
    merge_INFO['REJECT'] = reject if reject != "" else None
    return merge_INFO

def merge_samples(candidate, FORMAT_list, Sample_list, iMin, similar):
    origin_format = candidate[iMin].FORMAT.split(':')
    similar_format = [candidate[s].FORMAT.split(':') for s in similar]
    merge_Sample = candidate[iMin].samples.copy()
    for idx, sample in enumerate(merge_Sample):
        new_list= []
        for item in FORMAT_list:
            sim = []
            for ids, s in enumerate(similar):
                if item in similar_format[ids]:
                    sim.append(candidate[s].samples[idx][item])
            if len(sim) == 0:
                if item not in origin_format:
                    new_list.append(None)
                else:
                    new_list.append(sample[item])
            else:
                if item in origin_format:
                    sim.append(sample[item])
                    if any([s != sim[0] for s in sim]):
                        new_list.append(None)
                    else:
                        new_list.append(sample[item])
                else:
                    if any([s != sim[0] for s in sim]):
                        new_list.append(None)
                    else:
                        new_list.append(sim[0])
        sample.data = tuple(new_list)
    return merge_Sample

def vcf_combination(sample_content_list, caller_list, output_file):
    Read_list = [read_vcf(x) for x in sample_content_list[2]]
    Read_header = generate_header(Read_list, caller_list, [sample_content_list[0], sample_content_list[1]])
    vcf_writer = vcf.Writer(open(output_file,"w"), Read_header)
    INFO_list = [infos for infos in Read_header.infos]
    FORMAT_list = ['GT', 'AD', 'DP', 'AF']
    FORMAT_list.extend([x for x in Read_header.formats if x not in FORMAT_list])
    Sample_list = [sample_content_list[0], sample_content_list[1]]
    NEW_FORMAT_STR = ":".join(FORMAT_list)
    Somatic_list, num_list = get_somatic_list(Read_list, caller_list)
    Merged_list = merge_and_sort(Somatic_list, NEW_FORMAT_STR, INFO_list, caller_list, FORMAT_list, Sample_list, num_list)
    print(colored("\nWriting to ...."+output_file+"\n",'yellow'))
    for idx, chrom in enumerate(Merged_list):
        for record in chrom:
            vcf_writer.write_record(record)
    vcf_writer.close()
    print(colored(("=> Finish combining file: "+output_file+'\n'), 'green'))
