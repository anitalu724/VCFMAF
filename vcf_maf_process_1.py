############################################################################################
# FileName     [ vcf_maf.py ]
# PackageName  [ My_Lab_Tool]
# Synopsis     [ Control all functions ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020/7 ]
############################################################################################

from src.vcf_filter import *
from src.vcf_combination import vcf_combination
from src.vcf2maf import vcf2vep2maf
from src.maf_filter import *
from src.maf_filter import fast_read_maf
from termcolor import colored
from tabulate import tabulate
import numpy as np
import pandas as pd
import argparse, textwrap
import sys, os, json, time
from tqdm import tqdm, trange
import shutil
from enum import Enum

class ErrorPrint(Enum):
    TSV_format = "\nERROR: TSV file's format is wrong!\nEnd of implementation!\n"
    VCF_command = "\nERROR: Command -c, -v2m must required.\nEnd of implementation!\n"
    MAF_command = "\nERROR: Command -vf, -c, -v2m must not required.\nEnd of implementation!\n"
    FILTER_command = "\nERROR: Wrong format for command -vf or -mf!\n"

# Read args.file to determine whether implementing VCF or MAF
def read_tsv(tsv_file):
    file = pd.read_csv(tsv_file, sep='\t')
    category, category_print, category_caller = [], [], []
    flag = ""
    if file.shape[1] == 1:      # MAF only
        for i in range(file.shape[0]):
            if file.loc[i, 'MAF'][-3:]!= "maf":
                print(colored(ErrorPrint.TSV_format.value, 'red'))
                return False, False, False
            category.append(file.loc[i, 'MAF'])
        length = [[str(len(category))]]
        print(colored("\nThe input tsv file\'s content:\n", "green"))
        print(tabulate(length, headers=['# of MAF files'], tablefmt='orgtbl'))
        print("\n")
        flag = "maf"
        return flag, category, category_caller
    elif file.shape[1] == 9:    # VCF and MAF
        for i in range(file.shape[0]):
            sample, vcf_file, caller_list = [], [], []
            sample.append(file.loc[i,'NORMAL'])
            sample.append(file.loc[i,'TUMOR'])
            variant_list = ['MuSe', 'Mutect2', 'SomaticSniper', 'Strelka2', 'VarScan2']
            for j in variant_list:
                if not isinstance(file.loc[i, j], float):
                    if file.loc[i, j][-3:] != "vcf":
                        print(colored(ErrorPrint.TSV_format.value, 'red'))
                        return False, False, False  
                    vcf_file.append(file.loc[i, j])
                    caller_list.append(j)
            sample.append(vcf_file)
            sample.append(file.loc[i,'At Least # CALLS'])
            sample.append(file.loc[i,'At Most # REJECT'])
            category.append(sample)
            category_caller.append(caller_list)
        for i in category:
            if i[3] == 0 or i[3] > len(i[2]):
                print(colored(ErrorPrint.TSV_format.value, 'red'))
                return False, False, False
            sample = []
            sample.extend((i[0], i[1], len(i[2]), i[3], i[4]))
            category_print.append(sample)
        print(colored("\nThe input tsv file\'s content:\n", "green"))
        print(tabulate(category_print, headers=['NORMAL', 'TUMOR', '# of VCF files', 'At least # of variants', 'At most # of REJECT'], tablefmt='orgtbl'))
        flag = "vcf"
        return flag, category, category_caller
    else:
        print(colored(ErrorPrint.TSV_format.value, 'red'))
        return False, False, False

def get_filter_data(vcf_flt, num=4):
    flt_list = [False]*num
    if len(vcf_flt)%2 != 0:
        print(colored(ErrorPrint.FILTER_command.value, 'red'))
    for idx, data in enumerate(vcf_flt):
        if idx%2 == 0:
            if data == "GI":
                info = vcf_flt[idx+1]
                if "{" not in info and "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    new_info = []
                    for i in info:
                        if ":" not in i:
                            new_info.append(int(i))
                        else:
                            s, e = i.split(":")
                            for no in range(int(s),int(e)+1):
                                new_info.append(no)
                    flt_list[0] = new_info
                elif "{" in info:
                    info = info.strip('{').strip('}').split("],")
                    new_dict = {}
                    for i in info:
                        i = i.split(":")
                        s, e = i[1].split(",")
                        s, e = int(s.strip('[')), int(e.strip(']'))
                        new_dict[str(i[0])] = [s, e]
                    flt_list[0] = new_dict
            elif data == "CI":
                flt_list[1] = [float(x) for x in (vcf_flt[idx+1].replace(" ", "").replace("[", "").replace("]", "")).split(",")]
            elif data == "P":
                flt_list[2] = bool(int(vcf_flt[idx+1]))
            elif data == "FFPE":
                flt_list[3] = float(vcf_flt[idx+1])
    return flt_list

def get_maf_filter_data(maf_flt, num = 5):
    flt_list = [False]*num
    for idx, data in enumerate(maf_flt):
        if idx%2 == 0:
            if data == "GI":
                info = maf_flt[idx+1]
                if "{" not in info and "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    new_info = []
                    for i in info:
                        if ":" not in i:
                            new_info.append(int(i))
                        else:
                            s, e = i.split(":")
                            for no in range(int(s),int(e)+1):
                                new_info.append(no)
                    flt_list[0] = new_info
                elif "{" in info:
                    info = info.strip('{').strip('}').split("],")
                    new_dict = {}
                    for i in info:
                        i = i.split(":")
                        s, e = i[1].split(",")
                        s, e = int(s.strip('[')), int(e.strip(']'))
                        new_dict[str(i[0])] = [s, e]
                    flt_list[0] = new_dict
            elif data == "CI":
                info = maf_flt[idx+1]
                if "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    new_info = []
                    for i in info:
                        new_info.append(int(i))
                    flt_list[1] = new_info
            elif data == "TE":
                info = maf_flt[idx+1]
                if "[" in info:
                    info = info.strip("[").strip("]").split(",")
                    flt_list[2] = info
            elif data == "PF":
                flt_list[3] = bool(int(maf_flt[idx+1]))
            elif data == "H":
                flt_list[4] = int(maf_flt[idx+1])
    return flt_list

def nt_test(vcf_read, col_number, caller):
    col_0, col_1 = [], []
    if caller == "SomaticSniper" or caller == "MuSe":
        for record in vcf_read:
            col_0.append(record.samples[0]['SS'])
            col_1.append(record.samples[1]['SS'])
        if all((i == 2) for i in col_0) and not all ((i == 2) for i in col_1):
            return True
        else:
            return False
    elif caller == "Strelka2":
        return False
    else:
        for record in vcf_read:
            col_0.append(record.samples[0]['GT'])
            col_1.append(record.samples[1]['GT'])
        if all((i == "0/0" or i == "0|0") for i in col_1) and not all ((i == "0/0" or i == "0|0") for i in col_0):
            return True
        else:
            return False


def main():
    parser = argparse.ArgumentParser(description='VCF and MAF files processing', \
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f","--file", help="Input the tsv file.\n\n", required = True)
    parser.add_argument("-vf", "--vcf_filter",nargs='*', \
                                          help=textwrap.dedent("GI: Genome Interval\n"
                                                               "CI: Caller Information\n"
                                                               "P: Keep or exclude non-PASS tag\n"
                                                               "FFPE: Artifact variant filter: FFPE filter (Only for caller = Mutect2)\n\n"))
    parser.add_argument("-c","--combine", action='store_true',help="Input the path for combination VCF files. The path must end with a folder.\n\n")
    parser.add_argument("-v2m","--vcf2maf", nargs=1, help="The first value is the path for transformed MAF files. This path must end with a folder.\n"
                                                             "The second value is \"max_filter_ac\" which is an important parameter(int) when transforming files to MAF.\n\n")
    parser.add_argument("-mf", "--maf_filter", nargs="*", \
                        help = textwrap.dedent("GI: Genome Interval\n"
                                               "CI: Caller Information\n"
                                               "TE: Tissue Expression\n"
                                               "PF: Population Frequency\n"
                                               "H: Hypermutation or Sample Exclusion\n\n"))
    parser.add_argument("-o", "--output", required = True, help="The path for storing output files.\nThis path must end with a folder.\n")
    parser.add_argument("-m","--metafile", required = True, help="The path for storing metafiles.\nThis path must end with a folder.\n")
    args = parser.parse_args()

    # Read TSV file
    print(colored("\nReading TSV file....", "yellow"))
    flag, category, category_caller = read_tsv(args.file)
    maf_output_list = []
    folder, meta = args.output, args.metafile
    if folder[-1:] != "/":
        folder += "/"
    if meta[-1:] != "/":
        meta += "/"

# Solution 1: VCF
    if flag == "vcf":  
        if not args.combine or not args.vcf2maf:
            print(colored(ErrorPrint.VCF_command.value, "red"))
            return
    # # VCF Formalizing (Required) and Filtering (Optional)
    #     filter_list = []
    #     if args.vcf_filter:
    #         print(colored("\n\nFormalizing and filtering VCF files....\n", "yellow"))
    #         filter_list = get_filter_data(args.vcf_filter, 4)
    #     else:
    #         print(colored("\n\nFormalizing VCF files....\n", "yellow"))
    #     for s_idx, sample in enumerate(category):
    #         formalized_data_list = []
    #         vcf_read_set = [read_vcf(x) for x in sample[2]]
    #         for vcf_idx, vcf_read in enumerate(vcf_read_set):
    #             new_filter_list = list(filter_list)
    #             if args.vcf_filter:
    #                 FFPE = (category_caller[s_idx][vcf_idx] != "Mutect2")
    #                 if FFPE and new_filter_list[3]:
    #                     new_filter_list[3] = None
    #                     print(colored(("Warning: FFPE filter does not apply to variant = "+ category_caller[s_idx][vcf_idx]), 'yellow'))
    #                     print("Skip the FFPE filter for " + sample[2][vcf_idx])
    #             del_count, change = 0, False
    #             output_path = meta+(sample[2][vcf_idx].split("/")[-1])[:-4]+"_formalized.vcf" if not args.vcf_filter else meta+(sample[2][vcf_idx].split("/")[-1])[:-4]+"_formalized_filter.vcf"
    #             formalized_data_list.append(output_path)
    #             if len(vcf_read.samples) == 2:
    #                 change = nt_test(vcf_read, len(vcf_read.samples), category_caller[s_idx][vcf_idx])
    #             vcf_read = read_vcf(sample[2][vcf_idx])
    #             vcf_read.samples = [sample[0], sample[1]] if len(vcf_read.samples) == 2 else [sample[1]]
    #             vcf_writer = vcf.Writer(open(output_path,"w"), vcf_read)
    #             for record in vcf_read:
    #                 PASS = True
    #                 if change and PASS:
    #                     x = record.samples[0]
    #                     record.samples[0] = record.samples[1]
    #                     record.samples[1] = x
    #                 if "chr" in record.CHROM and PASS:
    #                     record.CHROM = record.CHROM[3:]
    #                 if len(record.ALT) != 1 and PASS:
    #                     del_count += 1
    #                     PASS = False
    #                 if args.vcf_filter and PASS:
    #                     filter_score = np.zeros(4, dtype = bool)
    #                     for i in range(len(new_filter_list)):
    #                         if new_filter_list[i]:
    #                             score = True
    #                             if i == 0 and type(new_filter_list[i]) is list: # 1-1
    #                                 score = chromosome_select(new_filter_list[i], record)
    #                             elif i == 0 and type(new_filter_list[i]) is dict: # 1-2
    #                                 score = interval_select(new_filter_list[i], record)
    #                             elif i == 1:
    #                                 call = category_caller[s_idx][vcf_idx]
    #                                 [DP_N, DP_T, AD_N, AD_T, AF_N, AF_T, NLOD, TLOD] = new_filter_list[i]
    #                                 score = caller_info(call, record, DP_N, DP_T, AD_N, AD_T, AF_N, AF_T, NLOD, TLOD)
    #                             elif i == 2:
    #                                 score = find_pass(record)
    #                             elif i == 3:
    #                                 call = category_caller[s_idx][vcf_idx]
    #                                 score = FFPE_filter(call, record, new_filter_list[i])
    #                             filter_score[i] = score
    #                         else:
    #                             filter_score[i] = True
    #                     if not all(x == True for x in filter_score):
    #                         PASS = False
    #                         del_count += 1
    #                 if PASS:
    #                     vcf_writer.write_record(record)
    #             if del_count!= 0:
    #                 print("NOTICE: "+str(del_count)+" data have been removed from "+sample[2][vcf_idx])
    #             vcf_writer.close()
    #         listToStr = ', '.join([str(elem) for elem in formalized_data_list]) 
    #         sample[2] = formalized_data_list
    #         print(colored("\n=> Finish formalizing files: \n[ "+ listToStr+" ]\n", "green"))
    
    # # VCF Combination (Required)
    #     if args.combine:
    #         print(colored("\nStart VCF combination....\n", "yellow"))
    #         combine_output_list = []
    #         for idx in range(len(category)):
    #             path_name = meta+category[idx][0]+"_"+category[idx][1]+"_combination.vcf"
    #             combine_output_list.append(path_name)
    #         for idx, sample in enumerate(category):
    #             vcf_combination(category[idx], category_caller[idx], combine_output_list[idx])

    #         ## Deal with "At Least...." condition
    #         combine_filter_filelist = [file[:-4]+"_f.vcf" for file in combine_output_list]
    #         for idx, file in enumerate(combine_output_list):
    #             del_cou = 0
    #             vcf_read = read_vcf(file)
    #             vcf_writer = vcf.Writer(open(combine_filter_filelist[idx],"w"), vcf_read)
    #             for record in vcf_read:
    #                 calls = len(record.INFO['CALLS'].split('_')) if record.INFO['CALLS']!= None else 0
    #                 reject = len(record.INFO['REJECT'].split('_')) if record.INFO['REJECT'] != None else 0
    #                 if calls >= category[idx][3] and reject <= category[idx][4]:
    #                     vcf_writer.write_record(record)
    #                 else:
    #                     del_cou+=1
    #             vcf_writer.close()
    #             print("NOTICE: "+str(del_cou)+" data have been removed from "+combine_output_list[idx])
    #         print("\n")
    
    # VCF2MAF (Required)
        if args.vcf2maf:
            combine_filter_filelist = ['examples/ms_vcf/colon1-sample.vcf',
                                        'examples/ms_vcf/colon2-sample.vcf',
                                        'examples/ms_vcf/colon3-sample.vcf',
                                        'examples/ms_vcf/liver1-sample.vcf',
                                        'examples/ms_vcf/liver2-sample.vcf',
                                        'examples/ms_vcf/liver3-sample.vcf',
                                        'examples/ms_vcf/intestine1-sample.vcf',
                                        'examples/ms_vcf/intestine2-sample.vcf',
                                        'examples/ms_vcf/intestine3-sample.vcf']

            

            print(colored("Start transforming VCF to MAF....\n", "yellow"))
            print("WARNING: This transformation tool must be implemented in the direction of \"mskcc-vcf2maf-bbe39fe\"!\n")
            for idx, file in enumerate(combine_filter_filelist):
                fileName = file[:-4]+"_2maf.maf"
                maf_output_list.append(fileName)
            vcf2vep2maf(combine_filter_filelist, maf_output_list, folder, category, args.vcf2maf[0])

    # MAF Filtering (Optional)
        if args.maf_filter:
            print(colored("Start MAF filtering....\n", "yellow"))
            maf_flt_list = get_maf_filter_data(args.maf_filter)
            maf_filtered_list = [x[:-4]+"_filtered.maf" for x in maf_output_list]
            ALL_DICT = {}
            if maf_flt_list[2] != False:
                if not os.path.isfile('src/auxiliary_file/rna_tissue_consensus.json'):
                    ALL_DICT = read_TSV("src/auxiliary_file/rna_tissue_consensus.tsv")
                    with open("src/auxiliary_file/rna_tissue_consensus.json", "w") as jsonfile:  
                        json.dump(ALL_DICT, jsonfile) 
                else:
                    with open("src/auxiliary_file/rna_tissue_consensus.json", "r") as jsonfile:  
                        ALL_DICT = json.load(jsonfile)
            for idx, maf in enumerate(maf_output_list):
                maf_filter(maf, maf_flt_list, ALL_DICT, maf_filtered_list[idx])
                print(colored(("\n=> Finish filtering MAF file: "+maf_filtered_list[idx]+"\n"), 'green'))
            
            # MAF combination
            if len(maf_filtered_list) > 1:
                print(colored("Start MAF combination....\n", "yellow"))
                maf_df, head = pd.DataFrame(), ""
                for maf_file in maf_filtered_list:
                    head, maf = fast_read_maf(maf_file)
                    maf_df = pd.concat([maf_df, maf])
                maf_df.to_csv(folder+"maf_combination.maf", sep="\t", index= False, header = list(maf_df.columns.values))
                if head != "":
                    with open(folder+"maf_combination.maf", "r+") as filtered_file:
                        lines = filtered_file.readlines()
                        filtered_file.seek(0)
                        filtered_file.write(head)
                        filtered_file.writelines(lines)
                print(colored(("=> Finish combining MAF files to "+folder+"maf_combination.maf"+"\n"), 'green'))
        else:
            if len(maf_output_list) > 1:
                print(colored("Start MAF combination....\n", "yellow"))
                maf_df, head = pd.DataFrame(), ""
                for maf_file in maf_output_list:
                    head, maf = fast_read_maf(maf_file)
                    maf_df = pd.concat([maf_df, maf])
                maf_df.to_csv(folder+"maf_combination.maf", sep="\t", index= False, header = list(maf_df.columns.values))
                if head != "":
                    with open(folder+"maf_combination.maf", "r+") as filtered_file:
                        lines = filtered_file.readlines()
                        filtered_file.seek(0)
                        filtered_file.write(head)
                        filtered_file.writelines(lines)
                print(colored(("=> Finish combining MAF files to " + folder + "maf_combination.maf" + "\n"), 'green'))

    elif flag == "maf":
        if args.vcf_filter or args.combine or args.vcf2maf:
            print(colored(ErrorPrint.MAF_command.value, 'red'))
            return

    # MAF Filtering (Optional)
        if args.maf_filter:
            print(colored("Start MAF filtering....\n", "yellow"))
            maf_output_list = list(category)
            maf_flt_list = get_maf_filter_data(args.maf_filter)
            maf_filtered_list = [x[:-4]+"_filtered.maf" for x in maf_output_list]
            ALL_DICT = {}
            if maf_flt_list[2] != False:
                if not os.path.isfile('src/auxiliary_file/rna_tissue_consensus.json'):
                    ALL_DICT = read_TSV("src/auxiliary_file/rna_tissue_consensus.tsv")
                    with open("src/auxiliary_file/rna_tissue_consensus.json", "w") as jsonfile:  
                        json.dump(ALL_DICT, jsonfile) 
                else:
                    with open("src/auxiliary_file/rna_tissue_consensus.json", "r") as jsonfile:  
                        ALL_DICT = json.load(jsonfile)
            for idx, maf in enumerate(maf_output_list):
                maf_filter(maf, maf_flt_list, ALL_DICT, maf_filtered_list[idx])
                print(colored(("\n=> Finish filtering MAF file: " + maf_filtered_list[idx]+"\n"), 'green'))
            # MAF combination
            if len(maf_filtered_list) > 1:
                print(colored("Start MAF combination....\n", "yellow"))
                maf_df, head = pd.DataFrame(), ""
                for maf_file in maf_filtered_list:
                    head, maf = fast_read_maf(maf_file)
                    maf_df = pd.concat([maf_df, maf])
                maf_df.to_csv(folder+"maf_combination.maf", sep="\t", index= False, header = list(maf_df.columns.values))
                if head != "":
                    with open(folder+"maf_combination.maf", "r+") as filtered_file:
                        lines = filtered_file.readlines()
                        filtered_file.seek(0)
                        filtered_file.write(head)
                        filtered_file.writelines(lines)
                print(colored(("=> Finish combining MAF files to " + folder + "maf_combination.maf" + "\n"), 'green'))
        else:
            if len(category) > 1:
                print(colored("Start MAF combination....\n", "yellow"))
                maf_df, head = pd.DataFrame(), ""
                for maf_file in category:
                    head, maf = fast_read_maf(maf_file)
                    maf_df = pd.concat([maf_df, maf])
                maf_df.to_csv(folder+"maf_combination.maf", sep="\t", index= False, header = list(maf_df.columns.values))
                if head != "":
                    with open(folder+"maf_combination.maf", "r+") as filtered_file:
                        lines = filtered_file.readlines()
                        filtered_file.seek(0)
                        filtered_file.write(head)
                        filtered_file.writelines(lines)
                print(colored(("=> Finish combining MAF files to " + folder + "maf_combination.maf" + "\n"), 'green'))
            

if __name__ == '__main__':
    main()
