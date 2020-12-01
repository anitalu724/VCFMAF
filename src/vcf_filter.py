############################################################################################
# FileName     [ vcf_filter.py ]
# PackageName  [ src ]
# Synopsis     [ Sifting data ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020 4 ]
############################################################################################

import vcf
from .vcf_tool import *
from termcolor import colored

# 1.
# Select the chromosome(s)
# <select> must be a list (ex: [1,3,4])
def chromosome_select(select,record):            
    selected_list = []
    for i in select:
        selected_list.append(str(i))
    if record.CHROM in selected_list:
        return True
    else:
        return False

# Select the chromosomes' specific interval
# <select> must be a dictionary (ex: {"1":[123456, 234567], "2":[345678, 456789]})
# if want to touch the head or end, please enter "-1" at interval (ex: {"1":[-1, 2345679], "2":[1234567, -1]})
def interval_select(select,record):
    for i in select:
        if record.CHROM == (str(i)):
            if select[i][0] == -1 and select[i][1] == -1:
                return True
            elif select[i][0] == -1:
                if record.POS <= select[i][1]:
                    return True
                else:
                    return False
            elif select[i][1] == -1:
                if record.POS >= select[i][0]:
                    return True
                else:
                    return False
            else:
                if record.POS >= select[i][0] and record.POS <= select[i][1]:
                    return True
                else:
                    return False
        else:
            return False

# 2.
# Filter from caller information
# Caller: MuSe, Mutect2, SomaticSniper, Strelka2, VarScan2
def caller_info(call, record, DP_N, DP_T, AD_N, AD_T, AF_N, AF_T, NLOD, TLOD):
    if call == "Mutect2":
        PASS = True
        if record.INFO['NLOD'][0] < NLOD or record.INFO['TLOD'][0] < TLOD:
            return False
        for sample in record.samples:
            if "-N" in str(sample):
                if sample['DP'] < DP_N or sample['AD'][1] < AD_N or not(type(sample['AF']) is float and sample['AF'] >= float(AF_N)):
                    PASS = False
                else:
                    PASS = True
            elif "-T" in str(sample):
                if sample['DP'] < DP_T or sample['AD'][1] < AD_T or not(type(sample['AF']) is float and sample['AF'] >= float(AF_T)):
                    PASS = False
                else:
                    PASS = True
        return PASS
    elif call == "MuSe":
        PASS = True
        for sample in record.samples:
            if sample.sample == "NORMAL":
                DP = (sample['DP'] < DP_N)
                AD = (sample['AD'][1] < AD_N)
                AF = (round(sample['AD'][1]/sample['DP'],3) < AF_N)
                if DP or AD or AF or len(sample['AD'])> 2:
                    PASS = False
            elif sample.sample == "TUMOR":
                DP = (sample['DP'] < DP_T)
                AD = (sample['AD'][1] < AD_T)
                AF = (round(sample['AD'][1]/sample['DP'],3) < AF_T)
                if DP or AD or AF or len(sample['AD'])> 2:
                    PASS = False
            else:
                PASS = False
        return PASS
    elif call == "SomaticSniper":
        PASS = True
        for sample in record.samples:
            AD[0] = sample['DP4'][0]+sample['DP4'][1]
            AD[1] = sample['DP4'][2]+sample['DP4'][3]
            if "-N" in str(sample):
                DP = (sample['DP'] < DP_N)
                AD = (sample['AD'][1] < AD_N)
                AF = (round(sample['AD'][1]/sample['DP'],3) < AF_N)
                if DP or AD or AF:
                    PASS = False
            elif "-T" in str(sample):
                DP = (sample['DP'] < DP_T)
                AD = (sample['AD'][1] < AD_T)
                AF = (round(sample['AD'][1]/sample['DP'],3) < AF_T)
                if DP or AD or AF:
                    PASS = False
            else:
                PASS = False
        return PASS
    elif call == "VarScan2":
        PASS = True
        for sample in record.samples:
            if sample.sample == "NORMAL":
                DP = (sample['DP'] < DP_N)
                AD = (sample['AD'] < AD_N)
                AF = (round(sample['AD']/sample['DP'],3) < AF_N)
                if DP or AD or AF:
                    PASS = False
            elif sample.sample == "TUMOR":
                DP = (sample['DP'] < DP_T)
                AD = (sample['AD'] < AD_T)
                AF = (round(sample['AD']/sample['DP'],3) < AF_T)
                if DP or AD or AF:
                    PASS = False
            else:
                PASS = False
        return PASS
    elif call == "Strelka2":
        PASS = True
        for sample in record.samples:
            _AD = []
            if sample["AU"][0] > 0:
                _AD.append(sample["AU"][0]) 
            if sample["CU"][0] > 0:
                _AD.append(sample["CU"][0]) 
            if sample["GU"][0] > 0:
                _AD.append(sample["GU"][0]) 
            if sample["TU"][0] > 0:
                _AD.append(sample["TU"][0])
            if len(_AD) == 2:
                if _AD[1] > _AD[0]:
                    tmp = _AD[1]
                    _AD[1] = _AD[0]
                    _AD[0] = tmp
            
            if sample['DP'] == 0 or len(_AD)!= 2:
                return False
            if sample.sample == "NORMAL":
                DP = (sample['DP'] < DP_N)
                AD = (_AD[1] < AD_N)
                AF = (round(_AD[1]/sample['DP'],3) < AF_N)
                # print(round(_AD[1]/sample['DP'],3))
                if DP or AD or AF:
                    PASS = False
            elif sample.sample == "TUMOR":
                DP = (sample['DP'] < DP_T)
                AD = (_AD[1] < AD_T)
                AF = (round(_AD[1]/sample['DP'],3) < AF_T)
                if DP or AD or AF:
                    PASS = False
            else:
                PASS = False
        return PASS

# 3.
# Select the data that "pass" all the filters
def find_pass(record):
    if len(record.FILTER) == 0:
        return True
    else:
        return False

# 4.
# Artifact variant filter: FFPE filter
# Only for Mutect2
def FFPE_filter(call,record,delta):
    if call == "Mutect2":
        F1R2 = record.samples[1]['F1R2'][1]
        F2R1 = record.samples[1]['F2R1'][1]
        d = abs((F1R2-F2R1)/(F1R2+F2R1))
        if d <= delta:
            return True
        else:
            return False
    else:
        return False


            

########
# main #
########        
# vcf_read = read_vcf("../vcf_ex/MuSe/CBCP00111-T1.muse.variants.vcf")
# vcf_read = read_vcf("../vcf_ex/Mutect2/CBCP00111-T1_somatic_oncefiltered.vcf")
# vcf_read = read_vcf("../vcf_ex/SomaticSniper/CBCP00111-T1.somaticsniper_variants.vcf")
# vcf_read = read_vcf("../vcf_ex/Strelka2/CBCP00111-T1.strelka.somatic.snvs.vcf")
# vcf_read = read_vcf("../vcf_ex/VarScan2/CBCP00111-T1.snv.Somatic.hc.vcf")
# call = check_caller(vcf_read)
# vcf_writer = vcf.Writer(open("result.vcf","w"), vcf_read)
# chromosome_select([2,3,4],"result.vcf", vcf_read)
# interval_select({"1":[1234567, 2345678], "2":[3456787, 4567898]},"result.vcf", vcf_read)
# find_pass("result.vcf", vcf_read)
# caller_info(call, "result.vcf", vcf_read, 15, 15, 0, 0, 0, 0.05, 8, 8)

# for record in vcf_read:
#     if FFPE_filter(call, record, 0.9):
#         vcf_writer.write_record(record)

# for record in vcf_read:
#     if caller_info(call, record, 15, 15, 0, 0, 0, 0.01, 8, 8):
#         vcf_writer.write_record(record)

# for record in vcf_read:
#     for sample in record.samples:
#         if sample['DP'] == 0 and len(record.FILTER) == 0:
#             print(sample)

