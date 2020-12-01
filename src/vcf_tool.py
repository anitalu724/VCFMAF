############################################################################################
# FileName     [ vcf_tool.py ]
# PackageName  [ src ]
# Synopsis     [ Assisting vcf filters ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020 4 ]
############################################################################################

import vcf
# Read vcf file
def read_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file),'r')
    return vcf_reader
    
# Check the caller of the vcf file
def check_caller(vcf_reader):
    caller = ""
    if "source" in vcf_reader.metadata:
        if vcf_reader.metadata["source"][0] == 'strelka':
            caller = "Strelka2"
        elif vcf_reader.metadata["source"][0] == 'VarScan2':
            caller = "VarScan2"
        elif vcf_reader.metadata["source"][1] == 'Mutect2':
            caller = "Mutect2"
    else:
        if "MuSE_version" in vcf_reader.metadata:
            caller = "MuSe"
        elif "BCOUNT" in vcf_reader.formats:
            caller = "SomaticSniper"
    return caller
