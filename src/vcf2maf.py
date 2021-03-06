############################################################################################
# FileName     [ vcf2maf.py ]
# PackageName  [ src ]
# Synopsis     [ Use vcf2maf module to convert VCF file to MAF ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020.7 ]
############################################################################################

import vcfpy
from vcfpy import header
import os
from termcolor import colored
import datetime
import string
import multiprocessing


def vcf2vep2maf(vcf_file_list, maf_file_list, path, category, max_filter_ac):
    run_file = open(path+"run.sh", "w")
    run_file.write("YELLOW='\\033[0;33m'\n")
    run_file.write("GREEN='\\033[0;32m'\n")
    run_file.write("NC='\\033[0m'\n")

    perl_path = input("Please enter vcf2maf.pl's path: ")
    if perl_path == "":
        perl_path = "../vcf2maf.pl"
    fork = str(multiprocessing.cpu_count())
    if isinstance(vcf_file_list, list) and isinstance(maf_file_list, list) and len(vcf_file_list) == len(maf_file_list):
        for index, vcf_file in enumerate(vcf_file_list):
            run_file.write("printf \"${YELLOW}\nStart transforming file:\nVCF: "+vcf_file+"\nMAF: "+maf_file_list[index]+"\n\n${NC}\"\n\n")
            run_file.write("printf \"VCF file's size: "+str(os.stat(vcf_file).st_size/10**6)+" MB\"\n\n")
            run_file.write("printf \"\n\"\n")
            run_file.write("perl "+perl_path+" \\\n"+ \
                            "--tumor-id "+category[index][1]+" \\\n"+\
                            "--normal-id "+category[index][0]+" \\\n"+\
                            "--vcf-tumor-id "+category[index][1]+" \\\n"+\
                            "--vcf-normal-id "+category[index][0]+" \\\n"+\
                            "--max-filter-ac "+max_filter_ac+" \\\n"+\
                            "--input-vcf "+vcf_file+" \\\n"+\
                            "--output-maf "+maf_file_list[index]+" \\\n")
            if fork != "":
                run_file.write("--vep-forks "+fork+" \\\n")
            run_file.write("\n")
            run_file.write("printf \"${GREEN}\n=> Finish transforming file: "+vcf_file+" to "+maf_file_list[index]+"${NC}\n\n\"\n\n")
    else:
        print(len(vcf_file_list),len(maf_file_list))
        print(colored("ERROR: Different number of VCF files and MAF files!\n", "red"))
    run_file.close()
    os.system("sh "+path+"/run.sh\n")
    os.system("rm "+path+"/run.sh")
