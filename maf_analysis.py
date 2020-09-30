############################################################################################
# FileName     [ maf_analysis.py ]
# PackageName  [ VCFMAF ]
# Synopsis     [ Control all functions ]
# Author       [ Cheng-Hua (Anita) Lu ]
# Copyright    [ 2020/8 ]
############################################################################################

import numpy as np
import pandas as pd
import argparse, textwrap
from src.maf_analysis import CoMutAnalysis, SigMutatedGeneDetection, KnownCancerGeneAnnotation, MutationalSignature, TotalMutationBurden

def main():
    parser = argparse.ArgumentParser(description="MAF analysis", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file", nargs=1, metavar="MAF file", required = True)
    parser.add_argument("-comut", "--comut_analysis", action="store_true")
    parser.add_argument("-smg", "--significantly_mutated_gene", action="store_true")
    parser.add_argument("-kcga", "--known_cancer_gene_annotaiton", action="store_true")
    parser.add_argument('-ms', '--mutational_signature', nargs=1, metavar="figures_folder")
    parser.add_argument("-tmb", "--total_mutation_burden", nargs=2, help="One item must be entered:\n \
                                                                          1. Sequencing Length\n" )
    parser.add_argument("-oncokb", "--oncokb_annotator", nargs=3, help="Three items must be entered:\n \
                                                                        1. The relative path of the folder \"oncokb-annotator\".\n \
                                                                        2. The token of your OncoKB account.\n \
                                                                        3. The level of the drug.")
    parser.add_argument("-o", "--output", required = True, metavar="OUTPUT folder",help="The path for storing every generated file.\nThis path must end with a folder.\n")
    args = parser.parse_args()

    folder = args.output
    if folder[-1:] != "/":
        folder += "/"

# 1. CoMut Plot Analysis
    if args.comut_analysis:
        df = CoMutAnalysis(args.file[0])
        df.data_analysis(folder)
# 2. Significantly mutated gene detection (Only oncodriveCLUST works until now)
    if args.significantly_mutated_gene:
        df = SigMutatedGeneDetection(args.file[0])
        df.oncodriveCLUST(folder)
        # df.dNdScv(folder)
# 3. Known cancer gene annotation
    if args.known_cancer_gene_annotaiton:
        df = KnownCancerGeneAnnotation(args.file[0])
        df.annotation(folder)   
# 4. Mutational signature
    if args.mutational_signature:
        df = MutationalSignature(args.file[0])
        if args.mutational_signature[0][-1:] != "/":
            args.mutational_signature[0] += "/"
        df.data_analysis(folder, args.mutational_signature[0])
# 5. Mutation burden statistics
    if args.total_mutation_burden:
        df = TotalMutationBurden(args.file[0])
        df.data_analysis(folder, args.total_mutation_burden[0], int(args.total_mutation_burden[1]))
# 6. OncoKB Annotator
    if args.oncokb_annotator:
        return



if __name__ == '__main__':
    main()

