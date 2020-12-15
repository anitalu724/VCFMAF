############################################################################################
# FileName     [ maf_analysis.py ]
# PackageName  [ VCFMAF ]
# Synopsis     [ Control all functions ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020/8 ]
############################################################################################

import numpy as np
import pandas as pd
import argparse, textwrap
import ast
from src.maf_analysis import (
    CoMutAnalysis,
    SigMutatedGeneDetection,
    KnownCancerGeneAnnotation,
    MutationalSignature,
    TotalMutationBurden,
    OncoKBAnnotator,
    HRDScore,
    WGDnCIN
)




def main():
    parser = argparse.ArgumentParser(
        description="MAF analysis", formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-f", "--file", nargs=1, metavar="MAF file", required=True)
    parser.add_argument("-comut", "--comut_analysis", action="store_true")
    parser.add_argument("-smg", "--significantly_mutated_gene", action="store_true")
    parser.add_argument("-kcga", "--known_cancer_gene_annotaiton", action="store_true")
    parser.add_argument("-ms", "--mutational_signature", nargs=2, metavar=('step_#', 'list_of_params'),
                        help=textwrap.dedent('''\
                        Two steps are included: Estimation(0) and Plotting(1).
                        For each step, two values are required:
                        * Estimation
                          1. 1
                          2. [<Rank start number>, <Rank end number>, <Epoch number>]
                        * Plotting
                          1. 2
                          2. [Signature Number]'''))
                        
    parser.add_argument("-tmb","--total_mutation_burden",nargs=2,help="One item must be entered:\n \
                                                                       1. Sequencing Length\n",)
    parser.add_argument("-oncokb","--oncokb_annotator",nargs=5,help='Three items must be entered:\n \
                                                                     1. The relative path of the folder "oncokb-annotator".\n \
                                                                     2. The token of your OncoKB account.\n \
                                                                     3. The level of the drug. \n\
                                                                     4. The file path of clinical input.\n',)
    parser.add_argument("-hrd","--hrd_score",nargs=2,help="Two items must be entered:\n\
                                                           1. The CSV_input file.\n\
                                                           2. The reference for HRD Score.\n")
    parser.add_argument("-wgdcin", "--wgd_cin", nargs=1)
    parser.add_argument("-o","--output",required=True,metavar="OUTPUT folder",help="The path for storing every generated file.\n\
                                                                                    This path must end with a folder.\n")
    parser.add_argument("-p", "--picture", required=True,metavar="picture folder" ,help="The path for storing every picture.\n")                                                                        
    
    args = parser.parse_args()

    folder, pic = args.output, args.picture
    if folder[-1:] != "/":
        folder += "/"
    if pic[-1:] != "/":
        pic += "/"

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
        params = ast.literal_eval(args.mutational_signature[1])
        if args.mutational_signature[0] == '1':
            df.data_analysis(folder, pic, params[0], params[1], params[2])
        elif args.mutational_signature[0] == '2':
            df.plotting(folder, pic, params[0])
        
    # 5. Mutation burden statistics
    if args.total_mutation_burden:
        df = TotalMutationBurden(args.file[0])
        df.data_analysis(folder, args.total_mutation_burden[0], int(args.total_mutation_burden[1]))
    # 6. OncoKB Annotator
    if args.oncokb_annotator:
        df = OncoKBAnnotator(args.file[0])
        # df.data_analysis(folder,args.oncokb_annotator[0],args.oncokb_annotator[1],args.oncokb_annotator[3],args.oncokb_annotator[4])
        df.plotting(folder,pic, args.oncokb_annotator[2])
    # 7. HRD score
    if args.hrd_score:
        df = HRDScore(args.hrd_score[0])
        df.data_analysis(folder, args.hrd_score[1])
        df.plotting(folder, pic)
    #8. WGD and CIN
    if args.wgd_cin:
        df = WGDnCIN(args.wgd_cin[0])
        df.data_analysis(folder)
        df.plotting(folder, pic)

if __name__ == "__main__":
    main()
