############################################################################################
# FileName     [ maf_plotting.py ]
# PackageName  [ VCFMAF ]
# Synopsis     [ Visualization of MAF analysis ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020/9 ]
############################################################################################

import numpy as np
import pandas as pd
import argparse, textwrap
from src.maf_plot import CoMutPlot, CosineSimilarity, SBSPlot, DonutPlot


def main():
    parser = argparse.ArgumentParser(description="MAF plotting", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-cmp', '--comut_plot',nargs=4, help="Three items need to be entered:\n\
                                                              1. TSV file for data paths.\n\
                                                              2. TSV file which contain all information for image.\n\
                                                              3. plot color theme(0: cold, 1: warm)\n\
                                                              4. CoMut_plot picture file name")
    parser.add_argument("-cs", '--cosine_similarity',nargs=1, help="Enter the S2S_matrix.csv outputed from NMF analysis.")
    parser.add_argument("-sbs", "--sbs_plot", nargs=1)
    parser.add_argument("-dn", '--donut_plot', nargs=1)
    parser.add_argument("-o", "--output", required = True, help="The path for storing every generated file.\nThis path must end with a folder.\n")
    args = parser.parse_args()


    folder = args.output
    if folder[-1:] != "/":
        folder += "/"

    if args.comut_plot:
        plot1 = CoMutPlot(args.comut_plot[0], args.comut_plot[1])
        plot1.plot(folder, args.comut_plot[2], args.comut_plot[3])
    if args.cosine_similarity:
        plot = CosineSimilarity(args.cosine_similarity[0])
        plot.plot(folder)
    if args.sbs_plot:
        plot = SBSPlot(args.sbs_plot[0], folder)
        plot.plot(folder)
    if args.donut_plot:
        plot = DonutPlot(args.donut_plot[0])
        plot.plot(folder)
        



if __name__ == '__main__':
    main()