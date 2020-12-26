############################################################################################
# FileName     [ comut_analysis.py ]
# PackageName  [ src ]
# Synopsis     [ Use the entered MAF files to do some analysis ]
# Author       [ LU, CHENG-HUA ]
# Copyright    [ 2020 12 ]
############################################################################################

from .maf_filter import *
import pandas as pd
import numpy as np
from numpy.lib.scimath import logn
from math import e
from termcolor import colored
import time
import os
import math
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as mtick
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from comut import comut
from comut import fileparsers

COLOR_MAP = ['#266199','#b7d5ea','#acc6aa','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22','#E08B69']

# 4.1 CoMut Plot Analysis
class CoMutAnalysis:
    """CoMut plot analysis

    Arguments:
        file        {string}    -- A MAF file path
        folder      {string}    -- The path for every output file

    Outputs:
        mutation_data.tsv
        mutation_burden.tsv
    
    """    
    def __init__(self, file):
        print(colored(("\nStart CoMut_Plot_Analysis...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder):
        def mutation_type():
            maf = self.df
            chosen_col = maf[['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification']]
            chosen_col = chosen_col.rename({'Tumor_Sample_Barcode':'sample', 'Hugo_Symbol':'category', 'Variant_Classification':'value'}, axis=1)
           
            value_old_list = ['Missense_Mutation', 'Nonsense_Mutation','In_Frame_Del', 'In_Frame_Ins','Splice_Site',
                              'Silent','Frame_Shift_Del','Frame_Shift_Ins','Nonstop_Mutation','Translation_Start_Site']
            remove_idx = []
            for i in range(len(chosen_col['value'])):
                if chosen_col['value'][i] not in value_old_list:
                    remove_idx.append(i)
                else:
                    if chosen_col['value'][i] == 'Missense_Mutation':
                        chosen_col['value'][i] = 'Missense'
                    elif chosen_col['value'][i] == 'Nonsense_Mutation':
                        chosen_col['value'][i] = 'Nonsense'
                    elif chosen_col['value'][i] == 'Nonstop_Mutation':
                        chosen_col['value'][i] = 'Nonstop'
                    elif chosen_col['value'][i] == 'In_Frame_Del' or chosen_col['value'][i] == 'In_Frame_Ins':
                        chosen_col['value'][i] = 'In frame indel'
                    elif chosen_col['value'][i] == 'Frame_Shift_Del' or chosen_col['value'][i] == 'Frame_Shift_Ins':
                        chosen_col['value'][i] = 'Frameshift indel'
                    elif chosen_col['value'][i] == 'Splice_Site':
                        chosen_col['value'][i] = 'Splice site'
                    elif chosen_col['value'][i] == 'Translation_Start_Site':
                        chosen_col['value'][i] = 'Translation start site'
            chosen_col = chosen_col.drop(chosen_col.index[remove_idx])
            # print(chosen_col)
            unique_chosen_col = chosen_col.drop_duplicates()
            # print(unique_chosen_col)
            # os._exit()
            unique_chosen_col.to_csv(folder+"mutation_data.tsv", sep = "\t", index = False)
            print(colored(("=> Generate CoMut_Analysis output files:"), 'green'))
            print(colored(("   "+folder+"mutation_data.tsv"), 'green'))
        def mutation_clonality():
            mutation_type_file = pd.read_csv(folder+"mutation_data.tsv", sep='\t', header=0)
            sample_dict = dict()
            for idx, data in mutation_type_file.iterrows():
                if data['sample'] not in sample_dict:
                    sample_dict[data['sample']] = [0, 0]
                if data['value'] == "Silent":
                    sample_dict[data['sample']][1]+=1
                else:
                    sample_dict[data['sample']][0]+=1

            mutation_clone = pd.DataFrame.from_dict(sample_dict, orient='index')
            mutation_clone.reset_index(level=0, inplace=True)
            mutation_clone.to_csv(folder+"mutation_burden.tsv", sep='\t', header=['sample', 'Nonsynonymous', 'Synonymous'], index=False)
            print(colored(("   "+folder+"mutation_burden.tsv"), 'green'))
        mutation_type()
        mutation_clonality() 
