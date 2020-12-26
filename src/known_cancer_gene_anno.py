############################################################################################
# FileName     [ known_cancer_gene_anno.py ]
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

# 2. Known cancer gene annotation
class KnownCancerGeneAnnotation:
    """Known cancer gene annotation

    Arguments:
        file        {string}    -- A MAF file path
        folder      {string}    -- The path for every output file

    Outputs:
        kcga.output.maf
            
    """ 
    def __init__(self, file):
        print(colored(("\nStart Known_Cancer_Gene_Annotation...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def annotation(self, folder):
        anno_file = "src/auxiliary_file/cancerGeneList.txt"
        origin_file = self.df
        output_file = folder+"kcga.output.maf"
        anno = pd.read_csv(anno_file, header = 0, sep = "\t")
        origin_file['#_of_occurrence_within_resources'] = "No"
        origin_file['OncoKB_Annotated'] = "No"
        origin_file['Is_Oncogene'] = "No"
        origin_file['Is_Tumor_Suppressor_Gene'] = "No"
        origin_file['MSK-IMPACT'] = "No"
        origin_file['MSK-HEME'] = "No"
        origin_file['FOUNDATION_ONE'] = "No"
        origin_file['FOUNDATION_ONE_HEME'] = "No"
        origin_file['Vogelstein'] = "No"
        origin_file['SANGER_CGC(05/30/2017)'] = "No"
        List = ['#_of_occurrence_within_resources', 'OncoKB_Annotated', 'Is_Oncogene','Is_Tumor_Suppressor_Gene','MSK-IMPACT','MSK-HEME','FOUNDATION_ONE','FOUNDATION_ONE_HEME','Vogelstein','SANGER_CGC(05/30/2017)']
        for item in range(len(anno)):
            origin_file.loc[origin_file['Hugo_Symbol'] == anno["Hugo Symbol"][item],List] = anno['# of occurrence within resources (Column D-J)'][item], anno['OncoKB Annotated'][item],anno['Is Oncogene'][item], anno['Is Tumor Suppressor Gene'][item], anno['MSK-IMPACT'][item], anno['MSK-HEME'][item], anno['FOUNDATION ONE'][item], anno['FOUNDATION ONE HEME'][item],anno['Vogelstein'][item],anno['SANGER CGC(05/30/2017)'][item]
        origin_file.to_csv(output_file, sep = '\t', index = False)
        print(colored("=> Generate output file: ", 'green'))
        print(colored(("   "+output_file), 'green'))
