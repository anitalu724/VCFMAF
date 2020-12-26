############################################################################################
# FileName     [ sig_mutated_gene_detect.py ]
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

# 1. Significantly mutated gene detection (Only oncodriveCLUST works until now)
class SigMutatedGeneDetection:
    """Significantly mutated gene detection (Only oncodriveCLUST works until now)

    Arguments:
        file        {string}    -- A MAF file path
        folder      {string}    -- The path for every output file

    Outputs:
        oncodriveCLUST inputs:
            oncodriveCLUST.nonsyn.txt
            oncodriveCLUST.syn.txt
        oncodriveCLUST output:
            oncodriveclust_results.tsv

    """    
    def __init__(self, file):
        print(colored(("\nStart Significantly_Mutated_Gene_Detection...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def oncodriveCLUST(self, folder):
        def get_input():
            selected_col = self.df[['Hugo_Symbol', 'Variant_Classification','Tumor_Sample_Barcode', 'Transcript_ID', 'Gene', 'Protein_position']]
            selected_col.columns = ['symbol','Sample','Tumor_Sample_Barcode', 'transcript', 'gene', 'position']
            list1, list2 = [],[]
            for idx in range(selected_col.shape[0]):
                if type(selected_col.iloc[idx]['position']) != float:
                    pos = selected_col.iloc[idx]['position'].split("/")[0]
                    if "-" not in pos:
                        selected_col.iloc[idx]['position'] = pos
                        if selected_col.iloc[idx]['Sample'] == "Silent":
                            list2.append(idx)
                        else:
                            list1.append(idx)
            file_1, file_2 = selected_col.iloc[list1, :], selected_col.iloc[list2, :]
            for idx in range(file_1.shape[0]):
                if file_1.iloc[idx]['Sample'] == 'Nonsense_Mutation':
                    file_1.iloc[idx]['Sample'] = "stop"
                else:
                    file_1.iloc[idx]['Sample'] = "non-synonymous"
            for idx in range(file_2.shape[0]):
                file_2.iloc[idx]['Sample'] = "synonymous"
            file_1 = file_1.loc[:,['symbol', 'gene', 'transcript', 'Tumor_Sample_Barcode', 'Sample', 'position']]
            file_1.columns = ['symbol', 'gene', 'transcript', 'Sample', 'ct', 'position']
            file_2 = file_2.loc[:,['symbol', 'gene', 'transcript', 'Tumor_Sample_Barcode', 'Sample', 'position']]
            file_2.columns = ['symbol', 'gene', 'transcript', 'Sample', 'ct', 'position']
            file_1.to_csv(folder+"oncodriveCLUST.nonsyn.txt", sep = '\t', index = False)
            file_2.to_csv(folder+"oncodriveCLUST.syn.txt", sep = '\t', index = False)
            print(colored(("=> Generate OncodriveCLUST's input files:"), 'green'))
            print(colored(("   "+folder+"oncodriveCLUST.nonsyn.txt"), 'green'))
            print(colored(("   "+folder+"oncodriveCLUST.syn.txt\n"), 'green'))
        def implement():
            os.system("oncodriveclust -m 3 --cgc src/auxiliary_file/CGC_phenotype.tsv "+\
                      folder+"oncodriveCLUST.nonsyn.txt "+folder+"oncodriveCLUST.syn.txt "+\
                      "src/auxiliary_file/gene_transcripts.tsv -o "+folder+"oncodriveclust_results.tsv\n")
        get_input()
        implement()
        print("\n")
    def dNdScv(self, folder):
        def get_input():          
            selected_col = self.df[['Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']]
            selected_col.columns = [ 'sampleID','chr', 'pos', 'ref', 'mut']
            selected_col.to_csv(folder+"dNdScv.input.txt", sep = '\t', index = False)
            print(colored(("\n=> Generate dNdScv's input file:"), 'green'))
            print(colored(("   "+folder+"dNdScv.input.txt"), 'green'))
        def implement():
            code = open(folder+"dNdScv_run.sh", "w")
            path = os.path.abspath(folder+"dNdScv.input.txt")
            code.write("file_path <- "+"\""+path+"\"\n")
            code.write("mutation <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)\n")
            code.write("library(\"dndscv\")\n")
            code.write("data(\"dataset_simbreast\", package=\"dndscv\")\n")
            code.write("dndsout = dndscv(mutation)\n")
            code.write("sel_cv = dndsout$sel_cv\n")
            code.write("signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c(\"gene_name\",\"qglobal_cv\")]\n")
            code.write("rownames(signif_genes) = NULL\n")
            code.write("write.table(signif_genes, file = \""+os.path.abspath(folder+"dNdScv.output.txt")+"\",row.names = FALSE,sep=\"\t\")\n")
            os.system("Rscript "+ folder+"dNdScv_run.sh"+" | head\n")
        get_input()
        implement()
        print("\n")
    def get_oncodriveFML_input(self, origin_file, output_file):
        tStart = time.time()
        print("start oncodriveFML...")
        selected_col = origin_file[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']]
        selected_col.columns = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']
        for i in range(len(selected_col['CHROMOSOME'])):
            new_str = selected_col['CHROMOSOME'][i][3:]
            selected_col['CHROMOSOME'][i] = new_str
        selected_col.to_csv(output_file, sep = '\t', index = False)
        tEnd = time.time()
        print("FML cost %f sec" % (tEnd - tStart))
    def get_twenty_plus_input(self, origin_file, output_file):
        tStart = time.time()
        print("start 20/20+...")
        selected_col = origin_file[['Hugo_Symbol','Tumor_Sample_Barcode','Entrez_Gene_Id', 'Chromosome', 'Start_Position', 'End_Position',  'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2', 'HGVSp_Short', 'Center']]
        for i in range(len(selected_col['Entrez_Gene_Id'])):
            selected_col.loc[i, 'Entrez_Gene_Id'] = "Carcinoma"
            selected_col.loc[i, 'Center'] = "c.?"
        selected_col.columns = ["Gene","Tumor_Sample","Tumor_Type","Chromosome","Start_Position","End_Position","Variant_Classification","Reference_Allele","Tumor_Allele","Protein_Change","DNA_Change"]
        selected_col.to_csv(output_file, sep = '\t', index = False)
        tEnd = time.time()
        print("20/20+ cost %f sec" % (tEnd - tStart))
    def implement_dNdScv(self, origin_file, input_file, code_file, output_file):
        tStart = time.time()
        print("start dNdScv...")
        selected_col = origin_file[['Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']]
        selected_col.columns = [ 'sampleID','chr', 'pos', 'ref', 'mut']
        for i in range(len(selected_col['chr'])):
            selected_col['chr'][i] = selected_col['chr'][i][3:]
        selected_col.to_csv(input_file, sep = '\t', index = False)

        code = open(code_file, "w")
        path = os.path.abspath(input_file)
        code.write("file_path <- "+"\""+path+"\"\n")
        code.write("mutation <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)\n")
        code.write("library(\"dndscv\")\n")
        code.write("data(\"dataset_simbreast\", package=\"dndscv\")\n")
        code.write("dndsout = dndscv(mutation)\n")
        code.write("sel_cv = dndsout$sel_cv\n")
        code.write("signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c(\"gene_name\",\"qglobal_cv\")]\n")
        code.write("rownames(signif_genes) = NULL\n")
        code.write("write.table(signif_genes, file = \""+os.path.abspath(output_file)+"\",row.names = FALSE,sep=\"\t\")\n")
        os.popen("Rscript "+ code_file+" | head\n")
        tEnd = time.time()
        print("dNdScv cost %f sec" % (tEnd - tStart))
