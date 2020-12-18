############################################################################################
# FileName     [ maf_analysis.py ]
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

# 3. Mutation burden statistics
class TotalMutationBurden:
    """Mutation burden statistics

    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for every output file
        length          {int}       -- The length of genome (WGS = 60456963)

    Outputs:
        TMB_analysis.tsv
        TMB_statistic.tsv
            
    """ 
    def __init__(self, file):
        print(colored(("\nStart Total Mutation Burden...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder, length):
        select_df = self.df[['Variant_Classification', 'Tumor_Sample_Barcode']]
        non = ["Missense_Mutation","Splice_Site", "Translation_Start_Site","Nonstop_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins", "Nonsense_Mutation"]
        sample = select_df['Tumor_Sample_Barcode'].unique()
        sample_dict = {s:[0,0] for s in sample}     #list[all, nonsynonmous]
        for i in range(select_df.shape[0]):
            barcode, variant = select_df.iloc[i,:]['Tumor_Sample_Barcode'], select_df.iloc[i,:]['Variant_Classification']
            sample_dict[barcode][0] += 1
            if variant in non:
                sample_dict[barcode][1] += 1
            
        for s in sample_dict:
            sample_dict[s].append(sample_dict[s][0]*1000000/length)
            sample_dict[s].append(sample_dict[s][1]*1000000/length)
        d = pd.DataFrame(sample_dict).T.reset_index()
        d.columns = ['sample', 'All', 'nonsynonmous', 'TMB_All', 'TMB_nonsynonmous']
        stat_dict = {'mean' : [d['nonsynonmous'].sum()/len(sample), d['TMB_All'].sum()/len(sample), d['TMB_nonsynonmous'].sum()/len(sample)], 
                     'median' : [d['nonsynonmous'].median(), d['TMB_All'].median(), d['TMB_nonsynonmous'].median()], 
                     'max' : [d['nonsynonmous'].max(), d['TMB_All'].max(), d['TMB_nonsynonmous'].max()], 
                     'min' : [d['nonsynonmous'].min(), d['TMB_All'].min(), d['TMB_nonsynonmous'].min()],
                     'Q1' : [d['nonsynonmous'].quantile(q=0.25), d['TMB_All'].quantile(q=0.25), d['TMB_nonsynonmous'].quantile(q=0.25)],
                     'Q3' : [d['nonsynonmous'].quantile(q=0.75), d['TMB_All'].quantile(q=0.75), d['TMB_nonsynonmous'].quantile(q=0.75)]}
        stat_dict_df = pd.DataFrame(stat_dict).T
        stat_dict_df.columns = ['nonsynonmous', 'TMB_All', 'TMB_nonsynonmous']
        d.to_csv(folder+"TMB_analysis.tsv", sep="\t")
        stat_dict_df.to_csv(folder+"TMB_statistic.tsv", sep="\t")
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   "+folder+"TMB_analysis.tsv"), 'green'))
        print(colored(("   "+folder+"TMB_statistic.tsv"), 'green'))



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
            print(unique_chosen_col)
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

# 4.2 CoMutPlotting
class CoMutPlot:
    def __init__(self, file, info):
        print(colored(("\nStart plotting CoMut Plot...."), 'yellow'))
        self.fd = (pd.read_csv(file, sep="\t")).to_dict('list')
        for item in self.fd:
            cleanedList = [x for x in self.fd[item] if str(x) != 'nan']
            self.fd[item] = cleanedList
        self.info = (pd.read_csv(info, sep="\t")).to_dict('list')
        for item in self.info:
            cleanedList = [x for x in self.info[item] if str(x) != 'nan']
            self.info[item] = cleanedList
    def plot(self, folder, theme, name):
        fixed_category = ['Same Patient','Copy Number Alteration','Mutation Type','Purity','Mutation Signature','Mutation Classification','Frequency', 'Whole Genome Doubling']
        feature_category = [item for item in self.fd if item not in fixed_category]
        fixed_path_list = [self.fd[item][0] if len(self.fd[item]) != 0 else "" for item in fixed_category ]
        
        feature_path_list = [self.fd[item][0] for item in feature_category ]
      
        print(colored(("*   Reading Data...."), 'yellow'))
        same_patient_df = pd.read_csv(fixed_path_list[0], sep = '\t') if fixed_path_list[0] != "" else pd.DataFrame(data={})
        copy_num_df     = pd.read_csv(fixed_path_list[1], sep = '\t') if fixed_path_list[1] != "" else pd.DataFrame(data={})
        mut_type_df     = pd.read_csv(fixed_path_list[2], sep = '\t') if fixed_path_list[2] != "" else pd.DataFrame(data={})
        purity_df       = pd.read_csv(fixed_path_list[3], sep = '\t') if fixed_path_list[3] != "" else pd.DataFrame(data={})
        mut_sig_df      = pd.read_csv(fixed_path_list[4], sep = '\t') if fixed_path_list[4] != "" else pd.DataFrame(data={})
        mut_clone_df    = pd.read_csv(fixed_path_list[5], sep = '\t') if fixed_path_list[5] != "" else pd.DataFrame(data={})
        freq_df         = pd.read_csv(fixed_path_list[6], sep = '\t') if fixed_path_list[6] != "" else pd.DataFrame(data={})
        wgd_df          = pd.read_csv(fixed_path_list[7], sep = '\t') if fixed_path_list[7] != "" else pd.DataFrame(data={})
        feature_df = [pd.read_csv(i, sep = '\t') for i in feature_path_list]
        
        # primary_df      = pd.read_csv(path_list[4], sep = '\t') if path_list[4] != "" else pd.DataFrame(data={})
        # best_res_df     = pd.read_csv(path_list[5], sep = '\t') if path_list[5] != "" else pd.DataFrame(data={})

        cold_set = {'mut':['#b7d5ea','#acc6aa','#266199','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22'],
                    'purity':'Blues',
                    'burden':['#266199','#b7d5ea'],
                    'sig':['#b7d5ea','#acc6aa','#266199','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22'],
                    'copy_num':['#E5E5E5','#b7d5ea','#266199', {'facecolor': 'black', 'alpha': 0.5}, '#7566AA','#000060', '#71a0a5'],
                    'feature':['#DDDDDD','#acc6aa','#71a0a5','#77628c', '#00033F', '#E0CADB','#b7d5ea','888888'],
                    'primary':{'Skin':'#E5E5E5', 'Acral': '#8ec6c5', 'Occult':'#6983aa', 'Mucosal': '#8566aa'},
                    'response':{'SD': '#DDDDDD', 'MR': '#acc6aa', 'PR': '#71a0a5', 'CR':'#77628c', 'PD':  '#00033F'}}

        warm_set = {'mut':['#fbc687','#D4E2C3','#ea907a','#AA1D22','#F4D764','#EBEDAA','#BBBBBB','#aacdbe','#AA7733'],
                    'purity':'Oranges',
                    'burden':['#ea907a','#fbc687'],
                    'sig':['#ea907a','#fbc687','#f4f7c5','#aacdbe'],
                    'copy_num':['#eeeeee','#EEBBBB','#AA1C19', {'facecolor': 'black', 'alpha': 0.2}, '#CCDD99','#575033'],
                    'primary':{'Skin':'#f7e7bd', 'Acral': '#d9c6a5', 'Occult':'#a35d6a', 'Mucosal': '#c26565'},
                    'response':{'SD': '#ea907a', 'MR': '#fbc687', 'PR': '#f4f7c5', 'CR':'#aacdbe', 'PD':  '#ACBF89'}}

        used_set = cold_set if theme == "0" else warm_set

    # info
        #wgd
        wgd_mapping = {'Yes': {'facecolor': 'black', 'alpha': 0.5}, 'No': 'white'}
        #same patient
        indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 5}
        #copy number
        cna = ['Baseline', 'Allelic amplification', 'Allelic deletion','aCN = 0', 'CN-LOH', 'Complex', 'Multiple']
        cna_mapping = dict(zip(cna, used_set['copy_num']))
        cna_order = self.info['Copy Number Alteration']
        #mut
        category_order = self.info['Mutation Type']
        item = ['Multiple','Silent','Missense','Frameshift indel','In frame indel','Translation start site','Nonstop','Splice site','Nonsense']
        mut_mapping = dict(zip(item, used_set['mut']))
        #signature
        sig_kwargs = {'edgecolor': 'white', 'width': 0.96, 'linewidth': 0.5}
        sig = self.info['Mutation Signature']
        if len(sig) >= 10:
            print(colored('ERROR: The number of mutation signature must be less than 10!\n', 'red'))
            print(colored('CoMut Plotting Halt\n', 'red'))
            return
        sig_mapping = {sig[i]:used_set['sig'][i] for i in range(len(sig))}
        #burden(clonality)
        burden = ['Nonsynonymous','Synonymous']
        bar_mapping = dict(zip(burden, used_set['burden']))
        bar_kwargs = {'width': 0.8}
    
    # features
        feature_mapping = []
        for feature in feature_category:
            fea = self.info[feature]
            if len(fea) > 8:
                print(colored(('ERROR: The number of type of \"'+feature+"\" must be less than 8!\n"), 'red'))
                print(colored('CoMut Plotting Halt\n', 'red'))
                return
            fea_mapping = {fea[i]:used_set['feature'][i] for i in range(len(fea))}
            feature_mapping.append(fea_mapping)
        from matplotlib import rcParams
        custom_rcParams = {'font.family': 'arial','font.size': 14}   
        rcParams.update(custom_rcParams)


    # Plot
        toy_comut = comut.CoMut()
        toy_comut.samples = list(wgd_df['sample'])
        width, fig_info = 0, 0
        FigCol = 2 if width > 10 else 1
        height = int(width/10)*8
        # print(width)
        FigSize = (10, 10)
        
        if not same_patient_df.empty:
        #     width += 1
        #     fig_info += 1
            toy_comut.add_sample_indicators(same_patient_df, name = fixed_category[0], plot_kwargs = indicator_kwargs) 
        if not wgd_df.empty:
            toy_comut.add_categorical_data(wgd_df, name=fixed_category[7], mapping = wgd_mapping)

        if not copy_num_df.empty:
        #     width += 2
        #     fig_info += 7
            toy_comut.add_categorical_data(copy_num_df, name=fixed_category[1], mapping = cna_mapping, category_order = cna_order) 
        if not mut_type_df.empty:
            # width += 5*int(len(self.info['Mutation Type'])/10)
            # fig_info += 9
            toy_comut.add_categorical_data(mut_type_df, name = fixed_category[2], category_order = category_order, mapping = mut_mapping)
        
        if not purity_df.empty:
        #     width += 1
            toy_comut.add_continuous_data(purity_df, name = fixed_category[3], mapping = used_set['purity'], value_range = (0, 1))
        
        # for i in range(len(feature_df)):
        #     width += 1
        #     toy_comut.add_categorical_data(feature_df[i], name=feature_category[i], mapping = feature_mapping[i])
        
        if not mut_sig_df.empty:
        #     width += 2
            toy_comut.add_bar_data(mut_sig_df, name = fixed_category[4], mapping = sig_mapping, stacked = True, ylabel = 'Mutational\nsignatures',bar_kwargs = sig_kwargs)
        if not mut_clone_df.empty:
        #     # width += 2
            toy_comut.add_bar_data(mut_clone_df, name = fixed_category[5], mapping = bar_mapping, stacked = True, bar_kwargs = bar_kwargs, ylabel = 'Muts/Mb')
        
        
        toy_comut.plot_comut(figsize = FigSize, x_padding = 0.04, y_padding = 0.04, tri_padding = 0.03)
        toy_comut.add_unified_legend(ncol = 1)
        toy_comut.axes[fixed_category[1]].set_xticklabels([])
        toy_comut.figure.savefig(folder+name, dpi = 300, bbox_inches = 'tight')
        os._exit()
        
        # print(FigSize)
        
        
        # patient_leg = toy_comut.add_axis_legend('Same Patient', bbox_to_anchor = (1.025, 1.2), frameon = False, numpoints = 2)
        toy_comut.plot_comut(figsize = FigSize, x_padding = 0.04, y_padding = 0.04, tri_padding = 0.03)
        toy_comut.add_unified_legend(ncol = 2)
        toy_comut.figure.savefig(folder+name, dpi = 300, bbox_inches = 'tight')
        print(colored(("=> Generate CoMut Plot: "+folder+name), 'green'))

    





# 5. Mutational signature (V)
class MutationalSignature:
    """Mutational signature

    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for output files
        pic             {string}    -- The path especially for output figures(.pdf)
        rank1, rank2    {int}       -- The range for estimate # signature
        epoch           {int}       -- # estimation running
        sig             {int}       -- The final factorization rank(# signature)
    
    Outputs:
        ms_input.tsv
        96_sig.csv
        sig_sample.csv
        SBS.tsv
        
    Pictures:
        Estimation.pdf
        SBS_96_plots.pdf
        S2S.pdf
        SigContribution.pdf
        SigSamHeatmap.pdf
        Donut_plot.pdf

    """
    def __init__(self, file):
        print(colored(("\nStart Mutational_Signature...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder, pic, rank1, rank2, epoch):
        def get_input_file():
            output_file = folder+"ms_input.tsv"
            selected_col = self.df[['Tumor_Sample_Barcode','flanking_bps', 'Reference_Allele', 'Tumor_Seq_Allele2']]
            selected_col.columns = ['SampleID', 'Three_Allele', 'Ref', 'Mut']
            sample_list = selected_col.SampleID.unique()
            grouped = selected_col.groupby(selected_col['SampleID'])
            df_list = [grouped.get_group(sample).reset_index(drop=True) for sample in sample_list]
            final_dict = {}
            for d, df in enumerate(df_list):
                # order: 'C>A','C>G','C>T','T>A','T>C','T>G'
                cata_list = [[],[],[],[],[],[]]
                for i in range(len(df)):
                    item = df.loc[i]
                    if (item['Ref'] == 'C' and item['Mut'] == 'A') or (item['Ref'] == 'G' and item['Mut'] == 'T'):
                        cata_list[0].append(item)
                    elif (item['Ref'] == 'C' and item['Mut'] == 'G') or (item['Ref'] == 'G' and item['Mut'] == 'C'):
                        cata_list[1].append(item)
                    elif (item['Ref'] == 'C' and item['Mut'] == 'T') or (item['Ref'] == 'G' and item['Mut'] == 'A'):
                        cata_list[2].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'A') or (item['Ref'] == 'A' and item['Mut'] == 'T'):
                        cata_list[3].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'C') or (item['Ref'] == 'A' and item['Mut'] == 'G'):                       
                        cata_list[4].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'G') or (item['Ref'] == 'A' and item['Mut'] == 'C'):
                        cata_list[5].append(item)
                list_96 = []
                for cata in range(len(cata_list)):
                    cata_sum_list = [int(0)]*16
                    if cata in [0,1,2]:
                        three_allele_dict={"ACA":0,     "TGT":0,    "ACC":1,    "GGT":1,    "ACG":2,    "CGT":2,    "ACT":3,    "AGT":3, \
                                           "CCA":4,     "TGG":4,    "CCC":5,    "GGG":5,    "CCG":6,    "CGG":6,    "CCT":7,    "AGG":7, \
                                           "GCA":8,     "TGC":8,    "GCC":9,    "GGC":9,    "GCG":10,   "CGC":10,   "GCT":11,   "AGC":11,\
                                           "TCA":12,    "TGA":12,   "TCC":13,   "GGA":13,   "TCG":14,   "CGA":14,   "TCT":15,   "AGA":15 }   
                    elif cata in [3,4,5]:
                        three_allele_dict={"ATA":0,     "TAT":0,    "ATC":1,    "GAT":1,    "ATG":2,    "CAT":2,    "ATT":3,    "AAT":3, \
                                           "CTA":4,     "TAG":4,    "CTC":5,    "GAG":5,    "CTG":6,    "CAG":6,    "CTT":7,    "AAG":7, \
                                           "GTA":8,     "TAC":8,    "GTC":9,    "GAC":9,    "GTG":10,   "CAC":10,   "GTT":11,   "AAC":11,\
                                           "TTA":12,    "TAA":12,   "TTC":13,   "GAA":13,   "TTG":14,   "CAA":14,   "TTT":15,   "AAA":15 }  

                    for j in range(len(cata_list[cata])):
                        if (cata_list[cata][j])['Three_Allele'] in three_allele_dict:
                            cata_sum_list[three_allele_dict[(cata_list[cata][j])['Three_Allele']]] += 1;
                    list_96 += cata_sum_list
                final_dict[sample_list[d]] = list_96

            new_df = pd.DataFrame.from_dict(final_dict)
            list_a = ["A.A", "A.C", "A.G", "A.T", "C.A", "C.C", "C.G", "C.T",\
                      "G.A", "G.C", "G.G", "G.T", "T.A", "T.C", "T.G", "T.T"]
            list_b = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
            new_row_name = []
            for item in list_b:
                for allele in list_a:
                    new_str = allele[0]+"["+item+"]"+allele[2]
                    new_row_name.append(new_str)
            new_df.index = new_row_name
            new_df.to_csv(output_file, sep = '\t', index = True)
            print(colored("=> Generate input file: ", 'green'))
            print(colored(("   "+output_file), 'green'))
        def estimation():
            os.system("git clone https://github.com/mims-harvard/nimfa.git\n")
            os.chdir("nimfa")
            os.system("python3 setup.py install --user")
            code = open("nimfa.py", "w")
            code.write("import nimfa\nfrom collections import defaultdict, Counter\nimport urllib\nimport numpy as np\nfrom matplotlib import pyplot as plt\nimport matplotlib.gridspec as gridspec\nfrom sklearn import preprocessing\nimport scipy.cluster.hierarchy as sch\nimport pandas as pd\n")
            code.write("df = (pd.read_csv(\"../"+folder+"ms_input.tsv\", sep=\"\t\")).T\n")
            code.write("data = (df.to_numpy())[1:]\n")
            code.write("rank_cands = range("+str(rank1)+","+ str(rank2)+", 1)\n")
            code.write("snmf = nimfa.Snmf(data, seed='random_vcol', max_iter=100)\n")
            code.write("summary = snmf.estimate_rank(rank_range=rank_cands, n_run="+str(epoch)+", what='all')\n")
            
            code.write("rss = [summary[rank]['rss'] for rank in rank_cands]\n")
            code.write("coph = [summary[rank]['cophenetic'] for rank in rank_cands]\n")
            code.write("disp = [summary[rank]['dispersion'] for rank in rank_cands]\n")
            code.write("spar = [summary[rank]['sparseness'] for rank in rank_cands]\n")
            code.write("spar_w, spar_h = zip(*spar)\n")
            code.write("evar = [summary[rank]['evar'] for rank in rank_cands]\n")
            code.write("fig, axs = plt.subplots(2, 3, figsize=(12,8))\n")
            code.write("axs[0,0].plot(rank_cands, rss, 'o-', color='#266199', label='RSS', linewidth=3)\n")
            code.write("axs[0,0].set_title('RSS', fontsize=16,fontweight='bold')\n")
            code.write("axs[0,0].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[0,0].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[0,1].plot(rank_cands, coph, 'o-', color='#695D73', label='Cophenetic correlation', linewidth=3)\n")
            code.write("axs[0,1].set_title('Cophenetic', fontsize=16,fontweight='bold')\n")
            code.write("axs[0,1].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[0,1].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[0,2].plot(rank_cands, disp,'o-', color='#71a0a5', label='Dispersion', linewidth=3)\n")
            code.write("axs[0,2].set_title('Dispersion', fontsize=16,fontweight='bold')\n")
            code.write("axs[0,2].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[0,2].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[1,0].plot(rank_cands, spar_w, 'o-', color='#B88655', label='Sparsity (Basis)', linewidth=3)\n")
            code.write("axs[1,0].set_title('Sparsity (Basis)', fontsize=16,fontweight='bold')\n")
            code.write("axs[1,0].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[1,0].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[1,1].plot(rank_cands, spar_h, 'o-', color='#E08B69', label='Sparsity (Mixture)', linewidth=3)\n")
            code.write("axs[1,1].set_title('Sparsity (Mixture)', fontsize=16,fontweight='bold')\n")
            code.write("axs[1,1].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[1,1].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("axs[1,2].plot(rank_cands, evar,  'o-', color='#841D22', label='Explained variance', linewidth=3)\n")
            code.write("axs[1,2].set_title('Explained variance', fontsize=16,fontweight='bold')\n")
            code.write("axs[1,2].tick_params(axis='both', labelsize=12)\n")
            code.write("axs[1,2].set_xticks(np.arange("+str(rank1)+", "+str(rank2)+", 1))\n")
            code.write("fig.tight_layout(pad=1.0)\n")
            code.write("plt.savefig(\"../"+pic+"Estimation.pdf\",dpi=300,bbox_inches = 'tight')\n")
            code.close()
            print(colored(("\nStart Estimation(may need a few minutes)...."), 'yellow'))
            p = os.popen("python3 nimfa.py\n")
            x = p.read()
            print(x)
            p.close()
            print(colored("=> Generate estimation figure: ", 'green'))
            print(colored(("   "+pic+"Estimation.pdf\n"), 'green'))
            os.chdir("..")
            os.system("rm -rf nimfa\n")
        get_input_file()
        estimation()  
    def plotting(self, folder, pic, sig):
        LABEL_SIZE, TITLE_SIZE = 24,30
        print(colored(("\nStart Mutational_Signature Plotting(signature number must be in the range of 2 to 9)...."), 'yellow'))
        def nmf():
            print(colored(("\nStart NMF...."), 'yellow'))
            from sklearn.decomposition import NMF
            df = (pd.read_csv(folder+"ms_input.tsv", sep="\t")).T
            sample_list = df.index[1:]
            index_96 = df.to_numpy()[0]
            data = (df.to_numpy())[1:]
            model = NMF(n_components=int(sig),init='random', random_state=0)
            W = model.fit_transform(data)
            H = model.components_
            Hdf, Wdf = pd.DataFrame(H.T), pd.DataFrame(W.T)
            Hdf.columns = ["Signature "+str(i+1) for i in range(int(sig))]
            Wdf.columns = sample_list
            Hdf.index = index_96
            Wdf.index = ["Signature "+str(i+1) for i in range(int(sig))]
            Hdf.to_csv(folder+"96_sig.csv")
            Wdf.to_csv(folder+"sig_sample.csv")
            print(colored("=> Generate file: ", 'green'))
            print(colored(("   "+folder+"96_sig.csv"), 'green'))
            print(colored(("   "+folder+"sig_sample.csv"), 'green'))
        def SBSPlot():
            df = (pd.read_csv(folder+"96_sig.csv"))
            df = df.set_index(list(df.columns[[0]]))
            fig_x = tuple([ ' '+i[0]+' '+i[6] for i in list(df.index)])
            y_pos = np.arange(len(fig_x))
            fig_name = list(df.columns)
            fig, axes = plt.subplots(df.shape[1], 1, figsize=(12,2*df.shape[1]))#
            if df.shape[1] == 1:
                return
            for r in range(df.shape[1]):
                color_set = ['#02bdee', '#010101','#e32925','#cac9c9', '#a1cf63', '#ecc7c4']
                color_96 = [ c for c in color_set for i in range(16)]
                all_data = df.iloc[:, r]
                all_data /= (all_data.sum())
                maximum = max(all_data)*1.25
                data_list = all_data.tolist()
                axes[r].text(0.01, 0.86, fig_name[r], horizontalalignment='left',verticalalignment='center', transform=axes[r].transAxes, fontweight='bold')
                axes[r].bar(y_pos, data_list, color=color_96, width=0.4)
                axes[r].spines['bottom'].set_color('#cac9c9')
                axes[r].spines['top'].set_color('#cac9c9') 
                axes[r].spines['right'].set_color('#cac9c9')
                axes[r].spines['left'].set_color('#cac9c9')
                if r != df.shape[1]-1:
                    axes[r].xaxis.set_visible(False)
                    axes[r].set_xticklabels([])
                axes[r].tick_params(axis='x',length=0)
                axes[r].set_xlim([-0.8,len(data_list)-.8])

                axes[r].tick_params(axis='y',direction='in', color='#cac9c9', labelsize=10)
                axes[r].set_ylabel('Percentage', fontweight='bold')
                axes[r].tick_params(axis='y', labelsize=10)
                axes[r].set_ylim(top = max(all_data)*1.25)
                axes[r].yaxis.set_major_locator(ticker.LinearLocator(5))
                axes[r].yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=1))
                for i in range(6):
                    axes[r].add_patch(matplotlib.patches.Rectangle((0+16*i ,maximum*0.95), 15.6 , 0.01, color=color_set[i],transform=axes[r].transData))
            mut_list = ['C>A','C>G','C>T','T>A','T>C','T>G']
            for i in range(6):
                plt.text(0.19+0.13*i,0.916-df.shape[1]*0.0029, mut_list[i], horizontalalignment='center',verticalalignment='center',transform=plt.gcf().transFigure, fontweight='bold', fontsize=14)
            plt.xticks(y_pos, fig_x, color='#999999',rotation=90, fontsize=9,horizontalalignment='center',verticalalignment='top',fontname='monospace')#verticalalignment='bottom',
            space = 0.008075
            y_scale = [0.072, 0.084, 0.09, 0.094, 0.097, 0.0987, 0.1, 0.1013, 0.1023]
            for i in range(6):
                for j in range(16):
                    if i < 3:
                        plt.text((0.131+space*16*i)+space*j, y_scale[df.shape[1]-2], 'C',horizontalalignment='center',verticalalignment='center',transform=plt.gcf().transFigure, color=color_set[i], fontsize=9, rotation=90,fontname='monospace', fontweight='bold')
                    else:
                        plt.text((0.131+space*16*i)+space*j, y_scale[df.shape[1]-2], 'T',horizontalalignment='center',verticalalignment='center',transform=plt.gcf().transFigure, color=color_set[i], fontsize=9, rotation=90,fontname='monospace', fontweight='bold')
            plt.savefig(pic+"SBS_96_plots.pdf",dpi=300, bbox_inches='tight')
            print(colored(("=> Generate SBS Plot: "+pic+"SBS_96_plots.pdf"), 'green'))
        def CosineSimilarity():
            from sklearn.metrics.pairwise import cosine_similarity
            my_file, aux_file = folder+"96_sig.csv", "src/auxiliary_file/COSMIC_72.tsv"
            my_df, aux_df = pd.read_csv(my_file, index_col=0), pd.read_csv(aux_file, sep="\t",index_col=0)
            my_list, aux_list = my_df.columns, aux_df.columns
            X = np.array(my_df.T.to_numpy())
            Y = np.array(aux_df.T.to_numpy())
            M = cosine_similarity(X, Y, dense_output=True)
            Mdf= pd.DataFrame(M)
            Mdf.index, Mdf.columns = my_list, aux_list
            Mdf.to_csv(folder+"SBS.tsv", sep="\t")
            print(colored("=> Generate file: ", 'green'))
            print(colored(("   "+folder+"SBS.tsv"), 'green'))
            
            height, length = len(my_list), len(aux_list)
            sns.set(font_scale=2)
            sns.set_style("white")
            grid_kws = {"height_ratios": (.9, .2),"hspace": 0.3}  
            f, (ax, cbar_ax) = plt.subplots(2,figsize=(20,6), gridspec_kw=grid_kws)
            ax = sns.heatmap(M, vmin=0, vmax=1, xticklabels =aux_list, yticklabels = my_list, square=False, linewidth=1, cbar_ax=cbar_ax,ax=ax,
                                cmap="Blues",cbar_kws={"orientation": "horizontal",'shrink':1, 'aspect':70})
            # ax.set_title('Cosine Similarity',fontsize=TITLE_SIZE,weight='bold',pad=0,verticalalignment='bottom')
            ax.set_xticklabels(ax.get_xticklabels(),rotation=90, horizontalalignment='center', fontsize=LABEL_SIZE-6, color='#222222')
            ax.tick_params(axis='both',length=0)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=LABEL_SIZE-6,color='#222222',verticalalignment='center')
            plt.ylim(bottom=0, top=height+0.5)
            plt.savefig(pic+"S2S.pdf",dpi=300,bbox_inches='tight')
            plt.clf()
            print(colored(("=> Generate Cosine Similarity Plot: "+pic+"S2S.pdf"), 'green'))  
        def SigDistribution():
            df = pd.read_csv(folder+"sig_sample.csv", index_col=0)
            sample_list, sig_list = list(df.columns),list(df.index)
            SUM = (df.sum(axis = 0, skipna = True)).tolist()
            df = df/SUM
            dft = df.T
            # dft.columns = ['sample']+dft.columns
            dft.to_csv(folder+"SigContribution.tsv",index_label='sample', sep='\t')
            print(colored(("   "+folder+"SigContribution.tsv"), 'green'))
            ind = np.arange(df.shape[1])
            data = []
            for i in range(df.shape[0]):
                d = tuple(df.iloc[i].tolist())
                data.append(d)
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_axes([0,0,1,1])
            
            

            for i in range(len(data)):
                if i == 0:
                    ax.bar(ind, data[i], 0.8, color = COLOR_MAP[i])
                else:
                    b = np.array(data[0])
                    for k in range(1,i):
                        b = b+np.array(data[k])
                    ax.bar(ind, data[i], 0.8, bottom=b,color = COLOR_MAP[i])
            # ax.set_title('Relative Contribution',fontsize=TITLE_SIZE, fontweight='bold')
            ax.spines['bottom'].set_color('#cac9c9')
            ax.spines['top'].set_color('#FFFFFF') 
            ax.spines['right'].set_color('#FFFFFF')
            ax.spines['left'].set_color('#cac9c9')
            ax.set_xlim([-1,len(ind)])
            ax.tick_params(axis='y',direction='in', color='#cac9c9', labelsize=LABEL_SIZE-4)
            ax.tick_params(axis='x',direction='in', length=0)
            ax.xaxis.set_visible(False)
            ax.set_yticks(np.arange(0, 1+0.1, 0.25))
            ax.legend(title="",labels=sig_list,loc='lower center',ncol=3, fontsize=LABEL_SIZE-4, edgecolor='white',
                      labelspacing=0.5, bbox_to_anchor=(0.5, (-0.1-(math.ceil(len(sig_list)/3)*0.065))))
            plt.savefig(pic+"SigContribution.pdf", dpi=300,bbox_inches='tight')
            print(colored(("=> Generate Bar Plot: " + pic+"SigContribution.pdf"), 'green')) 
            
            height, length = len(sig_list), len(sample_list)  
            h_data = np.array(df.to_numpy())
            sns.set(font_scale=2)
            f,ax = plt.subplots(figsize=(9+length/20,2+height*0.3))
            ax = sns.heatmap(data, vmin=0, vmax=1, yticklabels = sig_list, linewidths=1,
                             square=False, cmap="Blues",cbar_kws={"orientation": "horizontal",'shrink':1, 'aspect':50})
            # ax.set_title('Signature Sample Heatmap', fontsize=TITLE_SIZE,weight='bold',va='bottom')
            ax.xaxis.set_visible(False)
            ax.set_xticklabels([])
            ax.tick_params(axis='both',length=0)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=LABEL_SIZE-4,color='#222222')
            plt.savefig(pic+"SigSamHeatmap.pdf",dpi=300,bbox_inches='tight')
            print(colored(("=> Generate Heatmap: "+pic+"SigSamHeatmap.pdf\n"), 'green'))
        def DonutPlot():
            df = pd.read_csv(folder+"sig_sample.csv", index_col=0)
            raw_data = df.sum(axis=1)/df.shape[1]
            SUM = raw_data.sum(axis=0)
            raw_data = raw_data/SUM
            names, sizes = list(raw_data.index), list(raw_data.iloc[:])
            names = [names[i]+": "+'{:.1%}'.format(sizes[i]) for i in range(len(sizes))]
            fig, ax = plt.subplots(figsize=(6, 8), subplot_kw=dict(aspect='equal'))
            wedges, texts = ax.pie(sizes, colors=COLOR_MAP[:len(names)],wedgeprops=dict(width=0.6,edgecolor='w',linewidth=2), startangle=-40) #,normalize=False

            bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0)
            kw = dict(arrowprops=dict(arrowstyle="-"),bbox=bbox_props, zorder=0, va="center")

            for i, p in enumerate(wedges):
                ang = (p.theta2 - p.theta1)/2. + p.theta1
                y = np.sin(np.deg2rad(ang))
                x = np.cos(np.deg2rad(ang))
                horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
                connectionstyle = "angle,angleA=0,angleB={}".format(ang)
                kw["arrowprops"].update({"connectionstyle": connectionstyle})
                ax.annotate(names[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),horizontalalignment=horizontalalignment, **kw, fontsize=LABEL_SIZE)
            plt.savefig(pic+"Donut_plot.pdf", dpi=300, bbox_inches='tight')
            print(colored(("=> Generate Donut Plot: "+pic+"Donut_plot.pdf"), 'green'))
        nmf()
        SBSPlot()
        DonutPlot()
        CosineSimilarity()
        SigDistribution()
        
# 6. HRD score (V)
class HRDScore:
    """HRD score

    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for output files
        ref             {string}    -- The reference genome used, grch38 or grch37 or mouse (default: grch38)
        pic             {string}    -- The path especially for output figures(.pdf)

    Outputs:
        all_HRDresults.csv

    Pictures:
        HRD_Score.pdf
        high_HRD_pie.pdf

    """
    def __init__(self, file):
        print(colored(("\nStart analysing HRD Score...."), 'yellow'))
        self.list = ((pd.read_csv(file, sep="\t"))[['CNV_input']].values.T)[0]
    def data_analysis(self, folder, ref):
        scar_r = open(folder + "scar.r", "a")
        scar_r.write("library(\"scarHRD\")\n")
        meta_list = []
        for i in self.list:
            scar_r.write("scar_score(\"" + i + "\", reference = \""+ref+"\", seqz = FALSE, outputdir = \"" + folder[:-1] + "\")\n")
        scar_r.close()
        os.system("Rscript " + folder + "scar.r\n")
        os.system("rm "+ folder + "scar.r\n")
        
        for file in os.listdir(folder):
            if file.endswith("_HRDresults.txt"):
                meta_list.append(file)
        meta_list.sort(reverse = True)
        final_df = pd.DataFrame()
        for meta in meta_list:
            df = pd.read_csv(folder+meta, sep="\t", index_col=False)
            final_df = pd.concat([df, final_df]) if not final_df.empty else df
            os.system("rm " + folder + meta + "\n")
        final_df.columns = [['Sample_id','HRD','Telomeric_AI','LST','HRD-sum']]
        final_df.to_csv(folder + "all_HRDresults.csv",  index=False)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + folder + "all_HRDresults.csv"), 'green'))
    def plotting(self, folder, pic):
        LABEL_SIZE, TITLE_SIZE = 24,30
        #Bar Plot
        df = pd.read_csv(folder+"all_HRDresults.csv")
        size = df.shape[0]
        HRD = tuple(list(df['HRD']))
        TAI = tuple(list(df['Telomeric_AI']))
        LST = tuple(list(df['LST']))
        SUM = list(df["HRD-sum"])
        Sample = tuple(list(df['Sample_id']))
        ind = np.arange(size)
        print(ind)
        width = 0.7
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_axes([0,0,1,1])
        ax.bar(ind, HRD, width, color=COLOR_MAP[7])
        ax.bar(ind, TAI, width, bottom=HRD, color=COLOR_MAP[2])
        ax.bar(ind, LST, width, bottom=np.array(TAI)+np.array(HRD), color=COLOR_MAP[6])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_color('#cac9c9')
        ax.spines['left'].set_color('#cac9c9')
        ax.set_ylabel('Scores', fontsize=LABEL_SIZE, fontweight='bold')
        ax.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax.tick_params(axis='y',direction='in', color='#cac9c9')
        ax.set_ylim(top = max(SUM)*1.25)
        # ax.set_title('HRD Scores',fontsize=TITLE_SIZE, fontweight='bold')
        # plt.xticks(ind, Sample,rotation=45,horizontalalignment='right',fontweight='light', fontsize=12)
        ax.set_xlim([-1,len(ind)])
        ax.xaxis.set_visible(False)
        plt.yticks(fontsize=LABEL_SIZE-4)
        ax.set_yticks(np.arange(0, max(SUM)*1.25+3, 10))
        ax.legend(labels=['HRD','Telomeric_AI','LST'], fontsize=LABEL_SIZE-4, edgecolor='white')
        plt.savefig(pic+"HRD_Score.pdf", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + pic + "HRD_Score.pdf"), 'green'))
        #Pie Plot
        over = len([i for i in SUM if i >=42])
        data = [over, size-over]
        labels='≧ 42','< 42'
        fig1, ax1 = plt.subplots()
        ax1.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, colors=[COLOR_MAP[7],COLOR_MAP[2]] ,textprops={'fontsize': LABEL_SIZE, 'size': LABEL_SIZE})
        ax1.axis('equal')
        # plt.title("Propotion of Samples with HRD Phenotype", fontsize=TITLE_SIZE, fontweight='bold')
        plt.savefig(pic+"high_HRD_pie.pdf", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + pic + "high_HRD_pie.pdf"), 'green'))

# 7. WGD and CIN (V)
class WGDnCIN:
    """HRD score

    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for output files
        pic             {string}    -- The path especially for output figures(.pdf)

    Outputs:
        WGD_result.csv
        CIN_result.csv

    Pictures:
        CIN_Score.pdf
        WGD_pie.pdf

    """
    def __init__(self, file):
        print(colored(("\nStart analysing WGD and CIN...."), 'yellow'))
        self.list = ((pd.read_csv(file, sep="\t"))[['CNV_input']].values.T)[0]
    def data_analysis(self, folder):
        genome_length = [249250621, 243199373, 198022430, 191154276, 180915260,
                         171115067, 159138663, 146364022, 141213431, 135534747,
                         135006516, 133851895, 115169878, 107349540, 102531392,
                         90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
        whole_length = sum(genome_length)
        WGD_list, CIN_list, sample_list = [], [], []
        for sample in self.list:
            df = pd.read_csv(sample, sep="\t")
            sample_name = pd.unique(df['SampleID'])[0]
            sample_list.append(sample_name)
            WGD_selected = df.loc[(df['A_cn'] >= 2)]
            CIN_selected = df.loc[(df['A_cn'] != 1) & (df['B_cn'] != 1)]
            WGD_SUM, CIN_SUM = 0, 0
            for i in range(WGD_selected.shape[0]):
                WGD_SUM += WGD_selected.iloc[i]['End_position']-WGD_selected.iloc[i]['Start_position']
            for i in range(CIN_selected.shape[0]):
                CIN_SUM += CIN_selected.iloc[i]['End_position']-CIN_selected.iloc[i]['Start_position']
            WGD_list.append(WGD_SUM >= 0.5*whole_length)
            CIN_list.append(CIN_SUM/whole_length)
        WGD_df = (pd.DataFrame([sample_list, WGD_list])).T
        CIN_df = (pd.DataFrame([sample_list, CIN_list])).T
        WGD_df.columns,CIN_df.columns = [['SampleID','WGD']], [['SampleID','CIN']]
        WGD_df.to_csv(folder + "WGD_result.csv",  index=False)
        CIN_df.to_csv(folder + "CIN_result.csv",  index=False)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + folder + "WGD_result.csv"), 'green'))
        print(colored(("   " + folder + "CIN_result.csv"), 'green'))
    def plotting(self, folder, pic):
        # WGD Pie Plot
        LABEL_SIZE, TITLE_SIZE = 24,30
        wgd_df = pd.read_csv(folder+"WGD_result.csv")
        data = [0,0]
        for i in range(wgd_df[['WGD']].shape[0]):
            if wgd_df[['WGD']].iloc[i]['WGD'] == False:
                data[1]+=1
            elif wgd_df[['WGD']].iloc[i]['WGD'] == True:
                data[0]+=1
        labels = 'WGD','Non-WGD'
        fig1, ax1 = plt.subplots()
        _, _, autotexts = ax1.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, colors=[COLOR_MAP[1],COLOR_MAP[0]] ,textprops={'fontsize': LABEL_SIZE}) #, 'size': LABEL_SIZE
        ax1.axis('equal')
        autotexts[0].set_color('black')
        autotexts[1].set_color('white')
        plt.savefig(pic+"WGD_pie.pdf", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + pic + "WGD_pie.pdf"), 'green'))
        
        # CIN Bar plot
        CIN_df = pd.read_csv(folder+"CIN_result.csv")
        size = CIN_df.shape[0]
        CIN = tuple(list(CIN_df['CIN']))
        Sample = tuple(list(CIN_df['SampleID']))
        ind = np.arange(size)
        width = 0.7
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_axes([0,0,1,1])
        ax.bar(ind, CIN, width, color=COLOR_MAP[0])
        ax.set_ylabel('Scores', fontsize=LABEL_SIZE, fontweight='bold')
        # ax.set_title('CIN Scores',fontsize=24, fontweight='bold')
        # plt.xticks(ind, Sample,rotation=45,horizontalalignment='right',fontweight='light', fontsize=12)
        plt.yticks(fontsize=LABEL_SIZE-4)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_color('#cac9c9')
        ax.spines['left'].set_color('#cac9c9')
        ax.set_xlim([-1,len(ind)])
        ax.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax.tick_params(axis='y',direction='in', color='#cac9c9')
        ax.set_yticks(np.arange(0, 1, 0.2))
        ax.xaxis.set_visible(False)
        plt.savefig(pic+"CIN_Score.pdf", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + pic + "CIN_Score.pdf"), 'green'))

# 8. OncoKB annotator (V)
class OncoKBAnnotator:
    """OncoKB annotator

    Arguments:
        file            {string}    -- A MAF file path
        folder          {string}    -- The path for output files
        path            {string}    -- The path for oncokb-annotator
        token           {string}    -- The personal token provided from OncoKB
        clinical        {string}    -- The path for clinical data
        pic             {string}    -- The path especially for output figures(.pdf)
        cna             {string}    -- (Optional) The path for cna data
        level           {string}    -- (Optional) The level the user chooses (default = 4)

    Outputs:
        maf_oncokb_output.txt
        clinical_oncokb_output.txt

    Pictures:
        oncokb_total_pie.pdf
        oncokb_freq_actionable_genes.pdf

    """
    def __init__(self, file):
        print(colored(("\nStart OncoKB annotator(drug)...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder, path, token, clinical, cna = ''):
        selected_df = (self.df[['NCBI_Build','Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'HGVSp',  'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']]).set_index("Hugo_Symbol")
        selected_df.to_csv(folder + "maf_oncokb_input.txt", sep="\t")
        ## Generate clinical and cna data automatically
        # sample_list = selected_df['Tumor_Sample_Barcode'].unique()
        # print(sample_list)
        # print(len(sample_list))
        # oncotree_code = ['BRCA']*len(sample_list)
        # clin_d = {'SAMPLE_ID':sample_list, 'ONCOTREE_CODE':oncotree_code}
        # clin_df = pd.DataFrame(clin_d)
        # clin_df.to_csv("examples/TCGA/clinical_input.txt",sep="\t",index=False)
        # dfa = pd.DataFrame(np.random.randint(-1,1,size=(23058,len(sample_list))), columns=sample_list)
        # first_3_col = (pd.read_csv("examples/TCGA/TCGA_Mutect_v10_white_rerun_cna_part.txt",sep="\t")).iloc[:,[0,1,2]]
        # cna_df = pd.concat([first_3_col, dfa], axis=1)
        # cna_df.to_csv("examples/TCGA/cna.txt",sep="\t",index=False)
        # os._exit()
        #
        # os.system("pwd")
        # os._exit()
        
        os.system("git clone https://github.com/oncokb/oncokb-annotator.git\n")
        os.system('cp src/auxiliary_file/autoChange.py oncokb-annotator\n')
        os.chdir("oncokb-annotator")
        os.system('python3 autoChange.py\n')
        os.system('pip3 install requests\n')
        
        p = os.popen("python3 MafAnnotator.py -i ../"+folder + "maf_oncokb_input.txt -o ../" + folder + "maf_oncokb_output.txt -c ../"+clinical+" -b " + token + "\n")
        print(p.read())
        p.close()
        p = os.popen("python3 ClinicalDataAnnotator.py -i ../"+clinical+" -o ../"+ folder +"clinical_oncokb_output.txt -a ../"+folder+"maf_oncokb_output.txt\n")
        print(p.read())
        p.close()
        if cna !='':
            p = os.popen("python3 CnaAnnotator.py -i ../"+cna+" -o ../"+folder+"cna_oncokb_output.txt -c ../"+clinical+" -b "+ token + "\n")
            print(p.read())
            p.close()
        os.chdir("..")
        os.system("rm -rf oncokb-annotator\n")
        # os.system("rm "+folder+"maf_oncokb_input.txt\n")
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + folder + "maf_oncokb_output.txt"), 'green'))
        print(colored(("   " + folder + "clinical_oncokb_output.txt"), 'green'))
        # print(colored(("   "+folder+"cna_oncokb_output.txt"), 'green'))
    def plotting(self, folder, pic, level='4'):
        
        LABEL_SIZE, TITLE_SIZE = 24,30
        self.file = folder + "clinical_oncokb_output.txt"
        df = pd.read_csv(self.file, sep="\t")
        df_level = df[['HIGHEST_LEVEL']]
        level_list = ['LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4']
        level_dict = dict.fromkeys(level_list,0)
        sample_size = df.shape[0]
        for i in range(sample_size):
            if df_level.iloc[i]['HIGHEST_LEVEL'] in level_list:
                level_dict[df_level.iloc[i]['HIGHEST_LEVEL']] += 1

        true_num = 0
        if level == '4':
            true_num = sum(level_dict.values())
        elif level == '3':
            true_num = sum(level_dict.values()) - level_dict['LEVEL_4']
            level_list = ['LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B']
        elif level == '2':
            true_num = level_dict['LEVEL_1'] + level_dict['LEVEL_2']
            level_list = ['LEVEL_1', 'LEVEL_2']
        elif level == '1':
            true_num = level_dict['LEVEL_1']
            level_list = ['LEVEL_1']
        # Pie Plot( Total pie plot )
        size = [true_num, sample_size - true_num]
        labels = "Actionable\nbiomarkers"," Current absence\n of new actionable\n biomarkers"
        fig1, ax1 = plt.subplots()
        _, _, autotexts = ax1.pie(size, labels=labels, autopct='%1.1f%%', startangle=90, colors=[COLOR_MAP[3],COLOR_MAP[4]] ,textprops={'fontsize': LABEL_SIZE})
        autotexts[1].set_color('white')
        ax1.axis('equal')
        # plt.title("Total", fontsize=18, fontweight='bold')
        plt.savefig(pic+"oncokb_total_pie.pdf", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + pic + "oncokb_total_pie.pdf"), 'green'))
        
        # Bar Plot( Frequency of Actionable Genes )
        df_drug_count = df[level_list]
        drug_total_dict = {}
        for i in range(sample_size):
            drug_list = [[], [], [], [], [], []]    # [total, level1, level2, level3A, level3B, level4]
            for idx, item in enumerate(level_list):
                data = df_drug_count.iloc[i][item]
                if not pd.isna(data):
                    new_drug_list = data.split("(")
                    new_drug_list.pop(0)
                    new_drug_list = [drug.split(" ")[0] for drug in new_drug_list]
                    for j in new_drug_list:
                        if j not in drug_list[0]:
                            drug_list[0].append(j)
                            drug_list[idx+1].append(j)
            for j in range(5):
                for item in drug_list[j+1]:
                    if item not in drug_total_dict:
                        drug_total_dict[item] = [0,0,0,0,0]
                        drug_total_dict[item][j] += 1
                    else:
                        drug_total_dict[item][j] += 1
        
        drug_df = pd.DataFrame(drug_total_dict)
        ALL_DRUG_LIST = drug_df.columns
        
        SUM = (drug_df.sum()).tolist()

        List1 = tuple(list(np.asarray(drug_df.iloc[0,:].tolist())/sample_size))
        List2 = tuple(list(np.asarray(drug_df.iloc[1,:].tolist())/sample_size))
        List3A = tuple(list(np.asarray(drug_df.iloc[2,:].tolist())/sample_size))
        List3B = tuple(list(np.asarray(drug_df.iloc[3,:].tolist())/sample_size))
        List4 = tuple(list(np.asarray(drug_df.iloc[4,:].tolist())/sample_size))
        LEVEL_COLOR = ['#359744', '#286ea0', '#8a4a8d', '#b490b5','#3d3d3b']
        width = 0.7
        fig2 = plt.figure(figsize=(10,5))
        ax2 = fig2.add_axes([0,0,1,1])
        
        ax2.bar(ALL_DRUG_LIST, List1 ,  width, color=LEVEL_COLOR[0])
        ax2.bar(ALL_DRUG_LIST, List2 ,  width, bottom=List1,color=LEVEL_COLOR[1])
        ax2.bar(ALL_DRUG_LIST, List3A , width, bottom = np.array(List1)+np.array(List2),color=LEVEL_COLOR[2])
        ax2.bar(ALL_DRUG_LIST, List3B , width, bottom = np.array(List1)+np.array(List2)+np.array(List3A),color=LEVEL_COLOR[3])
        ax2.bar(ALL_DRUG_LIST, List4 ,  width, bottom = np.array(List1)+np.array(List2)+np.array(List3A)+np.array(List3B), color=LEVEL_COLOR[4])
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_color('#cac9c9')
        ax2.spines['left'].set_color('#cac9c9')
        plt.ylabel("Percentage", fontsize=LABEL_SIZE, fontweight='bold')
        ax2.tick_params(axis='x',direction='in', color='#cac9c9', length=0)
        ax2.tick_params(axis='y',direction='in', color='#cac9c9')
        plt.yticks(fontsize=LABEL_SIZE-4)
        ax2.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
        ax2.set_yticks(np.arange(0, max(SUM)/sample_size*1.25, 0.2))
        plt.xticks(color='#222222',rotation=90, fontsize=LABEL_SIZE-4,horizontalalignment='center',verticalalignment='top')#verticalalignment='bottom',
        ax2.legend(labels=level_list, fontsize=LABEL_SIZE-4, edgecolor='white')
        # ax2.set_yticks(np.arange(0, max(freq)*1.25, 0.1))

        # plt.title("Frequency of Actionable Genes", fontsize=20, fontweight='bold')
        plt.savefig(pic+"oncokb_freq_actionable_genes.pdf", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + pic + "oncokb_freq_actionable_genes.pdf"), 'green'))
