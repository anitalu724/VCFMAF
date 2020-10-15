############################################################################################
# FileName     [ maf_plot.py ]
# PackageName  [ VCFMAF ]
# Synopsis     [ Functions for visualization of MAF analysis ]
# Author       [ Cheng-Hua (Anita) Lu ]
# Copyright    [ 2020/9 ]
############################################################################################


from comut import comut
from comut import fileparsers
import sigProfilerPlotting as sigPlt

import pandas as pd
import numpy as np
import math
from termcolor import colored

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


class CoMutPlot:
    def __init__(self, file, info):
        print(colored(("\nStart plotting CoMut Plot...."), 'yellow'))
        self.fd = (pd.read_csv(file, sep="\t")).to_dict('list')
        for item in self.fd:
            cleanedList = [x for x in self.fd[item] if str(x) != 'nan']
            self.fd[item] = cleanedList
        self.info = (pd.read_csv(info, sep="\t")).to_dict('list')
        for item in self.info:
            print(item)
            cleanedList = [x for x in self.info[item] if str(x) != 'nan']
            self.info[item] = cleanedList
    def plot(self, folder, theme, name):
        fixed_category = ['Same Patient','Copy Number Alteration','Mutation Type','Purity','Mutation Signature','Mutation Clonality','Frequency']
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
            


        # Plot
        toy_comut = comut.CoMut()
        width, fig_info = 0, 0
        if not same_patient_df.empty:
            width += 1
            fig_info += 1
            toy_comut.add_sample_indicators(same_patient_df, name = fixed_category[0], plot_kwargs = indicator_kwargs) 
        if not copy_num_df.empty:
            width += 2
            fig_info += 7
            toy_comut.add_categorical_data(copy_num_df, name=fixed_category[1], mapping = cna_mapping, category_order = cna_order) 
        if not mut_type_df.empty:
            width += 5*int(len(self.info['Mutation Type'])/10)
            fig_info += 9
            toy_comut.add_categorical_data(mut_type_df, name = fixed_category[2], category_order = category_order, mapping = mut_mapping)
        if not purity_df.empty:
            width += 1
            toy_comut.add_continuous_data(purity_df, name = fixed_category[3], mapping = used_set['purity'], value_range = (0, 1))
        
        for i in range(len(feature_df)):
            width += 1
            toy_comut.add_categorical_data(feature_df[i], name=feature_category[i], mapping = feature_mapping[i])
        
        if not mut_sig_df.empty:
            width += 2
            toy_comut.add_bar_data(mut_sig_df, name = fixed_category[4], mapping = sig_mapping, stacked = True, ylabel = 'Mutational\nsignatures',bar_kwargs = sig_kwargs)
        if not mut_clone_df.empty:
            width += 2
            toy_comut.add_bar_data(mut_clone_df, name = fixed_category[5], mapping = bar_mapping, stacked = True, bar_kwargs = bar_kwargs, ylabel = 'Muts/Mb')
        

        FigCol = 2 if width > 10 else 1
        height = int(width/10)*8
        # print(width)
        FigSize = (10, height)
        # print(FigSize)
        
        # # purity_ax = toy_comut.figure.add_axes([1.07, 0.46, 0.08, 0.014])
        # norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        # purity_colorbar = toy_comut.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=used_set['purity']),
                                                #  cax=purity_ax, orientation='horizontal')
        
        # patient_leg = toy_comut.add_axis_legend('Same Patient', bbox_to_anchor = (1.025, 1.2), frameon = False, numpoints = 2)
        toy_comut.plot_comut(figsize = FigSize, x_padding = 0.04, y_padding = 0.04, tri_padding = 0.03)
        toy_comut.add_unified_legend(ncol = 2)
        toy_comut.figure.savefig(folder+name, dpi = 300, bbox_inches = 'tight')
        print(colored(("=> Generate CoMut Plot: "+folder+name), 'green'))


class CosineSimilarity:
    def __init__(self, file):
        print(colored(("\nStart plotting Cosine Similarity Heatmap...."), 'yellow'))
        self.df = (pd.read_csv(file, header = 0, index_col=0))
        self.x_label = list(self.df.columns)
        self.y_label = list(self.df.index)
    def plot(self, folder):
        height, length = len(self.y_label), len(self.x_label)
        h , l = int(height/4 + 1)*3, int(length/6 + 1)*3
        data = np.array(self.df.to_numpy())
        f, ax = plt.subplots(figsize=(l,h))
        ax = sns.heatmap(data, vmin=0, vmax=1, xticklabels =self.x_label, yticklabels = self.y_label, 
                               square=True, cmap="Blues",
                               cbar_kws={"orientation": "vertical",'shrink':0.5})
        ax.set_title('Cosine Similarity Plot', fontsize=20,weight='bold')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right',fontweight='light')
        plt.savefig(folder+"S2S.png",dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Cosine Similarity Plot: "+folder+"S2S.png"), 'green'))


class SBSPlot:
    def __init__(self, file, folder):
        print(colored(("\nStart plotting SBS Plot...."), 'yellow'))
        self.df = pd.read_csv(file)
        col = list(self.df.columns)
        col[0]='MutationType'
        self.df.columns = col
        (self.df).to_csv(folder+"SBS.tsv", sep="\t", index=False)     
    def plot(self, folder):
        sigPlt.plotSBS(folder+"SBS.tsv", folder,"", "96", percentage=True)
        print(colored(("=> Generate SBS Plot: "+folder+"SBS_96_plots_.pdf"), 'green'))


class DonutPlot:
    def __init__(self, file):
        print(colored(("\nStart plotting Donut Plot...."), 'yellow'))
        self.df = pd.read_csv(file, index_col=0)
    def plot(self, folder):
        raw_data = self.df.sum(axis=1)/self.df.shape[1]
        SUM = raw_data.sum(axis=0)
        raw_data = raw_data/SUM
        names, sizes = list(raw_data.index), list(raw_data.iloc[:])
        color_map = ['#b7d5ea','#acc6aa','#266199','#E0CADB','#695D73','#B88655','#DDDDDD','#71a0a5','#841D22']
        my_circle = plt.Circle( (0,0), 0.6, color='white')
        plt.pie(sizes, labels=names, colors=color_map[:len(names)], autopct='%1.1f%%',textprops={'fontsize': 12},radius=1, pctdistance=1-0.4/2)
        
        p = plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig(folder+"donut_plot.png", dpi=300)
        print(colored(("=> Generate Donut Plot: "+folder+"donut_plot.png"), 'green'))



 
        