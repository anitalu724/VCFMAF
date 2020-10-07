############################################################################################
# FileName     [ maf_analysis.py ]
# PackageName  [ src ]
# Synopsis     [ Use the entered MAF files to do some analysis ]
# Author       [ Cheng-Hua (Anita) Lu ]
# Copyright    [ 2020 9 ]
############################################################################################

from .maf_filter import *
import pandas as pd
import numpy as np
from termcolor import colored

import time
import os
import math
import seaborn as sns
# sns.set(color_codes = True)
import matplotlib.pyplot as plt
import matplotlib.style
from matplotlib.ticker import MaxNLocator

# 1. CoMut Plot Analysis
class CoMutAnalysis:
    def __init__(self, file):
        print(colored(("\nStart CoMut_Plot_Analysis...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder):
        def mutation_type():
            maf = self.df
            chosen_col = maf[['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification']]
            chosen_col = chosen_col.rename({'Tumor_Sample_Barcode':'sample', 'Hugo_Symbol':'category', 'Variant_Classification':'value'}, axis=1)
            value_old_list = ['Missense_Mutation', 
                              'Nonsense_Mutation',
                              'In_Frame_Del', 
                              'In_Frame_Ins',
                              'Splice_Site',
                              'Silent',
                              'Frame_Shift_Del',
                              'Frame_Shift_Ins',
                              'Nonstop_Mutation',
                              'Translation_Start_Site']
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
            chosen_col.to_csv(folder+"mutation_data.tsv", sep = "\t", index = False)
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

# 2. Significantly mutated gene detection (Only oncodriveCLUST works until now)
class SigMutatedGeneDetection:
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

# 3. Known cancer gene annotation
class KnownCancerGeneAnnotation:
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

# 4. Mutational signature
class MutationalSignature:
    def __init__(self, file):
        print(colored(("\nStart Mutational_Signature...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder, fig_folder):
        def get_input_file():
            output_file = folder+"ms_input.tsv"
            selected_col = self.df[['Tumor_Sample_Barcode','flanking_bps', 'Reference_Allele', 'Tumor_Seq_Allele2']]
            selected_col.columns = ['SampleID', 'Three_Allele', 'Ref', 'Mut']
            new_dict = {}
            for i in range(len(selected_col['SampleID'])):
                if selected_col['SampleID'][i] not in new_dict:
                    new_dict[selected_col['SampleID'][i]] = [selected_col.loc[i, :]]
                else:
                    new_dict[selected_col['SampleID'][i]].append(selected_col.loc[i, :])
            # rearrange each data 
            for sampleID in new_dict:
                sample_dict = {'C>A':[],'C>G':[],'C>T':[],'T>A':[],'T>C':[],'T>G':[]}#,'OTHER':[]}
                for item in new_dict[sampleID]:
                    if (item['Ref'] == 'C' and item['Mut'] == 'A') or (item['Ref'] == 'G' and item['Mut'] == 'T'):
                        sample_dict['C>A'].append(item)
                    elif (item['Ref'] == 'C' and item['Mut'] == 'G') or (item['Ref'] == 'G' and item['Mut'] == 'C'):
                        sample_dict['C>G'].append(item)
                    elif (item['Ref'] == 'C' and item['Mut'] == 'T') or (item['Ref'] == 'G' and item['Mut'] == 'A'):
                        sample_dict['C>T'].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'A') or (item['Ref'] == 'A' and item['Mut'] == 'T'):
                        sample_dict['T>A'].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'C') or (item['Ref'] == 'A' and item['Mut'] == 'G'):
                        sample_dict['T>C'].append(item)
                    elif (item['Ref'] == 'T' and item['Mut'] == 'G') or (item['Ref'] == 'A' and item['Mut'] == 'C'):
                        sample_dict['T>G'].append(item)
                list_96 = []
                for data in sample_dict:
                    # data = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
                    if data != "OTHER":
                        new_list = [int(0)]*16
                        three_allele_dict = {}
                        if data in ['C>A', 'C>G', 'C>T']:
                            three_allele_dict={"ACA":0,     "TGT":0,    "ACC":1,    "TGG":1,    "ACG":2,    "TGC":2,    "ACT":3,    "TGA":3, \
                                               "CCA":4,     "GGT":4,    "CCC":5,    "GGG":5,    "CCG":6,    "GGC":6,    "CCT":7,    "GGA":7, \
                                               "GCA":8,     "CGT":8,    "GCC":9,    "CGG":9,    "GCG":10,   "CGC":10,   "GCT":11,   "CGA":11,\
                                               "TCA":12,    "AGT":12,   "TCC":13,   "AGG":13,   "TCG":14,   "AGC":14,   "TCT":15,   "AGA":15 }   
                        elif data in ['T>A', 'T>C', 'T>G']:
                            three_allele_dict={"ATA":0,     "TAT":0,    "ATC":1,    "TAG":1,    "ATG":2,    "TAC":2,    "ATT":3,    "TAA":3, \
                                               "CTA":4,     "GAT":4,    "CTC":5,    "GAG":5,    "CTG":6,    "GAC":6,    "CTT":7,    "GAA":7, \
                                               "GTA":8,     "CAT":8,    "GTC":9,    "CAG":9,    "GTG":10,   "CAC":10,   "GTT":11,   "CAA":11,\
                                               "TTA":12,    "AAT":12,   "TTC":13,   "AAG":13,   "TTG":14,   "AAC":14,   "TTT":15,   "AAA":15 }   
                        for j in sample_dict[data]:
                            # print(j['Three_Allele'])
                            if j['Three_Allele'] in three_allele_dict:
                                new_list[three_allele_dict[j['Three_Allele']]]+=1
                        for each_num in new_list:
                            list_96.append(each_num)
                new_dict[sampleID] = list_96
            new_df = pd.DataFrame.from_dict(new_dict)
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
        def run_r_estimation():
            input_file = folder+"ms_input.tsv"
            code_file = fig_folder+"estimation.r"
            img_file = fig_folder+"1_Estimation.png"
            path = os.path.abspath(input_file)
            img_path = os.path.abspath(img_file)
            code = open(code_file, "w")
            code.write("library(BSgenome)\nlibrary(MutationalPatterns)\nlibrary(NMF)\n")
            code.write("file_path <- "+"\""+path+"\"\n")
            code.write("data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)\n")
            # 1_Estimate
            code.write("png(filename=\""+img_path+"\",width = 1000, height = 1000)\n")
            code.write("data <- data+0.0001\n")
            code.write("estimate <- nmf(data, rank=2:5, method=\"brunet\", nrun=10, seed=123456)\n")
            code.write("plot(estimate)\ndev.off()\n")
            #  os.system("bash Install_script/install_r_estimation.sh")
            os.popen("Rscript "+ code_file+" | head\n")
            print(colored("=> Generate estimation figure: ", 'green'))
            print(colored(("   "+img_file), 'green'))
        def run_r_signature():
            input_file = folder+"ms_input.tsv"
            code_file = fig_folder+"signature.r"
            img_file = fig_folder+"1_Estimation.png"
            img_file1 = fig_folder+"2_Signature_96.png"
            img_file2 = fig_folder+"3_Signature_Distribution.png"
            img_heat = fig_folder+"4_Signature2sample_Heatmap.png" 
            img_compare = fig_folder+"5_Similarity_Comparison.png"
            img_sig_heat = fig_folder+"6_Signature_Similarity_heatmap.png"
            c_path = fig_folder+"nmf_res_contribution.csv"
            s_path = fig_folder+"nmf_res_signature.csv"
            m_path = fig_folder+"S2S_matrix.csv"
            rank = "4"#input("Input number of signatures: ")
            sig_list = [("Signature_"+str(i)) for i in range(int(rank))]
            new_str = ""
            for i in range(len(sig_list)):
                new_str+=("\""+sig_list[i]+"\"")
                if i != (len(sig_list)-1):
                    new_str+=","
            code = open(code_file, "w")
            code.write("library(BSgenome)\nlibrary(MutationalPatterns)\nlibrary(NMF)\nlibrary(\"gridExtra\")\n")
            code.write("file_path <- "+"\""+os.path.abspath(input_file)+"\"\n")
            code.write("data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)\n")
            code.write("data <- data+0.0001\n")
            # 2_Signature_96
            code.write("nmf_res <- extract_signatures(data, rank = "+rank+", nrun = 10)\n")
            code.write("colnames(nmf_res$signatures) <- c("+new_str+")\n")
            code.write("rownames(nmf_res$contribution) <- c("+new_str+")\n")
            code.write("png(filename=\""+os.path.abspath(img_file1)+"\",width = 1000, height = 1000)\n")
            code.write("plot_96_profile(nmf_res$signatures, condensed = TRUE)\ndev.off()\n")
            # 3_Signature_distribution
            code.write("pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = \"relative\")\n")
            code.write("pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = \"absolute\")\n")
            code.write("png(filename=\""+os.path.abspath(img_file2)+"\",width = 1000, height = 1000)\n")
            code.write("a <- grid.arrange(pc1, pc2)\nplot(a)\ndev.off()\n")
            # 4_Signature2sample_heatmap
            code.write("png(filename=\""+os.path.abspath(img_heat)+"\",width = 1000, height = 1000)\n")
            code.write("pch1 <- plot_contribution_heatmap(nmf_res$contribution,sig_order = c("+new_str+"))\n")
            code.write("pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)\n")
            code.write("grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6))\ndev.off()\n")
            # 5_Similarity_Comparison
            code.write("png(filename=\""+os.path.abspath(img_compare)+"\",width = 1000, height = 1000)\n")
            code.write("plot_compare_profiles(data[,1],nmf_res$reconstructed[,1],profile_names = c(\"Original\", \"Reconstructed\"),condensed = TRUE)\ndev.off()\n")
            # 6_Signature_Similarity_heatmap
            code.write("sp_url <- paste(\"https://cancer.sanger.ac.uk/cancergenome/assets/\",\"signatures_probabilities.txt\", sep = \"\")\n")
            code.write("cancer_signatures = read.table(sp_url, sep = \"\t\", header = TRUE)\n")
            code.write("new_order = match(row.names(data), cancer_signatures$Somatic.Mutation.Type)\n")
            code.write("cancer_signatures = cancer_signatures[as.vector(new_order),]\n")
            code.write("row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type\n")
            code.write("cancer_signatures = as.matrix(cancer_signatures[,4:33])\n")
            code.write("hclust_cosmic = cluster_signatures(cancer_signatures, method = \"average\")\n")
            code.write("cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]\n")
            code.write("cos_sim_samples_signatures = cos_sim_matrix(nmf_res$signature, cancer_signatures)\n")
            code.write("png(filename=\""+os.path.abspath(img_sig_heat)+"\",width = 800, height = 600)\n")
            code.write("plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)\ndev.off()\n")
            # back_up references
            code.write("write.csv(nmf_res$contribution,file=\""+os.path.abspath(c_path)+"\",sep=\",\",row.names=T, na = \"NA\")\n")
            code.write("write.csv(nmf_res$signature,file=\""+os.path.abspath(s_path)+"\",sep=\",\",row.names=T, na = \"NA\")\n")
            code.write("write.csv(cos_sim_matrix(nmf_res$signatures,cancer_signatures),file=\""+os.path.abspath(m_path)+"\",sep=\",\",row.names=T, na = \"NA\")\n")
            os.popen("Rscript "+ code_file+" | head\n")
            print(colored("=> Generate signature figures: ", 'green'))
            print(colored(("   "+img_file1), 'green'))
            print(colored(("   "+img_file2), 'green'))
            print(colored(("   "+img_heat), 'green'))
            print(colored(("   "+img_compare), 'green'))
            print(colored(("   "+img_sig_heat), 'green'))
            print(colored("=> Generate NMF references CSV files:", 'green'))
            print(colored(("   "+c_path), 'green'))
            print(colored(("   "+s_path), 'green'))
            print(colored(("   "+m_path), 'green'))
        get_input_file()
        run_r_estimation()
        run_r_signature()


# 5. Mutation burden statistics
class TotalMutationBurden:
    def __init__(self, file):
        print(colored(("\nStart Total Mutation Burden...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder, mode, length):
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


# 6. OncoKB annotator
class OncoKBAnnotator:
    def __init__(self, file):
        print(colored(("\nStart OncoKB annotator(drug)...."), 'yellow'))
        self.head, self.df = fast_read_maf(file)
    def data_analysis(self, folder, path, token, clinical):
        selected_df = (self.df[['Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'HGVSp']]).set_index("Hugo_Symbol")
        selected_df.to_csv(folder + "maf_5col_oncokb_input.txt", sep="\t")
        os.system("python3 " + path + "MafAnnotator.py -i " + folder + "maf_5col_oncokb_input.txt -o " + folder + "maf_5col_oncokb_output.txt -b " + token + "\n")
        os.system("python3 " + path + "ClinicalDataAnnotator.py -i "+clinical+" -o "+folder+"clinical_oncokb_output.txt -a "+folder+"maf_5col_oncokb_output.txt")
        os.system("rm "+folder+"maf_5col_oncokb_input.txt\n")
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   "+folder+"maf_5col_oncokb_output.txt"), 'green'))
        print(colored(("   " + folder + "clinical_oncokb_output.txt"), 'green'))
    def plotting(self, folder, level='4'):
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

        size = [true_num, sample_size - true_num]
        labels = "Actionable\nbiomarkers","Current absence\nof new actionable\nbiomarkers"
        fig1, ax1 = plt.subplots()
        ax1.pie(size, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#b7d5ea','#266199'] ,textprops={'fontsize': 11})
        ax1.axis('equal')
        plt.title("Total", fontsize=18, fontweight='bold')
        plt.savefig(folder+"oncokb_total_pie.png", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + folder + "oncokb_total_pie.png"), 'green'))
        df_drug_count = df[level_list]
        drug_total_dict = {}
        for i in range(sample_size):
            for item in level_list:
                data = df_drug_count.iloc[i][item]
                if not pd.isna(data):
                    new_drug_list = data.split("(")
                    new_drug_list.pop(0)
                    for drug in new_drug_list:
                        drug_name = drug.split(" ")[0]
                        if drug_name not in drug_total_dict:
                            drug_total_dict[drug_name] = 1
                        else:
                            drug_total_dict[drug_name] += 1
        langs = list(drug_total_dict.keys())
        count = list(drug_total_dict.values())
        fig2 = plt.figure(figsize=(8,5))
        ax2 = fig2.add_axes([0,0,1,1])
        ax2.bar(langs, count, color='#266199')
        plt.xticks(fontsize=14)
        ax2.set_yticks(np.arange(0, max(count) + 1, 5))
        plt.yticks(fontsize=14)
        plt.ylabel("Count")

        plt.title("Frequency of Actionable Genes", fontsize=18, fontweight='bold')
        plt.savefig(folder+"oncokb_actionable_genes.png", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + folder + "oncokb_actionable_genes.png"), 'green'))


class HRDScore:
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
    def plotting(self, folder):
        df = pd.read_csv(folder+"all_HRDresults.csv")
        size = df.shape[0]
        HRD = tuple(list(df['HRD']))
        TAI = tuple(list(df['Telomeric_AI']))
        LST = tuple(list(df['LST']))
        SUM = list(df["HRD-sum"])
        Sample = tuple(list(df['Sample_id']))
        ind = np.arange(size)
        width = 0.8
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_axes([0,0,1,1])
        ax.bar(ind, HRD, width, color='#266199')
        ax.bar(ind, TAI, width,bottom=HRD, color='#b7d5ea')
        ax.bar(ind, LST, width,bottom=np.array(TAI)+np.array(HRD), color='#acc6aa')
        ax.set_ylabel('Scores')
        ax.set_title('HRD Scores',fontsize=18, fontweight='bold')
        plt.xticks(ind, Sample,rotation=45,horizontalalignment='right',fontweight='light')
        ax.set_yticks(np.arange(0, max(SUM)+3, 10))
        ax.legend(labels=['HRD','Telomeric_AI','LST'])
        plt.savefig(folder+"HRD_Score.png", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + folder + "HRD_Score.png"), 'green'))
        over = len([i for i in SUM if i >=42])
        data = [over, size-over]
        labels='≧ 42','< 42'
        fig1, ax1 = plt.subplots()
        ax1.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#266199','#b7d5ea'] ,textprops={'fontsize': 11})
        ax1.axis('equal')
        plt.title("Propotion of Samples with HRD Phenotype", fontsize=18, fontweight='bold')
        plt.savefig(folder+"high_HRD_pie.png", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + folder + "high_HRD_pie.png"), 'green'))


class WGDnCIN:
    def __init__(self, file):
        print(colored(("\nStart analysing WGD and CIN...."), 'yellow'))
        self.list = ((pd.read_csv(file, sep="\t"))[['CNV_input']].values.T)[0]
    def data_analysis(self, folder):
        genome_length = [249250621, 243199373, 198022430, 191154276, 180915260,
                         171115067, 159138663, 146364022, 141213431, 135534747,
                         135006516, 133851895, 115169878, 107349540, 102531392,
                         90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
        whole_length = sum(genome_length)
        WGD_dict, CIN_dict = {}, {}
        for sample in self.list:
            df = pd.read_csv(sample, sep="\t")
            sample_name = pd.unique(df['SampleID'])[0]
            WGD_selected = df.loc[(df['A_cn'] >= 2)]
            CIN_selected = df.loc[(df['A_cn'] != 1) & (df['B_cn'] != 1)]
            WGD_SUM, CIN_SUM = 0, 0
            for i in range(WGD_selected.shape[0]):
                WGD_SUM += WGD_selected.iloc[i]['End_position']-WGD_selected.iloc[i]['Start_position']
            for i in range(CIN_selected.shape[0]):
                CIN_SUM += CIN_selected.iloc[i]['End_position']-CIN_selected.iloc[i]['Start_position']
            WGD_dict[sample_name] = (WGD_SUM >= 0.5*whole_length)
            CIN_dict[sample_name] = CIN_SUM/whole_length
        WGD_df = pd.DataFrame.from_dict(WGD_dict, orient='index').reset_index()
        CIN_df = pd.DataFrame.from_dict(CIN_dict, orient='index').reset_index()
        WGD_df.columns,CIN_df.columns = [['SampleID','WGD']], [['SampleID','CIN']]
        WGD_df.to_csv(folder + "WGD_result.csv",  index=False)
        CIN_df.to_csv(folder + "CIN_result.csv",  index=False)
        print(colored("=> Generate analysis files: ", 'green'))
        print(colored(("   " + folder + "WGD_result.csv"), 'green'))
        print(colored(("   " + folder + "CIN_result.csv"), 'green'))
    def plotting(self, folder):
        # WGD Pie Plot
        wgd_df = pd.read_csv(folder+"WGD_result.csv")
        data = [0,0]
        for i in range(wgd_df[['WGD']].shape[0]):
            if wgd_df[['WGD']].iloc[i]['WGD'] == False:
                data[1]+=1
            elif wgd_df[['WGD']].iloc[i]['WGD'] == True:
                data[0]+=1
        labels = 'True','False'
        fig1, ax1 = plt.subplots()
        ax1.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#266199','#b7d5ea'] ,textprops={'fontsize': 11})
        ax1.axis('equal')
        plt.title("Propotion of Samples with WGD", fontsize=18, fontweight='bold')
        plt.savefig(folder+"WGD_pie.png", dpi=300, bbox_inches='tight')
        print(colored(("=> Generate Pie Plot: " + folder + "WGD_pie.png"), 'green'))
        
        # CIN Bar plot
        CIN_df = pd.read_csv(folder+"CIN_result.csv")

        size = CIN_df.shape[0]
        CIN = tuple(list(CIN_df['CIN']))
        Sample = tuple(list(CIN_df['SampleID']))
        ind = np.arange(size)
        width = 0.8
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_axes([0,0,1,1])
        ax.bar(ind, CIN, width, color='#266199')
        ax.set_ylabel('Scores')
        ax.set_title('CIN Scores',fontsize=18, fontweight='bold')
        plt.xticks(ind, Sample,rotation=45,horizontalalignment='right',fontweight='light')
        ax.set_yticks(np.arange(0, 1, 0.2))
        plt.savefig(folder+"CIN_Score.png", dpi=300,bbox_inches='tight')
        print(colored(("=> Generate Bar Plot: " + folder + "CIN_Score.png"), 'green'))



            

