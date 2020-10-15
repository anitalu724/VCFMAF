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


        