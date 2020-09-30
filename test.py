
import pandas as pd
origin = "sample_data/20200929/maf_combination.maf"
five_col = "sample_data/20200929/maf_combination_5col.maf"
oncokb_input = "sample_data/20200928_test/maf_combination_oncokb_input.maf"
file = "sample_data/20200928_test/maf_combination.oncokb.maf"
output = 'sample_data/20200928_test/maf.oncokb.txt'


df = pd.read_csv(origin, sep="\t", skiprows=1, header=0, encoding='gb2312', low_memory=False)
# df = pd.read_csv(origin, sep="\t")
output_df = df[['Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'HGVSp']]
output_df = output_df.set_index("Hugo_Symbol")
output_df.to_csv(five_col, sep="\t")


# selected = df[['Hugo_Symbol','Tumor_Sample_Barcode','LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4', 'LEVEL_R1', 'LEVEL_R2', 'LEVEL_R3', 'HIGHEST_LEVEL']]
# selected = selected.set_index(['Hugo_Symbol','Tumor_Sample_Barcode'])
# selected = selected.dropna(how="all")
# selected.to_csv("sample_data/20200928_test/maf_oncokb.txt",sep="\t")


