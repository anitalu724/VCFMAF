import pandas as pd

file = "sample_data/20200930/CBCP00001-T1_HRDresults.txt"
df = pd.read_csv(file, sep="\t")
print(df)