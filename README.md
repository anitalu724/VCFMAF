# MyLabTool
MyLabTool is a Python tool for preprocessing variant call format(VCF) and mutation annotation format(MAF) files and 

## Prerequisite
```shell
pip install termcolor
pip install tqdm
pip install numpy
pip install pandas
pip install vcfpy
pip install seaborn
pip install oncodriveclust
pip install comut
pip install SigProfilerPlotting
```
## Preprocessing VCF files and MAF files
### Required input files
* a TSV file
   * for VCF files: 9 columns					
| NORMAL | TUMOR | MuSe | Mutect2 | SomaticSniper | Strelka2 | VarScan2 | At Least # CALLS | At Most # REJECT |
| ------ | ----- | ---- | ------- | ------------- | -------- | -------- | ---------------- | ---------------- |
| ...    | ...   | ...  | ...     | ...           | ...      | ...      | ...              | ...              |
  * for MAF files: 1 column 
| MAF |
| --- |
| ... |




# VCFMAF

## Data Preprocessing

####
* My oncokb token: ca398551-c549-49bd-80f0-0e68d9ca033c