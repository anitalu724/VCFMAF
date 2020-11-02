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
```shell
python3 vcf_maf_process.py ...
```
### Required input files(Required)
```shell
python3 vcf_maf_process.py \
-f [tsv file path]
```
- for VCF files: a 9 columns TSV file

| NORMAL | TUMOR | MuSe | Mutect2 | SomaticSniper | Strelka2 | VarScan2 | At Least # CALLS | At Most # REJECT |
| ------ | ----- | ---- | ------- | ------------- | -------- | -------- | ---------------- | ---------------- |
| ...    | ...   | ...  | ...     | ...           | ...      | ...      | ...              | ...              |
- for MAF files: a 1 column TSV file 

| MAF |
| --- |
| ... |

### For VCF preprocessing
#### VCF filtering
```shell
python3 vcf_maf_process.py \
...
-vf GI "[1:3,5]" \
...
```
#### VCF combination
```shell
python3 vcf_maf_process.py \
...
-c \
...
```
#### VCF transform to MAF 
```shell
python3 vcf_maf_process.py \
...
-v2m 48 \
...
```
### For MAF preprocessing
#### MAF filtering
```shell
python3 vcf_maf_process.py \
...
-mf GI "[1:3,5]"
...
```
### Output files and Meta files
```shell
python3 vcf_maf_process.py \
...
-o [output files path] \
-m [meta files path] \
...
```



# VCFMAF

## Data Preprocessing

####
* My oncokb token: ca398551-c549-49bd-80f0-0e68d9ca033c