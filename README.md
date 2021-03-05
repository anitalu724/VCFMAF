###### tags: `Lab`
# MutScape
MutScape is a Python tool for preprocessing variant call format(VCF) and mutation annotation format(MAF) files and 

## Prerequisite
```bash=
pip install termcolor
pip install tqdm
pip install numpy
pip install pandas
pip install vcfpy
pip install seaborn
pip install oncodriveclust
pip install comut
```
______

## Preprocessing VCF files and MAF files
```bash=
python3 vcf_maf_process.py ...
```
### Required input files(Required)
```bash=
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
```bash=
python3 vcf_maf_process.py \
...
-vf GI "[1:3,5]" CI "15,15,0,0,0,0.05,8,8" P 1 FFPE 0.9 \
...
```
#### VCF combination
```bash=
python3 vcf_maf_process.py \
...
-c \
...
```
#### VCF transform to MAF 
```bash=
python3 vcf_maf_process.py \
...
-v2m 48 \
...
```
### For MAF preprocessing
#### MAF filtering
```bash=
python3 vcf_maf_process.py \
...
-mf GI "[1:3,5]"
...
```
### Output files and Meta files
```bash=
python3 vcf_maf_process.py \
...
-o [output files path] \
-m [meta files path] \
...
```
______

## Data Analysis and Visualization
### CoMut Plot Analysis
### Mutational Signature

1. Preprocessing: Make all MAFs into one MAF file
```bash=
python3 vcf_maf_process.py \
-f examples/Tissue_samples/ms_maf/maf.tsv \
-m examples/Tissue_samples/ms_maf \
-o examples/Tissue_samples/ms_maf
```
2. Estimation
```bash=
python3 maf_analysis.py \
-f examples/Tissue_samples/ms_maf/maf_combination.maf \
-ms 1 "[2,9,10]" \
-o examples/Tissue_samples/outputs \
-p examples/Tissue_samples/pictures
```
3. Analysis and Visualization
```bash=
python3 maf_analysis.py \
-f examples/Tissue_samples/ms_maf/maf_combination.maf \
-ms 2 "[3]" \
-o examples/Tissue_samples/outputs \
-p examples/Tissue_samples/pictures
```
4. Outputs
![](https://i.imgur.com/mhWyWjf.png)


### HRD Score, Whole-Genome doubling(WGD) and Chromosome instability(CIN)

1. HRD Score
```bash=
python3 maf_analysis.py \
-f examples/CBCP_samples/CBCP_oncokb_maf.txt \
-hrd examples/CBCP_samples/HRD_input.tsv grch37 \
-o examples/CBCP_samples/outputs \
-p examples/CBCP_samples/pictures
```

2. WGD_CIN
```bash=
python3 maf_analysis.py \
-f examples/CBCP_samples/CBCP_oncokb_maf.txt \
-wgdcin examples/CBCP_samples/HRD_input.tsv \
-o examples/CBCP_samples/outputs \
-p examples/CBCP_samples/pictures
```

3. ONCOKB
```bash=
python3 maf_analysis.py \
-f examples/TCGA/TCGA_Mutect_v10_white_rerun.maf \
-oncokb ../oncokb-annotator/ ca398551-c549-49bd-80f0-0e68d9ca033c 4 examples/TCGA/clinical_input.txt \
-o examples/TCGA/outputs \
-p examples/TCGA/pictures
```

4. Outputs examples
![](https://i.imgur.com/hSv3VsM.png)




###

______
####
* My oncokb token: ca398551-c549-49bd-80f0-0e68d9ca033c
