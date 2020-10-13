library(BSgenome)
library(MutationalPatterns)
library(NMF)
file_path <- "/home/b06901011/mskcc-vcf2maf-bbe39fe/VCFMAF/sample_data/out_20201013/ms_input.tsv"
data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
png(filename="/home/b06901011/mskcc-vcf2maf-bbe39fe/VCFMAF/sample_data/pic/1_Estimation.png",width = 1000, height = 1000)
data <- data+0.0001
estimate <- nmf(data, rank=2:5, method="brunet", nrun=10, seed=123456)
plot(estimate)
dev.off()
