
######################################################
######## Volcano Plots from DESeq2 #######
######################################################
#### returning to default p-adjust of <1e-5 (not 1e-3)
DEseq_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/2-normalization/DESeq2/v1_SLR/DE_csvs"
wald_path <- file.path(DEseq_dir)


##########################
###### infection
##########################
csv <- "cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)



pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-5, FCcutoff = 1, widthConnectors = 0.25, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20, gridlines.minor = FALSE, labSize=3, boxedLabels=FALSE, drawConnectors = TRUE)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-5, FCcutoff = 1, widthConnectors = 0.25, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20, gridlines.minor = FALSE, labSize=3, boxedLabels=FALSE, drawConnectors = TRUE)
dev.off()



##########################
##### interaction infection x cell type
##########################
csv <- "cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)


####################
###### cell type
###################
csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)
