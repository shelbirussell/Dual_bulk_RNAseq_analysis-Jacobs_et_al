library(tximport)
library(DESeq2)
library(ggplot2)
library(rhdf5)
library(dplyr)
library(pheatmap)
library(EnhancedVolcano)
library(apeglm)
library("IHW")

## re-running the analysis on my kallisto output

kallisto_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/1-transcript_counting/kallisto/cds.gRNA_ref/pseudoaligned"
base_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/1-transcript_counting/kallisto/cds.gRNA_ref"
sample_id <- dir(kallisto_dir)

sample_dirs <- file.path(kallisto_dir, sample_id)

metadata_path <- file.path(base_dir,"samples.txt")
metadata <- read.table(metadata_path, header = TRUE)

##reassign "condition" from character to factor
metadata$condition <- as.factor(metadata$condition)

## append directories to table in a new column; label column "path"
metadata <- dplyr::mutate(metadata, path = sample_dirs)

################################################################################
## Transcripts need to be associated with gene IDs for gene-level summarization.
## If that information is present in the files, we can skip this step.
## For Salmon, Sailfish, and kallisto the files only provide the transcript ID.
## So, geneIDs need to be associated with transcript IDs here
################################################################################
## We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID.
## The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

ref_db <- file.path(base_dir, "Dmel6_cds.rna_from_genomic_gene_ids.tsv")
t2g <- read.table(ref_db, header=FALSE)

### create a named vector pointing to the (subset) Kallisto transcriptome quantification files
files <- file.path(kallisto_dir, metadata$run, "abundance_Drosophila.tsv")
names(files) <- metadata$sample
all(file.exists(files))

## run kallisto with TSV files - hd5 files are tricky to subset
## https://support.bioconductor.org/p/9135208/ needed the dropInReps flag
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = t2g, dropInfReps = TRUE, ignoreAfterBar = TRUE)

##Note: importing `abundance.h5` is typically faster than `abundance.tsv`
##reading in files with read_tsv
##1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
##transcripts missing from tx2gene: 27  ==> what are these???
##summarizing abundance
##summarizing counts
##summarizing length


###################################################
################# GO ON FOR DEseq2 DE analysis --- STOP HERE and take txi.kallisto.tsv to WGCNA.R for WGCNA
###################################################

################################################################################
################################################################################
### Create a DESeqDataset ####
## use the gene-level estimated counts from the quantification tools, and additionally to use the transcript-level abundance estimates
## to calculate a gene-level offset that corrects for changes to the average transcript length across samples
## the function DESeqDataSetFromTximport takes care of creation of the offset for you

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = metadata, design = ~celltype + infection + celltype:infection)

##Warning message:
##In DESeqDataSet(se, design = design, ignoreRank) :
##  some variables in design formula are characters, converting to factors

## The ddsTxi object here can then be used as dds in the following analysis steps.
#######################
## prefilter low-count/coverage genes across samples
## keep rows that have at least 10 reads total --> n=6 (all samples)

smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## normalization happens automatically, this command is not necessary (unless making a gct file) #dds <- estimateSizeFactors(dds)
## Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.

dds$group <- factor(paste0(dds$celltype, dds$infection))

############################################################################
#### set "controls" - JW18 vs S2 (control) and wMel vs uninfected (control)
############################################################################
## run relevel before running DESeq() - can also use "contrast" later

dds$infection <- relevel(dds$infection, ref = "uninfected")
dds$celltype <- relevel(dds$celltype, ref = "S2")

################################################################################################################################
################################################################################################################################
################################################################################################################################
########################################## DIFFERENTIAL EXPRESSION ANALYSIS ####################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
### Fit statistical model
## DESeq2 offers two kinds of hypothesis tests:
## the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero,
## and the likelihood ratio test (LRT).
###############################################
#########################################

## Model and interaction design
## two cell types: S2 and JW18
## two infections: uninfected and wMel

################################################################################
################################################################################
################################################################################
##################################################################
########################### Wald test ############################
##################################################################

#################################################################
########## WALD formatting: retain the terms to test ############
#################################################################

## interactions modeled - Wald test
design(dds) <- ~ celltype + infection + celltype:infection

dds <- DESeq(dds)
resultsNames(dds)

## before re-leveling
#[1] "Intercept"                    "celltype_JW18_vs_S2"
#[3] "infection_wMel_vs_uninfected" "celltypeJW18.infectionwMel"

## after re-leveling
#[1] "Intercept"                    "celltype_JW18_vs_S2"
#[3] "infection_wMel_vs_uninfected" "celltypeJW18.infectionwMel"


#################################################################
## infection effect
outputI <- results(dds, name="infection_wMel_vs_uninfected")
#outputI <- results(dds, contrast=c("infection", "wMel", "uninfected"))
write.csv(outputI, file="cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.csv")



## Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
outputIlfc <- lfcShrink(dds, coef="infection_wMel_vs_uninfected", type="apeglm")
write.csv(outputIlfc, file="cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.lfc.csv")

summary(outputI)
summary(outputIlfc)

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4178, 39%
LFC < 0 (down)     : 3623, 33%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(outputIlfc)

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4178, 39%
LFC < 0 (down)     : 3623, 33%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

########

## Independent hypothesis weighting = IHW
##  p value filtering is to weight hypotheses to optimize power. A Bioconductor package, IHW, is available that implements the method of Independent Hypothesis Weighting (Ignatiadis et al. 2016).
outputIihw <- results(dds, name="infection_wMel_vs_uninfected", filterFun=ihw)
summary(outputIihw)
sum(outputIihw$padj < 0.1, na.rm=TRUE)
metadata(outputIihw)$ihwResult

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4241, 39%
LFC < 0 (down)     : 3612, 33%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting

## sum ## [1] 7853

#ihwResult object with 10839 hypothesis tests
#Nominal FDR control level: 0.1
#Split into 7 bins, based on an ordinal covariate

write.csv(outputIihw, file="cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.ihw.csv")

###
## Independent hypothesis weighting (IHW) + Shrinkage of effect size (LFC estimates)
# Replace log2FoldChange, lfcSE, and p-values in outputIlfc with those from cell_model
outputIihwlfc <- outputIihw
outputIihwlfc$log2FoldChange <- outputIlfc$log2FoldChange
outputIihwlfc$lfcSE <- outputIlfc$lfcSE
outputIihwlfc$pvalue <- outputIlfc$pvalue
outputIihwlfc$padj <- outputIlfc$padj

write.csv(outputIihwlfc, file="cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.ihw.lfc.csv")

###
## In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
## Points will be colored red(?) if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

pdf("cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf_MAplot.ihw.pdf")
plotMA(outputIihw, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf_MAplot.pdf")
plotMA(outputI, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf_MAplot.lfc.pdf")
plotMA(outputIlfc, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf_MAplot.ihw.lfc.pdf")
plotMA(outputIihwlfc, ylim=c(-2,2))
dev.off()



#################################################################
## cell type effect     celltype_JW18_vs_S2

outputC <- results(dds, name="celltype_JW18_vs_S2")
write.csv(outputC, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv")

outputClfc <- lfcShrink(dds, coef="celltype_JW18_vs_S2", type="apeglm")
write.csv(outputClfc, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.lfc.csv")

summary(outputC)
summary(outputClfc)

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4621, 43%
LFC < 0 (down)     : 4494, 41%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(outputClfc)

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4621, 43%
LFC < 0 (down)     : 4494, 41%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


###############################
## ordered by contrast function - not needed because factors are reordered
#outputC <- results(dds, contrast=c("celltype", "JW18", "S2"))
#write.csv(outputC, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv")


########

## Independent hypothesis weighting = IHW
outputCihw <- results(dds, name="celltype_JW18_vs_S2", filterFun=ihw)
summary(outputCihw)
sum(outputCihw$padj < 0.1, na.rm=TRUE)
metadata(outputCihw)$ihwResult

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4638, 43%
LFC < 0 (down)     : 4499, 42%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting

## sum ## [1] 9137

#ihwResult object with 10839 hypothesis tests
#Nominal FDR control level: 0.1
#Split into 7 bins, based on an ordinal covariate

write.csv(outputCihw, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.ihw.csv")

###
## Independent hypothesis weighting (IHW) + Shrinkage of effect size (LFC estimates)
# Replace log2FoldChange, lfcSE, and p-values in outputClfc with those from cell_model
outputCihwlfc <- outputCihw
outputCihwlfc$log2FoldChange <- outputClfc$log2FoldChange
outputCihwlfc$lfcSE <- outputClfc$lfcSE
outputCihwlfc$pvalue <- outputClfc$pvalue
outputCihwlfc$padj <- outputClfc$padj

## these can't be merged until I can reorder the levels upstream
write.csv(outputCihwlfc, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.ihw.lfc.csv")

###
## In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
## Points will be colored red(?) if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

pdf("cells_Kallisto-DESeq2_WaldTest-JW18.S2.celltype_MAplot.ihw.pdf")
plotMA(outputCihw, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-JW18.S2.celltype_MAplot.pdf")
plotMA(outputC, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-JW18.S2.celltype_MAplot.lfc.pdf")
plotMA(outputClfc, ylim=c(-2,2))
dev.off()

#pdf("cells_Kallisto-DESeq2_WaldTest-celltype_MAplot.ihw.lfc.pdf")
#plotMA(outputCihwlfc, ylim=c(-2,2))
#dev.off()



########################################################

## infection*cell type interaction effect       celltypeJW18.infectionwMel

outputIC <- results(dds, name="celltypeJW18.infectionwMel")
write.csv(outputIC, file="cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.csv")

outputIClfc <- lfcShrink(dds, coef="celltypeJW18.infectionwMel", type="apeglm")
write.csv(outputIClfc, file="cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.lfc.csv")

summary(outputIC)
summary(outputIClfc)

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3865, 36%
LFC < 0 (down)     : 4020, 37%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> summary(outputIClfc)

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3863, 36%
LFC < 0 (down)     : 4022, 37%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


###################################
###### ORDERED by contrast function
####################################
#outputIC <- results(dds, contrast=c("group", "JW18wMel", "S2uninfected"))
#write.csv(outputIC, file="cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.csv")


######################################################
## Independent hypothesis weighting = IHW
##  p value filtering is to weight hypotheses to optimize power. A Bioconductor package, IHW, is available that implements the method of Independent Hypothesis Weighting (Ignatiadis et al. 2016).
outputICihw <- results(dds, name="celltypeJW18.infectionwMel", filterFun=ihw)
summary(outputICihw)
sum(outputICihw$padj < 0.1, na.rm=TRUE)
metadata(outputICihw)$ihwResult

out of 10839 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3906, 36%
LFC < 0 (down)     : 4056, 37%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting

## sum ## [1] 7962

#ihwResult object with 10839 hypothesis tests
#Nominal FDR control level: 0.1
#Split into 7 bins, based on an ordinal covariate

write.csv(outputICihw, file="cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.ihw.csv")

###
## Independent hypothesis weighting (IHW) + Shrinkage of effect size (LFC estimates)
# Replace log2FoldChange, lfcSE, and p-values in outputIClfc with those from cell_model
outputICihwlfc <- outputICihw
outputICihwlfc$log2FoldChange <- outputIClfc$log2FoldChange
outputICihwlfc$lfcSE <- outputIClfc$lfcSE
outputICihwlfc$pvalue <- outputIClfc$pvalue
outputICihwlfc$padj <- outputIClfc$padj

write.csv(outputICihwlfc, file="cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.ihw.lfc.csv")

###
## In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
## Points will be colored red(?) if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

pdf("cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel_MAplot.ihw.pdf")
plotMA(outputICihw, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel_MAplot.pdf")
plotMA(outputIC, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel_MAplot.lfc.pdf")
plotMA(outputIClfc, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel_MAplot.ihw.lfc.pdf")
plotMA(outputICihwlfc, ylim=c(-2,2))
dev.off()


################################################################################
################################################################################
################################################################################
#####################################################
######################################################
######## Volcano Plots from DESeq2 #######
######################################################
#### returning to default p-adjust of <1e-5 (not 1e-3)
DEseq_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/2-normalization/DESeq2/v1_SLR"
wald_path <- file.path(DEseq_dir)

csv <- "cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

## Warning message:
## One or more p-values is 0. Converting to 10^-1 * current lowest non-zero p-value...
###################
#######

csv <- "cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.lfc.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()


##

csv <- "cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.ihw.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()


##

csv <- "cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.ihw.lfc.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##################
## cell type
##############
csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.lfc.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.ihw.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

#########################
## interaction infection x cell type
##########################

csv <- "cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.lfc.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.ihw.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.ihw.lfc.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
library(ggVennDiagram)

x <- list(infection=row.names(outputI[which(outputI$padj<=0.01), ]), cell_type=row.names(outputC[which(outputC$padj<=0.01), ]), interaction=row.names(outputIC[which(outputIC$padj<=0.01), ]))

pdf("cells_Kallisto-DESeq2_WaldTest-Infection-vs-Celltype-vs-Interaction-WaldTest-VennDiagram.pdf")
ggVennDiagram(x)
dev.off()

# Create a Venn object explicitly
venn_obj <- Venn(x)

# Process region data to extract categories
regions <- process_region_data(venn_obj)

# View the regions
print(regions)

> print(regions)
# A tibble: 7 × 4
  id    name                            item          count
  <chr> <chr>                           <list>        <int>
1 1     infection                       <chr [444]>     444
2 2     cell_type                       <chr [1,563]>  1563
3 3     interaction                     <chr [396]>     396
4 1/2   infection/cell_type             <chr [1,273]>  1273
5 1/3   infection/interaction           <chr [750]>     750
6 2/3   cell_type/interaction           <chr [1,350]>  1350
7 1/2/3 infection/cell_type/interaction <chr [4,071]>  4071


region_data <- data.frame(
  Region = names(regions), # Region names (e.g., A, A&B, A&B&C)
  Elements = sapply(regions, function(x) paste(x, collapse = ", ")) # Elements in each region
)

# View the data
print(region_data)

# Save the data to a CSV file
write.csv(region_data, "venn_diagram_regions.p01.csv", row.names = FALSE)

# Output path
cat("Saved Venn diagram region data to 'venn_diagram_regions.p01.csv'\n")

################################################################################
################################################################################

## Note: use subset_DE_venn_regions.pl to create subsetted versions of
#cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.csv
#cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv
#cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.csv
## containing only the unique hits for each set (interaction, celltype, infection)
## input these subsetted csvs Volcano

## volcano plots of unique Wald Test hits

##
csv <- "subsetVenn_cells_Kallisto-DESeq2_WaldTest-wMel_vs_uninf.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "subsetVenn_cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

###

csv <- "subsetVenn_cells_Kallisto-DESeq2_WaldTest-celltypeJW18.infectionwMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 10e-32, FCcutoff = 2, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()


################################################################################
################################################################################
################################################################################
###################
##### genes by id
##################
gene <- "larval_protein"
ID <- "Dmel_CG3984"
###
gene <- "cilium_assembly"
ID <- "Dmel_CG4970"
####
gene <- "lncRNACR45304"
ID <- "Dmel_CR45304"
####
gene <- "unknown_protein"
ID <- "Dmel_CG43185"
####
gene <- "cookie_monster"
ID <- "Dmel_CG13493"
####
gene <- "unknown_protein"
ID <- "Dmel_CG46460"
####
gene <- "His2B"
ID <- "Dmel_CG17949"
####
gene <- "unknown_protein"
ID <- "Dmel_CG13622"
####
gene <- "asRNACR44340"
ID <- "Dmel_CR44340"
####
gene <- "unknown_protein"
ID <- "Dmel_CG34323"
####
gene <- "Ance3"
ID <- "Dmel_CG17988"
####
S2ids <- c("Dmel_CG8776", "Dmel_CG5958", "Dmel_CG4099", "Dmel_CG6953", "Dmel_CG1358", "Dmel_CG14715", "Dmel_CG1803", "Dmel_CG9372")
JW18ids <- c("Dmel_CG10580", "Dmel_CG4905", "Dmel_CG2759", "Dmel_CG4451", "Dmel_CG8250")
##
## crystal cells and plasmatocyte
##
gene <- "nemy"
ID <- "Dmel_CG8776"
####
gene <- "CG5958"
ID <- "Dmel_CG5958"
####
gene <- "Sr-Cl"
ID <- "Dmel_CG4099"
####
gene <- "fat-spondin"
ID <- "Dmel_CG6953"
####
##
## glycoproteins ###
##
gene <- "glt"
ID <- "Dmel_CG9280"
####
gene <- ""
ID <- ""
##
gene <- "fringe"
ID <- "Dmel_CG10580"
##
gene <- "regucalcin"
ID <- "Dmel_CG1803"
####


## fat body
##
####
gene <- "HisT"
ID <- "Dmel_CG1358"
####
gene <- "FKBP2"
ID <- "Dmel_CG14715"
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- "Syn2"
ID <- "Dmel_CG4905"
####
gene <- "SP34"
ID <- "Dmel_CG9372"
####
gene <- "w"
ID <- "Dmel_CG2759"
####
gene <- "Hs6st"
ID <- "Dmel_CG4451"
####
gene <- "Alk"
ID <- "Dmel_CG8250"
####
gene <- "MalA5"
ID <- "Dmel_CG30359"
####
gene <- "stj"
ID <- "Dmel_CG12295"
####
gene <- "argos"
ID <- "Dmel_CG4531"
####
gene <- "Tig"
ID <- "Dmel_CG11527"
####
gene <- "Idgf3"
ID <- "Dmel_CG4559"
####
gene <- "Idgf2"
ID <- "Dmel_CG4475"
####
gene <- "Idgf1"
ID <- "Dmel_CG4472"
####
gene <- "Mur2B"
ID <- "Dmel_CG14796"
####
gene <- "stumps"
ID <- "Dmel_CG31317"
####
gene <- "lncRNA:roX1"
ID <- "Dmel_CR32777"
####
gene <- "mira"
ID <- "Dmel_CG12249"
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
##################
## Infection ######
## padj0 lfc5
####################
gene <- "Ance"
ID <- "CG8827"
####
gene <- "superdeath"
ID <- "CG14516"
####

gene <- ""
ID <- "CG46339"
####
gene <- "yin"
ID <- ""
####
gene <- "Prip"
ID <- ""

###
## infection, unique
####
gene <- "plasmatocyte_protein"
ID <- "CG5397"
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####
gene <- ""
ID <- ""
####



####
pdf1 <- file.path(paste0("cell_KallistoDESeq2_", ID, "_", gene, "_counts.pdf"))
pdflog2 <- file.path(paste0("cell_KallistoDESeq2_", ID, "_", gene, "_log2.counts.pdf"))

obj <- plotCounts(dds, gene=ID, intgroup="group", returnData=TRUE)
pdf(pdf1)
ggplot(data=obj, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf(pdflog2)
ggplot(data=obj, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()



#######################################################################################
#######################################################################################
### Manual Plots

sordd2 <- plotCounts(dds, gene="Dmel_CG32581", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG32581_sordd2_counts.pdf")
ggplot(data=sordd2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG32581_sordd2_log2.counts.pdf")
ggplot(data=sordd2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()

##
#################################################################################
## lowest DE cell type (S2)     Dmel_CG8965
rau <- plotCounts(dds, gene="Dmel_CG8965", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG8965_rau_counts.pdf")
ggplot(data=rau, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG8965_rau_log2.counts.pdf")
ggplot(data=rau, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


## lowest DE cell type (JW18)     Dmel_CG1505
gd <- plotCounts(dds, gene="Dmel_CG1505", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG1505_gd_counts.pdf")
ggplot(data=gd, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG1505_gd_log2.counts.pdf")
ggplot(data=gd, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


mAChRA <- plotCounts(dds, gene="Dmel_CG4356", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG4356_mAChR-A_counts.pdf")
ggplot(data=mAChRA, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG4356_mAChR-A_log2.counts.pdf")
ggplot(data=mAChRA, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


asRNACR43767 <- plotCounts(dds, gene="Dmel_CR43767", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CR43767_asRNACR43767_counts.pdf")
ggplot(data=asRNACR43767, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CR43767_asRNACR43767_log2.counts.pdf")
ggplot(data=asRNACR43767, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


chitin_binding_protein <- plotCounts(dds, gene="Dmel_CG6933", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG6933_chitin_binding_protein_counts.pdf")
ggplot(data=chitin_binding_protein, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG6933_chitin_binding_protein_log2.counts.pdf")
ggplot(data=chitin_binding_protein, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()



#############
################################################################################

## lowest padj, lowest DE infection (wMel)
mew <- plotCounts(dds, gene="Dmel_CG1771", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG1771_mew_counts.pdf")
ggplot(data=mew, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG1771_mew_log2.counts.pdf")
ggplot(data=mew, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


## lowest padj, highest DE infection (wMel)
bowl <- plotCounts(dds, gene="Dmel_CG10021", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG10021_bowl_counts.pdf")
ggplot(data=bowl, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG10021_bowl_log2.counts.pdf")
ggplot(data=bowl, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()

## lowest padj, highest DE infection (wMel)
rt <- plotCounts(dds, gene="Dmel_CG6097", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG6097_rt_counts.pdf")
ggplot(data=rt, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG6097_rt_log2.counts.pdf")
ggplot(data=rt, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()

muc11A <- plotCounts(dds, gene="Dmel_CG32656", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG32656_muc11A_counts.pdf")
ggplot(data=muc11A, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG32656_muc11A_log2.counts.pdf")
ggplot(data=muc11A, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


##
################################################################################

## lowest DE infection x cell type
yin <- plotCounts(dds, gene="Dmel_CG44402", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG44402_yin_counts.pdf")
ggplot(data=yin, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG44402_yin_log2.counts.pdf")
ggplot(data=yin, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


## highest DE infection x cell type
cac <- plotCounts(dds, gene="Dmel_CG43368", intgroup="group", returnData=TRUE)
pdf("cell_KallistoDESeq2_Dmel_CG43368_cac_counts.pdf")
ggplot(data=cac, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf("cell_KallistoDESeq2_Dmel_CG43368_cac_log2.counts.pdf")
ggplot(data=cac, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
