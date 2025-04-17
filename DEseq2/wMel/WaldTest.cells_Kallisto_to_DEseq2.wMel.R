
############################################################################
#### set "controls" - JW18 vs S2 (control) and wMel vs uninfected (control)
############################################################################
## run relevel before running DESeq() - can also use "contrast" later

dds$celltype <- relevel(dds$celltype, ref = "S2")

################################################################################
################################################################################
################################################################################
##################################################################
########################### Wald test ############################
##################################################################

#################################################################
########## WALD formatting: retain the terms to test ############
#################################################################
design(dds) <- ~ celltype
dds <- DESeq(dds)
resultsNames(dds)

#[1] "Intercept"           "celltype_JW18_vs_S2"

#################################################################
## cell type effect     celltype_JW18_vs_S2

outputC <- results(dds, name="celltype_JW18_vs_S2")
write.csv(outputC, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.wMel.csv")

outputClfc <- lfcShrink(dds, coef="celltype_JW18_vs_S2", type="apeglm")
write.csv(outputClfc, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.lfc.wMel.csv")

summary(outputC)

#out of 1122 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 255, 23%
#LFC < 0 (down)     : 286, 25%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 9)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

summary(outputClfc)

#out of 1122 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 255, 23%
#LFC < 0 (down)     : 286, 25%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 9)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


###############################

## Independent hypothesis weighting = IHW
outputCihw <- results(dds, name="celltype_JW18_vs_S2", filterFun=ihw)
summary(outputCihw)
sum(outputCihw$padj < 0.1, na.rm=TRUE)
metadata(outputCihw)$ihwResult

#out of 1122 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 255, 23%
#LFC < 0 (down)     : 286, 25%
#outliers [1]       : 0, 0%
#[1] see 'cooksCutoff' argument of ?results
#see metadata(res)$ihwResult on hypothesis weighting

#> sum(outputCihw$padj < 0.1, na.rm=TRUE)
#[1] 541

#> metadata(outputCihw)$ihwResult
#ihwResult object with 1122 hypothesis tests
#Nominal FDR control level: 0.1
#Split into 1 bins, based on an ordinal covariate


write.csv(outputCihw, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.ihw.wMel.csv")

###
## Independent hypothesis weighting (IHW) + Shrinkage of effect size (LFC estimates)
# Replace log2FoldChange, lfcSE, and p-values in outputClfc with those from cell_model
outputCihwlfc <- outputCihw
outputCihwlfc$log2FoldChange <- outputClfc$log2FoldChange
outputCihwlfc$lfcSE <- outputClfc$lfcSE
outputCihwlfc$pvalue <- outputClfc$pvalue
outputCihwlfc$padj <- outputClfc$padj

## these can't be merged until I can reorder the levels upstream
write.csv(outputCihwlfc, file="cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.ihw.lfc.wMel.csv")

###
## In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
## Points will be colored red(?) if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

pdf("cells_Kallisto-DESeq2_WaldTest-JW18.S2.celltype_MAplot.ihw.wMel.pdf")
plotMA(outputCihw, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-JW18.S2.celltype_MAplot.wMel.pdf")
plotMA(outputC, ylim=c(-2,2))
dev.off()

pdf("cells_Kallisto-DESeq2_WaldTest-JW18.S2.celltype_MAplot.lfc.wMel.pdf")
plotMA(outputClfc, ylim=c(-2,2))
dev.off()



################################################################################
################################################################################
################################################################################
#####################################################
######################################################
######## Volcano Plots from DESeq2 #######
######################################################
#### returning to default p-adjust of <1e-5 (not 1e-3)
DEseq_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/2-normalization/DESeq2/wMel_SLR"
wald_path <- file.path(DEseq_dir)


##################
## cell type
##############
csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.wMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.wMel.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.wMel.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.wMel.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-3, FCcutoff = 0.5, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-3, FCcutoff = 0.5, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.lfc.wMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.wMel.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.wMel.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.wMel.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-3, FCcutoff = 0.5, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-3, FCcutoff = 0.5, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
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
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-3, FCcutoff = 0.5, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-3, FCcutoff = 0.5, widthConnectors = 0.75, col = c("grey30", "cyan3", "skyblue", "orangered1"), shape=20)
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

##

################################################################################
################################################################################
################################################################################
###################
##### genes by id
##################
gene <- "capsid"
ID <- "WD_RS02030"
####
gene <- "AAA_ATPase"
ID <- "WD_RS02730"
####
gene <- "subtilisin"
ID <- "WD_RS02575"
####
gene <- "AAA_ATPase"
ID <- "WD_RS02645"
####
gene <- "riboflabin_synthase"
ID <- "WD_RS01055"
####
gene <- "TrbC.VirB2"
ID <- "WD_RS02115"
####
gene <- "hypothetical"
ID <- "WD_RS04395"
####
gene <- "septal_ring_transglycosylase_RlpA"
ID <- "WD_RS02210"
####
gene <- "flavin_oxidoreductase"
ID <- "WD_RS05305"
####
gene <- "DUF4815_domain_protein"
ID <- "WD_RS02600"
####
gene <- "ankyrin_repeat_domain_protein"
ID <- "WD_RS01940"
####
gene <- "FtsK.SpoIIIE_family_DNA_translocase"
ID <- "WD_RS00540"
####
gene <- "ankyrin_repeat_domain_protein"
ID <- "WD_RS03475"
####
gene <- "hypothetical_protein"
ID <- "WD_RS03775"
####
gene <- "tRNA.Trp"
ID <- "WD_RS00090"
####
gene <- "tRNA.Pro"
ID <- "WD_RS05495"
####
gene <- "tRNA.Lys"
ID <- "WD_RS04410"
####
gene <- "tRNA.Leu"
ID <- "WD_RS01745"
####
gene <- "tRNA.Tyr"
ID <- "WD_RS03110"
####
gene <- "dihydroneopterin_aldolase"
ID <- "WD_RS03980"
####
gene <- "tRNA.Gly"
ID <- "WD_RS03105"
####
gene <- "elongation_factor_Tu"
ID <- "WD_RS03100"
####
gene <- "metK"
ID <- "WD_RS00610"
####
####
gene <- "tRNA.Asp"
ID <- "WD_RS02375"
####
gene <- "hypothetical"
ID <- "WD_RS00535"
####
gene <- "tRNA.Ala"
ID <- "WD_RS04405"
####
gene <- "tRNA.Val"
ID <- "WD_RS05450"
####
gene <- "hypothetical"
ID <- "WD_RS06605"
####
gene <- "peroxiredoxin"
ID <- "WD_RS03425"
####
####
gene <- "50S_ribosomal_L23"
ID <- "WD_RS03080"
####
gene <- "O-methyltransferase"
ID <- "WD_RS04790"
####
gene <- "transcription_termination_Rho"
ID <- "WD_RS03600"
####
gene <- "hypothetical"
ID <- "WD_RS04310"
####
gene <- "chaperone_groEL"
ID <- "WD_RS01360"
####
gene <- "thioredoxin-disulfide_reductase"
ID <- "WD_RS03420"
####
gene <- "6,7.dimethyl.8.ribityllumazine_synthase"
ID <- "WD_RS01055"
####
gene <- "50S_ribosomal_L3"
ID <- "WD_RS03090"
####
gene <- "Grx4_monothiol_glutaredoxin"
ID <- "WD_RS02185"
####
gene <- "PQQ.binding.like_beta.propeller_protein"
ID <- "WD_RS03400"
####
gene <- "polyribonucleotide_nucleotidyltransferase"
ID <- "WD_RS04085"
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
gene <- "HtpG"
ID <- "WD_RS05780"
####


####
pdf1 <- file.path(paste0("cell_KallistoDESeq2_", ID, "_", gene, "_counts.pdf"))
pdflog2 <- file.path(paste0("cell_KallistoDESeq2_", ID, "_", gene, "_log2.counts.pdf"))

obj <- plotCounts(dds, gene=ID, intgroup="celltype", returnData=TRUE)
pdf(pdf1)
ggplot(data=obj, aes(x=celltype, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()
pdf(pdflog2)
ggplot(data=obj, aes(x=celltype, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
dev.off()


#################
