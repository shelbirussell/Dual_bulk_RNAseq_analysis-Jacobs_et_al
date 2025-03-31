## Purpose: output WaldTest csvs

## Prereqs: load dds DEseq2 results 

##
#######################
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
