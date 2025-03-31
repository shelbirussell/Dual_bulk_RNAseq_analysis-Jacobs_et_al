# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)  # Use the appropriate organism package for gene ID mapping
library(DOSE)          # For disease ontology semantic and enrichment analysis
library(ggridges)
library(enrichplot)
library(ggupset)

## GSEA
## original GSEA paper, PNAS 2005 www.pnas.orgcgidoi10.1073pnas.0506580102
## ranked based on the correlation between their expression and the class distinction by using any suitable metric
## ==> microarrays used signal/noise
## ==> DE uses -log10(p) x sign(log2FC) to normalize pvalues by expression levels

## ==> upload DEseq csvs and rank by adjusted pvalue


## https://www.biostars.org/p/375584/
## https://crazyhottommy.blogspot.com/2016/08/gene-set-enrichment-analysis-gsea.html
## the first thing that GSEA does is to rank the genes in L based on "how well they divide the conditions"
## now you want to see whether the genes present in a gene set (S) are at the top or at the bottom of your list...or if they are just spread around randomly.
## to do that GSEA calculates the famous enrichment score, that becomes normalized enrichment score (NES) when correcting for multiple testing (FDR).
##  GSEA is a tool that uses **every** datapoint in its statistical algorithm.

## In the "classical" method genes are ranked by from most up-regulated to most down-regulated
## The rank metric itself varies, but two valid methods are to use signed p-value, or lower 90% confidence interval of the fold change.
#### assess whether members of a gene set appear enriched at one end of the profile
######## GSEA performs permutations of the profile, calculating the enrichment of the gene set a thousand or more times to estimate p-values empirically.
## The default metric for ranking genes is the signal-to-noise ratio.
## To use this metric, your phenotype file must define at least two categorical phenotypes [infected and cell type]
## and your expression dataset must contain at least three (3)samples for each phenotype. [6 each]

## internally GSEA is going to rank your list of genes first and then it will compare the public-curated gene sets with the ranked gene list you have.
## you need to supply ALL genes that are detected in the experiment

## https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=rank#running-score-and-preranked-list-of-gsea-result

## HOW TO RANK????
## The "stat" column from a DESeq2 analysis is a great ranking metric for Gene Set Enrichment Analysis (GSEA). It represents the Wald statistic, which is calculated as:
## Combines both effect size (log2FC) and significance: Unlike raw log2 fold change, which only measures magnitude, the stat column accounts for variability by incorporating the standard error.
## Better than using adjusted p-value alone: P-values can be misleading for ranking because small expression differences in highly expressed genes can produce extremely low p-values.
## Smooths over noise: Genes with high variability (large standard error) will get lower ranking even if they have high log2FC.

######################################################
######## GO/KEGG Plots from DESeq2 #######
######################################################
DEseq_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/2-normalization/DESeq2/v1_SLR/DE_csvs"
wald_path <- file.path(DEseq_dir)

###############################################################################################
############################## DE ~ cell type ###############################################
###############################################################################################
csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.csv"

# Load the DESeq2 results file
csv_path <-file.path(wald_path, csv)
deseq_results <- read.csv(csv_path)

# Check the column names to identify log2 fold change and p-value columns
head(deseq_results)

# Ensure the data has columns "log2FoldChange", "pvalue", and "padj"
# Filter to keep genes with valid log2 fold change and p-value
deseq_results <- deseq_results[!is.na(deseq_results$log2FoldChange) & !is.na(deseq_results$pvalue) & !is.na(deseq_results$padj), ]

## export top hit FBGNids to file for PANGEA GO analysis
expressed <- sub('^Dmel_', '', (deseq_results$X [which(deseq_results$padj<=1e-64)]))
expressed_ids<-bitr(expressed, fromType = "FLYBASECG", toType = "FLYBASE", OrgDb="org.Dm.eg.db")
write.table(expressed_ids, "celltype.DE_p1e-64_FBGNids.txt")


###############################################################################
########### KEGG uses Dmel_CG gene names ####################
# Rank genes by stat (high to low order)
gene_list <- deseq_results$stat
names(gene_list) <- deseq_results$X  # Ensure this column contains processed gene names
gene_list <- sort(gene_list, decreasing = TRUE)



########################### GSEA using KEGG pathways #####################################################
################## KEGG ##################################################################################
######### ALL pathways ###################################################################################
##########################################################################################################
## https://rdrr.io/bioc/clusterProfiler/man/gseKEGG.html
## [default] Expected input gene ID: Dmel_CG17654,Dmel_CG17324,Dmel_CG17333,Dmel_CG32072,Dmel_CG12229,Dmel_CG4001
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "dme",  # Replace with the appropriate organism code
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.5)


# View KEGG GSEA results
print(head(gsea_kegg))
length(gsea_kegg$ID)

# Save the results to CSV
write.csv(as.data.frame(gsea_kegg), "celltype_WaldTest-DE_GSEA_gseKEGG_results.p0.5.csv")

### plot results
# Check if gsea_kegg contains any results before plotting
if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
        pdfGseKEGG <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.dotplot.p0.5.pdf"))
        pdf(pdfGseKEGG)
        print(dotplot(gsea_kegg, showCategory = 20))
        dev.off()


        for (i in seq_along(gsea_kegg$ID)) {
                pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset", i, ".gseKEGG.plot2.p0.5.pdf"))
                pdf(pdfGse2KEGG)
                print(gseaplot2(gsea_kegg, geneSetID = i, title=gsea_kegg$Description[i]))
                dev.off()
        }


        ## relaxed p-value to get more output - only lysozyme is significant **
        ## https://www.biostars.org/p/433131/

        pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset.1-", length(gsea_kegg$ID), ".gseKEGG.plot2.p0.5.pdf"))
        pdf(pdfGse2KEGG)
        print(gseaplot2(gsea_kegg, geneSetID = 1:length(gsea_kegg$ID), title=gsea_kegg$Description[i], pvalue_table = TRUE))
        dev.off()

        ## https://rdrr.io/bioc/enrichplot/man/ridgeplot.html
        pdfGseKEGGridge <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.ridgeplot.p0.05.pdf"))
        pdf(pdfGseKEGGridge)
        print(ridgeplot(gsea_kegg))
        dev.off()

        gsea_kegg2 <- pairwise_termsim(gsea_kegg)
        pdfGseKEGGemap <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.emapplot.p0.05.pdf"))
        pdf(pdfGseKEGGemap)
        print(emapplot(gsea_kegg2, showCategory=20))
        dev.off()

}

##########################################################################################################
################## KEGG ##################################################################################
######### Significant pathways ###########################################################################
##########################################################################################################

## GSEA using KEGG pathways
## https://rdrr.io/bioc/clusterProfiler/man/gseKEGG.html
## [default] Expected input gene ID: Dmel_CG17654,Dmel_CG17324,Dmel_CG17333,Dmel_CG32072,Dmel_CG12229,Dmel_CG4001
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "dme",  # Replace with the appropriate organism code
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05)


# View KEGG GSEA results
print(head(gsea_kegg))
length(gsea_kegg$ID)

# Save the results to CSV
write.csv(as.data.frame(gsea_kegg), "celltype_WaldTest-DE_GSEA_gseKEGG_results.csv")

### plot results
# Check if gsea_kegg contains any results before plotting
if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 1) {
        pdfGseKEGG <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.dotplot.pdf"))
        pdf(pdfGseKEGG)
        print(dotplot(gsea_kegg, showCategory = 20))
        dev.off()

        for (i in seq_along(gsea_kegg$ID)) {
                pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset", i, ".gseKEGG.plot2.pdf"))
                pdf(pdfGse2KEGG)
                print(gseaplot2(gsea_kegg, geneSetID = i, title=gsea_kegg$Description[i]))
                dev.off()

        }

        ## https://rdrr.io/bioc/enrichplot/man/ridgeplot.html
        pdfGseKEGGridge <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.ridgeplot.pdf"))
        pdf(pdfGseKEGGridge)
        print(ridgeplot(gsea_kegg))
        dev.off()

        gsea_kegg2 <- pairwise_termsim(gsea_kegg)
        pdfGseKEGGemap <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.emapplot.pdf"))
        pdf(pdfGseKEGGemap)
        print(emapplot(gsea_kegg2, showCategory=20))
        dev.off()

        pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset.1-", length(gsea_kegg$ID), ".gseKEGG.plot2.pdf"))
        pdf(pdfGse2KEGG)
        print(gseaplot2(gsea_kegg, geneSetID = 1:length(gsea_kegg$ID), title=gsea_kegg$Description[i], pvalue_table = TRUE))
        dev.off()

}


if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) == 1) {
        pdfGseKEGG <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.dotplot.pdf"))
        pdf(pdfGseKEGG)
        print(dotplot(gsea_kegg, showCategory = 20))
        dev.off()

        for (i in seq_along(gsea_kegg$ID)) {
                pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset", i, ".gseKEGG.plot2.pdf"))
                pdf(pdfGse2KEGG)
                print(gseaplot2(gsea_kegg, geneSetID = i, title=gsea_kegg$Description[i]))
                dev.off()
        }

        ## https://rdrr.io/bioc/enrichplot/man/ridgeplot.html
        pdfGseKEGGridge <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.ridgeplot.pdf"))
        pdf(pdfGseKEGGridge)
        print(ridgeplot(gsea_kegg))
        dev.off()

        gsea_kegg2 <- pairwise_termsim(gsea_kegg)
        pdfGseKEGGemap <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.emapplot.pdf"))
        pdf(pdfGseKEGGemap)
        print(emapplot(gsea_kegg2, showCategory=20))
        dev.off()

}



################################################################################
# GSEA using GO terms

### Remove the "Dmel_" prefix from gene names - GSEA GO uses ENSEMBL names #######
deseq_results$gene <- gsub("^Dmel_", "", deseq_results$X)  # Replace 'gene' with the actual column name if different

g <- deseq_results$gene

## Map FLYBASECG to ENTREZID for gseGO()
## https://rdrr.io/bioc/clusterProfiler/man/gseGO.html
## convert FLYBASECGs to keyType: one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
## We will lose some genes here because not all IDs will be converted
ids<-bitr(g, fromType = "FLYBASECG", toType = "ENTREZID", OrgDb="org.Dm.eg.db")

# Ensure all FLYBASECG have corresponding ENTREZID
if (nrow(ids) == 0) {
        message("No valid ENTREZ IDs found for module: ", module)
        next
}

# Filter out NA or missing mappings
ids <- ids[!is.na(ids$ENTREZID), ]

# Match FLYBASECG IDs to matrix column names and replace with ENTREZID
deseq_results$ENTREZID <- ids$ENTREZID[match(deseq_results$gene, ids$FLYBASECG)]

gene_list <- deseq_results$stat
names(gene_list) <- deseq_results$ENTREZID  # Use ENTREZ IDs as names
gene_list <- gene_list[!is.na(names(gene_list))]  # Remove NAs from the gene list
gene_list <- sort(gene_list, decreasing = TRUE)


gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = 'org.Dm.eg.db',  # Replace with the correct annotation package
                 ont = "ALL",  # Biological Process; can also use "MF" or "CC"
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05)

# View GO GSEA results
head(gsea_go)
length(gsea_go$ID)


# Save the results to CSV
write.csv(as.data.frame(gsea_go), "celltype_WaldTest-DE_GSEA_gseGO_results.csv")

## Plot results
# Check if gsea_go contains any results before plotting
if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {

        pdfGseGO <- file.path(paste0("celltype_WaldTest-DE_", ".gseGO.dotplot.pdf"))
        pdf(pdfGseGO)
        print(dotplot(gsea_go, showCategory = 20))
        dev.off()

        pdfGseGOridge <- file.path(paste0("celltype_WaldTest-DE_", ".gseGO.ridgeplot.pdf"))
        pdf(pdfGseGOridge)
        print(ridgeplot(gsea_go))
        dev.off()

        pdfGseGOemap <- file.path(paste0("celltype_WaldTest-DE_", ".gseGO.emapplot.pdf"))
        pdf(pdfGseGOemap)
        gsea_go2 <- pairwise_termsim(gsea_go)
        print(emapplot(gsea_go2))
        dev.off()

        ## plot all genes in category
        for (i in seq_along(gsea_go$ID)) {
                pdfGse2GO <- file.path(paste0("celltype_WaldTest-DE", ".geneset", i, ".gseGO.plot2.pdf"))
                pdf(pdfGse2GO)
                print(gseaplot2(gsea_go, geneSetID = i, title=gsea_go$Description[i]))
                dev.off()

                pdfGseGOcowplot <- file.path(paste0("celltype_WaldTest-DE", ".geneset", i, ".gseGO.cowplot.pdf"))
                pdf(pdfGseGOcowplot)
                p1 <- gseaplot(gsea_go, geneSetID = i, by = "runningScore", title = gsea_go$Description[i])
                p2 <- gseaplot(gsea_go, geneSetID = i, by = "preranked")
                print(cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2]))
                dev.off()
        }


        pdfGse2GO <- file.path(paste0("celltype_WaldTest-DE", ".geneset.1-", length(gsea_go$ID), ".gseGO.plot2.pdf"))
        pdf(pdfGse2GO)
        print(gseaplot2(gsea_go, geneSetID = 1:length(gsea_go$ID), title=gsea_go$Description[i]))
        dev.off()


        gsea_go2 <- pairwise_termsim(gsea_go)


        # Enhanced checks for gsea_go2
        if (is.null(gsea_go2)) {
            message(paste("pairwise_termsim() returned NULL for module", module, "and term", term))
            next
        }

        # Check if the object contains a valid similarity matrix
        if (is.null(gsea_go2@termsim) || nrow(gsea_go2@termsim) < 2 || ncol(gsea_go2@termsim) < 2) {
            message(paste("Insufficient or missing pairwise similarity matrix in gsea_go2 for module", module, "and term", term))
            next
        }

        # Check if the similarity matrix has any meaningful values
        if (all(is.na(gsea_go2@termsim)) || sum(gsea_go2@termsim, na.rm = TRUE) == 0) {
            message(paste("pairwise_termsim matrix is empty or invalid for module", module, "and term", term))
            next
        }

        # Check if there are enough terms to cluster
        if (nrow(gsea_go2@result) < 2) {
            message(paste("Not enough enriched terms to perform clustering for module", module, "and term", term))
            next
        }

        # Generate plots if all checks pass
        tryCatch({
            pdfGOtreeplot <- file.path(paste0("celltype_WaldTest-DE", ".gseGO.treeplot.pdf"))
            pdf(pdfGOtreeplot)
            print(treeplot(gsea_go2, cluster.params = list(method = "average")))
            dev.off()

            pdfGseGOemap <- file.path(paste0("celltype_WaldTest-DE", ".gseGO.emapplot.pdf"))
            pdf(pdfGseGOemap)
            print(emapplot(gsea_go2, showCategory=20))
            dev.off()

            pdfGseGOupset <- file.path(paste0("celltype_WaldTest-DE", ".gseGO.upset.pdf"))
            pdf(pdfGseGOupset)
            print(upsetplot(gsea_go))
            dev.off()

        }, error = function(e) {
            message(paste("Error in treeplot for module", module, "and term", term, ":", e$message))
        })
}


## this works, but is included in the gseaplot(by="preranked") function
##gsearank(gsea_go, geneSetID = 1)
