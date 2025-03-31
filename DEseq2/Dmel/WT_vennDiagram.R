## Purpose: calculate overlap of Wald Test gene sets and create unique lists

## Prereqs: load dds DEseq2 results 

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
