library(tximport)
library(DESeq2)
library(ggplot2)
library(rhdf5)
library(dplyr)
library(pheatmap)
library(WGCNA)
library(genefilter)
library(ggforce)
library(ComplexHeatmap)
library(tidyr)
library(readr)
library(flashClust)
library(AnnotationDbi)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(enrichplot)
library(ggforce)
library(ggplot2)

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
## Transcripts need to be associated with gene IDs for gene-level summarization - geneIDs need to be associated with transcript IDs here
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
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = t2g, dropInfReps = TRUE, ignoreAfterBar = TRUE)

################################################################################
## use the gene-level estimated counts from the quantification tools, and additionally to use the transcript-level abundance estimates
## to calculate a gene-level offset that corrects for changes to the average transcript length across samples
## the function DESeqDataSetFromTximport takes care of creation of the offset for you

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = metadata, design = ~celltype + infection + celltype:infection)


## prefilter low-count/coverage genes across samples #==> higher cutoff for WGCNA (50 vs 10 for DE)
## keep rows that have at least 10 reads total --> n=6 (all samples)
## keep only genes with 10 or more reads in total across the samples
## normalized = TRUE ==> normalized by sample depth and compositional bias

smallestGroupSize <- 6
minreads <- 70

################################################################################
## normalize counts: divide each column of counts(dds) by sizeFactors(dds)
################################################################################
## set the normalized=TRUE flag because we are not running DESeq2() --> transformed data needs to be normalized by gene length
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE) >= minreads) >= smallestGroupSize
dds <- dds[keep,]

dds$group <- factor(paste0(dds$celltype, dds$infection))

############################################################################
#### set "controls" - JW18 vs S2 (control) and wMel vs uninfected (control)
############################################################################
## run relevel before running DESeq() - can also use "contrast" later

dds$infection <- relevel(dds$infection, ref = "uninfected")
dds$celltype <- relevel(dds$celltype, ref = "S2")


#################################################################
################################################################################
## transformation code adapted from : https://rpubs.com/mahima_bose/921350
## Use the variance-stabilizing transformation (VST) or rlog for downstream correlation analyses:
##  the rlog() function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size ... more robust in the case when the size factors vary widely
## The argument blind = FALSE ensures the transformation accounts for experimental design, making it suitable for visualizations or downstream analyses like WGCNA.
## blind = FALSE because experimental conditions affect the data
rld <- rlog(dds, blind=FALSE)
normalized_rld_counts <- assay(rld)

## Gene Filtering for WGCNA: Consider additional variance-based filtering (e.g., selecting the top 5000 most variable genes) after transformation to focus on biologically relevant signals for co-expression analysis.
# Ensure your normalized_counts is a matrix
wgcna_rld_data <- as.matrix(normalized_rld_counts)

# WGCNA expects samples as rows and genes as columns.
wgcna_rld_data <- t(wgcna_rld_data)

# Remove Missing or Invalid Values:
# WGCNA does not handle missing or infinite values. Clean your data beforehand:
wgcna_rld_data <- wgcna_rld_data[complete.cases(wgcna_rld_data), ]

# Filter Genes Based on Variance:
# WGCNA works best with the most variable genes to focus on biologically relevant patterns. Consider selecting the top variable genes:
## all genes:
top_rld_genes = length(colnames(wgcna_rld_data))

#top_var_rld_genes <- head(order(rowVars(wgcna_rld_data), decreasing = TRUE), top_rld_genes)
#wgcna_rld_data <- wgcna_rld_data[top_var_rld_genes, ]

# Scale the Data (Optional):
# For some datasets, scaling gene expression to zero mean and unit variance improves co-expression detection:
#wgcna_rld_data <- scale(wgcna_rld_data)

length(colnames(wgcna_rld_data))


powers = c(c(1:20), seq(from =1, to=100, by=5)) #choosing a set of soft-thresholding powers
sft_rld <- pickSoftThreshold(wgcna_rld_data, powerVector = powers, verbose = 5, networkType="signed") #call network topology analysis function
rld_chosen_power <- sft_rld$powerEstimate

rld_chosen_power


dir.create("wgcna_output", showWarnings = TRUE, recursive = TRUE)


pdfAutoclusters <- file.path("wgcna_output", paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "genes_soft_threshold_power", rld_chosen_power, "_analysis.ggplot.pdf"))
pdf(pdfAutoclusters)
ggplot(sft_rld$fitIndices, aes(x = Power, y = SFT.R.sq, label = Power)) +
# Plot the points
geom_point() +
# We'll put the Power labels slightly above the data points
geom_text(nudge_y = 0.1) +
# We will plot what WGCNA recommends as an R^2 cutoff
geom_hline(yintercept = 0.80, col = "red") +
# Just in case our values are low, we want to make sure we can still see the 0.80 level
ylim(c(min(sft_rld$fitIndices$SFT.R.sq), 1.05)) +
# We can add more sensible labels for our axis
xlab("Soft Threshold (power)") +
ylab("SFT.R.sq") +
# This adds some nicer aesthetics to our plot
theme_classic()
dev.off()

##############
pdfAutoclusters2 <- file.path("wgcna_output", paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "genes_soft_threshold_power_connectivity", rld_chosen_power, "_analysis.ggplot.pdf"))
pdf(pdfAutoclusters2)
ggplot(sft_rld$fitIndices, aes(x = Power, y = mean.k., label = Power)) +
# Plot the points
geom_point() +
# We'll put the Power labels slightly above the data points
geom_text(nudge_y = 100) +
# We can add more sensible labels for our axis
xlab("Soft Threshold (power)") +
ylab("mean connectivity") +
# This adds some nicer aesthetics to our plot
theme_classic()
dev.off()

pdfAutoclusters2 <- file.path("wgcna_output", paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "genes_soft_threshold_power_log10connectivity", rld_chosen_power, "_analysis.ggplot.pdf"))
pdf(pdfAutoclusters2)
ggplot(sft_rld$fitIndices, aes(x = Power, y = log10(mean.k.), label = Power)) +
# Plot the points
geom_point() +
# We'll put the Power labels slightly above the data points
geom_text(nudge_y = 0.1) +
# We can add more sensible labels for our axis
xlab("Soft Threshold (power)") +
ylab("mean connectivity") +
# This adds some nicer aesthetics to our plot
theme_classic()
dev.off()

#######################

################################################################################
################################################################################
######## AUTOMATE CLUSTERING with blockwiseModules()
################################################################################
################################################################################

net <- blockwiseModules(
  wgcna_rld_data,
  maxBlockSize = 10000, # What size chunks (how many genes) the calculations should be run in
  power = rld_chosen_power,
  networkType = "signed", # Adjacency function options
  TOMType = "signed",  # topological overlap matrix
  saveTOMs = TRUE, # save TOM for downstream steps. The full file names have "block.1.RData", "block.2.RData" etc. appended.
  saveTOMFileBase = paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "_power", rld_chosen_power), # These files are standard R data files and can be loaded using the load function.
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  randomSeed = 1234,
  pamRespectsDendro = FALSE,
  verbose = 3
)

# Visualize the Network:
# Plot module dendrograms and heatmaps to evaluate the module structure.
pdfAutoclusters <- file.path("wgcna_output", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, ".pdf"))

pdf(pdfAutoclusters)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors,
  "Module Colors",  addTextGuide = TRUE, cex.colorLabels = 0.5, cex.dendroLabels = 0.5, cex.rowText = 0.5
)
dev.off()


net$goodGenes
net$goodSamples
net$dendrograms

readr::write_rds(net, file = file.path("wgcna_output", paste0("cells_Dmel_wgcna_", minreads, "reads_", top_rld_genes, "genes_auto.results_power", rld_chosen_power, ".RDS")))


## get data frame of eigengene module data for each sample in the MEs slot.
## These represent the collapsed, combined, and normalized expression of the genes that make up each module.
module_eigengenes <- net$MEs

# Print out a preview
head(module_eigengenes)


## pull out list of modules
modules_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)

write.csv(modules_df, file = file.path("wgcna_output", paste0("gene_modules_", minreads, "reads_", top_rld_genes, "genes_auto.results_power", rld_chosen_power, ".txt")))


## get data frame of eigengene module data for each sample in the MEs slot.
## These represent the collapsed, combined, and normalized expression of the genes that make up each module.
module_eigengenes <- net$MEs


## pull out list of modules
modules_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)

## Calculating and Visualizing Module Eigengenes
## investigate the expression profiles in the modules and relationships among the modules using the eigengenes as representative profiles:

MEs = moduleEigengenes(wgcna_rld_data, net$colors)$eigengenes

###############################
## create a heatmap of adjacencies among eigengenes and a dendrogram for their relationship
pdfEigengeneHeatmap <- file.path("wgcna_output", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.eigengene.heatmap", rld_chosen_power, ".pdf"))
pdf(pdfEigengeneHeatmap)
plotEigengeneNetworks(MEs, names(MEs), marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

## see the size of each module color-coded, run the following command:
table(net$colors)

## export the expression profiles of eigengenes in a text table, run the following command:
write.table(MEs, file="eigengenes.txt", sep="\t")

## export genes and their log2FC in each module as a separate table file in the working directory, run the following command:

for (color in net$colors) {
        module=wgcna_rld_data[which(net$colors==color)]
        write.table(module, file.path("wgcna_output", paste0("module_", color, ".txt",sep="")), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
}


####################################################################################
## To visualize the resulting network module in Cytoscape and perform downstream analyses, the information of nodes, edges, and connection weights need to be extracted.
## The following commands will create both edge and node files:
#table(net$colors) # To see the modules
#modules = c(“lightcyan”); # To select the specific module. Users can choose different modules of interests.

adjacency= adjacency(wgcna_rld_data, type = "signed", power = rld_chosen_power)

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

dir.create("cytoscape_files", showWarnings = TRUE, recursive = TRUE)

for (module in names(table(net$colors))) {
        probes = colnames(wgcna_rld_data) # Select module probes
        inModule = is.finite(match(net$colors, module));
        modProbes = probes[inModule];
        modTOM = TOM[inModule, inModule]; # Select the corresponding Topological Overlap
        dimnames(modTOM) = list(modProbes, modProbes)
        cyt = exportNetworkToCytoscape(modTOM, edgeFile = file.path("cytoscape_files", paste0("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep="")),
         nodeFile = file.path("cytoscape_files", paste0("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep="")),
         weighted = TRUE,
         threshold = 0.2,
         nodeNames = modProbes,
         nodeAttr = net$colors[inModule]);
}

####################################################################################
################################################################################
## DO ANY OF THE MODULES HAVE SIGNIFICANT ASSOCIATIONS WITH TREATMENT GROUPS?
## From: https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
### ==> do the eigengenes relate to metadata labels?

## check that the samples are still in order
all.equal(metadata$sample, rownames(module_eigengenes))

## create file to save stats to
stats_txt <- file.path("wgcna_output", paste0("cells_Kallisto-DESeq2_WGCNA_", minreads, "reads_", top_rld_genes, "genes_auto.stats_power", rld_chosen_power, ".csv"))

######################################
### Run a linear model on each module
#####################################
########## INFECTION ###############
####################################

# Create the design matrix from the `infection` and 'celltype' variables
des_mat <- model.matrix(~ metadata$infection + metadata$celltype)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

#########
# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
head(stats_df)
write.csv(stats_df, file=stats_txt)

module_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::inner_join(
    metadata %>%
      dplyr::select(all_of(c("sample", "infection", "celltype"))),
    by = c("sample" = "sample")
  )

########################
###### PLOT Modules - treatment/type vs eigenmodule values (These represent the collapsed, combined, and normalized expression of the genes that make up each module.)
########################
# a boxplot with outlier points hidden (they will be in the sina plot)
# A sina plot to show all of the individual data points

### plot all modules M# ---> M#

dir.create("output_plots", showWarnings = TRUE, recursive = TRUE)

for (M in stats_df$module) {
    # Debug: Check if the column exists in module_df
    if (!(M %in% names(module_df))) {
        print(paste("Column", M, "not found in module_df"))
        next  # Skip to the next iteration if M is not a valid column
    }

    # Construct the output PDF file path
    pdfout <- file.path("output_plots", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, ".pdf"))

    # Debug: Check if the column exists in the data
    print(paste("Plotting for module:", M))
    print(head(module_df[, c("infection", "celltype", M)]))  # Show first few rows for verification

    # Open a PDF device to save the plot
    pdf(pdfout)

    # Create the plot
    plot <- ggplot(module_df, aes(x = infection, y = !!sym(M), color = celltype)) +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        ggforce::geom_sina(maxwidth = 0.3) +
        theme_classic()

    # Print the plot to ensure it's being rendered
    print(plot)

    # Close the PDF device
    dev.off()

    # Optionally, print the module name for progress tracking
    print(paste("Saved plot for module:", M))
}


###############################################################################
###############################################################################

## What genes are a part of each module?
## If you want to know which of your genes make up a modules, you can look at the $colors slot. This is a named list which associates the genes with the module they are a part of. We can turn this into a data frame for handy use.

# Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
gene_module_key <- tibble::enframe(net$colors, name = "gene", value = "module") %>%  dplyr::mutate(module = paste0("ME", module))
readr::write_tsv(gene_module_key, file = file.path("wgcna_output", paste0("cells_wgcna_", minreads, "reads_", top_rld_genes, "genes_.auto_gene_to_module.tsv")))

dir.create("output_tsvs", showWarnings = TRUE, recursive = TRUE)

for (M in stats_df$module) {
    # Debug: Check if the column exists in module_df
    if (!(M %in% names(module_df))) {
        print(paste("Column", M, "not found in module_df"))
        next  # Skip to the next iteration if M is not a valid column
    }

    # Construct the output tsv file path
    tsvout <- file.path("output_tsvs", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, "_genes.tsv"))

    # Debug: Check if the column exists in the data
    print(paste("Plotting for module:", M))
    print(head(module_df[, c("infection", "celltype", M)]))  # Show first few rows for verification

    # Create the file
    readr::write_tsv(gene_module_key %>% dplyr::filter(module==M), file = tsvout)

    # Optionally, print the module name for progress tracking
    print(paste("Saved tsv for module:", M))
}


################################################################################
################################################################################
################################################################################
################################## HEATMAPS ####################################
################################################################################
dir.create("output_heatmaps", showWarnings = TRUE, recursive = TRUE)

##################################
#### gene vs sample clusters
###################################

for (M in stats_df$module) {
    # Debug: Check if the column exists in module_df
    if (!(M %in% names(module_df))) {
        print(paste("Column", M, "not found in module_df"))
        next  # Skip to the next iteration if M is not a valid column
    }

    expression_mat = normalized_rld_counts
    metadata_df = metadata
    gene_module_key_df = gene_module_key
    module_eigengenes_df = module_eigengenes
    module_name = M

    module_eigengene <- module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("sample")


    # Set up column annotation from metadata
    col_annot_df <- metadata_df %>%
      # Only select the treatment and sample ID columns
      dplyr::select(sample, infection, celltype) %>%

      # Add on the eigengene expression by joining with sample IDs
      dplyr::inner_join(module_eigengene, by = "sample") %>%

      # Arrange by sample and infection
      dplyr::arrange(infection, celltype) %>%

      # Store sample
      tibble::column_to_rownames("sample")

     rownames(col_annot_df)

     # Get a vector of the Ensembl gene IDs that correspond to this module
     module_genes <- gene_module_key_df %>%
       dplyr::filter(module == module_name) %>%
       dplyr::pull(gene)

    # Create the ComplexHeatmap column annotation object
    col_annot <- ComplexHeatmap::HeatmapAnnotation(
      # Supply treatment labels
      sample = col_annot_df$sample,
      infection = col_annot_df$infection,
      celltype = col_annot_df$celltype,
      # Pick colors for each experimental group
      col = list(infection = c("uninfected" = "#f1a340", "infected" = "#998ec3"))
    )

    # Set up the gene expression data frame
    mod_mat <- as.data.frame(t(expression_mat))

    # Only keep genes from this module
    mod_mat <- mod_mat[, colnames(mod_mat) %in% module_genes]

    mod_mat <- as.matrix(mod_mat)

    # Normalize the gene expression values
    mod_mat <- mod_mat %>%
      # Scale can work on matrices, but it does it by column so we will need to
      # transpose first
      t() %>%
      scale() %>%
      # And now we need to transpose back
      t()

    # Create a color function based on standardized scale
    color_func <- circlize::colorRamp2(
      c(-2, 0, 2),
      c("#67a9cf", "#f7f7f7", "#ef8a62")
    )

    # Plot on a heatmap
    heatmap <- ComplexHeatmap::Heatmap(mod_mat,
      name = module_name,
      # Supply color function
      col = color_func,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6)
    )

    # Debug: Check if the column exists in the data
    print(paste("Heatmap plotting for module:", M))
    print(head(module_df[, c("infection", "celltype", M)]))  # Show first few rows for verification

    # Construct the output PDF file path
    # print heatmap to file
    pdfout <- file.path("output_heatmaps", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, "_heatmap.pdf"))

    # Open a PDF device to save the plot
    pdf(pdfout)

    # Print the plot to ensure it's being rendered
    print(heatmap)

    # Close the PDF device
    dev.off()

    # Optionally, print the module name for progress tracking
    print(paste("Saved plot for module:", M))
}

################################################################################
###### samples vs genes + eigengene plots
################################################################################
## double check that our samples are still in order.
if (all.equal(metadata_df$sample, module_eigengene$sample)) {
        for (M in stats_df$module) {
            # Debug: Check if the column exists in module_df
            if (!(M %in% names(module_df))) {
                print(paste("Column", M, "not found in module_df"))
                next  # Skip to the next iteration if M is not a valid column
            }

            expression_mat = normalized_rld_counts
            metadata_df = metadata
            gene_module_key_df = gene_module_key
            module_eigengenes_df = module_eigengenes
            module_name = M

            module_eigengene <- module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("sample")

            # Set up column annotation from metadata
            col_annot_df <- metadata_df %>%
              # Only select the treatment and sample ID columns
              dplyr::select(sample, infection, celltype) %>%

              # Add on the eigengene expression by joining with sample IDs
              dplyr::inner_join(module_eigengene, by = "sample") %>%

              # Arrange by sample and infection
              dplyr::arrange(infection, celltype) %>%

              # Store sample
              tibble::column_to_rownames("sample")

             # Get a vector of the Ensembl gene IDs that correspond to this module
             module_genes <- gene_module_key_df %>%
               dplyr::filter(module == module_name) %>%
               dplyr::pull(gene)

            # Create the ComplexHeatmap column annotation object
            col_annot <- ComplexHeatmap::HeatmapAnnotation(
              # Supply treatment labels
              sample = rownames(col_annot_df),
              infection = col_annot_df$infection,
              celltype = col_annot_df$celltype,

              # Add annotation barplot - does not work - Error: number of observations in bottom annotation should be as same as ncol of the matrix.
              module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(all_of(col_annot_df), all_of(module_name))),

              # Pick colors for each experimental group
              col = list(infection = c("uninfected" = "#f1a340", "wMel" = "#998ec3"))
            )

            # Set up the gene expression data frame
            mod_mat <- as.data.frame(t(expression_mat))

            # Only keep genes from this module
            mod_mat <- mod_mat[, colnames(mod_mat) %in% module_genes]

            mod_mat <- as.matrix(mod_mat)

            # Normalize the gene expression values
            mod_mat <- mod_mat %>%
              # Scale can work on matrices, but it does it by column so we will need to
              # transpose first
              t() %>%
              scale()

            # Create a color function based on standardized scale
            color_func <- circlize::colorRamp2(
              c(-2, 0, 2),
              c("#67a9cf", "#f7f7f7", "#ef8a62")
            )

            # Plot on a heatmap
            heatmap <- ComplexHeatmap::Heatmap(mod_mat,
              name = module_name,
              # Supply color function
              col = color_func,

              # We don't want to cluster samples
              cluster_columns = FALSE,
              # We don't need to show sample or gene labels
              show_row_names = FALSE,
              show_column_names = FALSE,

              # Supply column annotation
              bottom_annotation = col_annot
              #row_names_gp = gpar(fontsize = 6),

            )

            # Debug: Check if the column exists in the data
            print(paste("Heatmap plotting for module:", M))
            print(head(module_df[, c("infection", "celltype", M)]))  # Show first few rows for verification

            # Construct the output PDF file path
            # print heatmap to file
            pdfout <- file.path("output_heatmaps", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, "_heatmap.annot.pdf"))

            # Open a PDF device to save the plot
            pdf(pdfout)

            # Print the plot to ensure it's being rendered
            print(heatmap)

            # Close the PDF device
            dev.off()

            # Optionally, print the module name for progress tracking
            print(paste("Saved plot for module:", M))
        }
}



################################################################################
## clusterProfiler is now available for Drosophila (used in: https://academic.oup.com/g3journal/article/12/10/jkac212/6670625)
## https://www.cell.com/the-innovation/fulltext/S2666-6758(21)00066-7?_returnURL&
## https://guangchuangyu.github.io/software/clusterProfiler/
## https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
## Supported Analysis:
## Over-Representation Analysis
## Gene Set Enrichment Analysis
## Biological theme comparison

## PDF manual: https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf

###########################################################################################################################################
# ==> For over representation analysis (ORA), all we need is a gene vector, that is a vector of gene IDs. These gene IDs can be obtained by differential expression analysis (e.g. with the DESeq2 package).

#######################################################################################################################################
## Gene ID mapping
## https://medium.com/computational-biology/gene-id-mapping-using-r-14ff50eec9ba
## https://www.biostars.org/p/70821/ ## entrezs <- mapIds(org.Dm.eg.db, keys=names(genes), column="ENTREZID", keytype="FLYBASECG", multiVals="first")
## https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
## https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf
## https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist
## https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html#setReadable
#################################################################################
## GO Analysis
###############################################################################
## print out GO groups
dir.create("go_txt", showWarnings = TRUE, recursive = TRUE)
dir.create("go_pdf", showWarnings = TRUE, recursive = TRUE)
dir.create("module_fbgnids", showWarnings = TRUE, recursive = TRUE)

myGOs <- c("CC", "MF", "BP")
for (term in myGOs) {
        for (module in names(table(net$colors))) {
                genes=net$colors[which(net$colors==module)]
                ## get list of genes without Dmel_ prefix
                names(genes) <- sub('^Dmel_', '', names(genes))
                g <- names(genes)

                ## set readable=TRUE for gene symbols
                go <- groupGO(gene=g, OrgDb='org.Dm.eg.db', keyType='FLYBASECG', ont=term, level=3)

                #print(head(go))

                write.table(go, file=file.path("go_txt", paste0("cells_clusterProfiler_", minreads, "reads_", top_rld_genes, "_module", module, ".GOgroup.", term, ".txt")), sep="\t")
        }
}
#################################################################################
## perform GO enrichment for each module
################################################################################
## groupGO()
## Functional Profile of a gene set at specific GO level. Given a vector of genes, this function will return the GO profile at a specific level.
## enrichGO()
## GO Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment GO categories after FDR control.
## gseGO()
## Gene Set Enrichment Analysis of Gene Ontology

## go subontologies
## CC = cellular component
## MF = molecular function
## BP = biological process
## ALL = all three
################################################################################
## the background list is the full list of expressed transcripts (8322 genes)
expressed <- sub('^Dmel_', '', colnames(wgcna_rld_data))

expressed_ids<-bitr(expressed, fromType = "FLYBASECG", toType = "FLYBASE", OrgDb="org.Dm.eg.db")
write.table(expressed_ids, "background_expressed_FBGNids.txt")


enrichment_summary <- data.frame(
    module = character(),
    num_enriched_pathways = integer(),
    stringsAsFactors = FALSE
)

for (module in names(table(net$colors))) {
        genes=net$colors[which(net$colors==module)]
        ## get list of genes without Dmel_ prefix
        names(genes) <- sub('^Dmel_', '', names(genes))
        g <- names(genes)

        ## output FBGN ids for PANGEA
        ## Map FLYBASECG to FBGN
        ## We will lose some genes here because not all IDs will be converted
        ids<-bitr(g, fromType = "FLYBASECG", toType = "FLYBASE", OrgDb="org.Dm.eg.db")
        write.table(ids, file.path("module_fbgnids", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, "_gene_FBGNids.txt")))


        go <- enrichGO(gene=g, OrgDb='org.Dm.eg.db', universe=expressed, keyType='FLYBASECG', ont="ALL", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)

        #as.data.frame(go)
        #print(head(go@result))
        ## save results
        write.table(go, file=file.path("go_txt", paste0("cells_clusterProfiler_", minreads, "reads_", top_rld_genes, "_module", module, ".enrichGO.ALL.txt")), sep="\t")

        # Check for valid results
        if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
            enrichment_summary <- rbind(
                enrichment_summary,
                data.frame(
                    module = module,
                    num_enriched_pathways = nrow(as.data.frame(go))
                )
            )
        }
        else {
                ## else (is.null(go) || nrow(as.data.frame(go)) == 0)
                message(paste("No GO terms enriched for module", module))
                next
        }

        #print(class(go))
        #print(as.data.frame(go))       # Ensure the results are not empty
        #head(as.data.frame(go))        # View the first few rows
        #nrow(as.data.frame(go))        # Ensure it contains rows


        ## Set readable IDs
        rgo <- setReadable(go, OrgDb='org.Dm.eg.db')

        ## save results
        write.table(rgo, file=file.path("go_txt", paste0("cells_clusterProfiler_", minreads, "reads_", top_rld_genes, "_module", module, ".enrichGO.ALL.symbols.txt")), sep="\t")

        ## Dotplot
        pdfGOEnrichDP <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.", rld_chosen_power, ".module", module, ".GOEnrich.ALL.dotplot.pdf"))
        pdf(pdfGOEnrichDP, width=12, height=18)
        print(dotplot(rgo, showCategory=30))
        dev.off()
}


print(enrichment_summary)

pdfGObarplot <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.power", rld_chosen_power, ".GOenrich.barplot.pdf"))
pdf(pdfGObarplot)
print(ggplot(enrichment_summary, aes(x = module, y = num_enriched_pathways)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Number of Enriched Pathways per Module",
    x = "Module",
    y = "Number of Enriched Pathways"))
dev.off()

## save results
write.table(enrichment_summary, file="enriched_pathways.output")


#################################################################################


myGOs <- c("CC", "MF", "BP")
for (term in myGOs) {
        for (module in names(table(net$colors))) {
                genes=net$colors[which(net$colors==module)]
                ## get list of genes without Dmel_ prefix
                names(genes) <- sub('^Dmel_', '', names(genes))
                g <- names(genes)

                go <- enrichGO(gene=g, OrgDb='org.Dm.eg.db', universe=expressed, keyType='FLYBASECG', ont="CC", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)

                #as.data.frame(go)
                #print(head(go@result))

                # Check for valid results
                if (is.null(go) || nrow(as.data.frame(go)) == 0) {
                        message(paste("No GO terms enriched for module", module))
                        next
                }

                #print(class(go))
                #print(as.data.frame(go))       # Ensure the results are not empty
                #head(as.data.frame(go))        # View the first few rows
                #nrow(as.data.frame(go))        # Ensure it contains rows

                ## Set readable IDs
                rgo <- setReadable(go, OrgDb='org.Dm.eg.db')

                pdfGOEnrichGraph <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.", rld_chosen_power, "_module", module, ".GOgroup.", term, ".GOEnrich.graph.pdf"))
                pdf(pdfGOEnrichGraph)
                print(plotGOgraph(rgo))
                dev.off()

                ## Dotplot
                pdfGOEnrichDP <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.", rld_chosen_power, ".module", module, ".GOEnrich.ALL.dotplot.pdf"))
                pdf(pdfGOEnrichDP)
                print(dotplot(rgo, showCategory=30))
                dev.off()

                ## Gene-concept network
                pdfGOcnet <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.power", rld_chosen_power, "_module", module, ".GOgroup.", term, ".GOEnrich.cnet.pdf"))
                pdf(pdfGOcnet)
                print(cnetplot(rgo, foldChange=wgcna_rld_data))
                dev.off()

                ## Heatmap
                pdfGOheatmap <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.power", rld_chosen_power, "_module", module, ".GOgroup.", term, ".GOEnrich.heatmap.pdf"))
                pdf(pdfGOheatmap)
                print(heatplot(rgo, foldChange=wgcna_rld_data))
                dev.off()

                ## tree plot
                rgo2 <- pairwise_termsim(rgo)

                # Enhanced checks for rgo2
                if (is.null(rgo2)) {
                    message(paste("pairwise_termsim() returned NULL for module", module, "and term", term))
                    next
                }

                # Check if the object contains a valid similarity matrix
                if (is.null(rgo2@termsim) || nrow(rgo2@termsim) < 2 || ncol(rgo2@termsim) < 2) {
                    message(paste("Insufficient or missing pairwise similarity matrix in rgo2 for module", module, "and term", term))
                    next
                }

                # Check if the similarity matrix has any meaningful values
                if (all(is.na(rgo2@termsim)) || sum(rgo2@termsim, na.rm = TRUE) == 0) {
                    message(paste("pairwise_termsim matrix is empty or invalid for module", module, "and term", term))
                    next
                }

                # Check if there are enough terms to cluster
                if (nrow(rgo2@result) < 2) {
                    message(paste("Not enough enriched terms to perform clustering for module", module, "and term", term))
                    next
                }

                # Generate treeplot if all checks pass
                tryCatch({
                    pdfGOtreeplot <- file.path("go_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.power", rld_chosen_power, "_module", module, ".GOgroup.", term, ".GOEnrich.treeplot.pdf"))
                    pdf(pdfGOtreeplot)
                    print(treeplot(rgo2, cluster.params = list(method = "average")))
                    dev.off()
                }, error = function(e) {
                    message(paste("Error in treeplot for module", module, "and term", term, ":", e$message))
                })

        }
}

#################################################################################

#################################################################################
## perform KEGG enrichment for each module
################################################################################
## enrichKEGG()
## KEGG Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment KEGG categories after FDR control.
## gseKEGG()
## Gene Set Enrichment Analysis of KEGG pathways

## the scientific_name 'Drosophila melanogaster' returns wMel too.
#dmel <- search_kegg_organism('Drosophila melanogaster', by='scientific_name')
dmel <- search_kegg_organism('dme', by='kegg_code')
dim(dmel)


################################################################################
## ncbi-geneid == ENTREZID
## https://github.com/gencorefacility/r-notebooks/blob/master/ora.Rmd

## the background list is the full list of expressed transcripts (8322 genes)
expressed <- sub('^Dmel_', '', colnames(wgcna_rld_data))
expressed_ids<-bitr(expressed, fromType = "FLYBASECG", toType = "ENTREZID", OrgDb="org.Dm.eg.db")

write.table(expressed_ids, "background_expressed_entrezids.txt")

enrichment_summary <- data.frame(
    module = character(),
    num_enriched_pathways = integer(),
    stringsAsFactors = FALSE
)

mnrichment_summary <- data.frame(
    module = character(),
    num_enriched_pathways = integer(),
    stringsAsFactors = FALSE
)

dir.create("kegg_txt", showWarnings = TRUE, recursive = TRUE)
dir.create("kegg_pdf", showWarnings = TRUE, recursive = TRUE)
dir.create("module_entrezids", showWarnings = TRUE, recursive = TRUE)


for (module in names(table(net$colors))) {
        genes=net$colors[which(net$colors==module)]
        #print(length(names(genes)))
        #print(head(names(genes)))

        # Subset the matrix while retaining its structure - for rlog counts renamed by ENTREZIDs
        # Subset the matrix while retaining its structure
        wgcna_rld_subset <- wgcna_rld_data[, colnames(wgcna_rld_data) %in% names(genes), drop = FALSE]

        ## get list of genes without Dmel_ prefix
        colnames(wgcna_rld_subset) <- sub('^Dmel_', '', colnames(wgcna_rld_subset))
        g <- colnames(wgcna_rld_subset)

        ## Map FLYBASECG to ENTREZID
        ## Convert gene IDs for enrichKEGG function
        ## convert FLYBASECGs to keyType: one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
        ## We will lose some genes here because not all IDs will be converted
        ids<-bitr(g, fromType = "FLYBASECG", toType = "ENTREZID", OrgDb="org.Dm.eg.db")

        write.table(ids, file.path("module_entrezids", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, "_gene_entrezids.txt")))

        #print(head(ids))
        #print(class(ids))
        #print(mode(ids))
        #print(names(ids))

        # Ensure all FLYBASECG have corresponding ENTREZID
        if (nrow(ids) == 0) {
                message("No valid ENTREZ IDs found for module: ", module)
                next
        }

        # Filter out NA or missing mappings
        ids <- ids[!is.na(ids$ENTREZID), ]

        # Match FLYBASECG IDs to matrix column names and replace with ENTREZID
        entrez_ids <- ids$ENTREZID[match(colnames(wgcna_rld_subset), ids$FLYBASECG)]

        # Check if all columns have valid mappings
        if (any(is.na(entrez_ids))) {
                warning("Some FLYBASECG IDs could not be mapped to ENTREZIDs. Columns with NA will retain original names.")
        }

        # Replace column names with ENTREZIDs
        colnames(wgcna_rld_subset) <- entrez_ids

        # omit NA values
        kegg_genes <- na.omit(ids$ENTREZID)

        #print(length(ids$ENTREZID))
        print(length(kegg_genes))
        print(head(kegg_genes))

        print(length(expressed_ids$ENTREZID))
        print(head(expressed_ids$ENTREZID))

        kegg_organism = "dme"

        ################## KEGG pathway enrichment analysis ######################
        kk <- enrichKEGG(gene=kegg_genes, universe=expressed_ids$ENTREZID,organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
        #print(head(kk))

        if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
            enrichment_summary <- rbind(
                enrichment_summary,
                data.frame(
                    module = module,
                    num_enriched_pathways = nrow(as.data.frame(kk))
                )
            )
            output_file <- file.path("kegg_txt", paste0("module_", module, "_KEGG_results.csv"))
            write.csv(as.data.frame(kk), output_file, row.names = FALSE)

            ##Barplot
            pdfKEGGbarplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.barplot.pdf"))
            pdf(pdfKEGGbarplot)
            print(barplot(kk, showCategory = 10, font.size = 8))
            dev.off()

            ## Dotplot
            pdfKEGGdotplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.dotplot.pdf"))
            pdf(pdfKEGGdotplot)
            print(dotplot(kk, showCategory = 10, font.size = 8))
            dev.off()

            ## Category Netplot:
            ## Gene-concept network
            # categorySize can be either 'pvalue' or 'geneNum'
            pdfKEGGcnet <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.cnet.pdf"))
            pdf(pdfKEGGcnet)
            print(cnetplot(kk, categorySize="pvalue", foldChange=wgcna_rld_subset))
            dev.off()

            rkk <- setReadable(kk, OrgDb='org.Dm.eg.db', keyType='ENTREZID')
            pdfKEGGcnet <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.readable.cnet.pdf"))
            pdf(pdfKEGGcnet)
            print(cnetplot(rkk, categorySize="pvalue", foldChange=wgcna_rld_subset))
            dev.off()

            ## Heatmap
            pdfGOheatmap <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGEnrich.heatmap.pdf"))
            pdf(pdfGOheatmap)
            print(heatplot(rkk, foldChange=wgcna_rld_data))
            dev.off()


        } else {
            enrichment_summary <- rbind(
                enrichment_summary,
                data.frame(
                    module = module,
                    num_enriched_pathways = 0
                )
            )
            message("No enriched pathways found for module: ", module)
            next
        }

        ######################################################################

        ####### KEGG module over-representation analysis ####################
        mkg <- enrichMKEGG(gene=kegg_genes, universe=expressed_ids$ENTREZID,organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

        if (!is.null(mkg) && nrow(as.data.frame(mkg)) > 0) {
            mnrichment_summary <- rbind(
                mnrichment_summary,
                data.frame(
                    module = module,
                    num_enriched_pathways = nrow(as.data.frame(mkg))
                )
            )
            output_file <- file.path("kegg_txt", paste0("module_", module, "_mKEGG_results.csv"))
            write.csv(as.data.frame(mkg), output_file, row.names = FALSE)

            ##Barplot
            pdfKEGGbarplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".mKEGGenrich.barplot.pdf"))
            pdf(pdfKEGGbarplot)
            print(barplot(mkg, showCategory = 10, font.size = 8))
            dev.off()

            ## Dotplot
            pdfKEGGdotplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".mKEGGenrich.dotplot.pdf"))
            pdf(pdfKEGGdotplot)
            print(dotplot(mkg, showCategory = 10, font.size = 8))
            dev.off()

            ## Category Netplot:
            ## Gene-concept network
            # categorySize can be either 'pvalue' or 'geneNum'
            pdfKEGGcnet <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".mKEGGenrich.cnet.pdf"))
            pdf(pdfKEGGcnet)
            print(cnetplot(mkg, categorySize="pvalue", foldChange=wgcna_rld_subset))
            dev.off()

            ## Set readable IDs
            rmkg <- setReadable(mkg, OrgDb='org.Dm.eg.db', keyType='ENTREZID')
            pdfKEGGcnet <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".mKEGGenrich.readable.cnet.pdf"))
            pdf(pdfKEGGcnet)
            print(cnetplot(rmkg, categorySize="pvalue", foldChange=wgcna_rld_subset))
            dev.off()

            ## Heatmap
            pdfGOheatmap <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".mKEGGEnrich.heatmap.pdf"))
            pdf(pdfGOheatmap)
            print(heatplot(rmkg, foldChange=wgcna_rld_data))
            dev.off()


        } else {
            mnrichment_summary <- rbind(
                mnrichment_summary,
                data.frame(
                    module = module,
                    num_enriched_pathways = 0
                )
            )
            message("No enriched pathways found for module: ", module)
            next
        }

}

print(enrichment_summary)

pdfKEGGbarplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, ".KEGGenrich.barplot.pdf"))
pdf(pdfKEGGbarplot)
print(ggplot(enrichment_summary, aes(x = module, y = num_enriched_pathways)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Number of Enriched Pathways per Module",
    x = "Module",
    y = "Number of Enriched Pathways"))
dev.off()

## save results
write.table(enrichment_summary, file="enriched_pathways.output")


print(mnrichment_summary)
pdfmKEGGbarplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, ".mKEGGenrich.barplot.pdf"))
pdf(pdfmKEGGbarplot)
print(ggplot(mnrichment_summary, aes(x = module, y = num_enriched_pathways)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Module",
    y = "Number of Enriched Modules"))
dev.off()

## save results
write.table(mnrichment_summary, file="enriched_modules.output")

#################################################################################
