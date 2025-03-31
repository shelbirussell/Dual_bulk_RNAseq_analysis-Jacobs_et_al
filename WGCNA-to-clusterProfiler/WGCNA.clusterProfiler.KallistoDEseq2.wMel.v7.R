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

## re-running the analysis on my kallisto output

kallisto_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/1-transcript_counting/kallisto/cds.gRNA_ref/pseudoaligned"
base_dir <- "/Users/shelbilrussell/Library/CloudStorage/GoogleDrive-shelbilrussell@gmail.com/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/1-transcript_counting/kallisto/cds.gRNA_ref"
sample_id <- dir(kallisto_dir)

sample_dirs <- file.path(kallisto_dir, sample_id)
sample_dirs <- sample_dirs[grepl("wMel", sample_dirs)]

metadata_path <- file.path(base_dir,"wMel.samples.txt")
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

ref_db <- file.path(base_dir, "wMel2978.6_cds.rna_from_genomic_gene_ids.tsv")
t2g <- read.table(ref_db, header=FALSE)

### create a named vector pointing to the (subset) Kallisto transcriptome quantification files
files <- file.path(kallisto_dir, metadata$run, "abundance_Wolbachia.tsv")
names(files) <- metadata$run
all(file.exists(files))

## kallisto with TSV files
## avoid gene-level summations txOut=TRUE
txi.kallisto.tsv <- tximport(files, type = "kallisto", geneIdCol = t2g, ignoreAfterBar = FALSE, txOut=TRUE)


###################################################
################################################################################

### Create a DESeqDataset ####
## use the gene-level estimated counts from the quantification tools, and additionally to use the transcript-level abundance estimates
## to calculate a gene-level offset that corrects for changes to the average transcript length across samples
## the function DESeqDataSetFromTximport takes care of creation of the offset for you

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = metadata, design = ~celltype)


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

dds$group <- factor(paste0(dds$celltype))

############################################################################
#### set "controls" - JW18 vs S2 (control) and wMel vs uninfected (control)
############################################################################
## run relevel before running DESeq() - can also use "contrast" later

dds$celltype <- relevel(dds$celltype, ref = "S2")


#################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## transformation code adapted from : https://rpubs.com/mahima_bose/921350
## Use the variance-stabilizing transformation (VST) or rlog for downstream correlation analyses:
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


pdfAutoclusters <- file.path("wgcna_output", paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "genes_soft_threshold_power", rld_chosen_power, "_analysis.ggplot.wMel.pdf"))
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
pdfAutoclusters2 <- file.path("wgcna_output", paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "genes_soft_threshold_power_connectivity", rld_chosen_power, "_analysis.ggplot.wMel.pdf"))
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

pdfAutoclusters2 <- file.path("wgcna_output", paste0("cell_WGCNA_rld_", minreads, "reads_", top_rld_genes, "genes_soft_threshold_power_log10connectivity", rld_chosen_power, "_analysis.ggplot.wMel.pdf"))
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
pdfAutoclusters <- file.path("wgcna_output", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, ".wMel.pdf"))

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

#Cluster method   : average
#Number of objects: 763

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

write.csv(modules_df, file = file.path("wgcna_output", paste0("gene_modules_", minreads, "reads_", top_rld_genes, "genes_auto.results_power", rld_chosen_power, ".wMel.txt")))


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
pdfEigengeneHeatmap <- file.path("wgcna_output", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.full.eigengene.heatmap", rld_chosen_power, ".wMel.pdf"))
pdf(pdfEigengeneHeatmap)
plotEigengeneNetworks(MEs, names(MEs), marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

## see the size of each module color-coded, run the following command:
table(net$colors)

## export the expression profiles of eigengenes in a text table, run the following command:
write.table(MEs, file="eigengenes.wMel.txt", sep="\t")

## export genes and their log2FC in each module as a separate table file in the working directory, run the following command:

for (color in net$colors) {
        module=wgcna_rld_data[which(net$colors==color)]
        write.table(module, file.path("wgcna_output", paste0("module_", color, ".wMel.txt",sep="")), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
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
        cyt = exportNetworkToCytoscape(modTOM, edgeFile = file.path("cytoscape_files", paste0("CytoscapeInput-edges-", paste(module, collapse="-"), ".wMel.txt", sep="")),
         nodeFile = file.path("cytoscape_files", paste0("CytoscapeInput-nodes-", paste(module, collapse="-"), ".wMel.txt", sep="")),
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
stats_txt <- file.path("wgcna_output", paste0("cells_Kallisto-DESeq2_WGCNA_", minreads, "reads_", top_rld_genes, "genes_auto.stats_power", rld_chosen_power, ".wMel.csv"))

######################################
### Run a linear model on each module
#####################################
########## INFECTION ###############
####################################

# Create the design matrix from the `infection` and 'celltype' variables
des_mat <- model.matrix(~metadata$celltype)

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
      dplyr::select(all_of(c("sample", "run", "celltype"))),
    by = c("sample" = "run")
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
    pdfout <- file.path("output_plots", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, ".wMel.pdf"))

    # Debug: Check if the column exists in the data
    print(paste("Plotting for module:", M))
    print(head(module_df[, c("celltype", M)]))  # Show first few rows for verification
    print(paste("Module name:", M))
    print(paste("Module exists in data:", M %in% names(module_df)))
    print(paste("Available columns:", paste(names(module_df), collapse=", ")))

    # Open a PDF device to save the plot
    pdf(pdfout)

    # Create the plot
    plot <- ggplot(module_df, aes(x = celltype, y = !!sym(M), color = celltype)) +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        ggforce::geom_sina(maxwidth = 0.3) +
        theme_classic()

    # Print the plot to ensure it's being rendered
    print(plot)

    # Close the PDF device
    dev.off()

    # Optionally, print the module name for progress tracking
    print(paste("Module eigengene column names:", paste(colnames(module_eigengenes), collapse=", ")))
    print(paste("Saved plot for module:", M))
}


###############################################################################
###############################################################################

## What genes are a part of each module?
## If you want to know which of your genes make up a modules, you can look at the $colors slot. This is a named list which associates the genes with the module they are a part of. We can turn this into a data frame for handy use.

# Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
gene_module_key <- tibble::enframe(net$colors, name = "gene", value = "module") %>%  dplyr::mutate(module = paste0("ME", module))
readr::write_tsv(gene_module_key, file = file.path("wgcna_output", paste0("cells_wgcna_", minreads, "reads_", top_rld_genes, "genes_.auto_gene_to_module.wMel.tsv")))

dir.create("output_tsvs", showWarnings = TRUE, recursive = TRUE)

for (M in stats_df$module) {
    # Debug: Check if the column exists in module_df
    if (!(M %in% names(module_df))) {
        print(paste("Column", M, "not found in module_df"))
        next  # Skip to the next iteration if M is not a valid column
    }

    # Construct the output tsv file path
    tsvout <- file.path("output_tsvs", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, "_genes.wMel.tsv"))

    # Debug: Check if the column exists in the data
    print(paste("Plotting for module:", M))
    print(head(module_df[, c("celltype", M)]))  # Show first few rows for verification

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
      dplyr::select(sample, celltype) %>%

      # Add on the eigengene expression by joining with sample IDs
      dplyr::inner_join(module_eigengene, by = "sample") %>%

      # Arrange by sample and infection
      dplyr::arrange(celltype) %>%

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
      celltype = col_annot_df$celltype,
      # Pick colors for each experimental group
      col = list(celltype = c("S2" = "#f1a340", "JW18" = "#998ec3"))
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
    print(head(module_df[, c("celltype", M)]))  # Show first few rows for verification

    # Construct the output PDF file path
    # print heatmap to file
    pdfout <- file.path("output_heatmaps", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, "_heatmap.wMel.pdf"))

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
if (all.equal(metadata_df$run, module_eigengene$sample)) {
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
              dplyr::select(sample = run, celltype) %>%
              dplyr::inner_join(module_eigengene, by = c("sample" = "sample")) %>%
              dplyr::arrange(celltype) %>%
              tibble::column_to_rownames("sample")

             # Get a vector of the Ensembl gene IDs that correspond to this module
             module_genes <- gene_module_key_df %>%
               dplyr::filter(module == module_name) %>%
               dplyr::pull(gene)

            # Create the ComplexHeatmap column annotation object
            col_annot <- ComplexHeatmap::HeatmapAnnotation(
              # Supply treatment labels
              sample = rownames(col_annot_df),
              celltype = col_annot_df$celltype,

              # Add annotation barplot - does not work - Error: number of observations in bottom annotation should be as same as ncol of the matrix.
              module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(all_of(col_annot_df), all_of(module_name))),

              # Pick colors for each experimental group
              col = list(celltype = c("S2" = "#f1a340", "JW18" = "#998ec3"))
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
            print(head(module_df[, c("celltype", M)]))  # Show first few rows for verification

            # Construct the output PDF file path
            # print heatmap to file
            pdfout <- file.path("output_heatmaps", paste0("cell_WGCNA_rld", minreads, "reads_", top_rld_genes, "genes_blockwise.clusters_power", rld_chosen_power, "_", M, "_heatmap.annot.wMel.pdf"))

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

#################################################################################
## perform KEGG enrichment for each module
################################################################################
## enrichKEGG()
## KEGG Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment KEGG categories after FDR control.
## gseKEGG()
## Gene Set Enrichment Analysis of KEGG pathways

## the scientific_name 'Drosophila melanogaster' returns wMel too.
#dmel <- search_kegg_organism('Drosophila melanogaster', by='scientific_name')
wol <- search_kegg_organism('wol', by='kegg_code')
#wol <-  search_kegg_organism('weg', by='kegg_code')
dim(wol)


################################################################################
## ncbi-geneid == ENTREZID
## https://github.com/gencorefacility/r-notebooks/blob/master/ora.Rmd

## I may have an issue using the KEGG wol database
## https://github.com/YuLab-SMU/clusterProfiler/issues/572

## Or KEGG without a database?
## https://support.bioconductor.org/p/9155326/
## https://www.genome.jp/genome/wol
## https://www.genome.jp/kegg/tables/br08606.html
## https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/025/GCA_000008025.1_ASM802v1/

write.table(colnames(wgcna_rld_data), "background_expressed_WD_RSids.txt")

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
dir.create("module_geneids", showWarnings = TRUE, recursive = TRUE)

### written with help from Claude ###

# Function to convert WD_RS# to WD_# using the GTF file
convert_gene_ids <- function(gtf_file) {
  # Read GTF file
  gtf_data <- read.delim(gtf_file, header = FALSE, comment.char = "#", sep = "\t", fill = TRUE)

  # Extract gene information
  gene_info <- gtf_data[gtf_data$V3 == "gene", ]

  # Create a mapping dictionary
  gene_id_map <- list()

  for (i in 1:nrow(gene_info)) {
    info_fields <- strsplit(as.character(gene_info[i, 9]), ";")[[1]]

    # Extract gene_id (WD_RS#)
    gene_id_match <- grep("gene_id ", info_fields, value = TRUE)
    if (length(gene_id_match) > 0) {
      gene_id <- gsub('gene_id "?|"?', "", gene_id_match)

      # Look for old_locus_tag
      old_locus_tag_match <- grep("old_locus_tag ", info_fields, value = TRUE)
      if (length(old_locus_tag_match) > 0) {
        old_locus_tag <- gsub('old_locus_tag "?|"?', "", old_locus_tag_match)

        # Split if multiple old locus tags exist
        old_locus_tags <- unlist(strsplit(old_locus_tag, ","))

        # Use the second tag (WD_#) if it exists
        if (length(old_locus_tags) > 1) {
          gene_id_map[[gene_id]] <- old_locus_tags[2]
        } else {
          gene_id_map[[gene_id]] <- old_locus_tag
        }
      } else {
        # If no old locus tag, keep the original gene_id
        gene_id_map[[gene_id]] <- gene_id
      }
    }
  }

  return(gene_id_map)
}

# Function to fetch KEGG gene information
fetch_kegg_gene_info <- function(gene_ids) {
  gene_info <- list()

  for (gene_id in gene_ids) {
    # Construct KEGG API URL
    url <- paste0("http://rest.kegg.jp/get/wol:", gene_id)

    # Try to fetch gene information
    tryCatch({
      # Use httr for more robust HTTP requests
      response <- httr::GET(url)

      # Check if request was successful
      if (httr::status_code(response) == 200) {
        content <- httr::content(response, "text", encoding = "UTF-8")

        # Extract gene name and description
        name_match <- regexpr("NAME\\s+(.+)", content)
        if (name_match != -1) {
          gene_name <- regmatches(content, name_match)
          gene_name <- sub("NAME\\s+", "", gene_name)
        } else {
          gene_name <- "Unknown"
        }

        # Store information
        gene_info[[gene_id]] <- list(
          name = gene_name
        )
      }
    }, error = function(e) {
      # If there's an error, store NA
      gene_info[[gene_id]] <- list(name = "Not Found")
    })
  }

  return(gene_info)
}

# Perform KEGG enrichment with additional gene information
# Add this block after your existing KEGG enrichment code

# Ensure required libraries are loaded
library(httr)

# Convert WD_RS# to WD_#
# Use the GTF file path from earlier in the script
gtf_file <- "/Users/shelbilrussell/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/reference/gtfs/GCF_000008025.1_ASM802v1_genomic.gtf"
gene_id_mapping <- convert_gene_ids(gtf_file)

write.table(gene_id_mapping, "background_expressed_oldWD_ids.txt", col.names = FALSE, row.names = FALSE, sep = "\n")

# Prepare a more comprehensive enrichment summary
enrichment_summary_detailed <- data.frame(
  module = character(),
  num_enriched_pathways = integer(),
  enriched_pathways = character(),
  stringsAsFactors = FALSE
)

# Direct access function to handle the KEGG enrichment results properly
extract_kegg_results <- function(kegg_enrich, module) {
  if (is.null(kegg_enrich)) {
    return(NULL)
  }

  # Directly access the result slot, which contains all pathway information
  if (inherits(kegg_enrich, "enrichResult") && nrow(kegg_enrich@result) > 0) {
    result_df <- kegg_enrich@result
    # Add module information
    result_df$module <- module
    return(result_df)
  }

  return(NULL)
}

# Enhanced barplot function using direct result slot access
enhanced_barplot <- function(kegg_enrich, module, showCategory = 10, font.size = 8) {
  if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
    # Create empty plot with message
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Module", module, ": No significant pathways at p < 0.05")) +
      theme_void()
    return(p)
  }

  # Get results directly from the result slot
  result_df <- kegg_enrich@result

  # Limit to top pathways by p-value
  result_df <- result_df[order(result_df$pvalue), ]
  if (nrow(result_df) > showCategory) {
    result_df <- result_df[1:showCategory, ]
  }

  # Create barplot
  p <- ggplot(result_df, aes(x = reorder(Description, -pvalue), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Pathway", y = "Gene Count",
         title = paste("Top KEGG Pathways - Module", module),
         subtitle = paste("p-value < 0.05, q-value cutoff =", 0.2)) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = font.size),
          axis.text.x = element_text(size = font.size),
          axis.title = element_text(size = font.size + 2),
          plot.title = element_text(size = font.size + 4))

  return(p)
}

# Enhanced dotplot function using direct result slot access
enhanced_dotplot <- function(kegg_enrich, module, showCategory = 10, font.size = 8) {
  if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
    # Create empty plot with message
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Module", module, ": No significant pathways at q < 0.2")) +
      theme_void()
    return(p)
  }

  # Get results directly from the result slot
  result_df <- kegg_enrich@result

  # Process GeneRatio
  result_df$GeneRatio_num <- sapply(result_df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    if (length(parts) == 2) {
      return(parts[1] / parts[2])
    } else {
      return(NA)
    }
  })

  # Limit to top pathways by p-value
  result_df <- result_df[order(result_df$pvalue), ]
  if (nrow(result_df) > showCategory) {
    result_df <- result_df[1:showCategory, ]
  }

  # Create dotplot
  p <- ggplot(result_df, aes(x = reorder(Description, -pvalue), y = GeneRatio_num)) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradient(low = "red", high = "blue") +
    coord_flip() +
    labs(x = "Pathway", y = "Gene Ratio",
         title = paste("Top KEGG Pathways - Module", module),
         subtitle = paste("Size indicates gene count, color indicates p-value")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = font.size),
          axis.text.x = element_text(size = font.size),
          axis.title = element_text(size = font.size + 2),
          plot.title = element_text(size = font.size + 4))

  return(p)
}

# Initialize a list to collect all pathway results
all_module_results <- list()

# Modify your KEGG enrichment loop
for (module in names(table(net$colors))) {
  # Get genes in this module
  genes <- net$colors[which(net$colors == module)]
  gene_names <- names(genes)

  print(paste("Processing module", module, "with", length(gene_names), "genes"))

  # Convert WD_RS# to WD_# with improved error handling
  converted_genes <- sapply(gene_names, function(x) {
    # Look up the converted ID
    converted_id <- gene_id_mapping[[x]]
    if (is.null(converted_id)) {
      return(x)  # If no mapping exists, use the original ID
    }

    # Trim any leading/trailing spaces
    converted_id <- trimws(converted_id)

    return(converted_id)
  })

  # Perform KEGG enrichment
  if (length(converted_genes) > 0) {

  write.table(converted_genes, file.path("module_geneids", paste0("module", module, "_background_expressed_oldWD_ids.txt")), col.names = FALSE, row.names = FALSE, sep = "\n")


    tryCatch({
      # Try KEGG enrichment
      kegg_enrich <- enrichKEGG(
        gene = converted_genes,
        organism = "wol",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )

      # Save the full pathway results always, even if none pass significance cutoffs
      if (!is.null(kegg_enrich)) {
        # Save all results including those that don't pass significance thresholds
        # This helps document what was analyzed
        kegg_all_results <- kegg_enrich@result
        kegg_all_results$module <- module

        # Add to our collection of all results
        all_module_results[[module]] <- kegg_all_results

        enrichment_summary <- rbind(
            enrichment_summary,
            data.frame(
                module = module,
                num_enriched_pathways = nrow(as.data.frame(kegg_all_results))
            )
        )

        # Save the full results to file
        write.table(
          kegg_all_results,
          file = file.path("kegg_txt", paste0("cells_clusterProfiler_", minreads, "reads_", top_rld_genes, "_module", module, "_allPathways.wMel.txt")),
          sep = "\t", row.names = FALSE
        )

        # If there are significant results, also save those specifically
        if (nrow(kegg_enrich@result) > 0) {
          write.table(
            kegg_enrich@result,
            file = file.path("kegg_txt", paste0("cells_clusterProfiler_", minreads, "reads_", top_rld_genes, "_module", module, ".KEGGenrich.wMel.txt")),
            sep = "\t", row.names = FALSE
          )

          # Create visualization using enhanced functions
          tryCatch({
            # Barplot
            pdfKEGGbarplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.barplot.wMel.pdf"))
            pdf(pdfKEGGbarplot)
            print(enhanced_barplot(kegg_enrich, module, showCategory = min(10, nrow(kegg_enrich@result)), font.size = 8))
            dev.off()

            # Dotplot
            pdfKEGGdotplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.dotplot.wMel.pdf"))
            pdf(pdfKEGGdotplot)
            print(enhanced_dotplot(kegg_enrich, module, showCategory = min(10, nrow(kegg_enrich@result)), font.size = 8))
            dev.off()

            # If not one of the known problematic modules, try the complex plots
            if (!module %in% c("1", "4")) {
              tryCatch({
                # Try the standard complex plots
                pathway_genes <- unique(unlist(strsplit(kegg_enrich@result$geneID, "/")))

                # Create reverse mapping (KEGG ID -> original gene ID)
                reverse_mapping <- setNames(names(converted_genes), converted_genes)
                matching_genes <- reverse_mapping[pathway_genes]
                matching_genes <- matching_genes[!is.na(matching_genes)]

                if (length(matching_genes) > 0) {
                  # Subset the expression data
                  viz_data <- wgcna_rld_data[, colnames(wgcna_rld_data) %in% matching_genes, drop = FALSE]

                  if (ncol(viz_data) > 0) {
                    # cnetplot
                    pdfKEGGcnet <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.cnet.wMel.pdf"))
                    pdf(pdfKEGGcnet)
                    print(cnetplot(kegg_enrich, categorySize="pvalue", foldChange=viz_data))
                    dev.off()

                    # heatplot
                    pdfGOheatmap <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGEnrich.heatmap.wMel.pdf"))
                    pdf(pdfGOheatmap)
                    print(heatplot(kegg_enrich, foldChange=viz_data))
                    dev.off()
                  }
                }
              }, error = function(e) {
                print(paste("Error with complex plots for module", module, ":", e$message))
              })
            }
          }, error = function(e) {
            print(paste("Error creating visualizations for module", module, ":", e$message))
          })

          # Update the summary with significant pathways
          enrichment_summary_detailed <- rbind(
            enrichment_summary_detailed,
            data.frame(
              module = module,
              num_enriched_pathways = nrow(kegg_enrich@result),
              enriched_pathways = paste(head(kegg_enrich@result$Description, 3), collapse = "; "),
              min_pvalue = min(kegg_enrich@result$pvalue),
              min_qvalue = min(kegg_enrich@result$qvalue)
            )
          )
        } else {
          print(paste("No significant KEGG enrichment found for module", module))

          # Create empty plot with informative message
          pdfKEGGbarplot <- file.path("kegg_pdf", paste0("cell_clusterProfiler_rld", minreads, "reads_", "genes_blockwise.full.power", rld_chosen_power, top_rld_genes, "genes_module", module, ".KEGGenrich.barplot.wMel.pdf"))
          pdf(pdfKEGGbarplot)
          p <- ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                    label = paste("Module", module, ": No significant KEGG pathways found")) +
            theme_void()
          print(p)
          dev.off()

          # Add a row to the summary indicating no enrichment
          enrichment_summary_detailed <- rbind(
            enrichment_summary_detailed,
            data.frame(
              module = module,
              num_enriched_pathways = 0,
              enriched_pathways = "None",
              min_pvalue = NA,
              min_qvalue = NA
            )
          )
        }
      }
    }, error = function(e) {
      print(paste("Error processing KEGG enrichment for module", module, ":", e$message))

      # Add error information to the summary
      enrichment_summary_detailed <- rbind(
        enrichment_summary_detailed,
        data.frame(
          module = module,
          num_enriched_pathways = NA,
          enriched_pathways = paste("Error:", e$message),
          min_pvalue = NA,
          min_qvalue = NA
        )
      )
    })
  } else {
    print(paste("No genes to analyze for module", module))
  }
}

# Combine all results into a single data frame and save
if (length(all_module_results) > 0) {
  all_pathways_df <- do.call(rbind, all_module_results)

  # Save the combined results with significance information
  write.csv(
    all_pathways_df,
    file = file.path("kegg_txt", paste0("all_modules_KEGG_pathways_summary.wMel.csv")),
    row.names = FALSE
  )

  # Create a version with only pathways that pass significance thresholds
  sig_pathways_df <- all_pathways_df[all_pathways_df$pvalue < 0.05, ]
  write.csv(
    sig_pathways_df,
    file = file.path("kegg_txt", paste0("significant_KEGG_pathways_summary.wMel.csv")),
    row.names = FALSE
  )
}

# Save the enrichment summary
write.csv(
  enrichment_summary_detailed,
  file = file.path("kegg_txt", "kegg_enrichment_summary_detailed.csv"),
  row.names = FALSE
)

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
