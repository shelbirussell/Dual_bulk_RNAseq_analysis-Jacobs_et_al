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

## normalization happens automatically in DEseq(), below, so this command is not necessary (unless making a gct file) #dds <- estimateSizeFactors(dds)
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


################################################################################
################################################################################
################################################################################
################################################################################
