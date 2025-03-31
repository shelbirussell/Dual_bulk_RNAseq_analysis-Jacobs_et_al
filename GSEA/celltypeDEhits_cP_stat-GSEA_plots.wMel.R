# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)  # Use the appropriate organism package for gene ID mapping
library(DOSE)          # For disease ontology semantic and enrichment analysis
library(ggridges)
library(enrichplot)
library(ggupset)
# Ensure required libraries are loaded
library(httr)

# Convert WD_RS# to WD_#
# Use the GTF file path from earlier in the script
gtf_file <- "/Users/shelbilrussell/My Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/reference/gtfs/GCF_000008025.1_ASM802v1_genomic.gtf"
gene_id_mapping <- convert_gene_ids(gtf_file)


## GSEA
## original GSEA paper, PNAS 2005 www.pnas.orgcgidoi10.1073pnas.0506580102
## ranked based on the correlation between their expression and the class distinction by using any suitable metric
## ==> microarrays used signal/noise
## ==> DE uses -log10(p) x sign(log2FC) to normalize pvalues by expression levels

## ==> upload DEseq csvs and rank by adjusted pvalue

######################################################
######## GO/KEGG Plots from DESeq2 #######
######################################################
DEseq_dir <- "/Users/shelbilrussell/My\ Drive/lab/projects/infection-omics/infection_transcriptomics/bulkRNAseq/data/cell_culture_RNAseq/3-differential_expression/2-normalization/DESeq2/wMel_SLR/csvs"
wald_path <- file.path(DEseq_dir)

###############################################################################################
############################## DE ~ cell type ###############################################
###############################################################################################
csv <- "cells_Kallisto-DESeq2_WaldTest-celltype_JW18_vs_S2.wMel.csv"

# Load the DESeq2 results file
csv_path <-file.path(wald_path, csv)
deseq_results <- read.csv(csv_path)

# Check the column names to identify log2 fold change and p-value columns
head(deseq_results)

# Ensure the data has columns "log2FoldChange", "pvalue", and "padj"
# Filter to keep genes with valid log2 fold change and p-value
deseq_results <- deseq_results[!is.na(deseq_results$log2FoldChange) & !is.na(deseq_results$pvalue) & !is.na(deseq_results$padj), ]

####################################################
################################################
#######################################

# Rank genes by stat (high to low order)
# First, convert the gene IDs from WD_RS# to WD_#
deseq_results$converted_gene_id <- sapply(deseq_results$X, function(x) {
  # Look up the converted ID
  converted_id <- gene_id_mapping[[x]]
  if (is.null(converted_id)) {
    return(x)  # If no mapping exists, use the original ID
  }
  # Trim any leading/trailing spaces
  converted_id <- trimws(converted_id)
  return(converted_id)
})

# Now create the gene list with converted IDs
gene_list <- deseq_results$stat
names(gene_list) <- deseq_results$converted_gene_id  # Use converted gene names
gene_list <- sort(gene_list, decreasing = TRUE)


############################
###########################
############################
### Functions
#############################

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


#############################################
########################################

# Prepare a more comprehensive enrichment summary
enrichment_summary_detailed <- data.frame(
  module = character(),
  num_enriched_pathways = integer(),
  enriched_pathways = character(),
  stringsAsFactors = FALSE
)

#########################################

# Initialize a list to collect all pathway results
all_module_results <- list()

  # Perform KEGG enrichment
if (length(gene_list) > 0) {
        tryCatch({
                # Try KEGG enrichment
                        kegg_enrich <- gseKEGG(
                        geneList = gene_list,
                        organism = "wol",
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = 0.05
                )

                # Save the full pathway results always, even if none pass significance cutoffs
                if (!is.null(kegg_enrich)) {
                        # Save all results including those that don't pass significance thresholds
                        # This helps document what was analyzed
                        kegg_all_results <- kegg_enrich@result
                        kegg_all_results$module <- module

                        # Add to our collection of all results
                        all_module_results[[module]] <- kegg_all_results

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
}

# View KEGG GSEA results
print(head(kegg_enrich))
length(kegg_enrich$ID)

# Save the results to CSV
write.csv(as.data.frame(kegg_enrich), "celltype_WaldTest-DE_GSEA_gseKEGG_results.csv")

### plot results
# Check if kegg_enrich contains any results before plotting
if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
        pdfGseKEGG <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.dotplot.pdf"))
        pdf(pdfGseKEGG)
        print(dotplot(kegg_enrich, showCategory = 20))
        dev.off()


        for (i in seq_along(kegg_enrich$ID)) {
                pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset", i, ".gseKEGG.plot.pdf"))
                pdf(pdfGse2KEGG)
                print(gseaplot2(kegg_enrich, geneSetID = i, title=kegg_enrich$Description[i]))
                dev.off()
        }


        ## relaxed p-value to get more output - only lysozyme is significant **
        ## https://www.biostars.org/p/433131/

        pdfGse2KEGG <- file.path(paste0("celltype_WaldTest-DE", ".geneset.1-", length(kegg_enrich$ID), ".gseKEGG.plot.pdf"))
        pdf(pdfGse2KEGG)
        print(gseaplot2(kegg_enrich, geneSetID = 1:length(kegg_enrich$ID), title=kegg_enrich$Description[i], pvalue_table = TRUE))
        dev.off()

        ## https://rdrr.io/bioc/enrichplot/man/ridgeplot.html
        pdfGseKEGGridge <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.ridgeplot.p0.05.pdf"))
        pdf(pdfGseKEGGridge)
        print(ridgeplot(kegg_enrich))
        dev.off()

        kegg_enrich2 <- pairwise_termsim(kegg_enrich)
        pdfGseKEGGemap <- file.path(paste0("celltype_WaldTest-DE", ".gseKEGG.emapplot.p0.05.pdf"))
        pdf(pdfGseKEGGemap)
        print(emapplot(kegg_enrich2, showCategory=20))
        dev.off()

}
