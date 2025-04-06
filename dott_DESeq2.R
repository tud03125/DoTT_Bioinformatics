#!/usr/bin/env Rscript
library(argparse)
library(DESeq2)
library(EnhancedVolcano)

# Set up argument parsing
parser <- ArgumentParser()
parser$add_argument("--counts_file", required=TRUE, help="Cleaned counts file")
parser$add_argument("--output_dir", required=TRUE, help="Output directory")
parser$add_argument("--conditions", required=TRUE, help="Comma-separated list of condition labels")
parser$add_argument("--bootstrap", type="logical", default=FALSE, help="Whether to perform bootstrapping")
parser$add_argument("--n_boot", type="integer", default=100, help="Number of bootstrap iterations")
parser$add_argument("--consensus_threshold", type="double", default=0.5, help="Fraction threshold for consensus calls")
args <- parser$parse_args()

# Parse the conditions argument into a vector and trim whitespace
conditions <- unlist(strsplit(args$conditions, ","))
conditions <- trimws(conditions)

# Load the count data
counts <- read.delim(args$counts_file, row.names = 1)

# Check that the number of conditions matches the number of samples
if(length(conditions) != ncol(counts)){
  stop("The number of conditions provided does not match the number of samples in the counts file.")
}

# Create the coldata data frame using the provided conditions
coldata <- data.frame(row.names = colnames(counts), condition = factor(conditions))

# Create DESeq2 dataset and run the analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Write DESeq2 results to CSV (without quotation marks)
results_file <- file.path(args$output_dir, "3_UTR_extended_differential_analysis_results.csv")
write.csv(as.data.frame(res), file = results_file)
cat("DESeq2 results saved to", results_file, "\n")

# Generate MA Plot
ma_plot_file <- file.path(args$output_dir, "3_UTR_extended_MA_plot.svg")
svg(ma_plot_file, width = 8, height = 6)
plotMA(res, main = "MA Plot", ylim = c(-2, 2))
dev.off()

# Generate Volcano Plot
volcano_plot_file <- file.path(args$output_dir, "3_UTR_extended_Volcano_plot.svg")
svg(volcano_plot_file, width = 8, height = 6)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Differential 3UTR Expression',
                pCutoff = 0.05)
dev.off()

# Save significant genes based on DESeq2 results:
# Use absolute log2FoldChange > 1, padj < 0.05
signif_file <- file.path(args$output_dir, "significant_extended_genes.csv")
significant_genes <- as.data.frame(res)[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, ]
write.csv(significant_genes, file = signif_file, quote = FALSE)
cat("Significant genes saved to", signif_file, "\n")

# Normalize counts and save the normalized counts to a CSV file for unsupervised ML
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
norm_counts_file <- file.path(args$output_dir, "normalized_counts.csv")
write.csv(normalized_counts, file = norm_counts_file, quote = FALSE)
cat("Normalized counts saved to", norm_counts_file, "\n")

# Optional Bootstrapping step
if (args$bootstrap) {
  cat("Performing bootstrapping with", args$n_boot, "iterations...\n")
  gene_counts <- list()
  
  # Identify indices for each condition
  cond_levels <- levels(coldata$condition)
  if(length(cond_levels) != 2){
    stop("Bootstrapping currently supports exactly 2 conditions.")
  }
  control_idx <- which(coldata$condition == cond_levels[1])
  experimental_idx <- which(coldata$condition == cond_levels[2])
  
  for (i in 1:args$n_boot) {
    boot_control <- sample(control_idx, length(control_idx), replace=TRUE)
    boot_experimental <- sample(experimental_idx, length(experimental_idx), replace=TRUE)
    boot_idx <- c(boot_control, boot_experimental)
    
    boot_counts <- counts[, boot_idx, drop=FALSE]
    boot_coldata <- coldata[boot_idx, , drop=FALSE]
    
    # Run DESeq2 on bootstrapped data
    dds_boot <- DESeqDataSetFromMatrix(countData=boot_counts,
                                       colData=boot_coldata,
                                       design=~ condition)
    dds_boot <- DESeq(dds_boot)
    res_boot <- results(dds_boot)
    res_boot_df <- as.data.frame(res_boot)
    
    # Identify significant genes in bootstrap iteration
    sig_boot <- rownames(res_boot_df)[which(!is.na(res_boot_df$padj) & 
                                              res_boot_df$padj < 0.05 & 
                                              abs(res_boot_df$log2FoldChange) > 1)]
    for (gene in sig_boot) {
      gene_counts[[gene]] <- ifelse(is.null(gene_counts[[gene]]), 1, gene_counts[[gene]] + 1)
    }
  }
  # Determine consensus genes (present in at least consensus_threshold * n_boot iterations)
  consensus_genes <- names(gene_counts)[sapply(gene_counts, function(x) x >= args$consensus_threshold * args$n_boot)]
  consensus_file <- file.path(args$output_dir, "consensus_deseq2_genes_bootstrap.txt")
  write(consensus_genes, file=consensus_file)
  cat("Consensus DE genes from bootstrapping saved to", consensus_file, "\n")
}

# --- Begin additional filtering for DOTT-like genes ---
# Combine DESeq2 results with computed normalized mean counts
conds <- levels(coldata$condition)
if(length(conds) != 2) {
  stop("This script expects exactly 2 unique conditions for filtering DOTT-like genes.")
}
ctrl <- conds[1]
exp <- conds[2]

ctrl_samples <- colnames(normalized_counts)[coldata$condition == ctrl]
exp_samples <- colnames(normalized_counts)[coldata$condition == exp]

ctrl_mean <- rowMeans(normalized_counts[, ctrl_samples], na.rm = TRUE)
exp_mean <- rowMeans(normalized_counts[, exp_samples], na.rm = TRUE)

results_with_counts <- cbind(as.data.frame(res), ctrl_mean, exp_mean)

# Define predicted DOTT-like genes based solely on DESeq2 significance and fold change:
dott_genes <- results_with_counts[
  !is.na(results_with_counts$padj) & results_with_counts$padj < 0.05 &
  abs(results_with_counts$log2FoldChange) > 1, ]
  
# Save the DOTT-like genes with individual means to a CSV file
# (Updated file name to reflect absolute values)
dott_file <- file.path(args$output_dir, "absolute_significant_extended_genes_with_individual_means.csv")
write.csv(dott_genes, file = dott_file, quote = FALSE)
