#!/usr/bin/env Rscript
library(argparse)

# Set up argument parser
parser <- ArgumentParser()
parser$add_argument("--input_file", required=TRUE, help="Input file containing gene ranking metrics (e.g., gene, log2FoldChange)")
parser$add_argument("--output_dir", required=TRUE, help="Output directory for the pre-ranked list")
args <- parser$parse_args()

# Read the input gene data.
#gene_data <- read.csv(args$input_file, header=TRUE, stringsAsFactors=FALSE)
gene_data <- read.csv(args$input_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)

# Check if the 'gene' column exists; if not, create it from row names.
if (!("gene" %in% colnames(gene_data))) {
  cat("No 'gene' column found. Creating one from row names...\n")
  gene_data <- cbind(gene = rownames(gene_data), gene_data)
  rownames(gene_data) <- NULL
}

# Check if the expected 'log2FoldChange' column exists
if (!("log2FoldChange" %in% colnames(gene_data))) {
  stop("The input file must contain a 'log2FoldChange' column.")
}

# Create the ranked list (using only the two necessary columns)
gene_list <- gene_data[, c("gene", "log2FoldChange")]

# Save the pre-ranked list in .rnk format (no header, tab-delimited)
output_file <- file.path(args$output_dir, "gene_list_converted.rnk")
write.table(gene_list, file=output_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("GSEA pre-ranked list saved to:", output_file, "\n")
cat("Instructions for GenePattern upload:\n")
cat("  1) For each category, select 'all'.\n")
cat("  2) Choose the 'No Collapse' option.\n")
