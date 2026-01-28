#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

# Inputs (only these are required)
deseq2_file <- get_arg("--deseq2_results", NULL)
out_root    <- get_arg("--output_dir", NULL)
species_tag <- get_arg("--species", "hg38")  # passed from main.py

if (is.null(deseq2_file) || is.null(out_root)) {
  stop("Usage: Rscript dott_gsea_preranked.R --deseq2_results <csv> --output_dir <dir> [--species hg38|mm10|mm39|hg19]")
}

# Hard-coded “defaults” requested
padj_cut <- 0.05
topN     <- 80
esN      <- 80
minSize  <- 15
maxSize  <- 500

# Output structure under user's --output-dir
outdir <- file.path(out_root, "gsea_fullrnk_sig")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "ES_plots"), showWarnings = FALSE, recursive = TRUE)

# -------------------- Read full DESeq2 results --------------------
# Your DESeq2 module writes comma-separated CSV with a 'gene' column.
res <- suppressMessages(readr::read_csv(deseq2_file, show_col_types = FALSE))

if (!("gene" %in% names(res))) stop("DESeq2 results must have a 'gene' column.")

# Choose ranking column: stat preferred (GSEA classic), else log2FoldChange
rank_col <- if ("stat" %in% names(res)) "stat" else if ("log2FoldChange" %in% names(res)) "log2FoldChange" else NA_character_
if (is.na(rank_col)) stop("DESeq2 results must contain 'stat' or 'log2FoldChange' for ranking.")

rnk <- res %>%
  transmute(gene = as.character(gene),
            score = suppressWarnings(as.numeric(.data[[rank_col]]))) %>%
  filter(!is.na(gene), gene != "", !is.na(score), is.finite(score))

# Collapse duplicate genes (keep max abs score)
rnk <- rnk %>%
  group_by(gene) %>%
  summarise(score = score[which.max(abs(score))], .groups = "drop")

# Deterministic tie-breaking (optional but helps avoid fgsea “ties” warning)
# Sort by score then gene, then add a tiny monotonic jitter (does not change rank meaningfully)
rnk <- rnk %>% arrange(desc(score), gene)
rnk$score <- rnk$score + (seq_len(nrow(rnk)) * 1e-12)

# Save RNK used
rnk_path <- file.path(outdir, "gene_list_full.sorted.rnk")
write.table(rnk[, c("gene","score")], file = rnk_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

gene_ranks <- rnk$score
names(gene_ranks) <- rnk$gene
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# -------------------- MSigDB Hallmarks --------------------
is_mouse <- grepl("^mm", tolower(species_tag))

if (is_mouse) {
  msig <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
} else {
  msig <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = "H")
}

pathways <- split(msig$gene_symbol, msig$gs_name)

# -------------------- fgsea (multilevel) --------------------
set.seed(123)
fg <- fgseaMultilevel(
  pathways = pathways,
  stats    = gene_ranks,
  minSize  = minSize,
  maxSize  = maxSize
) %>%
  arrange(padj) %>%
  mutate(neglog10padj = -log10(padj))

# Flatten list-column for CSV
fg_out <- fg %>%
  mutate(
    leadingEdge = if ("leadingEdge" %in% colnames(.))
      vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    else NA_character_
  )

write.csv(fg_out, file = file.path(outdir, "fgsea_results_all_hallmarks.csv"),
          row.names = FALSE, quote = FALSE)

# -------------------- Bubble plot: significant only --------------------
fg_top <- fg %>%
  filter(!is.na(padj), is.finite(padj), padj <= padj_cut) %>%
  arrange(padj)

# Cap to topN significant if too many
fg_top <- fg_top %>%
  slice_head(n = topN) %>%
  mutate(Pathway = factor(pathway, levels = rev(pathway)))

if (nrow(fg_top) == 0) {
  message("No pathways passed padj_cut=", padj_cut, ". Skipping bubble/ES plots.")
  quit(save = "no", status = 0)
}

bubble_title <- if (is_mouse) {
  paste0("DoTT-region (RNK) GSEA — Mouse Hallmark (MH), FDR ≤ ", padj_cut, " (rank=", rank_col, ")")
} else {
  paste0("DoTT-region (RNK) GSEA — Human Hallmark (H), FDR ≤ ", padj_cut, " (rank=", rank_col, ")")
}

bubble_plot <- ggplot(fg_top, aes(x = NES, y = Pathway, color = abs(NES), size = neglog10padj)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "blue", high = "red", name = "NES Magnitude") +
  scale_size_continuous(name = "-log10(padj)", range = c(2, 16)) +
  labs(title = bubble_title, x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# Auto-height so labels don’t crush
plot_h <- max(8, 0.22 * nrow(fg_top) + 2)
ggsave(file.path(outdir, "bubble_sig_hallmarks.svg"), bubble_plot, width = 11, height = plot_h, dpi = 300)

# -------------------- ES plots: significant only --------------------
fg_es <- fg %>%
  filter(!is.na(padj), is.finite(padj), padj <= padj_cut) %>%
  arrange(padj) %>%
  slice_head(n = esN)

for (i in seq_len(nrow(fg_es))) {
  pw <- fg_es$pathway[i]
  p <- plotEnrichment(pathways[[pw]], gene_ranks) +
    labs(title = paste0(pw, " | NES=", round(fg_es$NES[i], 3), " | padj=", signif(fg_es$padj[i], 3))) +
    theme_bw(base_size = 12)
  safe_pw <- gsub("[^A-Za-z0-9_]+", "_", pw)
  ggsave(file.path(outdir, "ES_plots", paste0("ES_", i, "_", safe_pw, ".svg")),
         plot = p, width = 9, height = 5, dpi = 300)
}

message("Done. Outputs written to: ", outdir)
message("RNK used: ", rnk_path)
