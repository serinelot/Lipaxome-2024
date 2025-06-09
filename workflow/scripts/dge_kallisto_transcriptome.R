#!/usr/bin/env Rscript
#
# Differential expression (DESeq2) on the whole transcriptome
# (coding + non‑coding) from Kallisto HDF5 quantifications.
#

suppressPackageStartupMessages({
  library(readr)     # read_tsv
  library(tximport)  # tximport
  library(DESeq2)    # DESeq2
})

# — Snakemake inputs —
quant_files   <- snakemake@input[["quant"]]
design_file   <- snakemake@input[["samples"]]
comp_file     <- snakemake@input[["comparisons"]]
filter_thr    <- as.integer(snakemake@params[["filter_count_threshold"]])
out_dir       <- snakemake@output[["results"]]

# Create output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Load sample table
design_df <- read_tsv(design_file, col_types = cols())
sampleTable <- data.frame(
  condition = factor(design_df$condition, levels = unique(design_df$condition)),
  row.names = design_df$sample
)

# 2) Load comparisons
comparisons_df <- read_tsv(comp_file, col_types = cols())

# 3) Prepare abundance.h5 files
names(quant_files) <- basename(dirname(quant_files))
if (!all(file.exists(quant_files))) {
  stop("At least one abundance.h5 is missing under results/quant/kallisto/")
}

# 4) Import transcript-level counts (no collapsing to genes)
txi <- tximport(
  quant_files,
  type  = "kallisto",
  txOut = TRUE
)

# 5) Build DESeq2 dataset
dds <- DESeqDataSetFromTximport(
  txi,
  colData = sampleTable,
  design  = ~ condition
)

# 6) Filter out low-count transcripts
keep <- rowSums(counts(dds)) >= filter_thr
dds  <- dds[keep, ]

# 7) Run DESeq2
dds <- DESeq(dds)

# 8) Loop over comparisons and write results
for (i in seq_len(nrow(comparisons_df))) {
  cdn1    <- comparisons_df$cdn1[i]
  cdn2    <- comparisons_df$cdn2[i]
  comp_nm <- paste0(cdn1, "-", cdn2)
  res     <- results(dds, contrast = c("condition", cdn2, cdn1))
  res_df  <- as.data.frame(res)
  out_file <- file.path(
    out_dir,
    paste0(comp_nm, "_DESeq2_transcripts.csv")
  )
  write.csv(
    res_df,
    file      = out_file,
    quote     = FALSE,
    row.names = TRUE
  )
}
