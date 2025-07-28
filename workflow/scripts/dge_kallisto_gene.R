#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)    # read_tsv
  library(tximport) # tximport
  library(DESeq2)   # DESeq2
  library(dplyr)    # select, %>%
  library(tibble)   # rownames_to_column
})

# — Snakemake inputs & params —
quant_files      <- snakemake@input[["quant"]]
design_file      <- snakemake@input[["samples"]]
comparisons_file <- snakemake@input[["comparisons"]]
tx2gene_file     <- snakemake@input[["tx2gene"]]
filter_thr       <- as.integer(snakemake@params[["filter_count_threshold"]])
out_files        <- snakemake@output[["out_files"]]

# 1) Création des dossiers parents
for (f in out_files) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
}

# 2) Lecture du design
design_df <- read_tsv(design_file, col_types = cols())
sampleTable <- data.frame(
  condition = factor(design_df$condition, levels = unique(design_df$condition)),
  row.names = design_df$sample
)

# 3) Lecture des comparaisons
comparisons_df <- read_tsv(comparisons_file, col_types = cols())

# 4) Lecture du mapping transcript→gene (tous biotypes)
tx2gene_df <- read_tsv(tx2gene_file, col_types = cols()) %>%
  select(transcript_id, gene_id)

# 5) Préparation des fichiers Kallisto pour tximport
names(quant_files) <- basename(dirname(quant_files))
if (!all(file.exists(quant_files))) {
  stop("Au moins un fichier abundance.h5 est introuvable sous results/quant/kallisto/")
}

# 6) Import & agrégation au niveau gène
txi <- tximport(
  files           = quant_files,
  type            = "kallisto",
  tx2gene         = tx2gene_df,
  ignoreTxVersion = TRUE
)

# 7) Construction du DESeq2DataSet
dds <- DESeqDataSetFromTximport(
  txi,
  colData = sampleTable,
  design  = ~ condition
)

# 8) Filtrage des gènes à faible comptage
keep <- rowSums(counts(dds)) >= filter_thr
dds  <- dds[keep, ]
if (nrow(dds) == 0) stop("Aucun gène conservé après filtrage.")

# 9) Exécution de DESeq2
dds <- DESeq(dds)

# 10) Boucle sur les comparaisons et écriture des CSV
for (i in seq_len(nrow(comparisons_df))) {
  cdn1   <- comparisons_df$cdn1[i]
  cdn2   <- comparisons_df$cdn2[i]
  comp_nm <- paste0(cdn1, "-", cdn2)

  res     <- results(dds, contrast = c("condition", cdn2, cdn1))
  res_df  <- as.data.frame(res) %>% rownames_to_column("gene_id")

  write.csv(
    res_df,
    file      = out_files[i],
    quote     = FALSE,
    row.names = FALSE
  )
}

cat("✓ Gene-level DESeq2 terminé, résultats écrits dans :\n", 
    paste(out_files, collapse = "\n"), "\n")
