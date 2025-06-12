#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# --------------- INPUTS -------------------
quant_files      <- snakemake@input[["quant"]]
design_file      <- snakemake@input[["samples"]]
comparisons_file <- snakemake@input[["comparisons"]]
tx2gene_file     <- snakemake@input[["tx2gene"]]
filter_thr       <- as.integer(snakemake@params[["filter_count_threshold"]])
out_dir          <- snakemake@output[["results"]]
vst_mat_file     <- snakemake@output[["vst_mat"]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- LECTURE DU DESIGN -------------
design_df <- read_tsv(design_file, col_types = cols())
sampleTable <- data.frame(
  condition = factor(design_df$condition, levels = c("Control", "FXS")),
  row.names = design_df$sample
)
comparisons_df <- read_tsv(comparisons_file, col_types = cols())

# ----------- FILTRAGE protein_coding ---------
tx2gene <- read_tsv(tx2gene_file, col_types = cols())
protein_coding_genes <- unique(tx2gene$GENEID)

# ---- CONSTRUCTION MATRICE DE COUNTS ------
count_list <- lapply(quant_files, function(path) {
  df <- read_tsv(path, col_types = cols())
  if (!("gene_id" %in% colnames(df)) || !("count" %in% colnames(df))) {
    stop(sprintf("Fichier %s doit contenir 'gene_id' et 'count'", path))
  }
  if (!is.numeric(df$count)) {
    stop(sprintf("La colonne 'count' dans %s n'est pas numérique.", path))
  }
  samp <- tools::file_path_sans_ext(basename(path))
  df <- df %>% mutate(count = as.integer(round(count)))
  df <- df %>% filter(gene_id %in% protein_coding_genes)
  df %>% select(gene = gene_id, !!samp := count)
})

# Fusionne les matrices filtrées
merged_counts <- Reduce(function(x, y) full_join(x, y, by = "gene"), count_list)
gene_col <- merged_counts$gene

counts_mat <- merged_counts %>%
  select(-gene) %>%
  as.matrix()
rownames(counts_mat) <- gene_col

# Réordonne les colonnes selon sampleTable
counts_mat <- counts_mat[, rownames(sampleTable)]

# ------------- DESeq2 -------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sampleTable,
  design    = ~ condition
)
keep <- rowSums(counts(dds)) >= filter_thr
dds  <- dds[keep, ]
if (nrow(dds) == 0) stop("Aucun gène gardé après filtrage.")

dds <- DESeq(dds)
genes_kept <- rownames(dds)

# ----------- Export des résultats par comparaison -----------
for (i in seq_len(nrow(comparisons_df))) {
  cdn1  <- comparisons_df$cdn1[i]
  cdn2  <- comparisons_df$cdn2[i]
  comp  <- paste0(cdn1, "-", cdn2)
  res   <- results(dds, contrast = c("condition", cdn2, cdn1))
  resdf <- as.data.frame(res)
  resdf <- cbind(gene = genes_kept, resdf)
  out_f <- file.path(out_dir, paste0(comp, "_DESeq2_gene.csv"))
  write.csv(resdf, file = out_f, quote = FALSE, row.names = FALSE)
}

# ========== Sauvegarde matrice VST ==========
vsd <- vst(dds)
vst_mat <- assay(vsd)
vst_df <- as.data.frame(vst_mat)
vst_df$gene <- rownames(vst_df)
vst_df <- vst_df[, c("gene", setdiff(colnames(vst_df), "gene"))] # 'gene' colonne 1
write_csv(vst_df, vst_mat_file)
