#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(DESeq2)
})

# — Snakemake inputs —
quant_files      <- snakemake@input[["quant"]]
design_file      <- snakemake@input[["samples"]]
comparisons_file <- snakemake@input[["comparisons"]]
filter_thr       <- as.integer(snakemake@params[["filter_count_threshold"]])
out_dir          <- snakemake@output[["results"]]

# Créer le dossier de sortie
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# — Lecture du design et des comparaisons —
design_df <- read_tsv(design_file, col_types = cols())
sampleTable <- data.frame(
  condition = factor(design_df$condition, levels = unique(design_df$condition)),
  row.names = design_df$sample
)
comparisons_df <- read_tsv(comparisons_file, col_types = cols())

# — Chargement et fusion des counts COCO —
count_list <- lapply(quant_files, function(path) {
  df <- read_tsv(path, col_types = cols())
  # Vérifier la présence et le type des colonnes essentielles
  if (!all(c("transcript_id", "count") %in% colnames(df))) {
    stop(sprintf("Le fichier %s doit contenir les colonnes 'transcript_id' et 'count'", path))
  }
  if (!is.numeric(df$count)) {
    stop(sprintf("La colonne 'count' dans %s n'est pas numérique.", path))
  }
  samp <- tools::file_path_sans_ext(basename(path))
  # Arrondir et convertir en entier
  df <- df %>% mutate(count = as.integer(round(count)))
  # On sélectionne transcript_id (ENST), pas gene_id
  df %>% select(transcript_id, !!samp := count)
})

# Fusion en une seule table
merged_counts <- Reduce(function(x, y) full_join(x, y, by = "transcript_id"), count_list)

# Construire la matrice de counts
transcript_ids <- merged_counts$transcript_id
counts_mat <- merged_counts %>% select(-transcript_id) %>% as.matrix()
rownames(counts_mat) <- transcript_ids

# Réordonner les colonnes selon sampleTable
counts_mat <- counts_mat[, rownames(sampleTable), drop = FALSE]

# — DESeq2 —
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sampleTable,
  design    = ~ condition
)
# Filtrer transcripts à faible comptage
keep <- rowSums(counts(dds)) >= filter_thr
dds  <- dds[keep, ]
if (nrow(dds) == 0) stop("Aucun élément conservé après filtrage.")

dds <- DESeq(dds)

# Boucle sur chaque comparaison
for (i in seq_len(nrow(comparisons_df))) {
  cdn1   <- comparisons_df$cdn1[i]
  cdn2   <- comparisons_df$cdn2[i]
  comp_nm <- paste0(cdn1, "-", cdn2)
  res    <- results(dds, contrast = c("condition", cdn2, cdn1))
  res_df <- as.data.frame(res)
  # Réinjecter la colonne 'transcript_id'
  res_df <- cbind(transcript_id = rownames(res_df), res_df)
  out_f  <- file.path(out_dir, paste0(comp_nm, "_DESeq2_transcripts.csv"))
  write.csv(res_df, file = out_f, quote = FALSE, row.names = FALSE)
}
