#!/usr/bin/env Rscript
#
# Differential expression (DESeq2) on protein-coding genes
# from Kallisto HDF5 (abundance.h5) quantifications.
#

suppressPackageStartupMessages({
  library(readr)
  library(tximport)
  library(DESeq2)
})

# — Snakemake inputs —
quant_files   <- snakemake@input[["quant"]]
design_file   <- snakemake@input[["samples"]]
comp_file     <- snakemake@input[["comparisons"]]
tx2gene_file  <- snakemake@input[["tx2gene"]]
filter_thr    <- as.integer(snakemake@params[["filter_count_threshold"]])
out_dir       <- snakemake@output[["results"]]

# Création du dossier de sortie
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Lecture du design (échantillon → condition)
design_df <- read_tsv(design_file, col_types = cols())
sampleTable <- data.frame(
  condition = factor(design_df$condition, levels = unique(design_df$condition)),
  row.names  = design_df$sample
)

# Lecture des comparaisons
comparisons_df <- read_tsv(comp_file, col_types = cols())

# Lecture du mapping transcript -> gene
tx2gene <- read_tsv(
  tx2gene_file,
  col_names = c("TXNAME","GENEID"),
  col_types = cols()
)

# Préparation des fichiers abundance.h5
names(quant_files) <- basename(dirname(quant_files))
if(!all(file.exists(quant_files))){
  stop("Il manque au moins un fichier abundance.h5 dans results/quant/kallisto/*/")
}

# Import des données via tximport
txi <- tximport(
  quant_files,
  type    = "kallisto",
  tx2gene = tx2gene
)

# Construction du DESeq2 dataset
dds <- DESeqDataSetFromTximport(
  txi,
  colData = sampleTable,
  design  = ~ condition
)

# Filtrage des gènes à faible comptage
keep <- rowSums(counts(dds)) >= filter_thr
dds  <- dds[keep, ]

# Exécution de DESeq2
dds <- DESeq(dds)

# Boucle sur chaque comparaison et écriture des CSV
for(i in seq_len(nrow(comparisons_df))) {
  cdn1    <- comparisons_df$cdn1[i]
  cdn2    <- comparisons_df$cdn2[i]
  comp_nm <- paste0(cdn1, "-", cdn2)
  res     <- results(dds, contrast = c("condition", cdn2, cdn1))
  out_f   <- file.path(
    out_dir,
    paste0(comp_nm, "_DESeq2_gene.csv")
  )
  write.csv(
    as.data.frame(res),
    file      = out_f,
    quote     = FALSE,
    row.names = TRUE
  )
}
