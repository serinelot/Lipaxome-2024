#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tibble)
})

# — Inputs / Outputs Snakemake —
deseq2_csv <- as.character(snakemake@input[["deseq2"]])
gtf_files  <- as.character(unlist(snakemake@input[["gtf"]]))
output_csv <- as.character(snakemake@output[["stats"]])

# 1) Création du répertoire de sortie
out_dir <- dirname(output_csv)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2) Lecture des résultats DESeq2
dt <- fread(deseq2_csv)
setnames(dt, old = names(dt)[1], new = "gene_id")

# 3) Lecture et parsing du/des GTF(s)
gtf_list <- lapply(gtf_files, function(f) {
  read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE)
})
gtf <- bind_rows(gtf_list)
colnames(gtf)[3]  <- "feature"
colnames(gtf)[9]  <- "attributes"
gtf_gene <- filter(gtf, feature == "gene")

extract_attr <- function(attrs, key) {
  str_match(attrs, paste0(key, ' "([^"]+)"'))[,2]
}

gene_table <- tibble(
  gene_id      = extract_attr(gtf_gene$attributes, "gene_id"),
  gene_name    = extract_attr(gtf_gene$attributes, "gene_name"),
  gene_biotype = extract_attr(gtf_gene$attributes, "gene_biotype")
)

# fallback si besoin
if (all(is.na(gene_table$gene_biotype))) {
  gene_table$gene_biotype <- extract_attr(gtf_gene$attributes, "gene_type")
}

gene_table <- distinct(gene_table, gene_id, .keep_all = TRUE)

# 4) Fusion et réorganisation des colonnes
dt2 <- left_join(dt, gene_table, by = "gene_id") %>% as.data.table()
all_cols  <- names(dt2)
rest_cols <- setdiff(all_cols, c("gene_name","gene_biotype"))
pos       <- match("gene_id", rest_cols)
final_cols<- append(rest_cols, c("gene_name","gene_biotype"), after = pos)
dt2       <- dt2[, ..final_cols]

# 5) Écriture du fichier de stats
fwrite(dt2, file = output_csv, sep = ",", na = "NA", quote = FALSE)

cat("✓ Gene-stats written to:", output_csv, "\n")
