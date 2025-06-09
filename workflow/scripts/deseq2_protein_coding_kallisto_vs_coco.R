#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Entrées Snakemake
coco_file     <- snakemake@input[["coco_stat"]]
kallisto_file <- snakemake@input[["kallisto_stat"]]
gtf_file      <- snakemake@input[["gtf"]]
out_file      <- snakemake@output[["summary"]]

# Lecture des fichiers stats COCO et KALLISTO
coco      <- read.csv(coco_file,      stringsAsFactors = FALSE)
kallisto  <- read.csv(kallisto_file,  stringsAsFactors = FALSE)

# On garde SEULEMENT les gènes avec padj <= 0.05
coco_sig     <- coco     %>% filter(!is.na(padj), padj <= 0.05)
kallisto_sig <- kallisto %>% filter(!is.na(padj), padj <= 0.05)

# Sélection et renommage pour la jointure finale
coco_sig <- coco_sig %>%
  select(gene, gene_symbol, baseMean, log2FoldChange, padj) %>%
  distinct() %>%
  rename(
    baseMean_coco        = baseMean,
    log2FoldChange_coco  = log2FoldChange,
    padj_coco            = padj
  )

kallisto_sig <- kallisto_sig %>%
  select(gene, gene_symbol, baseMean, log2FoldChange, padj) %>%
  distinct() %>%
  rename(
    baseMean_kallisto        = baseMean,
    log2FoldChange_kallisto  = log2FoldChange,
    padj_kallisto            = padj
  )

# Jointure complète sur "gene" et "gene_symbol"
summary_tab <- full_join(
  kallisto_sig,
  coco_sig,
  by = c("gene", "gene_symbol")
)

# Renommer la colonne 'gene' en 'ensembl'
colnames(summary_tab)[colnames(summary_tab) == "gene"] <- "ensembl"

# Ajout des fold change linéaires
summary_tab$foldChange_kallisto <- 2^(summary_tab$log2FoldChange_kallisto)
summary_tab$foldChange_coco     <- 2^(summary_tab$log2FoldChange_coco)

# Ajout des pourcentages d’expression différentielle
# Formule : (foldChange - 1) * 100 (% d'augmentation par rapport au contrôle)
summary_tab$expression_diff_percent_kallisto <- NA
summary_tab$expression_diff_percent_coco     <- NA

not_na_kal <- !is.na(summary_tab$foldChange_kallisto)
summary_tab$expression_diff_percent_kallisto[not_na_kal] <-
  round((summary_tab$foldChange_kallisto[not_na_kal] - 1) * 100, 2)

not_na_coco <- !is.na(summary_tab$foldChange_coco)
summary_tab$expression_diff_percent_coco[not_na_coco] <-
  round((summary_tab$foldChange_coco[not_na_coco] - 1) * 100, 2)

# ------------------------------
# AJOUT DE LA COLONNE niveau_expression
# ------------------------------

summary_tab$niveau_expression <- NA_character_
for (i in seq_len(nrow(summary_tab))) {
  lfc_k <- summary_tab$log2FoldChange_kallisto[i]
  lfc_c <- summary_tab$log2FoldChange_coco[i]
  if (is.na(lfc_k) & is.na(lfc_c)) {
    summary_tab$niveau_expression[i] <- NA
  } else if (!is.na(lfc_k) & !is.na(lfc_c)) {
    if (lfc_k > 0 & lfc_c > 0) {
      summary_tab$niveau_expression[i] <- "surexpression"
    } else if (lfc_k < 0 & lfc_c < 0) {
      summary_tab$niveau_expression[i] <- "sousexpression"
    } else {
      summary_tab$niveau_expression[i] <- "mixte"
    }
  } else if (!is.na(lfc_k) & is.na(lfc_c)) {
    summary_tab$niveau_expression[i] <- ifelse(lfc_k > 0, "surexpression", ifelse(lfc_k < 0, "sousexpression", NA))
  } else if (is.na(lfc_k) & !is.na(lfc_c)) {
    summary_tab$niveau_expression[i] <- ifelse(lfc_c > 0, "surexpression", ifelse(lfc_c < 0, "sousexpression", NA))
  }
}

# Réorganisation des colonnes pour mettre niveau_expression juste après gene_symbol
summary_tab <- summary_tab %>%
  select(
    ensembl, gene_symbol, niveau_expression,
    baseMean_kallisto,
    log2FoldChange_kallisto, foldChange_kallisto, expression_diff_percent_kallisto, padj_kallisto,
    baseMean_coco,
    log2FoldChange_coco, foldChange_coco, expression_diff_percent_coco, padj_coco,
    everything()
  )

# Lecture du GTF pour extraire gene_biotype
gtf_lines <- readLines(gtf_file)
gtf_lines <- gtf_lines[!startsWith(gtf_lines, "#") & nzchar(gtf_lines)]
gtf_df <- str_split_fixed(gtf_lines, "\t", 9) %>%
  as.data.frame(stringsAsFactors = FALSE)
colnames(gtf_df) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

extract_attr <- function(attr, key) {
  m <- str_match(attr, paste0(key, ' "([^"]+)"'))
  ifelse(is.na(m[,2]), NA, m[,2])
}
gtf_df$gene_id      <- extract_attr(gtf_df$attribute, "gene_id")
gtf_df$gene_biotype <- extract_attr(gtf_df$attribute, "gene_biotype")

# Mapping unique (un seul gene_biotype par gene_id)
biotype_map <- gtf_df %>%
  filter(!is.na(gene_id), !is.na(gene_biotype)) %>%
  distinct(gene_id, gene_biotype)

# Merge sur ensembl
summary_tab <- summary_tab %>%
  left_join(biotype_map, by = c("ensembl" = "gene_id"))

# Ajout de la comparaison log2FC CoCo vs Kallisto en pourcentage
summary_tab$log2FC_CoCo_vs_Kallisto_percent <- NA
is_not_na <- !is.na(summary_tab$log2FoldChange_coco) & !is.na(summary_tab$log2FoldChange_kallisto) & summary_tab$log2FoldChange_kallisto != 0
summary_tab$log2FC_CoCo_vs_Kallisto_percent[is_not_na] <-
  100 * (summary_tab$log2FoldChange_coco[is_not_na] - summary_tab$log2FoldChange_kallisto[is_not_na]) / abs(summary_tab$log2FoldChange_kallisto[is_not_na])
summary_tab$log2FC_CoCo_vs_Kallisto_percent <- round(summary_tab$log2FC_CoCo_vs_Kallisto_percent, 2)

# Sauvegarde du tableau final
write.csv(
  summary_tab,
  file = out_file,
  row.names = FALSE,
  quote = FALSE
)
