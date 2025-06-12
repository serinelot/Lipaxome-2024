# ========== Librairies ==========
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ========== Arguments Snakemake ==========
input_file      <- snakemake@input[["stats"]]
output_total    <- snakemake@output[["enrich_total_csv"]]
output_up       <- snakemake@output[["enrich_up_csv"]]
output_down     <- snakemake@output[["enrich_down_csv"]]
plot_total      <- snakemake@output[["barplot_total"]]
plot_up         <- snakemake@output[["barplot_up"]]
plot_down       <- snakemake@output[["barplot_down"]]
organism        <- snakemake@params[["organism"]] # "human" ou "mouse"
top_n           <- snakemake@params[["top_n"]]

# ========== Lecture du fichier ==========
df <- read_csv(input_file, show_col_types = FALSE)

# ========== Sélection des listes ==========
genes_total <- unique(df$gene)
genes_up    <- unique(df$gene[df$log2FoldChange > 0])
genes_down  <- unique(df$gene[df$log2FoldChange < 0])

# ========== Conversion ENSEMBL -> ENTREZID ==========
convert_ensembl_to_entrez <- function(gene_list) {
  gene_df <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  entrez_ids[!is.na(entrez_ids)]
}

entrez_total <- convert_ensembl_to_entrez(genes_total)
entrez_up    <- convert_ensembl_to_entrez(genes_up)
entrez_down  <- convert_ensembl_to_entrez(genes_down)

# ========== Fonction enrichissement Reactome ==========
run_enrichPathway <- function(entrez_ids, organism) {
  if(length(entrez_ids) == 0) return(NULL)
  enrichPathway(
    gene         = entrez_ids,
    organism     = organism,  # "human" ou "mouse"
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable     = TRUE
  )
}

# ========== Fonction de plot ==========
plot_enrichment <- function(res, top_n, title, file) {
  if (!is.null(res) && nrow(res) > 0) {
    top_terms <- res %>%
      arrange(p.adjust) %>%
      head(top_n)
    p <- ggplot(top_terms, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = title,
           x = "Reactome pathway",
           y = "-log10(adjusted p-value)") +
      theme_minimal(base_size = 11) +
      theme(axis.text.y = element_text(size = 10)) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 50))
    ggsave(file, p, width = 10, height = 7)
  }
}

# ========== Enrichissement et export ==========
# TOTAL
epath_total <- run_enrichPathway(entrez_total, organism)
if (!is.null(epath_total)) {
  res_total <- as.data.frame(epath_total)
  write_csv(res_total, output_total)
  plot_enrichment(res_total, top_n, paste("Top", top_n, "Reactome pathways TOTAL"), plot_total)
} else {
  write_csv(tibble(), output_total)
}

# UP
epath_up <- run_enrichPathway(entrez_up, organism)
if (!is.null(epath_up)) {
  res_up <- as.data.frame(epath_up)
  write_csv(res_up, output_up)
  plot_enrichment(res_up, top_n, paste("Top", top_n, "Reactome pathways UP"), plot_up)
} else {
  write_csv(tibble(), output_up)
}

# DOWN
epath_down <- run_enrichPathway(entrez_down, organism)
if (!is.null(epath_down)) {
  res_down <- as.data.frame(epath_down)
  write_csv(res_down, output_down)
  plot_enrichment(res_down, top_n, paste("Top", top_n, "Reactome pathways DOWN"), plot_down)
} else {
  write_csv(tibble(), output_down)
}
