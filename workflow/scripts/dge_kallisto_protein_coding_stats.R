#!/usr/bin/env Rscript
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 0) Chargement des packages et options
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)       # fread, fwrite, setnames
  library(dplyr)            # left_join, mutate, filter
  library(ggplot2)          # ggplot, ggsave
  library(ggrepel)          # geom_text_repel
  library(DESeq2)           # DESeq2 core
  library(clusterProfiler)  # enrichGO, dotplot
  library(org.Hs.eg.db)     # OrgDb pour humain
  library(AnnotationDbi)    # select
  library(patchwork)        # combiner plots
})

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1) Lecture du résultat DESeq2 brut (un seul fichier .csv)
input_csv <- snakemake@input[["deseq2"]]
résultats <- fread(input_csv, sep = ",", header = TRUE)

# Renommer la première colonne (Ensembl IDs) en "gene"
setnames(résultats,
         old = names(résultats)[1],
         new = "gene")

# 2) Annotation des IDs Ensembl en symboles HGNC
gene_info <- AnnotationDbi::select(
  x       = org.Hs.eg.db,
  keys    = résultats$gene,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
) %>%
  rename(gene       = ENSEMBL,
         gene_symbol = SYMBOL)

résultats <- résultats %>%
  left_join(gene_info, by = "gene") %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), gene, gene_symbol))

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3) Fonction pour créer un volcano plot
create_volcano_plot <- function(df, title = "", xlims = NULL, ylims = NULL) {
  df2 <- df %>%
    mutate(status = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "overexpressed",
      padj < 0.05 & log2FoldChange < -1 ~ "underexpressed",
      TRUE                              ~ "not_signif"
    ))
  p <- ggplot(df2, aes(x = log2FoldChange, y = -log10(pvalue), color = status)) +
    geom_point(size = 2) +
    scale_color_manual(values = c(
      not_signif    = "black",
      overexpressed = "red",
      underexpressed= "blue"
    )) +
    labs(title = title,
         x     = "Log2 Fold Change",
         y     = "-Log10(p-value)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 16, hjust = 0.5),
      axis.title      = element_text(size = 14),
      axis.text       = element_text(size = 12)
    ) +
    geom_text_repel(
      data = filter(df2, padj < 0.05 & abs(log2FoldChange) > 1),
      aes(label = gene_symbol),
      size         = 4,
      box.padding  = unit(0.35, "lines"),
      point.padding= unit(0.5,  "lines")
    )
  if (!is.null(xlims)) p <- p + xlim(xlims)
  if (!is.null(ylims)) p <- p + ylim(ylims)
  return(p)
}

# 4) Génération et sauvegarde des volcano plots
volcano_total <- create_volcano_plot(résultats, "Volcano Plot – Global")
volcano_zoom  <- create_volcano_plot(
  résultats,
  "Volcano Plot – Zoom",
  xlims = c(-5, 5),
  ylims = c(0, 15)
)

ggsave(
  filename = snakemake@output[["volcano_total"]],
  plot     = volcano_total,
  width    = 8, height = 6, dpi = 300
)
ggsave(
  filename = snakemake@output[["volcano_zoom"]],
  plot     = volcano_zoom,
  width    = 8, height = 6, dpi = 300
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 5) Fonction générique pour enrichissement GO
run_go <- function(filter_expr, output_key, showCategory = 30) {
  genes <- résultats %>%
    filter(eval(filter_expr)) %>%
    pull(gene) %>%
    na.omit()
  if (length(genes) == 0) {
    message("Aucun gène pour ", output_key)
    return(NULL)
  }
  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = "ALL",
    pAdjustMethod = "fdr",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (!is.null(ego) && nrow(ego@result) > 0) {
    p <- dotplot(ego, showCategory = showCategory) +
         ggtitle(output_key) +
         theme_minimal()
    ggsave(
      filename = snakemake@output[[output_key]],
      plot     = p,
      width    = 8, height = 6, dpi = 300
    )
  } else {
    message("Pas de termes GO significatifs pour ", output_key)
  }
}

# 6) Enrichissements GO global / up / down
run_go(TRUE,                                   "go_total", showCategory = 30)
run_go(quote(padj < 0.05 & log2FoldChange >  0), "go_enrich_up",   showCategory = 25)
run_go(quote(padj < 0.05 & log2FoldChange <  0), "go_enrich_down", showCategory = 25)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 7) Sauvegarde du tableau annoté complet (statistiques + symboles)
fwrite(
  résultats,
  file      = snakemake@output[["stat"]],
  sep       = ",",
  quote     = FALSE,
  row.names = FALSE
)
