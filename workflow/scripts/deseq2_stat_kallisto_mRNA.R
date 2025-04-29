#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 0) Packages et options
options(stringsAsFactors = FALSE)
library(data.table)       # fread
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1) Lecture des résultats DESeq2
résultats <- fread(
  snakemake@input[[1]],
  sep    = ",",
  header = TRUE
) %>%
  rename(ensembl = 1)

# 2) Annotation locale via org.Hs.eg.db
gene_info <- AnnotationDbi::select(
  x       = org.Hs.eg.db,
  keys    = résultats$ensembl,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
) %>%
  rename(
    ensembl     = ENSEMBL,
    gene_symbol = SYMBOL
  )

# 3) Jointure et fallback à l’ID Ensembl
résultats <- résultats %>%
  left_join(gene_info, by = "ensembl") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol), ensembl, gene_symbol)
  )

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4) Fonction pour volcano plot
create_volcano_plot <- function(df, title = "", xlims = NULL, ylims = NULL) {
  df2 <- df %>%
    mutate(status = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "surexprimé",
      padj < 0.05 & log2FoldChange < -1 ~ "sousexprimé",
      TRUE                              ~ "non_signif"
    ))
  
  p <- ggplot(df2, aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      color = status
    )) +
    geom_point(size = 2) +
    scale_color_manual(values = c(
      non_signif = "black",
      surexprimé = "red",
      sousexprimé = "blue"
    )) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(p-value)") +
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

# 5) Création et sauvegarde des volcano plots
volcano_plot_total <- create_volcano_plot(résultats, "Volcano Plot - global")
volcano_plot_zoom  <- create_volcano_plot(
  résultats,
  "Volcano Plot - zoom",
  xlims = c(-10,10),
  ylims = c(0,15)
)

ggsave(snakemake@output[["volcano_plot_total"]],
       plot   = volcano_plot_total,
       width  = 8, height = 6, dpi = 300)
ggsave(snakemake@output[["volcano_plot_zoom"]],
       plot   = volcano_plot_zoom,
       width  = 8, height = 6, dpi = 300)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 6) Fonction générique pour enrichment GO
run_go <- function(filter_expr, output_name, showCategory = 30) {
  genes <- résultats %>%
    filter(eval(filter_expr)) %>%
    pull(ensembl) %>%
    na.omit()
  
  if (length(genes) == 0) {
    message("Aucun gène pour ", output_name)
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
         ggtitle(output_name) +
         theme_minimal()
    ggsave(
      filename = snakemake@output[[output_name]],
      plot     = p,
      width    = 8, height = 6, dpi = 300
    )
  } else {
    message("Pas de termes GO significatifs pour ", output_name)
  }
}

# Enrichissements GO global, up, down
run_go(TRUE,                         "go_enrichment_total", showCategory = 30)
run_go(quote(padj < 0.05 & log2FoldChange >  0),
       "go_enrichment_up",    showCategory = 30)
run_go(quote(padj < 0.05 & log2FoldChange <  0),
       "go_enrichment_down",  showCategory = 25)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 7) Sauvegarde du tableau annoté en CSV
fwrite(
  résultats,
  file      = snakemake@output[["deseq2_stat"]],
  sep       = ",",
  quote     = FALSE,
  row.names = FALSE
)
