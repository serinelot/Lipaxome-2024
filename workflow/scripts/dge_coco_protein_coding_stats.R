#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(data.table)
})

dge_csv <- snakemake@input[["deseq2"]]

résultats <- read.csv(
  dge_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"gene" %in% colnames(résultats)) stop("Colonne 'gene' absente du CSV!")

résultats <- résultats[, c("gene", setdiff(names(résultats), "gene"))]

# --- Annoter SYMBOL ---
gene_info <- AnnotationDbi::select(
  x       = org.Hs.eg.db,
  keys    = résultats$gene,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

colnames(gene_info)[colnames(gene_info) == "ENSEMBL"] <- "gene"
colnames(gene_info)[colnames(gene_info) == "SYMBOL"]  <- "gene_symbol"
gene_info <- gene_info[, c("gene", "gene_symbol")]

résultats <- résultats %>%
  left_join(gene_info, by = "gene") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol), gene, gene_symbol)
  )

# ----- Filtrer padj <= 0.05 -----
résultats <- résultats %>% filter(!is.na(padj) & padj <= 0.05)

# --- Volcano global seulement ---
create_volcano_plot <- function(df, title = "") {
  df2 <- df %>%
    mutate(status = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "overexpressed",
      padj < 0.05 & log2FoldChange < -1 ~ "underexpressed",
      TRUE                              ~ "not_signif"
    ))
  p <- ggplot(df2, aes(
        x = log2FoldChange,
        y = -log10(pvalue),
        color = status
      )) +
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
  p
}

# Volcano total avec les gènes filtrés
volcano_total <- create_volcano_plot(résultats, "Volcano Plot – Global")
ggsave(
  filename = snakemake@output[["volcano_total"]],
  plot     = volcano_total,
  width    = 8, height = 6, dpi = 300
)

# Écrire le tableau filtré
fwrite(
  résultats,
  file      = snakemake@output[["stat"]],
  sep       = ",",
  quote     = FALSE,
  row.names = FALSE
)
