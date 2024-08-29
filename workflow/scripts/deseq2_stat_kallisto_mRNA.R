# Chargement des bibliothèques nécessaires
library(ggplot2)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(ggrepel)
library(data.table)

# Lecture des données
résultats <- read.csv(snakemake@input[[1]], header = TRUE, sep = ",")
colnames(résultats)[1] <- "ensembl"

# Configuration de biomaRt pour le mapping des symboles des gènes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = résultats$ensembl, 
                   mart = ensembl)

# Jointure pour ajouter les symboles de gènes
résultats <- résultats %>%
  left_join(gene_info, by = c('ensembl' = 'ensembl_gene_id')) %>%
  rename(gene_symbol = hgnc_symbol) %>%
  mutate(gene_symbol = as.character(gene_symbol))

# Création de volcano plots avec texte agrandi
create_volcano_plot <- function(data, title, xlims = NULL, ylims = NULL) {
  colors <- c("non_signif"="black", "surexprimé"="red", "sousexprimé"="blue")
  plot <- ggplot(data, aes(x=log2FoldChange, y=-log10(pvalue), color = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "surexprimé",
    padj < 0.05 & log2FoldChange < -1 ~ "sousexprimé",
    TRUE ~ "non_signif"))) +
    geom_point(size = 2) +
    scale_color_manual(values = colors) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 p-value") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16),  # Augmenter la taille du titre
      axis.title.x = element_text(size = 14),  # Augmenter la taille du titre de l'axe X
      axis.title.y = element_text(size = 14),  # Augmenter la taille du titre de l'axe Y
      axis.text.x = element_text(size = 12),  # Augmenter la taille du texte de l'axe X
      axis.text.y = element_text(size = 12)   # Augmenter la taille du texte de l'axe Y
    )
  
  if (!is.null(xlims)) plot <- plot + xlim(xlims)
  if (!is.null(ylims)) plot <- plot + ylim(ylims)
  
  plot + geom_text_repel(
    data = subset(data, padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene_symbol),
    size = 4,  # Augmenter la taille des étiquettes de gènes
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
}

# Création et sauvegarde des plots
volcano_plot_total <- create_volcano_plot(résultats, "Volcano Plot")
volcano_plot_zoom <- create_volcano_plot(résultats, NULL, c(-10, 10), c(0, 15))
ggsave(snakemake@output[["volcano_plot_total"]], plot = volcano_plot_total)
ggsave(snakemake@output[["volcano_plot_zoom"]], plot = volcano_plot_zoom)

# Analyse d'enrichissement GO
genes_of_interest <- filter(résultats, padj <= 0.05 & !is.na(ensembl) & ensembl != "")
if (nrow(genes_of_interest) > 0) {
  ego_ensembl <- enrichGO(gene = genes_of_interest$ensembl, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                          ont = "ALL", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_ensembl) && 'result' %in% slotNames(ego_ensembl) && nrow(ego_ensembl@result) > 0) {
    plot_total <- dotplot(ego_ensembl, showCategory = 30) + ggtitle("Analyse d'enrichissement - Gene Ontology") + theme_minimal()
    ggsave(snakemake@output[["go_enrichment_total"]], plot = plot_total)
  } else {
    message("No significant GO terms found.")
  }
}

# Analyse d'enrichissement GO and log2FoldChange > 0 (SUREXPRIMÉ)
genes_of_interest <- filter(résultats, padj < 0.05 & log2FoldChange > 0 & !is.na(ensembl) & ensembl != "")
if (nrow(genes_of_interest) > 0) {
  ego_ensembl <- enrichGO(gene = genes_of_interest$ensembl, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                          ont = "ALL", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_ensembl) && 'result' %in% slotNames(ego_ensembl) && nrow(ego_ensembl@result) > 0) {
    plot_up <- dotplot(ego_ensembl, showCategory = 30) + theme_minimal()
    ggsave(snakemake@output[["go_enrichment_up"]], plot = plot_up)
  } else {
    message("No significant GO terms found.")
  }
}

# Analyse d'enrichissement GO and log2FoldChange < 0 (SOUSEXPRIMÉ)
genes_of_interest <- filter(résultats, padj < 0.05 & log2FoldChange < 0 & !is.na(ensembl) & ensembl != "")
if (nrow(genes_of_interest) > 0) {
  ego_ensembl <- enrichGO(gene = genes_of_interest$ensembl, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                          ont = "ALL", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_ensembl) && 'result' %in% slotNames(ego_ensembl) && nrow(ego_ensembl@result) > 0) {
    plot_down <- dotplot(ego_ensembl, showCategory = 25) + theme_minimal()
    ggsave(snakemake@output[["go_enrichment_down"]], plot = plot_down)
  } else {
    message("No significant GO terms found.")
  }
}

# Sauvegarde du fichier CSV
write.csv(résultats, snakemake@output[["deseq2_stat"]], row.names = FALSE)
