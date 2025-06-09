# ========== Librairies ==========
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db) # À adapter selon l'espèce (ex : org.Mm.eg.db pour mouse)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ========== Arguments Snakemake ==========
input_file  <- snakemake@input[["stats"]]
output_csv  <- snakemake@output[["enrich_csv"]]
output_plot <- snakemake@output[["barplot"]]
org_db      <- snakemake@params[["org_db"]]
ont         <- snakemake@params[["ont"]]
top_n       <- snakemake@params[["top_n"]]

# ========== Lecture du fichier ==========
df <- read_csv(input_file)

# Adapter la colonne selon ton identifiant de gène
# On suppose ici "gene_id" correspond à des symboles (ex: "BRCA1"), sinon adapter "keyType"
gene_list <- unique(df$gene_id)

# ========== Enrichissement GO ==========
ego <- enrichGO(
  gene          = gene_list,
  OrgDb         = get(org_db),
  keyType       = "SYMBOL",   # Adapter à "ENSEMBL", "ENTREZID", etc selon tes IDs
  ont           = ont,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# ========== Export des résultats ==========
res <- as.data.frame(ego)
write_csv(res, output_csv)

# ========== Barplot ==========
if (nrow(res) > 0) {
  top_terms <- res %>%
    arrange(p.adjust) %>%
    head(top_n)
  
  p <- ggplot(top_terms, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top", top_n, "GO terms (", ont, ")"),
         x = "GO term",
         y = "-log10(adjusted p-value)") +
    theme_minimal()
  
  ggsave(output_plot, p, width = 8, height = 5)
}
