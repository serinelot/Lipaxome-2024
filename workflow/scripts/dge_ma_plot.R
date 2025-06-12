# ========== Librairies ==========
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ========== Arguments Snakemake ==========
input_file <- snakemake@input[["stats"]]
output_file <- snakemake@output[["ma_plot"]]
padj_thr <- as.numeric(snakemake@params[["padj_threshold"]])

# ========== Lecture du CSV ==========
df <- read_csv(input_file, show_col_types = FALSE)

# ========== Vérification ==========
if (!all(c("baseMean", "log2FoldChange", "padj") %in% colnames(df))) {
  stop("Le fichier d'entrée doit contenir les colonnes 'baseMean', 'log2FoldChange', 'padj'.")
}

# ========== Préparation ==========
df <- df %>%
  mutate(
    log10baseMean = log10(baseMean + 1), # +1 pour éviter log(0)
    Significant = ifelse(!is.na(padj) & padj <= padj_thr, "DEG", "non DEG")
  )

# ========== Plot ==========
p <- ggplot(df, aes(x = log10baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("DEG" = "#d73027", "non DEG" = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = "MA plot",
    x = "log10(baseMean + 1)",
    y = "log2(Fold Change)",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave(output_file, p, width = 7, height = 5)
