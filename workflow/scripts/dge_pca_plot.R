suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ----- Arguments Snakemake -----
vst_file    <- snakemake@input[["vst_mat"]]
design_file <- snakemake@input[["samples"]]
out_file    <- snakemake@output[["pca_plot"]]
ntop        <- as.integer(snakemake@params[["ntop"]])

# ----- Lecture -----
vst <- read_csv(vst_file, show_col_types = FALSE)
design <- read_tsv(design_file, show_col_types = FALSE)

# ----- Préparation de la matrice -----
rownames(vst) <- vst$gene
vst <- vst[, setdiff(colnames(vst), "gene")]
vst <- as.matrix(vst)

# ----- Sélection des gènes les plus variables -----
vars <- apply(vst, 1, var)
vst_top <- vst[order(vars, decreasing = TRUE)[1:min(ntop, nrow(vst))], , drop = FALSE]

# ----- PCA -----
pca <- prcomp(t(vst_top), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, design, by = c("sample" = "sample"))

# ----- Pourcentage de variance -----
percentVar <- round(100 * (pca$sdev)^2 / sum(pca$sdev^2), 1)

# ----- Plot -----
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1.2, size = 4) +
  labs(
    title = "PCA des échantillons (top gènes variables)",
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    color = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

ggsave(out_file, p, width = 8, height = 6)
