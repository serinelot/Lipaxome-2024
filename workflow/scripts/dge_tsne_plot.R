suppressPackageStartupMessages({
  library(Rtsne)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ----- Arguments Snakemake -----
vst_file    <- snakemake@input[["vst_mat"]]
design_file <- snakemake@input[["samples"]]
out_file    <- snakemake@output[["tsne_plot"]]
ntop        <- as.integer(snakemake@params[["ntop"]])
perplexity  <- as.numeric(snakemake@params[["perplexity"]])

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

# ----- t-SNE -----
set.seed(42)
tsne <- Rtsne(t(vst_top), perplexity = perplexity, verbose = TRUE, max_iter = 1000)

tsne_df <- as.data.frame(tsne$Y)
tsne_df$sample <- rownames(vst_top)
tsne_df$sample <- colnames(vst_top)  # Correction: les colonnes sont les samples
tsne_df <- left_join(tsne_df, design, by = c("sample" = "sample"))

# ----- Plot -----
p <- ggplot(tsne_df, aes(x = V1, y = V2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1.2, size = 4) +
  labs(
    title = "t-SNE des échantillons (top gènes variables)",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

ggsave(out_file, p, width = 8, height = 6)
