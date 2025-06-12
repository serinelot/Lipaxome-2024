suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(readr)
})

# Entrées
vst_file    <- snakemake@input[["vst_mat"]]
degs_file   <- snakemake@input[["degs"]]
design_file <- snakemake@input[["samples"]]
heatmap_file <- snakemake@output[["heatmap"]]
top_n       <- as.integer(snakemake@params[["top_n"]])

# Lecture et nettoyage
vst_mat <- read_csv(vst_file, show_col_types = FALSE)
degs <- read_csv(degs_file, show_col_types = FALSE)
design <- read_tsv(design_file, show_col_types = FALSE)

# Correction IDs
degs$gene <- sub("\\..*", "", degs$gene)
vst_mat$gene <- sub("\\..*", "", vst_mat$gene)
degs$gene <- trimws(degs$gene)
vst_mat$gene <- trimws(vst_mat$gene)

# Sélection top N DEGs
degs_sel <- degs %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n)

genes_heat <- degs_sel$gene

# Extraction, fusion, conversion en matrice, rownames = symbols
mat_heat <- vst_mat %>%
  filter(gene %in% genes_heat) %>%
  arrange(match(gene, genes_heat))

symbols <- degs_sel$gene_symbol[match(mat_heat$gene, degs_sel$gene)]
symbols[is.na(symbols) | symbols == ""] <- mat_heat$gene[is.na(symbols) | symbols == ""]

# ==> On fait une matrice avec rownames = symbols, plus de colonne gene !
mat_heat2 <- as.matrix(mat_heat %>% select(-gene))
rownames(mat_heat2) <- make.unique(symbols)

# Annotation (même principe que splicing)
annotation_col <- design %>%
  select(sample, condition) %>%
  as.data.frame()
rownames(annotation_col) <- annotation_col$sample
annotation_col <- annotation_col[colnames(mat_heat2), , drop = FALSE]

# LARGEUR spéciale pour la lisibilité
largeur_img <- 15
hauteur_img <- max(10, nrow(mat_heat2) * 0.6)

# 🎨 HEATMAP avec gene symbols, style identique à splicing
pheatmap(
  mat_heat2,
  annotation_col = annotation_col["condition", drop=FALSE],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 12,
  fontsize_row = 10,
  fontsize_col = 11,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  width = largeur_img,
  height = hauteur_img,
  filename = heatmap_file,
  main = "Heatmap des top DEGs (gene symbols visibles)"
)
