#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(biomaRt)
})

# Fichiers en entrée/sortie
summary_file <- snakemake@input[["summary"]]
out_file     <- snakemake@output[["enriched"]]

# Lecture du tableau
tab <- read.csv(summary_file, stringsAsFactors = FALSE)

# Nom de colonne contenant les ENSG
ens_col <- c("ensembl", "gene", "GeneID", "ensembl_gene_id")
col_used <- intersect(ens_col, colnames(tab))[1]
if (is.na(col_used)) stop("Aucune colonne Ensembl trouvée dans le tableau !")

ensembl_ids <- unique(tab[[col_used]][!is.na(tab[[col_used]])])

# Initialisation biomart (catch retry)
for (try_i in 1:3) {
  mart <- tryCatch(
    useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
    error = function(e) NULL
  )
  if (!is.null(mart)) break
  Sys.sleep(2)
}
if (is.null(mart)) stop("Impossible de contacter Ensembl biomart.")

# Récupérer annotation GO BP et MF
go_annot <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = mart
)

go_bp <- go_annot %>%
  filter(namespace_1003 == "biological_process" & !is.na(name_1006)) %>%
  group_by(ensembl_gene_id) %>%
  summarize(processus_biologiques = paste(unique(name_1006), collapse = "; "), .groups = "drop")

go_mf <- go_annot %>%
  filter(namespace_1003 == "molecular_function" & !is.na(name_1006)) %>%
  group_by(ensembl_gene_id) %>%
  summarize(fonctions_moleculaires = paste(unique(name_1006), collapse = "; "), .groups = "drop")

# Récupérer les voies Reactome
pw_annot <- getBM(
  attributes = c("ensembl_gene_id", "reactome"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = mart
)

reactome <- pw_annot %>%
  filter(!is.na(reactome)) %>%
  group_by(ensembl_gene_id) %>%
  summarize(voies_metaboliques = paste(unique(reactome), collapse = "; "), .groups = "drop")

# Fusion sur le tableau principal
tab2 <- tab %>%
  left_join(go_bp,      by = setNames("ensembl_gene_id", col_used)) %>%
  left_join(go_mf,      by = setNames("ensembl_gene_id", col_used)) %>%
  left_join(reactome,   by = setNames("ensembl_gene_id", col_used))

# Sauvegarde du tableau annoté
write.csv(
  tab2,
  file      = out_file,
  row.names = FALSE,
  quote     = FALSE
)
