# Charger les packages nécessaires
library(dplyr)

# Importer les variables depuis Snakemake
gtf_file <- snakemake@input[["gtf"]]
counts_file <- snakemake@input[["counts"]]
output_file <- snakemake@output[["counts_with_biotypes"]]

# Lire le fichier GTF
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, comment.char="#", stringsAsFactors=FALSE)

# Lire le fichier de comptage
counts_data <- read.delim(counts_file, sep="\t")

# Filtrer pour les lignes qui contiennent des transcrits seulement
gtf_data <- gtf_data %>%
  filter(V3 == "transcript")

# Décomposer la colonne des attributs en une liste
attributes_df <- strsplit(gtf_data$V9, ";")  # Séparer chaque attribut par ";"

# Extraire les valeurs pour 'gene_id', 'gene_name', et 'gene_biotype'
extracted_attributes <- data.frame(
  gene_id = sapply(attributes_df, function(x) {
    gene_id_entry <- grep("gene_id", x, value = TRUE)
    if (length(gene_id_entry) > 0) {
      return(gsub("gene_id \"|\"", "", gene_id_entry))
    } else {
      return(NA)
    }
  }),
  gene_name = sapply(attributes_df, function(x) {
    gene_name_entry <- grep("gene_name", x, value = TRUE)
    if (length(gene_name_entry) > 0) {
      return(gsub("gene_name \"|\"", "", gene_name_entry))
    } else {
      return(NA)
    }
  }),
  gene_biotype = sapply(attributes_df, function(x) {
    gene_biotype_entry <- grep("gene_biotype", x, value = TRUE)
    if (length(gene_biotype_entry) > 0) {
      return(gsub("gene_biotype \"|\"", "", gene_biotype_entry))
    } else {
      return(NA)
    }
  })
)

# Joindre les attributs extraits au gtf_data
gtf_data <- cbind(gtf_data, extracted_attributes)

# Garder uniquement les colonnes nécessaires et réorganiser
gtf_data <- gtf_data %>%
  select(gene_id, gene_name, gene_biotype)

# Nettoyer les colonnes
gtf_data <- gtf_data %>%
  mutate(
    gene_id = gsub("gene_id ", "", gene_id),
    gene_name = gsub("gene_name ", "", gene_name),
    gene_biotype = gsub("gene_biotype ", "", gene_biotype)
  )

# Agréger les données GTF par 'gene_id'
gtf_data_aggregated <- gtf_data %>%
  group_by(gene_id) %>%
  summarize(
    gene_name = first(gene_name),
    gene_biotype = first(gene_biotype)
  )

# Agréger les données de comptage par 'gene_id'
counts_data_aggregated <- counts_data %>%
  group_by(gene_id) %>%
  summarize(
    count = sum(count),
    cpm = sum(cpm),
    tpm = sum(tpm)
  )

# Fusionner les données agrégées
merged_data <- counts_data_aggregated %>%
  left_join(gtf_data_aggregated, by = "gene_id")

# Écrire le fichier de sortie
write.table(merged_data, file=output_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
