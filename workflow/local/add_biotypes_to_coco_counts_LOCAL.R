library(dplyr)

# Lire le fichier GTF en ignorant les lignes de métadonnées
gtf_file <- "C:/Users/bens3307/Documents/Globus/Homo_sapiens.GRCh38.110_snoRNAs_tRNAs_correct_annotation.gtf"
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, comment.char="#", stringsAsFactors=FALSE)

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

# Agréger les données par 'gene_id' pour éviter les duplications dans le GTF
gtf_data_aggregated <- gtf_data %>%
  group_by(gene_id) %>%
  summarize(
    gene_name = first(gene_name),
    gene_biotype = first(gene_biotype)
  )

# Charger votre fichier de données (en TSV)
data_file <- "C:/Users/bens3307/Documents/Globus/LipC03.tsv"
data <- read.delim(data_file, sep = "\t")

# Agréger les données de comptage par 'gene_id'
data_aggregated <- data %>%
  group_by(gene_id) %>%
  summarize(
    count = sum(count),
    cpm = sum(cpm),
    tpm = sum(tpm)
  )

# Joindre les données GTF agrégées avec les données de comptage agrégées
merged_data <- data_aggregated %>%
  left_join(gtf_data_aggregated, by = "gene_id")

# Afficher les premières lignes du dataframe fusionné
head(merged_data)

# Sauvegarder le résultat dans un nouveau fichier au format TSV
write.table(merged_data, "output_with_gene_name_and_biotype.tsv", sep="\t", quote=FALSE, row.names=FALSE)