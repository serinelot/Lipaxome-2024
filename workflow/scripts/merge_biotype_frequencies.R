# Charger les packages nécessaires
library(dplyr)
library(tidyr)

# Importer les variables depuis Snakemake
input_files <- snakemake@input[["biotype_files"]]
output_frequencies <- snakemake@output[["merged_frequencies"]]
output_abundance <- snakemake@output[["merged_abundance"]]

# Fonction pour lire et sélectionner les colonnes `gene_biotype`, `frequency`, et `abundance`
read_and_select <- function(file, id) {
  data <- read.delim(file, sep="\t")
  frequency_data <- data %>%
    select(gene_biotype, frequency) %>%
    rename(!!id := frequency)  # Renommer la colonne frequency avec le nom du participant (id)
  abundance_data <- data %>%
    select(gene_biotype, abundance) %>%
    rename(!!id := abundance)  # Renommer la colonne abundance avec le nom du participant (id)
  list(frequency = frequency_data, abundance = abundance_data)
}

# Initialiser des dataframes vides pour stocker les résultats fusionnés
merged_frequencies <- NULL
merged_abundances <- NULL

# Fonction d'arrondi spécial pour les petits pourcentages
round_special <- function(x) {
  ifelse(x < 0.01 & x > 0, format(round(x, 5), nsmall = 5), format(round(x, 2), nsmall = 2))
}

# Parcourir chaque fichier d'entrée et fusionner les résultats
for (file in input_files) {
  id <- gsub("_biotype_proportions.tsv", "", basename(file))  # Extraire le `id` à partir du nom de fichier
  data <- read_and_select(file, id)
  
  # Fusionner les données de fréquence
  if (is.null(merged_frequencies)) {
    merged_frequencies <- data$frequency
  } else {
    merged_frequencies <- full_join(merged_frequencies, data$frequency, by = "gene_biotype")
  }
  
  # Fusionner les données d'abondance
  if (is.null(merged_abundances)) {
    merged_abundances <- data$abundance
  } else {
    merged_abundances <- full_join(merged_abundances, data$abundance, by = "gene_biotype")
  }
}

# Supprimer la ligne "Total" si elle existe
merged_frequencies <- merged_frequencies %>%
  filter(gene_biotype != "Total")
merged_abundances <- merged_abundances %>%
  filter(gene_biotype != "Total")

# Calculer la moyenne, l'écart-type et le coefficient de variation (CV) pour les fréquences
merged_frequencies <- merged_frequencies %>%
  rowwise() %>%
  mutate(
    moyenne = mean(c_across(-gene_biotype), na.rm = TRUE),  # Calculer la moyenne
    ecart_type = sd(c_across(-gene_biotype), na.rm = TRUE),  # Calculer l'écart-type
    CV = ifelse(moyenne != 0, (ecart_type / moyenne) * 100, NA)  # Calculer le CV
  ) %>%
  ungroup() %>%
  mutate(
    moyenne = round_special(moyenne * 100),
    ecart_type = round_special(ecart_type * 100),
    CV = round_special(CV),
    across(c(-gene_biotype, -moyenne, -ecart_type, -CV), ~round_special(. * 100))  # Convertir en pourcentage et arrondir
  ) %>%
  arrange(desc(as.numeric(moyenne)))  # Trier par la moyenne du plus grand au plus petit

# Calculer la moyenne, l'écart-type et le coefficient de variation (CV) pour les abondances
merged_abundances <- merged_abundances %>%
  rowwise() %>%
  mutate(
    moyenne = mean(c_across(-gene_biotype), na.rm = TRUE),  # Calculer la moyenne
    ecart_type = sd(c_across(-gene_biotype), na.rm = TRUE),  # Calculer l'écart-type
    CV = ifelse(moyenne != 0, (ecart_type / moyenne) * 100, NA)  # Calculer le CV
  ) %>%
  ungroup() %>%
  mutate(
    moyenne = round_special(moyenne * 100),
    ecart_type = round_special(ecart_type * 100),
    CV = round_special(CV),
    across(c(-gene_biotype, -moyenne, -ecart_type, -CV), ~round_special(. * 100))  # Convertir en pourcentage et arrondir
  ) %>%
  arrange(desc(as.numeric(moyenne)))  # Trier par la moyenne du plus grand au plus petit

# Gérer les NA en remplaçant par 0 si nécessaire
merged_frequencies <- merged_frequencies %>%
  mutate(across(c(moyenne, ecart_type, CV), ~replace_na(.x, 0)))
merged_abundances <- merged_abundances %>%
  mutate(across(c(moyenne, ecart_type, CV), ~replace_na(.x, 0)))

# Sauvegarder les fichiers fusionnés
write.table(merged_frequencies, file=output_frequencies, sep="\t", quote=FALSE, row.names=FALSE)
write.table(merged_abundances, file=output_abundance, sep="\t", quote=FALSE, row.names=FALSE)
