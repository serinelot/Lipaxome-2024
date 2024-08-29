# Charger les packages nécessaires
library(dplyr)

# Importer les variables depuis Snakemake
input_file <- snakemake@input[["counts_with_biotypes"]]
output_file <- snakemake@output[["biotype_proportions"]]

# Charger le fichier TSV
data <- read.delim(input_file, sep="\t")

# Filtrer les lignes où count, cpm, et tpm ne sont pas tous égaux à 0
filtered_data <- data %>%
  filter(!(count == 0 & cpm == 0 & tpm == 0))

# Calculer la somme des counts pour chaque type de biotype
biotype_counts <- filtered_data %>%
  group_by(gene_biotype) %>%
  summarize(
    gene_count = n(),  # Nombre de gènes par biotype
    total_count = sum(count)  # Somme des counts par biotype
  )

# Calculer les proportions en pourcentage pour les gènes et les counts
total_gene_counts <- sum(biotype_counts$gene_count)  # Total des gènes
total_counts_sum <- sum(biotype_counts$total_count)  # Total des counts

biotype_counts <- biotype_counts %>%
  mutate(
    perc_frequency = paste0(round((gene_count / total_gene_counts) * 100, 2), "%"),  # Proportion formatée en pourcentage (frequency)
    perc_abundance = paste0(round((total_count / total_counts_sum) * 100, 2), "%"),  # Proportion formatée en pourcentage (abundance)
    frequency = (gene_count / total_gene_counts),  # Proportion sur une échelle de 0 à 1 (nombre de gènes)
    abundance = (total_count / total_counts_sum)  # Proportion sur une échelle de 0 à 1 (counts)
  )

# Réorganiser les colonnes dans l'ordre souhaité
biotype_counts <- biotype_counts %>%
  select(gene_biotype, gene_count, total_count, perc_frequency, perc_abundance, frequency, abundance)

# Ajouter une ligne avec le total
biotype_counts <- biotype_counts %>%
  add_row(
    gene_biotype = "Total",
    gene_count = total_gene_counts,
    total_count = total_counts_sum,
    perc_frequency = "100%",
    perc_abundance = "100%",
    frequency = 1,
    abundance = 1
  )

# Sauvegarder les résultats dans un fichier TSV
write.table(biotype_counts, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
