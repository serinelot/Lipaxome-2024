suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(stringr)
  library(dplyr)
})

# Chemins via snakemake
deseq2_csv <- snakemake@input[["deseq2"]]
gtf_file   <- snakemake@input[["gtf_corr"]]
output_csv <- snakemake@output[["stat"]]

# 1. Lire le CSV DESeq2
dt <- fread(deseq2_csv)
colnames(dt)[1] <- "transcript_id" # sécurise le nom

# 2. Lire le GTF corrigé, extraire seulement les lignes "transcript"
gtf <- read_tsv(
  gtf_file, 
  comment = "#", 
  col_names = FALSE,
  show_col_types = FALSE
)

colnames(gtf)[9] <- "attributes"
colnames(gtf)[3] <- "feature"

# On ne garde que les "transcript"
gtf_trans <- gtf %>% filter(feature == "transcript")

# Extraire transcript_id et transcript_name des attributs
extract_attr <- function(attr, key) {
  stringr::str_match(attr, paste0(key, ' "([^"]+)"'))[,2]
}

gtf_table <- tibble(
  transcript_id = extract_attr(gtf_trans$attributes, "transcript_id"),
  transcript_name = extract_attr(gtf_trans$attributes, "transcript_name")
) %>% distinct(transcript_id, .keep_all = TRUE)

# 3. Fusion avec le CSV DESeq2 (ajoute transcript_name mais conserve toutes les colonnes)
dt2 <- left_join(dt, gtf_table, by = "transcript_id") %>% as.data.table()

# 4. Réorganiser pour placer transcript_name juste après transcript_id, sans perdre les autres colonnes
first_col <- "transcript_id"
insert_after <- "transcript_id"
all_cols <- colnames(dt2)
all_cols <- all_cols[all_cols != "transcript_name"]
insert_pos <- match(insert_after, all_cols)
final_cols <- append(all_cols, "transcript_name", after = insert_pos)
dt2 <- dt2[, ..final_cols]

# 5. Écrire la sortie
fwrite(dt2, file = output_csv, sep = ",", na = "NA", quote = FALSE)

cat("Terminé : transcript_name ajouté, colonnes conservées !\n")
