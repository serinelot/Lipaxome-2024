library(readr)
library(DESeq2)

# Variables venant de Snakemake
quant_dir <- snakemake@params[["quant_dir"]]
output_dir <- snakemake@output[["results"]]
design_file <- snakemake@input[["samples"]]
comparison_file <- snakemake@input[["comparisons"]]
transcript_gene_file <- snakemake@input[["gene_id"]]
filter_count_threshold <- as.numeric(snakemake@params[["filter_count_threshold"]])

# Créer le répertoire qui contiendra les résultats
dir.create(output_dir, showWarnings=FALSE)

# Charger les données
sampleTable <- read.table(design_file, header=TRUE, row.names="sample", check.names=FALSE)
conditions <- unique(sampleTable$condition)
samples <- rownames(sampleTable)
sampleTable$condition <- factor(sampleTable$condition, levels=conditions)

comparisons <- read.table(comparison_file, header=TRUE)
tx2gene <- read_tsv(transcript_gene_file, col_names=c('TXNAME', 'GENEID'))

# Lire les fichiers de quantification et vérifier la structure
files <- file.path(quant_dir, paste0(samples, ".tsv"))
quant_data <- lapply(files, function(x) {
    data <- read_tsv(x, col_types = cols())
    if (!"count" %in% colnames(data)) {
        stop("La colonne 'count' n'existe pas dans le fichier de quantification.")
    }
    return(data)
})

# Convertir les données en un format utilisable par DESeq2
counts <- sapply(quant_data, function(x) round(x$count))
colnames(counts) <- samples
rownames(counts) <- quant_data[[1]]$gene_id

# Préparer l'analyse avec DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~condition)
keep <- rowSums(counts(dds)) >= filter_count_threshold
dds <- dds[keep,]
dds$condition <- factor(dds$condition, conditions)
dds <- DESeq(dds)

# Boucler sur toutes les comparaisons pour calculer et enregistrer les résultats
for (row in 1:nrow(comparisons)) {
    condition1 <- comparisons[row, "cdn1"]
    condition2 <- comparisons[row, "cdn2"]
    exp <- sprintf("%s-%s", condition1, condition2)
    res <- results(dds, contrast = c("condition", condition2, condition1))
    res_df <- as.data.frame(res)
    fname <- paste(output_dir, paste(exp, "csv", sep='.'), sep='/')
    write.csv(res_df, file=fname, quote=FALSE)
}
