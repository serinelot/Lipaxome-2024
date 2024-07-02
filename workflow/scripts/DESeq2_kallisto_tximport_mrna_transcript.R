library(readr)
library(tximport)
library(DESeq2)

# Variables coming from the snakemake
kallisto_dir <- snakemake@params[["kallisto_dir"]]
output_dir <- snakemake@output[["results"]]
design_file <- snakemake@input[["samples"]]
comparison_file <- snakemake@input[["comparisons"]]
transcript_gene_file <- snakemake@input[["gene_id"]]
filter_count_threshold <- as.numeric(snakemake@params[["filter_count_threshold"]]) 

# Create the directory that will contain the results
dir.create(output_dir, showWarnings=FALSE)

# Loading data
sampleTable <- read.table(design_file, header=TRUE, row.names="sample", check.names=FALSE)
conditions <- unique(sampleTable$condition)
samples <- rownames(sampleTable)
sampleTable$condition <- factor(sampleTable$condition, levels=conditions)

comparisons <- read.table(comparison_file, header=TRUE)
tx2gene <- read_tsv(transcript_gene_file, col_names=c('TXNAME', 'GENEID'))

files <- file.path(kallisto_dir, samples, "abundance.h5")
names(files) <- samples
txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE)

# DESeq2 analysis setup
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= filter_count_threshold
dds <- dds[keep,]
dds$condition <- factor(dds$condition, conditions)
dds <- DESeq(dds)

# Looping over all comparisons
for (row in 1:nrow(comparisons)) {
    condition1 <- comparisons[row, "cdn1"]
    condition2 <- comparisons[row, "cdn2"]
    exp <- sprintf("%s-%s", condition1, condition2)
    res <- results(dds, contrast = c("condition", condition2, condition1))
    res_df <- as.data.frame(res)
    fname <- paste(output_dir, paste(exp, "csv", sep='.'), sep='/')
    write.csv(res_df, file=fname, quote=FALSE)
}
