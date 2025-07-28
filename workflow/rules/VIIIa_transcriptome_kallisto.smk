rule transcriptome_deseq2_kallisto:
    input:
        quant = expand("results/quant/kallisto/{id}/abundance.h5", id=id_list),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
    output:
        results = directory("results/transcriptome/kallisto/deseq2"),
        out_files = expand("results/transcriptome/kallisto/deseq2/{comp}_DESeq2_transcripts.csv", comp=comparisons)
    params:
        filter_count_threshold = config["dge"]["filter_count_threshold"]
    log:
        "logs/dge/kallisto_transcriptome.log"
    message:
        "Transcriptome-wide DESeq2 (Kallisto, complet)"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_kallisto_transcriptome.R"


rule transcriptome_stats_kallisto:
    input:
        deseq2 = "results/transcriptome/kallisto/deseq2/{comp}_DESeq2_transcripts.csv",
        gtf = config["download"]["human_gtf"]
    output:
        stats = "results/transcriptome/kallisto/deseq2/{comp}_transcriptome_stats.csv"
    log:
        "logs/dge/kallisto_transcriptome_stats_simple_{comp}.log"
    message:
        "Ajout annotation et stats (transcript_name, biotype) — {wildcards.comp}"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_kallisto_transcriptome_stats_simple.R"


rule gene_deseq2_kallisto:
    """
    Gene-level DESeq2 (Kallisto, codant + non-codant) sans directory()
    Sommation des isoformes via tx2gene_all
    """
    input:
        quant       = expand(
                        "results/quant/kallisto/{id}/abundance.h5",
                        id = id_list
                      ),
        samples     = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        tx2gene     = rules.build_tx2gene_all.output.tx2gene_all
    output:
        out_files   = expand(
                        "results/transcriptome/kallisto/deseq2/gene/{comp}_DESeq2_genes.csv",
                        comp = comparisons
                      )
    params:
        filter_count_threshold = config["dge"]["filter_count_threshold"]
    log:
        "logs/dge/kallisto_gene.log"
    message:
        "Running gene-level DGE (Kallisto) with DESeq2"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_kallisto_gene.R"


rule gene_stats_kallisto:
    """
    Ajout d'annotations et stats (gene_name, gene_biotype) —
    gene-level DESeq2 Kallisto
    """
    input:
        deseq2 = "results/transcriptome/kallisto/deseq2/gene/{comp}_DESeq2_genes.csv",
        gtf    = config["download"]["human_gtf"]      # peut être un vecteur de GTF
    output:
        stats  = "results/transcriptome/kallisto/deseq2/gene/{comp}_gene_stats.csv"
    log:
        "logs/dge/kallisto_gene_stats_simple_{comp}.log"
    message:
        "Ajout annotation et stats (gene_name, gene_biotype) — {wildcards.comp}"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_kallisto_gene_stats_simple.R"
