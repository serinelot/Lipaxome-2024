#### FONCTIONNE POUR FXS-Ccontrol SEULEMENT ##########

rule deseq2:
    input:
        quant = expand("results/dge/kallisto/{id}/abundance.tsv", id = id_list),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = "data/references/tx2gene.tsv"
    output:
        results = directory("results/dge/deseq2"),
        deseq2_output = "results/dge/deseq2/FXS-Control.csv"
    params:
        kallisto_dir = "results/dge/kallisto"
    log:
        "logs/deseq2.log"
    message:
        "Perform differential expression analysis for various conditions."
    script:
        "../scripts/DESeq2_kallisto_tximport.R"
