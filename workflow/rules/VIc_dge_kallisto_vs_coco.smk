rule deseq2_protein_coding_kallisto_vs_coco:
    """
    Comparaison des gènes différentiellement exprimés (padj <= 0.05) entre Kallisto et COCO.
    Produit un tableau: gene, gene_symbol, log2FoldChange_kallisto, padj_kallisto, log2FoldChange_coco, padj_coco.
    """
    input:
        coco_stat     = "results/dge/coco_protein_coding/{comp}_deseq2_stats.csv",
        kallisto_stat = "results/dge/kallisto_protein_coding/{comp}_deseq2_stats.csv",
        gtf           = "data/references/gtf/Homo_sapiens.GRCh38.110_snoRNAs_tRNAs_correct_annotation.gtf"
    output:
        summary = "results/dge/kallisto_vs_coco/deseq2_protein_coding_kallisto_vs_coco_{comp}.csv"
    log:
        "logs/dge/compare_deseq2_kallisto_vs_coco_{comp}.log"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/deseq2_protein_coding_kallisto_vs_coco.R"

