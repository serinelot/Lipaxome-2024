rule transcriptome_deseq2_coco:
    """
    Transcriptome‑wide DESeq2 (codants + non‑codants)
    à partir des TSV COCO cc.
    """
    input:
        quant       = expand("results/quant/coco/{id}.tsv", id=id_list),
        samples     = "data/design.tsv",
        comparisons = "data/comparisons.tsv"
    output:
        results   = directory("results/transcriptome/coco/deseq2"),
        out_files = expand(
            "results/transcriptome/coco/deseq2/{comp}_DESeq2_transcripts.csv",
            comp = comparisons
        )
    params:
        coco_dir               = "results/quant/coco",
        filter_count_threshold = config["dge"]["filter_count_threshold"]
    log:
        "logs/dge/coco_transcriptome.log"
    message:
        "Running transcriptome‑wide DGE (COCO cc) with DESeq2"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_coco_transcriptome.R"


rule transcriptome_stats_coco:
    """
    Statistiques post‑DESeq2 pour le transcriptome complet (COCO cc) :
      - tableau récapitulatif
      - volcano plots (total et zoom)
      - enrichissement GO (global, up, down)
    """
    input:
        deseq2       = "results/transcriptome/coco/deseq2/{comp}_DESeq2_transcripts.csv",
        gtf    = config["download"]["human_gtf"]
    output:
        stat         = "results/transcriptome/coco/deseq2/{comp}_transcriptome_stats.csv",

    log:
        "logs/dge/coco_transcriptome_stats_{comp}.log"
    message:
        "Calculating post‑DESeq2 stats for transcriptome COCO — comparison {wildcards.comp}"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_coco_transcriptome_stats.R"
