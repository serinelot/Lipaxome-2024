rule filter_protein_coding_gtf_coco:
    """
    Extraire toutes les lignes du GTF corrigé par CoCo pour les gènes
    dont gene_biotype est "protein_coding", en conservant aussi les
    commentaires (#). Utile pour se concentrer strictement sur les
    gènes protéinocodants dans l’annotation CoCo.
    """
    input:
        gtf_corr = rules.coco_correct_annotation.output.gtf_corr
    output:
        pc_gtf_corr = "data/references/gtf/Homo_sapiens.GRCh38.110_snoRNAs_tRNAs_correct_annotation_protein_coding.gtf"
    log:
        "logs/coco/filter_pc_gtf_coco.log"
    message:
        "Filtering protein-coding genes from CoCo‑corrected GTF"
    shell:
        r"""
        mkdir -p $(dirname {output.pc_gtf_corr})
        grep -E '^#|gene_biotype "protein_coding"' {input.gtf_corr} \
            > {output.pc_gtf_corr} 2> {log}
        """

rule build_tx2gene_pc_coco:
    """
    Construction de la map transcript→gène pour les transcripts protéinocodants
    à partir du GTF corrigé par CoCo et filtré (only protein_coding).
    """
    input:
        gtf_corr_pc = rules.filter_protein_coding_gtf_coco.output.pc_gtf_corr
    output:
        tx2gene_pc_coco = "data/references/tx2gene_protein_coding_coco.tsv"
    log:
        "logs/dge/build_tx2gene_pc_coco.log"
    message:
        "Building transcript→gene map for CoCo‑corrected protein‑coding transcripts"
    run:
        import pandas as pd
        import re

        df = pd.read_csv(
            input.gtf_corr_pc,
            sep="\t",
            comment="#",
            header=None,
            usecols=[8],
            names=["attrs"]
        )

        mapping = {}
        for attr in df["attrs"]:
            if "gene_id" in attr and "transcript_id" in attr:
                gid = re.search(r'gene_id "([^"]+)"', attr).group(1)
                tid = re.search(r'transcript_id "([^"]+)"', attr).group(1)
                mapping[tid] = gid

        pd.DataFrame(
            mapping.items(),
            columns=["TXNAME","GENEID"]
        ).to_csv(
            output.tx2gene_pc_coco,
            sep="\t",
            index=False
        )


rule dge_coco_protein_coding:
    """
    Analyse différentielle d’expression (DESeq2) sur les gènes protéinocodants,
    à partir des comptages COCO cc (TSV).
    """
    input:
        quant       = expand("results/quant/coco/{id}.tsv", id=id_list),
        samples     = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        tx2gene     = rules.build_tx2gene_pc_coco.output.tx2gene_pc_coco
    output:
        results   = directory("results/dge/coco_protein_coding"),
        out_files = expand(
            "results/dge/coco_protein_coding/{comp}_DESeq2_gene.csv",
            comp = comparisons
        )
    params:
        coco_dir               = "results/quant/coco",
        filter_count_threshold = config["dge"]["filter_count_threshold"]
    log:
        "logs/dge/coco_dge_protein_coding.log"
    message:
        "Running gene‑level DGE (COCO cc, protein‑coding) with DESeq2"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_coco_protein_coding.R"


rule dge_coco_protein_coding_stats:
    """
    Statistiques post‑DESeq2 pour les gènes protéinocodants (COCO cc) :
      - calcul d’un tableau récapitulatif
      - tracé du volcano plot global
    """
    input:
        deseq2          = "results/dge/coco_protein_coding/{comp}_DESeq2_gene.csv"
    output:
        stat            = "results/dge/coco_protein_coding/{comp}_deseq2_stats.csv",
        volcano_total   = "results/dge/coco_protein_coding/{comp}_volcano_total.png"
    log:
        "logs/dge/coco_dge_protein_coding_stats_{comp}.log"
    message:
        "Calculating post-DESeq2 stats for COCO (protein‑coding) — comparison {wildcards.comp}"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_coco_protein_coding_stats.R"


rule dge_go_enrichment:
    """
    Analyse d'enrichissement GO sur les gènes différentiellement exprimés (padj <= 0.05).
    Génère un tableau de résultats et un barplot des top termes enrichis.
    """
    input:
        stats = rules.dge_coco_protein_coding_stats.output.stat
    output:
        enrich_csv = "results/dge/coco_protein_coding/go/{comp}_go_enrichment.csv",
        barplot    = "results/dge/coco_protein_coding/go/{comp}_go_enrichment_barplot.png"
    params:
        org_db = "org.Hs.eg.db",     # Adapter si nécessaire
        ont = "BP",                  # "BP" (biological process), "MF", "CC"
        top_n = 15                   # Nombre de termes à afficher dans le barplot
    script:
        "../scripts/dge_go_enrichment.R"

