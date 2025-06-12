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
        ),
        vst_mat   = expand(
            "results/dge/coco_protein_coding/{comp}_vst_matrix.csv",
            comp = comparisons)  
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
    Analyse d'enrichissement GO (total, up, down) dans des dossiers séparés.
    """
    input:
        stats = "results/dge/coco_protein_coding/{comp}_deseq2_stats.csv"
    output:
        enrich_total_csv = "results/dge/coco_protein_coding/go/total/{comp}_go_enrichment_total.csv",
        enrich_up_csv    = "results/dge/coco_protein_coding/go/up/{comp}_go_enrichment_up.csv",
        enrich_down_csv  = "results/dge/coco_protein_coding/go/down/{comp}_go_enrichment_down.csv",
        barplot_total    = "results/dge/coco_protein_coding/go/total/{comp}_go_enrichment_barplot_total.png",
        barplot_up       = "results/dge/coco_protein_coding/go/up/{comp}_go_enrichment_barplot_up.png",
        barplot_down     = "results/dge/coco_protein_coding/go/down/{comp}_go_enrichment_barplot_down.png"
    params:
        org_db = "org.Hs.eg.db",
        ont = "BP",
        top_n = 15
    script:
        "../scripts/dge_go_enrichment.R"


rule dge_kegg_enrichment:
    """
    Enrichissement KEGG pour tous les DEGs, up et down, avec sorties dans dossiers dédiés.
    """
    input:
        stats = "results/dge/coco_protein_coding/{comp}_deseq2_stats.csv"
    output:
        enrich_total_csv = "results/dge/coco_protein_coding/kegg/total/{comp}_kegg_enrichment_total.csv",
        enrich_up_csv    = "results/dge/coco_protein_coding/kegg/up/{comp}_kegg_enrichment_up.csv",
        enrich_down_csv  = "results/dge/coco_protein_coding/kegg/down/{comp}_kegg_enrichment_down.csv",
        barplot_total    = "results/dge/coco_protein_coding/kegg/total/{comp}_kegg_enrichment_barplot_total.png",
        barplot_up       = "results/dge/coco_protein_coding/kegg/up/{comp}_kegg_enrichment_barplot_up.png",
        barplot_down     = "results/dge/coco_protein_coding/kegg/down/{comp}_kegg_enrichment_barplot_down.png"
    params:
        species = "hsa", # Code KEGG ("hsa" humain, "mmu" mouse, etc)
        top_n = 15
    script:
        "../scripts/dge_enrichKEGG.R"


rule dge_reactome_enrichment:
    """
    Enrichissement Reactome pour tous les DEGs, up et down, avec sorties dans dossiers dédiés.
    """
    input:
        stats = "results/dge/coco_protein_coding/{comp}_deseq2_stats.csv"
    output:
        enrich_total_csv = "results/dge/coco_protein_coding/reactome/total/{comp}_reactome_enrichment_total.csv",
        enrich_up_csv    = "results/dge/coco_protein_coding/reactome/up/{comp}_reactome_enrichment_up.csv",
        enrich_down_csv  = "results/dge/coco_protein_coding/reactome/down/{comp}_reactome_enrichment_down.csv",
        barplot_total    = "results/dge/coco_protein_coding/reactome/total/{comp}_reactome_enrichment_barplot_total.png",
        barplot_up       = "results/dge/coco_protein_coding/reactome/up/{comp}_reactome_enrichment_barplot_up.png",
        barplot_down     = "results/dge/coco_protein_coding/reactome/down/{comp}_reactome_enrichment_barplot_down.png"
    params:
        organism = "human", # "human" ou "mouse"
        top_n = 15
    conda:
        "../envs/reactome.yml"
    script:
        "../scripts/dge_reactome_enrichment.R"


rule dge_ma_plot:
    """
    Génère un MA plot (log2FC vs expression moyenne) à partir du CSV DESeq2.
    """
    input:
        stats = "results/dge/coco_protein_coding/{comp}_DESeq2_gene.csv"
    output:
        ma_plot = "results/dge/coco_protein_coding/{comp}_MAplot.png"
    params:
        padj_threshold = 0.05
    script:
        "../scripts/dge_ma_plot.R"


rule dge_heatmap_degs:
    """
    Heatmap des gènes différentiellement exprimés (padj < 0.05, top 50 par |log2FC|).
    """
    input:
        vst_mat = "results/dge/coco_protein_coding/{comp}_vst_matrix.csv",
        degs    = "results/dge/coco_protein_coding/{comp}_deseq2_stats.csv",
        samples = "data/design.tsv"
    output:
        heatmap = "results/dge/coco_protein_coding/heatmap/{comp}_heatmap_topDEG.png"
    params:
        padj_threshold = 0.05,
        top_n = 60
    conda:
        "../envs/pheatmap.yml"
    script:
        "../scripts/dge_heatmap_degs.R"


rule dge_pca_plot:
    """
    Génère un PCA plot des échantillons à partir de la matrice VST.
    """
    input:
        vst_mat = "results/dge/coco_protein_coding/{comp}_vst_matrix.csv",
        samples = "data/design.tsv"
    output:
        pca_plot = "results/dge/coco_protein_coding/{comp}_PCAplot.png"
    params:
        ntop = 500  # Nombre de gènes les plus variables à utiliser
    conda:
        "../envs/DESeq2.yml"   # Ou ton env DESeq2 si tu veux ggplot2
    script:
        "../scripts/dge_pca_plot.R"


rule dge_tsne_plot:
    """
    Génère un t-SNE plot des échantillons à partir de la matrice VST.
    """
    input:
        vst_mat = "results/dge/coco_protein_coding/{comp}_vst_matrix.csv",
        samples = "data/design.tsv"
    output:
        tsne_plot = "results/dge/coco_protein_coding/{comp}_tSNEplot.png"
    params:
        ntop = 500,           # Nombre de gènes les plus variables à utiliser
        perplexity = 2        # À adapter (<= n/3 pour n échantillons)
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_tsne_plot.R"




