rule filter_protein_coding_gtf:
    """
    Extraire toutes les lignes du GTF pour les gènes dont gene_biotype est "protein_coding",
    en conservant aussi les commentaires (#).
    Utile pour se concentrer strictement sur les gènes protéinocodants.
    """
    input:
        gtf = config["download"]["human_gtf"]
    output:
        pc_gtf = "data/references/gtf/Homo_sapiens.GRCh38.110_protein_coding.gtf"
    log:
        "logs/kallisto/filter_pc_gtf.log"
    message:
        "Filtering GTF for protein-coding genes only"
    shell:
        r"""
        mkdir -p $(dirname {output.pc_gtf})
        grep -E '^#|gene_biotype "protein_coding"' {input.gtf} \
            > {output.pc_gtf} 2> {log}
        """

rule build_tx2gene_pc:
    input:
        gtf = rules.filter_protein_coding_gtf.output.pc_gtf
    output:
        tx2gene_pc = "data/references/tx2gene_protein_coding.tsv"
    log:
        "logs/dge/build_tx2gene_pc.log"
    message:
        "Building transcript→gene map for protein‑coding transcripts"
    run:
        import pandas as pd
        import re

        # Lire le GTF complet (commentaires ignorés par comment="#")
        df = pd.read_csv(
            input.gtf,
            sep="\t",
            comment="#",
            header=None,
            usecols=[8],  # on ne lit que la colonne des attributs
            names=["attrs"]
        )

        mapping = {}
        for attr in df["attrs"]:
            # si on a à la fois gene_id et transcript_id dans cette ligne
            if "gene_id" in attr and "transcript_id" in attr:
                gid = re.search(r'gene_id "([^"]+)"', attr).group(1)
                tid = re.search(r'transcript_id "([^"]+)"', attr).group(1)
                mapping[tid] = gid

        # On écrit le tsv final
        pd.DataFrame(
            mapping.items(),
            columns=["TXNAME","GENEID"]
        ).to_csv(
            output.tx2gene_pc,
            sep="\t",
            index=False
        )

rule extract_cdn_values:
    """
    Extraire condition1 et condition2 depuis data/comparisons.tsv
    """
    output:
        cdn1_file = "data/cdn1.txt",
        cdn2_file = "data/cdn2.txt"
    run:
        with open("data/comparisons.tsv") as f, \
             open(output.cdn1_file, 'w') as f1, \
             open(output.cdn2_file, 'w') as f2:
            lines = f.read().strip().split("\n")
            # on prend la 2e ligne
            cdn1, cdn2 = lines[1].split()  
            f1.write(cdn1)
            f2.write(cdn2)


rule dge_kallisto_protein_coding:
    """
    Analyse différentielle d’expression (DESeq2) sur les gènes
    protéinocodants, à partir des abundance.h5 de Kallisto.
    """
    input:
        quant      = expand("results/quant/kallisto/{id}/abundance.h5",   id=id_list),
        samples    = "data/design.tsv",
        comparisons= "data/comparisons.tsv",
        tx2gene    = rules.build_tx2gene_pc.output.tx2gene_pc
    output:
        results    = directory("results/dge/kallisto_protein_coding"),
        out_files  = expand(
            "results/dge/kallisto_protein_coding/{comp}_DESeq2_gene.csv",
            comp = comparisons
        )
    params:
        kallisto_dir           = "results/quant/kallisto",
        filter_count_threshold = config["dge"]["filter_count_threshold"]
    log:
        "logs/dge/kallisto_dge_protein_coding.log"
    message:
        "Running gene-level DGE (kallisto, protein-coding) with DESeq2"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_kallisto_protein_coding.R"


rule dge_kallisto_protein_coding_stats:
    """
    Statistiques post‑DESeq2 pour les gènes protéinocodants :
    - calcul d’un tableau récapitulatif
    - tracé de volcano plots (total et zoom)
    - analyses d’enrichissement GO (total, up‑régulés, down‑régulés)
    """
    input:
        # un seul fichier DESeq2 par wildcard {comp}
        deseq2 = "results/dge/kallisto_protein_coding/{comp}_DESeq2_gene.csv"
    output:
        stat           = "results/dge/kallisto_protein_coding/{comp}_deseq2_stats.csv",
        volcano_total  = "results/dge/kallisto_protein_coding/{comp}_volcano_total.png",
        volcano_zoom   = "results/dge/kallisto_protein_coding/{comp}_volcano_zoom.png",
        go_total       = "results/dge/kallisto_protein_coding/{comp}_go_enrichment_total.png",
        go_enrich_up   = "results/dge/kallisto_protein_coding/{comp}_go_enrichment_up.png",
        go_enrich_down = "results/dge/kallisto_protein_coding/{comp}_go_enrichment_down.png"
    log:
        "logs/dge/kallisto_deseq2_stats_{comp}.log"
    message:
        "Calculating post-DESeq2 stats for comparison {wildcards.comp}"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/dge_kallisto_protein_coding_stats.R"

