rule build_transcriptome:
    """
    Génère le FASTA du transcriptome à partir du GTF + génome.
    Utile pour fournir à kallisto les séquences de transcriptes.
    """
    input:
        genome = rules.download_human_genome.output.genome,
        gtf    = config["download"]["human_gtf"]
    output:
        tx_fa  = config["path"]["transcriptome"]
    conda:
        "../envs/gffread.yml"
    message:
        "Building reference transcriptome with gffread"
    log:
        "logs/build_transcriptome/build_transcriptome.log"
    shell:
        r"""
        mkdir -p $(dirname {output.tx_fa})
        gffread {input.gtf} \
          -g {input.genome} \
          -w {output.tx_fa} \
        &> {log}
        """


rule build_kallisto_index:
    """
    Construit l'index kallisto à partir du transcriptome FASTA généré.
    """
    input:
        transcriptome = rules.build_transcriptome.output.tx_fa
    output:
        index = "data/references/kallisto.idx"
    message:
        "Building kallisto index from {input.transcriptome}"
    conda:
        "../envs/kallisto.yml"
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.index})
        kallisto index \
          -i {output.index} \
          {input.transcriptome}
        """

rule quant_kallisto:
    """
    Estime les abondances des transcrits par pseudo‑alignement avec kallisto,
    produit abundance.tsv pour chaque échantillon.
    """
    input:
        idx  = rules.build_kallisto_index.output.index,
        fq1  = "results/preprocess/fastp/{id}_1_trimmed.fastq.gz",
        fq2  = "results/preprocess/fastp/{id}_2_trimmed.fastq.gz"
    output:
        tsv  = "results/quant/kallisto/{id}/abundance.tsv"
    params:
        outdir = "results/quant/kallisto/{id}"
    threads: 8
    message: "Quantifying sample {wildcards.id} with kallisto"
    conda:
        "../envs/kallisto.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        kallisto quant \
          -i {input.idx} \
          -o {params.outdir} \
          -t {threads} \
          {input.fq1} {input.fq2}
        """

