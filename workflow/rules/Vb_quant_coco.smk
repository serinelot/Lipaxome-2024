rule coco_correct_annotation:
    """
    Corrige la GTF de référence avec COCO pour éliminer les chevauchements
    et préparer l’annotation pour la quantification.
    """
    input:
        gtf      = config["download"]["human_gtf"],
        coco_dir = rules.clone_coco.output.repo_dir
    output:
        gtf_corr = "data/references/gtf/Homo_sapiens.GRCh38.110_snoRNAs_tRNAs_correct_annotation.gtf"
    message: "Correcting GTF annotation by COCO"
    conda:
        "../envs/coco.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.gtf_corr})
        python {input.coco_dir}/bin/coco.py ca \
          {input.gtf} \
          -o {output.gtf_corr}
        """

rule coco_quant:
    """
    Quantifie les reads alignés (BAM STAR) avec COCO,
    produit counts, CPM et TPM pour chaque échantillon.
    """
    input:
        gtf_corr = "data/references/gtf/Homo_sapiens.GRCh38.110_snoRNAs_tRNAs_correct_annotation.gtf",
        bam      = "results/alignment/star/{id}/Aligned.sortedByCoord.out.bam"
    output:
        counts   = "results/quant/coco/{id}.tsv"
    params:
        coco_script = rules.clone_coco.output.repo_dir + "/bin/coco.py"
    threads: 32
    message: "Quantifying sample {wildcards.id} with COCO"
    conda:
        "../envs/coco.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.counts})
        python {params.coco_script} cc \
          --countType both \
          --thread {threads} \
          --strand 1 \
          --paired \
          {input.gtf_corr} \
          {input.bam} \
          {output.counts}
        """

rule merge_coco_cc:
    """
    Fusionne les résultats CoCo (count, CPM, TPM) de tous les échantillons
    en trois matrices consolidées.
    """
    input:
        coco_tsv = expand("results/quant/coco/{id}.tsv", id=id_list)
    output:
        merged_counts = "results/quant/coco_merged_counts.tsv",
        merged_cpm    = "results/quant/coco_merged_cpm.tsv",
        merged_tpm    = "results/quant/coco_merged_tpm.tsv"
    conda:
        "../envs/coco.yml"
    message:
        "Merging CoCo quantification outputs (count / CPM / TPM) across all samples"
    script:
        "../scripts/merge_coco_quant.py"
