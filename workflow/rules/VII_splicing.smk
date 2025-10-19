rule extract_primary_alignments:
    """
    Conserver uniquement les alignements primaires (flag -F 256) pour l’échantillon {wildcards.id}.
    """
    input:
        bam = rules.star_alignReads.output.bam
    output:
        primary_bam = "results/splicing/star/{id}/{id}_Aligned.sortedByCoord.out.primary.bam"
    log:
        "logs/splicing/star/{id}_primary.log"
    conda:
        "../envs/genomecov.yml"
    message:
        "Keep primary alignments only for sample {wildcards.id}"
    shell:
        r"""
        mkdir -p $(dirname {output.primary_bam})
        samtools view -b -F 256 {input.bam} > {output.primary_bam} 2> {log}
        """

rule index_primary_bam:
    """
    Créer l'index BAI pour le BAM primaire de l'échantillon {wildcards.id}.
    """
    input:
        primary_bam = rules.extract_primary_alignments.output.primary_bam
    output:
        bai = "results/splicing/star/{id}/{id}_Aligned.sortedByCoord.out.primary.bam.bai"
    log:
        "logs/splicing/star/{id}_primary_index.log"
    conda:
        "../envs/genomecov.yml"
    message:
        "Indexation du BAM primaire pour l'échantillon {wildcards.id}"
    shell:
        """
        samtools index {input.primary_bam} &> {log}
        """

rule coverage_bedgraph:
    """
    Calculer la couverture génomique en BedGraph pour {wildcards.id}.
    """
    input:
        primary_bam = rules.extract_primary_alignments.output.primary_bam
    output:
        bedgraph = "results/splicing/genomecov/{id}.bedgraph"
    conda:
        "../envs/genomecov.yml"
    message:
        "Génération du BedGraph de couverture pour {wildcards.id}"
    shell:
        """
        bedtools genomecov -bg -split -ibam {input.primary_bam} \
        | sort -k1,1 -k2,2n > {output.bedgraph}
        """

rule bedgraph_to_bigwig:
    """
    Convertit le fichier bedGraph en bigWig, basé sur les tailles chromosomiques.
    """
    input:
        bedgraph = "results/splicing/genomecov/{id}.bedgraph",
        chrom_sizes = rules.fai_to_chromsizes.output.chrom_sizes
    output:
        bigwig = "results/splicing/genomecov/{id}.bw"
    conda:
        "../envs/genomecov.yml"
    message:
        "Conversion BedGraph -> BigWig pour {wildcards.id}"
    shell:
        """
        bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bigwig}
        """

rule index_gtf:
    input:
        gtf="data/references/gtf/Homo_sapiens.GRCh38.115.gtf"
    output:
        gtf_gz="data/references/gtf/Homo_sapiens.GRCh38.115.gtf.sorted.gz",
        tbi="data/references/gtf/Homo_sapiens.GRCh38.115.gtf.sorted.gz.tbi"
    conda:
        "../envs/genomecov.yml"
    shell:
        """
        sort -k1,1 -k4,4n {input.gtf} > {input.gtf}.sorted
        bgzip -c {input.gtf}.sorted > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        rm {input.gtf}.sorted
        """



rule make_split_script:
    """
    Générer le script Bash pour séparer les BAM par condition.
    """
    input:
        design = "data/design.tsv"
    output:
        script = "scripts/split_bams_by_cond.sh"
    message:
        "Création du script de séparation des BAM par condition"
    shell:
        """
        python scripts/split_bams_by_cond.py
        chmod +x {output.script}
        """

rule run_split_script:
    """
    Exécuter le script de séparation des BAM par condition.
    """
    input:
        script = rules.make_split_script.output.script
    output:
        done = "results/splicing/star/split_bam_done.txt"
    message:
        "Séparation des BAM par condition terminée"
    shell:
        """
        bash {input.script}
        touch {output.done}
        """

rule list_bams_for_rmats:
    """
    Générer les listes de BAM (FXS & Control) pour rMATS.
    """
    input:
        done = rules.run_split_script.output.done
    output:
        fxs_list     = "results/splicing/rmats/fxs_bam_list.txt",
        control_list = "results/splicing/rmats/control_bam_list.txt"
    message:
        "Création des listes de BAM pour rMATS"
    shell:
        """
        find results/splicing/star/FXS/ -name '*primary.bam' \
            | sort | paste -sd, - > {output.fxs_list}
        find results/splicing/star/Control/ -name '*primary.bam' \
            | sort | paste -sd, - > {output.control_list}
        """


rule merge_kallisto_tpm:
    """
    Fusionne les fichiers abundance.tsv de tous les échantillons
    en une matrice TPM pour les gènes protéinocodants.
    """
    input:
        abundances = expand("results/quant/kallisto/{id}/abundance.tsv", id=id_list),
        gtf        = rules.filter_protein_coding_gtf.output.pc_gtf,
        tx2gene    = rules.build_tx2gene_pc.output.tx2gene_pc
    output:
        tpm_matrix = "results/quant/kallisto_merged_tpm.tsv"
    log:
        "logs/kallisto/merge_tpm.log"
    message:
        "Merging kallisto abundance files into TPM matrix"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/merge_kallisto_tpm_quant.py"


rule rmats_original:
    """
    Exécuter rMATS (GTF original) pour la comparaison {wildcards.comp}.
    """
    input:
        bams    = expand(rules.extract_primary_alignments.output.primary_bam, id=id_list),
        group1  = rules.list_bams_for_rmats.output.fxs_list,
        group2  = rules.list_bams_for_rmats.output.control_list,
        gtf     = config["download"]["human_gtf"]
    output:
        raw_dir = directory("results/splicing/rmats/{comp}/raw"),
        tmp_dir = directory("results/splicing/rmats/{comp}/tmp"),
        summary = "results/splicing/rmats/{comp}/raw/summary.txt"
    params:
        readlength = 80
    log:
        "logs/splicing/rmats/{comp}.log"
    conda:
        "../envs/rmats.yml"
    message:
        "Exécution de rMATS (GTF original) pour {wildcards.comp}"
    shell:
        """
        rmats.py \
          --b1 {input.group1} --b2 {input.group2} \
          --gtf {input.gtf} -t paired \
          --readLength {params.readlength} \
          --variable-read-length \
          --nthread 4 \
          --od {output.raw_dir} --tmp {output.tmp_dir} \
          &> {log}
        """

rule rmats_coco:
    """
    Exécuter rMATS (GTF CoCo‑corrigée) pour la comparaison {wildcards.comp}.
    """
    input:
        bams    = expand(rules.extract_primary_alignments.output.primary_bam, id=id_list),
        group1  = rules.list_bams_for_rmats.output.fxs_list,
        group2  = rules.list_bams_for_rmats.output.control_list,
        gtf_corr= rules.filter_protein_coding_gtf_coco.output.pc_gtf_corr
    output:
        raw_dir = directory("results/splicing/rmats_coco/{comp}/raw"),
        tmp_dir = directory("results/splicing/rmats_coco/{comp}/tmp"),
        summary = "results/splicing/rmats_coco/{comp}/raw/summary.txt"
    params:
        readlength = 80
    log:
        "logs/splicing/rmats_coco/{comp}.log"
    conda:
        "../envs/rmats.yml"
    message:
        "Exécution de rMATS (CoCo‑corrigée) pour {wildcards.comp}"
    shell:
        """
        rmats.py \
          --b1 {input.group1} --b2 {input.group2} \
          --gtf {input.gtf_corr} -t paired \
          --readLength {params.readlength} \
          --variable-read-length \
          --nthread 4 \
          --od {output.raw_dir} --tmp {output.tmp_dir} \
          &> {log}
        """

rule filter_rmats_original:
    """
    Filtrer la sortie brute de rMATS (GTF original) pour {wildcards.comp}.
    """
    input:
        summary    = rules.rmats_original.output.summary,
        tpm_matrix = rules.merge_kallisto_tpm.output.tpm_matrix,
        gtf        = config["download"]["human_gtf"]
    output:
        filtered   = "results/splicing/rmats/{comp}/filtered/SE.tsv"
    params:
        fdr     = 0.05,
        dpsi    = 0.10,
        min_tpm = 0.75
    log:
        "logs/splicing/rmats/filter_{comp}.log"
    conda:
        "../envs/python.yml"
    message:
        "Filtrage de rMATS (GTF original) pour {wildcards.comp}"
    script:
        "../scripts/filter_rmats.py"


rule filter_rmats_coco:
    """
    Filtrage de la sortie brute de rMATS (GTF CoCo‑corrigée) pour {wildcards.comp}.
    """
    input:
        summary     = rules.rmats_coco.output.summary,
        tpm_matrix  = rules.merge_coco_cc.output.merged_tpm,
        gtf_corr    = rules.filter_protein_coding_gtf_coco.output.pc_gtf_corr
    output:
        filtered    = "results/splicing/rmats_coco/{comp}/filtered/SE.tsv"
    params:
        fdr         = 0.05,
        dpsi        = 0.10,
        min_tpm     = 0.75
    log:
        "logs/splicing/rmats_coco/filter_{comp}.log"
    conda:
        "../envs/python.yml"
    message:
        "Filtrage de rMATS (CoCo‑corrigée) pour {wildcards.comp}"
    script:
        "../scripts/filter_rmats_coco.py"
