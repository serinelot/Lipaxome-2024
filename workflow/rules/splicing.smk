rule primary_alignments:
    # More robust results with multimapped reads
    input:
        rules.star_alignReads.output.bam
    output:
        "results/splicing/star/{id}/{id}_Aligned.sortedByCoord.out.primary.bam"
    log:
        "logs/STAR/{id}_primary.log"
    conda:
        "../envs/genomecov.yml"
    message:
        "Keep primary alignments only for {wildcards.id}."
    shell:
        "samtools view -b -F 256 -o {output} {input} "
        "&> {log}"


rule bam_index:
    input:
        rules.primary_alignments.output
    output:
        "results/splicing/star/{id}/{id}_Aligned.sortedByCoord.out.primary.bam.bai"
    log:
        "logs/star/{id}_primary_index.log"
    conda:
        "../envs/genomecov.yml"
    message:
        "Create a BAI index for {wildcards.id}."
    shell:
        "samtools index {input} "
        "&> {log}"


rule genomecov:
    input:
        rules.primary_alignments.output
    output:
        "results/splicing/genomecov/{id}.bedgraph"
    conda:
        "../envs/genomecov.yml"
    message:
        "Report {wildcards.id} genome coverage in BEDGRAPH format."
    shell:
        "bedtools genomecov -bg -split -ibam {input} | sort -k1,1 -k2,2n > {output}"


rule split_bams_files:
    input:
        design="data/design.tsv"
    output:
        "scripts/split_bams_files.sh"
    shell:
        """
        python scripts/split_bams_by_cond.py
        """

rule execute_split_bams_files:
    input:
        "scripts/split_bams_files.sh"
    output:
        "results/splicing/star/split_bam_done.txt"
    shell:
        """
        bash {input}
        touch results/splicing/star/split_bam_done.txt
        """


rule create_bam_lists_for_rmats:
    output:
        fxs_bam_list = "results/splicing/rmats/fxs_bam_list.txt",
        control_bam_list = "results/splicing/rmats/control_bam_list.txt"
    shell:
        """
        find results/splicing/star_FXS/ -name '*_Aligned.sortedByCoord.out.primary.bam' -print0 | sort -z | xargs -0 printf '%s,' | sed 's/,$/\\n/' > {output.fxs_bam_list}
        find results/splicing/star_Control/ -name '*_Aligned.sortedByCoord.out.primary.bam' -print0 | sort -z | xargs -0 printf '%s,' | sed 's/,$/\\n/' > {output.control_bam_list}
        """


rule run_rmats:
    input:
        bam = expand(rules.primary_alignments.output, id = id_list),
        group1 = "results/splicing/rmats/fxs_bam_list.txt",
        group2 = "results/splicing/rmats/control_bam_list.txt",
        gtf = config['download']['human_gtf']
    output:
        outdir = directory("results/splicing/rmats/{comp}/raw"),
        tmpdir = directory("results/splicing/rmats/{comp}/tmp"),
        summary = "results/splicing/rmats/{comp}/raw/summary.txt"
    params:
        readlength = 80
    conda:
        "../envs/rmats.yml"
    log:
        "logs/rmats/{comp}.log"
    message:
        "Run rMATS for {wildcards.comp}."
    shell:
        "rmats.py --b1 {input.group1} --b2 {input.group2} "
        "--gtf {input.gtf} -t paired --readLength {params.readlength} --variable-read-length "
        "--nthread 4 --od {output.outdir} --tmp {output.tmpdir}"
        "&> {log}"


rule filter_rmats:
    input:
        summary = rules.run_rmats.output.summary
    output:
        result = 'results/splicing/rmats/{comp}/filtered/SE.tsv'
    params:
        dir = directory("results/splicing/rmats/{comp}"),
        tpm = rules.merge_kallisto_quant.output.tpm,
        gtf = config['download']['human_gtf'],
        fdr = 0.05,
        deltapsi = 0.10,
        min_tpm = 0.75
    conda:
        "../envs/python.yml"
    log:
        "logs/rmats/filter_{comp}.log"
    message:
        "Filter raw rMATS output for {wildcards.comp}."
    script:
        "../scripts/filter_rmats.py"
