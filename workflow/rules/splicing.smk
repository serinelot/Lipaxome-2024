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



rule majiq_build:
    input:
        gff3 = config["path"]["human_gff3"],
        bai = expand(rules.bam_index.output, id = id_list)
    output:
        splicegraph = "results/splicing/majiq/build/splicegraph.sql",
        majiq_files = expand('results/splicing/majiq/build/{id}_Aligned.sortedByCoord.out.primary.majiq', id = id_list, allow_missing=True)
    params:
        config = "data/majiq.conf",
        outdir = directory("results/splicing/majiq/build/"),
        majiq = config["tools"]["majiq_voila"]
    log:
        "logs/majiq/build.log"
    message:
        "Analyze RNA-seq data to detect LSV candidates using MAJIQ. "
        "Preinstallation of MAJIQ is required prior to running this rule."
    shell:
        "source {params.majiq} && "
        "majiq build {input.gff3} -c {params.config} -j 8 -o {params.outdir} "
        "&> {log} && deactivate"



rule majiq_deltapsi_quant:
    input:
        build_dir = expand(rules.majiq_build.output.majiq_files, id = id_list),
        comp = lambda wildcards: get_majiq_deltapsi_group_id(wildcards.id)
    output:
        voila = "results/splicing/majiq/deltapsi_quant.deltapsi.voila"
    params:
        outdir = directory("results/splicing/majiq/deltapsi_quant"),
        group = "FXS-Control",
        majiq_tool = config["tools"]["majiq_voila"]
    log:
        "logs/majiq/deltapsi_quant.log"
    message:
        "Quantify differential splicing between two different groups: FXS-Control."
    shell:
        "source {params.majiq_tool} && "
        "majiq deltapsi -grp1 {input.comp[0]} -grp2 {input.comp[1]} "
        "-j 8 -o {params.outdir} --name {params.group} "
        "&> {log} && deactivate"

